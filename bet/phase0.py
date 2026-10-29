"""Phase 0: per-query magic demand profiles of routed block encodings.

Encodings (Pauli units, H = sum_b J_b (XX + YY + ZZ)):
  square J1-J2, J2/J1 = 0.5, periodic L x L
    swap     routed SWAP-LCU: per bond J (2 SWAP - 1); target = 1 CSWAP; 4 swap trees; 6 fix-up CSWAPs
    triangle routed spectral triangle (lcvcbe clusters_v2 'tri_0.5'); 4 trees; 6 fix-up CSWAPs
  kagome Heisenberg, periodic L x L cells (N = 3 L^2)
    swap     routed SWAP-LCU on the 6 bonds per cell; 3 trees; 4 fix-up CSWAPs (lcvcbe fixup.json 'kagome bond')
    triangle routed up/down triangles (isotropic, the same normalized block as 'tri_0.5'); 3 trees, no fix-up
             (routing.kagome_check), addressing + 2 controlled modular decrements of the address

Demand categories per query: routing CSWAPs (tree, fix-up, layer structure), addressing ANDs, PREPARE
(and its inverse), the local target, and the rotation angles of PREPARE and the local block. Rotation
precision is not fixed here (Phase 2); `lcvcbe_rule` reproduces lcvcbe routed_cost2's T per query as a check.

Query counts for an absolute energy error eps (intensive task): N_q = ceil(pi X / (2 eps)) with X = alpha
(worst case) or X = lambda_eff = sqrt(delta (2 alpha - delta)), delta = E0 - (c - alpha).

Usage: python -I -m bet.phase0  ->  results/phase0_profiles.{json,md}
"""

import json
import math
from pathlib import Path

from bet import lcvcbe_io as IO
from bet import routing as R

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"

# ground-state energies per site, Pauli units (4x the S.S convention)
E0_SITE = {
    "square J1-J2 0.5": -1.9872,  # SU(2) DMRG 2D limit, Gong et al. PRL 113, 027201 (2014), as in IonQ Table I
    "kagome": -1.7544,  # 4 x (-0.4386), Depenbrock, McCulloch, Schollwoeck PRL 109, 067201 (2012)
}
EPS = [1e-2, 1e-3]
SQUARE_L = [8, 10, 16, 20, 30, 32, 64]
KAGOME_L = [4, 6, 8, 10, 16, 18]


def qcounts(bloq):
    from qualtran.resource_counting import QECGatesCost, get_cost_value

    g = get_cost_value(bloq, QECGatesCost())
    return {"and": int(g.and_bloq), "tof": int(g.toffoli), "cswap": int(g.cswap), "t": int(g.t)}


def uniform(n):
    """PREPARE of a uniform superposition over n items: counts and rotation angles."""
    if n <= 1:
        return {}, []
    from qualtran.bloqs.state_preparation import PrepareUniformSuperposition

    b = PrepareUniformSuperposition(n)
    return qcounts(b), IO.synth().bloq_angles(b)


def add(*parts):
    c, a = {}, []
    for pc, pa in parts:
        for k, v in pc.items():
            c[k] = c.get(k, 0) + v
        a = a + list(pa)
    return c, a


def local_block(key):
    cl = IO.clusters()[key]
    circ = IO.block_circuit(cl)[0]
    c = {k: v for k, v in cl["counts"].items() if k not in ("clifford", "rot")}
    ang = IO.synth().circuit_angles(circ)
    assert len(ang) == cl["counts"]["rot"]
    return c, ang, cl


def lambda_eff(alpha, center, e0):
    delta = e0 - (center - alpha)
    return math.sqrt(max(delta, 0.0) * (2 * alpha - delta)), delta


def lcvcbe_rule(counts, angles, alpha, eps, budget):
    """lcvcbe routed_cost2 T per query (gridsynth, precision budget * eps / (alpha R) per rotation)."""
    S = IO.synth()
    Rn = len(angles)
    er = budget * eps / (alpha * Rn) if Rn else 1.0
    trot = sum(S.tcount(t, er) for t in angles)
    return counts.get("t", 0) + 4 * (counts.get("and", 0) + counts.get("tof", 0) + counts.get("cswap", 0)) + trot


def profile(lattice, method, L):
    fx = IO.fixups()
    if lattice == "square J1-J2 0.5":
        N, chi = L * L, 4
        e0 = E0_SITE[lattice] * N
        addressing = {"and": 4 * math.ceil(math.log2(L))}
        if method == "swap":
            fix = fx["square J1-J2 bond"]
            p = 0.5 / 1.5  # NNN / total bond weight
            lab = ({}, [2 * math.asin(math.sqrt(p))])  # one rotation splits NN / NNN; directions by H
            pos = uniform(N)
            # as lcvcbe routed_cost2: the cheaper of (position x label) and a joint uniform over N x 6 entries
            # (items 2, 2, 1, 1 -> 6 entries, decoded with 2 ceil(log2 6) ANDs)
            opts = [add(pos, pos, lab, lab), add(uniform(6 * N), uniform(6 * N), ({"and": 2 * 3}, []))]
            prep = min(opts, key=lambda o: lcvcbe_rule(o[0], o[1], 6.0 * N, 1e-3, 0.1))
            local = ({"cswap": 1}, [])
            alpha, center, budget = 6.0 * N, -3.0 * N, 0.1
        else:
            fix = fx["square cluster"]
            prep = add(uniform(4 * N), uniform(4 * N))
            c, ang, cl = local_block("tri_0.5")
            local = (c, ang)
            alpha, center, budget = 4 * N * cl["hw"], 4 * N * cl["shift"], 0.05
    else:  # kagome
        N, chi = 3 * L * L, 3
        e0 = E0_SITE["kagome"] * N
        nbits = math.ceil(math.log2(L))
        if method == "swap":
            fix = fx["kagome bond"]
            addressing = {"and": 4 * nbits}
            prep = add(uniform(6 * L * L), uniform(6 * L * L))
            local = ({"cswap": 1}, [])
            alpha, center, budget = 4.0 * N, -2.0 * N, 0.1  # 2N bonds, 2J each; shift -J per bond
        else:
            fix = 0
            # address register + two controlled modular decrements (compute and uncompute) of nbits each
            addressing = {"and": 4 * nbits + 2 * 2 * nbits}
            prep = add(uniform(L * L), uniform(L * L))  # label up/down by H
            c, ang, cl = local_block("tri_0.5")
            local = (c, ang)
            alpha, center, budget = 2.0 * N, 0.0, 0.05  # 2N/3 triangles, half-width 3

    rt = R.routing_profile(N, chi, fix, power_of_two=(N // chi) & (N // chi - 1) == 0)
    routing_counts = {"cswap": rt["cswap_tree"] + rt["cswap_fixup"]}
    total = add((routing_counts, []), (addressing, []), prep, local)
    leff, delta = lambda_eff(alpha, center, e0)
    row = {
        "lattice": lattice,
        "method": method,
        "L": L,
        "N": N,
        "alpha": alpha,
        "center": center,
        "E0": e0,
        "delta": delta,
        "lambda_eff": leff,
        "routing": rt,
        "addressing": addressing,
        "prep": {"counts": prep[0], "angles": prep[1]},
        "local": {"counts": local[0], "angles": local[1]},
        "total": {"counts": total[0], "n_rotations": len(total[1])},
        "queries": {
            f"{eps:g}": {"alpha": math.ceil(math.pi * alpha / (2 * eps)), "lambda_eff": math.ceil(math.pi * leff / (2 * eps))}
            for eps in EPS
        },
        "check_T_query_lcvcbe_rule_eps1e-3": lcvcbe_rule(total[0], total[1], alpha, 1e-3, budget),
    }
    return row


def main():
    OUT.mkdir(exist_ok=True)
    rows = [profile("square J1-J2 0.5", m, L) for m in ("swap", "triangle") for L in SQUARE_L]
    rows += [profile("kagome", m, L) for m in ("swap", "triangle") for L in KAGOME_L]
    (OUT / "phase0_profiles.json").write_text(json.dumps(rows, indent=1))

    # cross-check against lcvcbe results/routed_cost2.json (square, gridsynth, eps 1e-3)
    ref = json.loads((IO.LCVCBE / "results" / "routed_cost2.json").read_text())
    lines = ["# Phase 0: per-query magic demand profiles", ""]
    lines.append("Counts per query (both routing directions, PREPARE and its inverse). CSWAP = 1 Toffoli-class gate.")
    lines.append("")
    lines.append("| lattice | method | L | N | CSWAP routing | AND addr | PREP (AND, T, rot) | local (AND+Tof, T, rot) | "
                 "routing share of non-Clifford gates | alpha/N | lambda_eff/N | T/query (lcvcbe rule) | lcvcbe ref |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        pc, lc = r["prep"]["counts"], r["local"]["counts"]
        nonclif = r["total"]["counts"].get("cswap", 0) + r["total"]["counts"].get("and", 0) + r["total"]["counts"].get("tof", 0)
        nonclif += r["total"]["counts"].get("t", 0) + r["total"]["n_rotations"]
        share = (r["routing"]["cswap_tree"] + r["routing"]["cswap_fixup"]) / nonclif
        refv = ""
        if r["lattice"].startswith("square"):
            m = {"swap": "swap", "triangle": "cluster"}[r["method"]]
            hit = [x for x in ref if x["model"] == "J1-J2 J2=0.5" and x["L"] == r["L"] and x["eps"] == 1e-3 and x["mode"] == "gridsynth"]
            refv = str(hit[0][m]["T_query"]) if hit else "-"
        lines.append(
            f"| {r['lattice']} | {r['method']} | {r['L']} | {r['N']} | {r['routing']['cswap_tree'] + r['routing']['cswap_fixup']} | "
            f"{r['addressing']['and']} | {pc.get('and', 0)}, {pc.get('t', 0)}, {len(r['prep']['angles'])} | "
            f"{lc.get('and', 0) + lc.get('tof', 0) + lc.get('cswap', 0)}, {lc.get('t', 0)}, {len(r['local']['angles'])} | "
            f"{share:.3f} | {r['alpha'] / r['N']:.3f} | {r['lambda_eff'] / r['N']:.3f} | {r['check_T_query_lcvcbe_rule_eps1e-3']} | {refv} |"
        )
    lines += ["", "Routing layers (power-of-two trees): each layer j of a tree holds 2^j CSWAPs sharing one address bit;",
              "all trees run in parallel. Example rows:"]
    for r in rows:
        if r["L"] in (16, 8) and "layers_per_direction" in r["routing"]:
            lay = ", ".join(str(x["width"]) for x in r["routing"]["layers_per_direction"])
            lines.append(f"- {r['lattice']} {r['method']} L={r['L']}: widths per direction [{lay}], "
                         f"{r['routing']['depth_cswap_layers']} CSWAP layers per query")
    (OUT / "phase0_profiles.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

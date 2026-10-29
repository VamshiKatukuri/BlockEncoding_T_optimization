"""Phase 5: end-to-end comparison in factory space-time (physical qubitcycles per run).

Task: ground-state energy to an absolute error eps (Pauli units), lambda_eff query counts, failure budget 0.1.
Methods:
  walk triangle  routed spectral triangle walk;  walk swap  routed SWAP-LCU walk
  trotter        2nd- or 4th-order Trotter phase estimation, xyz fragments, Hamming-weight phasing (lcvcbe
                 trotter/trotter_cost.py, empirical kappa from <= 24 sites, read only); cheaper of chunk 31 and N/2
Pricing:
  T4-flat   every Toffoli-class gate = 4 T, rotations at 1e-10 (gridsynth); T from the cheapest adequate factory
  T4        as T4-flat with the first-order rotation rule (lcvcbe convention)
  BEST      CCZ factory or 4 T per Toffoli-class gate, catalysed CCZ->2T or T factory for T states; second-order
            rotation rule for the triangle walk (Phase 2), first-order otherwise
Wall clock: run volume / (factory qubits) cycles, at 1e6 physical qubits of factories and 1 us per cycle; a
lower bound when the factories are the bottleneck (Phase 1 and 4 show they are for routing).

Usage: python -m bet.compare  ->  results/phase5_compare.{json,md}
"""

import json
import sys
from pathlib import Path

from bet import factories as F
from bet import lcvcbe_io as IO

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"
FACTORY_QUBITS, CYCLE_S = 1e6, 1e-6


def trotter_module():
    p = IO.LCVCBE
    sys.path[:0] = [str(p), str(p / "trotter"), str(p / "simplex")]
    import trotter_cost as TC

    return TC


def run_price(n3, n1, policy, queries, p=1e-3):
    v, fac = F.price(n3, n1, policy, p, queries=queries)
    return None if v is None else v * queries, fac


def main():
    S = IO.synth()
    TC = trotter_module()
    kap = TC.load_kappa()
    rows = json.loads((OUT / "phase0_profiles.json").read_text())
    fid = {(o["lattice"], o["method"], o["L"], o["eps"]): o for o in json.loads((OUT / "phase2_fidelity.json").read_text())}
    res = []
    targets = [("kagome", L) for L in (4, 6, 8, 10, 16, 18)] + [("square J1-J2 0.5", L) for L in (8, 10, 16, 20, 30, 32)]
    for lat, L in targets:
        for eps in (1e-2, 1e-3):
            out = {"lattice": lat, "L": L, "eps": eps}
            for meth in ("triangle", "swap"):
                r = next(x for x in rows if x["lattice"] == lat and x["method"] == meth and x["L"] == L)
                o = fid[(lat, meth, L, eps)]
                out["N"] = r["N"]
                c = r["total"]["counts"]
                n3 = c.get("and", 0) + c.get("tof", 0) + c.get("cswap", 0)
                nq = o["queries"]
                tA = o["rule_A"]["T_rot_gridsynth"] if "rule_A" in o else 0
                tB = o["rule_B"]["T_rot_gridsynth"] if "rule_B" in o else 0
                tC = o["rule_C"]["T_rot_gridsynth"] if "rule_C" in o else tB
                for name, pol, trot in (("T4-flat", "T4", tA), ("T4", "T4", tB), ("BEST", "BEST", tC)):
                    v, fac = run_price(n3, c.get("t", 0) + trot, pol, nq)
                    out[f"walk {meth} {name}"] = v
            N = out["N"]
            bps = 2 if lat == "kagome" else 4
            case = "kagome" if lat == "kagome" else "square J1-J2 J2=0.5"
            best = {}
            for order in (2, 4):
                if (case, order) not in kap:
                    continue
                for chunk in (31, max(1, N // 2)):
                    t = TC.trotter(kap[(case, order)][0], order, N, bps, eps, chunk)
                    n3t, n1t = t["toffoli"], t["t"]
                    for name in ("T4", "BEST"):
                        v, _ = run_price(n3t, n1t, name, 1)
                        if v is not None and (name not in best or v < best[name][0]):
                            best[name] = (v, order, chunk)
            for name, (v, order, chunk) in best.items():
                out[f"trotter {name}"] = v
                out[f"trotter {name} variant"] = f"order {order}, chunk {chunk}"
            res.append(out)
    (OUT / "phase5_compare.json").write_text(json.dumps(res, indent=1))

    def h(v):
        return "n/a" if v is None else f"{v:.2e} ({v / FACTORY_QUBITS * CYCLE_S / 3600:.2g} h)"

    L = ["# Phase 5: end-to-end comparison in factory space-time", "",
         f"Physical qubitcycles per run (p = 1e-3, Litinski factories); in brackets the wall clock with {FACTORY_QUBITS:.0e} "
         "physical qubits of factories at 1 us per cycle (factory-bound lower bound).", "",
         "| lattice | L | N | eps | walk tri T4-flat | walk tri T4 | walk tri BEST | walk swap T4 | walk swap BEST | Trotter T4 | Trotter BEST | "
         "tri BEST gain vs T4-flat | best walk / best Trotter |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for o in res:
        g = o["walk triangle T4-flat"] / o["walk triangle BEST"] if o.get("walk triangle BEST") and o.get("walk triangle T4-flat") else None
        walks = [o[k] for k in ("walk triangle BEST", "walk swap BEST") if o.get(k)]
        tro = o.get("trotter BEST")
        ratio = min(walks) / tro if walks and tro else None
        L.append(f"| {o['lattice']} | {o['L']} | {o['N']} | {o['eps']:g} | {h(o.get('walk triangle T4-flat'))} | {h(o.get('walk triangle T4'))} | "
                 f"{h(o.get('walk triangle BEST'))} | {h(o.get('walk swap T4'))} | {h(o.get('walk swap BEST'))} | {h(o.get('trotter T4'))} | "
                 f"{h(tro)} | {'-' if g is None else f'{g:.2f}'} | {'-' if ratio is None else f'{ratio:.2f}'} |")
    (OUT / "phase5_compare.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

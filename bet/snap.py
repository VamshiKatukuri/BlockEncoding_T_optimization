"""Phase 3: patch designs matched to magic-state supply.

1. Exactness obstruction. A Clifford+T circuit (with ancillas prepared and post-selected in |0>) has matrix
   entries in Z[1/sqrt2, i] (Giles & Selinger, PRA 87, 032332 (2013)), so the block it encodes has entries in
   that ring. The optimal-normalization patch block A = (h - delta)/lambda of an isotropic triangle (lambda = 3,
   Pauli units) is the reflection 2 P(S=3/2) - 1, with entries -1/3 and 2/3: no exact Clifford+T circuit encodes
   it, so its basis change needs synthesized rotations at any precision.
2. Ring-compatible normalization. A/lambda' with integer h needs lambda' in {2^(k/2)} (times units). The smallest
   admissible lambda' >= half-width, with the patch re-centred so its lower edge stays put (the lower edge sets
   lambda_eff), gives a lambda_eff penalty computed below. This is a necessary condition only: an exact
   dilation at lambda' is not constructed here.
3. Synthillation upper bound. Distilling the fixed local block as one gate can at most remove its share of the
   per-query cost: bound 1 / (1 - local share).
4. Nested reflections. Any patch with class eigenvalues e_0 < ... < e_r equals c0 + sum_k b_k (2 P_{>=k} - 1),
   b_k = (e_k - e_{k-1}) / 2, sum b_k = half-width: optimal normalization with no eigenvalue rotations, but the
   weights b_k move into PREPARE (r - 1 rotations each way), so the rotation count does not drop.

Usage: python -m bet.snap  ->  results/phase3_snap.{json,md}
"""

import json
import math
from pathlib import Path

import numpy as np

from bet import factories as F
from bet import lcvcbe_io as IO
from bet.bias import reference_triangle

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"


def ring_lambda(half_width):
    """Smallest 2^(k/2) >= half_width (k integer)."""
    k = math.ceil(2 * math.log2(half_width) - 1e-12)
    return 2 ** (k / 2)


def lam_eff(alpha, delta):
    return math.sqrt(max(delta, 0) * (2 * alpha - delta))


def tetra_classes():
    """Isotropic tetrahedron sum s_i.s_j over 6 edges (Pauli units): S = 0, 1, 2 -> -6, -2, 6."""
    return [-6.0, -2.0, 6.0]


def main():
    rows = json.loads((OUT / "phase0_profiles.json").read_text())
    fid = {(o["lattice"], o["method"], o["L"], o["eps"]): o for o in json.loads((OUT / "phase2_fidelity.json").read_text())}
    out = {"obstruction": {}, "ring": [], "synthillation": [], "nested": []}

    # 1. obstruction: entries of the optimal reflection block of the isotropic triangle
    h = reference_triangle(1.0, 1.0).real
    A = h / 3.0
    ent = sorted({round(float(v), 6) for v in A.flatten()})
    out["obstruction"] = {"block": "isotropic triangle, lambda = 3", "distinct_entries": ent,
                          "in_Z[1/sqrt2,i]": False, "reason": "entries 1/3 and 2/3 are not dyadic"}

    # 2. ring-compatible normalization and its lambda_eff penalty (triangle encodings)
    for r in rows:
        if r["method"] != "triangle":
            continue
        n_p = 4 * r["N"] if r["lattice"].startswith("square") else 2 * r["N"] // 3
        hw = r["alpha"] / n_p
        scale = 4.0 if r["lattice"].startswith("square") else 1.0  # integer Pauli form: 4 h (square, weights 1/4) or h
        lam_int = hw * scale
        lam_ring = ring_lambda(lam_int) / scale
        alpha2 = n_p * lam_ring
        # keep the lower edge: center moves up by (alpha2 - alpha); delta unchanged
        le0, le2 = r["lambda_eff"], lam_eff(alpha2, r["delta"])
        le_nocenter = lam_eff(alpha2, r["delta"] + (alpha2 - r["alpha"]))
        out["ring"].append({"lattice": r["lattice"], "L": r["L"], "N": r["N"], "half_width": hw, "ring_lambda": lam_ring,
                            "alpha_penalty": lam_ring / hw, "lambda_eff_penalty_recentred": le2 / le0,
                            "lambda_eff_penalty_same_center": le_nocenter / le0})

    # 3. synthillation upper bound and 2.+3. break-even with the rotation share (BEST pricing, p = 1e-3, rule C)
    for r in rows:
        for eps in (1e-2, 1e-3):
            o = fid[(r["lattice"], r["method"], r["L"], eps)]
            c, lc = r["total"]["counts"], r["local"]["counts"]
            n3 = c.get("and", 0) + c.get("tof", 0) + c.get("cswap", 0)
            rule = "rule_C" if "rule_C" in o else ("rule_B" if "rule_B" in o else None)
            t_rot = o[rule]["T_rot_gridsynth"] if rule else 0
            n1 = c.get("t", 0) + t_rot
            v, fac = F.price(n3, n1, "BEST", 1e-3, queries=o["queries"])
            if v is None:
                continue
            fc = F.cheapest("CCZ", 1e-3, 1.0)
            v3 = 4 * F.cheapest("T", 1e-3, 1.0)["qubitcycles"] if fac["Toffoli-class"] == "4 T" else fc["qubitcycles"]
            vt = v / max(n1, 1) if n3 == 0 else None
            # unit prices actually used
            name_t = fac["T"]
            vt = fc["qubitcycles"] / 2 if name_t.startswith("catalysed") else next(f["qubitcycles"] for f in F.factories(1e-3) if f["name"] == name_t)
            if fac["Toffoli-class"] == "4 T":
                v3 = 4 * vt
            n3_loc = lc.get("and", 0) + lc.get("tof", 0) + lc.get("cswap", 0)
            loc_angles = r["local"]["angles"]
            all_angles = r["prep"]["angles"] + loc_angles
            t_rot_loc = t_rot * len(loc_angles) / len(all_angles) if all_angles else 0
            v_loc = n3_loc * v3 + (lc.get("t", 0) + t_rot_loc) * vt
            v_rot = t_rot * vt
            share_loc = v_loc / v
            out["synthillation"].append({"lattice": r["lattice"], "method": r["method"], "L": r["L"], "eps": eps,
                                         "local_share": share_loc, "synthillation_bound": 1 / (1 - share_loc),
                                         "rotation_share": v_rot / v})

    # 4. nested reflections for 3-class patches
    for name, ev in (("J1-J2 triangle r=0.25 (lcvcbe tri_0.25)", None), ("isotropic tetrahedron", tetra_classes())):
        if ev is None:
            cl = IO.clusters()["tri_0.25"]
            ev = sorted({round(cl["hw"] * mu + cl["shift"], 9) for mu, _ in cl["classes"]})
        b = [(ev[k] - ev[k - 1]) / 2 for k in range(1, len(ev))]
        hw = (ev[-1] - ev[0]) / 2
        out["nested"].append({"patch": name, "class_eigenvalues": ev, "weights_b": b, "sum_b": sum(b), "half_width": hw,
                              "eigen_rotations_spectral": 2 * sum(1 for e in ev if abs(abs((e - (ev[0] + ev[-1]) / 2) / hw) - 1) > 1e-9),
                              "prep_rotations_nested": 2 * (len(b) - 1)})

    (OUT / "phase3_snap.json").write_text(json.dumps(out, indent=1))
    L = ["# Phase 3: patch designs and magic-state supply", "",
         "## 1. Exactness obstruction", "",
         f"Optimal block of the isotropic triangle (lambda = 3): distinct entries {out['obstruction']['distinct_entries']}. "
         "Entries 1/3, 2/3 are outside Z[1/sqrt2, i], so no exact Clifford+T circuit encodes it (Giles-Selinger).", "",
         "## 2. Smallest ring-compatible normalization (necessary condition) and its cost", "",
         "| lattice | L | half-width | ring lambda | alpha penalty | lambda_eff penalty (re-centred) | (same center) |",
         "|---|---|---|---|---|---|---|"]
    for x in out["ring"]:
        L.append(f"| {x['lattice']} | {x['L']} | {x['half_width']:.3f} | {x['ring_lambda']:.3f} | {x['alpha_penalty']:.3f} | "
                 f"{x['lambda_eff_penalty_recentred']:.3f} | {x['lambda_eff_penalty_same_center']:.3f} |")
    L += ["", "## 3. Synthillation upper bound and rotation share (BEST pricing, p = 1e-3, rotations by rule C/B)", "",
          "Removing every rotation pays for a lambda_eff penalty P only if rotation share > 1 - 1/P.", "",
          "| lattice | method | L | eps | local-block share | synthillation bound | rotation share |", "|---|---|---|---|---|---|---|"]
    for x in out["synthillation"]:
        L.append(f"| {x['lattice']} | {x['method']} | {x['L']} | {x['eps']:g} | {x['local_share']:.3f} | {x['synthillation_bound']:.3f} | "
                 f"{x['rotation_share']:.3f} |")
    L += ["", "## 4. Nested-reflection decomposition", "", "| patch | class eigenvalues | b_k | sum b = half-width | eigen rotations (spectral) | PREP rotations (nested) |",
          "|---|---|---|---|---|---|"]
    for x in out["nested"]:
        L.append(f"| {x['patch']} | {x['class_eigenvalues']} | {[round(v, 4) for v in x['weights_b']]} | {x['sum_b']:.4f} = {x['half_width']:.4f} | "
                 f"{x['eigen_rotations_spectral']} | {x['prep_rotations_nested']} |")
    (OUT / "phase3_snap.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

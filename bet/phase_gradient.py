"""Phase 6: rotations as phase-gradient additions, re-priced with CCZ states.

lcvcbe (LAB_NOTES entry 4, coarse_angles.py) rejected phase-gradient rotations at 4(b - 2) T per rotation against
gridsynth: both are ~50-60 T at a block error of 1e-4. That comparison prices a Toffoli at 4 T. Under the Phase 1
pricing a Toffoli-class gate is one CCZ state, 2.6-3.3x cheaper at matched fidelity, so the comparison may flip.

Model per rotation at precision eps_r (the Phase 2 rule of the row: C for triangle walks, B otherwise):
  synth   gridsynth T count at eps_r (lcvcbe synth.tcount), priced as T states
  PG      angle rounded to the grid 2 pi k / 2^b and applied by adding k into a b-bit phase-gradient register:
          b - 2 Toffoli-class gates (Gidney, Quantum 2, 74 (2018); Sanders et al., PRX Quantum 1, 020312 (2020)).
          Rounding error |delta| <= pi / 2^b, operator error |delta| / 2, so b = ceil(log2(pi / eps_r) - 1).
          The rounded angle is a deterministic unitary, so the exact-inverse-pair argument of rule C holds for it.
  hybrid  the cheaper of synth and PG for the whole row
One-time cost: the phase-gradient state (b rotations) and 2b logical qubits (register plus adder workspace); not
included in the per-query price, listed in the table.

Prices: BEST (factories.price, Litinski Table 1, tier-honest: falls back to 4 T where the single CCZ factory is not
good enough) and matched-fidelity T-equivalents (Toffoli-class = CCZ/T ratio of phase1_factories.md: 1.53 at
p = 1e-3, 1.20 at p = 1e-4), which are free of tier jumps. T4 = every Toffoli-class gate as 4 T.

Go/no-go 6: at least 1.5x on the whole walk (factory space-time per run) beyond Phase 5's BEST.

Usage: python -m bet.phase_gradient  ->  results/phase6_phase_gradient.{json,md}
"""

import json
import math
from pathlib import Path

from bet import factories as F
from bet import lcvcbe_io as IO

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"


def bits(eps_r):
    return max(3, math.ceil(math.log2(math.pi / eps_r) - 1))


def ccz_over_t(p):
    fc = F.cheapest("CCZ", p, 1.0)
    return fc["qubitcycles"] / F.cheapest("T", p, fc["p_out"])["qubitcycles"]


def main():
    S = IO.synth()
    prof = {(r["lattice"], r["method"], r["L"]): r for r in json.loads((OUT / "phase0_profiles.json").read_text())}
    fid = json.loads((OUT / "phase2_fidelity.json").read_text())
    ratio = {p: ccz_over_t(p) for p in (1e-3, 1e-4)}
    res = []
    for o in fid:
        r = prof[(o["lattice"], o["method"], o["L"])]
        angles = r["prep"]["angles"] + r["local"]["angles"]
        if not angles:
            continue
        rule = "C" if "rule_C" in o else "B"
        er = o[f"rule_{rule}"]["eps_r"]
        c = r["total"]["counts"]
        n3 = c.get("and", 0) + c.get("tof", 0) + c.get("cswap", 0)
        t0 = c.get("t", 0)
        t_rot = sum(S.tcount(a, er) for a in angles)
        b = bits(er)
        pg3 = len(angles) * (b - 2)
        nq = o["queries"]
        out = {"lattice": o["lattice"], "method": o["method"], "L": o["L"], "N": o["N"], "eps": o["eps"], "rule": rule,
               "eps_r": er, "R": len(angles), "b": b, "n3": n3, "t": t0, "T_rot_synth": t_rot,
               "T_per_rot_synth": t_rot / len(angles), "toffoli_per_rot_PG": b - 2,
               "one_time": {"logical_qubits": 2 * b, "gradient_state_rotations": b}}
        out["T4"] = {"synth": 4 * n3 + t0 + t_rot, "PG": 4 * (n3 + pg3) + t0}
        for p in (1e-3, 1e-4):
            k = f"{p:g}"
            x = ratio[p]
            m = {"synth": x * n3 + t0 + t_rot, "PG": x * (n3 + pg3) + t0, "no rotations": x * n3 + t0}
            m["rotation share (synth)"] = t_rot / m["synth"]
            m["gain PG"] = m["synth"] / m["PG"]
            m["gain bound"] = m["synth"] / m["no rotations"]
            out[f"matched {k}"] = m
            vs, _ = F.price(n3, t0 + t_rot, "BEST", p, queries=nq)
            vp, fp = F.price(n3 + pg3, t0, "BEST", p, queries=nq)
            out[f"BEST {k}"] = {"synth_run": None if vs is None else vs * nq, "PG_run": None if vp is None else vp * nq,
                                "PG_factories": fp}
            if vs and vp:
                out[f"BEST {k}"]["gain hybrid"] = vs / min(vs, vp)
                out[f"BEST {k}"]["gain PG"] = vs / vp
        res.append(out)
    (OUT / "phase6_phase_gradient.json").write_text(json.dumps(res, indent=1))

    f2 = lambda v: "n/a" if v is None else f"{v:.2f}"  # noqa: E731
    L = ["# Phase 6: rotations as phase-gradient additions", "",
         f"Per rotation: gridsynth T at the row's eps_r (rule C for triangles, B for swap) against b - 2 Toffoli-class gates "
         f"on a b-bit phase-gradient register, b = ceil(log2(pi / eps_r) - 1). Matched-fidelity CCZ/T = "
         f"{ratio[1e-3]:.2f} (p = 1e-3), {ratio[1e-4]:.2f} (p = 1e-4). Gains are synth / PG on the whole query or run "
         "(> 1 means phase gradient wins); 'bound' is the gain if rotations were free.", "",
         "| lattice | method | L | N | eps | R/q | b | T/rot synth | Tof/rot PG | T4: synth / PG | rot. share 1e-3 | "
         "matched 1e-3: gain PG / bound | matched 1e-4: gain PG / bound | BEST run 1e-3: gain PG (hybrid) | "
         "BEST run 1e-4: gain PG (hybrid) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for o in res:
        m3, m4, b3, b4 = o["matched 0.001"], o["matched 0.0001"], o["BEST 0.001"], o["BEST 0.0001"]
        L.append(f"| {o['lattice']} | {o['method']} | {o['L']} | {o['N']} | {o['eps']:g} | {o['R']} | {o['b']} | "
                 f"{o['T_per_rot_synth']:.1f} | {o['toffoli_per_rot_PG']} | {o['T4']['synth'] / o['T4']['PG']:.2f} | "
                 f"{m3['rotation share (synth)']:.2f} | {m3['gain PG']:.2f} / {m3['gain bound']:.2f} | "
                 f"{m4['gain PG']:.2f} / {m4['gain bound']:.2f} | {f2(b3.get('gain PG'))} ({f2(b3.get('gain hybrid'))}) | "
                 f"{f2(b4.get('gain PG'))} ({f2(b4.get('gain hybrid'))}) |")
    tri = [o for o in res if o["method"] == "triangle"]
    g = [o["matched 0.001"]["gain PG"] for o in tri]
    gb = [o["matched 0.001"]["gain bound"] for o in tri]
    L += ["", f"Triangle walk, matched fidelity p = 1e-3: gain {min(g):.2f}-{max(g):.2f}x, bound {min(gb):.2f}-{max(gb):.2f}x. "
          f"One-time overhead: 2b = {2 * min(o['b'] for o in res)}-{2 * max(o['b'] for o in res)} logical qubits."]
    (OUT / "phase6_phase_gradient.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

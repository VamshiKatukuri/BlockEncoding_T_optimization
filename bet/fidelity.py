"""Phase 2a: rotation-precision rules and the magic-state fidelity budget per circuit.

Rotation precision eps_r per rotation (energy error budget fraction b of the target eps):
  A  flat 1e-10 (IonQ arXiv:2610.32946)
  B  first order (eigenphase perturbation, Bauer-Fike for the unitary walk): Lambda_eff R eps_r <= b eps
  C  second order, for SU(2)-invariant reflection patches (isotropic triangles) built as exact inverse pairs
     (W' ... W'^dag, PREP' ... PREP'^dag) and a singlet, translation-invariant ground state: the first-order bias
     vanishes (the SU(2) twirl of i[K, A] is i[twirl(K), A] = 0, A central in the permutation algebra; uniform
     PREPARE weight errors sum to zero). Measured (bet.bias, 16-18 sites): bias <= kappa n_p hw dA^2 with
     kappa <= 0.06; we use kappa = 0.1. dA <= 2 R_W eps_r (R_W rotations per W). So
     eps_r = sqrt(b eps / (kappa n_p hw)) / (2 R).
T per rotation: gridsynth (deterministic, unitary; lcvcbe synth.tcount). For reference, the expected count of
mixed-fallback repeat-until-success synthesis used by IonQ, 0.53 log2(1/eps) + 4.86 (a channel, so rule C's
exact-inverse-pair argument does not apply to it).

Fidelity budget: a circuit with n magic states needs p_state <= delta / n (delta = 0.1). Cultivation gives 2e-9
per T state at p = 1e-3. Splitting coherent phase estimation (K_full = N_q queries) into circuits of K queries
costs a repetition overhead K_full / K in the standard-quantum-limit model (energy precision from shots); this is
a conservative model, Heisenberg-limited short-depth schemes (e.g. QCELS) pay less.

Usage: python -m bet.fidelity  ->  results/phase2_fidelity.{json,md}
"""

import json
import math
from pathlib import Path

from bet import factories as F
from bet import lcvcbe_io as IO

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"
KAPPA, BUDGET, DELTA = 0.1, 0.1, 0.1
P_CULT = 2e-9


def rus_t(eps_r):
    return 0.53 * math.log2(1 / eps_r) + 4.86


def rules(row, eps):
    angles = row["prep"]["angles"] + row["local"]["angles"]
    R = len(angles)
    if R == 0:
        return {}
    out = {"A": 1e-10, "B": BUDGET * eps / (row["lambda_eff"] * R)}
    if row["method"] == "triangle":
        n_p = (4 * row["N"]) if row["lattice"].startswith("square") else (2 * row["N"] // 3)
        hw = row["alpha"] / n_p
        out["C"] = math.sqrt(BUDGET * eps / (KAPPA * n_p * hw)) / (2 * R)
    return out


def main():
    S = IO.synth()
    rows = json.loads((OUT / "phase0_profiles.json").read_text())
    res = []
    for r in rows:
        angles = r["prep"]["angles"] + r["local"]["angles"]
        c = r["total"]["counts"]
        n3 = c.get("and", 0) + c.get("tof", 0) + c.get("cswap", 0)
        for eps in (1e-2, 1e-3):
            o = {"lattice": r["lattice"], "method": r["method"], "L": r["L"], "N": r["N"], "eps": eps, "R_query": len(angles)}
            nq = r["queries"][f"{eps:g}"]["lambda_eff"]
            o["queries"] = nq
            for k, er in rules(r, eps).items():
                tg = sum(S.tcount(t, er) for t in angles)
                o[f"rule_{k}"] = {"eps_r": er, "T_rot_gridsynth": tg, "T_rot_rus_expected": len(angles) * rus_t(er)}
            best_rule = "C" if "rule_C" in o else ("B" if "rule_B" in o else None)
            n1 = c.get("t", 0) + (o[f"rule_{best_rule}"]["T_rot_gridsynth"] if best_rule else 0)
            sq = 4 * n3 + n1  # T states per query when Toffoli-class gates use 4 T
            o["T_states_query"] = sq
            o["p_needed_coherent"] = DELTA / (sq * nq)
            kmax = DELTA / (P_CULT * sq)
            o["cultivation_K_max"] = kmax
            o["split_overhead_SQL"] = max(1.0, nq / kmax)
            v_best, fac = F.price(n3, n1, "BEST", 1e-3, queries=nq)
            v_cult = sq * F.CULTIVATION["qubitcycles"]
            o["BEST_qubitcycles_query"] = v_best
            o["BEST_factories"] = fac
            o["CULT_qubitcycles_query"] = v_cult
            if v_best:
                o["cultivation_saving_per_query"] = v_best / v_cult
                o["net_gain_split_cultivation"] = v_best / v_cult / o["split_overhead_SQL"]
            res.append(o)
    (OUT / "phase2_fidelity.json").write_text(json.dumps(res, indent=1))

    L = ["# Phase 2a: rotation precision and magic-state fidelity budget", "",
         f"Rules: A flat 1e-10; B first order; C second order (triangle encodings, kappa = {KAPPA}). Budget {BUDGET} of eps for "
         f"rotations; failure budget delta = {DELTA} per circuit; cultivation {P_CULT:g} per T at p = 1e-3 (volume estimate "
         f"{F.CULTIVATION['qubitcycles']:.2g} qubitcycles).", "",
         "| lattice | method | L | eps | R/q | gridsynth T/q A | B | C | RUS T/q A | C | T states/q | queries | p needed (coherent) | "
         "cultivation K_max | split overhead | cult. saving/q | net |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    g = lambda o, k, f: f"{o[k][f]:.0f}" if k in o else "-"  # noqa: E731
    for o in res:
        net = o.get("net_gain_split_cultivation")
        L.append(f"| {o['lattice']} | {o['method']} | {o['L']} | {o['eps']:g} | {o['R_query']} | {g(o, 'rule_A', 'T_rot_gridsynth')} | "
                 f"{g(o, 'rule_B', 'T_rot_gridsynth')} | {g(o, 'rule_C', 'T_rot_gridsynth')} | {g(o, 'rule_A', 'T_rot_rus_expected')} | "
                 f"{g(o, 'rule_C', 'T_rot_rus_expected')} | {o['T_states_query']} | {o['queries']} | {o['p_needed_coherent']:.1e} | "
                 f"{o['cultivation_K_max']:.3g} | {o['split_overhead_SQL']:.3g} | "
                 f"{o.get('cultivation_saving_per_query', float('nan')):.3g} | {'-' if net is None else f'{net:.3g}'} |")
    (OUT / "phase2_fidelity.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

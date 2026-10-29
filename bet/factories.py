"""Phase 1: magic-state factory cost model and re-pricing of the Phase 0 demand profiles.

Cost unit: physical qubitcycles (physical qubits x surface-code cycles) per output state, as in Litinski,
Quantum 3, 205 (2019), Table 1 (arXiv:1905.06903v3). Output errors are per state.

Policies for one block-encoding query:
  T4      every Toffoli-class gate (AND, Toffoli, CSWAP) as 4 T states, T from the cheapest T factory that meets
          the per-state error target (the convention of IonQ arXiv:2610.32946 and Summer & Jamet arXiv:2610.03686)
  BEST    each Toffoli-class gate from the cheaper of a CCZ factory (Litinski 8-to-CCZ) and 4 T states; T states
          from the cheaper of a T factory and catalysed CCZ -> 2 T (Gidney & Fowler, Quantum 3, 135 (2019)).
          Litinski lists one CCZ factory per physical error rate, so BEST falls back to 4 T when that factory is
          overkill or insufficient for the target.
  CULT    magic state cultivation (Gidney, Shutty, Jones arXiv:2409.17595) for T; Toffoli-class = 4 T
          (temporary AND with measurement-based uncomputation). Only reaches 2e-9 per state at p = 1e-3.

Rotations: T count per rotation from gridsynth (lcvcbe synth.tcount) at precision eps_r, priced as T states.

Cultivation volume is not tabulated in the paper's text; it states the volume is "nearly four times" that of a
distance-15 lattice-surgery CNOT and an order of magnitude below prior factories at ~1e-9. We use
1.5e5 qubitcycles per T state (band 5e4 .. 3e5), labelled as an estimate.

Usage: python -m bet.factories  ->  results/phase1_factories.{json,md}
"""

import json
import math
from pathlib import Path

from bet import lcvcbe_io as IO

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"

# (name, kind, p_phys, p_out, qubits, cycles, qubitcycles per output state)  -- Litinski Table 1
LITINSKI = [
    ("(15-to-1)_{7,3,3}", "T", 1e-4, 4.4e-8, 810, 18.1, 14_600),
    ("(15-to-1)_{9,3,3}", "T", 1e-4, 9.3e-10, 1_150, 18.1, 20_700),
    ("(15-to-1)_{11,5,5}", "T", 1e-4, 1.9e-11, 2_070, 30.0, 62_000),
    ("(15-to-1)^4_{9,3,3} x (20-to-4)_{15,7,9}", "T", 1e-4, 2.4e-15, 16_400, 90.3, 371_000),
    ("(15-to-1)^4_{7,3,3} x (8-to-CCZ)_{15,7,9}", "CCZ", 1e-4, 7.2e-14, 12_400, 36.1, 447_000),
    ("(15-to-1)_{17,7,7}", "T", 1e-3, 4.5e-8, 4_620, 42.6, 197_000),
    ("(15-to-1)^6_{13,5,5} x (20-to-4)_{23,11,13}", "T", 1e-3, 1.4e-10, 43_300, 130, 1_410_000),
    ("(15-to-1)^4_{13,5,5} x (20-to-4)_{27,13,15}", "T", 1e-3, 2.6e-11, 46_800, 157, 1_840_000),
    ("(15-to-1)^6_{11,5,5} x (15-to-1)_{25,11,11}", "T", 1e-3, 2.7e-12, 30_700, 82.5, 2_540_000),
    ("(15-to-1)^6_{13,5,5} x (15-to-1)_{29,11,13}", "T", 1e-3, 3.3e-14, 39_100, 97.5, 3_810_000),
    ("(15-to-1)^6_{13,7,7} x (8-to-CCZ)_{25,15,15}", "CCZ", 1e-3, 5.2e-11, 47_000, 60.0, 2_820_000),
]
CULTIVATION = {"name": "cultivation d1=5 (estimate)", "kind": "T", "p_phys": 1e-3, "p_out": 2e-9,
               "qubitcycles": 1.5e5, "band": (5e4, 3e5)}


def factories(p_phys):
    return [dict(zip(("name", "kind", "p_phys", "p_out", "qubits", "cycles", "qubitcycles"), f))
            for f in LITINSKI if f[2] == p_phys]


def cheapest(kind, p_phys, p_target):
    ok = [f for f in factories(p_phys) if f["kind"] == kind and f["p_out"] <= p_target]
    return min(ok, key=lambda f: f["qubitcycles"]) if ok else None


def rotation_tcounts(angles, eps_r):
    S = IO.synth()
    return sum(S.tcount(t, eps_r) for t in angles)


def demand(row, eps, budget=0.1, queries="lambda_eff"):
    """States per query and per run for one Phase 0 row: Toffoli-class n3, T n1 (incl. synthesized rotations)."""
    c = row["total"]["counts"]
    angles = row["prep"]["angles"] + row["local"]["angles"]
    nq = row["queries"][f"{eps:g}"][queries]
    X = row["lambda_eff"] if queries == "lambda_eff" else row["alpha"]
    eps_r = budget * eps / (X * len(angles)) if angles else 1.0
    n3 = c.get("and", 0) + c.get("tof", 0) + c.get("cswap", 0)
    n1 = c.get("t", 0) + rotation_tcounts(angles, eps_r)
    return {"n3_query": n3, "n1_query": n1, "queries": nq, "eps_r": eps_r}


def price(n3, n1, policy, p_phys, delta_fail=0.1, queries=1):
    """Qualified price per query (qubitcycles) and the factories used, for a run of `queries` queries."""
    states_T4 = (4 * n3 + n1) * queries
    if policy == "T4":
        f = cheapest("T", p_phys, delta_fail / states_T4)
        return (None, None) if f is None else ((4 * n3 + n1) * f["qubitcycles"], {"T": f["name"]})
    if policy == "BEST":
        for n_states in ((n3 + n1) * queries, (4 * n3 + n1) * queries):
            pt = delta_fail / n_states
            fc, ft = cheapest("CCZ", p_phys, pt), cheapest("T", p_phys, pt)
            if fc is None and ft is None:
                continue
            vt_opts = [(ft["qubitcycles"], ft["name"])] if ft else []
            if fc:
                vt_opts.append((fc["qubitcycles"] / 2, "catalysed CCZ->2T"))
            vt, tname = min(vt_opts)
            v3_opts = [(4 * vt, "4 T")] + ([(fc["qubitcycles"], fc["name"])] if fc else [])
            v3, cname = min(v3_opts)
            return n3 * v3 + n1 * vt, {"Toffoli-class": cname, "T": tname}
        return None, None
    if policy == "CULT":
        if p_phys != CULTIVATION["p_phys"] or CULTIVATION["p_out"] > delta_fail / states_T4:
            return None, None
        return (4 * n3 + n1) * CULTIVATION["qubitcycles"], {"T": CULTIVATION["name"]}
    raise ValueError(policy)


def burst(row, fc_cycles, factory_qubits_budget, fc_qubits):
    """Cycles per query when the CCZ supply is the bottleneck: n3 states from F parallel CCZ factories."""
    F = max(1, int(factory_qubits_budget // fc_qubits))
    n3 = row["total"]["counts"].get("cswap", 0) + row["total"]["counts"].get("and", 0) + row["total"]["counts"].get("tof", 0)
    widest = max((x["width"] for x in row["routing"].get("layers_per_direction", [])), default=None)
    return {"factories": F, "cycles_per_query_supply_bound": n3 * fc_cycles / F,
            "widest_layer": widest, "factories_for_one_layer_per_25_cycles": None if widest is None else math.ceil(widest * fc_cycles / 25)}


def main():
    rows = json.loads((OUT / "phase0_profiles.json").read_text())
    res = []
    for r in rows:
        for eps in (1e-2, 1e-3):
            d = demand(r, eps)
            for p in (1e-3, 1e-4):
                out = {"lattice": r["lattice"], "method": r["method"], "L": r["L"], "N": r["N"], "eps": eps, "p_phys": p, **d}
                for pol in ("T4", "BEST", "CULT"):
                    v, fs = price(d["n3_query"], d["n1_query"], pol, p, queries=d["queries"])
                    out[pol] = {"qubitcycles_query": v, "factories": fs,
                                "qubitcycles_run": None if v is None else v * d["queries"]}
                if out["T4"]["qubitcycles_query"] and out["BEST"]["qubitcycles_query"]:
                    out["T4/BEST"] = out["T4"]["qubitcycles_query"] / out["BEST"]["qubitcycles_query"]
                if out["T4"]["qubitcycles_query"] and out["CULT"]["qubitcycles_query"]:
                    out["T4/CULT"] = out["T4"]["qubitcycles_query"] / out["CULT"]["qubitcycles_query"]
                fc = cheapest("CCZ", p, 1.0)
                out["burst"] = burst(r, fc["cycles"], 1e6, fc["qubits"])
                res.append(out)
    (OUT / "phase1_factories.json").write_text(json.dumps(res, indent=1))

    L = ["# Phase 1: factory-priced cost per query", "",
         "Qubitcycles per query (physical qubits x code cycles), lambda_eff query count for the per-state error target "
         "(failure budget 0.1 per run). T4 = all Toffoli-class gates as 4 T; "
         "BEST = cheaper of CCZ factory and 4 T per Toffoli-class gate; CULT = cultivation (p = 1e-3 only, estimated volume). Rotations: gridsynth at eps_r = 0.1 eps / (lambda_eff R).", "",
         "| lattice | method | L | N | eps | p | Toffoli-class/q | T/q | queries | T4 | BEST | CULT | T4/BEST | T4/CULT | Toffoli-class source |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    fmt = lambda v: "n/a" if v is None else f"{v:.3g}"  # noqa: E731
    for o in res:
        L.append(f"| {o['lattice']} | {o['method']} | {o['L']} | {o['N']} | {o['eps']:g} | {o['p_phys']:g} | {o['n3_query']} | "
                 f"{o['n1_query']} | {o['queries']} | {fmt(o['T4']['qubitcycles_query'])} | {fmt(o['BEST']['qubitcycles_query'])} | "
                 f"{fmt(o['CULT']['qubitcycles_query'])} | {fmt(o.get('T4/BEST'))} | {fmt(o.get('T4/CULT'))} | "
                 f"{(o['BEST']['factories'] or {}).get('Toffoli-class', 'n/a')} |")
    L += ["", "## Unit prices at matched fidelity (Litinski Table 1)", "",
          "| p | CCZ factory (p_out, qubitcycles) | nearest T factory with p_out <= CCZ's (p_out, qubitcycles) | CCZ / T | "
          "price of a Toffoli-class gate under T4 / under CCZ |", "|---|---|---|---|---|"]
    for p in (1e-3, 1e-4):
        fc = cheapest("CCZ", p, 1.0)
        ft = cheapest("T", p, fc["p_out"])
        L.append(f"| {p:g} | {fc['p_out']:.2g}, {fc['qubitcycles']:.3g} | {ft['p_out']:.2g}, {ft['qubitcycles']:.3g} | "
                 f"{fc['qubitcycles'] / ft['qubitcycles']:.2f} | {4 * ft['qubitcycles'] / fc['qubitcycles']:.2f}x |")
    L += ["", "## Matched-fidelity ratio per encoding: (4 n3 + n1) v_T / (n3 v_CCZ + n1 v_T), eps = 1e-3", "",
          "Free of factory-tier jumps: both policies use the CCZ factory and the T factory of the row above.", "",
          "| lattice | method | L | Toffoli-class share of states (T4 count) | p = 1e-3 | p = 1e-4 |", "|---|---|---|---|---|---|"]
    for o in res:
        if o["eps"] == 1e-3 and o["p_phys"] == 1e-3:
            n3, n1 = o["n3_query"], o["n1_query"]
            rat = []
            for p in (1e-3, 1e-4):
                fc = cheapest("CCZ", p, 1.0)
                vt = cheapest("T", p, fc["p_out"])["qubitcycles"]
                rat.append((4 * n3 + n1) * vt / (n3 * fc["qubitcycles"] + n1 * vt))
            L.append(f"| {o['lattice']} | {o['method']} | {o['L']} | {4 * n3 / (4 * n3 + n1):.3f} | {rat[0]:.2f} | {rat[1]:.2f} |")
    L += ["", "## CCZ supply versus routing-layer width (1e6 physical qubits of CCZ factories)", "",
          "| lattice | method | L | p | factories | widest layer (CCZ) | factories to finish it in 25 cycles | supply-bound cycles/query |",
          "|---|---|---|---|---|---|---|---|"]
    for o in res:
        if o["eps"] == 1e-3 and o["burst"]["widest_layer"]:
            b = o["burst"]
            L.append(f"| {o['lattice']} | {o['method']} | {o['L']} | {o['p_phys']:g} | {b['factories']} | {b['widest_layer']} | "
                     f"{b['factories_for_one_layer_per_25_cycles']} | {b['cycles_per_query_supply_bound']:.3g} |")
    (OUT / "phase1_factories.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

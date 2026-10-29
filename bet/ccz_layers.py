"""Phase 4: routing layers as batches of shared-control CCZ gates.

11. SWAP = V^dag CZ V with V = (H x I) CNOT (Bell basis -> computational basis), so a layer of m CSWAPs that share
    the control c equals (prod V_i)^dag [prod_i CCZ(c, a_i, b_i)] (prod V_i). Checked by dense unitaries, m <= 3.
12. T count of prod_i CCZ(c, a_i, b_i):
    - phase-polynomial (unitary, no ancilla) form: each CCZ is 7 pi/4 parity rotations; the m copies of the
      rotation on c merge, giving 6m + (m mod 2). PyZX full_reduce is run as an independent optimizer.
    - measurement-based form (temporary AND, Jones 2013 / Gidney 2018): 4 T per CCZ, i.e. 4m. The function
      c * IP(a, b) has multiplicative complexity m (inner product), so m AND gates are needed in any
      AND-based circuit; with 4 T per AND the batch stays at 4m T, i.e. one CCZ state per CSWAP.
13. Parallel supply: CCZ states per layer and the factory count needed to finish the widest layer within one
    logical step, against Litinski's CCZ factory. Magic-state-fountain constants (arXiv:2501.19375,
    arXiv:2601.06736) are asymptotic only, so no factory-equivalent is assigned to them.

Usage: python -m bet.ccz_layers  ->  results/phase4_ccz_layers.{json,md}
"""

import json
import math
from pathlib import Path

import numpy as np

from bet import factories as F

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"


def kron(*ops):
    U = np.array([[1.0 + 0j]])
    for o in ops:
        U = np.kron(U, o)
    return U


def perm_unitary(n, f):
    """Unitary of a classical reversible map f on n-bit basis states (bit 0 most significant)."""
    dim = 2**n
    U = np.zeros((dim, dim), dtype=complex)
    for b in range(dim):
        bits = [(b >> (n - 1 - i)) & 1 for i in range(n)]
        nb = f(bits)
        U[sum(v << (n - 1 - i) for i, v in enumerate(nb)), b] = 1
    return U


def check_bell_rewrite(m):
    n = 1 + 2 * m  # control, then pairs (a_i, b_i)

    def cswaps(bits):
        bits = list(bits)
        if bits[0]:
            for i in range(m):
                a, b = 1 + 2 * i, 2 + 2 * i
                bits[a], bits[b] = bits[b], bits[a]
        return bits

    target = perm_unitary(n, cswaps)
    H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
    I2 = np.eye(2)
    cnot = perm_unitary(2, lambda x: [x[0], x[0] ^ x[1]])
    V1 = kron(H, I2) @ cnot
    V = kron(I2, *([V1] * m))
    diag = np.ones(2**n, dtype=complex)
    for b in range(2**n):
        bits = [(b >> (n - 1 - i)) & 1 for i in range(n)]
        par = sum(bits[0] * bits[1 + 2 * i] * bits[2 + 2 * i] for i in range(m)) % 2
        diag[b] = -1 if par else 1
    got = V.conj().T @ np.diag(diag) @ V
    return float(np.linalg.norm(got - target))


def pyzx_tcount(m):
    import pyzx as zx

    n = 1 + 2 * m
    c = zx.Circuit(n)
    for i in range(m):
        c.add_gate("CCZ", 0, 1 + 2 * i, 2 + 2 * i)
    base = zx.Circuit.from_graph(c.to_graph()).to_basic_gates()
    t0 = zx.tcount(base)
    g = c.to_graph()
    zx.full_reduce(g)
    t1 = zx.tcount(g)
    return t0, t1


def main():
    res = {"bell_rewrite_error": {m: check_bell_rewrite(m) for m in (1, 2, 3)}, "tcount": [], "supply": []}
    for m in range(1, 9):
        t_naive, t_pyzx = pyzx_tcount(m)
        res["tcount"].append({"m": m, "7m": 7 * m, "merged_phase_polynomial": 6 * m + (m % 2), "pyzx_input": t_naive,
                              "pyzx_full_reduce": t_pyzx, "measurement_based": 4 * m})
    rows = json.loads((OUT / "phase0_profiles.json").read_text())
    for p in (1e-3, 1e-4):
        fc = F.cheapest("CCZ", p, 1.0)
        for r in rows:
            lay = r["routing"].get("layers_per_direction")
            if not lay or r["L"] not in (8, 16, 32, 64, 4, 18):
                continue
            w = max(x["width"] for x in lay)
            res["supply"].append({"lattice": r["lattice"], "method": r["method"], "L": r["L"], "p": p, "widest_layer": w,
                                  "layers_per_query": 2 * len(lay),
                                  "factories_one_layer_per_d_cycles": {d: math.ceil(w * fc["cycles"] / d) for d in (15, 25, 35)},
                                  "factory_qubits_for_25": math.ceil(w * fc["cycles"] / 25) * fc["qubits"]})
    (OUT / "phase4_ccz_layers.json").write_text(json.dumps(res, indent=1))
    L = ["# Phase 4: routing layers as shared-control CCZ batches", "",
         "## 11. Bell-basis rewrite (dense check)", "",
         "|| (prod V_i)^dag prod_i CCZ(c, a_i, b_i) (prod V_i) - prod_i CSWAP(c; a_i, b_i) || = " +
         ", ".join(f"{v:.1e} (m = {m})" for m, v in res["bell_rewrite_error"].items()), "",
         "## 12. T count of m CCZ gates sharing one control", "",
         "| m | 7m | merged phase polynomial 6m + (m mod 2) | PyZX input | PyZX full_reduce | measurement-based 4m |", "|---|---|---|---|---|---|"]
    for x in res["tcount"]:
        L.append(f"| {x['m']} | {x['7m']} | {x['merged_phase_polynomial']} | {x['pyzx_input']} | {x['pyzx_full_reduce']} | {x['measurement_based']} |")
    L += ["", "## 13. CCZ supply needed to run one routing layer per logical step", "",
          "| lattice | method | L | p | widest layer (CCZ) | CSWAP layers / query | factories for one layer per 15 / 25 / 35 cycles | "
          "factory qubits (25 cycles) |", "|---|---|---|---|---|---|---|---|"]
    for x in res["supply"]:
        f = x["factories_one_layer_per_d_cycles"]
        L.append(f"| {x['lattice']} | {x['method']} | {x['L']} | {x['p']:g} | {x['widest_layer']} | {x['layers_per_query']} | "
                 f"{f[15]} / {f[25]} / {f[35]} | {x['factory_qubits_for_25']:.2e} |")
    (OUT / "phase4_ccz_layers.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

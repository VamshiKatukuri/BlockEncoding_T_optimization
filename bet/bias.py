"""Phase 2b: exact dyadic rotations keep the spectral patch an exact reflection; measure the energy bias.

The spectral patch block is B = W (sum_j R(mu_j) (x) Q_j) W^dag with R(mu) = Ry(2 arccos mu) Z. For any angles,
R(mu) is a Hermitian unitary, so B stays an exact reflection when the rotations of W or the eigenvalue rotations
are replaced by exact rotations through nearby dyadic angles (phase-gradient addition). The walk then
block-encodes a nearby H' = sum_p (hw A'_p + shift) exactly, with
    |E0(H') - E0(H)| <= ||H' - H|| <= sum_p hw ||A'_p - A_p||.

Here the continuous angles of the lcvcbe triangle blocks (clusters_v2) are rounded to b bits (multiples of
2 pi / 2^b), and the ground-state energy of H' is computed by Lanczos on small periodic lattices:
  square J1-J2, 4 x 4 (16 sites): 4 triangles per plaquette, block 'tri_0.5' (J2 = 0.5) or 'tri_0.25'
  kagome 2 x 2 and 2 x 3 cells (12, 18 sites): up and down triangles, block 'tri_0.5' with hw = 3

Usage: python -m bet.bias  ->  results/phase2_bias.{json,md}
"""

import json
import math
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as sla

from bet import lcvcbe_io as IO

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"
BITS = [6, 8, 10, 12, 14, 16]


def snapped_cluster(cl, bits):
    """Copy of a cluster with its continuous W angles rounded to multiples of 2 pi / 2^bits.

    Eigenvalue classes mu with |mu| < 1 are moved so that 2 arccos(mu) is on the same grid."""
    g = 2 * math.pi / 2**bits
    x = [round(t / g) * g if not s else t for t, s in zip(cl["x"], cl["snapped"])]
    cls = []
    for mu, size in cl["classes"]:
        if abs(abs(mu) - 1) > 1e-9 and abs(mu) > 1e-9:
            th = 2 * math.acos(max(-1.0, min(1.0, mu)))
            mu = math.cos(round(th / g) * g / 2)
        cls.append([mu, size])
    return {**cl, "x": x, "classes": cls}


def apply_local(psi, op, sites, n):
    """Apply a 2^k x 2^k operator to sites (qubit 0 = most significant) of an n-qubit state."""
    k = len(sites)
    t = psi.reshape([2] * n)
    t = np.moveaxis(t, sites, list(range(k)))
    sh = t.shape
    t = (op @ t.reshape(2**k, -1)).reshape(sh)
    return np.moveaxis(t, list(range(k)), sites).reshape(-1)


def ground_energy(terms, n, const):
    """Lowest eigenvalue of sum_t op_t on sites_t + const (matrix-free Lanczos)."""
    dim = 2**n

    def mv(v):
        v = v.astype(complex, copy=False)
        out = np.zeros(dim, dtype=complex)
        for op, s in terms:
            out += apply_local(v, op, s, n)
        return out

    A = sla.LinearOperator((dim, dim), matvec=mv, dtype=complex)
    w = sla.eigsh(A, k=1, which="SA", tol=1e-10, maxiter=5000)[0]
    return float(w[0].real) + const


def square_triangles(L):
    """4 triangles per plaquette: (a, m, b) with m the right-angle corner, a and b its two neighbours."""
    site = lambda x, y: (x % L) * L + (y % L)  # noqa: E731
    tris = []
    for x in range(L):
        for y in range(L):
            c = [site(x, y), site(x + 1, y), site(x + 1, y + 1), site(x, y + 1)]
            for m in range(4):
                tris.append((c[(m - 1) % 4], c[m], c[(m + 1) % 4]))
    return tris


def kagome_triangles(Lx, Ly):
    """Up triangles A(x,y), B(x,y), C(x,y); down triangles A(x,y), B(x-1,y), C(x,y-1)."""
    idx = lambda s, x, y: 3 * ((x % Lx) * Ly + (y % Ly)) + s  # noqa: E731
    tris = []
    for x in range(Lx):
        for y in range(Ly):
            tris.append((idx(0, x, y), idx(1, x, y), idx(2, x, y)))
            tris.append((idx(0, x, y), idx(1, x - 1, y), idx(2, x, y - 1)))
    return tris


def reference_triangle(J_nn, J_nnn):
    """Pauli-unit triangle J_nn (s_a.s_m + s_m.s_b) + J_nnn s_a.s_b on (a, m, b)."""
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1.0 + 0j, -1.0])
    I2 = np.eye(2)

    def ss(i, j):
        tot = np.zeros((8, 8), dtype=complex)
        for P in (X, Y, Z):
            ops = [I2, I2, I2]
            ops[i], ops[j] = P, P
            tot += np.kron(np.kron(ops[0], ops[1]), ops[2])
        return tot

    return J_nn * (ss(0, 1) + ss(1, 2)) + J_nnn * ss(0, 2)


def case(name, tris, n, cl, hw, shift, href):
    A = IO.cluster_matrix(cl)
    rec = float(np.linalg.norm(hw * A + shift * np.eye(8) - href, 2))
    e0 = ground_energy([(hw * A, t) for t in tris], n, shift * len(tris))
    rows = []
    for b in BITS:
        Ab = IO.cluster_matrix(snapped_cluster(cl, b))
        dA = float(np.linalg.norm(Ab - A, 2))
        herm = float(np.linalg.norm(Ab - Ab.conj().T, 2))
        refl = float(np.linalg.norm(Ab @ Ab - np.eye(8), 2)) if all(abs(abs(m) - 1) < 1e-9 for m, _ in cl["classes"]) else None
        eb = ground_energy([(hw * Ab, t) for t in tris], n, shift * len(tris))
        rows.append({"bits": b, "dA": dA, "hermiticity_err": herm, "reflection_err": refl,
                     "E0_snapped": eb, "bias": eb - e0, "bound": len(tris) * hw * dA,
                     "bias_per_site": (eb - e0) / n})
    return {"case": name, "n": n, "triangles": len(tris), "hw": hw, "shift": shift,
            "decomposition_err": rec, "E0": e0, "E0_per_site": e0 / n, "rows": rows}


GATES = {
    "H": np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2),
    "S": np.diag([1, 1j]),
    "T": np.diag([1, np.exp(1j * math.pi / 4)]),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]]),
    "Z": np.diag([1.0 + 0j, -1.0]),
    "I": np.eye(2, dtype=complex),
    "W": np.exp(1j * math.pi / 4) * np.eye(2),  # global phase symbol of pygridsynth
}


def gridsynth_ry(theta, eps):
    """Unitary of the gridsynth Clifford+T approximation of Ry(theta) = S H Rz(theta) H S^dag (up to phase)."""
    import mpmath
    from pygridsynth.gridsynth import gridsynth_gates

    g = str(gridsynth_gates(theta=mpmath.mpf(theta), epsilon=mpmath.mpf(eps), up_to_phase=True))
    U = np.eye(2, dtype=complex)
    for ch in g:  # gate string in circuit order
        U = GATES[ch] @ U
    return GATES["S"] @ GATES["H"] @ U @ GATES["H"] @ GATES["S"].conj().T


def gridsynth_cluster_matrix(cl, eps):
    """A' = W' D W'^dag with every continuous W rotation replaced by its gridsynth approximation (exact inverse
    pair: the same approximate W' and its exact adjoint)."""
    n = int(cl.get("n", 3))
    x, sn, L = cl["x"], cl["snapped"], cl["L"]

    def ry(t):
        c, s = math.cos(t / 2), math.sin(t / 2)
        return np.array([[c, -s], [s, c]], dtype=complex)

    def layer(ops):
        U = ops[0]
        for o in ops[1:]:
            U = np.kron(U, o)
        return U

    def cnot(q, n):
        dim = 2**n
        U = np.zeros((dim, dim), dtype=complex)
        for b in range(dim):
            bits = [(b >> (n - 1 - i)) & 1 for i in range(n)]
            if bits[q]:
                bits[q + 1] ^= 1
            U[sum(v << (n - 1 - i) for i, v in enumerate(bits)), b] = 1
        return U

    W, k = np.eye(2**n, dtype=complex), 0
    for _ in range(L):
        W = layer([ry(x[k + q]) if sn[k + q] else gridsynth_ry(x[k + q], eps) for q in range(n)]) @ W
        k += n
        for q in range(n - 1):
            W = cnot(q, n) @ W
    W = layer([ry(x[k + q]) if sn[k + q] else gridsynth_ry(x[k + q], eps) for q in range(n)]) @ W
    d = []
    for mu, size in cl["classes"]:
        d += [mu] * size
    return W @ np.diag(d) @ W.conj().T


def extra_tests(cls):
    """(1) gridsynth-approximated W (exact inverse pair) on a reflection patch; (2) PREPARE weight errors."""
    out = []
    cl = cls["tri_0.5"]
    A = IO.cluster_matrix(cl)
    tris = kagome_triangles(2, 3)
    hw = 4 * cl["hw"]
    e0 = ground_energy([(hw * A, t) for t in tris], 18, 0.0)
    for eps in (1e-2, 1e-3, 1e-4):
        Ag = gridsynth_cluster_matrix(cl, eps)
        dA = float(np.linalg.norm(Ag - A, 2))
        eb = ground_energy([(hw * Ag, t) for t in tris], 18, 0.0)
        out.append({"test": "kagome 2x3, gridsynth W (exact inverse pair)", "eps": eps, "dA": dA,
                    "reflection_err": float(np.linalg.norm(Ag @ Ag - np.eye(8), 2)), "bias": eb - e0,
                    "bound": len(tris) * hw * dA})
    rng = np.random.default_rng(7)
    for eta in (1e-2, 1e-3, 1e-4):
        dw = rng.normal(size=len(tris))
        dw -= dw.mean()
        dw *= eta / np.abs(dw).max()
        w = 1 + dw  # relative weights; sum preserved (PREPARE stays normalized)
        eb = ground_energy([(hw * wi * A, t) for wi, t in zip(w, tris)], 18, 0.0)
        out.append({"test": "kagome 2x3, PREPARE weights 1 + dw, sum dw = 0", "eps": eta, "dA": eta,
                    "reflection_err": None, "bias": eb - e0, "bound": hw * float(np.abs(dw).sum())})
    return out


def main():
    cls = IO.clusters()
    res = []
    for key, J2 in (("tri_0.5", 0.5), ("tri_0.25", 0.25)):
        cl = cls[key]
        res.append(case(f"square J1-J2 J2={J2} 4x4 ({key})", square_triangles(4), 16, cl, cl["hw"], cl["shift"],
                        reference_triangle(0.25, J2 / 2)))
    cl = cls["tri_0.5"]
    for Lx, Ly in ((2, 2), (2, 3)):
        res.append(case(f"kagome {Lx}x{Ly} cells (tri_0.5 x 4)", kagome_triangles(Lx, Ly), 3 * Lx * Ly, cl,
                        4 * cl["hw"], 0.0, reference_triangle(1.0, 1.0)))
    extra = extra_tests(cls)
    (OUT / "phase2_bias.json").write_text(json.dumps({"dyadic": res, "extra": extra}, indent=1))
    L = ["# Phase 2b: energy bias of exact dyadic rotations (Lanczos)", ""]
    L += ["## Other error sources on a reflection patch (kagome 2x3 cells, 18 sites)", "",
          "| test | eps / eta | ||A'-A|| or max dw | A'^2 = 1 err | E0' - E0 | first-order bound |", "|---|---|---|---|---|---|"]
    for x in extra:
        refl = "-" if x["reflection_err"] is None else f"{x['reflection_err']:.1e}"
        L.append(f"| {x['test']} | {x['eps']:g} | {x['dA']:.2e} | {refl} | {x['bias']:.2e} | {x['bound']:.2e} |")
    L.append("")
    for r in res:
        L += [f"## {r['case']}", "",
              f"N = {r['n']}, triangles = {r['triangles']}, decomposition check ||hw A + shift - h|| = {r['decomposition_err']:.1e}, "
              f"E0/N = {r['E0_per_site']:.6f}", "",
              "| bits | max angle error | ||A'-A|| | A' Hermitian err | A'^2 = 1 err | E0' - E0 | bound sum hw ||A'-A|| | bound / actual |",
              "|---|---|---|---|---|---|---|---|"]
        for x in r["rows"]:
            ratio = x["bound"] / abs(x["bias"]) if abs(x["bias"]) > 1e-12 else float("inf")
            refl = "-" if x["reflection_err"] is None else f"{x['reflection_err']:.1e}"
            L.append(f"| {x['bits']} | {math.pi / 2**x['bits']:.1e} | {x['dA']:.2e} | {x['hermiticity_err']:.1e} | {refl} | "
                     f"{x['bias']:.2e} | {x['bound']:.2e} | {ratio:.1f} |")
        L.append("")
    (OUT / "phase2_bias.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

"""Swap-tree routing (Summer & Jamet arXiv:2610.03686; IonQ arXiv:2610.32946): structure, counts and a
classical permutation check.

A tree over M = 2^m positions brings position a to position 0 with M - 1 controlled SWAPs: for bit j = m-1 .. 0,
controlled on bit j of a, swap position k with k + 2^j for k = 0 .. 2^j - 1. All 2^j CSWAPs of layer j share
the control bit j (one fan-out). Routing is applied forward and inverse in every query.
"""

import math


def tree_layers(M):
    """Layer widths of one tree, in application order (bit m-1 first). Requires M a power of two."""
    m = int(round(math.log2(M)))
    assert 2**m == M, "layer structure needs M = 2^m"
    return [2**j for j in range(m - 1, -1, -1)]


def route(perm_positions, a, M):
    """Apply the tree classically: perm_positions[k] = item at position k. Returns the new list."""
    p = list(perm_positions)
    m = int(round(math.log2(M)))
    for j in range(m - 1, -1, -1):
        if (a >> j) & 1:
            for k in range(2**j):
                p[k], p[k + 2**j] = p[k + 2**j], p[k]
    return p


def check_tree(M):
    """Every address a brings item a to position 0."""
    return all(route(range(M), a, M)[0] == a for a in range(M))


def kagome_check(L):
    """Kagome on L x L cells, sublattices A, B, C, one tree per sublattice over the L^2 cells.

    Up triangle at cell (x, y): A(x, y), B(x, y), C(x, y). Down triangle: A(x, y), B(x - 1, y), C(x, y - 1)
    (any fixed convention works). The address of tree B is ((x - u) mod L, y), of tree C (x, (y - u) mod L), with
    u = 1 for down triangles: a controlled modular decrement of the address, no fix-up swaps (the Heisenberg
    triangle is invariant under permutations of its sites, so slot order does not matter).
    Returns True if every triangle lands in slots (A0, B0, C0)."""
    M = L * L
    idx = lambda x, y: (x % L) * L + (y % L)  # noqa: E731
    for x in range(L):
        for y in range(L):
            for u in (0, 1):
                want = {("A", idx(x, y)), ("B", idx(x - u, y)), ("C", idx(x, y - u))}
                got = {
                    ("A", route(range(M), idx(x, y), M)[0]),
                    ("B", route(range(M), idx(x - u, y), M)[0]),
                    ("C", route(range(M), idx(x, y - u), M)[0]),
                }
                if want != got:
                    return False
    return True


def routing_profile(N, chi, fixup_cswaps=0, power_of_two=True):
    """Routing demand per query (both directions).

    chi trees over M = N / chi positions; 2 (M - 1) CSWAPs per tree and query. Layers of all trees run in
    parallel; each layer j holds chi * 2^j CSWAPs in chi shared-control groups of size 2^j."""
    M = N // chi
    out = {"N": N, "chi": chi, "M": M, "cswap_tree": 2 * chi * (M - 1), "cswap_fixup": fixup_cswaps}
    if power_of_two and M & (M - 1) == 0:
        lay = tree_layers(M)
        out["layers_per_direction"] = [{"width": chi * w, "group": w, "groups": chi} for w in lay]
        out["depth_cswap_layers"] = 2 * len(lay)
    return out

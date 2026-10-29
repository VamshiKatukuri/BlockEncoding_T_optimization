"""Read-only access to the verified circuits and helpers of the lcvcbe repository.

The location is taken from the environment variable LCVCBE_DIR (default: a sibling checkout ../lcvcbe).
Nothing is written to that repository.
"""

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LCVCBE = Path(os.environ.get("LCVCBE_DIR", ROOT.parent / "lcvcbe")).resolve()


def _path():
    if not LCVCBE.is_dir():
        raise FileNotFoundError(f"lcvcbe checkout not found at {LCVCBE}; set LCVCBE_DIR")
    if str(LCVCBE) not in sys.path:
        sys.path.append(str(LCVCBE))


@lru_cache(maxsize=None)
def clusters():
    """Verified spectral 3-site clusters (lcvcbe results/clusters_v2.json)."""
    return json.loads((LCVCBE / "results" / "clusters_v2.json").read_text())


@lru_cache(maxsize=None)
def fixups():
    """Fix-up CSWAPs per query, both directions (lcvcbe results/fixup.json)."""
    return json.loads((LCVCBE / "results" / "fixup.json").read_text())


def block_circuit(cl):
    """cirq circuit of a controlled spectral cluster block (lcvcbe spectral_circuit.block_circuit)."""
    _path()
    import spectral_circuit as SC

    return SC.block_circuit(int(cl.get("n", 3)), cl["x"], cl["snapped"], cl["L"], cl["classes"])


def verify_block(cl, A=None):
    """Re-run lcvcbe's unitary check of a cluster block; returns its error dict."""
    _path()
    import spectral_circuit as SC

    circ, c, a, s, f = block_circuit(cl)
    if A is None:
        A = cluster_matrix(cl)
    return SC.verify(circ, c, a, s, f, A)


def cluster_matrix(cl):
    """Normalized cluster operator A = sum_j mu_j P_j reconstructed from the circuit's W and classes."""
    import numpy as np

    _path()
    import spectral_circuit as SC

    n = int(cl.get("n", 3))
    import cirq

    s = [cirq.NamedQubit(f"s{i}") for i in range(n)]
    W = SC.w_circuit(s, cl["x"], cl["snapped"], cl["L"]).unitary(qubit_order=s)
    d = []
    for mu, size in cl["classes"]:
        d += [mu] * size
    # block_circuit applies W^-1, the class rotations, then W: as a matrix W D W^dag
    return W @ np.diag(d) @ W.conj().T


def synth():
    _path()
    import synth as S

    return S

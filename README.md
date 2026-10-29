# BlockEncoding_T_optimization

Co-design of magic-state preparation (distillation, cultivation, CCZ factories) with routed, low-frustration
block encodings of lattice Hamiltonians. The plan and go/no-go criteria are in `docs/PLAN.md`.

## Setup

    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt

Versions are pinned because Qualtran gate counts can change between releases.

## Layout

- `docs/PLAN.md`: phases 0-5, go/no-go criteria, risks, compute needs.
- `RESULTS.md`: outcome of every phase (non-HPC scope) with confidence tags and caveats.
- `bet/`: one module per phase, each writing `results/phase*.{json,md}`:
  - `phase0.py`: demand profiles per query (routing, PREPARE, local block, rotation angles, swap-tree layers)
  - `routing.py`: swap-tree structure and classical routing checks (square sublattices, kagome)
  - `factories.py`: Phase 1 factory model (Litinski Table 1, catalysed CCZ->2T, cultivation estimate)
  - `bias.py`: Phase 2b Lanczos check of exact-reflection bias (dyadic, gridsynth, PREPARE weights)
  - `fidelity.py`: Phase 2a rotation-precision rules and magic-state fidelity budget
  - `snap.py`: Phase 3 exactness obstruction, ring-compatible normalization, synthillation bound
  - `ccz_layers.py`: Phase 4 Bell rewrite, CCZ-batch T counts (PyZX), parallel CCZ supply
  - `compare.py`: Phase 5 end-to-end comparison with Trotter + Hamming-weight phasing
  - `lcvcbe_io.py`: read-only access to `lcvcbe` (set `LCVCBE_DIR`, default `../lcvcbe`)

## Running

    export LCVCBE_DIR=/path/to/lcvcbe
    for m in phase0 factories bias fidelity snap ccz_layers compare; do .venv/bin/python -m bet.$m; done

Order matters: later phases read `results/phase0_profiles.json` and `results/phase2_fidelity.json`.

## Related repositories

- `VamshiKatukuri/lcvcbe`: source of the verified routed spectral-cluster circuits (read only).
- `VamshiKatukuri/BeyondBlockEncoding`: surveys and directions beyond global block encodings.

# BlockEncoding_T_optimization

Co-design of magic-state preparation (distillation, cultivation, CCZ factories) with routed, low-frustration
block encodings of lattice Hamiltonians. The plan and go/no-go criteria are in `docs/PLAN.md`.

## Setup

    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt

Versions are pinned because Qualtran gate counts can change between releases.

## Layout

- `docs/PLAN.md`: phases 0-5, go/no-go criteria, risks, compute needs.

## Related repositories

- `VamshiKatukuri/lcvcbe`: source of the verified routed spectral-cluster circuits (read only).
- `VamshiKatukuri/BeyondBlockEncoding`: surveys and directions beyond global block encodings.

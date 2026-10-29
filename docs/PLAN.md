# Plan: magic-state-aware block encodings

Status 2026-10-29: phases 0-5 done at non-HPC scope, see `RESULTS.md`. Go/no-go: 1 pass, 2 fail at moderate N
(2b adds a second-order precision rule), 3 fail, 4 fail.

Question: can magic-state preparation (distillation, cultivation, CCZ factories, parallel "fountains") be
co-designed with routed, low-frustration block encodings so that the real cost (factory space-time and
wall-clock) falls well below what T-count accounting suggests?

Confidence tags: **[Certain]** = read in a paper or derived directly, **[Likely]** = strong inference,
**[Guessing]** = to be tested.

## Starting point

- **[Certain]** In routed encodings (Summer & Jamet arXiv:2610.03686; IonQ arXiv:2610.32946) the routing is
  most of each query. IonQ's 3x3 J1-J2 patch at 120x120 spends about 1.15e5 of 1.29e5 T per query on routing.
- **[Certain]** Routing is controlled-SWAPs, i.e. Toffoli/CCZ-class gates, not T gates. Both papers count each
  one as 4 T.
- **[Likely]** At extensive accuracy and large N, phase estimation only returns the trial state's energy. So
  the benchmarks here are intensive-precision energies on moderate lattices (10x10 to 30x30) and time
  evolution.

## Four ways to combine magic-state preparation with block encoding

1. **CCZ supply for routing.** Use CCZ factories and catalysed CCZ -> 2T. Price routing in factory
   space-time, not in T gates.
2. **Synthillation for the fixed local patch** (Campbell & Howard, PRL 118, 060501). Translation invariance
   means one fixed patch circuit runs in every query. Choose patch shift, normalisation and weights so its
   spectral angles are multiples of pi/2^k.
3. **Routing layers as batched shared-control CCZs.** Feed them from parallel sources, such as magic state
   fountains on QLDPC codes or a code with a transversal CCZ.
4. **Fidelity co-design.**
   - Cut the non-Clifford count per circuit, using statistical or single-ancilla phase estimation or
     light-cone-local encodings, until cultivation or STAR-type rotations suffice without distillation.
   - Use exact dyadic rotations (phase-gradient addition) so that the spectral patch stays an exact reflection.
     The walk then encodes a nearby H' with a known bias.

## Phases

### Phase 0: Task and cost currency (about 1 week, laptop)
1. **Fix the tasks:** intensive-precision energy on 10x10 to 30x30 lattices, and time evolution.
2. **Fix the encodings:**
   - routed SWAP bonds;
   - routed triangles;
   - IonQ's routed 3x3 patch, rebuilt from their appendices;
   - kagome triangles, which are reflections with no rotations.
3. **Build a per-query magic demand profile** for each encoding:
   - controlled-SWAPs per routing layer, with shared controls marked;
   - fixed-angle rotations, listed by angle;
   - variable rotations;
   - T gates;
   - depth and parallelism per layer.

   Use Qualtran counts on the verified circuits in `lcvcbe` (read-only).

**Deliverable:** `results/demand_profiles.json` plus a table.

### Phase 1: Factory cost model (about 2 weeks)
4. **Space-time per output state at a target error** for each source:
   - 15-to-1 T factories;
   - Gidney-Fowler CCZ factories, including catalysed CCZ -> 2T;
   - magic state cultivation;
   - QLDPC (bicycle architecture, fountains).
5. **Re-price the Phase 0 profiles,** including factory count against the burst demand of routing layers.

**Go/no-go 1:** if CCZ pricing changes the total by less than 1.3x versus 4 T per controlled-SWAP, option 1 is
minor.

### Phase 2: Fidelity budget and error structure (about 2-3 weeks)
6. **Required magic-state error per circuit,** p <~ delta / (non-Cliffords per circuit). Compute it for:
   - standard qubitized phase estimation;
   - single-ancilla or statistical phase estimation;
   - split or light-cone-local runs.
7. **Exact dyadic rotations.**
   - Prove the spectral patch stays an exact reflection.
   - Bound the energy bias by Lambda * max|delta a_j|.
   - Check it with exact diagonalisation on 16-24 sites.
   - Read off the angle bits needed per accuracy target.

**Go/no-go 2:** the no-distillation regime is viable if each circuit has at most about 1e5-1e6
non-Cliffords and the extra repetitions cost at most 2x in total.

### Phase 3: Patch designs that suit the factories (about 3-4 weeks, after Phase 1)
8. **Optimise the free parameters** (shift, lambda_p >= half-width, overlap weights) so spectral angles land
   on the pi/2^k grid. Trade the Lambda_eff penalty against the rotations saved.
9. **Synthillation of the fixed diagonal part D_ref.** Find its Clifford-hierarchy level, pick a protocol and
   compute the factory saving.
10. **Patches that cannot hit the grid:** QET remapping, or splitting into a small LCU.

**Go/no-go 3:** at most 10% Lambda_eff penalty, and at least 50% of patch rotation cost removed. This mainly
matters at moderate N, where the patch is a larger share of each query.

### Phase 4: Routing layers as batched CCZs (about 3-4 weeks, after Phase 1)
11. **Rewrite each controlled-SWAP layer** as shared-control CCZs by conjugating with Bell-basis Cliffords.
    Verify on small lattices.
12. **T-count optimisation** of the phase polynomial `c * sum_i a_i b_i` (TODD, Amy-Mosca, PyZX). Look for
    fewer than 4 T per controlled-SWAP.
13. **Map onto parallel sources:** fountain throughput against layer width, or the cost of switching into a
    transversal-CCZ code.

**Go/no-go 4:** fewer than 4 T per controlled-SWAP, or a wall-clock reduction about equal to the layer width.

### Phase 5: End-to-end comparison (about 2 weeks)
14. **Compare the combination of the options that passed** against:
    - T-count-priced routed encodings (IonQ-style);
    - Trotter with Hamming-weight phasing.

    Report factory space-time and wall-clock time for surface code with cultivation, and for bicycle QLDPC.

## Risks

- **[Guessing]** Fountain and transversal-CCZ overheads for real codes are not well established.
- **[Likely]** The extra repetitions of statistical phase estimation can cancel the saving from skipping
  distillation.
- **[Guessing]** Patches larger than about 6 spins may leave too little freedom to snap angles.

## Compute

- **Phases 0-1:** laptop scale.
- **Phase 2:** exact diagonalisation up to 24 sites is local; 36 sites needs a cluster.
- **Phase 4:** T-count optimisers scale steeply; large batches need a cluster.

## Related repositories

- `VamshiKatukuri/lcvcbe`: verified routed spectral-cluster circuits and cost pipeline. Read only from here.
- `VamshiKatukuri/BeyondBlockEncoding`: literature surveys and directions that drop the global block encoding.

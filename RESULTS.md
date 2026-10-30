# Results: magic-state-aware block encodings (phases 0-6, non-HPC scope)

Written 2026-10-29. Plan and go/no-go criteria: `docs/PLAN.md`. Every number below comes from `results/` and can be
re-run with `python -m bet.<module>` (about 15 minutes in total on 4 cores, most of it the Lanczos runs in
`bet.bias`).

Confidence tags: **[Certain]** = computed or checked here, **[Likely]** = follows from the computation under the
stated model, **[Guessing]** = depends on an unverified input.

## Scope

- **Encodings:**
  - routed SWAP-LCU and routed spectral-triangle walks;
  - square J1-J2 at J2/J1 = 0.5 (L = 8-64);
  - kagome Heisenberg (L = 4-18 cells, N up to 972).
- **Task:** ground-state energy to an absolute error eps in {1e-2, 1e-3} (Pauli units), with lambda_eff query
  counts.
- **Inputs from `lcvcbe`, read only:** the verified triangle-block circuits, fix-up counts, the gridsynth helper
  and Trotter error constants. Square-lattice T per query reproduces `lcvcbe` `routed_cost2` exactly. **[Certain]**
- **Factories:** Litinski's Table 1 (Quantum 3, 205), including 8-to-CCZ; Gidney-Fowler catalysed CCZ -> 2T; and
  magic state cultivation at p = 1e-3. The cultivation volume is estimated, see caveats.

## Outcome by phase

| Phase | Question | Result | Go/no-go |
|---|---|---|---|
| 0 | Where does the non-Clifford demand sit? | **[Certain]** Routing CSWAPs are 74-98% of non-Clifford gates per query and rise with N. Swap-tree layer j holds 2^j CSWAPs sharing one address bit. | done |
| 1 | Does pricing Toffoli-class gates as CCZ states matter? | **[Certain]** At matched fidelity a CCZ state costs 1.53 (p = 1e-3) or 1.20 (p = 1e-4) T states. A Toffoli-class gate is therefore 2.6x or 3.3x cheaper than "4 T". Per query this gives 1.5-2.6x (p = 1e-3) and 1.6-3.3x (p = 1e-4). | **pass** (>= 1.3x) |
| 2a | Can cultivation replace distillation by splitting phase estimation? | **[Likely]** Coherent runs need 1e-10 to 1e-13 per state from L = 16 (eps = 1e-3); cultivation gives 2e-9. Splitting costs 10-1e4x more repetitions (shot-noise model) against a 5-17x per-state saving. The net gain is above 1 only for the smallest lattices or eps = 1e-2. | **fail** at moderate N |
| 2b | Does exactness relax the rotation precision? | **[Certain]** Exact dyadic angles, or any synthesis applied as an exact inverse pair, keep the patch an exact reflection (to 1e-14). For isotropic patches with a singlet, translation-invariant ground state the first-order energy bias vanishes. Lanczos on 12-18 sites shows a quadratic bias for basis-change, gridsynth and PREPARE-weight errors. **[Likely]** The resulting second-order rule cuts rotation T per query by 1.7-1.9x against the first-order rule and 2.2-2.6x against a flat 1e-10. | pass (new result) |
| 3 | Do patches chosen to suit the factories pay? | **[Certain]** The optimal triangle reflection has entries 1/3, 2/3, outside Z[1/sqrt2, i], so it has no exact Clifford+T circuit. The smallest admissible normalisation (4/3 larger) costs 1.16-1.18x in lambda_eff even after re-centring. **[Likely]** That pays only where rotations exceed ~15% of the query cost (L <= 10). Synthillation of the local block is bounded by <= 1.2x from L = 16 up. Nested reflections reach optimal normalisation without eigenvalue rotations but move the same number of rotations into PREPARE. | **fail** (needs <= 10% penalty) |
| 4 | Can routing layers be cheaper than one CCZ per CSWAP? | **[Certain]** The Bell rewrite of CSWAP layers into shared-control CCZ batches checks to 1e-14. PyZX finds only 6m + (m mod 2) T for m CCZs, never below the measurement-based 4m. **[Likely]** One CCZ per CSWAP is the floor, because c * IP(a, b) needs m AND gates. Running one layer per logical step needs 3.6e6-2.3e8 factory qubits, so routing is magic-supply bound. | **fail** |
| 5 | End to end | **[Certain]** Combining CCZ pricing (Phase 1) with the second-order rule (Phase 2b) gives 2.1-2.7x less factory space-time than the T-count convention for the routed triangle walk. Where the run outgrows the single published CCZ factory at p = 1e-3 (5.2e-11 per state), the gain falls to 1.03-1.5x. The triangle walk beats Trotter with Hamming-weight phasing by 1.2-100x here. On kagome at eps = 1e-2 the gap closes with N (0.82 at N = 972). | - |
| 6 | Do rotations pay as phase-gradient additions once Toffolis are CCZ-priced? | **[Certain]** Under 4 T per Toffoli, phase gradient ties synthesis (0.94-1.01x), as `lcvcbe` found. At matched-fidelity CCZ pricing a rotation costs b - 2 = 9-14 CCZ states (14-21 T-equivalents) against 28-46 T for gridsynth under rule C. **[Certain]** On the whole triangle walk this gives 1.01-1.34x at p = 1e-3 (1.01-1.51x at p = 1e-4). The upper bound, if rotations were free, is 1.01-1.85x. Rotations are 5-46% of a query's cost, and their share falls with N as routing grows. Tier-honest BEST pricing gives at most 1.2x, apart from factory-tier jumps such as kagome L = 16 (2.6x). There, removing rotation T states keeps the run inside the single CCZ factory's capacity. **[Certain]** QROM angle lookup has nothing to act on: the 4-8 rotations per query are fixed angles applied in every query, not selected by an index. | **fail** (needs >= 1.5x) |

## What this means

1. **[Likely]** For routed block encodings the useful pairing with magic-state preparation is supply-side:
   - price and produce Toffoli-class states (CCZ) rather than T states;
   - build patches so that rotation errors are second order.

   Together these give about 2-2.7x in factory space-time, without changing the algorithm.
2. **[Likely]** The three demand-side routes fail at these sizes:
   - skipping distillation;
   - patches matched to the factories, or synthillation of the local block;
   - T-count optimisation of the CCZ batches.

   The reason is that routing dominates, and routing is already at one CCZ per CSWAP.
3. **[Likely]** The log-depth parallelism of swap-tree routing cannot be used on factory-limited hardware. The
   widest layer needs hundreds to thousands of CCZ factories, so wall-clock time is set by run volume divided by
   factory qubits.
4. **[Guessing]** The biggest remaining lever on the factory side is a CCZ factory menu that is denser in
   fidelity. Litinski lists one CCZ protocol per physical error rate, and the CCZ gain vanishes where that one
   protocol is not good enough (square L >= 30 at eps = 1e-3).

5. **[Likely]** Rotations are not the lever. Pricing Toffolis as CCZ states turns the phase-gradient trade-off in
   favour of additions, but rotations are too small a share of a routed query for this to matter beyond L ~ 10
   (Phase 6). Only the routing cost moves the total.

## Caveats

- **Factory menu.** Litinski's single CCZ factory per p makes the BEST policy jump between tiers; the
  matched-fidelity ratios in `phase1_factories.md` remove this. Kagome L = 4 at eps = 1e-2 (10.9x) is such a jump,
  not a real gain.
- **Cultivation volume.** 1.5e5 qubitcycles per T state is estimated from the paper's statements (about 4x a d =
  15 lattice-surgery CNOT, and an order of magnitude below prior factories), not read from its figure. Band:
  5e4-3e5. **[Guessing]**
- **Splitting model.** The overhead of splitting phase estimation uses a shot-noise model. Heisenberg-limited
  short-depth schemes (e.g. QCELS) pay less, so Phase 2a is conservative.
- **Second-order rule.** kappa = 0.1 is about 2x above the largest value measured on 16-18 sites. The cancellation
  needs a non-degenerate singlet ground state and isotropic patches. It does not apply to the swap baseline's
  NN/NNN label rotation, or to non-isotropic 3-class patches such as `tri_0.25`.
- **Trotter.** Error constants come from <= 24 sites (`lcvcbe`), and its controlled steps are priced as
  uncontrolled, which favours Trotter.
- **Ground-state energies.** These are literature values: J1-J2 0.5 from DMRG (Gong et al. 2014), kagome from
  Depenbrock et al. 2012.
- **Excluded:** state preparation and overlap.

## What HPC would add

- Lanczos on 24-36 sites, to confirm that the second-order bias constant converges and to cover the pyrochlore.
- T-count optimisation (TODD-class) of larger CCZ batches and full routed circuits.
- A CCZ factory model denser in fidelity, plus cultivation at p = 5e-4, from simulation rather than published
  points.

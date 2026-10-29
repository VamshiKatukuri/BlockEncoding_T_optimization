# What the magic-state study adds to the existing manuscripts

Written 2026-10-29 after phases 0-5 (`RESULTS.md`), and checked against the `lcvcbe` notes, which already contain
physical estimates. The manuscripts concerned are the triangle paper (square J1-J2) and the clique paper (3D
frustrated magnets).

Confidence tags: **[Certain]** computed or checked, **[Likely]** follows under the stated model, **[Guessing]**
not verified.

## Bottom line

**[Certain]** About half of what this study found is already in the `lcvcbe` physical estimates and should not be
presented as new:
- CCZ -> 2T factories;
- "the walk is Toffoli-heavy, Trotter is T-heavy";
- runtimes bound by magic-state production.

**[Likely]** Two additions are worth a section or appendix: a second-order rotation-precision result, and a
quantified rebuttal of IonQ's routing-parallelism outlook. Three more fit as a remark or footnote. None of them
needs HPC.

## Already covered: do not re-add

| Finding here | Where it already is |
|---|---|
| A CCZ factory makes a Toffoli for about 2 T, not 4 T, so physical ratios favour the walk | `trotter/NOTES.md` step 3; `NATCOMM_PLAN.md` step 3; `gap/NOTES.md` |
| Runtime is bound by magic-state production | `trotter/NOTES.md` step 3, stated as an assumption (single factory, no reaction limit) |
| Rotations can be synthesised up to a global phase because each occurs with its inverse | `synth.py` header, which points to the manuscript's error appendix |
| Exact phase-gradient angles cost as much as synthesis at the needed accuracy | README follow-up 2 (coarse angles) |
| Kagome triangle routing needs no fix-up swaps | `results/fixup.json` ("kagome cluster": 0) |

## Additions, ranked

### A1. Second-order synthesis error for SU(2)-central blocks (appendix, both papers)

**Statement.**
- Let a patch block be built as W' D W'^dag. Here W' is any synthesised basis change, applied with its exact
  adjoint, and D holds eigenvalue rotations that are exact (+-1 classes, or exact dyadic angles).
- Then the block is an exact reflection and the walk exactly encodes a nearby H' (**[Certain]**, checked to 1e-14).
- Now suppose the patch operator lies in the centre of the permutation algebra, as an isotropic triangle or K_n
  clique does, and the ground state is a non-degenerate SU(2) singlet.
- Then the first-order energy shift vanishes: the SU(2) twirl of i[K, A] equals i[twirl(K), A], which is 0 because
  A commutes with every permutation.
- The same holds for PREPARE-weight errors on a translation-invariant ground state, because the weight errors sum
  to zero.

**Evidence** (`results/phase2_bias.md`, Lanczos). **[Certain]**
- 4x4 J1-J2 (J2/J1 = 0.5) and 12- and 18-site kagome.
- Basis-change angle errors, gridsynth approximations and PREPARE weights all give a quadratic bias. For example,
  a block error of 1.1e-3 gives an E0 shift of 2.5e-6 on 18 sites.
- The 3-class block tri_0.25 shows a first-order bias, as expected.

**Consequence.**
- The precision per rotation can follow eps_r ~ sqrt(eps / (kappa n_p hw)) / R instead of eps / (Lambda_eff R),
  with kappa <= 0.06 measured and 0.1 used.
- **[Certain]** Rotation T per query falls 1.7-1.9x. Total T per query falls 1.16-1.21x at N <= 108 and 1.02-1.09x
  at N >= 256 (`results/phase2_fidelity.md`). The gain is modest, and largest in the near-term window the clique
  paper targets.

**Conditions to state.**
- Non-degenerate singlet ground state.
- Isotropic patches.
- The swap baseline's NN/NNN label rotation and non-isotropic 3-class patches stay first order.
- Pyrochlore tetrahedra (S = 0, 1, 2) are mixed: the basis-change rotations are second order, the S = 1 eigenvalue
  rotation is first order.
- **[Likely]** Mixed (random-compiling) synthesis breaks the exact adjoint pairing unless one sample is drawn per
  query and reused for W' and W'^dag.

**Prior art to cite.**
- Tran, Su, Carney, Taylor, PRX Quantum 2, 010323 (2021): symmetry protection of Trotter errors.
- Campbell, PRA 95, 042306 (2017), and Hastings, arXiv:1612.01011: random compiling, quadratic suppression without
  symmetry.

**[Likely]** I found no prior statement for block encodings in a literature search on 2026-10-29.

**Before inclusion (laptop):**
- one 24-site run (about 3 GB, a few hours);
- a pyrochlore 16-site check of the mixed case.

### A2. Routing parallelism is not available on factory-limited hardware (triangle paper, discussion)

- **[Certain]** IonQ's outlook says their routing networks "allow substantial parallelism", with T depth as the
  next step.
- **[Certain]** To run one swap-tree layer per logical step (about 25 cycles), the widest layer needs this many
  Litinski CCZ factories at p = 1e-3 (`results/phase4_ccz_layers.md`):

  | Lattice | Widest layer (CCZ states) | Factories | Physical qubits |
  |---|---|---|---|
  | square L = 8 | 32 | 77 | 3.6e6 |
  | square L = 16 | 128 | 308 | 1.5e7 |
  | square L = 32 | 512 | 1229 | 5.8e7 |
  | kagome L = 8 | 96 | 231 | 1.1e7 |

  With realistic factory budgets, wall-clock time is run volume divided by factory qubits, independent of the
  log-depth structure.
- This turns your single-factory, production-bound assumption into a result, and it qualifies IonQ's outlook
  directly.
- **[Certain]** One CCZ state per controlled-SWAP is also the floor: the Bell-basis rewrite plus the
  multiplicative-complexity argument; PyZX finds nothing below it.
- Size: one paragraph and a small table.

### A3. Align comparison conventions with IonQ (triangle paper, comparison section)

**[Certain]** IonQ's headline 15-17x combines conventions that differ from yours:
- the product-wise Pauli LCU baseline (Yoshioka et al.);
- 4 T per Toffoli;
- a flat rotation error of 1e-10 with repeat-until-success expected counts (0.53 log2(1/eps) + 4.86, about 22.5 T);
- leading-order PREPARE only;
- energy error 0.01 J1 N, which grows with N.

Decomposed:
- routing versus product-wise SELECT: 2.98x;
- SWAP shift plus Lambda_eff: 1.48x;
- triangles: 2.63x;
- 3x3 with QET: 1.42x (`BeyondBlockEncoding/docs/block-encoding-improvements-vs-ionq.md`).

**Proposed addition:** a convention table and one row in IonQ's setting computed with your pipeline. **[Guessing]**
That row is not computed yet. It is laptop-scale and the most useful single number for referees.

**[Likely]** Also state that at 0.01 J1 N and N = 14 400, phase estimation cannot resolve the ground state, so
their benchmark measures the trial-state energy.

### A4. Exactness obstruction (footnote, both papers)

- **[Certain]** The optimal isotropic-triangle reflection 2 P(S=3/2) - 1 has entries 1/3 and 2/3, outside
  Z[1/sqrt2, i]. By Giles-Selinger, no exact Clifford+T circuit, even with ancillas, encodes it, so its basis
  change needs synthesised rotations. This explains the coarse-angle result.
- **[Likely]** The smallest ring-compatible normalisation is 4/3 of optimal. After re-centring it costs 1.16-1.18x
  in Lambda_eff, which pays only where rotations exceed ~15% of the query cost (L <= 10).

### A5. Distillation remains necessary at the target sizes (one sentence, physical-resources section)

- **[Likely]** Coherent runs at N >= 256 and 1e-3 need 1e-10 to 1e-13 error per magic state.
- Cultivation at p = 1e-3 reaches 2e-9 and is sufficient only for N <~ 100 at 1e-2.
- Splitting phase estimation to fit costs 10-1e4x more repetitions (shot-noise model) against a 5-17x per-state
  saving.
- This supports the CCZ2T distillation model the papers already use.

## Optional, low value

- **Nested-reflection identity.** Any patch with class eigenvalues e_0 < ... < e_r equals
  c_0 + sum_k b_k (2 P_{>=k} - 1), with sum b_k equal to the half-width. This reaches optimal normalisation without
  eigenvalue rotations, but moves the same number of rotations into PREPARE. Use it at most as a footnote in the
  clique paper's K_n section.

## Not worth adding

- The CCZ-versus-4T pricing as a result: it is already in your physical estimates and is standard practice.
- Synthillation of the local block: bounded by 1.2x from L = 16 up.
- The PyZX T counts of CCZ batches: trivial merges.

## Order of work, all laptop-scale

1. Row in IonQ's setting (A3).
2. 24-site and pyrochlore checks for A1.
3. Write A1 into the error appendix and A2 into the discussion.
4. Footnotes A4 and A5.

Then move on to light-cone-local block encodings (`BeyondBlockEncoding/docs/radical-directions.md`, direction A).

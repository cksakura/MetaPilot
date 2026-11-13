# Identification v2: design

Status: design, 2026-10-24. Implements in MetaPilot what the audit harness established
(`scripts/identification_audit/`; audit report, revision 2). Principles: `identification-principles.md`
(P1-P9). Selected with `identification_version = v2`; `v1` (the current joint identification) stays
available unchanged for comparison.

## 1. What v2 changes, and why

| Id | Change | Where | Evidence (audit report, revision 2) |
|---|---|---|---|
| C1 | One query per (spectrum, precursor); closed and open candidates merge only for the same precursor; FDR groups split by precursor class | `JointCandidate`, `JointIdentification` | Precursor queries +5.7-8.1% peptides at matched FDP (R4) |
| C2 | RESET: training / estimating decoys by peptide; decoy weight 2 at every level | `JointPsmClassifier`, `JointIdentification`, `GenomeEvidenceTest` | P4; tree calibration (R3) |
| C3 | Gradient-boosted trees (XGBoost, monotone constraints) as the joint model | `SemiSupervisedPsmClassifier.ModelFactory` | Trees +3-5% peptides, +5-10% proteins, calibrated with C4 (R3, R4) |
| C4 | Pseudo-reversed sample decoys by default in v2 | `sample_decoy_method` | Closes the twin gap (R2, R5) |
| C5 | Protein filter (P8): one-peptide proteins walked to 1%, two-peptide proteins kept, peptides of rejected proteins dropped | `GenomeEvidenceTest`, razor assignment | One-peptide proteins 8.5-22% FDP before, all proteins 0.17-0.32% after |
| C6 | DIA-NN predictions to open rank 3 | `joint_ms2_max_open_rank` default 3 | Rank 3 had no prediction features; candidates are reported to rank 3 |

## 2. Theory per change

**C1, precursor queries (P3).**
- A PSM is the claim that one peptide explains one precursor. Crystal-C may assign a scan's open candidates to a different, co-isolated precursor. 98% of closed/open disagreements were such precursors.
- Two candidates for different precursors are therefore not rivals. One winner per precursor query, not per spectrum.
- A reassigned precursor carries a different error rate (about twice the overall rate), so it gets its own FDR group. Groups are the evidence class (`EXACT`, `OPEN_UNMODIFIED`, `OPEN_KNOWN_SHIFT`, `OPEN_UNEXPLAINED`) crossed with the precursor class (same, reassigned).
- The same peptide found by both searches merges into one candidate only when both explain the same precursor.
- A query key holds one reassigned precursor per scan. Crystal-C writes one per scan; the code refuses input with more.

**C2, RESET (P4).**
- Decoys that train the model are scored lower than decoys it never saw, so using them to estimate error is optimistic. Flexible models (C3) make this worse.
- Each decoy peptide (I = L) is assigned once, by a stable hash, to training or estimating.
- In training, training decoys are the only negatives. Positives are chosen among pseudo-targets (targets and estimating decoys) at q ≤ `trainQValue` against training decoys.
- In reporting, a query whose winner is a training decoy reports nothing. Every estimate counts each estimating decoy twice: group q-values, the peptide walk, the genome test (`decoyWeight = 2`) and the C5 walk.
- Ties: a tied block is accepted or rejected as one, decoys counted. Production's `qValuesConservative` already does this; the weight is added to it.

**C3, trees.**
- With realistic decoys (C4) and RESET (C2), monotone gradient-boosted trees rank better than logistic regression and stay calibrated at every level on both entrapment benchmarks.
- With reversed decoys they were optimistic: they exploited the decoy/false-target difference that C4 removes.
- Monotone increasing constraints on the score features (hyperscores, −log10 expect, hyperscore margin, predicted-spectrum entropy and matched fraction) keep the model from learning that a better match is worse.
- XGBoost (`ml.dmlc:xgboost4j`) is used through a weka `Classifier` adapter, so the semi-supervised loop is unchanged. Settings match the audit harness: depth 4, 200 rounds, learning rate 0.1, L2 1, one thread per fold, fixed seed.

**C4, pseudo-reversed decoys.**
- Whole-protein reversal moves cleavage sites, so no decoy has the mass of a present peptide. Same-mass neighbour errors are therefore invisible to it.
- Pseudo-reversal keeps each tryptic segment's length, composition and ends.
- Implemented and off by default (`sample_decoy_method`); v2 turns it on.
- **Gate:** the 12mix genome-test calibration check fails under pseudo-reversal (p = 0.016 to 5×10⁻⁵). This is investigated before v2 becomes the product default, not before the v2 test run.

**C5, protein filter (P8).**
- A false peptide rarely joins a supported protein; it usually creates a one-peptide protein. So one-peptide proteins carry nearly all protein-level error.
- `GenomeEvidenceTest` already assigns every evidence peptide, target and decoy alike, to the protein of an accepted genome that holds most evidence peptides (its picked-protein step).
- v2 adds, on that symmetric assignment:
  - proteins with at least two assigned peptides pass;
  - one-peptide target and decoy proteins are walked together by their peptide's score to 1% (decoy weight as C2), and those at or under it pass.
- Razor assignment then considers only passing proteins. A peptide with no passing protein is excluded from protein and genome outputs, as a peptide without a confident genome is today.

**C6, rank-3 predictions.** The open search reports up to three candidates per precursor, and all compete. Without predictions, rank 3 lacked the strongest features. `joint_ms2_max_open_rank` defaults to 3 in v2 (it stays a parameter).

## 3. What stays

- The closed and open searches, Crystal-C, the delta-mass peak gate and the evidence classes.
- The joint model's features (`JointFeatures`) and its cross-validation by spectrum (P5).
- The genome test's statistics (Poisson, near-neighbour exclusion, BH).
- `v1` output for the same run, through `identification_version = v1`.

## 4. Acceptance test

On both entrapment runs (`Qiao_faecal_pseudorev`, `Qiao_12mix_pseudorev`), v2 production output is compared with the audit harness (B0/Btm, RESET, pseudo-reversed decoys). Both readings are compared at every level: PSMs, final peptides, genomes, one-peptide proteins, proteins, retained peptides, and additional IDs.

v2 passes when:
1. each level's entrapment FDP is within its decoy estimate's 95% interval or below it;
2. yields at matched FDP are within 2% of the harness Btm arm, or above it;
3. all 12 mix genomes are accepted.

Differences in razor assignment (intensity in place of PSM counts) are reported, not treated as failures.

## 5. Implementation order

C6 → C2 → C1 → C3 → C5 → C4 default under v2, each with unit tests. Then the HPC test run on both entrapment runs, identification only, reusing the searches.

# Joint closed+open identification

Status: implemented (first build), not yet run on server data. Replaces the cascade as the default
when the closed stage is on (`closed_open_integration = joint`); the cascade stays selectable
(`cascade`) and is also computed inside a joint run as a comparison arm.

## 1. Goal

One identification method that uses the closed and the open search as evidence for each other, and
that beats, at the same measured error:

- the closed search alone,
- the open search alone,
- their plain PSM-level union,
- the cascade (closed first, open on the remaining spectra).

"Beats" is judged on yield against error measured without decoys (12mix: share outside the spiked
genomes; PXD050781: share never accepted in another replicate run), not on counts at each method's
nominal 1%.

## 2. Why the cascade could not win

Measured on the 12mix cascade run (`DDA/mzML/MetaPilot`, 2026-10-20):

| Observation | Value |
|---|---|
| Open winner differs from the closed answer on closed-accepted target spectra | 213,149 / 580,247 (37%) |
| Peptides the open stage adds (accepted, not accepted by the closed stage) | 5,707 (+6.7%) |
| Their decoy/target ratio | 266 / 5,707 = 4.66% |
| Their share outside the spiked genomes (protein-mapping approximation) | 7.39% |
| Peptide-union check: union FDR before tightening, factor applied to every class | 1.38%, 0.656 |
| Closed-stage peptides lost to that factor | 2,980 |

The closed stage decides each spectrum alone, the open ranking never sees the closed answer, and the
two stages score on two scales (the genome test's `pi_g` window then mixed them: chi-square
p = 1e-5). These are properties of the architecture.

## 3. Design

Code: `core/opensearch/JointCandidate`, `JointFeatures`, `JointPsmClassifier`,
`PredictedSpectrumFeatures`, `JointIdentification`, `IdentificationArms`; runner
`MagOpenSearchHeadlessRunner.runJoint`.

1. **Candidates.** Per spectrum: closed top-N (`closed_stage_topn`, default 3) and open top-N (both
   engines). Merged by (file, scan, peptide): a peptide both searches propose is one candidate with
   both bodies of evidence.
2. **Group** (decided before scoring, identically for targets and decoys):
   `EXACT` (proposed by the closed search) > `OPEN_UNMODIFIED` (open only, accepted zero family) >
   `OPEN_KNOWN_SHIFT` (open only, accepted allow-list peak) > `OPEN_UNEXPLAINED` (anything else;
   walked and reported, not an identification unless `joint_accept_unexplained`).
3. **Features** (`JointFeatures`), one row per candidate, computed on demand:
   - shared: charge, precursor mass, length, termini, missed cleavages;
   - presence: found by closed / MSFragger open / Sage open;
   - each search's own block (rank, matched ions and fraction, hyperscore, hyperscore - next,
     -log10 expect; closed: mass error after the isotope offset, isotope offset, variable mods;
     open MSFragger: delta-mass score gain, localization delta);
   - open peak context: |normalized delta|, zero family, known shift, distance in sigmas, local
     density, peak size, isotope offset;
   - spectrum competition: does the spectrum have any closed candidate, -log10 expect of the best
     *other* closed candidate, hyperscore margin over the best other candidate;
   - predicted spectrum (DIA-NN 1.8.x, `joint_ms2_features`): entropy similarity, matched predicted
     intensity, entropy margin over the spectrum's other candidates, calibrated RT error.
   Columns with no information on the data (e.g. Sage without Sage, MS2 without DIA-NN) are
   withheld automatically; `joint_withheld_features` withholds more by name.
   **Not used: same-peptide and same-protein counts.** They aggregate targets and decoys
   separately, and a false target lands on peptides/proteins that true identifications also hit:
   the count marks target identity, not correctness. On a synthetic run with known truth they took
   the exact group's measured FDP to 7.2% at a decoy estimate of 0.3%; without them 0.38%. Note:
   the existing open-search, combined and closed-stage classifiers (`OpenSearchFeatures`,
   `CombinedPsmFeatures`, `ClosedOpenCascade.CLOSED_FEATURES`) still use them - see section 7.
4. **Predicted spectrum** (`PredictedSpectrumFeatures`). The unmodified backbone is predicted (C
   carbamidomethyl only, so no modification needs a DIA-NN name) and masses are placed on the
   fragments: a closed candidate's variable modifications where the closed search put them; an open
   delta on MSFragger's best-localized residue, or, unlocalized, on the residue that best explains
   the spectrum. Top 24 predicted fragments, `joint_fragment_tolerance_ppm`. Candidates beyond open
   rank `joint_ms2_max_open_rank` (default 2) and not closed get none. RT: per file, DIA-NN iRT mapped
   to observed RT by an isotonic piecewise-linear curve through binned medians of the 3,000
   strongest rank-1 candidates (chosen by score only, never by label), error in robust SDs. DIA-NN is
   called per file in chunks of 100k precursors; a chunk DIA-NN fails on loses its own features only.
5. **Scoring** (`JointPsmClassifier`): Percolator scheme, the same helpers as the existing
   classifiers - positives re-seeded each iteration from target spectrum winners at
   the classifier training q-value, every training-partition decoy a negative, folds by spectrum, held-out rows
   scored only by models that never saw them. Logistic regression on standardised inputs. Positives
   are not restricted to validated peaks (the groups keep unexplained shifts apart already). The
   score is the log-odds from both class probabilities, `log(p_target) - log(p_decoy)`: weka keeps
   the decoy probability precise after the target probability has rounded to 1.0, so this is the
   linear predictor up to about +-700. The earlier logit of the target probability alone clamped at
   +27.63 and tied 251,164 of 1,132,941 closed winners on the 12mix cascade run.
6. **Winner:** the top-scoring candidate per spectrum, whichever search proposed it.
7. **PSM FDR:** target-decoy walk within each group, `(decoys + 1) / targets`, tied scores as one
   block (`TargetDecoyFdr.qValuesConservative`; the open groups can be small); accepted at q <= 1%.
   Every row of a candidate carries the joint q-value afterwards (the winner's group q, 1.0 for the
   others), so no row keeps a value a comparison arm wrote.
8. **Peptide FDR:** each peptide (base sequence, I = L) belongs to the highest-priority group in which
   it has an accepted PSM, scored by its best accepted PSM there, and is walked in that group alone
   (same estimator as the PSM walk).
   A peptide confirmed at exact mass anywhere is an exact peptide; a peptide only the open search
   identifies has to pass the open group's own walk. Identifications are accepted PSMs whose peptide
   passes. Each group at 1% keeps the union at 1% or below - no union tightening factor.

## 4. Outputs (combined folder)

| File | Content |
|---|---|
| `validated_psms.tsv` | as before, plus `id_group`; `q_value_scope` = `JOINT_GROUP_<group>` |
| `scored_psms.tsv` | one row per candidate: searches, each search's scores, joint score, winner, group, PSM and peptide q, MS2 columns. Read by the genome test (`id_group`) |
| `joint_identification_summary.tsv` | per group: winners, PSMs and peptides at 1%, kept; union decoy ratios; open-only peptides; features used/withheld; coefficients per fold |
| `joint_feature_report.tsv` | target/decoy AUC of every feature among winners, per group |
| `identification_arms.tsv.gz`, `identification_arm_peptides.tsv.gz`, `identification_arms_summary.tsv` | closed-only, open-only, cascade and joint winners with q <= 0.05 (`joint_compare_arms`, default on) |
| `closed_open_summary.tsv` | the cascade arm's own summary (compare mode) |

Genome test: an `EXACT` winner is the closed class and its decoys join the `pi_g` window (one score
scale); open groups map to unmodified / known shift / unexplained.

## 5. Parameters (`additional_parameters`)

`closed_open_integration` (joint | cascade), `closed_stage_topn` (3, every mode; the open search's
`round2_output_report_topn` is 3 too), `joint_ms2_features` (true),
`joint_ms2_max_open_rank` (2), `joint_fragment_tolerance_ppm` (20), `joint_accept_unexplained`
(false), `joint_compare_arms` (true), `joint_withheld_features` (""). Changing any of them, or the
mode, rewrites `identification_mode.txt` and re-runs validation; so does a change in what the MS2
features can use (DIA-NN executable present, number of runs with an mzML), so installing DIA-NN
after a run makes that run's validation stale. `closed_stage_topn` changes the closed search
parameters, so the closed searches re-run once. An unrecognised `closed_open_integration` value or
withheld feature name is warned about.

## 6. Evaluation

```
python docs/developer/scripts/evaluate_identification_arms.py <combined> --spiked 12genomes.txt
python docs/developer/scripts/evaluate_identification_arms.py <combined> --replication
```

`arms_evaluation.tsv` (each method at its own acceptance), `arms_curves.tsv` (yield vs measured
error over q), `arms_matched.tsv` (yield at the error closed-only reaches at its 1%). Pass:
joint at or above every other arm on PSMs and peptides at matched error, its open-only peptides
cleaner than the cascade's 4.66% / 7.39%, 12 of 12 spiked genomes accepted.

## 7. Open items

- **Same-peptide/protein count features in the existing classifiers.** The closed-stage, open,
  combined and stage-1 classifiers use them. On synthetic data they make the decoy estimate
  optimistic (closed classifier: 1.0% estimated, 8.8-9.2% measured). On the 12mix closed search
  (2026-10-20) no such effect is visible: the 89,086 PSMs the counts add at 1% are 2.93% outside the
  spiked genomes, against a 3.08% floor for unambiguous matches and about 32% for random target
  matches. That measure cannot see a wrong match inside a spiked genome, which is where the leak
  would put it; an entrapment test decides. Kept until then.
- **Existing classifiers still clamp the logit** (`SemiSupervisedPsmClassifier.toLogit`, used by the
  closed, open and combined classifiers): only the joint model uses `logitOf`. Switching them changes
  the cascade and open-only results and their comparison arms.
- **Known-shift gate** decides a peak from its own rank-1 decoy fraction, so a peak that just passes
  has a selected-low decoy count. On 12mix all six accepted known shifts pass by a wide margin.
- **Model family.** Logistic regression only. Tree ensembles give tied probabilities, which a
  target-decoy walk breaks arbitrarily; any nonlinear model must give continuous scores.
- **Group-walk** (Freestone, Keich & Noble 2022) is more powerful than separate walks per group; the
  separate walks are the conservative first build.
- **Feature selection** is by information content plus the AUC report; no automatic search yet.
- The open-only arm is validated with the cascade's peak settings (rejected peaks are not folded),
  so it differs slightly from a true open-only run.
- Not run on server data. Memory: the candidate pool adds one small object per (spectrum, peptide)
  to the existing open-search rows; DIA-NN predictions are held per 100k-precursor chunk.

## 8. Tests

`JointIdentificationTest`: merge and groups; separate group walks; peptide group priority with I/L
merging; winner choice; fragment masses, entropy similarity, isotonic fit; predicted-spectrum score
with localized, wrongly localized and closed variable modifications; synthetic end-to-end run with
known truth (per-group measured FDP within 3%, open-only identifications present).

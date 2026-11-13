# Unified rescoring with batch-wise FDR control (design)

> **Superseded in part (2026-10-21)** by `stagewise-exact-first-identification-design.md`. Exact and
> shifted candidates no longer compete for a spectrum in one model; they are decided in two
> sequential stages (exact first). The components below (annotator, DIA-NN block, peaks, listwise
> model, Group-walk) are reused there.

Status: design, revision 2 (2026-10-20). Supersedes the joint design in
`joint-closed-open-identification-plan.md` (kept for history). Nothing here is implemented yet.

Fixed requirements: rescoring, machine learning, spectrum prediction, peak detection on the
mass-difference axis, and FDR control per batch.

Revision 2 decisions:
- Peaks are detected and filtered first, without any database. Only the retained peaks are mapped
  to UniMod, as annotation.
- Prediction is DIA-NN only: the unmodified backbone, matched without the mass difference.
- The training instance is a spectrum with all of its closed and open candidates.

## 1. What the joint build taught us (12mix, measured 2026-10-20)

| Observation | Number | Consequence for the design |
|---|---|---|
| Largest joint coefficient: `closed_neg_log10_expect` | +22.5 (standardised) | Expect values measure the *search space*, not the match. Engine scores are not used. |
| Open-only candidates get the training mean for every closed column | all of them | No per-search feature blocks, no imputation. |
| Open **decoys** beating closed targets accepted at 1% | 9,883 spectra | One score across searches was not calibrated. |
| Closed-accepted spectra lost to an unexplained winner | 51,195 | An unreportable candidate must never delete a spectrum's identification. |
| Joint vs cascade at matched error (outside-spiked share) | 708,968 vs 787,880 PSMs | The joint build is worse than the cascade. |
| Scores tied at the clamp (old logit) | 251,164 / 1,132,941 closed winners | Use raw margins. |

From the literature:

- **Separate target-decoy competition (TDC) per group does not rigorously control the FDR.**
  Group-walk does, and gains 3.6-26% discoveries over ignoring groups (Freestone et al.,
  Bioinformatics 2022). Our cascade and joint builds both use separate walks.
- **Percolator and PeptideProphet can fail to control the FDR** (entrapment, Freestone et al. JPR
  2024). Percolator-RESET fixes the training side: train against a random half of the decoys,
  estimate with the other half, doubled (Freestone, Käll, Noble, Keich, JPR 2025).
- **CONGA** combines narrow and open searches with neighbour filtering, mass-difference groups
  (minimum 2K, adaptively merged by a KS test) and Group-walk, with entrapment-verified FDR. It
  keeps narrow and open hits in different groups so their raw scores never compete.
- **Entrapment** (FDRBench) is the valid check of FDR control (Wen et al., Nature Methods 2025).
- **Nonlinear rescorers can overfit the target database** (HUPO 2024, entrapment-tested). Any
  nonlinear model must be validated by entrapment.

## 2. Principles

1. **The search engines only propose candidates.** Every score the model sees is computed by
   MetaPilot from the spectrum and the hypothesis, the same way for every candidate. No expect, no
   engine score or rank, no per-search block, no imputed column.
2. **A hypothesis is (peptide, mass explanation).** The explanation decides its batch, not the
   search that found it.
3. **Batches come from mass alone, before any label is read.** Peak detection and filtering never
   look at target/decoy labels or at a modification database.
4. **Scores are compared for FDR only within a batch.** Group-walk combines the batches.
5. **A competing candidate is evidence, never a deletion.** A strong unexplained alternative
   lowers the confidence of a spectrum's reportable winner (a margin feature). It can't remove the
   identification the way the joint build did.
6. **The unit of learning is the spectrum.** The model decides which of a spectrum's candidates is
   right, or that none is.
7. **Training never uses the decoys that estimate the FDR** (RESET).
8. **Entrapment decides the default.**

## 3. Workflow

```
 closed search (top-3)   open search (top-3)       [existing MSFragger/Sage, same refined DB]
            \                /
  1  precursor recalibration; delta-mass peak DETECTION + label-blind FILTERING
  2  batches from retained peaks; UniMod mapping of retained peaks (annotation only)
  3  per spectrum: candidate set = closed + open hypotheses, neighbour filtering
  4  MetaPilot features per candidate: own fragment annotation + DIA-NN backbone + RT
  5  spectrum-level listwise model with a "none" option, semi-supervised, RESET, CV
  6  per spectrum: reportable winner, its score = margin over its best alternative
  7  batch-wise FDR: Group-walk over batches (PSM), then peptide units, Group-walk again
  8  outputs, genome test, proteins
```

### 3.1 Searches (unchanged)
- Closed: +-20 ppm, isotope 0/1/2, M oxidation and protein N-terminal acetyl, **top 3**.
- Open: -150/+500 Da, **top 3** (the same depth as the closed search, 2026-10-21), Crystal-C. Sage optional.
- Phase 2 option: peptide-paired decoys (each tryptic peptide reversed in place, C-terminal K/R
  kept). This gives paired peptide-level competition, at the cost of one database rebuild and one
  re-search.

### 3.2 Peak detection and filtering (no labels, no modification database)
**Input:** delta = precursor - peptide mass of every rank-1 open hit, targets and decoys pooled,
after per-run precursor recalibration (score-weighted, label-blind).

**Detect:** kernel density on the delta axis with a bandwidth from the measured precursor error,
then the existing Gaussian peak fit for centroid, sigma and area.

**Filter, in order, all label-blind:**

| # | Filter | Removes |
|---|---|---|
| F1 | Prominence over the local background density, and a minimum area (spectra) | noise bumps |
| F2 | Run support: present in >= m runs (default 2, or all runs if fewer) | single-run artefacts |
| F3 | Isotope folding: a peak at parent + k * 1.00336 (k = 1, 2) whose width matches is folded into its parent | precursor isotope errors |
| F4 | Zero folding: inside the zero window -> the exact explanation | near-zero shifts the closed search already covers |
| F5 | Fragment support: share of the peak's rank-1 hits where placing the delta on a residue gains shifted fragment ions, against the share at random off-peak deltas (one-sided binomial). Pooled over targets and decoys. | precursor-mass artefacts that don't sit on the peptide |
| F6 | Minimum batch size 2K (K = 40); smaller peaks are merged with the nearest retained peak only if within 3 sigma, otherwise they go to `UNEXPLAINED` | peaks too small to walk |

None of these reads a label, a score from the rescoring model, or a database. Every hit on a
retained peak belongs to that peak's batch. Hits on filtered-out peaks, or on no peak, are
`UNEXPLAINED`.

**Map** (after filtering, annotation only): each retained peak against:
- UniMod, single and pairwise combinations;
- isotope-offset variants;
- amino-acid substitutions.

The mapping is written to the peak report and the PSM tables. It changes no batch, score or
q-value.

### 3.3 Candidate sets
- **Candidate set of a spectrum:** the union of its closed top-3 and open top-3 hypotheses. One
  hypothesis per (spectrum, peptide); a peptide both searches propose is one hypothesis.
- **Explanation**, from mass alone:
  - `EXACT`: proposed by the closed search, or |delta| inside the zero family;
  - `PEAK_p`: delta on retained peak p;
  - `UNEXPLAINED`.
- **Neighbour filtering (CONGA):** within a spectrum, drop a hypothesis whose theoretical b/y ions
  share more than tau = 5% with a higher-ranked one. Ranking is by the label-free annotation score
  of 3.4, which removes hitchhikers (I/L variants, close homologues).

### 3.4 Features of a candidate
Every candidate is annotated by MetaPilot against the observed spectrum.

| Block | Features |
|---|---|
| Fragment match, unshifted | b/y ions of the backbone (charge <= min(z-1, 2)): matched peaks; matched intensity / TIC; matched intensity / top-50 intensity; b and y coverage; longest consecutive b and y run; mean and SD of fragment ppm error; top-10 peaks unexplained |
| Fragment match, shifted (`PEAK_p`, `UNEXPLAINED`) | best site for the delta: shifted ions gained, best - second-best site |
| Precursor | explanation error in sigmas (exact: ppm / run sigma; peak: distance to centroid / peak sigma); isotope offset; log(1 + abs(delta)) |
| DIA-NN prediction | backbone only (C carbamidomethyl; variable modifications and the delta not placed). Matched to the observed spectrum without the delta: spectral-entropy similarity, Pearson on matched top predicted fragments, predicted intensity observed. RT error in robust SD after per-run isotonic calibration on score-chosen anchors. |
| Peptide | length, charge, missed cleavages, termini |
| Explanation | type one-hot; log batch size (pooled targets + decoys) |
| Spectrum context | see 3.5 |

**Ignoring the delta in the prediction match is deliberate for this stage.** Fragments carrying
the shifted residue don't match, so a shifted hypothesis gets a systematically lower similarity.
That is acceptable because FDR is decided within a batch, where every hypothesis carries the same
kind of shift. Placing the delta is a later improvement, measured on its own.

**Not used:** expect, engine scores and ranks, same-peptide/protein counts, anything a decoy
cannot also have. A failed prediction gives a neutral value plus an indicator, identically for
targets and decoys.

### 3.5 Learning: one instance = one spectrum and all its candidates
**Model: listwise choice with a "none" option.** For spectrum sigma with candidate set C_sigma,
each candidate j gets a utility u_j = f(x_j, z_sigma,j), where:
- x_j is the candidate's own features (3.4);
- z_sigma,j is its **context in the spectrum**, computed without labels:
  - rank of j by annotation score within the spectrum;
  - margin to the best other candidate, overall and within each explanation type;
  - number of candidates within 10% of the best;
  - whether the closed and the open search agree on the top peptide;
  - the best other candidate's explanation type.

A virtual candidate "none" has fixed utility u_0. The model gives

  P(j is correct | sigma) = exp(u_j) / (exp(u_0) + sum_k exp(u_k)),

which is learned as XGBoost pairwise ranking (`rank:pairwise`, query = spectrum), with "none" as
an extra row of constant features in each query.

**Semi-supervised labels per iteration**, with Percolator-style re-seeding and the RESET split:
- Each decoy is assigned once, at random with probability 1/2, to *training* or *estimating*.
- **Positive spectra:** the current winner is a target at pseudo-target q <= train q. Target > none
  > every other candidate of that spectrum.
- **Negative spectra:** the current winner is a training decoy. None > every candidate. These
  spectra teach that the best candidate of an unidentifiable spectrum, target or decoy, is still
  below "none".
- Every other spectrum, including all estimating-decoy spectra, is left out of training.

**Folds by spectrum:** every candidate of a spectrum is scored by a model that never saw that
spectrum. Monotone constraints on the matched-intensity, coverage and prediction-similarity
features.

**Why this shape:**
- The model learns what makes one candidate beat the others in the same spectrum. It sees closed
  and open candidates side by side, as the spectrum actually presents them.
- The "none" row anchors the utilities, so they are comparable across spectra, which the
  threshold needs.
- Because context is a feature, a strong unexplained rival lowers the winner's confidence instead
  of replacing it.

### 3.6 The spectrum's winner and its score
- **Winner** = the candidate with the highest utility among the **reportable** explanations
  (`EXACT`, `PEAK_p`).
- `UNEXPLAINED` wins only a spectrum that has no reportable candidate. It is then reported in its
  own batch, if that batch is reported at all.
- **Score:** the winner's log-odds against its best alternative, *any* explanation or "none":

  s_sigma = u_win - log(exp(u_0) + sum_{k != win} exp(u_k)).

  An unexplained rival therefore lowers s; it never deletes the winner.
- The winner's batch is that of its explanation.
- **Target-decoy symmetry:** every step above depends only on label-free features, so an
  incorrect winner is equally likely to be a target or a decoy. That is the TDC assumption.
  Entrapment checks it.

### 3.7 Batch-wise FDR
**PSM level:**
- Group-walk over the reportable batches (`EXACT`, each retained-peak batch, optionally
  `UNEXPLAINED`) on s_sigma, with the RESET estimate (2 * estimating decoys + 1) / targets.
- The frontier moves to the batch with the most decoy wins in the last K (K = 40).
- A batch is in or out *before* the walk. Dropping a batch afterwards is not valid.

**Peptide level:**
- Unit = base sequence (I = L), scored by its best accepted PSM.
- A unit's batch is the highest-priority batch where it has an accepted PSM (`EXACT` first).
- Group-walk again over the batches.
- Phase 2 replaces this with paired target-decoy competition at the stem level.

### 3.8 Outputs and downstream
- PSM and peptide tables with batch, peak id, UniMod mapping of the peak (annotation), q-values,
  and the winner's best alternative.
- Peak report: detected, each filter's verdict, retained, mapping.
- The genome test reads the discovered peptide units and their batches, with the `pi_g` window
  drawn from the batch-wise decoys.

## 4. Validation (decides the default)
1. **Entrapment on 12mix** (FDRBench, combined + paired estimators), with a foreign proteome at 1:1
   in the refined databases (one re-search). The FDP per batch and overall must be <= alpha + noise
   at 1% and 5%.
2. **Matched-error comparison** against closed-only, open-only, union and cascade
   (`evaluate_identification_arms.py`), at PSM and peptide level.
   - Pass: at or above every arm on both.
   - Its new peptides must be cleaner than the cascade's 7.5% outside the spiked genomes.
3. Replicate runs (PXD050781): unreplicated share.
4. Ablations, each judged by entrapment:
   - no DIA-NN block;
   - pointwise instead of listwise model;
   - linear instead of trees;
   - no RESET;
   - separate walks instead of Group-walk;
   - no context features.

## 5. Build order
| Phase | Content | Needs a re-search |
|---|---|---|
| 1 | Peak detection + filters F1-F6 + post-hoc UniMod mapping; MetaPilot annotator; DIA-NN backbone block; listwise XGBoost with "none", RESET and CV; Group-walk (Java port of the CRAN `groupwalk` logic); unpaired peptide walk. Offline harness on the cached 12mix searches. | no |
| 2 | Peptide-paired decoys and stem-level competition, closed top-3, entrapment database. | once |
| 3 | Placing the delta on predicted fragments, chimeric second peptides, RESET two-decoy mode. | no |

## 6. Scale (12mix)
- ~1.25M spectra, ~6.6M candidates (about 5 per spectrum) plus one "none" row each.
- **Annotation:** ~50 fragments per candidate by binary search; minutes with threads. mzML is read
  one file at a time (about 1-2 GB of scans).
- **DIA-NN:** ~3-4M backbone precursors in 100k chunks. Predictions are per peptide, not per
  hypothesis, so open and closed candidates share them.
- **XGBoost ranking:** about 1M training queries per fold, capped by sampling whole spectra.

## 7. Open questions
- The RESET correction inside Group-walk (decoy weight 2) follows the logic of both papers but
  isn't published as one procedure; confirm by entrapment or use the two-decoy mode.
- Whether `UNEXPLAINED` is reported is a product decision.
- Phase 1's peptide FDR is the unpaired walk; phase 2 adds pairs.

## Sources
- Group-walk: https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii82/6701992 ;
  https://search.r-project.org/CRAN/refmans/groupwalk/html/group_walk.html
- CONGA: https://www.biorxiv.org/content/10.1101/2023.05.02.539167v3.full
- Percolator-RESET: https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00742 ;
  https://github.com/freejstone/Percolator-RESET
- Decoy-based FDR re-investigated: https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00902
- Entrapment / FDRBench: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12240826/
- Nonlinear rescoring overfitting (HUPO 2024):
  https://programme.conventus.de/hupo-2024/posters/67a3062a-89b3-4590-aa55-91c2be9d9075
- Spectrum + sequence rescoring: https://www.biorxiv.org/content/10.1101/2024.08.19.606805v2.full.pdf

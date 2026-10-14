# DDA MAG Open-Search Workflow — Design, Rationale and Open Problems

**Status:** working prototype, GUI-only, under active evaluation.
**Audience:** reviewers evaluating the design. This document is written to be
criticised — the "Open problems" section is the point, not an appendix.

**Purpose of this document:** MetaPilot's other DDA MAG workflows delegate
search and post-processing to FragPipe/Philosopher or pFind. This workflow
drives MSFragger directly and implements its own PTM validation. Everything
below the search engine is new code and therefore new risk.

---

## 1. What problem this solves

MetaPilot analyses metaproteomics data against MAG (metagenome-assembled
genome) catalogues such as MGnify's MGYG. The existing FragPipe-based DDA MAG
workflow supports open search nominally, but **open-search results cannot be
quantified through it**: FragPipe's downstream quantification (IonQuant via
`msfragger.mass_diff_to_variable_mod`) does not carry open-search delta masses
into the quantified peptide list, so open-search peptides are silently dropped
from final output. The open-search mode is therefore not genuinely available.

This workflow replaces the FragPipe/Philosopher stack with:

- MSFragger invoked directly (`java -jar MSFragger.jar params.params *.mzML`)
- a custom delta-mass peak model
- a Percolator-style semi-supervised classifier
- target-decoy FDR and a per-modification statistical gate

The MAG-specific services (genome selection, protein inference, razor
assignment, taxonomy) are **shared, unmodified** with the other DDA MAG
workflows.

---

## 2. Pipeline overview

Two iterative search rounds, as in every MetaPilot MAG workflow. Both rounds
are **open** searches. That is a design decision, not a setting.

```
Step A  Build marker target-decoy FASTA (1:1, rev_)
Step B  Round 1: open search vs marker DB, per sample (one batch invocation)
          -> fit ONLY the zero peak
          -> take unmodified PSMs at 1% FDR
          -> genome selection (shared MAG service)
          -> build refined per-sample FASTA (1:1 decoys)
Step C  Round 2: open search vs each sample's refined DB (per-sample invocation)
Step D  Pooled validation across all samples:
          -> read all pepXMLs
          -> discover multi-scale KDE modes from spectrum-balanced target excess
          -> require cross-file support and infer conservative isotope families
          -> preserve observed mass/mode; normalize family members to their root
          -> semi-supervised classifier scores every candidate
          -> collapse to best candidate per spectrum
          -> target-decoy FDR
          -> fitted peaks: per-peak gate + local FDR
          -> no fitted peak: identification-only at group-local target-decoy FDR
Step G  Post-analysis: shared MAG services -> final_peptides / razor_proteins /
        all_proteins / id.tsv
```

### Result layout

```
<result>/opensearch_MAG_search/
    search_results/<sample>/
        marker/     <sample>.pepXML, fragger_open.params,
                    <sample>_unmodified_psms.tsv,
                    marker_genome_{evidence,ranked,selected}.tsv,
                    <sample>_refined*.fasta
        sample_db/  <sample>.pepXML, fragger_open.params
    combined/       delta_mass_peaks.tsv, scored_psms.tsv, validated_psms.tsv,
                    peak_validated_psms.tsv, unexplained_psms.tsv,
                    modification_verdicts.tsv, id.tsv
<result>/final_peptides.tsv, razor_proteins.tsv, all_proteins.tsv
```

Folder names match the other MAG engines deliberately, so result folders look
the same whichever engine produced them and the cleaner's existing rules apply
without special cases. Each search folder holds its own `.params` file so it is
self-describing.

---

## 3. Design decisions and their reasoning

### 3.1 Round 1 is an open search, but only its zero peak is used

Round 1 exists to pick genomes, and only unmodified peptides are used for that
— so a closed search looks sufficient. It is not, and the difference is about
purity, not yield.

In a **closed** search, a spectrum from a modified peptide either fails to
match or gets force-matched to whatever wrong peptide falls inside the narrow
precursor window. That wrong match can pull an incorrect genome into the
refined database. In an **open** search, the same spectrum matches its true
peptide carrying a mass shift, so it is diverted *out* of the zero peak rather
than contaminating it. The surviving zero-peak set is therefore purer.

Round 1 also reads **rank-1 hits only**: it never arbitrates between competing
candidates, so lower ranks are pure overhead there. Round 2 reads all ranks.

### 3.2 1:1 target:decoy, `rev_` only

Earlier MetaPilot workflows used `shuf_` decoys because FragPipe does not emit
decoy peptides in its quantified list, and genome-level FDR needs them. Driving
MSFragger directly removes that constraint. 1:1 reversed decoys, one scheme,
both rounds. **The FragPipe MAG workflow's combined stage uses 1:3 decoys** —
copying it here would have made FDR ~3× wrong.

### 3.3 Isotope modes are evidence; family members normalize to one root

A precursor triggered on another isotope yields a delta mass displaced by
approximately `k * 1.00335484` Da. These are precursor-selection modes, not
independent modifications. The raw mass difference and observed mode are kept
for audit, while validation and downstream quantification use a normalized
mass and one canonical family root. Round 1 therefore includes confidently
linked zero-family satellites after normalization instead of calling them
separate modifications or discarding their peptide evidence.

This depends on `isotope_error = 0` in the MSFragger parameters. With isotope
correction enabled MSFragger folds those spectra into the zero peak internally
and the distinction becomes unrecoverable. The runtime check
`countNearIsotopeOffsets` warns if no PSMs appear near ±1.003/±2.006, which is
the signature of that misconfiguration.

Isotope families are bounded stars, not connected components. Connected
spacing graphs percolated through dense, unrelated integer-spaced
modifications and produced impossible offsets as large as +21. A satellite can
never recruit another mode, offsets are limited to ±3, and the zero peak is
always its own root.

**Round 1 uses a tighter bound than round 2: ±2, not ±3**
(`round1MaxIsotopeOffset`, clamped to `maxIsotopeOffset`). The two stages fold
on different evidence. Round 2 links a satellite only after consecutive-satellite,
amplitude-monotonicity and spacing-residual tests, so a generous bound is
cheap. Round 1 folds on mass proximity alone, so each additional offset is
another chance to absorb an unrelated shift. Measured on the 12mix marker
stage: monoisotopic 93.9% correct, +1 94.1%, +2 91.8% — declining. ±3 had never
been measured; it was inherited from `maxIsotopeOffset` rather than chosen, and
contributed ~2% of the round-1 set. Widening it again requires a correct-ratio
measurement for the offsets being admitted, not a plausibility argument. For a non-zero family, centroid spacing alone is
insufficient: at least two consecutive satellites must be present and KDE
amplitude may increase by at most 25% as distance from the root grows. The
spacing residual must also be within
`min(10 mDa, max(1 mDa, 4 * hypot(width1,width2)))`. This explicitly keeps
+0.984016 Da deamidation separate from zero +1 isotope (19.3 mDa apart).

### 3.4 No hard filters in round 2 — every condition is a classifier feature

The governing principle: the classifier, not the engine and not a threshold,
decides. Peak membership, distance to centroid, local density, isotope offset
are **features**, not filters. A PSM outside every peak gets `peak_id = -1`,
which is a feature value, not a rejection.

Three deliberate exceptions:

- **Round 1 selection is a hard cut.** Its job is genome mapping, not PTM
  discovery, and running the classifier there would be circular — a
  well-trained classifier wants the refined database that round 1 exists to
  produce.
- **Spectrum-level collapse before FDR** (§3.6) is not a filter on evidence;
  it sets the unit of analysis.
- **Final reporting is tiered.** Peak membership does not prevent classifier
  scoring, but only an accepted fitted peak can support a modification claim.
  An out-of-peak winner may remain an identification with unexplained mass.

### 3.5 Semi-supervised classifier

Percolator-style: seed confident positives from targets passing a q-value
threshold, decoys as negatives, iterate. Cross-validated by fold, **split by
spectrum** so a spectrum's rank-1 and rank-2 candidates cannot straddle the
train/test boundary (near-duplicate features would leak the answer).

This deliberately avoids four defects in the pre-existing
`core.ml.SemiSuperClassifier`:

1. **Inverted bootstrap** — it re-seeds positives from targets with *low*
   classifier scores, training on the worst targets as positives.
2. **Non-determinism** — folds assigned with `Math.random()`.
3. **Evaluation on training data** — every PSM scored by a model that saw it.
4. **A non-ensemble ensemble** — every decoy added to every fold, so the models
   share almost identical negatives.

36 features: search scores (hyperscore, nextscore, expect), raw and
isotope-normalized delta mass,
localization scores, peak context, charge, missed cleavages, aggregate counts.

### 3.6 Spectrum-level collapse before FDR

With `output_report_topN = 5`, one spectrum contributes up to five rows.
Computing FDR over all of them is wrong: decoy counts stop being independent
false discoveries, accepted counts report *candidates* rather than
identifications, and two candidates from one spectrum could both be accepted.

The classifier still scores every candidate and is free to promote a rank-3 hit
over the engine's rank-1. Only afterwards is the set collapsed to the best
candidate per `file:scan`, and FDR computed on that. Same step Percolator and
MS²Rescore perform.

Simulation (20k spectra × topN 5, correlated lower-rank targets): per-candidate
FDR reported 1% while the true FDR was **2.65%**; per-spectrum reported 1% at a
true **1.59%**.

### 3.7 Per-peak gate + local FDR

Each delta-mass peak is a hypothesis: *is this mass shift a real modification,
or a cluster of bad matches?* Under the null with 1:1 decoys, targets and
decoys are exchangeable, so an exact binomial test on (targets, decoys) is the
natural test. Benjamini-Hochberg across peaks. Peaks with ≤ N decoys are
auto-accepted (absence of decoys is itself the evidence).

Within an accepted peak, a monotonised local-FDR walk on classifier score
decides how far down to cut.

`peak_id = -1` is explicitly excluded from the peak gate. It is not a peak and
therefore not one mass-shift hypothesis. A separate target-decoy q-value within
that population controls identification error only; passing it is not evidence
that the masses share a mechanism.

**The gate p-values are screening statistics, not calibrated.** Peaks are
defined by fitting a target-dominated histogram, so conditioning on "is a peak"
pre-selects for target enrichment. FDR control comes from the local-FDR walk,
not the gate. This is stated in the output file header.

Peak discovery itself uses spectrum-balanced weights: all target candidates
from one spectrum share one vote. This preserves lower-ranked evidence without
letting `output_report_topN` multiply the apparent support for a mass mode.

### 3.8 Pooled, not per-sample, validation

Decided on three grounds:

- **FDR** is a property of the collection reported. `final_peptides.tsv` spans
  all samples. Filtering ten samples at 1% each and unioning gives a union
  well above 1%: true peptides repeat across samples and collapse on union,
  false ones are random and accumulate.
- **The peak gate is a counting test.** 40 targets / 2 decoys pooled gives
  p ≈ 1e-9; split ten ways that is 4 / 0.2 per sample and the test has no
  power at all.
- **The classifier needs confident positives** to seed on.

The one real argument for per-sample fitting — per-run mass calibration offsets
smearing a pooled histogram — is largely removed by `calibrate_mass = 2`.
**This is an assumption that should be verified**, see §6.

---

## 4. Performance work (and why it was also a correctness fix)

The fitting step originally ran for hours. Three separate causes:

**Linear scans over a sorted list.** Nearest-peak assignment was O(peaks × PSMs)
≈ 10¹⁰ operations; isotope linking was O(peaks²). Both are now binary searches
over the centroid-sorted peak array, ~780× fewer operations.

The subtlety: proximity is measured in **sigmas**, but binary search orders by
**daltons**, and peaks have different widths — so the dalton-nearest peak is
often not the sigma-nearest. Checking only adjacent peaks would silently change
results; measured, the true nearest peak is not a dalton-adjacent neighbour in
**33.7%** of cases. The search walks outward with a provably safe bound: once
dalton-distance ÷ *widest* sigma exceeds the best sigma-distance so far, nothing
further on that side can win. Verified identical to the linear scan on 90,000
randomised assignments and 345,252 isotope lookups.

**Fitting noise.** `minPeakCount = 3` gates the *apex bin*. Over 650,000 bins
holding ~1M PSMs (~1.5 per bin), noise clears that constantly — one run
detected **24,491 "peaks"**.

This is not just slow, it is **statistically destructive**. With 1:1 decoys a
peak of n PSMs has best-case binomial p = 0.5ⁿ — the p-value it gets if *every*
PSM in it is a target. So:

| PSMs in peak | best-case p | |
|---|---|---|
| 3 | 1.3e-01 | cannot reach p<0.05 ever |
| 7 | 7.8e-03 | cannot survive BH |
| 16 | 1.5e-05 | testable |

A peak below ~7 PSMs cannot be reported no matter what. Worse, each one
inflates *m*: at m = 24,491 a genuine peak needs ≥22 all-target PSMs to stay
callable; at m = 200 it needs 15. **Fitting noise suppresses real
modifications.**

Hence `minPeakArea = 15`, filtering on segment *area* before fitting, defaulted
from the gate's own detection floor rather than by taste. Simulation on a
realistic histogram: 27,564 segments → 60 fits, **12/12** planted modifications
recovered including the weakest at 60 PSMs, BH burden 27,564 → 60.

Plus `maxPeaksToFit = 500` (rank by area, largest kept) to bound the worst case,
`maxFitIterations = 200`, and parallel fitting.

---

## 5. Empirical results from the first complete run

Three samples, MGYG, ~294k validated PSMs, 68 fitted peaks.

### 5.1 Peak shape is not Gaussian

```
zero peak, 1 mDa bins
  -3 mDa   3,386  ######
  -2 mDa   8,916  ################
  -1 mDa  28,802  ####################################################
  +0 mDa  32,997  ############################################################
  +1 mDa   8,446  ###############
  +2 mDa   3,725  ######
  +3 mDa   1,809  ###
 +14 mDa     171
```

A sharp ~3 mDa FWHM core with heavy tails. Empirical sigma ≈ 1.27 mDa; fitted
sigma 0.9–1.06 mDa — **the fits are correct**. But `±3σ` then cuts membership at
±3 mDa and discards the tails, which are the same peptides measured with
slightly worse mass accuracy.

**Fix implemented:** membership is now `max(3σ, membershipFloorDa)` with a
0.010 Da default. Rationale: the relevant scale for membership is instrument
mass accuracy, not fit geometry.

| floor | PSMs reabsorbed | T:D of reabsorbed |
|---|---|---|
| 5 mDa | 8,420 (5.0%) | 221:1 |
| 10 mDa | 16,160 (9.5%) | 204:1 |
| 20 mDa | 26,323 (15.5%) | 129:1 |

Reference: fitted peaks 299:1, unassigned pool 100:1. At 5–10 mDa the
reabsorbed PSMs look like peak members; past ~15 mDa the ratio decays toward
background.

### 5.2 The `peak_id = -1` population is a mixture, not one modification

In the 2026-09 three-file run, **156,826 of 293,580 validated target PSMs
(53.42%)** had `peak_id = -1`. The pre-local-FDR group contained 223,560
targets and 45,557 decoys. Pooling all of those rows into a single binomial
test made the group pass because it was large; it did not establish one common
mass-shift mechanism.

The evidence says the population is neither all random noise nor all missed
Gaussian peaks:

| observation | result | interpretation |
|---|---:|---|
| unassigned retained targets / decoys | 156,826 / 1,568 | many peptide IDs are real at the stated FDR |
| target histogram correlation across files (10 mDa bins) | 0.991 | strong reproducible structure, not independent noise |
| farther than 50 mDa from every fitted centroid | 80.91% | most are not merely fitted-peak tails |
| conservative reproducible narrow residual modes | 57 modes / 15,006 target winners | a minority are plausible missed peaks |

Examples include +2.0055 Da (likely a two-isotope precursor-selection
artefact), +21.9795 Da (compatible with Na-H), and -18.0105 Da (compatible
with water loss). These examples also showed why a zero-baseline Gaussian is
too strict: the -18.0105 Da mode improved from R²=0.507 to R²=0.905 after a
constant local background was included.

**Decision implemented:** `peak_id = -1` is no longer sent through the
modification peak gate. It may survive its identification-level target-decoy cutoff
as an identification with unexplained mass, but it is never reported as a
validated modification hypothesis. The two populations are written separately
to `peak_validated_psms.tsv` and `unexplained_psms.tsv`; their backwards-compatible
union remains `validated_psms.tsv` for genome/protein coverage.

### 5.3 Evidence-driven upgrade plan

The immediate upgrade is deliberately narrower than claiming that every
repeatable mass mode is chemistry:

**Current default discovery path:** a Gaussian-convolved binned KDE is
evaluated at 1, 2 and 4 mDa bandwidths. Modes must persist across at least two
scales, recur in at least two files, have at least three spectrum-equivalent
target votes per contributing file, and exceed local decoy density. Empirical
KDE valleys define membership, capped at 50 mDa. Gaussian fitting remains a
compatibility mode; the Gaussian-specific safeguards below still apply there.

1. **Implemented — unbiased support unit.** Each spectrum contributes a total
   vote per target/decoy class to peak discovery, divided across candidates in
   that class. Candidate ranks 4–5 remain available to the classifier, but a
   spectrum with five correlated candidates can no longer create five votes
   for a peak or for its decoy background.
2. **Implemented — signal plus background.** Fit
   `background + amplitude * Gaussian(mean, sigma)` to the raw local histogram.
   Retain the original zero-background fit as a numerical fallback. Export the
   fitted background so acceptance can be audited.
3. **Implemented — conservative isotope-pattern families.** Observed mode id,
   raw mass, canonical family id, normalized mass and isotope offset are stored
   separately. Families are bounded, non-transitive stars. Non-zero roots need
   at least two coherent satellites with non-increasing amplitude; the known
   zero root is the only single-satellite exception. Validation tests the
   family once instead of testing each satellite as a separate modification.
4. **Implemented — honest reporting tiers.** Only fitted peaks receive peak
   verdicts/local FDR. Unassigned rows use a separate identification-only
   target-decoy q-value. `scored_psms.tsv` now separates peak
   hypothesis acceptance, modification validation, identification-only status,
   and the reporting union.
5. **Implemented after the first rerun — identifiable peak shapes.** Require at
   least seven raw bins per four-parameter fit, `sigma >= 0.5 * binWidth`, and
   `amplitude / sqrt(background + 1) >= 5`. The rerun produced 31 accepted fits
   narrower than 0.5 mDa even though the histogram bins are 1 mDa, plus 8 fits
   whose fitted background exceeded their peak height. These fits can have high
   R² but are not resolved measurements. The defaults retain the zero, +1/+2
   isotope, -18.0105 and oxidation modes in the observed run.
6. **Implemented after the first rerun — resolution-aware consolidation.** Two
   alternative fits are treated as unresolved when their centroid separation is
   below `max(2 * binWidth, 3 * hypot(sigma1, sigma2))`; the stronger integrated
   signal proxy is kept. The 10 mDa membership floor is deliberately not used
   here because assignment tolerance is not instrumental resolving power.
7. **Implemented after the uncapped diagnostic — no low-decoy bypass.** Every
   peak must pass the exact binomial/BH gate, including peaks with zero or one
   decoy. The exact test already treats an absence of decoys as evidence while
   still requiring enough targets; auto-accepting low-decoy groups admitted
   examples as small as 3 targets / 1 decoy in the uncapped run.
8. **Implemented after unknown-mass auditing — interval-complete assignment.**
   Assignment chooses the smallest sigma-distance among every peak whose
   `max(3 sigma, membershipFloorDa)` interval contains the PSM. Choosing the
   global sigma-nearest peak first and checking only its interval left 382
   retained target PSMs marked unexplained even though another fitted peak was
   within the 10 mDa floor.
9. **Independent validation before publication.** Discover peak hypotheses on
   two files and test them on the held-out file, rotating all files. Add an
   entrapment database benchmark. Gaussian fit, reproducibility, localization,
   and target-decoy evidence answer different questions and must all be shown.
10. **Quantification contract.** Build `quant_psms.tsv` only from stable canonical
   analyte keys. Unmodified/isotope-corrected peptides and uniquely localized
   modifications may be site-level analytes; ambiguous localization remains
   mass-shift-level; unexplained-mass IDs do not become modification analytes.

### 5.4 Closed-loop result after KDE and isotope normalization

On the same three pepXML files used for iteration 04, the conservative
500-mode configuration changed the validated target count only from 293,670 to
294,177 while moving 87,204 targets from the unexplained tier into validated
mass families. Unexplained targets fell from 145,771 (49.6%) to 59,074 (20.1%).
The overall decoy/target percentage changed from 0.712% to 0.722%; no accepted
family had fewer than three target winners in any file. The zero family
contained 112,597 target PSMs (38.3% of validated targets), including its
normalized ±1/±2/±3 modes.

An exploratory 1,500-mode run reduced unexplained targets further to 45,844
(15.6%), but total validated targets fell by 601 versus the 500-mode run and 53
accepted roots had fewer than three target winners in at least one file. The
default therefore remains 500: the extra weak hypotheses improve assignment
coverage but fail the current robustness guardrail. Both runs remain as audit
artifacts under `test/open/closed_loop/`.

---

## 6. Open problems

**(a) Resolved operationally: the `-1` mixture is not a modification
hypothesis.** It is retained separately at identification-level FDR for peptide,
protein and genome coverage. It cannot enter a modification verdict or a
site-specific quantification key. See §5.2.

**(b) Residual modes still need independent confirmation.** The background-aware
fit is expected to recover some narrow modes missed by the zero-baseline model,
but a narrow reproducible peak can still be an isotope-selection error, adduct,
in-source loss or search artefact rather than a biological PTM. Held-out-file
recovery and entrapment remain publication blockers, not optional polish.

**(c) Resolved: delta-mass localization evidence is read from `<ptm_result>`.**
MSFragger writes this evidence in a sibling element, not in `<search_score>`.
The reader now consumes both the legacy compact attributes
(`best_score_with_ptm`, `score_without_ptm`, ...) and the expanded MSFragger
4.1+ attributes (`score_best_position`, `score_all_unshifted`, ...), and logs
the actual `ptm_result` vocabulary through the progress callback.

**(d) Per-sample search spaces differ.** Percolator pools files that searched
the *same* database; here each sample has its own refined FASTA. A sample with
300 genomes has more chances for a high-scoring random match than one with 50,
so a single pooled threshold is mildly conservative for small databases and
anti-conservative for large ones. Unquantified. Standard treatment would be to
give the classifier search-space size as a feature.

**(e) `minPeakArea = 15` is derived from the gate but untested against real
low-abundance modifications.** It cannot be validated without ground truth.

**(f) Resolved operationally, not scientifically: sub-bin fits are rejected.**
The first upgraded rerun yielded 185 fits (169 accepted), with median sigma
0.673 mDa and 31 accepted fits below 0.5 mDa. A width below half of the 1 mDa
histogram bin cannot be estimated reliably. The fitter now enforces a 0.5-bin
sigma floor and at least seven input bins. The exact thresholds remain
provisional until external/entrapment validation; changing `binWidth` must scale
the sigma floor with it.

**(g) The classifier is now the bottleneck.** Up to 3 folds × 10 iterations = 30
logistic regression fits over ~667k instances × 36 features, feature table
~625 MB. Each iteration also scores *all* instances including the held-out fold,
which is not needed until the fold ends — ~⅓ waste at 3 folds.

**(h) Resolved: `scored_psms.tsv` is now written.** It contains every engine
candidate plus classifier scores, spectrum-winner status, spectrum q-value,
peak-hypothesis status, local q-value, modification-validation status,
identification-only status and final reporting status. Decoys remain in this
audit table and the PSM audit outputs, but are excluded from `id.tsv` and all
biological peptide/protein/taxon reports.

**(i) Sage as a second engine** is parameterised but not implemented. Intended
to be measured against the single-engine baseline, not assumed beneficial.

---

## 7. Key parameters

| parameter | default | rationale |
|---|---|---|
| `precursor_mass_lower/upper` | −150 / +500 Da | matches FragPipe MAG open search |
| `isotope_error` | 0 | **load-bearing**, see §3.3 |
| `calibrate_mass` | 2 | makes pooled fitting defensible, §3.8 |
| `output_report_topN` | 5 | classifier arbitrates, not the engine |
| `localize_delta_mass` | 1 | mandatory; see §6(c) |
| `membershipSigmas` | 3.0 | |
| `membershipFloorDa` | 0.010 | §5.1 |
| `useKdePeakDiscovery` | true | non-parametric mode acceptance, §5.3 |
| `kdeBandwidthDa` / scales | 0.001 / 1,2,4 | multi-scale persistence |
| `kdeMinPersistentScales` | 2 | reject bandwidth-specific fluctuations |
| `kdeMinReproducibleFiles` | 2 | independent-file support |
| `kdeMinFileSupport` | 3 | spectrum-equivalent votes per file |
| `kdeMaxMembershipDa` | 0.050 | cap empirical KDE basins |
| `kdeMaxModes` | 500 | real-data guardrail, §5.4 |
| isotope spacing link | min(10 mDa, max(1 mDa, 4 combined widths)) | §3.3 |
| `maxIsotopeOffset` | 3 | bounded family, no transitive growth |
| `minIsotopeFamilySatellites` | 2 | spacing alone is insufficient |
| `isotopeMaxAbundanceIncrease` | 1.25 | coherent outward amplitude pattern |
| `minPeakArea` | 15 | derived from binomial floor, §4 |
| `maxPeaksToFit` | 500 | bounds legacy Gaussian fitting only |
| `peakMinFitBins` | 7 | ≥3 residual degrees for four-parameter fit, §5.3 |
| `peakMinSigmaBins` | 0.5 | reject sub-bin-width fits, §5.3 |
| `peakMinSnr` | 5.0 | measurable prominence above fitted background, §5.3 |
| `peakMinSeparationBins` | 2.0 | absolute resolution floor, §5.3 |
| `peakMinSeparationSigmas` | 3.0 | consolidate width-unresolved alternatives, §5.3 |
| `classifierFolds` | 3 | |
| `localFdr` | 0.01 | |

---

## 8. Where the code lives

```
metapilot-core/.../core/opensearch/
    DeltaMassKdePeakFinder     target-excess multi-scale KDE + empirical valleys
    DeltaMassPeakFitter        KDE/Gaussian orchestration, isotope families,
                               normalization, binary-searched assignment
    UnmodifiedPsmSelector      round-1 zero-family normalization + FDR
    OpenSearchPepxmlReader     streaming pepXML reader (dom4j ElementHandler + detach)
    OpenSearchFeatures         36 classifier features
    SemiSupervisedPsmClassifier
    SpectrumLevelCollapse      best candidate per spectrum
    TargetDecoyFdr             monotonised q-values
    ModificationValidator      binomial gate + BH + local FDR
    OpenSearchLayout           folder layout + resume sentinels
    OpenSearchMetaPsm/Peptide  adapters to shared MAG services

metapilot-core/.../core/msfragger/MSFraggerHeadlessRunner
metapilot-core/.../core/runner/MagOpenSearchHeadlessRunner
metapilot-core/.../core/runner/OpenSearchPostAnalysisTsvExporter
metapilot-core/src/main/resources/msfragger/open_fragger_standalone.params
metapilot-windows-gui/.../dda/mag/gui/MetaPilotDdaMagOpenSearchPanel
```

---

## 9. Questions for a reviewer

1. Is the three-way split of the `-1` population right, and what should the
   continuum's treatment be? (§6a, §6b)
2. Are the KDE bandwidth/persistence and empirical-valley rules stable across
   instruments and acquisition methods, or should they be calibrated per run?
3. Does the pooled-validation argument hold given per-sample databases? (§6d)
4. Is deriving `minPeakArea` from the binomial floor sound, or does it bake the
   test's assumptions into the data it will be applied to?
5. Is the spectrum-level collapse placed correctly — after scoring, before all
   FDR — and should the classifier's *training* seed also be collapsed? (it
   currently is)

# Plan: identification, closed search first and open search on the remaining spectra

Status: built (switch `closed_stage_enabled`, on by default), not yet validated against section 6.
Differences from sections 3-5 as built: the D1 classifier uses search features only (no MSBooster
features for closed PSMs); the D2 walks over the remaining spectra are pooled per class rather
than per peak; the per-peak D1 agree and contradicted counts and the allow-list exemption go to
`combined/closed_open_summary.tsv` rather than `substitution_peaks.tsv` and
`modification_verdicts.tsv`; the peptide-level check filters `validated_psms.tsv` (so FlashLFQ and
`final_peptides.tsv` follow it) and its run-info is the same summary file. With the switch off the
open-only path runs unchanged, folding included. Version 8 of what was section 8 of
`genome-catalogue-resolution-plan.md`, now a document of its own; that plan keeps the genome test
(its sections 1-7), and "step 1" to "step 5" and E1-E3 below refer to it. Version 8 applies the
review of version 7 (all in section 3): the variant outputs get the unexplained-shift PSMs through
a separate diagnostic set once that tier leaves the identifications, no mask of closed-contradicted
candidates in peak detection, allow-list masses exempt from the substitution withdrawal instead of
a closed-agreement overrule, a peptide-level check of the union, and the reading of the genome
chi-square with the closed class.
After the build, E1 and E2 of the genome plan are run again on its output (section 6).

Comes out of E2 extended to the whole result (PSMs, peptides, proteins,
genomes), measured against the closed searches on 12mix. Scripts in the agent store
`...\files\mix12_eval\`: `50_open_vs_closed.py`, `51_spectrum_level_closed_vs_open.py`,
`52_cascade_simulation.py`, `53_weak_match_value.py`, `54_peak_picking_closed_informed.py`,
`55_contradicted_negatives.py`; tables in `out/c_*.tsv`. PXD050781: `...\files\pxd050781_eval\39_replication_ground_truth.py`.

## 1. Criterion

Once statistics that do not hold are removed, the open search is good enough when its accuracy is
comparable with the closed search and its sensitivity is higher. No further refinement after that.

Measured the same way for every workflow on 12mix:
- Sensitivity: target PSMs, distinct base peptides (I = L), proteins at q <= 0.01, genomes accepted.
- Accuracy: share of distinct peptides and of PSMs whose base peptide is not a substring of any
  spiked genome's proteome (I = L); proteins at q <= 0.01 and accepted genomes on non-spiked
  genomes, by relation to the spiked genomes (same genus, same family, other).

Peptides that the closed and the open search both identify are 3.0% outside the spiked proteomes:
real peptides absent from the UHGG representatives (strain differences). The shares are therefore
compared between workflows, not read as false discovery proportions. They are also a lower bound on
errors: most of each run's refined database is spiked proteome, so a wrong peptide is usually still
"spiked". Where the open and closed searches name unrelated peptides for the same precursor (finding
5), both sides are about 3.5% outside, although at most one of them is right.

## 2. Evidence

Workflows as on the server (closed searches of 2026-10-11, open searches of 2026-10-12):

| | FragPipe closed | pFind | Open, MSFragger | Open, MSFragger + Sage |
|---|---|---|---|---|
| Target PSMs | 564,712 | 896,063 | 850,230 | 884,280 |
| Distinct peptides | 93,978 | 104,280 | 107,644 | 107,971 |
| Peptides not in spiked | 4.61% | 4.76% | 6.05% | 5.47% |
| PSMs not in spiked | 3.35% | 1.87% | 2.70% | 2.51% |
| Proteins q <= 0.01: spiked / non-spiked | 10,110 / 3,038 | 9,612 / 2,522 | 9,805 / 1,541 | 9,578 / 1,504 |
| Genomes: spiked / same genus / same family / other | 12 / 11 / 15 / 5 | 12 / 8 / 13 / 2 | 12 / 8 / 4 / 2 | 12 / 9 / 4 / 2 |

(FragPipe PSMs from `psm.tsv`; peptides and shares from `final_peptides.tsv`. The closed workflows
still use the resampling genome score, so part of the genome and protein difference is the genome
test, not the search mode.)

1. **The unexplained-shift tier does not hold its FDR.** `validated_psms.tsv` claims 1% PSM FDR in
   every tier (decoys / targets 0.99-1.00%), but the share outside the spiked proteomes is 2.16% of
   PSMs and 4.70% of peptides for `PEAK_VALIDATED`, against 5.62% / 8.03% (`NO_PEAK_GROUP`) and
   4.40% / 7.28% (`REJECTED_PEAK_GROUP`) for `IDENTIFICATION_ONLY`. Open-search peptides with an
   unmodified form are 4.66% outside (closed 4.61%); peptides seen only with a shift 11.1%.
   Probable mechanism: an unexplained shift absorbs a sequence difference, so a relative's peptide
   plus a compensating shift beats the true peptide; a reversed decoy cannot model such a
   near-miss. Decision: unexplained-shift PSMs are not identifications (section 3).
2. **The open search loses closed answers in its ranking.** Of 564,712 closed-search target PSMs,
   the open-search spectrum winner (same file and scan; retention time and charge agree) has the
   same peptide in 367,342. Of the 197,272 disagreements: 93,980 are a different unmodified,
   validated peptide (98% of them in a spiked proteome; 91,936 of them have a precursor mass more
   than 0.02 Da from the closed search's, i.e. a co-isolated precursor where both answers may be
   right, and only 2,044 the same precursor); 25,376 an open-search decoy; 63,874 an unexplained shift (33,044
   identification-only, 30,830 not validated). In these groups 95-97% of the closed peptides are in
   a spiked proteome. 99.7% of the closed PSMs' peptides are in the open search's refined per-run
   database, so the loss is the ranking over a +-500 Da window, not the database.
3. **The open search finds what the closed search misses.** Open-validated spectra the closed
   search did not accept: 185,280 unmodified (97.6% in spiked) and 32,889 known-shift (97.4%).
   "Known shift" here means an accepted peak at a mass that is not 0-2 isotopes, which is not what it
   should mean: 12,689 are in accepted known-shift peaks (1.9% outside the spiked proteomes), 17,883
   were folded onto the zero peak by the fragment-evidence gate (finding 5), and 2,317 are zero-peak
   isotope members at offset -1 or +3 (8.6% outside).
4. **A cascade keeps both** (simulation, `52_cascade_simulation.py`). Stage 1: the FragPipe closed
   PSMs. Stage 2: open-search spectrum winners on the remaining spectra, q recomputed per class
   within those spectra, unmodified and known shift at q <= 0.01:

| | Target PSMs | Distinct peptides | Peptides not in spiked | PSMs not in spiked |
|---|---|---|---|---|
| Closed search | 564,712 | 92,148 | 4.26% | 2.50% |
| Open search, peak-validated tier | 686,036 | 86,067 | 4.70% | 2.16% |
| Open search, all validated tiers | 850,230 | 107,644 | 6.05% | 2.70% |
| Cascade (unmodified, known shift) | 776,371 | 99,597 | 4.61% | 2.35% |
| Cascade plus unexplained shift (reference) | 900,732 | 117,651 | 5.50% | 2.70% |

   The cascade gives +37% PSMs and +8% peptides at a PSM share below the closed search's and a
   peptide share 0.35 points above it. The 7,449 peptides it adds are 8.9% outside the spiked
   proteomes; the closed search's own peptides that the open search lacks are 10.2%, so the margin
   of each search is about equally clean. Stage 2 accepted 179,543 unmodified (1,795 decoys) and
   32,116 known-shift (321 decoys) spectra. The known-shift class includes the gate-folded PSMs of
   finding 3; without them the cascade has up to about 18,000 fewer PSMs, still about 1.34x the
   closed search.

Findings 5-7 compare each open rank-1 candidate with the closed search's accepted PSM on the same
spectrum (`54_peak_picking_closed_informed.py`). A candidate is *contradicted* when it is a decoy or
a different base peptide on the same precursor (observed masses within 0.02 Da): 102,173 of
1,253,437 rank-1 candidates (76,492 targets, 25,681 decoys); 367,292 agree; 95,149 differ on another
precursor.

5. **The fragment-evidence gate folds wrong peptides onto the zero peak.** `DeltaMassEvidenceGate`
   treats a peak whose members gain nothing from placing the shift on fragments as unmodified
   peptides with a wrong precursor mass, and moves them to the zero family. For the +3.5 to +5.6 Da
   bands: 18,250 spectrum winners folded; the closed search accepted 7,196 of these spectra and
   agrees with none; 7,134 name a different peptide on the same precursor. The two peptides are
   unrelated (median 5% of residues shared at the ends; 31% the same length), and 98% of the open
   candidates score worse with the shift placed. Example: `HTFEYLR` (+4.05 Da) against the closed
   `LVPEGVEAR`. No gain means the shift is not on that peptide, which a wrong peptide explains as
   well as a precursor error does, and where the closed search can test it, it is the wrong peptide.
   Folding turns them into validated unmodified PSMs.
6. **Oxidation was withdrawn as a substitution, but the closed search confirms it.** The
   substitution analysis withdrew the oxidation (+15.9948) and carbamidomethyl (+57.0214) peaks.
   84% of oxidation's rank-1 candidates have a closed answer (the closed search has variable Met
   oxidation), and 93% of those are the same peptide. Other peaks with a closed-confirmed majority:
   +17.027 (70%), +42.010 (55%, acetyl), +16.031 (53%). The substitution-shaped peaks the analysis
   withdrew (+-14.016, +30.011, -58.005, ...) agree in 0-6% and are contradicted in 93-100%, with the
   closed peptides cleaner (7.8% outside the spiked proteomes, open 18.1%): those withdrawals are
   right.
7. **Removing contradicted candidates does not change the peak set that matters.** Peak discovery
   (`DeltaMassKdePeakFinder`, ported; it recovers 1,996 of the workflow's 2,000 modes) on all rank-1
   candidates minus the contradicted ones: modes after the separation gate 5,687 -> 4,222. The 684
   modes that disappear are junk (96% of their closed-answered members contradicted; 14.6% of
   their peptides on unanswered spectra outside the spiked proteomes, 5.9% for the modes kept), but
   `kdeMaxModes = 2000` binds in both cases, so 684 equally poor modes (15.3%) take their place. The
   known-shift screen gives the same verdict for all eight known-shift peaks (decoy fractions move
   by under 0.6 points). Discovering on unanswered spectra only (the cascade's open-stage
   population) also loses 696 modes, but half of their answered members agree with the closed
   search: it removes support from real peaks (oxidation, deamidation: 84% and 90% of their
   candidates have a closed answer).

   Where a peak is a modification the closed search cannot express, a contradiction says nothing
   about the open answer: in accepted known-shift peaks (pyro-Glu, Na, Ca, deamidation) the
   contradicted open peptides are 5.5% outside the spiked proteomes and the closed ones 7.7%. In
   substitution-shaped peaks (5,791 contradicted targets), other unexplained peaks (23,402) and out
   of peak (32,872) the closed answer is the better one (3.5-7.8% outside against 13-18%).
8. **PXD050781: replication across 21 replicate runs** (`pxd050781_eval/39_replication_ground_truth.py`;
   single strain, MSFragger-only open search against FragPipe closed and pFind). Replication is used
   as extra ground truth here, not as a method.
   - *Random matches replicate too.* Of the open search's accepted decoy PSMs, 70.6% carry a decoy
     peptide accepted again in another run (pFind 96%: 11,500 PSMs on 1,086 decoy peptides); 80% of
     recurring decoy peptides sit on the same precursor mass (within 0.02 Da) in every run. A wrong
     match is a function of the spectrum, and replicates reproduce the spectrum. Replication
     therefore marks part of the random matches only: a PSM whose peptide is accepted in no other
     run is suspect, a replicated one is not thereby correct. Per peptide it separates better
     (decoy peptides: median 1 run, 90th percentile 4; target peptides: median 8).
   - *By class* (share of PSMs whose peptide the open search accepts in no other run; share never
     identified unmodified by any workflow in another run): unmodified on closed-accepted spectra
     0.46% / 0.25%; unmodified on spectra the closed search did not accept (107,561 PSMs) 1.81% /
     1.56%, by classifier-score quintile 7.5, 0.8, 0.4, 0.2, 0.1% for the first measure; known shift
     0.7-1.4% / 10-21%; unexplained shift (identification-only, 152,792) 2.8-3.4% / 34-38%. The
     unexplained tier is the weakest by both measures, as on 12mix. The open search's extra
     unmodified PSMs are as replicated as the reference except in the lowest score fifth.
   - *Same-precursor contradictions with the closed search*, which peptide is identified unmodified
     in another run: unexplained shift (40,816 pairs) only the closed one 13,325, only the open one
     2,071; known shift (3,704) 874 against 427; gate-folded (3,251) 190 against 295; unmodified
     (3,643) 149 against 437. Unlike 12mix, the closed answer also wins against accepted known
     shifts here.
   - *Fragment-evidence gate.* 16,592 of the 25,972 folded PSMs agree with the closed search, all at
     +0.5 to +2.5 Da: precursor isotope errors the closed search accepts too (`isotope_error =
     0/1/2`). In every other band none agrees and all 1,570 with a closed answer are contradicted on
     the same precursor (12mix: the +3.5 to +5.6 Da bands, none of 7,196). Folding is right for
     isotope-like deltas and wrong elsewhere; with the closed stage first, the isotope cases are
     accepted by D1 (section 3).
   - *Cascade* (closed, then open unmodified and known shift on the remaining spectra, folding
     stopped; q not recomputed within the remaining spectra): 775,094 PSMs (1.21x the closed
     search's 638,886), 36,040 peptides (+13%); PSMs whose peptide is accepted in no other run 0.76%
     against 0.62%.
9. **Closed-contradicted PSMs as classifier negatives, offline** (`mix12_eval/55_contradicted_negatives.py`,
   both datasets, MSFragger-only). The workflow's classifier rebuilt from `scored_psms.tsv` (logistic
   regression on standardised features, semi-supervised, three folds by run, fold undone), trained on
   decoys only and on decoys plus the unexplained-shift targets the closed search contradicts on the
   same precursor (12mix 73,204, PXD050781 88,583). Judged on spectra the closed search did not
   accept, q per class:

| Unexplained shift, q <= 0.01 | 12mix PSMs | not in spiked | PXD050781 PSMs | no exact-mass ID in another run |
|---|---|---|---|---|
| Decoys only | 157,820 | 4.20% | 115,511 | 33.7% |
| Decoys and contradicted | 115,860 | 2.79% | 88,981 | 26.5% |
| Open unmodified (reference) | about 180,000 | 1.9% | about 103,000 | 0.7-0.9% |
| Open known shift (reference) | about 13,700 | 1.7-1.8% | about 28,300 | 9.1% |

   At equal numbers the extra negatives halve the error share at the top (12mix top 25,000: 1.78% ->
   0.89%; PXD050781: 19.5% -> 12.4%), the curves meet at about 125,000 (12mix) and 100,000
   (PXD050781), and deeper the new score is worse: near-misses are moved down, not removed. Under
   the decoys-only score not one decoy ranks among the top 75,000 unexplained-shift targets on either
   dataset, so the decoys do not see these errors; the extra negatives make the 1% cut stricter only
   because decoys then mix with near-miss-like targets. On closed-answered spectra of held-out runs
   the closed search still contradicts 47% (12mix) and 63% (PXD050781) of the accepted
   unexplained-shift PSMs (decoys only: 51% and 67%). Unmodified and known-shift classes are
   unchanged.

## 3. Design

- **Stage C, closed search.** A closed MSFragger search per sample on the same refined database as
  the open search, settings as the FragPipe closed workflow: precursor +-20 ppm, `isotope_error =
  0/1/2`, strict trypsin, 2 missed cleavages, length 7-50, Met clipping, fixed carbamidomethyl C,
  variable M oxidation and protein N-terminal acetyl, `output_report_topN = 1`, mass calibration on.
  Its own folder and pepXML beside `sample_db`; staleness by parameters text, as the open search.
- **Stage D1, closed validation.** Closed PSMs scored by the PSM classifier trained on closed PSMs
  only (search scores, mass error, charge, length, missed cleavages; MSBooster features when the
  workflow runs MSBooster); target-decoy q per spectrum; accepted at q <= 0.01. Class: unmodified.
- **Stage D2, open validation.** Peak fitting, fragment-evidence gate, substitution analysis and the
  peak screen as now, score-independent, on all rank-1 open candidates. Candidates D1 contradicts
  are not masked out of peak detection: finding 7 shows the mask changes no accepted peak while
  `kdeMaxModes` binds, and without it peak detection does not depend on the closed stage. The
  classifier is trained as now. After the spectrum-level collapse, q-values are
  computed only over spectra D1 did not accept, per class: unmodified, known shift (accepted peak).
  Unexplained-shift winners are not accepted as identifications; they stay in `scored_psms.tsv` and
  `unexplained_psms.tsv` with their q for diagnosis.
- **Fragment-evidence gate: reject, do not fold.** A peak without shifted-fragment support is still
  rejected as a modification, but its PSMs keep their observed delta and peak and are unexplained
  shifts (not identifications), instead of being moved to the zero family (findings 5 and 8). Where
  the peptide really is unmodified with a wrong precursor mass, D1 finds it at the right mass: on
  PXD050781 that is the +0.5 to +2.5 Da folds, and only the 1,223 of them without a closed answer
  are lost, against 4,865 removed from the other bands. Allow-list masses stay exempt, as now.
- **Allow-list masses are exempt from the substitution withdrawal.** `SubstitutionReassigner` does
  not withdraw a peak within tolerance of an allow-list mass (`ExpectedPeakSelector.DEFAULT_SHIFTS`),
  the same exemption the fragment-evidence gate already applies; the peak then goes through the
  known-shift screen like any allow-list peak. On 12mix this keeps oxidation (+15.9949) and
  carbamidomethyl (+57.0215) in the screen (finding 6). For every substitution-shaped peak the
  rank-1 targets with a D1 answer on the same precursor are counted as agreeing (same base peptide)
  or contradicted, and both counts are reported; they decide nothing, because they exist only for
  modifications D1 can express (finding 7) and only when the closed stage runs.
- **Variant outputs keep their input.** `callSequenceVariants`, `writeSubstitutionForms` and the
  pangenome variant check (`OpenSearchVariantCandidateSource.collect` in `runPangenomeVariantCheck`)
  read `validation.validated` today, which holds the identification-only tier, where the PSMs of
  withdrawn substitution peaks and of rejected peaks land. Once that tier leaves `validated`, these
  three receive accepted plus unexplained-shift PSMs as a separate diagnostic set; otherwise
  variant calling would lose its candidates without notice. Variant outputs are diagnostics and are
  not identifications.
- **Contradicted candidates as known negatives: tested offline, not adopted** (finding 9). They
  improve the ranking of unexplained-shift PSMs near the top but not the error estimate, so they
  cannot make that class an identification class. Unexplained shifts stay out of the
  identifications.
- **One identification per spectrum.** A spectrum accepted in D1 is not re-identified in D2; the
  number of D1 spectra whose open winner differs is reported.
- **Peptide-level FDR.** Per stage and class, unit base peptide (I = L) with its best score, at
  q <= 0.01, for `final_peptides.tsv` (FragPipe applies a 1% peptide filter too). Same units and
  classes as step 1 of the genome test.
- **FDR statement.** Each stage controls its spectrum-level FDR at 1% within its own spectra, so the
  union is at most 1% (Kertesz-Farkas, Keich and Noble, J Proteome Res 2015, cascaded search).
  The argument holds for spectra, which belong to one stage; peptides do not. A true peptide is
  often accepted in several stages and classes and counts once in the union, while decoy peptides
  rarely repeat, so the union of per-class peptide lists can exceed 1%.
- **Peptide-union check.** After the per-class peptide walks, the union of accepted peptides (each
  base peptide once, target or decoy) gives decoys / targets. Above 1%, every class's peptide
  threshold is lowered by the same factor until the union is at or below 1%; the union FDR, the
  factor and the per-class counts before and after go to the run-info.
- **Genome test, step 1.** Classes become closed (D1), open unmodified, open known shift and
  unexplained (never evidence); each walks its own score. `pi_g` keeps one score scale: the open
  search's scored decoys of any class, as now; F counts unique decoy evidence peptides of every
  evidence class, the closed class included.
- **Reading the genome chi-square.** `F` then includes closed-stage decoys while `pi_g` comes from
  open-search decoys. If the two searches place decoys differently across genomes (the closed
  search sees only candidates within +-20 ppm), the chi-square of accepted decoys against
  `F * pi_g` can fail because of the two sources rather than a wrong `pi_g`. The run-info gives the
  chi-square for all evidence decoys and for the open-stage decoys alone; only a failure of the
  second points at `pi_g`.

## 4. Outputs

- `validated_psms.tsv`: closed and open accepted PSMs only, with a `stage` column (`closed`,
  `open`); `validation_tier` loses `IDENTIFICATION_ONLY`.
- `scored_psms.tsv`: unchanged rows plus `stage` and `closed_accepted`.
- Run-info of the validation: accepted spectra and decoys per stage and class, the D1/D2
  disagreement count, and the peptide-union check (union decoys / targets, the tightening factor,
  per-class peptides before and after).
- `substitution_peaks.tsv` and `modification_verdicts.tsv`: per substitution-shaped peak the D1
  agree and contradicted counts, and whether the allow-list exempted it.
- `genome_scores_runinfo.txt`: the decoy chi-square for all evidence decoys and for open-stage
  decoys alone.
- Variant outputs (`sequence_variants.*`, `substitution_forms.tsv`, `variant_*.tsv`): from accepted
  plus unexplained-shift PSMs, as now.
- FlashLFQ, `final_peptides.tsv`, protein tables: from `validated_psms.tsv` as now.

## 5. Code touch points (expected)

- `MagOpenSearchHeadlessRunner`: step C adds the closed search (`runSampleDbClosedSearches`); step D
  reads its pepXML, runs D1, passes the accepted spectra to `validate`; `writeValidatedPsms` drops
  the identification-only tier from the validated output.
- `MSFraggerHeadlessRunner`: closed-search settings next to `OpenSearchSettings`.
- `OpenSearchLayout`: closed-search folder, pepXML and params.
- `ModificationValidator.validate`: acceptance restricted to spectra outside D1; unexplained
  identifications no longer in `validated`.
- `MagOpenSearchHeadlessRunner`, variant callers: `callSequenceVariants`, `writeSubstitutionForms`
  and `runPangenomeVariantCheck` (`OpenSearchVariantCandidateSource.collect`) take the diagnostic
  set (accepted plus unexplained-shift PSMs) instead of `validation.validated`.
- `DeltaMassEvidenceGate.apply`: rejected peaks no longer reassigned to the zero family
  (`REASSIGNMENT_TAG` retired); `Result.psmsReassigned` becomes the count left as unexplained.
- `SubstitutionReassigner` (or `Gate.withoutPeaks`): allow-list masses exempt from withdrawal;
  per-peak D1 agree / contradicted counts, reported in `substitution_peaks.tsv` and
  `modification_verdicts.tsv`.
- Peptide-level walks (`final_peptides.tsv` writer): the peptide-union check and its run-info.
- `GenomeEvidenceTest` calibration: chi-square for all evidence decoys and for open-stage decoys.
- `SemiSupervisedPsmClassifier` / `PerEngineClassifier`: a closed-only training set.
- `OpenSearchGenomeEvidence.readHits`, `GenomeEvidenceTest.PeptideClass`: the closed class.
- `DDAMagOpenSearchParameters` and its IO: switch for the closed stage (default on).
- Tests: cascade acceptance (a D1 spectrum never re-identified; D2 q over remaining spectra only),
  identification-only PSMs absent from `validated_psms.tsv`, closed class in step 1; variant
  outputs unchanged when only the identification tier changes (same diagnostic set in, same
  `sequence_variants.tsv`, `substitution_forms.tsv` and variant candidates out); a no-gain peak's
  PSMs keep their delta; an allow-list peak of substitution shape not withdrawn, a
  non-allow-list one withdrawn as before; the peptide-union check tightening all classes by one
  factor when the union exceeds 1%.

## 6. Validation and stop rule

- The PSM gain depends on the dataset (cascade simulations, folding stopped: 12mix about 1.34x,
  PXD050781 1.21x), so the target is per dataset: target PSMs above the closed search by at least
  90% of the simulated gain (12mix at least +30%, PXD050781 at least +19%), distinct peptides at
  least equal to the closed search.
- 12mix MSFragger-only, against FragPipe closed: peptide share outside the spiked proteomes at most
  0.5 points above, PSM share not above; 12 of 12 spiked genomes and 0 reversed genomes accepted;
  protein and genome rows of the table in section 2.
- PXD050781 MSFragger-only: share of PSMs whose peptide is accepted in no other run at most the
  simulation's 0.76% (closed 0.62%); in same-precursor contradictions, the D2 winners D1 overrules
  are the side less often identified unmodified in other runs (finding 8); target genome accepted
  first.
- 12mix MSFragger + Sage: same measures as 12mix MSFragger-only.
- The spiked-proteome shares are a lower bound (section 1), so two direct checks go with them: on
  D1-accepted spectra, the share of D2's would-be winners that D1 contradicts on the same precursor,
  per class (it should fall for the zero family once folding stops), and the oxidation peak
  accepted.
- Runtime added by the closed search, per run.
- E1 and E2 of the genome plan are run again on the build's output (12mix, both open workflows),
  since the evidence classes change: species, genus and reversed-genome acceptance, the chi-square
  (both readings, section 3), and the non-spiked proteins by relation.
- When the 12mix criteria hold, identification work stops; remaining effort goes to step 4.

## 7. Risks

- The simulation's stage 1 is FragPipe's (MSBooster and Percolator). A closed stage scored by the
  MetaPilot classifier may accept fewer spectra; the open stage then takes most of them back as
  unmodified, so the loss is bounded, but the 12mix numbers decide.
- Co-isolated precursors: one identification per spectrum keeps either the closed or the open
  peptide when both may be right (part of the 93,980 unmodified disagreements). Searching chimeric
  spectra is a separate change.
- The peptide share of the cascade is 0.35 points above the closed search before peptide-level
  FDR; the peptide filter should close part of that, the 12mix run shows how much.
- Stopping the fold removes up to about 18,000 open PSMs on 12mix and 6,088 on PXD050781 on spectra
  D1 does not answer. Outside the isotope-like +0.5 to +2.5 Da bands the closed search never agrees
  with a folded PSM on either dataset; inside them (PXD050781) it agrees, and 1,223 unanswered ones
  there are lost. The simulations in section 6 already exclude folded PSMs.
- Replication across runs does not see systematic errors (finding 8), and the spiked share is a
  lower bound (section 1); no single accuracy measure is complete, hence the several in section 6.
- `kdeMaxModes = 2000` decides the peak set on 12mix (7,416 modes pass support, persistence and
  file recurrence; 5,687 after separation). Not changed: modes outside accepted peaks make no
  claim, and their PSMs are not identifications under section 3.
- Zero-family isotope members at offset -1 or +3 are 8.6% outside the spiked proteomes (2,317
  open-only PSMs); left as they are unless section 6 fails.

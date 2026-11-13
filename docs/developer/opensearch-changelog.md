# Open search changelog (unreleased)

Changes to the MetaPilot open search workflow (`DDA_MAG_OPENSEARCH`), held back from the public
`CHANGELOG.md` while the workflow is hidden from the GUI (`ReleaseFeatures.OPEN_SEARCH_WORKFLOW`).
Move them back into `CHANGELOG.md` when the workflow is released. Entries shared with public
workflows (pangenome variant check, resumed runs, parquet export, enzyme names) stay in
`CHANGELOG.md` without their open-search parts; the full text before the split is in tag
`archive/2026-10-14-opensearch`.

### Changed

- **One substitution call (2026-11-10; entrapment run pending, not yet adopted).** A shifted
  spectrum is a modification, a substitution or undecided by one decision, the identification's
  own (`ShiftLocaliser`), and every output reads it. Counts are from the round-5 pseudo-reversed
  entrapment runs (`search_results/combined_peaktest5/`, Qiao faecal and 12mix, 9 runs each).
  - An accepted-shift winner explained on a residue a substitution fits moves to
    `OPEN_SUBSTITUTION_SITE` (reported, not identified, no genome evidence), as off-peak winners
    already did. `PeakSiteShare` withdraws a whole peak only at a 0.6 share of such sites; below
    it 508 (faecal) and 4,336 (12mix) PSMs stayed identified, +40.031 above all. These PSMs still
    count in their peak's fragment test.
  - UNIMOD's amino acid substitution entries (`Ala->Thr`, `Val->Xle`) are no modification
    explanation (`UnimodCatalog.Entry#isAminoAcidSubstitutionOnly`). They were the only UNIMOD
    explanation of 80% (faecal) and 83% (12mix) of the substitution-site PSMs, so a substitution
    alone read as "substitution or modification".
  - New column `shift_class` after `shift_explanation` in `scored_psms.tsv`,
    `validated_psms.tsv`, `peak_validated_psms.tsv` and `unexplained_psms.tsv`: `modification`,
    `substitution` or `modification_or_substitution`; empty unless the shift is explained.
  - With the joint identification, a PSM carries a substitution verdict exactly when it is
    `OPEN_SUBSTITUTION_SITE` and a substitution alone explains its site
    (`SubstitutionReassigner#alignWithIdentification`). `substitution_forms.tsv`,
    `sequence_variants.tsv` and the side quantification follow. Before, the verdict was the
    peak-level one made before scoring: 639 of 6,400 (faecal) and 1,289 of 5,169 (12mix) PSMs
    with a verdict were also in the main result, and 489 and 495 set-aside substitution-site PSMs
    had none. `substitution_peaks.tsv` and `substitution_candidates.fasta` are still the
    peak-level analysis.
  - An accepted-shift PSM whose fragments do not place the shift is no longer exported as the
    unmodified peptide (FlashLFQ was then given a mass its precursor did not have). It takes the
    UNIMOD modification of its peak's mass on a residue or terminus of the peptide where UNIMOD
    lists one, a common specificity before a rare one (`ShiftLocaliser#assign`, new
    `shift_state = assigned`, no site probability), and is exported and quantified as that
    modified peptide, as a closed search places a variable modification. Without such a site it
    moves to `OPEN_UNEXPLAINED`. So does a PSM whose fragments place the shift on a residue nothing
    explains. No mass range is special: the +1.95-1.98, +4.0-4.1 and +4.9-5.0 Da bands have no
    UNIMOD entry and leave the main result. Counted on round 5: of 6,219 (faecal) and 51,251
    (12mix) unplaced PSMs, 425 and 8,628 have a UNIMOD site (0 and 8 entrapment PSMs); with the
    1,627 and 4,685 contradicted ones, 7,421 (1.4%) and 47,308 (5.5%) PSMs leave the main result,
    holding 274 of 64,986 and 662 of 90,996 peptides.
  - The peak fragment test now runs before these per-PSM decisions, on each peak's winners
    accepted with the whole peak as one group, so the decisions cannot change which peaks pass;
    the q-values that are reported are walked afterwards, on the final groups.
  - `sequence_variants.tsv`: "reference also observed" is read from the whole identified result
    (identified PSMs exported as the unmodified peptide, closed search included), not from the
    open rows alone. On round 5, 45 of 483 (faecal) and 44 of 286 (12mix) `SEQUENCE_VARIANT`
    calls had their reference identified unmodified, 41 and 27 of them by the closed search only.
  - Every level now shows the decision (`open-search-result-files.md` lists the files):
    `final_peptides.tsv` gains `Shift state` and `Shift class` after `Shift explanation`;
    `combined/set_aside_peptides.tsv` and `combined/set_aside_genomes.tsv` (`SetAsideReport`) give
    the set-aside PSMs per peptide and reason, and per genome beside its identified PSMs.
  - The catalogue (UHGP-100) lookup takes a peak-level proposal only where it is also the
    identification's call for that PSM (`OPEN_SUBSTITUTION_SITE`, the fragments' site, a residue
    its explanation names); other identifications keep every proposal.
  - Not changed: the off-peak site-probability threshold, and the 0.005 Da tolerance of an
    explanation on a single PSM's own mass.

- **Variant outputs in `sequence_variants/`, with their own quantification (2026-11-09).** Every
  variant file now sits in `sequence_variants/`; `substitution_peaks.tsv` and
  `substitution_candidates.fasta` moved there too, and an earlier run's copy in `combined/` is
  moved on first use. Substitution PSMs that are not identifications are quantified by a separate
  FlashLFQ run (`SubstitutionQuantification`), so no main table changes:
  - `substitution_psms.tsv`: the target PSMs with a substitution verdict that are not in
    `validated_psms.tsv` (`quant_role=variant`, as base peptide plus delta, protein `UNASSIGNED`,
    with `variant_id`), and the validated unmodified target PSMs of the same base peptides
    (`quant_role=reference`), so a variant and its reference are quantified in one run.
  - No match-between-runs, since these PSMs are not identifications. FlashLFQ 2.0.0.308 cannot
    switch it off (`--mbr` is a switch that defaults to true), so its `MBR` values are set to 0
    and marked `MBR_removed`; MBR fills only files without an MS/MS detection, so the MS/MS
    values are those of a run without it. Normalisation off: FlashLFQ would fit per-file factors
    on a few hundred peptides, unstable and on a different scale from the main run's.
  - Fixed with it, for every workflow: `FlashLfqTask` passed `--mbr false` and `--nor false`,
    which FlashLFQ reads as switches turned on, so a run asked not to normalise (the pFind
    label-free option, this side run) normalised. The switches are now passed only when on. The
    open-search and pFind MAG runs ask for both, so their results do not change.
  - `substitution_quantified_peptides.tsv`: FlashLFQ's peptide table plus `quant_role`,
    `variant_ids`, `reference_sequence`, `Ratio_<sample>` (variant over reference, empty unless
    both are above zero) and `shares_main_feature`: samples where the variant's MS1 peak is also
    a main-run peak of another peptide (same file and charge, m/z within 10 ppm, FlashLFQ's own
    tolerance, apex inside the main peak). Two runs cannot resolve a shared feature, so it is
    flagged. `substitution_quantified_peaks.tsv` keeps the peak table; the protein table, one
    `UNASSIGNED` protein, is not kept.
  - FlashLFQ rewrites `ExperimentalDesign.tsv` in the spectrum folder; the side run restores it,
    or removes it when there was none, and leaves every other file alone (another job may be
    writing there).

- **`OPEN_KNOWN_SHIFT` is now `OPEN_ACCEPTED_SHIFT` (2026-11-09).** "Known shift" meant a mass on
  a predefined list; the group holds open-search PSMs on a delta-mass peak accepted by the
  evidence rule, and no list is left. Renamed with it: the cascade class (`open_accepted_shift`),
  the genome-test peptide class `ACCEPTED_SHIFT`, the classifier feature `open_accepted_shift`
  (features are read by position, so the model is unchanged), and
  `ShiftLocaliser.checkAllowListShift`, now `checkAcceptedShift`. Its "allowed residues" are
  called explained residues: residues where a UNIMOD entry explains the shift. The
  `known_shift` column of `modification_verdicts.tsv`, empty since peaks are accepted on
  evidence, is removed. Only labels change; identifications do not. Results written before
  this date still say `OPEN_KNOWN_SHIFT`.

- **No hand-picked mass list is left (2026-11-08).** The seven-mass list (oxidation, deamidation,
  pyroGlu/ammonia loss, water loss, Na-H, Ca-2H, carbamidomethyl) had two remaining uses:
  - Stage 1 admitted PSMs at its masses into genome selection when they passed four checks per
    run. With the default closed marker search there is nothing to test: oxidation and N-terminal
    acetyl are variable modifications at mass difference zero, and no shifted peak forms in a ppm
    window. In all 158 closed stage-1 runs on disk nothing was admitted; the only populated window,
    0.984 Da isotope-error matches, sat at 0.77-1.04 decoys per target. The step is now skipped
    there and `expected_peaks.tsv` is not written. With the open marker search
    (`round1_search_mode = open`) the candidates are now every UNIMOD mass except isotopic labels,
    which exist only in a sample labelled on purpose; the four checks are unchanged. That mode was
    not re-measured.
  - The fragment-gain gate never folded a peak at a list mass. It now exempts only peaks inside
    the engine's `delta_mass_exclude_ranges` (-1.5 to +3.5 Da), where MSFragger scores no shifted
    fragment and no gain is no evidence, which was the reason for the exemption. UNIMOD cannot
    serve here: even without isotopic labels it matches 5 of the 73 fitted peaks in the +4 and +5
    Da artefact bands on 12mix and 5 of 52 on faecal. The gate folds only in the open-only
    identification; in the joint and cascade identifications it changes nothing.
  `unimod_compact.tsv` gains a `classifications` column, from the same `unimod.xml`. The same
  column takes isotopic labels out of the two UNIMOD lookups that explain observed shifts, since a
  label mass in a label-free sample is an isotope peak or a wrong precursor mass:
  - Shift labels: about 32 closed isotope-offset PSMs per run at 0.992 and 1.996 Da were
    explained only as `Label:15N(1)`; they stay identified and are now localised but unexplained.
    One faecal open PSM at +7.039 was identified only through `Label:2H(6)15N(1)` and no longer is.
  - The isotope-link guard: labels alone kept some +4/+5 Da artefact peaks (12mix +4.0104 at
    `Label:18O(2)`, faecal +4.0023, +4.0094, +5.0303) and faecal +44.0176 (`Acetyl:13C(2)`) out of
    isotope families; they may now link. Checked by entrapment re-run (round 5).

- **Mass shifts are accepted on evidence, not from a list (2026-11-06).** The stage-2 allow-list
  of seven masses is gone. Every non-zero peak is judged by one rule on its rank-1 PSMs of 13 or
  more residues: at least 50 of them, at most 5% decoys, and a long-peptide share of at least 0.75
  of that of all other non-zero peaks. A peak that passes is withdrawn when 60% or more of its
  localised PSMs (at least 20 of them) sit on a site where residue plus shift is an amino acid: a
  localisation-only test with no list or database check, written to `peak_site_shares.tsv`. The
  in-database test and the N->D/Q->E list no longer affect peak acceptance; they remain in the
  substitution reports. Likewise a winner localised and explained on a site where a substitution
  fits is no identification: it moves to the new group `OPEN_SUBSTITUTION_SITE`, walked and written
  to `unexplained_psms.tsv` and read by the variant outputs. `OPEN_LOCALISED` now holds only sites
  no substitution could produce; on 12mix the two parts had 4.0% and 17.7% of PSMs outside the true
  genomes, at 0.1% FDP or below. Checked again with the current localiser (2026-11-08): 11.7% in
  the withdrawn peaks (111 PSMs) and 18.8% elsewhere (1,943), against 3.5% for `OPEN_LOCALISED` and
  2.3% for exact-mass PSMs; the one-residue variant was a true peptide in 5 of 378 of the wrong ones,
  so fragment support on a substitution site confirms where the mass sits, not the peptide. With
  that localiser deamidation (0.984 Da) is withdrawn as well, at site shares 0.80 (12mix) and 0.87
  (faecal): its residues N and Q are the substitution sites of the same mass, so the rule cannot
  separate the two and makes no claim. Its exact-mass PSMs (closed-search isotope-offset matches)
  stay identified with a label that lists every explanation; of its open winners, those the
  fragments place on N or Q (10 on 12mix) and those they do not place (about 1,900) are not
  identified. Kept peaks reach a share of at most 0.50 and withdrawn ones start at 0.72, so any
  threshold from 0.51 to 0.72 gives the same result. `mass_class` in `modification_verdicts.tsv` now reads `ESTABLISHED`
  instead of `KNOWN`, and the table gains `long_peptide_psms`, `long_peptide_decoy_fraction` and
  `long_peptide_share_ratio`. Stage-1 genome selection and the fragment-evidence gate kept their
  list until 2026-11-08 (entry above). The rule decides which mass shifts are established, not what they are: names still
  come from localisation with fragment support. Evaluated offline on the entrapment runs
  (`Claude outputs/major_peak_test/`); on 12mix 69,455 PSMs at 1.21% FDP against 17,378 at 1.09%
  for the list, on faecal 11,673 at 1.23% against 4,334 at 2.03%.

### Fixed

- **Shift labels follow the fragments at peak level, for closed PSMs too (2026-11-07).** Two gaps
  in the rule that names a shift only with fragment support:
  - Each peak is now tested as a whole (`PeakFragmentTest`, written to `peak_fragment_test.tsv`):
    r1 is the share of its accepted winners whose shift the fragments carry, r0 the same share for
    a fixed sample of 2,000 accepted unmodified zero-peak winners scored at the peak's mass, taken at
    its one-sided 95% upper bound. Both rates come from accepted winners, the same population.
    Ions whose shifted position the 20 ppm tolerance cannot tell from the unshifted ion or its 13C
    peaks are no evidence either way: peak lists often hold a large fragment's 13C peak without its
    monoisotopic one, which read as +0.984 on the backbone in 8.6% of accepted unmodified 12mix
    PSMs (2.9% with the rule; closed 0.984 PSMs 81% -> 42%). Where deamidation still cannot be
    shown at this tolerance, its PSMs stay unlabelled. When r0 / r1 is above 0.10 (more than one localisation in ten
    expected by chance) the shift is not on the fragments and every winner in the peak is exported
    as the unmodified peptide. Before, a per-PSM check alone named shifts by chance: on 12mix about
    490 PSMs of the +4/+5 Da precursor-error band were localised, some "explained" as
    `Label:13C(3)15N(1)` or kynurenine. Accepted winners only, because below the FDR cut wrong
    peptides win by explaining peaks with the free shift; the bound, because a null of 0-3 in 2,000
    is not "no chance". -17.027 Da passes at 0.052-0.059 (its null is high: unmodified spectra
    carry b-NH3 and y-NH3 ions); +17.027 Da fails (48 of 1,448 PSMs supported).
  - Exact-mass winners that the closed search accepted with an isotope offset, and no variable
    modification, are checked like open shifts when their mass difference lies in a fitted peak
    other than an isotope peak of the zero family, whether or not that peak is accepted as an
    identification group: deamidation is withdrawn as substitution-shaped, yet its closed PSMs
    still need their label. On 12mix, 18,849 such PSMs sat at 0.984 Da
    (deamidation, 19 mDa from the +1 isotope step) and their fragments carried +0.984 on N/Q in 79%
    of them against 3% for unmodified peptides; they were exported as unmodified, and FlashLFQ
    searched for the unmodified mass. Their groups, scores and q-values are unchanged.
- **The host proteome is searched (2026-11-03).** `append_host_db` and `host_database` were
  ignored by the open-search workflow, so runs that set them searched the selected genomes alone.
  The host now joins each sample-specific database that stage 2 and its closed stage search, with
  one decoy per host protein made like the genomes' decoys, as the closed FragPipe workflow adds
  it to its final search. The marker search still holds the marker database only, and entrapment
  databases stay genomes-only. A run that asks for the host and cannot find it stops.
  - In the genome test, host proteins are no candidate genomes, and a peptide a host protein
    contains is no genome evidence (counted as `host` in the uniqueness step of
    `genome_scores_runinfo.txt`). Host proteins are reported as proteins with no genome, as in
    the other DDA MAG workflows.
  - On resuming, a changed host setting rebuilds the sample databases from the genomes already
    selected, without repeating genome selection; their stage-2 searches then run again.
- **Modifications are no longer folded into another peak as isotopes (2026-11-03, Unimod since
  2026-11-06).** A peak within 2 mDa of a Unimod mass is not linked as an isotope satellite when
  the isotope step misses by more than 1.5 mDa per step and by more than twice the Unimod match; it
  stays a peak of its own and is judged like any other. On PXD013427, an acetylome, acetyl
  (42.0105) had been linked as the -1 isotope of carbamyl (43.0058), 7.6 mDa off a true isotope
  step, and left the modification table. Replayed on the 1,443 isotope links of seven earlier
  runs, the guard breaks 75, none in the zero family.

### Changed

- **Defaults are the frozen workflow (2026-11-02).** A parameter file without these keys now
  runs the workflow that was validated by entrapment and used for the publication runs:
  - `identification_version` v2 (was v1); a blank value also means v2;
  - `sample_decoy_method` `pseudo_reverse` (was `reverse`);
  - `joint_ms2_max_open_rank` 3 (was 2);
  - `genome_selector_mode` Marker4 (was Marker2), in the GUI as well;
  - `joint_compare_arms` off (was on; an evaluation option);
  - `kde_max_modes` 500 (was 2000).

  Parameter files that set these keys behave as before. A resumed run whose file predates v2 and
  omits them now uses v2 with pseudo-reversed decoys, so its stage-2 searches and identification
  run again.
- **MSFragger + Sage is removed; results sit in `opensearch_MAG_search/` like
  the other DDA MAG workflows.** Closed and open MSFragger searches of the
  selected genomes replace a second engine. The GUI no longer offers Sage, and
  a saved `search_engines: "MSFragger + Sage"` is coerced to MSFragger only.
  Pooled validation, final tables, taxonomy, function and `report/` now live
  directly under `opensearch_MAG_search/` beside `search_results/` and
  `spectra_summary/`, matching `fragpipe_MAG_search/` and `pfind_MAG_search/`.
  There is no nested `msfragger_only/` or `msfragger_sage/`. An existing
  `msfragger_only/` folder is lifted to the engine root on the next run if
  `combined/validated_psms.tsv` is not already there.

### Added

- **Thermo .raw and timsTOF .d inputs.** Vendor files are searched by MSFragger
  as they are. The marker search also writes MSFragger's uncalibrated mzML
  (`write_uncalibrated_mzml = 1`), which is moved to
  `opensearch_MAG_search/spectra_mzml/<run>.mzML`. MSBooster, Crystal-C, the MS2
  features, shift localisation, the variant ranking, FlashLFQ and the spectra
  summary read that file, whose scan numbers match MSFragger's. For timsTOF
  ddaPASEF it holds one MS2 spectrum per precursor, summed over its PASEF frames,
  which msconvert does not produce. mzML inputs are read directly, as before, and
  their search parameters are unchanged. A run takes all mzML or all vendor
  inputs. A run whose converted mzML is missing repeats its marker search to
  write it again (`OpenSearchSpectra`).

- **GUI paired-shuffle entrapment switch** (MetaPilot open search). Off by
  default. Same processing as `mag-opensearch-entrapment`: after genome
  selection the sample databases also contain shuffled copies of the selected
  proteins, results go to `opensearch_entrapment_MAG_search/`, and the run stops
  after identification. A saved `additional_parameters.entrapment_enabled`
  JSON key selects the same path from `mag-opensearch`.

- **MSFragger `.pepindex` files are removed after the searches** (MetaPilot open
  search). Same job as FragPipe MAG's `cleanupMarkerSearchIntermediates`: once
  every sample's pepXML is on disk, fragment indexes under `search_results/`
  (and any fallback `marker_target_decoy.fasta` index in the engine folder) are
  deleted. Resume still skips from the pepXML. The shared `Marker4_td.fasta`
  index in the MAG database is not touched.

- **`mag-opensearch-report`** (MetaPilot open search). Rebuilds `final_peptides.tsv`,
  `razor_proteins.tsv`, `all_proteins.tsv`, `id.tsv`, `final_summary.tsv`,
  `taxonomy_analysis/` and `excluded_evidence/` from an existing
  `validated_psms.tsv` and `QuantifiedPeptides.tsv`. No search, Crystal-C,
  classification or FlashLFQ. Applies `MagResultFinalisationService` and the
  current export columns (shift explanation, spectrum count, mzML identification
  rate). Use `--out` when the parameter file still points at an HPC path.

- **Identification rate for runs started from mzML** (MetaPilot open search).
  `final_summary.tsv` reported MS/MS = MS/MS Identified (100%) when there was
  no spectra-conversion folder with `ms2Count.txt`, as on HPC. The run now
  writes `ms2Count.txt` and `rt.txt` from the mzML files to
  `opensearch_MAG_search/spectra_summary/`, never beside the inputs.
  `final_summary.tsv` then uses the real MS2 count. An existing
  `ms2Count.txt` beside the inputs (written by the Windows conversion) is
  still used as it is. Nested-layout migration used to move that folder under
  `search_results/` (via `sampleDbParams` before the summary was written),
  which hid the counts again. `spectra_summary/` is now a fixed top-level
  MAG folder; a previously nested copy is restored; the summary writer finds
  either location. Same path for FragPipe MAG CLI/HPC.

- **Spectrum-based identification rate, PSMs reported apart** (all DDA MAG
  exports).
  - `final_summary.tsv`: "MS/MS Identified" and "MS/MS Identified [%]" count
    distinct identified spectra (each scan once). A new last column, "PSMs",
    counts every identification.
  - `id.tsv` gains "Spectrum count".
  - Under identification v2, a scan can carry one identification per
    precursor, so the open search's PSMs can exceed its identified spectra.
  - Workflows with one identification per spectrum report equal counts.

- **Identification v2** (MetaPilot open search, `identification_version`,
  default `v1`; parameter file only). With `v2`, the joint closed+open
  identification follows `docs/developer/identification-v2-design.md`:
  - **Precursor queries.** An open row whose precursor Crystal-C reassigned
    (another charge, or a mass not within 0.02 Da of the closed precursor after
    whole isotope spacings) is its own query and FDR group. It never merges with
    a closed row. A scan with more than one reassigned precursor stops the run.
  - **RESET.** Decoys are split by peptide hash. Training decoys are the only
    negatives and are never reported. Estimating decoys count twice in every
    walk: PSM, peptide, genome test and protein filter.
  - **XGBoost trees** (`hist`, depth 4, 200 rounds) with monotone constraints on
    the search and predicted-spectrum scores. These replace logistic regression
    in the joint model.
  - **Protein filter** in the genome test. Proteins with two or more peptides
    are kept. One-peptide proteins are walked to 1%. Peptides whose razor
    protein is rejected are dropped.
  - **Defaults.** `joint_ms2_max_open_rank` 3 (DIA-NN predictions for open
    ranks 1-3) and `sample_decoy_method` `pseudo_reverse`, unless the file sets
    them.

  - **Shift localisation.** After the model has chosen the winners, each
    unexplained-shift winner is localised from its fragments
    (`ShiftLocaliser`): the fragment support, the best site and its
    probability.
    - Explained winners move to a new identification group, `OPEN_LOCALISED`,
      walked as its own FDR group. A winner is explained when it has support
      ≥ 2, site probability ≥ 0.75, and a substitution or UNIMOD entry at the
      site.
    - Every other unexplained-shift winner stays unexplained, and only its
      mass is claimed.
    - The model is unchanged.
    - For the genome test, `OPEN_LOCALISED` is not evidence.
    - The pangenome variant check now also takes `OPEN_LOCALISED` PSMs whose
      explanation at the site includes a substitution, whether on a peak or
      not. Each becomes a candidate at the localised site, with
      `front_end_call` = `localised` and `peak_substitution_shaped` NA.
      `variant_candidates.tsv` gains `localisation_site_probability` and
      `localisation_support`. A peak proposal naming the same PSM, site and
      residue gains these two values instead of becoming a second candidate.
    - Export: an `OPEN_LOCALISED` PSM's `Full Sequence` carries the shift on
      the localised site, with the explanation's theoretical mass. The PSMs of
      one peptide and explanation therefore form one peptide row, and the
      FlashLFQ mass is the peptide mass plus that explained mass.
      `validated_psms.tsv` gains `shift_site`, `shift_site_probability` and
      `shift_explanation`; `final_peptides.tsv` gains `Shift explanation`.
      Other workflows' peptide tables are unchanged.
    - Allow-list shifts are checked by their fragments (identification v2;
      design doc, section 8). A known-shift winner is exported in one of
      three ways:
      - with no fragment support, as the unmodified peptide (the observed
        shift stays in the mass-difference columns);
      - with fragment support and an allowed residue within a factor of 10 of
        the best site, as the modification on that residue;
      - when the fragments contradict every allowed residue, as the shift on
        its own site with no explanation.

      The identification and its FDR group are unchanged.
      `validated_psms.tsv` gains `shift_state`. On the entrapment runs this
      changed what is claimed for 26-35% (unsupported) and 18-26%
      (contradicted) of known-shift PSMs; for example, +0.984 on peptides
      without N or Q is no longer reported as deamidation.
    - See `docs/developer/shift-localisation-explanation-design.md`.

  A decoy wins a tied query. `scored_psms.tsv` gains these columns:
  `query_reassigned`, `reset_role`, `shift_state`, `shift_site` (1-based),
  `shift_site_probability`, `shift_support` and `shift_explanation`. Under
  evaluation by entrapment; v1 is unchanged.

- **Pseudo-reversed decoys for the stage-2 sample databases** (MetaPilot open
  search, `sample_decoy_method`, default `reverse`; parameter file only). With
  `pseudo_reverse`, each decoy protein is cut after every K/R and, in each
  segment, the residues between the first and the last are reversed, so every
  tryptic decoy peptide has its target's length, composition and mass. A
  palindromic segment is rotated by one instead. This applies to the normal and
  the entrapment runs. The shared marker database keeps whole-protein reversal.
  A database built with the other method is rebuilt together with its
  MSFragger `.pepindex` files, and the stage-2 searches then re-run. Under
  evaluation by entrapment; see
  `docs/developer/peptide-fdr-and-precursor-aware-design.md`.

- **Closed search first, open search on the remaining spectra** (MetaPilot
  open search, `closed_stage_enabled`, on by default; GUI: "Closed search
  first, open search on the remaining spectra"). Step C adds a closed MSFragger
  search per sample on the refined database (+-20 ppm, isotope errors 0/1/2,
  variable M oxidation and protein N-terminal acetyl) into `sample_db_closed/`.
  Step D scores the closed PSMs with a classifier trained on them alone and
  accepts spectra at 1% FDR; the open search then identifies only the spectra
  the closed stage left, as unmodified or allow-list shifts. Unexplained-shift
  PSMs are no longer identifications: they stay in `unexplained_psms.tsv`,
  `scored_psms.tsv` and the variant outputs. Peaks without shifted-fragment
  support are rejected without being folded onto the zero family, and
  allow-list masses are exempt from the substitution withdrawal. A peptide-level
  1% filter per class, tightened together when the union of classes exceeds 1%,
  decides which accepted PSMs reach `validated_psms.tsv`. New outputs:
  `combined/closed_open_summary.tsv` (spectra and decoys per stage, closed/open
  disagreements, the peptide check, per substitution peak the closed-stage
  agreement), `combined/identification_mode.txt`, and a `stage` column in
  `validated_psms.tsv` and `scored_psms.tsv` (plus `closed_accepted` in the
  latter). The genome test counts closed PSMs as their own evidence class and
  reports the decoy chi-square for open-stage decoys alone as well. Existing
  results rerun step C (the closed search) and step D. Turning the switch off
  restores the open-only identification exactly.

- **Thermo .raw input for the MetaPilot open search.** Loaded `.raw` files
  are converted to mzML with the msconvert configured under Resources before
  Step A, into the result folder's `spectra/`, where the other workflows put
  converted spectra, and every step reads the mzML. Before this, Crystal-C failed on `.raw` input, and MSBooster, Sage
  and the variant MS2 ranking could not read it. mzML input is read as loaded,
  as before. The command is `ProteoWizard`'s
  `--mzML --zlib --filter "peakPicking true 1-"`, the same command MetaPilot's
  other workflows use. A conversion is reused while it is complete (ends with
  `</indexedmzML>`) and no older than its raw file; a failed conversion's
  file is deleted. msconvert's messages go to `logging.log` with the ProteoWizard
  release. A file the installed release cannot read, such as an instrument
  newer than the release, is reported there as a version problem, and the
  run stops after converting the remaining files. Mixing `.raw` with other
  spectra in one run is refused.

- **MS2 form ranking for the pangenome variant check** (`variant_ranking`, on
  by default; MetaPilot open search only). An open-search substitution is the
  reassigner's reading of a mass shift. The same shift could also be a
  modification on any residue, or another substitution of the same mass at a
  different site. A modification at the substitution site produces exactly the
  same b/y m/z as the substitution, so only fragment intensities can separate
  them. For each candidate the ranking scores every form within its own
  spectrum: the substitution, every other same-mass substitution
  (`SubstitutionMassTable`), and the delta mass placed on each residue of the
  assigned peptide. The fragment intensities come from the DIA-NN executable
  already configured for MSBooster: its 1.8.x `--predict` interface, read the
  way MSBooster reads it. A modification form uses the assigned peptide's
  prediction, with the delta added to the fragments that contain the site. The
  score is the unweighted spectral-entropy similarity between the predicted
  intensities and the most intense observed peak within
  `variant_ranking_fragment_tolerance_ppm` (20) of each predicted fragment. The
  peaks are read from `<sample>.mzML`, either the input itself or the file
  beside it.

  Scores are only compared between forms of the same spectrum. The verdict
  compares only the forms at the candidate's own site: the substitution, the
  delta as a modification of that residue, and any other same-mass
  substitution of that residue. When the top two are within
  `variant_ranking_margin` the verdict is `AMBIGUOUS`; otherwise the top form
  wins. The default margin of 0.05 is provisional until Kai sets it. Forms at
  other sites are not in the verdict: with one placement per residue, the best
  of them would beat the single substitution form by selection alone. Their
  best score is reported as `ms2_score_elsewhere`, and
  `ms2_localises_elsewhere` is true when it beats every same-site form by more
  than the margin, meaning the delta mass seems to sit on another residue.
  `variant_candidates.tsv` gains `ms2_verdict`, the four scores,
  `ms2_localises_elsewhere` and `ms2_best_form`. DIA-NN's input and output are
  kept in `variant_ranking/` beside the two files.

  **Secondary evidence.** On PXD013427, 86% of the verdicts were `AMBIGUOUS`,
  and G->P (+40.031) got many `SUBSTITUTION` verdicts although its variants
  are no more often in UHGP-100 than its controls. The ranking is therefore a
  column only: it does not set `pangenome_class` and does not split
  `variant_pangenome_summary.tsv`, which rest on the database evidence.

  **What this does not do.** It is report-only, like the check it belongs to:
  it does not touch the FDR, the `SubstitutionReassigner` gate or
  `validated_psms.tsv`. It does not re-score pFind, whose calls are trusted
  (`NOT_APPLICABLE`). It needs DIA-NN 1.8.x, since a 2.x build has no
  `--predict`, and centroided mzML. MS-Numpress spectra and raw files without
  an mzML beside them are not read. Every candidate it cannot rank is
  `NOT_RANKED` with the reason in the log.

- **Delta-mass evidence gate** (`delta_mass_evidence_gate`, on by default). Because
  the searches set `localize_delta_mass = 1`, MSFragger scores every spectrum
  twice - with the delta mass on the fragment ions and without - and the
  difference already reaches this workflow as `delta_mass_score_gain`. A real
  modification gains from the shifted fragments; a wrong precursor mass does not.
  Peaks whose members overwhelmingly gain nothing are precursor errors, and their
  PSMs are now folded onto the zero family as the unmodified peptides they are -
  reassigned, not demoted, because the evidence says positively what they are.
  `mass_diff` is untouched; only `normalized_mass_diff` moves.

  The threshold is not tuned. Measured on 12mix (2026-09-14), the per-peak
  fraction of members with no gain, peaks over 200 PSMs, zero family excluded:
  `3 6 6 8 8 9 10 15 18 23 28 37 47 58 65 66 | 86 87 88 91 92 ... 99 100`. There
  is a twenty-point gap; gates at 0.70, 0.80 and 0.85 reject the identical 42
  peaks / 40,014 PSMs. The rejected side is every peak in the +3.9..+4.3 and
  +4.8..+5.4 Da bands (94-96% no-gain, median gain negative). The kept side is
  pyroGlu (8.3%), Ca-2H (5.9%), oxidation (18.2%), Na-H (37.2%),
  carbamidomethyl (58.4%).

- **Substitution analysis** (`substitution_analysis`, on by default). A delta mass
  can mean the SEQUENCE is wrong rather than the mass. The two accepted peaks at
  +/-14.0157 Da - CH2, i.e. G/A, S/T, V/L, V/I, N/Q, D/E - passed the per-peak
  local-FDR walk at 0.94% and 0.96% decoy while their true error against the
  twelve spiked genomes was 18.35% and 17.85%: a nineteen-fold under-estimate.
  Target-decoy is blind to it by construction, because the wrong peptide is
  itself a target sequence that no reversed or shuffled decoy ever competes with;
  and the evidence gate above passes it, because a substituted residue genuinely
  shifts part of the ion series.

  `SubstitutionMassTable` indexes every single-residue substitution by mass;
  `SubstitutionReassigner` proposes, for each PSM on a matching peak, the one
  sequence its delta AND its localization both support; `TrypticPeptideIndex`
  says whether that sequence was ever in the searched database. A peak that is
  mostly substitution-consistent is withdrawn from the modification verdict -
  not rejected on FDR grounds, but on the grounds that the question was wrong.
  New columns `delta_mass_reassignment`, `substitution_verdict` and
  `substitution_candidate` in `validated_psms.tsv`; new
  `substitution_peaks.tsv` and `substitution_candidates.fasta`.

  Mass match alone would be wrong: N->D collides exactly with deamidation
  (0.98401 vs 0.98402), A->S and F->Y with oxidation, A->Q and G->N with
  carbamidomethyl. On the same run those peaks carried 22,186 / 4,790 / 1,913
  PSMs at 2.80% / 1.68% / 1.17% entrapment - overwhelmingly the modification.
  The localization evidence is what separates them: if +15.9949 sits on an M it
  is oxidation, and A->S is never proposed.

  **What this does not do.** It does not rewrite a PSM's peptide. Choosing
  between the shifted assignment and the substituted one needs the spectrum
  scored against the candidate, and MSFragger is the scorer. That matters because
  the majority case cannot be settled by lookup at all: of the 2,369 PSMs on the
  CH2 peaks, only 25.4% propose a sequence already in the database and
  identified, 4.1% one in the database but unidentified, and **70.5% a sequence
  absent from every database this run searched** - genome selection picked one
  strain and the peptide came from a relative. `substitution_candidates.fasta` is
  the input to the bounded re-search that can settle them.


- **Crystal-C precursor refinement in the open-search workflow**, between the
  round-two MSFragger search and everything that reads a delta mass
  (`crystalc_enabled`, on by default). An open search reports
  `massdiff = precursor_neutral_mass - calc_neutral_pep_mass`; when the
  instrument picked the wrong peak of the precursor's isotope cluster, or a
  co-eluting neighbour inside the isolation window, that difference is a
  precursor-assignment error a few Da from zero, not a modification. Crystal-C
  re-reads the MS1 scan, re-picks the monoisotopic peak and rewrites the
  precursor mass. Validation, fitting and the classifier now read
  `sample_db/<sample>_c.pepXML`.

  Why it was missing and what it cost, measured on 12mix (2026-09-14): this
  workflow's open search and FragPipe's find the *same* number of rank-1
  unmodified target PSMs in the raw pepXML - 352,484 against 352,360, pooled
  over nine files - so the searches were never the difference. FragPipe runs
  Crystal-C and this workflow did not. That took FragPipe's zero-delta rank-1
  target count from 317,000 to 591,641, while this workflow's validated output
  left 316,120 of 785,180 target PSMs sitting at 0.05-5.5 Da and spent 335 of
  its 500 fitted peaks inside +/-5.5 Da on those artifacts. In `YS-194-1-1`
  alone Crystal-C moved 29,485 of 106,340 rank-1 target scans to |delta| < 0.02
  Da, keeping the peptide in 97.6% of them; 92.8% of those had a pre-correction
  delta within +/-5 Da, piled on integer values. `DeltaMassPeakFitter`'s own
  isotope-family normalisation recovers only about 45,000 of them, because just
  13% of the corrections are within 0.02 Da of a clean n x 1.00336 step - the
  rest are co-isolation artifacts that arithmetic isotope subtraction cannot
  reach.

  **It also runs on the round-one marker search**, before MSBooster and before
  genome selection reads a delta mass, because stage 1 is where the relocation
  actually pays: its admission rule is the zero family plus a fixed
  expected-shift allow-list, so a mis-picked precursor does not make a PSM
  slightly wrong there, it removes it. MSBooster is unaffected - it reads the
  `.pin`, which Crystal-C does not touch, and the rescorer joins its features
  back by `file:scan:rank`. `base_name` is preserved in the refined file, so
  per-file grouping downstream is unchanged.

  Verified on the real Crystal-C outputs for the nine stage-2 pepXMLs
  (2026-09-14), rank-1 PSMs at 1% FDR on hyperscore, purity = share whose
  protein set includes one of the 12 spiked genomes:

  | population | | accepted target PSMs | purity | peptides | pep purity |
  |---|---|---:|---:|---:|---:|
  | \|delta\| < 0.05 Da | before | 310,049 | 96.98% | 64,550 | 96.50% |
  | | after | 581,451 | 97.03% | 81,194 | 95.16% |
  | \|delta\| <= 5 Da | before | 659,561 | 96.98% | 83,861 | 94.91% |
  | | after | 675,149 | 96.98% | 85,721 | 94.91% |
  | all rank-1 | before | 577,440 | 96.35% | - | - |
  | | after | 577,440 | 96.35% | - | - |

  The last row is the honest headline: **Crystal-C creates no identifications.**
  Over the whole rank-1 set the accepted count and purity are identical to the
  digit. It changes no score and, at rank 1, essentially no peptide. What it
  does is move 271,402 correct PSMs and 16,644 peptides out of spurious
  delta-mass peaks into the zero peak at unchanged PSM purity. Not free:
  decoys inside the rank-1 zero window went from 1,192 to 6,465 (0.38% to 1.10%
  of the window), so about 1.9% of what it relocates is decoy and that window's
  internal FDR threshold tightens.

  Crystal-C is Apache-2.0, the same licence as MetaPilot, and ships inside the
  FragPipe installation Resources already points at, beside MSBooster; no new
  dependency. Cite Chang HY et al., J Proteome Res 2020. Its parameter keys,
  value syntax and defaults are taken from the shipped tool and from FragPipe's
  bundled `crystalc.params` rather than guessed, so a MetaPilot run refines the
  same way the FragPipe run above did. It deletes nothing - same
  `spectrum_query` and `search_hit` counts in and out, `alternative_protein`
  entries preserved - so it refines the candidate list rather than filtering it,
  and acceptance stays this workflow's own decision.

### Changed

- **Open-search GUI shows only study and machine settings.** Rank depth
  (candidates per spectrum), the open precursor window, round-1 mode, the
  closed-stage checkbox, engine choice, local FDR and peak-fitting knobs are
  no longer on the Parameters tab. A GUI run, and last-used files written when
  those spinners still existed, restore the tested defaults
  (`applyGuiLockedDefaults`): round 1 reports one candidate, round 2 and the
  closed stage report three, marker search is closed, identification is joint
  closed+open. Enzyme, missed cleavages, genome selector, strategy, threads
  and MSFragger memory stay. CLI and HPC still read the same keys from the
  parameter JSON.

- **The MetaPilot open search's marker search (Step A) is now closed by
  default.** It searches the marker database with a +/-20 ppm precursor
  window, isotope errors 0/1/2, and variable Met oxidation and protein
  N-terminal acetylation, and Crystal-C is skipped for it. On the 12-species
  mix the closed and open marker searches selected the same spiked genomes
  (11.3 and 11.4 per run), and on PXD050781 both selected the target in all 21
  runs. The MSFragger step takes about a fifth of the time (about 3.5 against
  18 min per run). The closed lists carried more non-spiked genomes (32 against
  20 distinct genomes across 9 runs), which genome scoring downstream has to
  remove. The closed lists in this comparison came from the FragPipe workflow
  (Percolator and its protein filter), so the first runs of this build are the
  first measurement under MetaPilot's own stage-1 rescoring. Set
  `round1_search_mode = open` to restore the open search
  (`round1_closed_precursor_ppm` sets the window). The marker search now reruns
  when its parameters change, so existing result folders rerun Step A and
  everything after it on the next run.

- **Open search: one folder per engine configuration.** MSFragger only and
  MSFragger + Sage now share only the per-sample searches, in
  `opensearch_MAG_search/search_results/`. Everything from pooled validation on
  goes to the configuration's own folder, `opensearch_MAG_search/msfragger_only/`
  or `opensearch_MAG_search/msfragger_sage/`:
  - `combined/`, with the validation tables, FlashLFQ output and `variant_*.tsv`;
  - the final tables, `metapilot.psm.parquet` and `metainfo.tsv`;
  - `taxonomy_analysis/`, `functional_annotation/` and `report/`;
  - a copy of `mag_opensearch_analysis_parameters.json`, plus the GUI run's
    `logging.log` and `system_output.log`.

  Previously only validation and FlashLFQ were split (`combined_<mode>/`). The
  final tables, taxonomy, function, report, parameters and log were in one
  shared place, so the two configurations overwrote each other there, and
  the GUI's `combined/` was an empty folder. The parameter file is still also
  saved to the chosen result folder. Existing result folders are not migrated:
  a run against one reuses its searches and writes a fresh mode folder.

- **`OpenSearchPepxmlReader` derives each hit's delta mass** as
  `precursor_neutral_mass - calc_neutral_pep_mass` instead of reading the
  `massdiff` attribute. This turned out to be required, not defensive: on a real
  top-5 Crystal-C output (`YS-194-1-1_c.pepXML`, 2026-09-14) the attribute is
  rewritten on **rank 1 only**. Rank 1 agrees with the derived value to 1e-5 Da,
  while 29,662 of 137,739 rank-2 hits disagree by more than 0.005 Da - up to
  499 Da - and ranks 3-5 likewise, because the spectrum's precursor moved
  underneath a `massdiff` that was left as searched. Reading the attribute would
  have scored corrected precursors against uncorrected deltas for about a fifth
  of the candidate pool. The identity was also checked on an uncorrected pepXML
  (659,521 hits, max deviation 4.7e-4 Da, float rounding), so this is a no-op
  when the step is off.


- **The stage-2 peak gate now tests on files that took no part in finding the
  peaks.** Peaks are located by looking for a target-over-decoy excess, so asking
  "do targets exceed decoys here?" at those locations, on the data that chose
  them, could not fail - and did not: in the 2026-09-10 run all 500 fitted peaks
  were accepted, which is what a test with no power looks like. Discovery now
  uses alternating files from the sorted list and the gate runs on the rest, so
  each side spans every batch. Below four files the split is skipped and the log
  says the p-values are descriptive rather than inferential.
- **The gate runs before scoring, on rank-1 candidates only.** One candidate per
  spectrum keeps the binomial draws independent - with five candidates a single
  spectrum contributes five correlated rows and the p-value is overstated - and
  using the engine's rank-1 rather than the classifier's winner is what makes the
  verdict usable as a classifier feature without circularity. `validate()` no
  longer re-gates on the collapsed set when a gate is supplied, and the verdict
  file now reports the counts and p-values the decision was actually made on.
- **New classifier feature `in_validated_peak`**, distinct from `in_peak`: a
  fitted peak is a location discovery proposed, a validated peak is one the
  held-out files confirmed. Both are kept, because "in a peak that failed the
  gate" is itself informative.
- **Confident positives are seeded only from PSMs inside a validated peak**
  (`restrict_positives_to_validated_peaks`, on by default), so the model learns
  what a correct *and explainable* PSM looks like rather than what a high-scoring
  one looks like. Training still runs over every candidate and every decoy stays
  a negative - narrowing the training set itself would make the peak-context
  features constant and leave PSMs outside every peak with no calibrated score.
- **Round two now reports one candidate per spectrum** (`round2_output_report_topn`,
  default 1, was 5). Measured on 12mix: at every matched entrapment budget,
  rank-1 alone recovered more true PSMs than top-5 (-1.0% at 50k entrapment,
  -6.8% at 259k). The extra ranks open up a population at roughly 50-60%
  genome-attribution purity; the classifier picks well within it (promoted
  rank>=2 winners are right 59.6% of the time against 49.7% for those spectra's
  own rank-1 candidates) but admitting the population at all is what costs the
  error budget. The old justification - peak-discovery votes and classifier
  negatives - does not hold for negatives either: 40.9% of candidates outside any
  peak are already decoys.

- **Stage-1 rescoring and MSBooster no longer degrade silently.** A run
  configured to rescore that quietly falls back is indistinguishable in its
  output from one that worked; a nine-run benchmark was completed on the
  engine-score cut before anyone noticed. Any failure now aborts: a missing
  `.pin`, a missing MSBooster jar or DIA-NN, a non-zero MSBooster exit, an absent
  edited pin, a pin without the expected feature columns, or a classifier error.
  Deliberate opt-out is `stage1_rescoring_enabled=false` /
  `stage1_msbooster_enabled=false`.
- MSBooster is launched the way FragPipe launches it -
  `-Djava.awt.headless=true -cp <MSBooster.jar>;<batmass-io.jar> mainsteps.MainClass`
  - rather than with `java -jar`. Its jar is not shaded, so `-jar` makes the JVM
  ignore every other classpath entry: it started, printed its version banner and
  died at `NoClassDefFoundError: org/slf4j/Logger` without writing anything.
  `batmass-io` is required, not incidental: MSBooster reads the mzML itself to
  compare observed against predicted spectra.
- MSBooster now uses the DIA-NN that FragPipe bundles under `tools/diann/`, not
  the standalone DIA-NN configured for the DIA workflows. MSBooster drives DIA-NN
  through the 1.8.x prediction interface, which DIA-NN 2.x does not provide.
- MSBooster is now given the marker search's `.pin` exactly as MSFragger wrote
  it, `marker/<sample>.pin`, which is what FragPipe does with its own
  `<run>/<run>.pin`. It was previously handed a narrowed copy containing only the
  PSMs stage 1 consumes, written as `<sample>_rescore.pin`; MSBooster pairs a pin
  with its spectra file by base name (`mainsteps.PinMzmlMatcher`), so that suffix
  sent it looking for `<sample>_rescore.mzML` and it aborted with
  `No matching mzml/mgf file` - an error naming the mzML directory, which was
  correct all along. The narrowing bought little in any case: MSBooster fits its
  retention-time calibration on a confident subset it selects itself, and
  features are read back per PSM by file/scan/rank, so rows this stage does not
  consume cost nothing. MSBooster's output now lands as
  `marker/<sample>_edited.pin` with predictions in `marker/MSBooster/`, mirroring
  FragPipe. The runner checks the pin/mzML base-name pairing before launching, so
  a future mismatch is reported against the pin's name rather than the spectra
  folder.
- MSBooster's parameter file now sets `FragmentationType = auto`, as FragPipe's
  own does; left unset it falls back to a fixed assumption. On the 12mix marker
  search this correctly detects LUMOS / HCD 30.
- **Stage-1 MSBooster now runs once over every sample, immediately after the
  marker searches**, rather than once per sample inside genome selection. DIA-NN
  loads its prediction model on every invocation, and the marker search is
  already a single batched MSFragger run, so there is a natural point at which
  every `.pin` exists at once - which is also how FragPipe invokes MSBooster over
  the runs in its manifest. The resulting layout mirrors FragPipe's:
  `msbooster_params.txt` and the shared `MSBooster/` folder of predictions and QC
  plots sit at the search-results root, and each `<sample>_edited.pin` lands
  beside the `<sample>.pin` it came from. The edited pin is that step's resume
  sentinel, so a re-run after genome selection does not re-predict.
- MSBooster is refused a path containing a space, with an explanation, instead of
  failing obscurely. It splits `mzmlDirectory` and `pinPepXMLDirectory` on spaces
  with no quoting or escaping, so such a path is torn into fragments and reported
  as `No matching mzml/mgf file` - which reads as a missing file rather than an
  unusable path.
- The bundled DIA-NN is now chosen as the highest version *below 2.0* rather than
  the highest version outright, falling back to the highest overall with a
  warning. MSBooster drives DIA-NN through the 1.8.x prediction interface, so a
  2.x build newly bundled by a FragPipe upgrade would otherwise be picked and
  fail - the same failure a standalone DIA-NN 2.5.1 caused on 2026-09-10.

- **Pooled validation (step D) no longer reuses PSMs from before the stage-2
  searches that supersede them.** `validated_psms.tsv`, `scored_psms.tsv`,
  `peak_validated_psms.tsv`, `modification_verdicts.tsv` and
  `unexplained_psms.tsv` were skipped on existence alone, so a run that correctly
  re-did stage 1 and re-searched all nine samples then reported the previous
  run's PSMs - and `final_peptides.tsv` and `razor_proteins.tsv` were rebuilt
  from them, with fresh timestamps, describing a search that had been replaced
  hours earlier. They are now re-run when any stage-2 pepXML or Sage result is
  newer.
- **Genome selection no longer reuses a cached result when its evidence has
  changed.** `MagGenomeSelectionService` has a resume check of its own, inside
  the step that calls it: if `*_evidence.tsv`, `*_ranked.tsv` and
  `*_selected.tsv` all exist it returned the previous run's genomes without
  looking at the PSMs it had just been handed. So the first run with stage-1
  rescoring active computed a 37% larger PSM set, wrote it to
  `unmodified_psms.tsv`, and then discarded it: the genomes came back from
  2026-09-10 and the refined FASTA was rebuilt byte-identical, while the
  selection audit files were left describing evidence that no longer existed. A
  caller can now pass `Request.inputs`, and a selection older than any of them is
  recomputed; the open-search workflow passes the PSM file it has just written.
  Callers that pass nothing keep the existence-only behaviour.
- **Resume sentinels now check staleness, not just existence.** Genome selection
  is skipped when the refined FASTA already exists; the stage-2 search is skipped
  when its pepXML already exists. Neither knew that the evidence behind those
  files could change. Switching stage-1 rescoring on after a completed run was
  exactly that: MSBooster wrote new edited pins, genome selection was skipped
  because the old FASTA was still there, and stage 2 re-searched a database built
  by the very cut the rescoring was meant to replace - with nothing in the log or
  the output saying so. Each of those two steps now re-runs when an input it was
  derived from is newer than it (the marker pepXML and the MSBooster edited pin
  for genome selection; the refined FASTA for the stage-2 search), and says which
  input made it stale. Inputs that do not exist are ignored, so an optional step
  that never ran does not force a rebuild.

- **Two of the seven MSBooster feature columns were named wrongly and were
  therefore never read.** `top_6_matched_intensity` and
  `delta_RT_loess_real_units` do not exist; MSBooster writes
  `top6_matched_intensity` and `delta_RT_loess_real`. The underscore-separated
  names come from the camel-case list MSBooster prints in its own log, not from a
  pin header. Until the missing-column check above was added this was silent: the
  classifier trained on five of the seven predicted features and said nothing.
- A sample whose `<sample>_edited.pin` already exists no longer fails when its
  raw `<sample>.pin` is gone. The pin has done its job once MSBooster has run and
  nothing downstream reads it again; at ~54 MB per sample it is exactly the kind
  of file a disk cleanup removes, and a finished sample should not fail over a
  deleted input. The resume check now runs before the input check.
- **MSBooster's predictions are now kept per stage**, in
  `search_results/MSBooster_marker/` rather than `search_results/MSBooster/`, and
  a future stage-2 batch will use `MSBooster_sample_db/`. MSBooster names its own
  output folder and always calls it `MSBooster`, placing it under the deepest
  common directory of the pins it was given - and a stage-one batch over
  `<sample>/marker/*.pin` and a stage-two batch over `<sample>/sample_db/*.pin`
  both resolve to the search-results root, so the second would have silently
  overwritten the first. It takes no output-directory setting, so the folder is
  moved after the run instead.
  - Each sample's QC images are also copied into its own search folder as
    `<sample>/marker/MSBooster_plots/`, so a search folder explains itself
    without a reader having to know where the batch put things.
  - The folder is merged rather than replaced, because a resumed run re-predicts
    only the samples that need it and replacing would discard the QC plots of
    every sample correctly skipped. `msbooster_runs.txt` records what each
    invocation covered, since the pooled `spectraRT*` library then describes only
    the most recent one.
- **Every remaining MSBooster fallback is gone.** With the step switched on it is
  mandatory, and anything missing stops the run before it starts rather than
  producing a result that cannot be compared with anything:
  - no absolute `java` resolvable from `java.home` (the launch would fall through
    to the PATH and give Windows exit 9009, which looks identical to an MSBooster
    failure);
  - no `batmass-io` beside `MSBooster.jar` (the run would print its version
    banner and die at `NoClassDefFoundError`);
  - no DIA-NN 1.8.x bundled with FragPipe - the standalone DIA-NN configured for
    the DIA workflows is no longer accepted as a substitute, because a 2.x build
    predicts nothing while reporting success;
  - an edited pin missing any of MSBooster's seven feature columns, rather than
    training on whichever subset survived;
  - any stage-1 candidate with no row in the edited pin. Round one runs at
    `output_report_topN = 1`, so the pepXML and the pin describe the same PSMs
    and the correspondence is exact; a partial match would score some candidates
    on search features and their neighbours on search plus prediction.
  The only way not to run MSBooster is `stage1_msbooster_enabled=false`.

### Measured

- Pooled peak detection does not broaden peaks. Across 154 peaks measurable in
  all nine files, the pooled-to-within-file width ratio is 1.000 median, 1.053 at
  p90, 1.167 at worst, and between-file drift is 9.3% of peak width (median) -
  which adds 0.4% in quadrature. Per-file calibration offsets span 3e-4 Da
  against a 2e-3 Da peak sigma, because MSFragger calibrates each file. Pooling
  is also what makes weak peaks findable at all: +40.03 Da has 6,468 PSMs pooled
  but roughly 700 per file.
- Peak membership is the strongest single signal in the stage-2 data: the decoy
  rate is 7.4% inside fitted peaks against 40.9% outside.

- On the 12mix marker search (`YS-194-1-1`, 139,496 stage-1 PSMs), MSBooster's
  predicted features are worth substantially more inside the zero family than
  over the pin as a whole - which is where this workflow uses them. Target/decoy
  AUC over the whole pin against the zero family: weighted spectral entropy
  0.643 -> 0.919, hypergeometric probability 0.636 -> 0.900, against 0.603 ->
  0.859 for hyperscore, the best search feature. The spike of PSMs whose
  predicted-spectrum similarity is exactly zero - the open-search caveat, where
  the prediction describes the unmodified peptide and the PSM carries a real
  delta mass - is 16.4% of the pin but 0.4% of the zero family.
- Semi-supervised rescoring of that zero family (3-fold cross-validated) accepts
  5,828 targets at 1% FDR on the engine score alone, 6,434 on search features,
  and 7,422 with MSBooster's features added (+27.4%). Against the 12mix ground
  truth, 90.7% of the PSMs gained map to a 12mix genome, and holding measured
  purity equal to the engine score's own 95.52% still leaves +18.2% true PSMs.
  Unique-peptide evidence for the true genomes rises by 15-25%.
- `unmodified_psms.tsv` records what decided the cut: `stage1_score`
  (`classifier` or `engine_expect`), the feature count, how many of those came
  from MSBooster, and the accepted target/decoy counts.

### Fixed

- **Open search rejected the unmodified peak and oxidation.** The peak filter
  admitted fitted peaks in order of decoy fraction until the admitted pool
  reached 1.5% decoys (`peak_filter_fdr`). On PXD050781 hundreds of small
  peaks with no decoys filled that budget first. The filter then rejected the
  zero peak, which had 1.93% decoys over 592,038 rank-1 PSMs for MSFragger
  and 1.59% for MSFragger + Sage, and never reached oxidation. The
  modifications reported were small peaks of unknown mass. The walk is
  replaced by a per-peak screen:
  - The peak nearest 0 Da is always accepted.
  - A peak matching an expected shift (oxidation, deamidation, pyroGlu,
    water loss, Na-H, Ca-2H, carbamidomethyl; 8 mDa) is accepted when the
    Wilson 95% upper bound on its decoy fraction is below the decoy fraction
    of all other non-zero peaks combined.
  - Any other mass is not accepted. Its PSMs are reported as identifications,
    with no modification claim.

  `peak_filter_fdr` is removed, from the GUI as well; old parameter files that
  contain it still load. `modification_verdicts.tsv` gains `centroid_da`,
  `mass_class`, `known_shift`, `decoy_upper_bound`,
  `background_decoy_fraction`, `fragment_supported` and `reason`, and loses
  `cumulative_fdr`. In `potential_modifications.tsv`, `cumulative_fdr` is
  replaced by `decoy_upper_bound`.

- **Identification-only PSMs had one pooled FDR.** PSMs in rejected peaks and
  PSMs outside every peak are exported under different scopes
  (`REJECTED_PEAK_GROUP`, `NO_PEAK_GROUP`) but shared one q-value walk, so the
  clean rejected-peak PSMs diluted the decoys of the background. Each
  population now has its own walk at the local FDR.

- **The delta-mass evidence gate folded deamidation onto zero.** MSFragger
  does not search shifted fragments between -1.5 and +3.5 Da
  (`delta_mass_exclude_ranges`), so peaks in that range report no fragment
  gain whether or not the shift is real. On PXD050781 the gate folded the
  +0.9836 Da peak (96.5% no-gain) onto the zero peak. Expected shifts are no
  longer rejected or folded by the gate; the peak screen judges them.
  With the gate fixed, the substitution check withdrew the same peak as N->D /
  Q->E (19,485 PSMs, 0.68% decoys): 10.2% of the substituted sequences existed
  in a database holding E. coli and several Blautia relatives, above the 5%
  threshold. N->D and Q->E are the same molecule as deamidated N and Q, so MS/MS
  cannot separate them; they are no longer tested
  (`SubstitutionReassigner.Config.modificationEquivalent`) and the peak is
  reported as deamidation.

- **Open-search folders were moved into `search_results/`.** The legacy
  flat-layout migration moves every unrecognised folder beside
  `search_results/` into it, and in earlier result folders it had moved
  `combined_msfragger_sage/`. The open search now exempts its engine
  configuration folders.

- **Open search with per-engine classification accepted too few genomes.** On
  the 2026-10-08 IPX0003851001 MSFragger-only run, 2 of 32 genomes passed the
  genome q-value gate (MSFragger + Sage passed all 32), and 76% of peptides had
  no confident razor protein. This had two causes:
  - **The calibrated score was being reused for genome scoring.**
    `PerEngineClassifier`, which runs whenever candidates are not combined
    (MSFragger only, or MSFragger + Sage with `combined_psm_instances = false`),
    overwrote `class_logit` with a calibrated score capped at 12, and most
    accepted PSMs reach that cap. Genome scoring read the capped value. The
    calibrated score now has its own field, and only the spectrum-level collapse
    and local FDR read it. `class_logit` stays the real classifier logit.
    Validation is unchanged: the same 858,300 PSMs and 132 peaks.
  - **The genome-scoring null was made of targets.** Open search hands
    post-analysis target PSMs only, so the decoy pool was the roughly 100
    target peptides that also map to a `rev_` protein, and those score like
    targets. Because each genome's null takes the top 5 of M draws, and M can
    reach 14,000 peptides, one high background value is enough to reject every
    large genome. The validated decoys (8,448 PSMs here, with their FlashLFQ
    intensities) are now the genome-scoring background. They are never counted
    or reported, and are not used for protein inference. This is what the pFind
    workflow already does.

  With the real logit alone, 31 of 32 genomes pass. With the real logit and the
  decoy background, all 32 pass in both engine modes, and every peptide has a
  confident razor protein.

- **Delta-mass peak detection was being decided by a cap, not by the data.** The
  12mix run returned exactly 500 peaks against `kdeMaxModes = 500`, with the
  weakest retained peak still carrying 161 target PSMs at SNR 6 - far above every
  evidence gate - so peaks were discarded by a quota and nothing in the task log
  said so. The default is now 2000, and the finder reports how many modes
  qualified, how many were retained, and explicitly says when the cap truncated
  the set.
- **One delta-mass feature was being reported as several peaks.** Five separate
  peaks were emitted inside 0.021 Da where the fitted sigma is 0.002 Da, and the
  same pattern repeated near +2.9 and +3.9 Da. `minPeakSeparationSigmas` existed
  for exactly this but only in the Gaussian-fitting branch, which the default KDE
  path bypasses, so it had never taken effect on a real run. It is now enforced
  on the KDE path: modes closer than the configured multiple of their own width
  are absorbed into the strongest of the group. On the existing peak table this
  merges 150 of 500 modes at 3 sigma.
- Peak **prominence** is now measured and can be required
  (`kde_min_prominence_ratio`, off by default). Between roughly +1.9 and +4 Da
  the delta-mass axis carries a broad elevated region that mode-finding alone
  carved into about twelve peaks per 0.06 Da - every ripple passing the absolute
  support test thousands of times over. Prominence separates a peak standing on a
  background from a ripple in it, and the finder now reports the decile
  distribution of prominence for the modes it kept, so a threshold can be chosen
  from the data instead of guessed.
- An optional width gate (`kde_max_sigma_multiple`, off by default) rejects modes
  much broader than the median mode width - an unresolved blend rather than a
  peak. The scale is taken from the data, so it assumes nothing about the
  instrument.

- Stage-1 rescoring aborted on every run because its feature leak-check flagged
  `num_tol_term`, which is constant on any fully-tryptic search
  (`num_enzyme_termini = 2`) and therefore uninformative rather than leaky. The
  check now flags only a feature constant within one class while varying in the
  other - the actual signature of a label marker - and features constant across
  both classes are dropped with a note.

### Known

- `peak_r2` is NaN for every peak on the default path: the KDE branch never fits
  a Gaussian, and it is one of the 36 features the stage-2 classifier trains on.
  The same branch makes `minR2`, `minFitBins`, `minSigmaBins`, `maxFitIterations`,
  `maxPeaksToFit`, `smoothHalfWidth`, `smoothPolyOrder` and
  `minPeakSeparationBins` inert - they remain exposed in the parameter file and
  have no effect.
- 137 of the 268 monoisotopic peaks in the 2026-09-10 run sit at mass defects no
  elemental composition can produce, yet are reproducible in 9 of 9 files and
  carry 15% of all PSMs assigned to a peak. Reproducibility, scale persistence
  and target-decoy enrichment do not separate real chemistry from systematic
  artefacts, because a systematic artefact is reproducible by construction.

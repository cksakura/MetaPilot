# Changelog

All notable user-facing changes to MetaPilot are documented here.

## 1.2.0 - Unreleased

Planned release after the open-search strategy is finalized.

### Added

- **PSM coverage option for DDA MAG genome selection (pFind, FragPipe, and open search).** Evidence and ranking are unchanged. Automatic (Focus or Comprehensive) remains the default and, with default settings, selects the same genomes as before this option. "At a PSM coverage" keeps genomes until a user-set share of the marker PSMs is explained, with the coverage and step-guard rule adopted from MetaLab-MAG (J Proteome Res 2023). It is not the MetaLab-MAG workflow: MetaPilot still ranks by marker-peptide gain. The control is labelled "When to stop adding genomes".
- Parameter keys `genome_selection_strategy` (`PSM coverage`), `genome_selection_psm_coverage` (default 0.95), `genome_selection_step_guard` (default 0.001; 0 switches the guard off; allowed 0 to 0.002) and `genome_selection_max_genomes` (default 300). The step-guard control shows a percentage (default 0.1%, range 0–0.2%, steps of 0.01 percentage points) and the parameter file stores the fraction. A parameter file without the new keys keeps automatic Focus or Comprehensive selection.
- CLI flags on `pfind-mag`: `--genome-selection-strategy`, `--psm-coverage`, `--step-guard`, `--max-genomes`.
- Selected-table footer rows `#SelectionStrategy`, `#PsmCoverageThreshold`, `#StepGuard`, `#MaxGenomes`, `#StopReason`, `#PsmCoverageReached`. Ranked tables gain `CutStep` and `InRanking`. Each run writes `pfind_MAG_search/genome_selection_summary.tsv`, one row per raw file. After a pFind run the GUI shows PSM coverage against genomes in rank order.

Sage uses the same genome-selection service, but its panel is not changed here. DIA evidence is precursors, so a coverage option for DIA would need its own label and testing.

- CLI-only `mag-opensearch-entrapment`: adds paired shuffled proteins after genome
  selection to the shared second-stage closed/open database, with decoys of both
  populations. Uses DIA Deep's reproducible shuffle, isolated evaluation outputs,
  hash-verified database manifests and all four identification arms. Stops after
  identification; the descriptive accession counts do not claim calibrated FDP.
  The MAG open-search GUI now has the same evaluation switch (off by default).

- DIA Marker Deep together GUI: off-by-default **Paired-shuffle entrapment
  (evaluation)** switch. Same processing as `dia-marker-deep-entrapment`; results
  go to `marker_deep_entrapment/` and do not overwrite a production Deep together
  run. Hidden unless the workflow is Deep together.

### Changed

- Open-search sequence-variant tables (`sequence_variants.tsv`, `sequence_variants.fasta`, the UHGP `variant_*.tsv` files, and new `substitution_forms.tsv`) are written to `opensearch_MAG_search/sequence_variants/`. A `substitution_forms.tsv` already in `combined/` is left there. `validated_psms.tsv` and the rest of the MAG layout stay where they are.
- Open search and FragPipe write `ms2Count.txt` and `rt.txt` in the result's `spectra/` folder, the same place pFind uses. A count left in `spectra_summary/` is moved there on the next run and is still read if the move cannot be done.
- `excluded_evidence/` stays at the search-folder root, beside `taxonomy_analysis/`. Layout migration no longer treats it as a sample folder, and a result clean keeps it with the other final folders.
- Open search writes `sequence_variants/substitution_psms.tsv` for substitution PSMs that are absent from `validated_psms.tsv`, and quantifies that file with FlashLFQ with match-between-runs off. The peptide table is `sequence_variants/substitution_quantified_peptides.tsv`. It is not read by `final_peptides.tsv` or the other final tables. A failure of this run does not fail the workflow, and the spectrum folder is restored afterwards, including `ExperimentalDesign.tsv`.

### Fixed

- **Open-search genome selection no longer discards the marker evidence of sparse samples.** Stage-1 rescoring fitted its logistic model with almost no regularisation. On marker searches with a few hundred confident PSMs the coefficients diverged, target probabilities rounded to exactly 1.0, decoys tied with targets, and rescoring kept 0 PSMs, so no sample database was built. On PXD066732 this zeroed files that the engine-score cut had given 120 to 890 marker PSMs. Stage 1 now uses ridge 0.1 and ranks on the logit; replayed on 20 PXD066732 marker searches, every collapsed file kept at least its engine-score count, and the largest files moved by under 1.5%. Identification (stage 2) is unchanged. An existing sample database is reused: delete `marker/<sample>.fasta` to select a sample again, and its stage-2 searches then rerun because they are older than the new database. A sample skipped for lack of stage-1 evidence has no database and is selected again on the next run.
- **Open-search identification reproducible across runs.** DIA-NN 1.8.1
  predicts slightly different intensities and iRT on every multithreaded run.
  These values reached the classifier as features, so two identifications of
  the same searches differed: 547,604 vs 547,770 faecal PSMs, enough to move
  genomes across the 1% gate. Predictions now run as fixed 5,000-precursor
  slices, each in its own single-threaded DIA-NN, several at a time. The
  output depends on the input alone, not on the thread count.
- **CLI log size.** Without a configuration file, logback (bundled through a
  dependency) logged at DEBUG, and the Parquet export wrote one line per field:
  49 GB for one faecal run. The CLI now ships a `logback.xml` that keeps
  Parquet and Hadoop at WARN.
- **DDA MAG identification rate after nested `search_results/` layout**
  (`final_summary.tsv`). `spectra_summary/` (`ms2Count.txt`) was treated as a
  leftover sample folder and moved under `search_results/`, so the summary
  writer fell back to identified spectra as the MS/MS count and reported 100%.
  The folder now stays at the engine root; the writer also finds a previously
  nested copy. FragPipe MAG started from mzML (no conversion folder) writes the
  counts there instead of skipping them. Open search and GUI MAG share the same
  writer.

- **DDA MAG FragPipe identification rate** (`final_summary.tsv`). FragPipe
  renames runs (`Qiao-FB1-DDA-1` becomes `Qiao_FB1_DDA_1_1`), so no MS2 count in
  `ms2Count.txt` matched and the rate fell back to 100%. Both name forms now
  match.

- **DDA MAG per-file "Peptide Sequences Identified"** (`final_summary.tsv`,
  `id.tsv` "Unique peptide count"). It also counted peptides that FlashLFQ
  quantified in a file by match-between-runs, so a file could report more
  peptides than identified spectra. It now counts MS/MS identifications only;
  totals and `final_peptides.tsv` are unchanged.

## 1.1.2 - 2026-10-14

### Added

- **Decoy genome calibration in `genome_rank_scores_runinfo.txt`** (DDA MAG
  genome scoring). Each reversed proteome's decoy peptides go through the same
  top-5 statistic and resampling as a target genome, and the run info reports
  how many decoy genomes reach p <= 0.01 next to how many target genomes do,
  with one row per reversed genome. Decoy genomes over target genomes estimates
  the genome-level FDR of the threshold. On PXD050781 no reversed genome of 35
  reaches it against 39 target genomes. Calibration only: decoy genomes are
  never selected. It cannot see a real peptide attributed to the wrong genome.

- **Pangenome variant check** (`pangenome_variant_check`, `auto` by default;
  pFind open search in the GUI only). pFind already proposes single-residue
  substitutions, as `Ser->Thr[S]`-style tokens in `pFind-Filtered.spectra`,
  but cannot say whether the substituted peptide exists in nature. The
  database it searched can't answer that either: genome selection kept one
  strain, and the peptide may come from a
  relative. This check asks the UHGP-100 protein catalogue. It looks up each
  variant peptide, its base peptide and a matched control in every UHGP-100
  protein, in tryptic context by default (`pangenome_require_tryptic_context`,
  I/L folded, initiator-Met loss allowed). It resolves hits to species
  representatives and reports per candidate whether the variant is in the
  searched database, anywhere in UHGP-100, or in a species this run detected.
  "Searched database" is every database the run searched, across all samples
  (`in_searched_db`, `SUBSTITUTION_IN_SEARCHED_DB`). A variant found only in
  another sample's database counts too, although the candidate's own spectrum
  could not have matched it directly.
  Detected species are the representatives of the genomes in
  `taxonomy_analysis/genome_final_selected.tsv`, so the check runs after
  taxonomy. `in_same_species` narrows it to the species of the PSM's own
  proteins (`assigned_species`). A variant that other genomes of that species
  carry is the most parsimonious explanation the check can offer: natural
  sequence variation within the species, not a modification. It is not a claim
  about which strain is present.

  When MGnify's cluster membership table `uhgp-100.tsv` is available, the
  check also counts genomes. It looks at `uhgp100_membership_path` if set.
  Otherwise it looks beside `uhgp-100.faa.gz`, then in the `uhgp-100/` folder
  that extracting the archive creates, either inside `protein_catalogue/` or
  next to it. A `.tsv.gz` copy works too. The FASTA holds
  only cluster representatives, so the lookup records which representatives
  contain each variant. One streaming pass over the table then collects
  their member genomes. The pass adds the time to read 25 GB, about 4 minutes
  at 100 MB/s. Per candidate, `variant_candidates.tsv` gains:
  `n_genomes` (any species), `n_genomes_same_species`, `n_genomes_in_species`
  (every catalogued genome of the assigned species), `n_isolates_same_species`
  (from `Genome_type`) and `example_genome_same_species`, which prefers an
  isolate. For example, "carried by 40 of 312 genomes of S, 5 of them
  isolates". Members can be shorter than their representative, because MGnify
  clustered with a coverage threshold, so the counts are an upper bound.
  Without the table these columns are NA.

  Two new files in the combined result folder: `variant_candidates.tsv`, one row
  per candidate with its `pangenome_class`, and `variant_pangenome_summary.tsv`,
  percentages over unique variant peptides per engine, substitution and
  MS2 decidability. Each summary row is set against a matched control, so a
  percentage can be read as enrichment over chance rather than taken at face
  value. The control is the same X->Y applied at a random X in a random
  identified peptide, seeded by `pangenome_control_seed`, with
  `pangenome_controls_per_candidate` controls per candidate. A candidate is
  `ms2_decidable = false` when a same-residue modification has the same mass
  (S->T against methylation of S, N->D against deamidation, ...). pFind's alias
  form `Methyl[S](Ser->Thr[S])` is always undecidable.

  The resource is the EBI UHGP-100 download (`uhgp-100.tar.gz`) plus
  `genomes-all_metadata.tsv`. It is set up in the MAG database configuration,
  human-gut only. One option is the opt-in checkbox "Also download the UHGP-100
  protein catalogue" next to Download/Update. The other is "Prepare protein
  catalogue" for a database that is already downloaded. On the CLI the
  equivalents are `mag-download --protein-catalogue` and
  `--protein-catalogue-only`. The 28 GB archive download resumes if
  interrupted. A one-off preparation of about 45 minutes then streams it into
  `uhgp-100.faa.gz` (about 23 GB), `uhgp-100.rep_species.tsv.gz`
  (multi-species clusters only) and a manifest, all under
  `<catalogue>/<version>/protein_catalogue/`, where the check looks by
  default. The archive is deleted afterwards unless "Delete the archive after
  preparation" is unticked (`--keep-protein-archive`). A failure there leaves
  the MAG database download successful. `uhgp100_path` points elsewhere if
  needed. At analysis time `auto` uses a prepared resource or logs why it
  skipped; `on` also prepares from an archive it finds.

  **What this does not do.** It is report-only. Nothing it finds reaches the FDR
  or any other existing output. With the check `off`, or the resource absent,
  every existing result file is byte-identical to before. It does not decide
  which form a spectrum carries: presence in UHGP-100 makes a variant
  plausible, not identified, and an undecidable candidate stays undecidable.
  pFind calls are taken as pFind made them.
  The headless pFind runner is untouched. The pFind GUI has no controls for the
  new keys yet; they are always written to the parameter JSON and edited there.
  Host-only PSMs are left out because UHGP-100 is microbial.

- **Per-species summary and assigned-protein function for the pangenome
  variant check.** The searched database and the genome report hold species
  representatives only, but a variant is often carried by other genomes of the
  same species. A third file, `variant_species_summary.tsv`, beside the other
  two, has one row per species representative that genome selection kept
  (`in_genome_final_selected`) or that a substitution PSM was assigned to. It gives the
  species' catalogued genomes and isolates; its candidates by
  `pangenome_class`, with `n_candidates_shared` for those also assigned to
  other species, under which they count too; and, over the distinct
  `SUBSTITUTION_IN_PANGENOME` variant peptides, how many UHGP-100 has in this
  species and how many only in others. With the membership table, it also gives
  how many of the species' genomes other than its representative carry at
  least one of the this-species variants, and how many of those are isolates.
  Like the per-candidate counts, these are upper bounds. They list compatible
  genomes, not strains present, and are not ranked.

  `variant_candidates.tsv` gains `assigned_kegg_ko`, `assigned_cog`,
  `assigned_cog_category`, `assigned_preferred_name` and
  `assigned_description`. They come from the eggNOG tables the functional
  annotation reads (`<genome>_eggNOG.tsv`), as a union over the PSM's assigned
  proteins, and use the same COG derivation as `Annotated_proteins.tsv`. This
  is the function the rest of the report already counts the PSM under, not an
  annotation of the variant's carrier proteins. A non-representative genome has
  no eggNOG table, and the only UHGP catalogue MGnify annotates is UHGP-90.
  The columns are empty when the tables are read but the protein is
  unannotated. They are NA when there is no eggNOG table for any assigned
  protein.

  Results from before this change lack the new file. The next run notices that
  and reruns only the check.

- **Pangenome variant check: sites, PSM identifiers, statistics and a column
  dictionary.** Two more files sit beside the other three.
  `variant_sites.tsv` has one row per protein site and substitution, with the
  site's protein position, peptides, PSM, file and class counts. The
  searched-database pass now records where each base peptide sits in its
  assigned proteins (first occurrence). A peptide not found in an assigned
  protein gets a per-peptide row with `located = false`.
  `variant_columns.tsv` says what every column of the four tables means.

  `variant_candidates.tsv` changes:
  - `psm_id` is replaced by `charge`, which with `file` and `scan` joins back
    to the workflow's PSM table.
  - New PSM columns: `precursor_mz` and `psm_q_value`.
  - New columns for the PSM's proteins: `assigned_proteins`, `protein_sites`
    (protein:position), and `residue_before`/`residue_after`.
  - `n_substitutions_in_variant_peptide` marks a pFind multi-site peptide,
    whose sites apply together.
  - New `carrier_fraction_in_species`.
  - `front_end_evidence` is split into `front_end_call` and
    `localized_sequence`.

  `variant_pangenome_summary.tsv` changes:
  - New `level` column.
  - Found counts appear beside every percentage.
  - New one-sided Fisher exact p-values, with Benjamini-Hochberg q-values
    within each level. Unique peptides are treated as independent.

  In `variant_species_summary.tsv`, `n_candidates_*` count candidates and
  `n_peptides_*` count unique variant peptides, now for every class.
  `selected` is renamed `in_genome_final_selected`: the post-analysis
  selected the species, which does not prove it is present.

  Results without `variant_sites.tsv` rerun the check, as above.

### Changed

- **Genome scoring counts each peptide once.** The unit was modified sequence
  plus charge, so one peptide seen at two charges or with an oxidation filled
  several of the top-5 slots; most genomes' top 5 were 2-3 distinct peptides.
  The unit is now the base peptide (I = L), weighted by its strongest form.
  On PXD050781 no genome's p-value changed.

- **Resumed runs redo only what depends on a missing or changed file**, in the
  pFind MAG workflow. An output is current when it exists and is no older than
  what it is derived from.
  - Step 7's post-analysis is skipped when its outputs are current against
    `combined/psms.tsv` and the quantified peptides.
  - HTML reports are skipped when every report in `report/reports/` is no
    older than its source table, `metainfo.tsv` is unchanged in content, and
    `manifest.json` exists.

  Deleting `variant_candidates.tsv` and `variant_pangenome_summary.tsv`
  therefore reruns the variant check and nothing else. Parameter changes are
  not detected: after changing post-analysis or report settings without
  changing their inputs, delete the affected outputs.

- **DIA-NN protein matrices use a 5% run-specific protein q-value**
  (`--matrix-spec-q 0.05`, DIA workflows that write matrices). It was set to
  the search q-value, usually 1%, which removed proteins from a run's matrix
  column whenever their evidence in that run alone was weaker than 1%. 5% is
  DIA-NN's own default and the run-specific protein filter its documentation
  recommends (`PG.Q.Value` at 0.01 to 0.05, "typically 0.05 is sufficient").
  The precursor q-value (`--qvalue`) and the global matrix q-value
  (`--matrix-qvalue`) still follow the search q-value.

- The default enzyme for DDA MAG FragPipe is "Trypsin/P"
  (MSFragger `stricttrypsin`): cleavage C-terminal to K/R including before
  proline. This matches MaxQuant's and FragPipe's own defaults, and is close to a
  superset of the proline-restricted digest - where the proline rule leaves a
  K-P-spanning peptide as the only candidate, strict digestion offers that
  peptide and both split forms. Select "Trypsin" for the conventional Keil-rule
  digest, which now genuinely performs one.

### Fixed

- **Target peptides that also map to a reversed protein entered the genome
  and protein null.** Genome and protein scoring counted a peptide as decoy
  when any of its proteins was `rev_`. Short target sequences such as
  `ARIPLAK` also occur in reversed proteins, so they joined the resampling
  pool with target-level scores, and genomes with a few thousand unique
  peptides drew them on every resample, which could push a present genome's
  p-value above the threshold. A peptide is now a decoy only when all of its
  proteins are decoys. This
  applies to every DDA MAG workflow that uses PSM-based genome scoring.

- **The GUI pFind and Sage MAG workflows never wrote `metapilot.psm.parquet`.**
  Their shared post-analysis task did not call the QPX PSM export; only the
  headless pFind export and FragPipe did. It now writes
  the file beside `final_peptides.tsv`. When pFind skips post-analysis because
  the final tables are current, it writes the parquet from `combined/psms.tsv`
  if the file is missing or older than that table.

- **`metapilot.psm.parquet` was never written on Windows** (every DDA workflow
  that exports it). Parquet wrote through Hadoop's local file system,
  which sets the new file's permissions with `winutils.exe`; MetaPilot does not
  ship it, so the writer failed right after creating an empty
  `metapilot.psm.parquet.tmp` and the export was skipped with a log line. The
  QPX PSM writer, `TypedParquetTableWriter` and `TsvToParquetExporter` now write
  through a plain-Java Parquet `OutputFile` (`LocalParquetOutputFile`). Verified
  from the shaded GUI jar, which reproduced the empty `.tmp` before the change.

- **The pFind workflow's pangenome variant check was silent when it did not
  run.** When the check was off, the selected MAG database had no UHGP-100
  resource, its outputs were already current, or it failed, the reason went to
  the log file only, never to the console. Every outcome is now also printed
  to the console, and a closed search says the check runs only for open search.
  `pangenome_variant_check = off` now logs a line instead of returning silently.

- **Enzyme selection: the "Trypsin"/"Trypsin/P" names were inverted in the
  MSFragger-family DDA MAG workflows.** "Trypsin" performed a Trypsin/P digest -
  cleaving after K/R even before proline - and "Trypsin/P" performed a
  conventional trypsin digest. "LysC"/"LysC/P" carried the same inversion. The
  cause was the enzyme dropdown being populated from `enzyme_default.ini`, which
  is pFind's configuration vocabulary and declares `Trypsin` with no proline
  restriction; the MSFragger mappings then reproduced that naming, contradicting
  MetaPilot's own `core.enzyme.Enzyme` enum and its DIA-NN mapping, both of which
  use the conventional meaning. The names now mean what they say.

  **The default digest is unchanged.** The workflow default is now "Trypsin/P",
  which is exactly what the old "Trypsin" performed, so a new run with default
  settings produces the same digest as before and existing results stay
  comparable.

  **Saved parameter files are migrated on load.** A file written before this
  change that says "Trypsin" is loaded as "Trypsin/P" - the enzyme it was
  actually written with - and a warning is logged. Re-saving makes it permanent.
  Files written from now on are marked and never re-interpreted.

  Affects the DDA MAG FragPipe workflow. The pFind workflows keep pFind's own
  names, where those semantics are the engine's own (see the next entry). The
  standalone Sage and IGC FragPipe workflows keep their own separate enzyme
  mappings and are unchanged.

- **The pFind workflows ignored the selected enzyme.** In the DDA MAG pFind and
  DDA IGC pFind workflows the enzyme was saved but never written to
  `pFind.cfg`, so every search used the enzyme in pFind's own parameter
  template, `Trypsin KR _ C`, whatever was selected. Samples digested with
  another enzyme, such as Chymotrypsin or Asp-N, were searched as tryptic. The
  selected enzyme is now written to `pFind.cfg`.

  **The default digest is unchanged.** The default is pFind's "Trypsin"
  (`KR _ C`, cleaving after K/R also before proline), which is what every pFind
  search used. pFind names the enzymes the other way round from the
  MSFragger-family workflows: pFind's "Trypsin_P" (`KR P C`) is conventional
  trypsin, and a name ending in "_P" does not cleave before proline. The
  dropdown tooltip says so. Parameter files that say "Trypsin_P" now get that
  digest, not the template's. An enzyme name pFind does not define stops the
  run before searching instead of falling back to trypsin.

- **The DDA MAG FragPipe workflow ignored the modifications and the digestion
  mode.** The modifications in the panel were saved but never written to the
  FragPipe workflow, so every search used the workflow template's own:
  carbamidomethyl C fixed, and oxidation on M and protein N-terminal
  acetylation variable. The digestion mode was also ignored and every search
  was fully specific. Both are now written to the workflow at every search
  stage, in the GUI and headless runs. The digestion choices are Specific, Semi
  specific and Unspecific; the older "Semi free N-term"/"Semi free C-term"
  values are searched as semi-specific, because MSFragger cannot restrict
  semi-specificity to one terminus.

  **Default searches are unchanged.** The default modifications are the
  template's three, and the default digestion mode is Specific. A headless
  parameter file saved before this change keeps searching the template's
  modifications and logs a warning; saving it again applies the listed ones.

  Two related fixes. Asp-N was searched as cleaving after D; it now cleaves
  before D. The shuffled decoys and host database used the cleavage residues
  of an outdated enzyme table in the panel; they now use the same enzyme
  definition as the search. Loading saved parameters also used to reset the
  enzyme, missed cleavages, digestion mode, MS data type and modifications in
  the panel; they are now restored.

- **DIA Target and DIA RealTime spectral libraries ignored the modification
  checkboxes.** The Target library and the RealTime per-file libraries were
  predicted with DIA-NN's built-in defaults (N-terminal Met excision and
  carbamidomethyl C) whatever was ticked. They now use the selected
  modifications. With the default ticks the libraries are unchanged.

- **DIA-NN: Chymotrypsin and Asp-N used the wrong cleavage rule.**
  "Chymotrypsin" was passed to DIA-NN as `--cut null`, because the name lookup
  was case-sensitive and expected "ChymoTrypsin"; the lookup now ignores case,
  and an enzyme name it does not know stops the run instead of writing `null`.
  Asp-N was `D*` (cleave after D); it is now `*D` (cleave before D), and
  MetaPilot's missed-cleavage count in the DIA-NN result readers understands
  the "cleave before" form.

- **The DDA MAG pFind workflow dropped the isobaric correction factors.** The
  values entered in the isobaric table were never passed to the search, so
  reporter-ion intensities were not corrected for isotope impurities. They are
  now applied, as the DDA IGC pFind workflow already did. On first launch,
  without saved parameters, the search mode also showed Open search; it now
  defaults to Closed search, which is what loading a saved file already used.

### Improved

- The enzyme mapping for MSFragger-family workflows now has a single definition
  (`MSFraggerHeadlessRunner.enzymeInfo`), shared by the FragPipe MAG runner and
  its dropdown. It previously existed as two copies with different
  normalisation rules that happened to agree.

## 1.1.1 - 2026-08-12

### Fixed

- Fixed the FragPipe DDA open-search workflow so MAG FragPipe open search runs correctly.

## 1.1.0 - 2026-07-23

### Added

- Representative-peptide clustering for high-throughput DIA deep-search workflows, with configurable presets and similarity thresholds to reduce redundant genome searches.
- A workflow-aware result-folder cleaner with preview, Basic, and Deep modes. The cleaner preserves resumability checkpoints and final deliverables according to the selected workflow.
- Expanded CLI and HPC support for DIA Marker and MAG workflows, including additional post-analysis and database-preparation utilities.

### Improved

- Progress reporting now follows completed workflow work more closely across long-running GUI and headless tasks.
- DIA and MAG result layouts, restart behavior, parameter handling, and post-analysis exports are more consistent.
- Report generation and peptide/PSM parsing are more robust when inputs are incomplete or differ between search engines.
- User-facing validation, status messages, and workflow configuration controls have been refined.

### Fixed

- Fixed progress bars that could advance independently of the real workflow state.
- Fixed resume and cleanup edge cases that could remove files still needed by later post-analysis steps.
- Fixed several DIA-NN, FragPipe, pFind, and report-generation compatibility issues.

## 1.0.0

- First formal MetaPilot release.

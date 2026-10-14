# Changelog

All notable user-facing changes to MetaPilot are documented here.

## Unreleased

### Added

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

- **Decoy genome calibration in `genome_rank_scores_runinfo.txt`** (DDA MAG
  genome scoring). Each reversed proteome's decoy peptides go through the same
  top-5 statistic and resampling as a target genome, and the run info reports
  how many decoy genomes reach p <= 0.01 next to how many target genomes do,
  with one row per reversed genome. Decoy genomes over target genomes estimates
  the genome-level FDR of the threshold. On PXD050781 no reversed genome of 35
  reaches it against 39 target genomes. Calibration only: decoy genomes are
  never selected. It cannot see a real peptide attributed to the wrong genome.

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

- **Pangenome variant check** (`pangenome_variant_check`, `auto` by default;
  pFind open search in the GUI and the MetaPilot open search only). Both
  engines already propose single-residue substitutions: pFind as
  `Ser->Thr[S]`-style tokens in `pFind-Filtered.spectra`, the open search as
  `SubstitutionReassigner` proposals. Neither can say whether the substituted
  peptide exists in nature. The database each one searched can't answer that
  either: genome selection kept one strain, and the peptide may come from a
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
  `taxonomy_analysis/genome_final_selected.tsv`, so the open-search check runs
  after taxonomy. `in_same_species` narrows it to the species of the PSM's own
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
  `pangenome_controls_per_candidate` controls per candidate. Open-search rows
  are further split by mass-shift peak and `peak_substitution_shaped`. The
  proposals the reassigner computed on peaks it rejected are reported too, so
  the modification-like peaks (+0.9841, +40.0313, ...) can be compared with the
  substitution-shaped ones on the same footing. A candidate is
  `ms2_decidable = false` when a same-residue modification has the same mass
  (S->T against methylation of S, N->D against deamidation, ...). The
  isobaric table behind `sequence_variants.tsv` is now the shared
  `IsobaricModificationTable`, and the run log prints it once. pFind's alias
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

  **What this does not do.** It is report-only. Nothing it finds reaches the FDR,
  the substitution verdict, `validated_psms.tsv` or any other existing output.
  With the check `off`, or the resource absent, every existing result file is
  byte-identical to before. It does not decide which form a spectrum carries:
  presence in UHGP-100 makes a variant plausible, not identified, and an
  undecidable candidate stays undecidable. For the open search, the MS2 form
  ranking below adds a secondary column on that question. pFind calls are
  taken as pFind made them.
  The headless pFind runner is untouched. The pFind GUI has no controls for the
  new keys yet; they are always written to the parameter JSON and edited there.
  Host-only PSMs are left out because UHGP-100 is microbial.

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
  and reruns the check: pFind reruns only the check, and open search reruns
  Step D. Step D is deterministic, so `validated_psms.tsv` keeps its timestamp
  and quantification and export are not redone.

- **Pangenome variant check: sites, PSM identifiers, statistics and a column
  dictionary.** Two more files sit beside the other three.
  `variant_sites.tsv` has one row per protein site and substitution, with the
  site's protein position, peptides, PSM, file and class counts. The
  searched-database pass now records where each base peptide sits in its
  assigned proteins (first occurrence). A peptide not found in an assigned
  protein gets a per-peptide row with `located = false`.
  `variant_columns.tsv` says what every column of the four tables means.

  `variant_candidates.tsv` changes:
  - `psm_id` is replaced by `hit_rank` (open search) and `charge`. Either one,
    with `file` and `scan`, joins back to the workflow's PSM table.
  - New PSM columns: `precursor_mz` and `psm_q_value`.
  - New columns for the PSM's proteins: `assigned_proteins`, `protein_sites`
    (protein:position), and `residue_before`/`residue_after`.
  - `n_substitutions_in_variant_peptide` and `n_candidates_for_psm` tell a
    pFind multi-site peptide (all sites apply together) from the open-search
    reassigner's alternative readings (only one can be right).
  - New `carrier_fraction_in_species`.
  - `front_end_evidence` is split into `front_end_call`,
    `localized_sequence` and `reassigner_in_db`.

  `variant_pangenome_summary.tsv` changes:
  - New `level` column.
  - For open search, rows pooled over peaks come before the per-peak rows.
  - Found counts appear beside every percentage.
  - New one-sided Fisher exact p-values, with Benjamini-Hochberg q-values
    within each level. Unique peptides are treated as independent.

  In `variant_species_summary.tsv`, `n_candidates_*` count candidates and
  `n_peptides_*` count unique variant peptides, now for every class.
  `selected` is renamed `in_genome_final_selected`: the post-analysis
  selected the species, which does not prove it is present.

  Results without `variant_sites.tsv` rerun the check, as above.

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

- **Genome scoring counts each peptide once.** The unit was modified sequence
  plus charge, so one peptide seen at two charges or with an oxidation filled
  several of the top-5 slots; most genomes' top 5 were 2-3 distinct peptides.
  The unit is now the base peptide (I = L), weighted by its strongest form.
  On PXD050781 no genome's p-value changed.

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
- **Resumed runs redo only what depends on a missing or changed file**, in the
  open-search (GUI and headless) and pFind MAG workflows. An output is current
  when it exists and is no older than what it is derived from.
  - Open search: when Step D reruns and writes `validated_psms.tsv` with the
    same bytes, the file keeps its previous timestamp, so FlashLFQ is skipped by
    its existing check. Steps G-I (final tables, taxonomy, functional
    annotation) are skipped together when `final_summary.tsv`,
    `final_peptides.tsv`, `razor_proteins.tsv`, `all_proteins.tsv`,
    `taxonomy_analysis/genome_final_selected.tsv`, `Taxa_report.tsv` and (with a
    function database) `functional_annotation/Annotated_proteins.tsv` are
    current against `validated_psms.tsv` and `QuantifiedPeptides.tsv`.
    They are reused only when `post_analysis_source.tsv`, written after Steps
    G-I, names those two files (paths relative to the mode folder), so an
    interrupted export is redone. Results from before this record existed are
    regenerated once.
  - pFind: step 7's post-analysis is skipped under the same rule, against
    `combined/psms.tsv` and the quantified peptides.
  - HTML reports (both workflows) are skipped when every report in
    `report/reports/` is no older than its source table, `metainfo.tsv` is
    unchanged in content, and `manifest.json` exists.

  Deleting `variant_candidates.tsv` and `variant_pangenome_summary.tsv` therefore
  reruns Step D (open search; its substitution result lives only in memory)
  and the variant check, and nothing else. Parameter changes are not detected:
  after changing post-analysis or report settings without changing their
  inputs, delete the affected outputs.

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

- **Target peptides that also map to a reversed protein entered the genome
  and protein null.** Genome and protein scoring counted a peptide as decoy
  when any of its proteins was `rev_`. Short target sequences such as
  `ARIPLAK` also occur in reversed proteins, so they joined the resampling
  pool with target-level scores. On PXD050781 the MSFragger + Sage pool held
  553 of them, 68 at target-level weight, so genomes with a few thousand
  unique peptides drew them on every resample. `MGYG000002478` went from
  p_emp 0.004 to 0.18, and 18 genomes that MSFragger open kept were dropped.
  A peptide is now a decoy only when all of its proteins are decoys. This
  applies to every DDA MAG workflow that uses PSM-based genome scoring.

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

- **The GUI pFind and Sage MAG workflows never wrote `metapilot.psm.parquet`.**
  Their shared post-analysis task did not call the QPX PSM export; only the
  headless pFind export, open search and FragPipe did. It now writes
  the file beside `final_peptides.tsv`. When pFind skips post-analysis because
  the final tables are current, it writes the parquet from `combined/psms.tsv`
  if the file is missing or older than that table.

- **`metapilot.psm.parquet` was never written on Windows** (every DDA workflow
  that exports it; noticed in open search). Parquet wrote through Hadoop's local file system,
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

### Fixed

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

  Affects the DDA MAG FragPipe and DDA MAG open-search workflows, including the
  Sage engine within the open search. The pFind workflow is unaffected: it keeps
  reading `enzyme_default.ini`, where those semantics are the engine's own. The
  standalone Sage and IGC FragPipe workflows keep their own separate enzyme
  mappings and are unchanged.

### Changed

- The default enzyme for DDA MAG FragPipe and DDA MAG open search is "Trypsin/P"
  (MSFragger `stricttrypsin`): cleavage C-terminal to K/R including before
  proline. This matches MaxQuant's and FragPipe's own defaults, and is close to a
  superset of the proline-restricted digest - where the proline rule leaves a
  K-P-spanning peptide as the only candidate, strict digestion offers that
  peptide and both split forms. Select "Trypsin" for the conventional Keil-rule
  digest, which now genuinely performs one.

### Fixed

- Stage-1 rescoring aborted on every run because its feature leak-check flagged
  `num_tol_term`, which is constant on any fully-tryptic search
  (`num_enzyme_termini = 2`) and therefore uninformative rather than leaky. The
  check now flags only a feature constant within one class while varying in the
  other - the actual signature of a label marker - and features constant across
  both classes are dropped with a note.

### Changed

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

### Improved

- The enzyme mapping for MSFragger-family workflows now has a single definition
  (`MSFraggerHeadlessRunner.enzymeInfo`), shared by the FragPipe MAG runner, the
  open-search runner and Sage-within-open-search, and by the dropdowns those
  panels show. It previously existed as two copies with different normalisation
  rules that happened to agree.

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

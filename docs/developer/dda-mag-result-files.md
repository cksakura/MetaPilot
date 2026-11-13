# DDA MAG (closed search) result files

What the standard DDA MAG workflow writes, step by step, and what each file is for. This is the
pFind workflow (marker search -> sample-specific database -> second search -> FlashLFQ ->
post-analysis); FragPipe differences are at the end. It does **not** cover the MetaPilot open
search (`MagOpenSearchHeadlessRunner`, `<result>/opensearch_*`), which has its own layout.

Code entry points:

- GUI (all steps): `metapilot-windows-gui/.../task/dda/mag/task/DDAMagPfindTask.java`, base class
  `metapilot-windows-gui/.../task/mag/MetaPilotMagTask.java`.
- Headless (step 7 only, resumes from `combined/`): `metapilot-core/.../core/runner/PfindMagHeadlessRunner.java`.
- Layout helper: `metapilot-core/.../core/mag/MagResultLayoutUtil.java`.

Notation: `R` = `<result>/pfind_MAG_search/`; `<raw>` = raw file name without extension; `<sfx>` =
MS2 mode suffix (`_HCDFT`, `_HCDIT`, `_CIDFT`, `_CIDIT` or empty); `<smp>` = `R/search_results/<raw>/`.

## Tree

```
<result>/
  logging.log, system_output.log            task logs
  spectra/                                  step 1
    <raw><sfx>.pf2, <raw>.ms2, <raw>.mgf, <raw>.pParse.*, rt.txt, ms2Count.txt
  pfind_MAG_search/                         = R
    search_results/<raw>/                   one folder per raw file
      marker/                               steps 2-3
        pFind.cfg, pFind.spectra, pFind-Filtered.spectra, pFind.protein, pFind.summary
        genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv   (prefix marker_genome_ with Marker4)
        <raw>.fasta                         sample-specific database
      sample_db/                            step 4
        pFind.cfg, pFind.spectra, pFind-Filtered.spectra, pFind.protein, pFind.summary
    combined/                               steps 5-6 (+ id.tsv in step 7)
      psms.tsv
      QuantifiedPeptides.tsv, QuantifiedProteins.tsv, QuantifiedPeaks.tsv, LFQ.bat
      [ms2_intensity.tsv, IsobaricQuantifiedPeptides.tsv]   isobaric only
      id.tsv
      metapilot.psm.parquet
    final_summary.tsv, final_peptides.tsv, razor_proteins.tsv, all_proteins.tsv   step 7a
    metainfo.tsv                            step 7d (GUI only)
    taxonomy_analysis/                      step 7b
    functional_annotation/                  step 7c
    report/                                 step 7d
<rawDir>/ExperimentalDesign.tsv             step 6, written next to the raw files
```

Fixed top-level folders under `R` are `combined`, `taxonomy_analysis`, `functional_annotation`,
`report`, `host`, `search_results`. Any other folder directly under `R` is treated as a legacy
per-sample folder and moved into `search_results/` when the GUI task starts.

## Step 1. Spectra conversion (`<result>/spectra/`)

pParse (or pParse2Plus) converts each raw file.

| File | Content | Used by |
|---|---|---|
| `<raw><sfx>.pf2` | pFind spectrum format | steps 2 and 4 (pFind input) |
| `<raw>.ms2` | MS2 spectra | MGF/RT writer; isobaric reporter ions (step 6) |
| `<raw>.mgf` | MGF from `.ms2` | not used in steps 2-6 (only by the disabled DBReducer branch) |
| `<raw>.pParse.bat`, `.pParse.cfg` / `.pParse.yaml` | converter command and config | record only |
| `rt.txt` | `FileName	Scannum	Rt` (minutes), one row per MS2 scan | step 5 (PSM retention times) |
| `ms2Count.txt` | `FileName	MS2 count` | `final_summary.tsv` (MS/MS column) |

Skip: a raw file is not reconverted when its `.pf2` exists and is non-empty.

## Step 2. Marker search (`<smp>marker/`)

pFind searches each `.pf2` against the MAG database's marker FASTA (`Marker2.fasta` by default,
`Marker4.fasta` in Marker4 mode), decoys `REV_`.

| File | Content | Used by |
|---|---|---|
| `pFind-Filtered.spectra` | pFind PSMs after its FDR | step 3 (genome selection) |
| `pFind.spectra` | all pFind PSMs | result check |
| `pFind.protein`, `pFind.summary` | pFind protein list and run summary | result check / record |
| `pFind.cfg` | search config written by MetaPilot | record |

Skip: `marker/pFind-Filtered.spectra` exists and is non-empty.

## Step 3. Genome selection and sample database (`<smp>marker/`)

`MagGenomeSelectionService` -> `MarkerGenomeSelector` reads the marker PSMs (q <= 0.01, length >= 7)
and selects genomes per sample (privileged by unique peptides, then greedy set cover).
Automatic Focus or Comprehensive applies the knee cut and keeps the 300-genome cap.
PSM coverage, when selected, walks that same ranking until the PSM-coverage target,
the step guard, or `genome_selection_max_genomes`. File prefix `genome` (Marker2) or `marker_genome` (Marker4).

| File | Content | Used by |
|---|---|---|
| `genome_evidence.tsv` | per genome: `Genome	TotalPeptides	UniquePeptides	TotalPSM` (Marker4: per marker category counts, `GenusKey`, `FamilyKey`) | diagnostics |
| `genome_ranked.tsv` | selected genomes first, then the rest: `Rank` (0-based), `Genome`, counts, `Reason`, greedy and knee diagnostics, cumulative coverage, `CutStep`, `InRanking` (`InRanking` 0 is after the greedy cap and is not part of the ranking) | diagnostics |
| `genome_selected.tsv` | `Genome	Reason`, then a blank line and `#TotalMarkerPeptides`, `#TotalMarkerPSMs`, `#SelectedGenomes`, `#SelectionStrategy`, `#PsmCoverageThreshold`, `#StepGuard`, `#MaxGenomes`, `#StopReason`, `#PsmCoverageReached` | builds `<raw>.fasta`; cached on rerun |
| `genome_selection_summary.tsv` | one row in `<raw>/marker/`, and one combined table at `pfind_MAG_search/genome_selection_summary.tsv` with one row per raw file | the coverage curve and the stop reason |
| `<raw>.fasta` | concatenated proteomes (`.faa`) of the selected genomes (+ host FASTA if "append host" is on) | step 4 database |

`Reason` values: `PRIVILEGE_UNIQUE_PEPTIDE_1PCT` (Marker2), `PRIVILEGE_UNIQUE_3TYPES` and
`*_REP` (Marker4), `GREEDY`, `GREEDY_COVERAGE` (PSM coverage cut), `GREEDY_FULL_COVERAGE`, `FALLBACK_MIN1`.

Skip: `<raw>.fasta` exists and is non-empty.

## Step 4. Sample-specific search (`<smp>sample_db/`)

pFind searches the same `.pf2` against `marker/<raw>.fasta`. Same file set as step 2.

| File | Used by |
|---|---|
| `pFind-Filtered.spectra` | step 5 (combined PSMs); sample-folder discovery in step 7 |
| `pFind.spectra` | isobaric quantification (step 6) |

Skip: `sample_db/pFind-Filtered.spectra` exists and is non-empty. (A DBReducer sub-branch,
`sample_db/DBReducer/`, and a separate host search, `R/host/`, exist in code but are disabled.)

## Step 5. Combine (`R/combined/psms.tsv`)

All `sample_db/pFind-Filtered.spectra` merged, with retention times from `rt.txt`.

Header: `File Name	Scan number	Base Sequence	Full Sequence	Peptide Monoisotopic Mass	Scan Retention Time	Precursor Charge	Protein Accession	Score	PEP	Peptide Mass	Mass difference	Missed cleavages`

- `Score` = pFind raw score; **`PEP` holds pFind's q-value**, not a posterior error probability.
- `Protein Accession` = `;`-joined proteins; decoys start with `REV_`.
- Rows: MAG targets, then host PSMs, then decoys.

Used by: FlashLFQ (step 6, identification input) and all of step 7 (PSM counts, peptide scores,
genome and protein scoring, which need the `REV_` decoy rows).

Skip: `psms.tsv` larger than 200 bytes.

## Step 6. Quantification (`R/combined/`)

Label-free: FlashLFQ with match-between-runs and normalisation.

| File | Content | Used by |
|---|---|---|
| `<rawDir>/ExperimentalDesign.tsv` | `FileName	Condition	Biorep	Fraction	Techrep` | FlashLFQ; step 7 sample order (GUI) |
| `QuantifiedPeptides.tsv` | FlashLFQ peptides: `Sequence`, ..., `Intensity_<raw>`, `Detection Type_<raw>` | **step 7 main input** |
| `QuantifiedProteins.tsv` | FlashLFQ protein intensities | not used (MetaPilot does its own protein inference) |
| `QuantifiedPeaks.tsv` | FlashLFQ MS1 peaks | isobaric step |
| `LFQ.bat` | FlashLFQ command | record |

Isobaric labelling adds `ms2_intensity.tsv` (reporter intensities per PSM) and
`IsobaricQuantifiedPeptides.tsv` (`Sequence	Base Sequence	Protein Groups`, `Intensity_<raw>_<label>`...,
`Detection Type_<raw>`...), which then replaces `QuantifiedPeptides.tsv` as the step 7 input.

Skip: `QuantifiedPeptides.tsv` exists and is non-empty.

## Step 7. Post-analysis

Inputs: `combined/QuantifiedPeptides.tsv` + `combined/psms.tsv`. Shared core:
`DdaMagPostAnalysisService.analyze` -> genome scoring -> razor protein inference -> protein scores.
The headless runner (`PfindMagHeadlessRunner`) starts here.

### 7a. Peptide and protein tables

**`combined/id.tsv`**: `Raw file	Unique peptide count	PSM count`, per raw file (distinct base
sequences and distinct scans). Feeds `final_summary.tsv`.

**`final_summary.tsv`**: `Raw file	Experiment	MS/MS	MS/MS Identified	MS/MS Identified [%]	Peptide Sequences Identified`,
one row per raw file plus `Total`. Feeds `report_ID_summary.html`. (Headless: `MS/MS` = PSM count.)

**`final_peptides.tsv`**: one row per quantified peptide.
`Sequence	Base sequence	Length	Missed cleavages	Mass	Razor protein	Mapped proteins	Charges	Score	PEP	`
`Identification type <file>`... `MS/MS count <file>`... `Intensity	Intensity <sample>`...
`Reverse	Potential contaminant	id	Protein group IDs	MS/MS Count`.
- `Score`/`PEP` from the peptide's best PSM (pFind raw score / pFind q-value).
- `Identification type`: `By MS/MS`, `MBR` or `NotDetected`.
- `Razor protein` empty when the peptide has no protein in a confident genome (still listed).
- `Protein group IDs` = `id`s in `razor_proteins.tsv`.

**`razor_proteins.tsv`**: one row per protein that is the razor protein of at least one peptide,
sorted by intensity.
`Protein IDs	Peptide counts (all)	Peptide counts (razor)	Number of PSMs	protein_score_S	protein_p_emp	protein_q_value	Intensity	Intensity <sample>`...
`Reverse	Potential contaminant	id	Peptide IDs	Peptide is razor`.
- Razor rule: candidate proteins are those from **confident genomes (genome q <= 0.01)**; the
  winner has the highest shared intensity.
- Intensity: each peptide's intensity split over its razor-protein candidates.
- `protein_score_S` = sum of evidence weights (best PSM score clipped to [0, 50], times an intensity
  factor); `protein_p_emp` = resampling p against `REV_` decoy evidence (B = 10,000);
  `protein_q_value` = Benjamini-Hochberg. **Reported only, not used as a filter.**
- `id` is referenced by `final_peptides.tsv` and by `Group_ID` in functional annotation.
- GUI: `Peptide IDs` / `Peptide is razor` are empty; headless fills them.

**`all_proteins.tsv`**: same columns; every target protein of any peptide that has a razor protein,
including proteins of non-confident genomes, intensity split evenly over all target proteins.
No downstream reader.

**`combined/metapilot.psm.parquet`**: PSMs in QPX layout (export; failures are non-fatal).

### 7b. Taxonomy (`R/taxonomy_analysis/`)

Genome scoring (`DdaMagGenomeScoringService`): unit = base peptide (I = L) unique to one genome;
weight as for proteins. `genome_score_S` = sum of the genome's top 5 weights; `genome_p_emp` =
resampling p against `REV_` decoy peptides (B = 10,000, minimum about 1e-4); `genome_q_value` =
Benjamini-Hochberg. **Confident genome = q <= 0.01**; only confident genomes receive razor peptides
and intensity, so only they appear in the files below (except `genome_rank_scores.tsv`).

| File | Content | Used by |
|---|---|---|
| `genome_final_evidence.tsv` | `Genome	Intensity	RelativeAbundance	TotalPeptides	RazorPeptides	UniquePeptides	TotalPSM` | diagnostics |
| `genome_final_ranked.tsv` | `Rank` (0-based), `Genome`, the evidence columns, `Reason`, greedy and knee diagnostics, `TotalCoveredPeptides	Coverage	Increment` | `Selection_label` and row order of `Genome.tsv` |
| `genome_final_selected.tsv` | `Genome	Reason` (`PRIVILEGE_UNIQUE_PEPTIDE_1PCT`, `GREEDY_core`, `GREEDY_mid`, `GREEDY_tail`), footer `#PrivilegeGenomes`, `#GreedyGenomes`, `#FocusKnee`, `#ComprehensiveKnee` | checkpoint; variant check |
| `Genome.tsv` | **main genome table**, see below | users |
| `Genome_unique_peptides.tsv` | same columns and rows as `Genome.tsv`; `Intensity` is summed from peptides unique to that genome only | users |
| `Genome.html` | genome table with lineage and MGnify links | report (GUI copies it to `report/reports/`) |
| `Taxa.tsv` | `Name	Score	Rank	domain	kingdom	phylum	class	order	family	genus	species	Intensity <sample>`...; one row per taxon per rank; `Score` = sum of `genome_score_S` | users |
| `Taxa_report.tsv` | as `Taxa.tsv` without `Score` | `report_taxonomy_summary.html` |
| `genome_rank_scores.tsv` | `Rank	Genome	genome_score_S	n_peptides	genome_p_emp	genome_q_value`, **all scored genomes including non-confident**, by S | genome-test audit |
| `genome_rank_scores_runinfo.txt` | genome test diagnostics: peptide counts, sharing, parameters (top_k 5, B 10000, cap 50), reversed-genome calibration table (GUI only) | genome-test audit |
| `taxon_abundance.parquet` | opt-in (`--export-metaatlas`, headless): `taxon_name, taxon_rank, taxonomy_source, sample_id, abundance, abundance_type, evidence_score` | MetaAtlas |

`Genome.tsv` header:
`Genome	genome_score_S	genome_p_emp	genome_q_value	Selection_label	PSM count	Peptide count	Razor peptide count	Unique peptide count	Peptide I/L ratio	Razor peptide I/L ratio	Unique peptide I/L ratio	GTDB_Name	GTDB_Rank	GTDB_<rank>`... `NCBI_Name	NCBI_Rank	NCBI_<rank>`... `Intensity <sample>`...
- One row per confident genome, in ranked-selection order; `Selection_label` = `core`, `mid` or `tail`.
- `Peptide count` / `Razor` / `Unique`: peptides mapping to the genome / assigned to it as razor /
  mapping only to it.
- GTDB lineage from the MAG database taxa file; NCBI from the bundled GTDB->NCBI mapping.
- Intensity = sum of the genome's razor-protein intensities (taxonomic abundance).
- `Genome_unique_peptides.tsv` repeats this table with Intensity = sum of peptides that map only to that genome (no host proteins, no other MAG genomes). Shared peptides contribute 0.

### 7c. Functional annotation (`R/functional_annotation/`)

Razor proteins (from `razor_proteins.tsv`) matched to their genome's eggNOG annotation; intensities
are razor-protein intensities, split equally when a protein has several IDs in one category.

| File | Key columns |
|---|---|
| `Annotated_proteins.tsv` | `Group_ID` (= `razor_proteins.tsv` `id`), `Name`, `Peptide count`, `PSM count`, `Intensity <sample>`..., then eggNOG fields (description, GO, EC, KEGG KO/pathway/module/reaction, BRITE, CAZy, PFAMs, COG, NOG) |
| `Gene_Ontology_id_report.tsv` | `Gene_Ontology_id	Gene_Ontology_name	Gene_Ontology_namespace	Intensity <sample>`... |
| `EC_id_report.tsv` | `EC_id	EC_description	EC_reaction	EC_synonyms	Intensity`... |
| `KEGG_ko_report.tsv` | `KEGG_ko	KEGG_ko_description	Intensity`... |
| `BRITE_report.tsv` | `BRITE	Intensity`... |
| `COG_accession_report.tsv`, `NOG_accession_report.tsv` | accession, category, name, intensities |
| `COG_category_report.tsv`, `NOG_category_report.tsv` | category, intensities |
| `functions.html` | browsable protein-function table |
| `function_abundance.parquet` | opt-in MetaAtlas export |

Used by `report_function_summary.html`. Headless skips this step when `func_def.db` is not configured.

### 7d. Report (`R/report/`, `R/metainfo.tsv`)

- `metainfo.tsv` (GUI only): SDRF-like sample sheet, `source name	comment[data file]	comment[fraction identifier]	comment[technical replicate]	characteristics[...]	factor value[...]`; used by the HTML reports for grouping.
- `report/index.html` + `assets/` + `manifest.*`: report shell.
- `report/reports/`: `report_ID_summary.html` (from `final_summary.tsv`),
  `report_peptides_summary.html` (`final_peptides.tsv`), `report_proteinGroups_summary.html`
  (`razor_proteins.tsv`), `report_taxonomy_summary.html` (`Taxa_report.tsv`),
  `report_function_summary.html` (`functional_annotation/`), GUI copies of `Genome.html` and
  `functions.html`.

GUI rerun: if all step 7 outputs are newer than `psms.tsv` and `QuantifiedPeptides.tsv`, only the
reports are refreshed.

## Which file answers which question

| Question | File |
|---|---|
| Which genomes were candidates in a sample? | `search_results/<raw>/marker/genome_selected.tsv` |
| Identification counts per raw file | `final_summary.tsv` |
| Peptide identifications and intensities | `final_peptides.tsv` |
| Protein abundances (razor) | `razor_proteins.tsv` |
| Genomes present, their evidence and abundance | `taxonomy_analysis/Genome.tsv` |
| Genome test for all candidates, incl. rejected | `taxonomy_analysis/genome_rank_scores.tsv` |
| Abundance per taxon at each rank | `taxonomy_analysis/Taxa.tsv` |
| Function abundances | `functional_annotation/*_report.tsv` |
| Raw PSMs with decoys | `combined/psms.tsv` |

## FragPipe variant

Root `<result>/fragpipe_MAG_search/`, same `search_results/`, `combined/`, step 7 files and headers
(same writers). Differences: search and quantification by FragPipe (input
`combined/combined_modified_peptide.tsv` and FragPipe `psm.tsv`; score = hyperscore); decoy prefix
for genome and protein scoring is `shuf_`; MS/MS counts from FragPipe spectral counts;
`Peptide IDs` columns filled; no `genome_rank_scores_runinfo.txt`; the GUI generates the HTML
report twice (headless runner, then GUI with `metainfo.tsv`).

## Known quirks

- `psms.tsv`: `PEP` column holds pFind q-values; decoy rows can be written twice.
- GUI `razor_proteins.tsv` / `all_proteins.tsv`: `Peptide IDs` and `Peptide is razor` are empty.
- Isobaric mode writes `spectra/ms2count.txt`, which on Windows is the same file as `ms2Count.txt`
  (different header).
- Headless pFind writes no `metainfo.tsv` and no `genome_rank_scores_runinfo.txt`; its reports have no
  sample metadata.

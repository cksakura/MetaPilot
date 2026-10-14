# DIA Marker Deep result files

What the DIA Marker workflow writes in its two deep modes, step by step, and what each file is for.
Both modes use DIA-NN, start from a Marker4 search, and walk the ranked genome list in rounds of
5 genomes until new peptides plateau; they differ in whether the rounds run on all DIA files at
once (**Deep together**) or on each file separately (**Deep per file**). The final quantifying
search and the MAG post-analysis are the same. For the DDA workflow see `dda-mag-result-files.md`.

Code entry points:

- Deep together: headless `metapilot-core/.../core/runner/DiaMarkerDeepTogetherHeadlessRunner.java`;
  GUI `metapilot-windows-gui/.../task/dia/DiaMarkerDeepTogetherTask.java` (runs its own
  implementation in `DiaMarkerDeepPerFileTask` / `DiaMarkerCoreTask`).
- Deep per file: headless `metapilot-core/.../core/runner/DiaMarkerDeepPerFileHeadlessRunner.java`;
  GUI `metapilot-windows-gui/.../task/dia/DiaMarkerDeepPerFileTask.java` (calls the core runner
  without post-analysis, then its own post-analysis).
- Round loop and stop rule: `metapilot-core/.../task/dia/DeepIterativeCommonSupport.java`.
- Post-analysis: headless `DiaMagPostAnalysisTsvExporter`; GUI `MetaDiaTask.exportReport` ->
  `DiaMarkerCoreTask.runFullDownstreamFromNewScoring`.
- Rerun post-analysis only: CLI `diann-mag-report` (`DiaMagReportHeadlessRunner`).

Notation: `R` = `<result>/marker_deep/` (both modes); `<raw>` = DIA file name without extension
(`.raw.dia`, `.raw.d`, `.mzml`, `.dia`, `.raw`, `.d` stripped); `N` = round number.

## Tree

```
<result>/
  dia_marker_analysis_parameters.json        GUI only
  marker_deep/                               = R
    marker/
      --- Deep together ---
      marker.parquet                         step 1
      genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv   step 2
      round1/ ... roundN/                    step 3: combined.fasta, combined.speclib, roundN.parquet
      iterative_rounds_info.tsv              step 3
      --- Deep per file: the same set, once per DIA file ---
      <raw>/marker.parquet, genome_*.tsv, round1/ ... roundN/, iterative_rounds_info.tsv
    representative_groups.tsv                Deep per file, clustering on
    representative_pairwise_similarity.tsv   Deep per file, clustering on
    representative_clustering_summary.tsv    Deep per file, clustering on
    combined/                                step 4: combined.fasta, combined.speclib, combined.parquet,
                                                     combined.pr_matrix.tsv, combined.pg_matrix.tsv
    host/                                    optional host search
    id.tsv, final_peptides.tsv, razor_proteins.tsv, all_proteins.tsv, final_summary.tsv   step 5a
    metainfo.tsv                             GUI only
    taxonomy_analysis/                       step 5b
    functional_annotation/                   step 5c
    report/                                  step 5d
<MAG DB current folder>/                     step 0, shared across analyses
  Marker4.fasta, genomes-taxa.tsv, libraries/Marker4/Marker4.speclib
```

Every DIA-NN call leaves `<output>.cfg` (the DIA-NN command line) next to its output, e.g.
`marker.parquet.cfg`, `combined.speclib.cfg`. DIA-NN may also leave its own side files (stats, logs,
first-pass parquet, gene matrices); MetaPilot does not read them.

## Step 0. Setup (outside `R`)

- In the MAG database folder: `Marker4.fasta` (marker proteins), `genomes-taxa.tsv` (genome ->
  lineage, no header), `libraries/Marker4/Marker4.speclib` (DIA-NN predicted library), plus
  `Marker4.parquet` / `Marker4_info.txt`. Built once, reused by every analysis.
- GUI only: `.raw` / `.wiff` converted to `<file>.dia` **next to the raw file**; `R/report/` template
  copied at task start; parameters saved to `<result>/dia_marker_analysis_parameters.json`.

## Step 1. Marker search

DIA-NN searches `Marker4.speclib`, no quantification (`--no-prot-inf --no-maxlfq --no-quant-files`,
q 0.01).

| Mode | File |
|---|---|
| Deep together | `R/marker/marker.parquet` (all files in one search) |
| Deep per file | `R/marker/<raw>/marker.parquet` (one search per file) |

DIA-NN report columns MetaPilot reads from any `.parquet` here and below: `Run`, `Protein.Ids`,
`Modified.Sequence`, `Stripped.Sequence`, `Precursor.Mz`, `Precursor.Charge`,
`Precursor.Normalised`, `Precursor.Quantity`, `Ms1.Normalised`, `Ms1.Area`, `Q.Value`, `PEP`,
`Global.Q.Value`, `PG.Q.Value`, `Global.PG.Q.Value`, `Decoy`, `CScore`. Precursors count as
identified at `Global.Q.Value < 0.01`.

Used by: step 2 (genome selection); per file with clustering, also step 2b and the step 4 peptide map.
Skip: the parquet exists.

## Step 2. Genome ranking (same folder as step 1)

`MagGenomeSelectionService` -> `MarkerGenomeSelector`, Marker4 categories, comprehensive strategy,
at most 300 selected.

| File | Content | Used by |
|---|---|---|
| `genome_evidence.tsv` | `Genome	Ribosomal	Replication	Elongation	Synthetase	UniqueRibosomal	UniqueReplication	UniqueElongation	UniqueSynthetase	TotalMarkers	TotalUniqueMarkers	TotalPSM	GenusKey	FamilyKey` | diagnostics, resume check |
| `genome_ranked.tsv` | `Rank` (0-based), `Genome`, `GenusKey`, `FamilyKey`, the marker counts, `nCat	nCatUnique	Reason`, greedy and knee diagnostics, cumulative peptide and PSM coverage | **drives the rounds** (column `Genome`, in order) |
| `genome_selected.tsv` | `Genome	Reason`, blank line, `#TotalMarkerPeptides`, `#TotalMarkerPSMs`, `#SelectedGenomes` | resume check |

The rounds walk `genome_ranked.tsv` from the top, 5 genomes per round (`genomesPerRound`), past the
300 selected if needed. Skip: all three files exist.

### Step 2b. Representative clustering (Deep per file only, optional)

Groups DIA files whose marker peptide sets are similar (Jaccard / containment thresholds from the
preset CONSERVATIVE, BALANCED, AGGRESSIVE or CUSTOM); only one representative per group runs the
rounds. Files in `R/`:

- `representative_groups.tsv`: one row per DIA file: `RawFile	RawName	GroupId	RepresentativeRawFile	RepresentativeRawName	Role	MarkerPeptides	RepresentativeMarkerPeptides	ClusterSize	PeptideJaccard	PeptideContainment	PeptideDistance	GenomeTopNOverlapDiagnostic	BestMatchRawFile	BestMatchRawName	BestPeptideJaccard	BestPeptideContainment	BestGenomeTopNOverlapDiagnostic	Reason`
  (`Role` = `representative` / `member`).
- `representative_pairwise_similarity.tsv`: one row per compared pair (peptide set sizes, overlap,
  Jaccard, containment, pass flag).
- `representative_clustering_summary.tsv`: `Metric	Value` (thresholds, algorithm, counts of files,
  pairs, clusters, representatives, singletons).

Reports only; rewritten every run. Dry run stops the workflow after writing them. Non-representative
files keep only their step 1-2 files.

Per file also has a gate: every file entering the rounds must have its genome files; the GUI asks
whether to continue without the failed files, headless stops.

## Step 3. Iterative rounds

Per round `N`, in `R/marker/roundN/` (together) or `R/marker/<raw>/roundN/` (per file):

| File | Content | Used by |
|---|---|---|
| `combined.fasta` | round 1: full proteomes of the next 5 ranked genomes; round 2+: the next 5 genomes plus the peptides identified in round N-1 (one entry per peptide) | DIA-NN library prediction |
| `combined.speclib` | DIA-NN predicted library (trypsin, 1 missed cleavage, length 7-30, charge 2-3) | the round search |
| `roundN.parquet` | DIA-NN no-quant search of that library | next round's FASTA; **peptide -> protein map** for step 4 and post-analysis |

Stop rule (on new unique peptides at `Global.Q.Value < 0.01`): end of ranked list; relative gain
< 1% for 2 consecutive rounds; from round 8, knee/plateau rule at 1.5%; or 100 rounds.

`iterative_rounds_info.tsv` (one per mode unit, next to the round folders): comment lines with the
stop rule (and its fitted parameters after 8+ rounds), then
`Round	Peptides_in_round	New_peptides	Cumulative_peptides	Stop_reason` (`Stop_reason` on the last
row only). Report only.

Skip: a round's search is skipped when `roundN.parquet` exists; library prediction when
`combined.speclib` exists.

## Step 4. Combined quantifying search (`R/combined/`)

| File | Content | Used by |
|---|---|---|
| `combined.fasta` | one entry per peptide identified in the last rounds (per file: union over all files; with clustering also all `marker.parquet` peptides); header = one protein of that peptide | library prediction |
| `combined.speclib` | predicted library | the search |
| `combined.parquet` | DIA-NN report, **all DIA files**, with decoys (`--report-decoys`), match-between-runs, `--relaxed-prot-inf` | **post-analysis**: precursors, scores, decoys |
| `combined.pr_matrix.tsv` | precursor x run intensities (`Stripped.Sequence`, `Modified.Sequence`, `Precursor.Charge`, `Protein.Group`, one column per run) | **post-analysis**: peptide intensities |
| `combined.pg_matrix.tsv` | DIA-NN protein-group matrix | not read (existence check only) |

Because `combined.fasta` headers carry a single protein, the post-analysis rebuilds each peptide's
full protein list from the round parquets (`DiaDeepProteinMappingService`) before scoring.

Skip: together headless needs `combined.parquet`, `pr_matrix` and `pg_matrix`; GUI and per file
check `combined.parquet` (per file then fails if the matrices are missing).

Optional host search (`appendHostDb`): `R/host/` with `host.speclib` (FASTA input only),
`host.parquet`, `host.pr_matrix.tsv`, `host.pg_matrix.tsv`. Not read by the post-analysis.

## Step 5. MAG post-analysis

Inputs: `combined/combined.parquet`, `combined/combined.pr_matrix.tsv`, the round parquets (protein
remap), the MAG database, optionally `func_def.db` and `Resources/gtdb_ncbi`.

Scoring (`DiaMagScoringService`), decoys = proteins prefixed `shuf_`:

- Peptide evidence `e` = min(6, CScore) (or min(6, -log10 PEP)), weighted by
  0.5 + 0.25 x normalised log2 precursor quantity + 0.25 x normalised log2 MS1 area.
- `genome_score_S` = sum of the top 5 weights among peptides unique to the genome;
  `genome_p_emp` = resampling p against decoy peptides (B = 10,000); `genome_q_value` =
  Benjamini-Hochberg. **Confident genome = q <= 0.01**; only confident genomes receive intensity.
- `protein_score_S` / `protein_p_emp` / `protein_q_value`: same scheme over all of a protein's
  peptides. Reported, not used as a filter.
- Precursors kept at `Global.Q.Value < 0.01`; proteins of genomes with p = 1 dropped; razor protein
  = protein with the highest S.

### 5a. Peptide and protein tables (`R/`)

**`id.tsv`**: `Raw file	Unique peptide count	PSM count`, per run. (GUI bug: the peptide count
counts peptides *not* detected in that run; headless counts detected ones.)

**`final_peptides.tsv`**: one row per modified sequence.
`Sequence	Base sequence	Length	Missed cleavages	Mass	Razor protein	Mapped proteins	Charges	Qvalue	PEP	`
`Identification type <run>`... `MS/MS count <run>`... `Intensity	Intensity <run>`...
`Reverse	Potential contaminant	id	Protein group IDs	MS/MS Count`.
- `Intensity <run>` from `pr_matrix`, summed over charges. `Qvalue` / `PEP` = best over precursors.
- `Identification type`: `MS/MS`, `MBR` or `unmatched`.
- Feeds `report_peptides_summary.html`.

**`razor_proteins.tsv`**: one row per protein that is razor for at least one peptide.
`Protein IDs	Peptide counts (all)	Peptide counts (razor)	Number of PSMs	protein_score_S	protein_p_emp	protein_q_value	Intensity	Intensity <run>`...
`Reverse	Potential contaminant	id	Peptide IDs	Peptide is razor`.
- Peptide intensity split equally over its razor proteins; sorted by S then intensity.
- `id` is referenced by `final_peptides.tsv` and by `Group_ID` in functional annotation.
- `Peptide is razor` = `true`/`false` (GUI) or `1`/`0` (headless).
- Feeds `report_proteinGroups_summary.html` and functional annotation.

**`all_proteins.tsv`**: same columns, every mapped protein, peptide intensity split over all its
proteins. No downstream reader.

**`final_summary.tsv`**: `Raw file	Experiment	MS/MS	MS/MS Identified	MS/MS Identified [%]	Peptide Sequences Identified`,
per run plus `Total`. For DIA, `MS/MS` = precursor rows (not spectra), `[%]` is always 1.0.
Feeds `report_ID_summary.html`.

**`metainfo.tsv`** (GUI only): SDRF-like sample sheet, `source name	comment[data file]	comment[fraction identifier]	comment[technical replicate]	characteristics[...]	factor value[...]`;
used by the GUI reports for sample grouping.

### 5b. Taxonomy (`R/taxonomy_analysis/`)

Same writers and headers as DDA MAG (`MagTaxonomyAnalysisService`):

| File | Content | Used by |
|---|---|---|
| `genome_final_evidence.tsv` | `Genome	Intensity	RelativeAbundance	TotalPeptides	RazorPeptides	UniquePeptides	TotalPSM` | diagnostics |
| `genome_final_ranked.tsv` | `Rank` (0-based), evidence columns, `Reason` (`PRIVILEGE_UNIQUE_PEPTIDE_1PCT`, `GREEDY`, empty), greedy and knee diagnostics | `Selection_label` and row order of `Genome.tsv` |
| `genome_final_selected.tsv` | `Genome	Reason` (`GREEDY_core` / `GREEDY_mid` / `GREEDY_tail`), footer counts and knees | diagnostics |
| `Genome.tsv` | **main genome table**: `Genome	genome_score_S	genome_p_emp	genome_q_value	Selection_label	PSM count	Peptide count	Razor peptide count	Unique peptide count`, three I/L ratios, `GTDB_Name	GTDB_Rank	GTDB_<rank>`..., `NCBI_Name	NCBI_Rank	NCBI_<rank>`..., `Intensity <run>`... ; one row per confident genome | users |
| `Genome.html` | genome table with lineage | report (GUI copies to `report/reports/`) |
| `Taxa.tsv` | `Name	Score	Rank	domain	kingdom	phylum	class	order	family	genus	species	Intensity <run>`...; `Score` = sum of genome S | users |
| `Taxa_report.tsv` | as `Taxa.tsv` without `Score` | `report_taxonomy_summary.html` |
| `genome_rank_scores.tsv` | `Rank	Genome	genome_score_S	n_peptides	genome_p_emp	genome_q_value`, **all scored genomes incl. non-confident** | genome-test audit |
| `taxon_abundance.parquet` | opt-in, headless `--export-metaatlas` | MetaAtlas |

`Selection_label`: `core` (privileged or greedy within the focus knee), `mid` (within the
comprehensive knee), `tail`. Genome intensity = sum of its razor-protein intensities.

### 5c. Functional annotation (`R/functional_annotation/`)

Razor proteins of MAG genomes matched to their genome's eggNOG annotation; razor intensities split
equally across multiple IDs in a category. Required in the GUI; headless skips it when
`func_def.db` is not configured.

- `Annotated_proteins.tsv`: `Group_ID` (= `razor_proteins.tsv` `id`), `Name`, `Peptide count`,
  `PSM count`, `Intensity <run>`..., then eggNOG fields (description, GO, EC, KEGG, BRITE, CAZy,
  PFAMs, COG, NOG).
- `Gene_Ontology_id_report.tsv`, `EC_id_report.tsv`, `KEGG_ko_report.tsv`, `BRITE_report.tsv`,
  `COG_accession_report.tsv`, `NOG_accession_report.tsv`, `COG_category_report.tsv`,
  `NOG_category_report.tsv`: one row per function ID or category, then `Intensity <run>`...
- `functions.html`; `function_abundance.parquet` (opt-in, headless).

Feeds `report_function_summary.html`.

### 5d. Report (`R/report/`)

`report/` holds the template; `report/reports/` the outputs: `report_ID_summary.html`
(`final_summary.tsv`), `report_peptides_summary.html` (`final_peptides.tsv`),
`report_proteinGroups_summary.html` (`razor_proteins.tsv`), `report_taxonomy_summary.html`
(`Taxa_report.tsv`), `report_function_summary.html` (`functional_annotation/`), `manifest.json`;
GUI also copies `Genome.html` and `functions.html`. GUI skips existing non-empty reports; headless
regenerates them.

## Which file answers which question

| Question | File |
|---|---|
| Which genomes the marker search ranked (per file in per-file mode) | `marker/[<raw>/]genome_ranked.tsv` |
| How many rounds ran and why they stopped | `marker/[<raw>/]iterative_rounds_info.tsv` |
| Which files were grouped (per file, clustering) | `representative_groups.tsv` |
| Final DIA-NN precursors with decoys | `combined/combined.parquet` |
| Identification counts per run | `final_summary.tsv` |
| Peptide identifications and intensities | `final_peptides.tsv` |
| Protein abundances (razor) | `razor_proteins.tsv` |
| Genomes present, their evidence and abundance | `taxonomy_analysis/Genome.tsv` |
| Genome test for all candidates, incl. rejected | `taxonomy_analysis/genome_rank_scores.tsv` |
| Abundance per taxon at each rank | `taxonomy_analysis/Taxa.tsv` |
| Function abundances | `functional_annotation/*_report.tsv` |

## GUI vs headless

- Deep together GUI runs its own implementation (`DiaMarkerCoreTask`), not the headless runner;
  Deep per file GUI calls the core runner for steps 1-4 and then its own post-analysis.
- GUI only: raw conversion to `.dia`, saved parameters JSON, `metainfo.tsv`, report copies of
  `Genome.html` / `functions.html`, dialog on partial marker results.
- Headless only: MetaAtlas parquet exports; functional annotation optional.
- `id.tsv` peptide counts differ (GUI counts undetected peptides).
- Moving runs between machines: `diann-mag-peptide-map` writes the round-parquet peptide -> protein
  map as one TSV (`Peptide	Proteins`), which `diann-mag-report --peptide-protein-map` can use in
  place of the round folders.

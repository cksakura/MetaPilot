# MetaPilot User Manual

MetaPilot is a metaproteomics workflow platform for DDA and DIA datasets. It combines established search engines with MAG-aware genome selection, refined search-space construction, taxonomy export and functional annotation.

This page targets the upcoming `1.2.0` release, pending finalization of the open-search strategy.

## Applications

| Module | Purpose | Current role |
| --- | --- | --- |
| `metapilot-core` | Shared workflow, parsing, MAG, taxonomy and reporting logic. | Main reusable implementation layer |
| `metapilot-windows-gui` | JavaFX desktop application. | Broadest workflow UI |
| `metapilot-cli` | Headless command-line runner. | Linux/HPC runner for selected workflows |
| `metapilot-web` | Spring Boot web wrapper. | Server/API wrapper that launches a JAR subprocess |

## Installation

Install MetaPilot from the release ZIP package. Public releases are distributed as ready-to-run application artifacts; building from source is not part of the supported user installation path.

Typical release contents:

```text
MetaPilot-1.2.0/
  metapilot-cli-1.2.0-shaded.jar
  metapilot
  config/
    meta_resources.params.template
    MAGDB_1_0.json
  docs/
  examples/
```

Extract the ZIP to a local folder on Windows or to shared storage on Linux/HPC. Use the included executable, launcher script or CLI JAR for the workflows described below.

## Required Tools

Configure external tools through the GUI resource panel or `meta_resources.params`.

| Tool | Used by |
| --- | --- |
| pFind | DDA IGC and DDA MAG pFind workflows |
| MaxQuant | DDA IGC MaxQuant workflow, Windows-oriented |
| FragPipe/MSFragger | DDA IGC and DDA MAG FragPipe workflows |
| Sage | DDA IGC and DDA MAG Sage workflows |
| DIA-NN | DIA Marker, DIA Seed, DIA Target and DIA Real-Time workflows |
| DeepDetect | Marker4/peptide detectability prediction and reduced DIA library generation |
| Python | Report generation and helper scripts |
| Functional database resources | Functional annotation export |

Example Windows resource file values:

```properties
PFIND=C:\Users\...\Resources\DBReducer\pFind3\bin\pFind.exe
MAXQUANT=C:\Users\...\Resources\mq_bin\MaxQuantCmd.exe
FRAGPIPE=C:\FragPipe\FragPipe-24.0\bin\fragpipe.bat
SAGE=C:\Users\...\Resources\sage-v0.14.7-x86_64-pc-windows-msvc\sage.exe
DIANN=C:\DIA-NN\2.3.0\diann.exe
DEEPDETECT=C:\Users\...\Resources\DeepDetect\DeepDetect.exe
```

## MAG Databases

MAG databases are catalogue entries in `MAGDB_1_0.json` and local database folders. A typical MAG database folder contains:

- `original_db/`: per-genome protein FASTA files.
- `genomes-taxa.tsv` or catalogue taxonomy metadata.
- `Marker4.fasta`: extracted marker proteins.
- `Marker4.predicted.tsv`: DeepDetect peptide detectability output.
- `Marker4_highpep.fasta`: reduced high-detectability peptide FASTA when Marker4 is large.
- `libraries/Marker4/Marker4.speclib`: Marker4 DIA-NN spectral library for DIA Marker.
- `predicted/`: per-genome predicted peptide TSV files used by Seed mapping and peptide-based rounds.

The current local catalogue configuration includes MGnify entries such as `human-gut` v2.0.2 and `mouse-gut` v1.0.

## Workflow Families

### DDA IGC

Integrated gene catalogue workflows:

- `DDA_IGC_MAXQUANT`
- `DDA_IGC_PFIND`
- `DDA_IGC_FRAGPIPE`
- `DDA_IGC_SAGE`

These workflows use a conventional protein database rather than MAG genome selection as the main reference strategy.

### DDA MAG

MAG-aware DDA workflows:

- `DDA_MAG_PFIND`
- `DDA_MAG_FRAGPIPE`
- `DDA_MAG_SAGE`

The workflow performs an initial marker/high-abundance search, ranks genomes, builds a sample-specific FASTA from selected genomes, then performs final identification, quantification and report export.

Marker4 is the current standard marker strategy. Marker2 is retained for DDA continuity with older MetaLab-style marker evidence.

### DIA Marker

`DIA_MARKER` starts from a Marker4 spectral library generated from the selected MAG database.

Supported modes:

- Core / Focus (`FOCUS`)
- Comprehensive (`COMPREHENSIVE`)
- Deep together (`DEEP`)
- Deep per file (`DEEP_PER_FILE`)

The workflow searches Marker4 evidence, ranks genomes, builds `combined.fasta` and `combined.speclib` from selected genomes, runs final DIA-NN quantification and exports unified taxonomy tables.

### DIA Seed

`DIA_SEED` starts from user-selected DIA libraries rather than Marker4.

Supported GUI modes:

- Direct (`DIRECT`)
- Core / Focus (`FOCUS`)
- Comprehensive (`COMPREHENSIVE`)
- Deep together (`DEEP`)

Seed Deep per file is part of the canonical mode model but is currently blocked by GUI validation as unsupported.

Seed iterative modes map initial library peptides back to MAG proteins through `predicted/*.predicted.tsv`, then use the same genome-ranking and refined-library construction backbone as DIA Marker.

### DIA Target

`DIA_TARGET` is a targeted follow-up workflow for selected taxa, functions or proteins. It is separate from the shared Marker/Seed MAG genome-selection backbone.

### DIA Real-Time

`DIA_REALTIME` monitors incoming DIA files and runs DIA analysis as files arrive. It is a GUI workflow and should not be treated as part of the shared Marker/Seed cleanup path.

## Genome Selection Summary

The current MAG selection strategy is documented in detail in [Workflow Design](workflows.md).

In brief:

1. Extract peptide-to-protein-to-genome evidence from the initial search.
2. For Marker4, track four categories: ribosomal, replication/transcription, elongation and aminoacyl-tRNA synthetase.
3. Apply the Marker4 hard filter when available: ribosomal + replication/transcription + elongation evidence.
4. Add privileged genomes with at least three unique marker categories and sufficient new peptide contribution.
5. Optionally add taxonomy representatives.
6. Greedily order remaining genomes by newly covered marker peptides.
7. Choose when to stop adding genomes. Automatic, the default, detects the tail of the coverage curve. Focus uses the median of kneedle, curvature and slope-change knees; Comprehensive uses the later of chord-distance and piecewise knees. Focus favours depth on the main genomes; Comprehensive adds more low-abundance genomes.
8. pFind, FragPipe and open search can instead stop at a PSM coverage: keep genomes until a chosen share of the marker PSMs is explained (default 0.95). A step guard, adopted from MetaLab-MAG (J Proteome Res 2023), stops when further genomes add almost nothing. Set the guard to 0 to reach the target exactly. This is the coverage and step-guard rule adopted from MetaLab-MAG. It is not the MetaLab-MAG workflow, which ranks by plain greedy PSM cover. Higher targets add more genomes, and a larger database can lose low-abundance proteins from the other genomes.
9. After a pFind run, read the curve as PSM coverage against genomes in rank order. The plot marks where selection stopped, where Focus and Comprehensive would have stopped, and the target when PSM coverage is selected. The run summary gives genomes selected, coverage reached and the stop reason. If the target was missed, the stop reason says why.

## Outputs

Common output artifacts include:

- `genome_evidence.tsv`
- `genome_ranked.tsv`
- `genome_selected.tsv`
- `genome_selection_summary.tsv` (DDA MAG pFind: one row per raw file)
- `combined/combined.parquet`
- `combined/combined.pr_matrix.tsv`
- `combined/combined.pg_matrix.tsv`
- `taxonomy_analysis/Genome.tsv`
- `taxonomy_analysis/Taxa_report.tsv`
- `functional_analysis/annotated_proteins.tsv`
- `functional_analysis/functions.tsv`

GUI workflows may also generate HTML reports.

## CLI Use

Use the CLI for currently supported headless workflows:

```bash
java -jar metapilot-cli-1.2.0-shaded.jar init
java -jar metapilot-cli-1.2.0-shaded.jar params-validate --params study.json
java -jar metapilot-cli-1.2.0-shaded.jar run --params study.json
```

See [CLI Reference](cli-reference.md) for command details and current support limits.

## Recommended Workflow Choice

| Goal | Recommended workflow |
| --- | --- |
| Conventional DDA catalogue search | DDA IGC |
| MAG-aware DDA search | DDA MAG pFind, FragPipe or Sage |
| DDA-independent DIA from a MAG database | DIA Marker |
| DIA using existing user libraries | DIA Seed |
| Maximum DIA discovery depth | DIA Marker Deep together or Deep per file |
| Targeted taxa/functions/proteins | DIA Target |
| Instrument-folder monitoring | DIA Real-Time |
| Linux/HPC headless execution today | DIA Marker or DDA MAG FragPipe |

## Troubleshooting

- Validate all tool paths before running.
- Check `logging.log` in the result directory for headless runs.
- For DIA Marker, confirm `Marker4.fasta`, `genomes-taxa.tsv` and `libraries/Marker4/Marker4.speclib` can be generated or found.
- For Seed workflows, confirm selected libraries and MAG `predicted/` peptide mapping files exist.
- For FragPipe headless runs, confirm the installed FragPipe workflow files are present.
- For Deep per file, budget memory conservatively; the shared support code can size worker pools using `METAPILOT_DEEP_PER_FILE_POOL` and `METAPILOT_DEEP_PER_FILE_RAM_PER_TASK_GB`.

# MetaPilot

MetaPilot is a metaproteomics workflow platform for DDA and DIA datasets. It helps users move from raw mass spectrometry files to MAG-aware genome selection, refined database or library searches, taxonomy reports, and functional annotation outputs.

This repository hosts the public MetaPilot release. The README is intended to be self-contained for release users, so it does not assume access to a separate documentation folder.

MetaPilot supports two main ways of working:

- A Windows desktop GUI for interactive setup and the broadest workflow coverage.
- A command-line runner for Linux and HPC execution of supported headless workflows.

This README is written for users who want to install MetaPilot, choose the right workflow, and understand the result folders.

## What MetaPilot Does

MetaPilot is designed for studies where the search space can be improved using metagenome-assembled genomes (MAGs). In MAG-aware workflows, MetaPilot usually follows this pattern:

1. Run a small initial search using marker proteins, high-abundance proteins, or existing DIA libraries.
2. Use the initial peptide evidence to rank candidate genomes.
3. Select a focused or comprehensive set of genomes.
4. Build a refined sample-specific FASTA or DIA spectral library.
5. Run the final identification and quantification search.
6. Export taxonomy, genome-level, peptide/protein, and functional annotation reports.

Conventional DDA workflows are also available when users want to search an integrated gene catalogue or another standard protein database.

## Installation And Setup

### Release Package

Download the latest MetaPilot release package from this repository's release page. A release package may include:

- A Windows desktop application or runnable GUI JAR.
- A Linux/HPC command-line JAR.
- Example configuration files.
- Optional MAG database registry templates.
- Optional bundled report templates and helper resources.

You do not need any additional repository to run a public release. Building from source is not part of the supported user installation path.

### Requirements

MetaPilot requires Java 17 or newer.

External tools are required depending on the workflow:

| Tool | Used for |
| --- | --- |
| DIA-NN | DIA Marker, DIA Seed, DIA Target, DIA Real-Time |
| FragPipe / MSFragger | DDA FragPipe workflows |
| pFind | DDA pFind workflows |
| MaxQuant | DDA MaxQuant workflow, mainly Windows-oriented |
| Sage | DDA Sage workflows |
| DeepDetect | Marker4 peptide detectability prediction |
| Python | Report generation and helper scripts |
| Functional annotation resources | Functional result export |

### Install Java

Install Java 17 or newer and confirm it is available:

```bash
java -version
```

The current formal release version is `1.0.0`.

### Windows Desktop Setup

Use the Windows GUI when you want the full workflow set and an interactive setup experience. Start the application from the executable or GUI JAR included in the release package.

In the GUI, configure the external tool paths and MAG database locations through the resource/settings panels. Typical resource values include paths to `pFind.exe`, `MaxQuantCmd.exe`, `fragpipe.bat`, `sage.exe`, `diann.exe`, `DeepDetect.exe`, Python, and the functional annotation database.

Before running a study:

1. Confirm all required tool paths are set.
2. Register or select a MAG database when using MAG-aware workflows.
3. Select raw files or existing DIA libraries.
4. Choose a workflow mode.
5. Choose an output folder with enough disk space.

### Linux And HPC Setup

Use the CLI JAR for supported headless workflows, especially on servers and clusters. Copy the release JAR to a stable location such as `/opt/metapilot/` or a shared project directory visible to compute nodes.

Initialize configuration:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar init
```

Set tool paths:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar resources-set \
  --set DIANN=/opt/tools/diann/diann-linux \
  --set FRAGPIPE=/opt/tools/fragpipe-24.0/bin/fragpipe \
  --set FRAGPIPE_TOOLS=/opt/tools/fragpipe-24.0/tools \
  --set PYTHON=/opt/conda/envs/metapilot/bin/python \
  --set FUNCDB=/data/metapilot/function \
  --print
```

Validate the setup:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar resources-check
```

List registered MAG databases:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar mag-list
```

Download and prepare a supported MAG catalogue before registration:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar mag-download \
  --catalogue human-gut \
  --version v2.0.2 \
  --out /data/metapilot/MAG \
  --resume
```

Register the downloaded or otherwise prepared local MAG database:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar mag-register \
  --catalogue human-gut \
  --version v2.0.2 \
  --identifier MGYG \
  --path /data/metapilot/MAG/human-gut/v2.0.2
```

Validate and run a parameter file:

```bash
java -jar metapilot-cli-1.0.0-shaded.jar params-validate \
  --params /data/studies/study_parameters.json

java -jar metapilot-cli-1.0.0-shaded.jar run \
  --params /data/studies/study_parameters.json
```

Current headless execution supports:

- `DIA_MARKER`: Core/Focus, Comprehensive, Deep together, and Deep per file.
- `DDA_MAG_FRAGPIPE`: MAG-aware DDA FragPipe workflow.

Other workflows are currently GUI workflows.

### MAG Database Preparation

MAG-aware workflows need a registered MAG database. A typical database folder contains:

- `original_db/`: per-genome protein FASTA files.
- `genomes-taxa.tsv`: genome-to-taxonomy mapping.
- `Marker4.fasta`: marker proteins used for marker-based evidence search.
- `Marker4.predicted.tsv`: DeepDetect peptide detectability output.
- `Marker4_highpep.fasta`: high-detectability marker peptides for large marker sets.
- `libraries/Marker4/Marker4.speclib`: Marker4 spectral library for DIA Marker.
- `predicted/`: per-genome predicted peptide files used by DIA Seed mapping.

For DIA Marker, make sure the Marker4 library is available or can be generated before starting the study.

## Choosing A Workflow

Use this table as the first decision point.

| Goal | Recommended workflow |
| --- | --- |
| Standard DDA search against an integrated gene catalogue or fixed protein database | DDA IGC |
| MAG-aware DDA search with genome selection and refined sample FASTA | DDA MAG |
| DIA search starting from a MAG database and Marker4 library | DIA Marker |
| DIA search starting from existing user-selected spectral libraries | DIA Seed |
| Broader DIA discovery with iterative refinement | DIA Marker Deep or DIA Seed Deep |
| Focused follow-up for selected taxa, functions, or proteins | DIA Target |
| Folder monitoring as DIA files arrive | DIA Real-Time |
| Linux/HPC execution today | DIA Marker or DDA MAG FragPipe |

### DDA IGC

Choose DDA IGC when you already have a conventional protein database or integrated gene catalogue and do not need MAG genome selection.

Available search engines include MaxQuant, pFind, FragPipe, and Sage. This is the most familiar route for users who want a standard DDA search and protein-level outputs.

### DDA MAG

Choose DDA MAG when you have DDA data and want MetaPilot to reduce the search space using MAG evidence.

The workflow performs an initial marker or high-abundance search, ranks candidate genomes, builds a sample-specific FASTA from selected genomes, and then performs the final DDA search and downstream taxonomy export.

Use Focus/Core mode when you want a smaller, stricter genome set. Use Comprehensive mode when you prefer broader genome inclusion and can accept a larger refined search database.

### DIA Marker

Choose DIA Marker when your starting point is a MAG database and you want MetaPilot to use the database's Marker4 spectral library as the initial DIA evidence source.

This is the main MAG-aware DIA workflow. It is a good default when users do not already have study-specific DIA libraries.

Available modes:

| Mode | When to use it |
| --- | --- |
| Core / Focus | Default mode; stricter genome selection and smaller refined library |
| Comprehensive | Broader genome selection and larger refined library |
| Deep together | Iterative discovery across all DIA files before final quantification |
| Deep per file | Iterative discovery independently per raw file; useful for heterogeneous sample sets but more resource-intensive |

### DIA Seed

Choose DIA Seed when you already have spectral libraries and want to use them as the initial search space.

Direct Seed searches the selected libraries directly. Seed Focus, Comprehensive, and Deep modes map identified peptides back to MAG proteins, rank genomes, build a refined library, and run a final quantifying search.

Use DIA Seed when the existing libraries are relevant to the sample type or when you want to combine standalone libraries with MAG-aware refinement.

### DIA Target

Choose DIA Target for targeted follow-up analyses. This workflow is intended for selected taxa, functions, or proteins rather than whole-community genome selection.

### DIA Real-Time

Choose DIA Real-Time when DIA files are being generated over time and you want MetaPilot to monitor an input folder. This is a GUI workflow for interactive or operational monitoring.

## Interpreting Results

MetaPilot result folders vary by workflow, but the most important outputs are consistent across MAG-aware workflows.

### Genome Selection Outputs

| File | Meaning |
| --- | --- |
| `genome_evidence.tsv` | Evidence used to rank genomes, including peptide support and marker evidence where available |
| `genome_ranked.tsv` | Candidate genomes ordered by evidence and incremental contribution |
| `genome_selected.tsv` | Genomes selected for the refined FASTA or refined DIA library |

Start with `genome_selected.tsv` to see which genomes MetaPilot kept for the final search. Then inspect `genome_ranked.tsv` to understand high-confidence genomes, borderline genomes, and the long tail of lower-evidence candidates.

Focus/Core mode usually selects fewer genomes. Comprehensive mode usually selects more genomes. Deep modes may add evidence through iterative rounds before building the final refined search space.

### Final DIA Outputs

Common DIA outputs include:

| File | Meaning |
| --- | --- |
| `combined/combined.parquet` | Final DIA-NN report used for downstream export |
| `combined/combined.pr_matrix.tsv` | Precursor or peptide-level quantitative matrix |
| `combined/combined.pg_matrix.tsv` | Protein-group quantitative matrix |
| `marker/combined.speclib` or `seed/combined.speclib` | Refined DIA library built from selected genomes |

Use the `.pr_matrix.tsv` file for peptide or precursor-level abundance patterns and the `.pg_matrix.tsv` file for protein-group abundance patterns.

### DDA MAG Outputs

DDA MAG workflows produce tool-specific identification and quantification files, plus a refined sample-specific FASTA for the final search.

The key interpretation point is that the final DDA search is run against the refined FASTA built from selected genomes. Use the genome selection files to understand why that refined search space was chosen.

### Taxonomy Outputs

| File | Meaning |
| --- | --- |
| `taxonomy_analysis/Genome.tsv` | Genome-level abundance and support table |
| `taxonomy_analysis/Taxa_report.tsv` | Aggregated taxonomic report across ranks |
| `taxonomy_analysis/Genome.html` | Interactive genome-level report when generated |

Use `Genome.tsv` when you need the most detailed genome-level view. Use `Taxa_report.tsv` when you want a higher-level summary by taxonomic rank.

When interpreting taxonomy results, prioritize genomes or taxa with both abundance and peptide support. A high abundance estimate with weak peptide evidence should be reviewed cautiously.

### Functional Annotation Outputs

Functional outputs are written under `functional_annotation/` when functional annotation is enabled.

| File | Meaning |
| --- | --- |
| `annotated_proteins.tsv` | Protein-level annotations with peptide or PSM support and intensities |
| `functions.tsv` | Compatibility copy of the annotated protein table |
| `functions.html` | Interactive functional annotation report |
| `KEGG_ko_report.tsv` | KEGG Orthology summary |
| `Gene_Ontology_id_report.tsv` | Gene Ontology summary |
| `EC_id_report.tsv` | Enzyme Commission summary |
| `COG_accession_report.tsv` | COG summary |
| `NOG_accession_report.tsv` | eggNOG/NOG summary |

Use functional reports to summarize biological activity after checking that the underlying protein and peptide evidence is strong enough for the question being asked.

### Logs And Troubleshooting

Check the workflow log when a run fails or produces unexpected output:

- GUI runs write logs in the selected result folder and GUI log locations.
- Headless runs write `logging.log` in the study output directory when the path can be resolved.

Common checks:

- Confirm the required external tool paths are valid.
- Confirm the MAG database entry resolves to the expected local folder.
- Confirm Marker4 files and libraries exist for DIA Marker.
- Confirm `predicted/` peptide mapping files exist for DIA Seed.
- Confirm the output directory is writable.
- For Deep per file, reduce concurrency if memory usage is too high.

## Notes For Release Users

This README is the main user-facing guide for the public release repository. It is intentionally self-contained so users can start from the release page alone. If future releases add expanded examples, screenshots, or parameter templates, the essential installation, workflow-choice, and result-interpretation information should still remain here.

When reporting a problem, include:

- MetaPilot version.
- Operating system.
- Workflow name and mode.
- External tool versions, especially DIA-NN, FragPipe, pFind, MaxQuant, Sage, and DeepDetect where relevant.
- The parameter JSON or a description of the GUI settings.
- The workflow log, such as `logging.log`.
- The last visible error message.

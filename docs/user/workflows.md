# MetaPilot workflow design (DDA + DIA)

This document describes **how MetaPilot is designed to run MAG-aware workflows** on **DDA** and **DIA** datasets.
It is written to support generating a **single workflow chart** for the project.

> **Copy/paste chart prompt (condensed, Mermaid-ready wording)**  
> MetaPilot MAG workflows follow one backbone with different “evidence search spaces”:
> 1) **Convert spectra** (DDA: to tool formats; DIA: to DIA-NN formats).  
> 2) **Initial evidence search (small marker/high-abundance search space)**:  
>    - DDA MAG: search a MAG-provided **HAP/marker FASTA** (“high abundant database”) to get first-pass PSM evidence.  
>    - DIA Marker: search **Marker4 spectral library** to get first-pass precursor evidence (parquet).  
>    - DIA Seed: search **user-selected libraries** to get first-pass precursor evidence (parquet).  
> 3) **Genome evidence building**: map identified peptides to candidate genomes via protein IDs and MAG metadata; optionally restrict to Marker4 peptides (4-category markers).  
> 4) **Genome ranking + selection**: seed genomes (privilege + optional taxonomic representatives), then greedy set-cover ordering by incremental new peptides, then choose a stop point using a “knee” on the incremental-gain curve (Focus vs Comprehensive differs).  
> 5) **Refined search space construction**:  
>    - DDA MAG: build **sample-specific FASTA** from selected genomes (optionally add host DB).  
>    - DIA Marker/Seed: build **refined microbiome library**: combined.fasta → combined.speclib (and convert to parquet library when multi-lib search requires it).  
>    - Optional peptide-based strategy: build combined.fasta from predicted peptides above a detectability threshold to reduce library size.  
> 6) **Refined quantifying search** on the refined DB/library (final IDs + intensities).  
> 7) **Taxonomy export** in a unified DDA-aligned format: `taxonomy_analysis/Genome.tsv` and `Taxa_report.tsv` (GUI workflows may also generate HTML reports).

The core design pattern across MAG workflows is:

1. **Initial evidence search** using a *small “marker / high-abundance” search space* (fast, broad signal).
2. **Genome ranking + selection** using the initial evidence.
3. **Refined search space construction** (sample-specific FASTA or refined DIA library).
4. **Refined identification + quantification** using the refined search space.
5. **Unified taxonomy export** (DDA-aligned table layout).

---

## Terminology and on-disk artifacts

MetaPilot uses a consistent set of folders and files across workflows.

- **MAG database folder** (`MagDbItem.getCurrentFile()`):
  - `Marker4.fasta`: protein FASTA containing “Marker4” proteins (4-marker categories).
  - `Marker4_highpep.fasta`: peptide-only FASTA derived from Marker4 detectability prediction (used to keep the Marker4 library manageable when Marker4 is large).
  - `Marker4.predicted.tsv`: peptide detectability predictions for Marker4 proteins.
  - `genomes-taxa.tsv`: genome → lineage mapping used for taxonomy-aware genome selection.
  - `predicted/`: folder of per-genome predicted peptide TSVs used by DIA Seed deep/iterative mapping and peptide-based strategies.

- **DIA-NN outputs** (typical):
  - `*.parquet`: DIA-NN report parquet for identification/quantification.
  - `*.pr_matrix.tsv`: peptide quant matrix.
  - `*.pg_matrix.tsv`: protein quant matrix.
  - `*.speclib`: DIA-NN spectral library.
  - `*.lib.parquet` (naming varies): a parquet-format library derived from a speclib for “multi-library” searching.

- **Genome selection outputs** (written by genome ranking/selection step):
  - `genome_evidence.tsv`: evidence table used for ranking.
  - `genome_ranked.tsv`: genomes sorted by evidence (and used later for labels like core/mid/tail).
  - `genome_selected.tsv`: selected genomes after applying a selection rule (Core/Focus vs Comprehensive).

- **Unified taxonomy export**:
  - `taxonomy_analysis/Genome.tsv`: **canonical** genome-level table format (DDA-aligned).
  - `taxonomy_analysis/Taxa_report.tsv`: rank-level aggregation.

---

## How the Marker4 inputs are generated (MAG database prerequisites)

Several workflows rely on Marker4 as the initial evidence search space (DDA “high abundant” stage and DIA Marker stage).
MetaPilot can generate missing Marker4 inputs for a MAG DB (GUI-triggered when a DB is selected).

### 1) What “Marker4” means (the 4 marker categories)

Marker4 is not “4 proteins”. It is a **4-category marker system**: proteins are treated as marker proteins if their FASTA headers match **one of four functional categories**.

MetaPilot categorizes a protein header into one of these categories:

- **Ribosomal** (also called “Anchor” in some naming): ribosomal proteins
- **Replication** (also called “HighRes”): high-resolution, conserved replication/transcription markers such as gyrase / RNA polymerase / DNA polymerase family markers
- **Elongation** (also called “Growth/Abundance”): elongation factor family markers (EF‑Tu/EF‑G/EF‑Ts/EF‑P etc.)
- **Synthetase** (also called “Stability”): aminoacyl‑tRNA ligases/synthetases (aaRS families)

#### How MetaPilot detects a marker protein (regex rules)

The categorization logic supports both:

- **Descriptive headers** (protein names in the FASTA header line), and
- **Gene-name-only headers** (short tokens like `rplA`, `gyrB`, `tuf`, `ileS`).

For descriptive headers, MetaPilot applies case-insensitive regex matching such as:

- Ribosomal: contains “ribosomal protein”
- Replication: contains gyrase / RNA polymerase / DNA polymerase / `rpoA|rpoB|rpoC` / `dnaE|dnaQ` / `gyrA|gyrB`
- Elongation: contains “elongation factor” or `ef-tu|efg|efts|efp` etc.
- Synthetase: contains “aminoacyl‑tRNA ligase/synthetase” (also matches aa‑tRNA spelling variants)

For gene-name-only headers, MetaPilot tokenizes the header and uses gene-pattern matching such as:

- Ribosomal: `rpl*`, `rps*`, `rpm*`
- Replication: `gyrA|gyrB|rpoA|rpoB|rpoC|dnaE|dnaQ`
- Elongation: `tuf|tufA|tufB|tsf|fusA|efp`
- Synthetase: aaRS genes (`ileS`, `leuS`, `lysS`, … plus `metG`, `pheS`, `pheT`, etc.)

This is why Marker4 is robust across databases with different header styles.

### 2) Build `Marker4.fasta` (Marker4 proteins)

- Source: **MAG database original protein FASTAs** in the “original DB folder”.
- Process: extract proteins belonging to four marker categories into a single FASTA.
- Output: `Marker4.fasta` inside the MAG DB folder.

### 2) Predict peptide detectability (`Marker4.predicted.tsv`)

- Tool: **DeepDetect**.
- Input: `Marker4.fasta`
- Output: `Marker4.predicted.tsv`
- Purpose: provide a detectability score per peptide so MetaPilot can keep only “high-detectability” peptides when needed.

### 3) Build `Marker4_highpep.fasta` (peptide-only, high-detectability subset)

- Trigger: if `Marker4.fasta` contains more than a configured protein count threshold (to keep the library size manageable).
- Process:
  - read `Marker4.predicted.tsv`
  - keep peptides with detectability ≥ a default threshold (e.g. 0.8 in the GUI)
  - write them as a peptide FASTA (one peptide per FASTA entry)
- Output: `Marker4_highpep.fasta`

### 4) Build `genomes-taxa.tsv`

- Source: MAG DB taxa table.
- Process: derive a two-column file `Genome\tLineage` used during genome selection.
- Output: `genomes-taxa.tsv`

### 5) Build the **Marker4 spectral library** (DIA Marker prerequisite)

Once Marker4 FASTA (either `Marker4.fasta` or `Marker4_highpep.fasta`) exists:

- Use DIA-NN to generate `libraries/Marker4/Marker4.speclib` from the effective Marker4 FASTA.
- Also generate:
  - `libraries/Marker4/Marker4.parquet` (library converted to parquet form)
  - `libraries/Marker4/Marker4_info.txt` (metadata like peptide count, enzyme, missed cleavages)

This “Marker4 library bundle” is what DIA Marker uses for its initial evidence search.

---

## DDA MAG workflows (pFind / Sage / FragPipe): two-stage MAG-aware search

These workflows share the same high-level design:

### Step A — Convert spectra

- Convert raw MS files into the tool’s required input format (e.g., MGF/PF2 for pFind, mzML for Sage).

### Step B — Initial evidence search (high-abundance / marker stage)

- Use the MAG database’s **high-abundance/marker FASTA** (often referred to as the “HAP FASTA” in code and logs).
- Run a fast initial search to get confident identifications that indicate which genomes are present.

Outputs: tool-specific identification result files (PSMs/proteins) inside per-sample subfolders.

### Step C — Genome ranking and selection (from initial evidence)

From the initial evidence:

- Parse peptide/protein evidence and attribute it to genomes (via protein naming conventions and MAG metadata).
- Rank genomes by evidence.
- Select genomes using one of the selection modes:
  - **Core / Focus**: fewer genomes, stricter selection. `Focus` remains the legacy persisted name.
  - **Comprehensive**: more genomes, broader selection

Outputs per sample:

- `*_evidence.tsv`, `*_ranked.tsv`, `*_selected.tsv`

### Step D — Build the refined sample-specific database (FASTA)

Using selected genomes:

- Concatenate protein FASTAs from the MAG original DB folder (`<genome>.faa` / `<genome>.fasta`) into a **sample-specific FASTA**.
- Optionally append host proteins if enabled.

Output: `<sample>.fasta` placed in the sample’s result folder (the refined search DB).

### Step E — Refined search (main identification pass)

- Run the full database search against the **sample-specific FASTA**.
- Combine results across samples/runs (workflow-specific).

### Step F — Quantification

- Run quantification (label-free and/or isobaric, depending on workflow settings).
- Produce peptide/protein intensities that can be aggregated per genome.

### Step G — Taxonomy analysis + export (unified format)

- Map refined identifications to genomes and taxonomy ranks.
- Export to:
  - `taxonomy_analysis/Genome.tsv` (canonical layout)
  - `taxonomy_analysis/Taxa_report.tsv`

---

## DIA Marker workflows (Core/Focus / Comprehensive): marker-evidence -> genome selection -> refined DIA library search

DIA Marker mirrors the DDA “marker initial search” strategy, but uses DIA-NN and spectral libraries.

### Step A — Convert spectra

- Convert raw DIA files into DIA-NN compatible format.

### Step B — Initial evidence search (Marker4 library)

- Search the **Marker4 spectral library** against DIA files to obtain confident precursor evidence.
- Output: a marker-stage DIA-NN report (parquet). If a combined parquet cannot be produced, MetaPilot can fall back to per-file parquet outputs and aggregate precursors from them.

### Step C — Genome ranking and selection

- Use the identified precursors (q-value filtered) as genome evidence.
- Rank genomes and write:
  - `genome_evidence.tsv`
  - `genome_ranked.tsv`
  - `genome_selected.tsv`
- Selection mode is controlled by Core/Focus vs Comprehensive.

### Step D — Build refined microbiome library (from selected genomes)

- Create `marker/combined.fasta` by concatenating proteins for selected genomes from the MAG original DB folder.
- Use DIA-NN to build `marker/combined.speclib`.

### Step E — Refined quantifying search (final)

- Run DIA-NN quantifying search using `marker/combined.speclib` to produce:
  - `combined/combined.parquet`
  - `combined/combined.pr_matrix.tsv`
  - `combined/combined.pg_matrix.tsv`

### Step F — Taxonomy export (unified)

- Export `taxonomy_analysis/Genome.tsv` and `Taxa_report.tsv` using the DDA-aligned layout (including unique peptide and I/L ratio columns).

---

## DIA Deep Marker workflows (Deep together / Deep per file): iterative refinement before the final quantifying search

Deep Marker extends DIA Marker by adding iterative rounds:

1. Marker4 evidence search → genome ranking (`genome_ranked.tsv`)
2. Iterative rounds:
   - each round takes the next top-N genomes from `genome_ranked.tsv`
   - default top-N is 5 genomes per round (`genomes_per_round`)
   - builds a temporary FASTA + speclib
   - runs a **non-quant** identification search to discover more peptides
3. Stop rule:
   - the shared deep stop rule evaluates the cumulative discovery curve after an initial minimum of 8 rounds
   - it computes chord-distance, piecewise and kneedle-style knee candidates
   - it stops when new peptides are below the relative-gain threshold (default 1.5% of cumulative peptides) for the configured plateau requirement
4. Final step: merge evidence and run one **quantifying** search to produce final matrices/parquet
5. Export taxonomy in unified format

`Deep together` is the legacy `DEEP` mode: iterative evidence is accumulated across all DIA files before the final quantifying search.

`Deep per file` uses the same core logic, but runs the iterative marker workflow independently for each DIA raw file under `marker/<rawname>/...` before building one global combined search.

---

## DIA Seed workflows (Direct / Core/Focus / Comprehensive / Deep together)

Seed workflows start from **user-selected spectral libraries** rather than Marker4.

### Direct Seed

- Search the selected libraries directly (no MAG genome selection/refinement).
- If all selected libraries are microbiome/MAG-matched, MetaPilot can continue through report export.
- If all selected libraries are standalone/non-MAG, MetaPilot stops after identification because there is no MAG mapping for taxonomy export.
- Mixed standalone + microbiome sets are handled by using microbiome-compatible libraries for MAG-aware refinement.

### Seed iterative (Core/Focus / Comprehensive)

1. Initial library search with user-selected libraries → produce seed parquet(s).
2. Map identified peptides to MAG proteins/genomes using MAG `predicted/*.predicted.tsv`.
3. Rank and select genomes (`genome_ranked.tsv`, `genome_selected.tsv`).
4. Build refined microbiome library under `seed/`:
   - default: proteins from selected genomes → `seed/combined.fasta` → `seed/combined.speclib`
   - optional: peptide-based FASTA from predicted peptides above a detectability threshold → `seed/combined.fasta` → `seed/combined.speclib`
5. Final quantifying search:
   - **mixed-library mode**: convert `seed/combined.speclib` → `seed/combined.lib.parquet` and combine with standalone parquet libs
   - otherwise: use `seed/combined.speclib` directly
6. Export unified taxonomy tables.

### Seed deep together (iterative rounds)

Seed deep adds iterative rounds similar to Deep Marker, but with Seed-specific peptide→MAG mapping and ranking.
It supports a peptide-budget mode where each round adds predicted peptides genome-by-genome until a peptide budget trigger is exceeded.

Seed deep per file is part of the canonical naming model, but it is still blocked by GUI validation with an explicit unsupported-mode error instead of routing it to another workflow.

---

## DIA workflow taxonomy and compatibility

DIA MAG workflows are modeled as:

- **Initial search strategy**: `Marker` or `Seed`
- **Refinement mode**: `Direct`, `Core`, `Comprehensive`, `DeepTogether`, or `DeepPerFile`

`Direct` is seed-only. `Focus` is a legacy alias for `Core`. Existing saved parameter JSON files keep their persisted values for compatibility.

| Legacy UI/task name | Canonical workflow name | Persisted value |
| --- | --- | --- |
| Marker Focus / `MetaDiaMarkerTask` | `DiaMarkerCoreTask` | `FOCUS` |
| Marker Comprehensive / `MetaDiaMarkerTask` | `DiaMarkerComprehensiveTask` | `COMPREHENSIVE` |
| Deep Marker / `MetaDiaDeepMarkerTask` | `DiaMarkerDeepTogetherTask` | `DEEP` |
| Deep Marker per raw file / `MetaDiaDeepMarkerTask` | `DiaMarkerDeepPerFileTask` | `DEEP_PER_FILE` |
| Seed Direct / `MetaDiaSeedDirectTask` | `DiaSeedDirectTask` | `DIRECT` |
| Seed Focus / `MetaDiaSeedIterativeTask` | `DiaSeedCoreTask` | `FOCUS` |
| Seed Comprehensive / `MetaDiaSeedIterativeTask` | `DiaSeedComprehensiveTask` | `COMPREHENSIVE` |
| Seed Deep Comprehensive / `MetaDiaDeepComprehensiveTask` | `DiaSeedDeepTogetherTask` | `DEEP` |
| Seed Deep Per File | `DiaSeedDeepPerFileTask` | `DEEP_PER_FILE` |

Target and Real-Time DIA workflows are separate DIA workflow families and are intentionally left out of the shared Marker/Seed MAG backbone.

---

## Shared MAG taxonomy exporter boundary

`MagTaxonomyTsvExporter` is the canonical writer for MAG taxonomy outputs:

- `taxonomy_analysis/Genome.html`
- `taxonomy_analysis/Genome.tsv`
- `taxonomy_analysis/Taxa.tsv` when DIA-style genome scoring is available
- `taxonomy_analysis/Taxa_report.tsv`
- `taxonomy_analysis/genome_rank_scores.tsv`

GUI task code should normally call `MagTaxonomyAnalysisService`, not `MagTaxonomyTsvExporter` directly. The service resolves the MAG taxonomy file, copies `genome_ranked.tsv` labels when present, and then delegates to the exporter.

Classes/workflows that should use this shared path:

- DDA MAG GUI tasks through `MetaPilotMagTask.performTaxonomyAnalysis`
- DIA Marker tasks through `MetaDiaTask.exportReport`
- DIA Seed Direct only when selected libraries are MAG-matched and the workflow reaches MAG report export
- DIA Seed Core/Comprehensive/Deep Together through the shared Marker/Seed report export path
- DDA/DIA headless post-analysis runners once they have genome intensity and peptide-support maps

Classes/workflows that should not use the shared Marker/Seed export boundary yet:

- DIA Target, because it is a targeted library workflow rather than the shared MAG genome-selection backbone
- DIA Real-Time, because it is a monitor/automation workflow rather than a shared Marker/Seed batch workflow
- standalone Seed Direct runs that stop after identification because no MAG mapping exists

---

## Shared MAG functional annotation outputs

Formal-release functional annotation uses an in-memory annotation result rather than `functions.xml` as the workflow data model.

Canonical outputs under `functional_annotation/`:

- `annotated_proteins.tsv`: one row per annotated protein, with peptide/PSM support, intensities, and functional annotations.
- `functions.tsv`: compatibility copy of `annotated_proteins.tsv` for older report consumers.
- `functions.html`: interactive functional annotation report.
- `KEGG_ko_report.tsv`, `Gene_Ontology_id_report.tsv`, `EC_id_report.tsv`, `COG_accession_report.tsv`, `NOG_accession_report.tsv`, and category reports.

Function summary reports should include definition/name columns from `func_def.db` where available. Protein confidence scores remain in the protein result files (`razor_proteins.tsv`, `all_proteins.tsv`) and are not duplicated into the functional annotation tables.

`functions.xml` is no longer a canonical workflow artifact. If encountered in older output folders, it should be treated as a pre-release cache/intermediate file.

---

## DIA Target workflow (targeted, not part of the MAG backbone)

`DIA Target` is a separate workflow from `DIA Seed`.

- It uses a **pre-generated targeted library** for specific taxa/functions/proteins of interest.
- It does **not** use Marker4-based genome ranking or the shared MAG sample-specific refinement backbone described above.
- It is best documented as a targeted follow-up workflow rather than as part of the shared Marker/Seed genome-selection design.

---

## Notes for workflow charts

If you need a single chart, a useful abstraction is:

- **Evidence search space**
  - DDA MAG: “high-abundance / marker FASTA”
  - DIA Marker: “Marker4 spectral library”
  - DIA Seed: “user-selected libraries”
  - DIA Target: targeted library generated for selected taxa/functions/proteins

- **Selection**
  - genome ranking + Focus/Comprehensive selection
  - DIA Target does not use this shared genome-ranking stage

- **Refinement**
  - DDA: sample-specific FASTA
  - DIA: refined microbiome library (speclib/parquet)

- **Final search**
  - DDA: refined DB search + quant
  - DIA: refined DIA-NN quant search

- **Unified export**
  - taxonomy tables in DDA-aligned format

---

## How genome ranking and selection works (conceptual algorithm, code-aligned)

This section explains *what the genome selection stage actually does*, without requiring knowledge of internal method names.

### Inputs to genome selection

Genome selection operates on **identification evidence** from the initial stage:

- **DDA**: PSMs / peptide sequences / protein IDs from the initial high-abundance (HAP) search.
- **DIA Marker**: precursor-level identifications from the Marker4 library search (parquet), filtered at q≤0.01.
- **DIA Seed**: precursor-level identifications from the user-library search (parquet), then mapped to MAG proteins via `predicted/*.predicted.tsv`.

The core requirement is: each identified peptide can be associated with **one or more genomes** through its matched protein(s).

### Evidence extraction: peptides → genomes (marker peptides)

When Marker4 is available, MetaPilot focuses on **marker peptides**:

- A peptide is considered a “marker peptide” if at least one of its matched proteins is categorized into one of the four Marker4 categories (Ribosomal/Replication/Elongation/Synthetase).
- For each genome, MetaPilot accumulates:
  - total marker peptides (unique sequences)
  - total unique marker peptides (peptides mapping to exactly one genome)
  - total marker PSM support (counts of PSM/precursor evidence IDs)
  - coverage per marker category (how many categories have evidence; and how many categories have unique evidence)

This produces a **genome evidence table** (`genome_evidence.tsv`).

### Hard filter (Marker4-only safety filter)

Before ranking, MetaPilot can apply a hard filter that keeps only genomes with strong marker-category diversity.
If the hard filter would remove everything, MetaPilot falls back to greedy-only selection and selects until full peptide coverage is reached.

### Seed genomes: “privilege” + (optional) representative selection

To avoid missing key clades and to stabilize selection, MetaPilot can select a “seed” set before greedy ranking:

1. **Privileged genomes**
   - Genomes with **≥3 unique marker categories** are candidates.
   - A genome must also add at least ~1% new peptides relative to what is already covered.

2. **Representative genomes (optional auto-representation)**
   - After privileged genomes, MetaPilot can add additional genomes to represent under-covered taxonomy groups (class/order/family/genus),
     again requiring a minimum incremental peptide gain.

These are added to the selected set first, and recorded with reasons such as `PRIVILEGE_UNIQUE_3TYPES`, `CLASS_REP`, `ORDER_REP`, etc.

### Greedy genome ordering (set-cover style)

After seeding, MetaPilot ranks remaining genomes using a greedy set-cover approach:

- Maintain a set of **covered peptides** (initially those covered by seed genomes).
- Repeatedly choose the genome with the **largest incremental gain**:
  - \(gain = \#\{\text{peptides in genome not yet covered}\}\)
- Tie-breakers favor genomes with higher overall marker evidence, then stable ID ordering.

This yields a **greedy order** and per-step statistics:

- **marginalGain** at each step
- **ratio** per step: incremental gain / total covered
- **moving average** of the ratio (window size 5)

The greedy order becomes the backbone of `genome_ranked.tsv`.

### Choosing the stop point (knee selection)

MetaPilot chooses how many greedy steps to keep using a “knee” on the moving-average incremental-gain curve.

It computes several knee candidates and then picks a final knee step depending on the workflow mode:

- **Chord-distance knee** (max distance from line connecting endpoints)
- **Kneedle-style knee**
- **Curvature-based knee**
- **Piecewise breakpoint knee**
- **Slope-change knee**

The final selected step differs by mode:

- **Focus**: stricter / earlier stopping (fewer genomes)
- **Comprehensive**: later stopping (more genomes)

Finally, MetaPilot selects:

- all seed genomes
- plus greedy genomes up to the chosen knee step

and writes them to `genome_selected.tsv` along with selection reasons.

### Output ordering and labels (core/mid/tail)

`genome_ranked.tsv` records the full ranking order (selected-first plus remaining ranked genomes) and includes:

- evidence columns (marker peptides, unique marker peptides, PSM support, category counts)
- the greedy step, gain, and knee annotations

Downstream export uses this ranked file to attach user-facing labels (e.g., “core/mid/tail”) in `Genome.tsv`.

---

## Focus vs Comprehensive: what actually changes

MetaPilot uses the same evidence extraction and greedy ordering for both modes. The difference is **how the stop step is chosen** on the greedy curve.

### Shared parts

- Both modes can apply:
  - a Marker4 hard filter (when Marker4 categories are available)
  - “privileged” seed genomes (≥3 unique marker categories + minimum incremental gain)
  - optional taxonomic representative genomes (class/order/family/genus)
- Both modes then build a greedy order by repeatedly selecting the genome that contributes the largest number of **previously-uncovered** marker peptides.

### Knee candidates computed from the same curve

After greedy ordering, MetaPilot computes a knee step using multiple knee estimators on the same smoothed incremental-gain curve:

- chord-distance knee (max distance to endpoint chord)
- kneedle-style knee
- curvature-based knee
- piecewise breakpoint knee
- slope-change knee

### How **Focus** chooses the stop step (stricter)

Focus mode chooses a conservative stop step by taking the **median of three** knee estimators:

- kneedle
- curvature
- slope-change

This tends to stop earlier (fewer genomes).

### How **Comprehensive** chooses the stop step (later)

Comprehensive mode chooses the stop step as:

- `max(chord-distance knee, piecewise knee)`

This tends to stop later (more genomes).

### Practical implication for workflow charts

You can represent this as:

- **Focus** = “median knee” (earlier stop)
- **Comprehensive** = “max(chord, piecewise)” (later stop)


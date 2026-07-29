# MetaPilot → MetaAtlas Parquet Export Layer

## Project Design and Upgrade Guide

**Document purpose:**  
This guide defines how to extend the current MetaPilot project with an additional result-export layer that produces QPX-aligned Parquet files for MetaAtlas.

The upgrade is designed to preserve the current MetaPilot analysis workflows while adding a standardized, reusable, and query-friendly output layer for cross-project metaproteomics data integration.

---

# 1. Background

MetaPilot currently performs metaproteomics analysis using:

- **FragPipe** for DDA workflows
- **DIA-NN** for DIA workflows
- MAG-aware genome selection
- Taxonomic assignment
- Functional annotation
- Genome-, taxon-, and function-level reporting

For long-term data reuse, the current result files should be complemented by a standardized export layer.

The proposed architecture separates the project into two connected systems:

## quantms Portal

Stores and serves mass-spectrometry evidence:

- PSMs
- peptide evidence
- protein evidence
- precursor features
- retention times
- spectra-linked evidence
- search provenance
- links to PRIDE raw data

## MetaAtlas

Stores and serves microbiome-specific interpretation:

- genomes
- taxa
- functions
- peptide-to-genome relationships
- protein-to-genome relationships
- protein-to-function relationships
- genome abundance
- taxon abundance
- function abundance
- genome-selection evidence
- cross-project microbiome queries

MetaAtlas should not duplicate detailed PSM or spectrum-level data. Instead, it should reference the corresponding QPX/quantms evidence using stable identifiers.

---

# 2. Design Goal

The goal is to add a new **MetaAtlas Export Layer** after the existing MetaPilot reporting stage.

The current MetaPilot workflows should remain unchanged.

The new layer should:

1. Read finalized MetaPilot outputs.
2. Normalize identifiers.
3. Build QPX-aligned entity and mapping tables.
4. Export the tables as Parquet files.
5. Add file-level metadata and provenance.
6. Validate schema compliance.
7. Package the files for MetaAtlas ingestion.
8. Preserve links to PRIDE and quantms/QPX records.

The export layer should work for both DDA and DIA workflows.

---

# 3. High-Level Architecture

```text
PRIDE raw files
      │
      ▼
MetaPilot analysis
      │
      ├── DDA: FragPipe
      │
      └── DIA: DIA-NN
      │
      ▼
Current MetaPilot result layer
      │
      ├── PSM / precursor evidence
      ├── peptide results
      ├── protein results
      ├── genome reports
      ├── taxonomy reports
      ├── function reports
      └── genome-selection results
      │
      ▼
New MetaAtlas Export Layer
      │
      ├── identifier normalization
      ├── schema transformation
      ├── Parquet writing
      ├── provenance generation
      ├── validation
      └── packaging
      │
      ├──────────────► QPX / quantms export
      │                PSM, peptide, protein, feature evidence
      │
      └──────────────► MetaAtlas export
                       genome, taxonomy, function, abundance,
                       mappings, genome-selection evidence
```

---

# 4. Core Design Principles

## 4.1 Preserve the current workflow

The export layer must not alter:

- search parameters
- identification results
- FDR control
- protein inference
- genome selection
- taxonomic aggregation
- functional aggregation
- existing TSV or HTML reports

The new layer is an additional downstream transformation only.

## 4.2 One entity or relationship per Parquet view

Each Parquet file should represent one logical entity or one many-to-many relationship.

Examples:

- one genome per row
- one taxon per row
- one function per row
- one peptide–genome relationship per row
- one protein–function relationship per row
- one genome abundance value per sample per row

## 4.3 Use stable identifiers

Display names must not be used as primary keys.

Preferred identifiers include:

- PRIDE accession for projects
- SDRF source name or a MetaPilot-generated stable sample ID
- raw filename or QPX run identifier for runs
- unmodified peptide sequence for peptide-level biological mappings
- stable protein accession from the search FASTA
- MGnify or MAG genome identifier
- GTDB or NCBI taxon identifier
- KEGG KO, COG, EC, GO, or eggNOG identifiers

## 4.4 Separate reference entities from observations

Reference tables describe what an entity is.

Observation tables describe where it was detected and how much evidence supports it.

For example:

- `genome.parquet` describes genome metadata.
- `genome_abundance.parquet` describes sample-specific genome detection and abundance.

## 4.5 Preserve provenance

Every exported result must record:

- MetaPilot version
- search-engine version
- workflow mode
- catalogue name and release
- taxonomy source and release
- functional annotation source and release
- FASTA checksums
- important analysis parameters
- export schema version
- creation time

## 4.6 Do not store PSM-level data in MetaAtlas

PSMs, spectra, fragment evidence, and detailed chromatographic information should remain in QPX/quantms.

MetaAtlas should store only the biological interpretation and compact summary evidence required for microbiome queries.

---

# 5. Proposed Output Structure

For each PRIDE project:

```text
MetaAtlasExport/
└── PXDxxxxxx/
    ├── PXDxxxxxx.dataset.parquet
    ├── PXDxxxxxx.sample.parquet
    ├── PXDxxxxxx.run.parquet
    ├── PXDxxxxxx.peptide.parquet
    ├── PXDxxxxxx.protein.parquet
    ├── PXDxxxxxx.genome.parquet
    ├── PXDxxxxxx.taxon.parquet
    ├── PXDxxxxxx.function.parquet
    ├── PXDxxxxxx.peptide_genome.parquet
    ├── PXDxxxxxx.protein_genome.parquet
    ├── PXDxxxxxx.protein_function.parquet
    ├── PXDxxxxxx.genome_taxon.parquet
    ├── PXDxxxxxx.genome_abundance.parquet
    ├── PXDxxxxxx.taxon_abundance.parquet
    ├── PXDxxxxxx.function_abundance.parquet
    ├── PXDxxxxxx.genome_selection.parquet
    ├── PXDxxxxxx.provenance.parquet
    ├── PXDxxxxxx.ontology.parquet
    ├── PXDxxxxxx.manifest.json
    └── PXDxxxxxx.sdrf.tsv
```

Not every file is required in the first implementation.

---

# 6. Minimum Viable Export

The first MetaAtlas export version should generate the following files:

1. `dataset.parquet`
2. `sample.parquet`
3. `genome.parquet`
4. `taxon.parquet`
5. `function.parquet`
6. `peptide_genome.parquet`
7. `genome_abundance.parquet`
8. `taxon_abundance.parquet`
9. `function_abundance.parquet`
10. `provenance.parquet`
11. `manifest.json`

The following files can be added later:

- `run.parquet`
- `peptide.parquet`
- `protein.parquet`
- `protein_genome.parquet`
- `protein_function.parquet`
- `genome_taxon.parquet`
- `genome_selection.parquet`
- `ontology.parquet`

---

# 7. File Naming Convention

Use the QPX-style naming pattern:

```text
{project_accession}.{view_name}.parquet
```

Example:

```text
PXD008738.genome.parquet
PXD008738.taxon_abundance.parquet
PXD008738.provenance.parquet
```

For internal or unpublished projects:

```text
{internal_project_id}.{view_name}.parquet
```

Example:

```text
METAATLAS_TEST_001.genome.parquet
```

---

# 8. Shared Fields

All project-specific observation tables should contain the following fields where applicable:

| Field | Type | Description |
|---|---|---|
| `project_accession` | string | PRIDE or internal project ID |
| `sample_accession` | string | Stable biological sample ID |
| `run_file_name` | string, nullable | Raw file or acquisition run |
| `analysis_unit_id` | string, nullable | Unit used for genome selection or aggregation |
| `metapilot_version` | string | MetaPilot version |
| `schema_version` | string | MetaAtlas schema version |

---

# 9. Core Schemas

## 9.1 `dataset.parquet`

One row per project.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `project_accession` | string | yes | PRIDE or internal accession |
| `project_title` | string | no | Project title |
| `project_description` | string | no | Brief project description |
| `publication_id` | string | no | DOI or PubMed identifier |
| `acquisition_mode` | string | yes | DDA, DIA, or mixed |
| `sample_count` | int32 | yes | Number of biological samples |
| `run_count` | int32 | yes | Number of raw files |
| `metapilot_version` | string | yes | MetaPilot version |
| `schema_version` | string | yes | MetaAtlas format version |
| `qpx_dataset_uri` | string | no | Link or identifier for QPX evidence |
| `pride_uri` | string | no | Link to source PRIDE project |
| `creation_date` | timestamp | yes | Export creation time |

Primary key:

```text
project_accession
```

---

## 9.2 `sample.parquet`

One row per biological sample.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `project_accession` | string | yes | Project accession |
| `sample_accession` | string | yes | Stable sample ID |
| `source_name` | string | yes | SDRF source name |
| `host_taxon_id` | string | no | Host taxonomy identifier |
| `body_site` | string | no | Host body site |
| `environment_biome` | string | no | ENVO biome |
| `environment_feature` | string | no | ENVO feature |
| `environment_material` | string | no | ENVO material |
| `condition` | list<string> | no | Experimental condition |
| `disease` | list<string> | no | Disease metadata |
| `matched_metagenome_available` | bool | no | Whether matched metagenomics exists |
| `sdrf_row_id` | string | no | Source SDRF row |

Primary key:

```text
project_accession + sample_accession
```

---

## 9.3 `genome.parquet`

One row per genome represented in the project.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `genome_id` | string | yes | Stable genome identifier |
| `catalogue_name` | string | yes | MGnify or other catalogue |
| `catalogue_release` | string | yes | Catalogue version |
| `genome_name` | string | no | Human-readable name |
| `taxon_id` | string | no | Primary taxonomic assignment |
| `taxonomy_source` | string | yes | GTDB or NCBI |
| `taxonomy_release` | string | yes | Taxonomy version |
| `completeness` | float32 | no | Genome completeness percentage |
| `contamination` | float32 | no | Genome contamination percentage |
| `protein_count` | int32 | no | Number of predicted proteins |
| `genome_length` | int64 | no | Genome length |
| `assembly_accession` | string | no | Source assembly |
| `genome_uri` | string | no | External resource link |
| `fasta_checksum` | string | no | Genome FASTA checksum |

Primary key:

```text
genome_id + catalogue_release
```

---

## 9.4 `taxon.parquet`

One row per taxonomic entity.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `taxon_id` | string | yes | Stable taxon identifier |
| `taxonomy_source` | string | yes | GTDB or NCBI |
| `taxonomy_release` | string | yes | Taxonomy version |
| `taxon_name` | string | yes | Scientific name |
| `taxon_rank` | string | yes | Species, genus, family, etc. |
| `parent_taxon_id` | string | no | Direct parent |
| `lineage` | list<struct> | no | Structured lineage |
| `synonyms` | list<string> | no | Alternative names |

Suggested lineage structure:

```text
rank
taxon_id
taxon_name
```

Primary key:

```text
taxon_id + taxonomy_source + taxonomy_release
```

---

## 9.5 `function.parquet`

One row per functional entity.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `function_database` | string | yes | KEGG_KO, COG, EC, GO, eggNOG |
| `function_accession` | string | yes | Stable functional ID |
| `function_name` | string | no | Display name |
| `function_level` | string | no | Gene, orthologue, pathway, category |
| `parent_accession` | string | no | Parent function |
| `category_accessions` | list<string> | no | Higher-level categories |
| `description` | string | no | Functional description |
| `database_release` | string | yes | Annotation database version |
| `function_uri` | string | no | External resource link |

Primary key:

```text
function_database + function_accession + database_release
```

---

## 9.6 `peptide_genome.parquet`

One row per peptide–genome relationship.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `sequence` | string | yes | Unmodified peptide sequence |
| `genome_id` | string | yes | Mapped genome |
| `catalogue_release` | string | yes | Genome catalogue release |
| `protein_accessions` | list<string> | yes | Proteins containing the peptide |
| `is_genome_unique` | bool | yes | True if mapped to one genome |
| `genome_mapping_count` | int32 | yes | Number of mapped genomes |
| `taxon_specific_rank` | string | no | Lowest supported taxonomic rank |
| `lowest_common_ancestor_id` | string | no | LCA taxon |
| `mapping_source` | string | yes | Search FASTA, predicted peptides, etc. |
| `qpx_peptide_id` | string | no | Reference to QPX peptide evidence |
| `qpx_dataset_uri` | string | no | QPX dataset link |

Primary key:

```text
sequence + genome_id + catalogue_release
```

---

## 9.7 `genome_abundance.parquet`

One row per genome per sample.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `project_accession` | string | yes | Project accession |
| `sample_accession` | string | yes | Sample accession |
| `genome_id` | string | yes | Genome identifier |
| `catalogue_release` | string | yes | Catalogue version |
| `abundance` | float64 | no | Quantitative value |
| `abundance_type` | string | yes | Raw, normalized, relative, unique-peptide |
| `normalization` | string | yes | Normalization method |
| `peptide_count` | int32 | yes | Supporting peptides |
| `unique_peptide_count` | int32 | yes | Genome-unique peptides |
| `protein_count` | int32 | no | Supporting proteins |
| `feature_count` | int32 | no | DIA features or DDA quantitative features |
| `evidence_score` | float64 | no | MetaPilot evidence score |
| `detected` | bool | yes | Detection flag |
| `quantification_method` | string | yes | Aggregation method |

Primary key:

```text
project_accession + sample_accession + genome_id + abundance_type
```

---

## 9.8 `taxon_abundance.parquet`

One row per taxon per sample.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `project_accession` | string | yes | Project accession |
| `sample_accession` | string | yes | Sample accession |
| `taxon_id` | string | yes | Taxon identifier |
| `taxon_rank` | string | yes | Taxonomic rank |
| `taxonomy_source` | string | yes | GTDB or NCBI |
| `taxonomy_release` | string | yes | Taxonomy release |
| `abundance` | float64 | no | Taxon abundance |
| `abundance_type` | string | yes | Raw, normalized, relative |
| `normalization` | string | yes | Normalization method |
| `supporting_peptide_count` | int32 | yes | Supporting peptides |
| `unique_peptide_count` | int32 | yes | Taxon-unique peptides |
| `supporting_genome_count` | int32 | no | Contributing genomes |
| `minimum_evidence_passed` | bool | yes | Taxon evidence filter |
| `quantification_method` | string | yes | Aggregation method |

Primary key:

```text
project_accession + sample_accession + taxon_id + taxon_rank + abundance_type
```

---

## 9.9 `function_abundance.parquet`

One row per function per sample.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `project_accession` | string | yes | Project accession |
| `sample_accession` | string | yes | Sample accession |
| `function_database` | string | yes | COG, KO, EC, GO, eggNOG |
| `function_accession` | string | yes | Functional identifier |
| `annotation_release` | string | yes | Annotation database release |
| `abundance` | float64 | no | Functional abundance |
| `abundance_type` | string | yes | Raw, normalized, relative |
| `normalization` | string | yes | Normalization method |
| `supporting_protein_count` | int32 | yes | Supporting proteins |
| `supporting_peptide_count` | int32 | no | Supporting peptides |
| `contributing_genome_count` | int32 | no | Contributing genomes |
| `quantification_method` | string | yes | Aggregation method |

Primary key:

```text
project_accession + sample_accession +
function_database + function_accession + abundance_type
```

---

## 9.10 `genome_selection.parquet`

One row per candidate genome per analysis unit.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `project_accession` | string | yes | Project accession |
| `sample_accession` | string | no | Sample accession |
| `analysis_unit_id` | string | yes | Genome-selection unit |
| `genome_id` | string | yes | Genome identifier |
| `catalogue_release` | string | yes | Catalogue release |
| `rank` | int32 | yes | Genome evidence rank |
| `selected` | bool | yes | Selection status |
| `selection_reason` | list<string> | no | Privilege, representative, greedy, etc. |
| `selection_mode` | string | yes | Focus, Comprehensive, Deep |
| `marker_peptide_count` | int32 | no | Marker evidence |
| `unique_marker_peptide_count` | int32 | no | Unique marker evidence |
| `marker_category_count` | int16 | no | Marker category coverage |
| `unique_marker_category_count` | int16 | no | Unique category coverage |
| `marginal_gain` | int32 | no | Newly covered peptides |
| `marginal_gain_ratio` | float64 | no | Relative gain |
| `genome_score` | float64 | no | Final genome score |
| `empirical_pvalue` | float64 | no | Empirical significance |
| `qvalue` | float64 | no | Multiple-testing adjusted value |
| `knee_step` | int32 | no | Selected stopping step |
| `selected_at_step` | int32 | no | Greedy selection step |

Primary key:

```text
analysis_unit_id + genome_id
```

---

## 9.11 `provenance.parquet`

One row per processing step.

| Field | Arrow type | Required | Description |
|---|---|---:|---|
| `step_order` | int32 | yes | Processing order |
| `step_category` | string | yes | Search, selection, annotation, export |
| `step_name` | string | yes | Human-readable step |
| `tool_name` | string | yes | Tool name |
| `tool_version` | string | no | Tool version |
| `tool_uri` | string | no | Container or software URI |
| `parameters` | list<struct> | no | Important key-value parameters |
| `config_json` | string | no | Full configuration |
| `input_views` | list<string> | no | Input tables/files |
| `output_views` | list<string> | no | Output views |
| `database_name` | string | no | Database or catalogue |
| `database_release` | string | no | Database version |
| `database_checksum` | string | no | Database checksum |

Suggested parameter structure:

```text
key
value
```

Example processing steps:

```text
1  raw_conversion
2  initial_marker_search
3  genome_evidence_building
4  adaptive_genome_selection
5  refined_database_construction
6  final_database_search
7  quantification
8  taxonomic_aggregation
9  functional_aggregation
10 metaatlas_export
```

---

# 10. MetaPilot Code Upgrade

## 10.1 Add a new export module

Suggested package:

```text
metapilot.export.metaatlas
```

Suggested classes:

```text
MetaAtlasExporter
MetaAtlasExportContext
MetaAtlasSchemaRegistry
MetaAtlasParquetWriter
MetaAtlasValidator
MetaAtlasManifestWriter
MetaAtlasProvenanceBuilder
MetaAtlasIdentifierResolver
```

## 10.2 Suggested responsibilities

### `MetaAtlasExporter`

Coordinates the complete export.

```java
public interface MetaAtlasExporter {
    MetaAtlasExportResult export(
        MetaAtlasExportContext context
    ) throws MetaAtlasExportException;
}
```

### `MetaAtlasExportContext`

Contains all inputs required for export.

Suggested fields:

```text
project accession
result directory
workflow type
sample metadata
SDRF path
MAG database path
catalogue name
catalogue release
taxonomy source
taxonomy release
annotation versions
MetaPilot version
QPX dataset URI
output directory
```

### `MetaAtlasSchemaRegistry`

Loads and exposes schema definitions.

Schemas should be stored outside Java code where possible:

```text
resources/metaatlas-schema/1.0/
├── dataset.yaml
├── sample.yaml
├── genome.yaml
├── taxon.yaml
├── function.yaml
├── peptide_genome.yaml
├── genome_abundance.yaml
├── taxon_abundance.yaml
├── function_abundance.yaml
├── genome_selection.yaml
└── provenance.yaml
```

### `MetaAtlasParquetWriter`

Writes Arrow records to Parquet.

Recommended properties:

```text
compression = ZSTD
dictionary encoding = enabled
statistics = enabled
row group size = configurable
schema version in footer metadata
project accession in footer metadata
```

### `MetaAtlasValidator`

Checks:

- required files exist
- required columns exist
- Arrow types match schema
- primary keys are unique
- required values are not null
- foreign keys resolve
- abundance values are valid
- taxonomy and functional identifiers are versioned
- project and sample identifiers are consistent
- QPX links are valid strings
- manifest checksums match output files

### `MetaAtlasManifestWriter`

Generates:

```text
PXDxxxxxx.manifest.json
```

Example:

```json
{
  "schema_version": "1.0",
  "project_accession": "PXD008738",
  "created_by": "MetaPilot",
  "metapilot_version": "1.0.0",
  "files": [
    {
      "name": "PXD008738.genome.parquet",
      "view": "genome",
      "rows": 432,
      "size_bytes": 182734,
      "sha256": "..."
    }
  ]
}
```

---

# 11. Input Mapping from Current MetaPilot Outputs

## DDA workflows

Potential input sources:

- FragPipe PSM tables
- FragPipe peptide tables
- FragPipe protein tables
- sample-specific FASTA files
- `genome_evidence.tsv`
- `genome_ranked.tsv`
- `genome_selected.tsv`
- `taxonomy_analysis/Genome.tsv`
- `taxonomy_analysis/Taxa_report.tsv`
- functional annotation reports
- quantitative protein or peptide matrices

## DIA workflows

Potential input sources:

- DIA-NN report Parquet or TSV
- `*.pr_matrix.tsv`
- `*.pg_matrix.tsv`
- refined spectral library
- `genome_evidence.tsv`
- `genome_ranked.tsv`
- `genome_selected.tsv`
- `taxonomy_analysis/Genome.tsv`
- `taxonomy_analysis/Taxa_report.tsv`
- functional annotation reports

## Common mapping logic

```text
MetaPilot genome ID
    └── genome.parquet

MetaPilot peptide-to-protein mapping
    └── peptide_genome.parquet

MetaPilot Genome.tsv
    └── genome_abundance.parquet

MetaPilot Taxa_report.tsv
    └── taxon_abundance.parquet

MetaPilot functional matrices
    └── function_abundance.parquet

MetaPilot genome evidence/ranking/selection
    └── genome_selection.parquet
```

---

# 12. Identifier Normalization

A central identifier resolver should be used by all exporters.

## Project identifiers

```text
project_accession = PRIDE accession where available
```

## Sample identifiers

Preferred order:

1. SDRF source name
2. existing stable MetaPilot sample ID
3. normalized raw-file basename
4. generated deterministic UUID

Do not use row numbers as permanent sample identifiers.

## Peptide identifiers

MetaAtlas biological mapping:

```text
sequence = unmodified amino-acid sequence
```

QPX evidence linking:

```text
qpx_peptide_id = QPX peptidoform or peptide identifier
```

## Protein identifiers

Use the exact stable accession from the search FASTA.

Do not generate new protein identifiers unless required.

## Genome identifiers

Store:

```text
genome_id
catalogue_name
catalogue_release
```

## Taxonomy identifiers

Store:

```text
taxon_id
taxonomy_source
taxonomy_release
```

## Functional identifiers

Store:

```text
function_database
function_accession
database_release
```

---

# 13. Parquet Metadata

Every Parquet file should contain footer metadata.

Recommended keys:

```text
metaatlas_schema_version
metapilot_version
project_accession
view_name
creation_date
software_provider
compression_format
catalogue_name
catalogue_release
taxonomy_source
taxonomy_release
annotation_release
qpx_dataset_uri
```

---

# 14. Validation Rules

## File-level validation

- file is readable
- Parquet footer exists
- expected schema version exists
- expected view name exists
- project accession matches filename
- checksum can be calculated

## Schema validation

- all required columns exist
- no unexpected incompatible columns
- Arrow types are correct
- non-nullable columns do not contain null
- controlled-value fields contain accepted values

## Primary-key validation

- primary-key columns exist
- primary-key combinations are unique
- primary-key values are not null

## Foreign-key validation

Examples:

```text
genome_abundance.genome_id
    must exist in genome.genome_id

taxon_abundance.taxon_id
    must exist in taxon.taxon_id

function_abundance.function_accession
    must exist in function.function_accession

sample_accession
    must exist in sample.sample_accession
```

## Biological validation

- abundance is not negative unless the field explicitly permits signed values
- peptide counts are non-negative
- unique peptide count is not greater than total peptide count
- selected genomes are present in the genome table
- taxonomic ranks are valid
- catalogue release is present for all genomes
- annotation release is present for all functions

---

# 15. Export Configuration

Add a MetaPilot export configuration section.

Example:

```yaml
metaatlas_export:
  enabled: true
  schema_version: "1.0"
  project_accession: "PXD008738"
  output_directory: "MetaAtlasExport"
  qpx_dataset_uri: null

  catalogue:
    name: "MGnify Human Gut Genome Catalogue"
    release: "v2.0"
    checksum: null

  taxonomy:
    source: "GTDB"
    release: "R232"

  annotation:
    eggnog_release: "5.0.2"
    kegg_release: null
    go_release: null

  files:
    dataset: true
    sample: true
    genome: true
    taxon: true
    function: true
    peptide_genome: true
    genome_abundance: true
    taxon_abundance: true
    function_abundance: true
    genome_selection: true
    provenance: true

  parquet:
    compression: "ZSTD"
    row_group_size: 100000
    dictionary_encoding: true

  validation:
    enabled: true
    fail_on_error: true
    fail_on_warning: false
```

---

# 16. User Interface and CLI Integration

## GUI

Add an optional output section:

```text
Export MetaAtlas-compatible Parquet files
```

Suggested fields:

- Enable MetaAtlas export
- Project accession
- Catalogue name
- Catalogue release
- Taxonomy source
- Taxonomy release
- QPX dataset link
- Output directory
- Include genome-selection evidence
- Validate exported files

## CLI

Suggested flags:

```bash
--export-metaatlas
--metaatlas-project PXD008738
--metaatlas-output /path/to/export
--metaatlas-schema-version 1.0
--metaatlas-qpx-uri <URI>
--metaatlas-catalogue-name <NAME>
--metaatlas-catalogue-release <VERSION>
--metaatlas-taxonomy-source GTDB
--metaatlas-taxonomy-release R232
--metaatlas-validate
```

---

# 17. Logging

Suggested log messages:

```text
[MetaAtlas] Starting export for PXD008738
[MetaAtlas] Loading sample metadata
[MetaAtlas] Loading genome catalogue metadata
[MetaAtlas] Exporting genome view
[MetaAtlas] Exporting peptide-genome mappings
[MetaAtlas] Exporting genome abundance
[MetaAtlas] Exporting taxon abundance
[MetaAtlas] Exporting function abundance
[MetaAtlas] Building provenance
[MetaAtlas] Validating primary keys
[MetaAtlas] Validating foreign keys
[MetaAtlas] Writing manifest
[MetaAtlas] Export completed successfully
```

Errors should include the view name and offending identifier.

---

# 18. Testing Strategy

## Unit tests

Test:

- identifier normalization
- schema loading
- Arrow type conversion
- primary-key generation
- lineage serialization
- list and struct fields
- footer metadata
- checksum generation
- abundance validation

## Integration tests

Use one small DDA dataset and one small DIA dataset.

Validate:

- current MetaPilot results remain unchanged
- all expected Parquet files are generated
- all files can be read by PyArrow
- all files can be queried by DuckDB
- foreign keys resolve
- QPX identifiers can be retained
- results match existing TSV summaries

## Regression tests

For each exporter:

```text
existing TSV value
    equals
exported Parquet value
```

Examples:

- selected genome count
- peptide count per genome
- taxon abundance
- function abundance
- sample count
- run count

---

# 19. Recommended Implementation Phases

## Phase 1: Core export

Implement:

- `dataset.parquet`
- `sample.parquet`
- `genome.parquet`
- `peptide_genome.parquet`
- `genome_abundance.parquet`
- `provenance.parquet`
- manifest
- validation

Use one small DDA project.

## Phase 2: Biological aggregation

Add:

- `taxon.parquet`
- `taxon_abundance.parquet`
- `function.parquet`
- `function_abundance.parquet`

Validate against current MetaPilot reports.

## Phase 3: Genome-selection evidence

Add:

- `genome_selection.parquet`
- ranking
- selection reasons
- marginal gain
- knee information
- significance statistics

## Phase 4: QPX linking

Add:

- QPX project URI
- QPX peptide identifier
- QPX protein identifier
- quantms evidence links
- PRIDE raw-data links

## Phase 5: MetaAtlas ingestion

Build the MetaAtlas loader:

```text
Parquet files
      ↓
validation
      ↓
DuckDB / analytical warehouse
      ↓
MetaAtlas API
      ↓
website and cross-project queries
```

---

# 20. Example Cross-Project Queries

## Find projects containing a genome

```sql
SELECT DISTINCT
    project_accession,
    genome_id
FROM genome_abundance
WHERE genome_id = 'MGYG000123456'
  AND detected = TRUE;
```

## Find samples containing a taxon

```sql
SELECT
    project_accession,
    sample_accession,
    abundance
FROM taxon_abundance
WHERE taxon_id = '...'
  AND minimum_evidence_passed = TRUE;
```

## Find taxa contributing to a function

```sql
SELECT DISTINCT
    gt.taxon_id,
    pf.function_accession
FROM protein_function pf
JOIN protein_genome pg
  ON pf.protein_accession = pg.protein_accession
JOIN genome_taxon gt
  ON pg.genome_id = gt.genome_id
WHERE pf.function_database = 'KEGG_KO'
  AND pf.function_accession = 'K00001';
```

## Find genome-unique peptides

```sql
SELECT
    sequence,
    genome_id
FROM peptide_genome
WHERE is_genome_unique = TRUE;
```

---

# 21. Storage Strategy

MetaAtlas should avoid storing:

- PSM tables
- fragment-ion peaks
- raw spectra
- chromatograms
- duplicated QPX feature tables
- large intermediate search files

MetaAtlas should store:

- compact biological entities
- relationship tables
- sample-level abundances
- evidence counts
- provenance
- links to QPX and PRIDE

This keeps the repository substantially smaller than a full proteomics evidence repository.

---

# 22. Relationship with QPX

The custom MetaAtlas files should initially be described as:

> QPX-aligned MetaAtlas Parquet views

They should not be called official QPX files unless they are adopted by the QPX project.

The shared interoperability layer should include:

- project accession
- sample identifiers
- run identifiers
- peptide identifiers
- protein identifiers
- provenance references
- QPX dataset URI
- PRIDE accession

The division of responsibilities is:

```text
quantms/QPX:
mass-spectrometry evidence

MetaAtlas:
microbiome interpretation and cross-project knowledge
```

---

# 23. Definition of Completion

The upgrade is complete when MetaPilot can:

1. Run the existing DDA or DIA workflow without changes.
2. Generate the current standard result reports.
3. Generate MetaAtlas-compatible Parquet files.
4. Validate all exported files.
5. Produce a manifest with row counts and checksums.
6. Preserve full workflow and database provenance.
7. Link peptide and protein entities to QPX evidence where available.
8. Load the files into DuckDB without manual transformation.
9. Reproduce the current genome, taxon, and function summaries.
10. Support cross-project MetaAtlas queries.

---

# 24. Recommended First Pilot

Use one small DDA metaproteomics dataset.

The first pilot should demonstrate:

- FragPipe analysis through MetaPilot
- existing MetaPilot result generation
- QPX conversion for PSM, peptide, and protein evidence
- MetaAtlas export for genome, taxon, and function information
- schema validation
- DuckDB query
- links between MetaAtlas entities and quantms evidence

The first pilot should remain small enough that errors can be inspected manually.

After the DDA implementation is stable, apply the same export layer to a small DIA-NN project.

---

# 25. Summary

The proposed upgrade adds a new result-export layer to MetaPilot rather than redesigning its existing workflows.

The implementation should:

- preserve all current analysis behaviour
- convert finalized results into standardized Parquet views
- use QPX-aligned identifiers and provenance principles
- separate biological entities from observations
- avoid duplicating PSM and spectrum-level evidence
- support direct MetaAtlas ingestion
- enable future cross-project genome, taxon, and function queries

The resulting architecture will allow MetaPilot to serve as the analysis engine, quantms/QPX to preserve mass-spectrometry evidence, and MetaAtlas to provide the microbiome-specific knowledge and interpretation layer.

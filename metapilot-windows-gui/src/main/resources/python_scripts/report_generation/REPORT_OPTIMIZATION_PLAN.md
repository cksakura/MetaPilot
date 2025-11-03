## MetaPilot Report Generation: Optimization Plan

This document consolidates proposed optimizations for metaproteomics report generation across analysis, visualization, robustness, performance, packaging, and reproducibility.

### High-impact analysis additions
- **Differential abundance**
  - Add group-aware statistics for proteins, peptides, taxa, and functions (two-sample t/Welch, Wilcoxon; multi-group ANOVA/Kruskal). Apply Benjamini–Hochberg FDR.
  - Visualize volcano plots, ranked heatmaps, and summary tables of significant features with effect sizes and FDR.
- **Normalization and missing data**
  - Per-sample normalization options: median, quantile, VSN; plot pre/post intensity distributions (boxplot/violin) and MA plots.
  - Handle MNAR with left-censoring imputation (e.g., downshifted gaussian) and MCAR/MAR imputation options; show impact on distributions.
- **Correlation and QC**
  - Sample–sample Pearson/Spearman correlation heatmap with hierarchical clustering; replicate CV per group.
  - ID rate (PSMs/MSMS), TIC/PSMs per run, missed-cleavage rate by sample; optional mass error histogram and RT drift if present.
- **Dimension reduction and compositionality**
  - For taxonomy/function profiles, support CLR transform + Aitchison distances. Add Bray–Curtis and Jaccard where appropriate.
  - Add PCoA/NMDS and PERMANOVA to test group separations; annotate p-values on ordination plots.
- **Enrichment**
  - Functional enrichment (GO/KEGG/COG) via Fisher exact or preranked GSEA; bubble plots with FDR, gene count.
  - Compositional DA for taxa (ALDEx2/ANCOM-style methods; optional dependency or Python analogs).

### Visualization improvements
- **Taxonomy**
  - Replace heuristic sunburst parent mapping with explicit parent linkage using rank columns (e.g., `superkingdom → kingdom → phylum → ...`).
  - Per-sample stacked bars by rank with top-N + “Other”; heatmap of top taxa relative abundance across samples.
- **Functions**
  - Bubble charts for enriched terms (x: enrichment score; y: term; size: gene count; color: FDR).
  - Heatmap of top variable functions across samples/groups with clustering.
- **Proteins/Peptides**
  - MA plots, intensity boxplots pre/post normalization; PTM counts and localization confidence distributions if available.
- **Consistency and UX**
  - Unify color palette (color-blind-friendly) and legends; show explained variance on PCA/PCoA axes; optional group centroids/ellipses.
  - Create an index HTML that links to each sub-report and shows metadata summary and run parameters.

### Robustness and code quality
- **Column inference**
  - Generalize `extract_sample_info` to support MaxQuant (`Intensity`, `LFQ intensity`, `iBAQ`), FragPipe/Philosopher, DIA-NN, MSstats. Allow user-provided patterns/mappings.
- **Shared analytics utilities**
  - Move shared logic (log/scale, PCA/UMAP/PCoA, dendrograms, correlations, FDR correction) into `BaseReportGenerator` helpers for consistency.
- **Compositional transforms**
  - Add safe CLR helper (pseudo-count handling), and distance functions (Aitchison, Bray–Curtis).
- **Metadata handling**
  - Case-insensitive, trimmed matching; diagnostics for unmatched samples; miniature “metadata health” panel in each report.
- **Palette**
  - Expand `group_colors` to distinct, color-blind-safe hues (e.g., Okabe–Ito); avoid repeating same three tones.

### Performance and scalability
- **IO and memory**
  - Use categorical dtypes for taxonomy/function labels; optionally `engine="pyarrow"` or Polars for large TSVs.
  - Read only required columns per report; cache intermediate aggregates (e.g., `Total_Abundance`).
- **Computation**
  - Vectorize abundance/diversity calculations; reuse shared transformed matrices; guard dendrogram/UMAP for large N (switch to PCoA or sample down).
- **Output**
  - Offer static exports (PNG/SVG via `kaleido`) alongside HTML. Use consistent Plotly JS delivery (CDN or bundled) across all generators.
  - Save derived tables (top-N taxa/functions, DE results, normalization factors) as TSVs next to HTML.

### Packaging and ergonomics
- **CLI**
  - Provide `metapilot-report` CLI with options: data dir, metadata path, group column, normalization, imputation, distance metric, top-N, export formats, enable/disable sections.
- **Consolidation**
  - Deprecate standalone `report_generator.py` or wrap it to use the new `MetaPilotReportGenerator` to avoid divergence.
- **Reproducibility**
  - Centralize random seeds; record config and package versions in report footer and save a JSON manifest per run.

### Testing
- Unit tests for: column detection, metadata mapping, CLR transform, PCA/PCoA wrappers, diversity metrics, FDR correction, enrichment calculators.
- Golden tests for figure data specs (trace counts, axis labels) and table outputs.

### Dependencies (optional/conditional)
- `scikit-bio` (beta diversity, PCoA), `statsmodels` (multiple testing), `kaleido` (static export), `pyarrow`/`polars` (IO), `umap-learn` (if desired), `scikit-learn` (already used).

### Proposed initial implementation sequence
1) Unify palette, HTML export, and figure helpers in `BaseReportGenerator`.
2) Add correlation heatmap, normalization options, and pre/post QC plots (proteins/peptides).
3) Taxonomy/function compositional pipeline: CLR, Aitchison/Bray–Curtis, PCoA, PERMANOVA, stacked bars with top-N+Other.
4) Differential abundance and volcano plots with FDR for proteins/peptides/functions/taxa; save TSVs.
5) Enrichment bubble charts; group-wise comparison heatmaps.
6) CLI entrypoint and index HTML; config/manifest writing.
7) Performance: selective column reads, categorical dtypes, static exports.
8) Tests for the shared utilities and key outputs.

### Output artifacts to add
- HTML: per-section reports and an index page.
- PNG/SVG: selected plots for publication.
- TSV: top-N summaries, DE results, normalization factors, enrichment results.
- JSON: run manifest (versions, parameters, file hashes).

### Risks and constraints
- Large datasets: limit traces and taxon/function counts dynamically; offer sampling or rank-based truncation.
- Optional dependencies: design graceful fallbacks and clear annotations when features are unavailable.
- Multiple input formats: require a small mapping config when auto-detection is ambiguous.

— End of plan —


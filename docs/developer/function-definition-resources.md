# Function Definition Resources

Current version target: `1.1.1`.

MetaPilot uses `Resources/function/func_def.db` as the runtime lookup database for
functional annotation names and definitions. The maintained rebuild script is:

`metapilot-core/src/main/java/uk/ac/quadram/metapilot/core/function/build_function_resources.py`

The script replaces the older one-off download helpers in the same package. It
builds editable TSV resources under `annotation_maps/` and can also build the
SQLite database used by report generation and functional annotation export.

## Current Bundled Database

The current bundled database is:

`Resources/function/func_def.db`

Current local snapshot:

| Field | Value |
| --- | --- |
| Database size | 456,597,504 bytes |
| File timestamp UTC | 2026-05-16 09:43:04 |
| Builder metadata | `build_function_resources.py` |
| Build time UTC | 2026-05-16T09:41:31Z |
| Schema version | 1 |
| eggNOG OG annotations included | yes |
| Merged TSV snapshot | `annotation_maps/annotation_lookup_merged.tsv` |
| Merged TSV size | 11,578,931 bytes |
| eggNOG OG source snapshot | `annotation_maps/eggnog_ogs/e5.og_annotations.tsv` |
| eggNOG OG source size | 205,646,777 bytes |

SQLite row counts in the current database:

| Object | Rows |
| --- | ---: |
| `annotation_term` | 183,499 |
| `annotation_alias` | 54,013 |
| `annotation_property` | 16,882 |
| `eggnog_og_annotation` | 4,415,548 |
| `go_obo` view | 38,560 |
| `enzyme` view | 8,441 |
| `cog_20_def_tab` view | 5,061 |
| `kegg_pathway` view | 585 |
| `nog_e5_og_annotations` view | 4,415,548 |

`annotation_term` source counts:

| Source | Terms |
| --- | ---: |
| `BRITE` | 156 |
| `BiGG_Reaction` | 28,302 |
| `CAZy` | 3,000 |
| `COG` | 5,061 |
| `COG_category` | 30 |
| `EC` | 8,441 |
| `GO` | 38,560 |
| `KEGG_Module` | 570 |
| `KEGG_Pathway` | 585 |
| `KEGG_Reaction` | 12,385 |
| `KEGG_ko` | 28,239 |
| `KEGG_rclass` | 3,208 |
| `PFAM` | 54,962 |

## Sources Built By The Script

The rebuild script downloads or generates these source resources:

| Resource | Source behavior | Main output |
| --- | --- | --- |
| KEGG KO, pathway, module, reaction, rclass, enzyme, BRITE | KEGG REST `list` endpoints | `annotation_maps/kegg/*.tsv` |
| GO | `go-basic.obo` from OBO Library | `annotation_maps/go/GO.tsv`, `GO_full.tsv` |
| EC | ExPASy `enzyme.dat` | `annotation_maps/ec/EC.tsv`, `EC_full.tsv` |
| COG | NCBI COG2024 definition and category tables | `annotation_maps/cog/*.tsv` |
| Pfam | EBI Pfam current release metadata | `annotation_maps/pfam/PFAM.tsv` |
| BiGG | BiGG universal reactions API | `annotation_maps/bigg/BiGG_Reaction.tsv` |
| CAZy | Generated family labels for GH, GT, PL, CE, AA, and CBM | `annotation_maps/cazy/CAZy.tsv` |
| KEGG_TC | Placeholder table only | `annotation_maps/kegg_tc/KEGG_TC.tsv` |
| eggNOG OG | Manual input; not downloaded automatically | `annotation_maps/eggnog_ogs/e5.og_annotations.tsv` |

The full eggNOG OG annotation table is intentionally manual because it is large.
Provide it with `--eggnog-og-annotations` when the bundled database should
include the legacy-compatible `nog_e5_og_annotations` view.

## Rebuild Commands

Use this from the repository root when all remote sources should be refreshed:

```powershell
python metapilot-core/src/main/java/uk/ac/quadram/metapilot/core/function/build_function_resources.py `
  --outdir annotation_maps `
  --sqlite-out Resources/function/func_def.db `
  --eggnog-og-annotations annotation_maps/eggnog_ogs/e5.og_annotations.tsv `
  --force
```

Use this when the TSV/raw files in `annotation_maps/` are already current and
only the merged lookup plus SQLite database need rebuilding:

```powershell
python metapilot-core/src/main/java/uk/ac/quadram/metapilot/core/function/build_function_resources.py `
  --outdir annotation_maps `
  --sqlite-out Resources/function/func_def.db `
  --eggnog-og-annotations annotation_maps/eggnog_ogs/e5.og_annotations.tsv `
  --skip-downloads
```

If eggNOG OG annotations are unavailable, omit `--eggnog-og-annotations`. The
database will still contain `annotation_term`, `annotation_alias`,
`annotation_property`, and the GO/EC/COG/KEGG views, but it will not contain
`eggnog_og_annotation` or `nog_e5_og_annotations`.

## SQLite Schema

The current schema is intentionally small and compatibility-focused:

| Object | Purpose |
| --- | --- |
| `metadata` | Builder name, build time, schema version, and eggNOG inclusion flag |
| `annotation_term` | Canonical function terms by source and identifier |
| `annotation_property` | Extra key/value properties, currently EC reactions and synonyms |
| `annotation_alias` | Source-specific aliases such as prefixed KEGG IDs |
| `eggnog_og_annotation` | Optional manual eggNOG OG annotations |
| `go_obo` | Legacy-compatible GO view |
| `enzyme` | Legacy-compatible EC/enzyme view |
| `cog_20_def_tab` | Legacy-compatible COG view |
| `kegg_pathway` | Legacy-compatible KEGG pathway view |
| `nog_e5_og_annotations` | Legacy-compatible eggNOG OG view when eggNOG is included |

`IGCFuncSqliteSearcher` detects the new `annotation_term` table and uses the
normalized schema first, while the compatibility views keep older query paths
working.

## Validation After Rebuild

After rebuilding, verify the database before packaging:

```powershell
python -c "import sqlite3; con=sqlite3.connect('Resources/function/func_def.db'); cur=con.cursor(); print(cur.execute('select count(*) from annotation_term').fetchone()[0]); print(cur.execute('select value from metadata where key=\"eggnog_ogs_included\"').fetchone()[0]); con.close()"
```

Expected for the current snapshot:

- `annotation_term` is about 183k rows.
- `eggnog_ogs_included` is `true` when the manual eggNOG source is supplied.
- `nog_e5_og_annotations` exists only when eggNOG is included.

Do not package a regenerated database if public source downloads failed unless
the failure is intentional and documented in this file.

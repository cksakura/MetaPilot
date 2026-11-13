# Open search: what is decided for each spectrum, and where each level is written

Status: 2026-11-10. Describes the code as it stands (uncommitted), not a released version. The
entrapment run and the end-to-end run of this build have not been done yet.

All paths are under the run's `opensearch_MAG_search/` folder.

## 1. The decision for one spectrum

The closed and the open search run against the same per-sample database. Their candidates are merged,
scored by one model, and each precursor gets one winner. The winner falls into one group, and each
group has its own 1% FDR walk at spectrum level and at peptide level.

| The winner | Group | Main result |
|---|---|---|
| Found by the closed search (exact mass, with its variable modifications) | `EXACT` | yes |
| Open search only, on the zero peak | `OPEN_UNMODIFIED` | yes |
| Open search only, on an accepted mass peak, and it passes the per-spectrum rule below | `OPEN_ACCEPTED_SHIFT` | yes |
| Open search only, off any accepted peak, fragments place the shift on one residue, UNIMOD explains it there and no substitution fits | `OPEN_LOCALISED` | yes |
| Fragments place the shift on a residue where an amino acid substitution fits | `OPEN_SUBSTITUTION_SITE` | no, set aside |
| Anything else | `OPEN_UNEXPLAINED` | no, set aside |

### The per-spectrum rule on an accepted peak

A peak is accepted as a whole first (evidence rule, site share, fragment test). Each of its winners is
then judged by itself:

| What the fragments show (`shift_state`) | What happens |
|---|---|
| Shift on a residue UNIMOD allows, no substitution fits (`localised_explained`, `shift_class = modification`) | Stays. Exported as the peptide with that modification on that residue. |
| Shift on a residue where a substitution fits (`localised_explained`, class `substitution` or `modification_or_substitution`) | Set aside as `OPEN_SUBSTITUTION_SITE`. |
| Shift on a residue nothing explains (`localised_unexplained`) | Set aside as `OPEN_UNEXPLAINED`. |
| Shift not on the fragments, and UNIMOD has a modification of the peak's mass on a residue or terminus of the peptide (`assigned`) | Stays. Exported and quantified as that modified peptide; the site is UNIMOD's (common specificity before rare), not measured. |
| Shift not on the fragments, no such UNIMOD site (`unsupported`) | Set aside as `OPEN_UNEXPLAINED`. |

No mass range is treated specially. UNIMOD's isotopic-label entries and its amino acid substitution
entries are not modifications for this purpose.

### The three classes

| Class (`shift_class`) | Meaning | Where it goes |
|---|---|---|
| `modification` | Only a UNIMOD modification explains the shift at its site | Main result |
| `substitution` | Only an amino acid substitution explains it | Set aside; variant outputs and side quantification |
| `modification_or_substitution` | Both fit the same residue; the fragments are identical | Set aside; no variant call |

## 2. Files by level

### Spectrum (PSM)

| File | Content |
|---|---|
| `combined/validated_psms.tsv` | The main result: every identified PSM, closed and open. Input of FlashLFQ. |
| `combined/unexplained_psms.tsv` | The set-aside PSMs (`OPEN_SUBSTITUTION_SITE`, `OPEN_UNEXPLAINED`) that pass their group's 1% walk. Not identifications. |
| `combined/peak_validated_psms.tsv` | The open part of the main result. |
| `combined/scored_psms.tsv` | Every candidate of every spectrum with its score, group, winner flag and q-values. Read by the genome test. |

Shift columns in these files: `id_group`, `shift_state`, `shift_site` (1-based), `shift_site_probability`
(empty for `assigned`), `shift_explanation`, `shift_class`.

### Peptide

| File | Content |
|---|---|
| `combined/QuantifiedPeptides.tsv` | FlashLFQ intensities of the main result. |
| `final_peptides.tsv` | Peptides of the main result with intensities and razor protein; last columns `Shift explanation`, `Shift state`, `Shift class`. |
| `combined/set_aside_peptides.tsv` | The set-aside PSMs per base peptide, reason and shift mass. |
| `excluded_evidence/peptides.tsv` | Main-result peptides with no confident genome (kept out of protein and genome tables). |

### Protein

| File | Content |
|---|---|
| `razor_proteins.tsv`, `all_proteins.tsv` | Built from the main result only. |

### Genome and species

| File | Content |
|---|---|
| `taxonomy_analysis/` | Built from the main result. The genome test counts `EXACT`, `OPEN_UNMODIFIED` and `OPEN_ACCEPTED_SHIFT` peptides as evidence; set-aside groups and `OPEN_LOCALISED` are not evidence. |
| `combined/set_aside_genomes.tsv` | Per genome: set-aside PSMs and peptides by reason beside the genome's identified PSMs and peptides. A genome whose strain differs from the reference loses more; this file shows how much. A peptide shared by several genomes counts for each. |

### Function

| File | Content |
|---|---|
| `functional_annotation/` | Built from the main-result proteins; written only when the function database is configured. |

### Substitutions (`sequence_variants/`)

| File | Content |
|---|---|
| `substitution_forms.tsv` | One row per PSM with a substitution verdict: set aside, and a substitution is the only explanation at its site. |
| `sequence_variants.tsv`, `.fasta` | Those PSMs whose substituted peptide is in no searched database, grouped into sites and classed by recurrence, codon reachability and whether the reference peptide is identified unmodified anywhere in the main result (closed search included). |
| `substitution_psms.tsv`, `substitution_quantified_peptides.tsv` | Side quantification of the substitution PSMs beside their reference peptides, in a FlashLFQ run of its own. |
| `variant_candidates.tsv`, `variant_sites.tsv`, `variant_pangenome_summary.tsv`, `variant_species_summary.tsv` | Lookup of candidate peptides in UHGP-100 (MGnify human-gut v2.0.2 protein catalogue, all genomes), i.e. in genomes not selected for the sample. A candidate comes from an `OPEN_SUBSTITUTION_SITE` PSM at the fragments' site; `ms2_decidable = false` marks one where a modification fits too. Report-only. |
| `substitution_peaks.tsv`, `substitution_candidates.fasta` | The older peak-level analysis made before scoring. Diagnostic; nothing downstream follows it. |

### Peak-level audit (`combined/`)

`delta_mass_peaks.tsv`, `modification_verdicts.tsv`, `peak_site_shares.tsv`, `peak_fragment_test.tsv`,
`potential_modifications.tsv`, `joint_identification_summary.tsv`, `identification_mode.txt`.

## 3. Known limits

- A shift the fragments do place is explained from that PSM's own measured mass within 0.005 Da. A PSM
  with a larger mass error can be set aside although its peak's mass is in UNIMOD.
- A shift made of two modifications, or of a modification plus an isotope step (the +1.95 to +1.98 Da
  peaks), has no single UNIMOD entry and is set aside.
- The +4 to +5 Da band is set aside like any mass UNIMOD does not have, although an MS1 check on 12mix
  (2026-11-07, 3 of 9 files) showed these are mostly unmodified peptides with a wrong precursor mass.
- `OPEN_LOCALISED` PSMs are quantified but are not genome evidence. On the round-5 entrapment runs the
  group added no peptide of its own.
- "Reference also observed" for a variant counts unmodified identifications only, in any sample. A
  reference peptide shared with another detected genome is not distinguished.
- The parquet export (`metapilot.psm.parquet`) holds the main result only and has no `shift_class`.

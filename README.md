# MetaPilot

Source repository for MetaPilot 1.2.0, a metaproteomics workflow for DDA and DIA. MAG-aware workflows select genomes from a catalogue, build a sample-specific database or library, and search that reduced space. The user install and the HPC getting-started guide live in the release repository. This file explains how the workflows operate.

The version under development is 1.2.0. User-facing changes are listed in [CHANGELOG.md](CHANGELOG.md). The open-search design history is in [docs/developer/opensearch-changelog.md](docs/developer/opensearch-changelog.md).

## Modules

| Module | Role |
| --- | --- |
| `metapilot-core` | Workflows, genome selection, open-search identification, result export |
| `metapilot-cli` | Headless runner used on Linux and HPC |
| `metapilot-windows-gui` | Windows desktop application |

External programs do the spectrum search and quantification. MetaPilot writes their parameter files, reads the results back, and decides which genomes, peptides, and mass shifts to keep. For open search those programs are MSFragger, Crystal-C, MSBooster, DIA-NN, and FlashLFQ.

## MAG genome selection

Every MAG workflow starts from a registered catalogue (formal runs use human-gut v2.0.2). A first search is intentionally small: marker proteins, or an existing DIA library. Peptides from that search are assigned to genomes. Genomes are ranked by the new marker peptides they explain, and selection stops by one of two rules.

- **Focus** keeps the high-gain genomes and stops at the knee of the coverage curve.
- **Comprehensive** keeps adding genomes further down the same ranking.
- **PSM coverage** keeps genomes until a chosen share of the marker PSMs is explained. The ranking itself is unchanged.

The selected genomes become one sample-specific FASTA, or one refined DIA library. A later search uses only that database. A reversed genome that passes the same test is reported beside the accepted target genomes, which is the genome-level false-discovery check.

## DDA open search

`mag-opensearch` is the DDA MAG workflow that does not require the mass shifts to be named in advance. It has two stages. Stage 1 chooses the genomes. Stage 2 identifies peptides, including modified and variant forms, against the selected genomes.

### Stage 1: marker search and rescoring

Each raw file is searched against the target-decoy marker database (`Marker4_td.fasta`). The default marker search is closed: carbamidomethyl cysteine is fixed, and oxidation of methionine and protein N-terminal acetylation are variable. MSBooster and DIA-NN add predicted-spectrum and retention-time features.

A semi-supervised logistic model then rescores the marker PSMs. The model is trained per file at 1% FDR. Its ridge is 0.1 and the ranking score is the logit. A near-zero ridge lets the fit diverge when a marker search has only a few hundred confident PSMs: target probabilities collapse to 1, decoys tie with targets, and the file keeps no PSMs, so no sample database is built. Ridge 0.1 keeps at least the search engine's own accepted set on those files.

Only unmodified PSMs, those whose mass difference sits in the zero-peak window, are passed to genome selection. The selected genomes are written to `marker/<sample>.fasta`. An optional host proteome can be appended.

### Stage 2: closed and open search of the sample database

Each sample FASTA is searched twice with MSFragger.

- The closed search allows the same variable modifications as the marker search.
- The open search uses a wide precursor window (by default −150 to +500 Da) and lets MSFragger localise the mass difference.

Crystal-C refines the open-search precursor masses. Candidates from both searches are then classified together.

### Joint identification

Joint identification (version `v2`) trains one XGBoost model on closed-search and open-search features, including the predicted-spectrum features. RESET walks the scored list at 1% FDR separately for each evidence group, so a modified population is not held to the error rate of the unmodified one. Decoys are pseudo-reversed.

A mass shift ends in one of three classes.

- **Zero.** The precursor matches the unmodified peptide inside the calibration window.
- **Established.** The shift matches a known modification, a residue substitution, or another accepted explanation, and the fragment evidence supports it.
- **Unexplained.** The shift is reproducible but has no accepted explanation. These PSMs are reported apart from the validated set.

An isotope guard stops a real modification from being absorbed as a satellite of a neighbouring peak. Acetyl (42.0106 Da) is the case that motivated it: it can sit on the −1 isotope position of carbamyl.

Sequence variants are called from the substitution mass shifts. When a UHGP-100 catalogue is configured, each variant peptide is checked against that pangenome and against matched controls. The check reports whether the variant is already known in the species, including outside the selected genome.

FlashLFQ quantifies the validated peptides. Taxonomy and functional annotation are exported from the accepted genomes and proteins. A genome confidence gate can move weak genomes, and the peptides that depend on them, into `excluded_evidence/`.

### What is reused when a run is repeated

Each step is skipped when its output already matches the current configuration.

| Existing output | What is skipped |
| --- | --- |
| `marker/<sample>.pepXML` whose `fragger_open.params` matches the text that would be written now | The marker search |
| MSBooster output newer than that pepXML | MSBooster |
| `marker/<sample>.fasta` | Genome selection |
| `sample_db/<sample>.pepXML` newer than the sample FASTA | The open stage-2 search |

The sample FASTA is the one to remove when stage 1 itself has changed. Deleting `marker/<sample>.fasta` makes genome selection run again, and the stage-2 searches then rerun because they are older than the new database. Leaving the FASTA in place keeps the genomes chosen by the previous stage 1.

MSFragger writes `<sample>.mzBIN_calibrated` beside the mzML. Two jobs that search the same spectra will overwrite that file.

### Entrapment evaluation

`mag-opensearch-entrapment` uses the same marker search and genome selection, then adds a paired shuffled-protein population to the stage-2 database. It reports separate closed, open, and joint results for the real and the shuffled proteins and stops before quantification and taxonomy. The summary is descriptive. It does not claim a calibrated false-discovery proportion. Results go to `opensearch_entrapment_MAG_search/` and do not replace a normal open-search run.

## Result folder

A normal open-search run writes `<output_directory>/opensearch_MAG_search/`.

| Path | Contents |
| --- | --- |
| `search_results/<sample>/marker/` | Marker pepXML, pins, rescoring output, selected-genome FASTA |
| `search_results/<sample>/sample_db/` | Open stage-2 search and Crystal-C output |
| `search_results/<sample>/sample_db_closed/` | Closed stage-2 search |
| `combined/validated_psms.tsv` | PSMs accepted by joint identification |
| `combined/unexplained_psms.tsv` | Reproducible shifts with no accepted explanation |
| `combined/joint_identification_summary.tsv` | Counts per evidence group |
| `final_summary.tsv` | Per-file MS/MS identification rate, peptide count, and PSM count |
| `final_peptides.tsv` | Peptide table, including the shift class |
| `razor_proteins.tsv` | Protein groups after the genome confidence gate |
| `sequence_variants/` | Substitution calls and the UHGP pangenome check |
| `taxonomy_analysis/` | Genome and taxon tables |
| `functional_annotation/` | eggNOG, KEGG, GO, EC, and COG summaries |

`final_summary.tsv` counts MS/MS identifications only. A peptide that FlashLFQ fills in by match-between-runs is quantified, and it is not counted as identified in that file.

## Other workflows

DIA Marker, DIA Seed, DIA Target, and DIA Real-Time follow the same genome-selection idea, with DIA-NN and a spectral library in place of the DDA stage-2 search. DDA MAG FragPipe and DDA MAG pFind are the closed-search forms of the DDA MAG workflow: marker search, genome selection, sample FASTA, then one closed search. DDA IGC searches a catalogue or a user database directly and does not select genomes.

Closed DDA and the DIA modes are the right comparison when the question is what the open search adds. They do not discover a mass shift that was not declared as a variable modification.

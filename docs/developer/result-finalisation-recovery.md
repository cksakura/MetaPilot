# Genome protein peptide finalisation and recovery

The F1 catalogue-wide remapping is confined to final export in the built-in pipelines. It does not feed subsequent search databases. Searches do not need repeating to correct F1–F4; reconstruct the final outputs from retained original search and mapping inputs.

## Feedback trace

- `DiaMarkerHeadlessRunner` builds the selected-genome FASTA/library and completes `combined.parquet` and its PR matrix before invoking `DiaMagPostAnalysisTsvExporter`.
- `DiaSeedHeadlessRunner` selects genomes from `seed.parquet`, constructs its selected-genome or Deep library, performs final quantification, and only then invokes final export. Direct searches its input Seed library before final export.
- `DiaMarkerDeepTogetherHeadlessRunner` and `DeepIterativeCommonSupport` build subsequent rounds from recorded genome selections and previous **round parquet** evidence. Their final combined library also comes from round evidence. They do not read `final_peptides.tsv`, `razor_proteins.tsv`, or `all_proteins.tsv`.
- GUI Seed/Marker export likewise follows identification and final quantification. Seed's initial catalogue mapping is intentional selection evidence, separate from the erroneous CLI export remapping.
- Repository consumers of final TSVs are reporting, entrapment evaluation, output checks and cleanup; there is no built-in search database construction from them.

This trace establishes absence of automatic feedback in the inspected source. It does not certify arbitrary external scripts or a user-supplied Seed library manually built from an earlier final export.

## Recovery scope

The audit's 92 archived result folders were inspected for standard search-output availability. Ninety contain their standard inputs. Two nested historical backups (`Elo_12mix_DIA/stale/marker_comprehensive/old` and `Elo_12mix_DIA/stale/seed_direct/old`) hold reports only. Their parent folders retain combined and first-pass search outputs; match the original generation before choosing either as the recovery source. File availability is not a hash/provenance validation. The per-folder inventory is in `result-finalisation-recovery-inventory.tsv`.

| Workflow | Recovery inputs | Mapping before finalisation |
|---|---|---|
| Marker Focus/Comprehensive | Original combined parquet and PR matrix | DIA-NN `Protein.Ids` |
| Marker Deep / Deep per file | Combined parquet, PR matrix, round parquets or exported round peptide/protein map | Recorded round associations |
| Seed Focus/Comprehensive | Combined parquet, PR matrix, original `genome_selected.tsv`, matching predicted TSVs | Only those selected genomes |
| Seed Deep | Combined parquet, PR matrix, completed visited-round mappings | Recorded visited rounds; final reporting is pruned again |
| Seed Direct | Combined parquet, PR matrix, matching catalogue predicted TSVs | Catalogue-wide candidates, then confidence gate and final ranking |
| DDA MAG | Original quantified peptide and PSM tables plus matching configuration | Original search associations |

All post-19-October affected DIA reruns, earlier DIA results affected by final selection/gating, Seed Direct, DDA association inconsistencies and the ongoing faecal Seed jobs fit post-analysis recovery. No running search was stopped or resubmitted by this change. Retain search artefacts until post-analysis is validated.

## Finalisation contract

Candidate scoring retains the existing statistical method. The genome confidence gate is q ≤ 0.01 (or an explicitly supplied accepted set for score-free peptide tables). Provisional protein/peptide assignments supply ranking inputs. Final selection uses the existing privileged/greedy ranking; zero-increment leftovers are excluded. The selected set is then frozen.

All microbial associations are restricted to that set. Razor proteins are reassigned within retained associations, and quantities, peptide/protein tables, summaries, taxonomy and functional inputs are rebuilt. Original genome statistics are preserved rather than recalibrated after pruning. Taxonomy receives the frozen decision and cannot select a second, different genome set.

GUI and CLI DIA exports now call the same core exporter. Seed mapping is explicit; an unrelated `predicted/` directory cannot change Marker output. Final-export Seed mapping preserves the original decoy background without changing upstream Seed selection. Missing Deep round mappings fail rather than reverting to an arbitrary combined-library accession.

`excluded_evidence/peptides.tsv` contains excluded sequence identities, associations and reasons, outside final result counts and quantities. `taxonomy_analysis/genome_final_decisions.tsv` records gate and ranking exclusions. Ranking files retain provisional evidence; final quantitative tables use rebuilt quantities. Empty final selections overwrite ranking and taxonomy tables with empty outputs, rather than retaining earlier results.

## Re-export commands

Use `diann-mag-report` with the original parameters and catalogue registry. Write to a fresh output directory; provide original inputs explicitly. This command performs no search.

```text
# Marker Focus/Comprehensive
java -jar <fixed-cli.jar> diann-mag-report --params <params.json> --magdb <registry.json> --out <new-output> --combined-parquet <original/combined.parquet> --pr-matrix <original/combined.pr_matrix.tsv>

# Deep (Marker or Seed): use the complete recorded map for that run
java -jar <fixed-cli.jar> diann-mag-report --params <params.json> --magdb <registry.json> --out <new-output> --combined-parquet <original/combined.parquet> --pr-matrix <original/combined.pr_matrix.tsv> --protein-mapping-policy ROUND_PROTEINS --peptide-protein-map <recorded-round-map.tsv>

# Seed Focus/Comprehensive
java -jar <fixed-cli.jar> diann-mag-report --params <params.json> --magdb <registry.json> --out <new-output> --combined-parquet <original/combined.parquet> --pr-matrix <original/combined.pr_matrix.tsv> --protein-mapping-policy SEED_SELECTED --selected-genomes <original/genome_selected.tsv>

# Seed Direct
java -jar <fixed-cli.jar> diann-mag-report --params <params.json> --magdb <registry.json> --out <new-output> --combined-parquet <original/combined.parquet> --pr-matrix <original/combined.pr_matrix.tsv> --protein-mapping-policy SEED_DIRECT
```

The recorded host library must remain reachable from the original combined directory's `host_library.txt` when used. Keep catalogue/prediction versions identical to the original inputs. A missing mapping resource requires recovering that resource, not automatically rerunning a search. The deprecated incomplete-Deep flag no longer bypasses required provenance.

## Strict audit

The version-controlled script is `scripts/12_release/list_consistency_audit.py`; the publication archive carries the same script. It checks every mapped protein, razor membership, both protein tables, reverse peptide links, unknown IDs and decoys. Legacy reports without a razor column are explicitly incomplete. It never invents a razor from the first mapped protein or assumes an unknown ID is host.

New exports record catalogue identity in `list_consistency_context.json`. For older exports supply the registry and catalogue ID. Host accessions must come from the recorded host FASTA(s) or a complete accession list recorded in that context. Mixtures use explicit `--no-host`. Without host policy the audit cannot report PASS.

```text
python scripts/12_release/list_consistency_audit.py --root <archive-root> --result-dir <result-folder> --catalogue-config <registry.json> --catalogue-id human-gut --no-host
python scripts/12_release/list_consistency_audit.py --root <archive-root> --result-dir <faecal-result-folder> --catalogue-config <registry.json> --catalogue-id human-gut --host-fasta <recorded-human.fasta>
```

Summary counts always cover all rows. The details file defaults to at most 10,000 examples per result to control memory/disk on catalogue-expanded exports; `details_truncated` marks this, and `--max-details-per-result 0` requests every violation. Non-PASS results return exit code 1. Use `--output` to preserve the original audit during comparisons.

## Validation

Regression coverage includes gate/ranking separation, shared-peptide reassignment, quantities after pruning, DDA exclusions and summaries, Seed selected versus catalogue scope, preservation of export decoys, strict Deep provenance, and DIA rejection of a negative genome with empirical p < 1 but q > 0.01. The shared DIA exporter is exercised against a synthetic PR matrix and taxonomy reference.

Read-only archive checks: Elo 12mix Marker Focus passes (45,880 peptides); Seed Direct reproduces 6,293 peptides with a razor genome absent from the final genome list. The stricter DDA check finds 28 peptides without a razor and 478 `all_proteins.tsv` entries from unreported genomes in Elo 12mix Marker4 Comprehensive. This expands the earlier audit's DDA scope beyond F4.

Internal consistency does not establish error-rate calibration. Recompute the entrapment readout from regenerated outputs and retain the revision plan's independent error-control validation.


Validation status for this implementation: 16 focused Java tests and 12 Python audit tests pass; the core, GUI and CLI compile together. A broader core run executed 261 tests with 7 failures, 8 errors and 1 skip. Its failing paths include sandbox temporary-file moves, Windows registry access, and an open-search canonical-key assertion outside finalisation. The full suite is not a clean release gate. No biological result folder has been re-exported, and publication figures remain unchanged.

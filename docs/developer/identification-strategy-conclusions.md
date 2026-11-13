# Closed + open identification: established conclusions

Status: settled as of 2026-10-23. These points are not re-tested; new work builds on them.

## Evaluation basis

- **Deciding benchmark:** the fecal open-search entrapment run (IPX0003851004, `opensearch_entrapment_MAG_search`).
  - 1 shuffled protein per real protein (MGYG8..., seed 1), with cleavage sites kept.
  - Labels are resolved by sequence against each run's searched database, not by accession.
  - FDP = 2E / (T + E).
- **12mix** is sanity only (a mock community; no tuning on it). Replication between runs is supporting evidence only.

## Conclusions

1. **Cascade beats joint.** It does better at PSM level at matched error, and at peptide level.

   | | PSMs at PSM q 1% | Final peptides |
   |---|---|---|
   | Cascade | 499,132 (FDP 1.13%) | 64,331 (FDP 1.20%) |
   | Joint | 486,433 (FDP 1.12%) | 58,998 (FDP 1.33%) |

2. **PSM-level FDR is calibrated in every arm.** Measured FDP is about 1.1-1.2 times the nominal q.

3. **Peptide-level FDR is the bottleneck.**
   - Cascade and joint both exceed 1% at peptide level. Their decoy estimates are 0.91% and 0.99%, against measured 1.20% and 1.33%.
   - Any PSM-level gain is lost there until the peptide step is fixed.

4. **Predicted spectra (DIA-NN) help the closed classifier.**
   - Closed stage: 379,331 PSMs at FDP 1.09%, against 368,068 at 1.14% (+3%).
   - In the full cascade the extra peptides come at about 4% marginal FDP, so peptide FDP rises to 1.32% (see 3).
   - Keep the predictions in the closed classifier; fix the peptide level.

5. **Re-ranking the closed top-3 has no effect.** The rescored winner differs from closed rank 1 in 998 of 760k spectra (0.1%). Drop it.

6. **Crystal-C improves the open search.** It stays.

7. **Closed and open "conflicts" are mostly different precursors.**
   - Crystal-C reassigns the precursor mass or charge, mostly to a co-isolated precursor. 98% of the open-unmodified answers that differ from an accepted closed answer are for such a reassigned precursor, not rivals for the same one.
   - Both sides measure similar error in these spectra: closed 1.66%, open 1.22%.
   - These open answers carry about 4 times the error of open answers that agree with closed (1.26% against 0.31%), yet they are scored as confidently. No feature tells the score that the precursor differs.
   - **Consequence:** compete only candidates that explain the same precursor. An open PSM for a reassigned precursor is a separate hypothesis, not an alternative to the closed answer.

8. **Unified and stage-wise methods.**
   - On 12mix both stay below the cascade (stage-wise: 76,628 peptides; cascade: 90,963).
   - The cause is the exact-stage model's own features, which are weaker than the engine features. On the closed-only check they gave 74,839 peptides against 94,046 with the engine features.
   - Recurrence counts in the PSM score about double the retention-time-inconsistent share, so they are off by default.
   - Neither method has been run on the entrapment data.

## Code state

All new options are off by default.

| Item | What it does |
|---|---|
| `cascade_ms2_features` | DIA-NN features in the cascade's closed and open classifiers |
| `cascade_closed_rerank` | closed re-ranking of the top-3 (see 5) |
| `SemiSupervisedPsmClassifier.Config.unclampedLogits` | unclamped log-odds for the shared classifier |
| `--stagewise`, `--psm-recurrence` | unified CLI flags |

## Next

1. Fix the peptide-level FDR (conclusion 3).
2. Make the competition precursor-aware (conclusion 7).

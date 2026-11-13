# Identification principles

Fixed principles for MetaPilot's identification. They were set by the project owner and are not
re-tested or argued again. A result that seems to conflict with one is reported as a result; the
principle stays.

## P1. Rescore (2026-10-23)

Every identification decision is made on a rescored value, never on a raw search score:
- which candidate wins a spectrum;
- which class a candidate falls into;
- whether it is accepted.

A machine-learning model that captures the relevant features scores better than any raw score
(hyperscore, expect, delta score, or a hand-weighted rule built from them).

- Raw scores and search properties are **features**, not decisions. These include the hyperscore,
  each search's own expect, the search-space size, the mass class or shift, and the
  predicted-spectrum agreement.
- Improvement comes from **better features and a better model**, not from hand-made rules. A
  pattern seen in the data, such as unmodified answers being more likely a priori, goes to the
  model as a feature. It is not hard-coded as a threshold or weight.

## P2. No q values across searches for the winner (2026-10-23)

q values from different searches are not comparable, so they are never used to decide between
candidates from different searches. The winner comes from the rescoring model (P1).

## P3. The precursor defines the competition (2026-10-23)

Candidates compete only within a precursor query.
- A Crystal-C reassigned precursor is its own query and its own FDR group.
- A scan may yield one identification per precursor.

## P4. Training decoys are never estimating decoys (2026-10-23)

Decoys used to train a model are never used to estimate its FDR. This follows RESET, with the
split made per decoy peptide. It holds whatever the model.

## P5. Folds by observed scan (2026-10-23)

Cross-validation folds are made by observed scan, and all of a scan's queries go into the same
fold.

## P6. Judge ranking and calibration, PSM and peptide (2026-10-23)

Every method is judged on the final exported set, at PSM and at peptide level, in two ways:
- **ranking:** yield at matched entrapment FDP;
- **calibration:** entrapment FDP at a decoy-estimated 1%.

Additional identifications per scan (beyond the first) are judged separately in both ways, as well
as within the combined sets (extended 2026-10-23).

Every level is reported with both readings side by side: PSM, peptide, genome, the one-peptide
protein group, all proteins, and the peptides retained after the protein filter (extended
2026-10-23).
- **FDR:** the decoy (`rev_`) estimate, D/T.
- **FDP:** the entrapment measurement, 2E/(R+E).
- Each comes with its counts (targets, entrapment, decoys) and a 95% interval. A difference
  between arms counts only when the intervals support it, not by eye.

## P7. Calibration is a gate, yield comes after (2026-10-23)

A model or method is used only when its decoy estimate is calibrated at every reported level. This
means the entrapment FDP agrees with the decoy FDR within its interval. A method that ranks better
but whose decoy estimate is optimistic is not used: its threshold cannot be set without entrapment.
Its matched-FDP advantage stays a reason to look for the cause, not a reason to adopt it.

## P8. Protein-level filter (2026-10-23)

- Proteins with one razor peptide form their own FDR group, walked to 1% decoy FDR on their
  peptide's score.
- Proteins with at least two razor peptides are kept.
- Peptides whose razor protein is rejected are dropped, so the peptide and protein outputs agree.

## P9. After a genome filter, entrapment is twin-matched (2026-10-23)

Entrapment genomes are never present, so a genome filter removes entrapment peptides whereas
false peptides on present real genomes survive. Measured naively after the filter, the FDP is
therefore far too low. Every entrapment measurement after a genome-level filter treats the
entrapment twin of each accepted real genome as accepted too, and the decoy side uses the reversal
of the same genomes.

## P10. Closed and open search are combined, not chosen between (2026-10-24)

Open search is not a replacement for closed search. It is complementary: each finds identifications
the other misses. The question is always how to combine the two, never which one to use. Both enter
one rescoring model (P1), and each candidate competes within its precursor query (P3).

## P11. Decoys must be able to make the errors targets make (2026-10-25)

A decoy has to be able to win a spectrum the way a false target can. Each decoy peptide therefore
has the length, residue composition and precursor mass of a target peptide. The stage-2 databases
use cleavage-preserving pseudo-reversal for this.

Whole-protein reversal moves the cleavage sites, so its decoys cannot be same-mass neighbours of a
present peptide, while false targets can. Under reversal, the entrapment-to-decoy ratio of twin
errors (2E/wD) was:

| Run | Reversed decoys | Pseudo-reversed decoys |
|---|---|---|
| fecal | 1.69 | 1.21 |
| 12mix | 2.06 | 1.23 |

## P12. Features are label-blind and exchangeable (2026-10-25)

- Every feature is computed in the same way for target, decoy and entrapment candidates.
- Any choice a feature depends on, such as a rival, a modification site or a peak window, is made
  by a label-blind score. It is never made by the model's own score.
- A new feature is checked for exchangeability: its distribution must be compared between
  exported decoy winners and exported entrapment winners. A clear difference is a calibration risk
  even when the calibration gate (P7) passes.

## P13. What gets adopted (2026-10-25)

- **Mock first.** A change is tried as a mock on cached outputs first, with the audit harness and
  the same export chain. A new search or HPC run follows only when the mock passes.
- **The deciding benchmark is the fecal entrapment run.** The 12mix run is a sanity check only:
  - no level may drop by more than 1%;
  - all 12 mix genomes must be accepted.

  Nothing is tuned on 12mix.
- **A change must earn its place.** After the calibration gate (P7), it is adopted only when it
  adds at least 1% real peptides at matched entrapment FDP on the deciding benchmark, without
  fewer proteins. A change that is harmless but gains less is not adopted.
- **Production is verified through the same chain.** A production change is accepted only after
  its own output (`scored_psms.tsv`, via `prod2win.py`) has gone through the same export chain
  and entrapment reading as the harness reference. It must then match that reference.

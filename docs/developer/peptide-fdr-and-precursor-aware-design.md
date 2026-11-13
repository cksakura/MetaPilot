# Peptide-level FDR and precursor-aware competition (design)

Status: diagnosis and proposal, 2026-10-23. Builds on `identification-strategy-conclusions.md`
(cascade; peptide-level FDR is the bottleneck; conflicts are mostly different precursors).
Every number below is mocked from the cached fecal entrapment searches (IPX0003851004), with no
rerun. Labels are resolved by sequence, and FDP = 2E/(T+E) unless stated otherwise.

## 1. Where the peptide-level excess comes from

The diagnosis uses the closed stage alone (production `validateClosed`, rank 1, no DIA-NN). It
reproduces the run exactly (372,255 targets and 3,722 decoys at q 1%).

| | real | entrapment | decoys | decoy estimate | FDP | 2E/D |
|---|---|---|---|---|---|---|
| PSMs, q 1% | 368,068 | 2,104 | 3,722 | 1.00% | 1.14% | 1.13 |
| Peptides, walk at 1% | 59,179 | 386 | 599 | 1.00% | 1.30% | 1.29 |

The 2E/D column is entrapment-implied false matches per decoy. It is 1 when decoys model the false
targets. The peptide walk adds a factor of 1.14 on top of the PSM level. Three causes were
measured.

### 1.1 Neighbour wins that decoys cannot copy

`EntrapmentShuffle` permutes residues within each cleavage segment. So every real peptide has an
entrapment twin with the same length, composition and precursor mass. The `rev_` decoys come from
whole-protein reversal, which moves the cleavage sites, so a decoy peptide is not a mass twin of a
target peptide. A twin can win a spectrum of a present peptide when that spectrum is weak; a
decoy cannot. These are neighbour errors, which decoys do not model (Lin et al., J Proteome Res
2022, "double competition").

| | 2E/D, PSM | 2E/D, peptide |
|---|---|---|
| E/D whose composition equals an accepted real peptide | 1.81 | 2.29 |
| all other E/D | 1.05 | 1.17 |
| rank 2 or 3 of the spectrum is an identified real peptide | 1.71 | 2.10 (majority of the unit's PSMs) |
| no identified peptide below it | 1.07 | 1.20 |

Twins are 16% of entrapment PSMs and 10% of decoy PSMs; the decoy share is chance composition
matches. Neighbour wins explain the PSM-level excess. On the real side the same mechanism exists
through homologous MAG peptides. Those are labelled real, and decoys do not model them either.

### 1.2 The peptide-count feature moves the walk toward low-count units

The closed classifier was retrained in the mock without `log_same_peptide_count`. 2E/D then
becomes 1.16 at the peptide level and 1.15 at the PSM level, so the extra peptide-level factor
disappears.

| Closed classifier | peptides at about 1.1% FDP (walk threshold) |
|---|---|
| production features | 58,568 at 1.09% (0.008) |
| without peptide count | 55,685 at 1.11% (0.009) |
| without both counts | 47,402 at 1.07% (0.009) |
| production + DIA-NN | 61,293 at 1.07% (0.008) |
| DIA-NN, without peptide count | 60,368 at 1.10% (0.009) |

The feature still ranks better at matched FDP, so dropping it is not the fix. It changes which units
pass:
- The units only it admits have about 1.5 PSMs each, with 2E/D 1.40.
- The units it drops have about 10 PSMs each ("sticky" sequences that win many weak spectra,
  decoys and targets alike), with 2E/D 1.09.

Neighbour wins sit in the low-count units, so the decoy estimate falls behind there. At matched
score, decoy and entrapment units have the same count distributions (mean 10.1 against 9.8 in the
top band). The counts themselves are symmetric; it is the low-count part of the ranking that is
not.

DIA-NN predictions help at the peptide level too when compared at matched FDP within the closed
stage (+4.7%). This sharpens conclusion 4: the loss in the cascade comes from the threshold, not
from the predictions.

### 1.3 Entrapment space is 1.117 times the real space

One searched database (Qiao-FB1-DDA-1, as searched: K/R, 2 missed cleavages, 7-50 residues,
I = L) has 14.05 M distinct real peptides, 15.69 M entrapment peptides and 29.70 M decoy
peptides:
- Decoys match targets (ratio 0.999).
- Entrapment does not: a peptide shared by several real proteins is shuffled differently in each
  copy.

With r = 1.117, the combined estimator is E(1 + 1/r)/(T + E). That lowers every FDP above by a
factor of 0.95. This is a constant scale, not a cause; the comparisons above use r = 1.

The 342 "other" target peptides are real protein N-termini after Met clipping (211 of 217 in one
run), so they do not affect the FDP.

## 2. Proposal for the peptide level

The decoys must be able to make the neighbour errors that the targets make. Theory: under TDC, the
estimate is unbiased when a false target and a decoy are exchangeable. A decoy that is never the
same mass as a present peptide is not exchangeable with a target neighbour of that peptide.

1. **Cleavage-preserving decoys.**
   - Generate each decoy peptide by reversing the residues of its target segment and keeping the
     C-terminal K/R (pseudo-reversal per cleavage segment, the same construction as the
     entrapment shuffle).
   - Every target peptide then has a decoy partner of the same mass and composition. A present
     peptide's weak spectra can be won by its decoy partner as often as by a false neighbour.
2. **Paired peptide competition** (Lin et al. 2022, "PSM-and-peptide").
   - After spectrum-level competition, each target peptide competes with its paired decoy. The
     unit with the higher best-PSM score survives, and the walk counts survivors.
   - This needs the pairing that (1) provides, and it applies per class before the union step.
3. **Keep the peptide-count and DIA-NN features.** At matched FDP they rank better; (1) and (2)
   correct the estimate rather than the ranking.

(1) changes the searched database, so it needs a re-search on HPC (closed and open). The entrapment
run must use the same decoy rule for the benchmark to stay comparable, and the shuffle is already
cleavage-preserving. (2) is a change to `ClosedOpenCascade.filterPeptides`, gated on a decoy rule
that pairs.

What to expect: the PSM-level excess is explained by neighbours (2E/D 1.05 without them). With
paired decoys, both levels should come within the PSM-level calibration. Until the re-search
exists, this is a prediction, not a measurement.

## 3. Precursor-aware competition (mocked)

An open winner explains the *same* precursor when two conditions hold: the Crystal-C charge equals
the closed charge, and the neutral masses agree within 0.02 Da after a whole number of isotope
spacings (|k| <= 3). Otherwise Crystal-C *reassigned* it. Both precursors are read per scan from
`scored_psms.tsv`.

- Of the accepted open-unmodified PSMs in the open-only arm, 137,295 of 386,708 (36%) are
  reassigned. No known-shift PSM is.
- In the cascade's open stage (grouped mode), 75,795 of the 92,870 open-unmodified PSMs kept
  (82%) are for a reassigned precursor. The open stage mostly identifies co-isolated precursors.

| Mode | PSMs (real / FDP) | peptides (real / FDP) |
|---|---|---|
| cascade as run (reproduces it exactly) | 455,599 / 0.56% | 64,340 / 1.20% |
| reassigned in its own FDR group | 455,760 / 0.56% | 64,431 / 1.24% |
| + reassigned on every scan (second ID) | 509,419 / 0.55% | 66,259 / 1.28% |

- The second identification per scan adds 53,820 real PSMs (+11.8%) at unchanged PSM error. The
  reassigned class alone has 128,841 real PSMs at FDP 0.59%.
- It adds only 1,919 peptides, at a marginal peptide FDP of about 4% (40 entrapment). Most second
  IDs are peptides already identified, and the new ones arrive through the peptide union, which
  is section 1's problem.

Proposal:
- Make reassigned precursors their own hypothesis and FDR group, and allow the second
  identification per scan.
- Evaluate peptides only after section 2. Until then the second ID is a PSM-level gain that the
  peptide step does not keep.

## 4. Scratch tools

These are scratch files, not in the repo. They live under
`C:\Users\daz25hak\AppData\Local\Temp\claude\C--Users-daz25hak-MetaPilot-MetaPilot\8ed8beb1-4c62-4f51-86c3-c41bebdffdee\scratchpad\pep\`:
- `ClosedMock2.java`: closed-stage feature ablations.
- `eval_closed.py`, `curve_closed.py`: decoy estimate, FDP, twin split and matched curves.
- `twin.py`, `parasite.py`, `counts.py`, `marginal.py`: neighbour and count diagnostics.
- `reff.py`: distinct-peptide ratios.
- `precursors.py`, `closed_ranks.py`: per-scan precursors and closed top 3 from `scored_psms.tsv`.
- `pa_cascade.py`: precursor-aware cascade mock.

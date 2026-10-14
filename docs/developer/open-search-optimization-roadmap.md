# Open-Search Workflow — Optimization Roadmap

**Basis:** measurements from raw pepXML on the 12mix dataset (sample YS-194-1-1,
marker stage) against `12mix_genomes_mapping.tsv` (12 known genomes). Every
number below is from an original engine output file, not a processed report.

**Status of prior results:** two bugs were fixed after these runs — the CRLF
params bug (open searched with *trypsin* while closed used *stricttrypsin*) and
the marker target-decoy path. Every open-vs-closed number predating the fix is
contaminated by the enzyme difference. Phase 0 exists to re-baseline.

---

## The three questions, answered from data

### Q1. Is unmodified-only genome selection reasonable?

**The purity argument is not measurable in the data; the yield cost is.**

At marker stage, 1% FDR, evidence reaching the 12 truth genomes:

| | PSMs | truth genomes found | off-target genomes | **off-target PSM share** |
|---|---|---|---|---|
| closed (±20 ppm) | 6,333 | **12/12** | 4,685 | **95.6%** |
| open, zero peak only (today) | 4,733 | **12/12** | 4,672 | **95.5%** |
| open, all hits (any shift) | 5,090 | **12/12** | 4,436 | **95.6%** |

Three findings:

1. **Recall is not the problem.** All three find 12/12. Genome selection is not
   failing to see the true organisms.
2. **The off-target share is identical (95.5–95.6%) in all three.** The
   zero-peak restriction does not measurably improve purity at the level that
   matters. Off-target evidence is dominated by *marker peptides being shared
   across the catalogue*, not by wrong-peptide matches — so restricting the mass
   window cannot fix it.
3. **It costs ~25–30% of the evidence.** Per-genome, zero-peak gets 1,196 PSMs
   where closed gets 1,649 (K. pneumoniae), 1,110 vs 1,499 (C. freundii). Across
   all 9 samples the zero peak retains only **5.5%** of marker-search PSMs.

That thinner evidence is what makes selection unstable: across 9 replicates of
the same mix, pairwise Jaccard of the selected genome sets is **median 0.35**,
only **8 of 77** genomes are selected in all 9 samples, and 22 are selected by
exactly one.

**Verdict:** the restriction is not wrong in principle, but it is paying a real
cost for a benefit that does not appear in the data. §3.1's purity claim needs
either a measurement that supports it or retirement.

**However** — open's unique peptides *are* cleaner than closed's:

| partition | peptides | **correct %** |
|---|---|---|
| overlap | 3,486 | **94.4%** |
| open only | 188 | **76.4%** |
| closed only | 676 | **57.6%** |

So the open window does avoid forced wrong matches — it just doesn't show up in
genome-level off-target share, because that metric is swamped by peptide
sharing.

### Q2. How does performance compare with closed search?

Marker stage, unmodified peptides, 1% FDR, correct = ≥1 protein in a truth genome:

| | peptides | correct | closed advantage |
|---|---|---|---|
| closed total | 4,156 | **3,639** | +246 correct |
| open total | 3,663 | 3,393 | |

**Closed finds more correct unmodified peptides.** Open's precision on its own
uniques is better (76.4% vs 57.6%) but it does not compensate for finding 493
fewer peptides overall.

This is the honest baseline: **for unmodified identification, the open workflow
is currently behind the closed one.** Its case has to rest on modifications,
which have not yet been evaluated against truth.

### Q3. Is the current validation method sound?

Three measured defects.

**(a) Decoy-based FDR is badly miscalibrated on a marker/MAG database.**

| partition | decoy-estimated FDR | true error rate |
|---|---|---|
| overlap | 1.07% | **5.6%** |
| open only | 5.20% | **23.6%** |
| closed only | 3.69% | **42.4%** |

Understated by 5–11×, and it **ranks the two unique partitions in the wrong
order** — decoys said open-only was dirtier, truth says it is 19 points cleaner.
Reversed decoys model random matches. They do not model a peptide matching a
real marker gene from a genome that is not in the sample, which is the dominant
error mode here. **Every threshold in the workflow is set against a statistic
that does not measure the error that actually occurs.**

**(b) The zero peak is fragmented.** Deamidation appears smeared across
0.981–0.991 as ~10 separate accepted peaks; near-zero shifts (−0.008…+0.008)
account for ~36% of validated PSMs as distinct "modifications". Only **5.3%** of
validated PSMs sit at a recognisable modification mass.

**(c) The `peak_id = -1` group is one hypothesis over a mixture** — 52–58% of
validated PSMs in earlier runs, passing the binomial gate on size alone
(223,560 targets vs 45,557 decoys → p = 0.0), with the local-FDR walk then
ranking on classifier score with no delta-mass reasoning.

---

## Roadmap

Ordered by dependency. Each phase has a measurable gate — we now have ground
truth, so no phase should be declared done on a plausibility argument.

### Phase 0 — Re-baseline (blocking)

The enzyme and database bugs invalidate every comparison so far.

- Rebuild; confirm `search_enzyme_nocut_1` is empty in the written params **and**
  that the pepXML reports `stricttrypsin`
- Confirm the marker search uses `<db>/Marker4_td.fasta`
- Re-run 12mix, all 9 samples

**Gate:** open and closed marker pepXMLs report identical enzyme parameters.
Re-measure the three-partition correct% table. Everything downstream is
measured against this baseline, not against the numbers above.

### Phase 1 — Make error rates measurable (highest value)

Nothing else can be tuned honestly while FDR is understated 5–11×.

- Add an **entrapment database**: append genomes known absent from the sample
  (for 12mix, any catalogue genome outside the 12) and count entrapment hits as
  empirical false positives
- Report entrapment FDR alongside decoy FDR everywhere a threshold is applied
- Recalibrate `localFdr`, `modificationGateAlpha` and `classifierTrainQValue`
  against entrapment, not decoys

**Gate:** on 12mix, entrapment-estimated FDR agrees with truth-based error rate
within a factor of 1.5 across all three partitions.

### Phase 2 — Round 1: restore genome-selection evidence and stability

Given Q1, the zero-peak restriction costs 25–30% of evidence for no measured
purity gain.

- Admit **accepted, explained** mass shifts into genome selection — a peptide
  with a confident deamidation identifies its genome as well as an unmodified
  one. Exclude only *unassigned* shifts (`peak_id = -1`), which is what the
  purity argument actually targets
- Include confidently linked isotope-family satellites after normalization
- Add a **cross-replicate reproducibility filter**: require selection in ≥k of n
  samples. Measured separation is clean — truth genomes median 7/9 samples,
  spurious genomes median 2/9, only 1/59 spurious in all 9

**Gate:** pairwise Jaccard of per-sample genome selections ≥0.70 (from 0.35);
12/12 truth genomes retained; off-target genome count reduced from 59.

### Phase 3 — Peak model: stop splitting one modification into ten

- Merge modes within instrument mass accuracy before gating — deamidation must
  be **one** peak, not ~10 across 0.981–0.991
- Set peak resolution from measured precursor accuracy, not from KDE bandwidth
- Re-check the membership floor (`max(3σ, 0.010 Da)`) after merging

**Gate:** deamidation, oxidation and the zero peak each appear exactly once in
`delta_mass_peaks.tsv`; the share of validated PSMs at a recognisable
modification mass rises from 5.3% to >40%.

### Phase 4 — Resolve the `-1` population

With Phases 1–3 done, re-measure it before choosing a policy — it should shrink
substantially once peaks are merged and round 1 is fixed.

- Split `-1` into *near a fitted peak* / *isolated* and gate each separately
- Decide, with entrapment numbers in hand, whether isolated PSMs are reportable
  identifications at all

**Gate:** each `-1` sub-population has its own entrapment-validated error rate;
no single group exceeds 10% of validated PSMs without one.

### Phase 5 — Modification claims against truth

Only now is the workflow's actual purpose testable.

- Evaluate modified peptides the way we evaluated unmodified ones: overlap /
  open-only / closed-only against the 12 genomes
- Report the correct% of modification claims, not just their count

**Gate:** open's modified peptides show a correct% within 10 points of its
unmodified overlap (94.4%). If they do not, the modification output is not
reportable regardless of yield.

### Phase 6 — Quantification, then second engine

- LFQ through the existing `FlashLfqTask` seam, keyed on
  `base sequence + localized site + peak-family ID + centroid mass`
- Sage: `combined_msfragger_sage/` is **empty** — the arm has never produced
  output. Get one run, then measure whether the second engine earns its cost

---

## What would make me abandon the approach

Worth stating explicitly, since the workflow is expensive:

If after Phases 0–3 the open workflow still finds fewer correct unmodified
peptides than closed (currently 3,393 vs 3,639) **and** Phase 5 shows
modification claims below ~85% correct, then open search is not paying for
itself on MAG data and the honest conclusion is to use the closed workflow for
identification and treat open search as a discovery-only side analysis.

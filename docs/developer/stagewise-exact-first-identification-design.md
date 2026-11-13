# Stage-wise, exact-first identification (design)

Status: design, revision 2 (2026-10-21, after review). To be built as an **experimental arm**, not
the default. It builds on the unified design (`unified-rescoring-batchwise-fdr-design.md`), whose
components it reuses, and replaces that design's single per-spectrum competition across all
explanations.

Revision 2, from the review:
- **(R1)** No recurrence in the PSM score by default; it is an ablation arm only (4.2).
- **(R2)** The FDR statements are qualified (section 3). Stage-1 blocking, RESET roles across stages
  and the validation of the union are now specified (4.3).
- **(R3)** RESET roles are assigned per decoy *peptide* and kept at the peptide level; the peptide
  model trains only on training-decoy peptides (4.3).
- **(R4)** Replication supports the evaluation but does not decide adoption. Entrapment is a
  prerequisite for promotion, the exported sets are evaluated, and results carry uncertainty
  intervals (section 7).
- **(R5)** Equivalent hypotheses are merged across searches, distinct peptidoforms kept (4.1). The
  shifted-rival context features are explicit (4.2).

Revision 2 also reverses an earlier decision (2026-10-21): the recurrence counts had been put into
the PSM classifier, and they are now off there by default (R1).

## 1. Question

How should the closed and the open search be combined so that:
- the error is controlled by target-decoy competition (TDC) on exactly what is reported;
- every spectrum gets the best candidate both searches offer;
- an exact-mass answer is never lost to a shifted one that is wrong?

## 2. Evidence

**12mix** (cascade run, 2026-10-20; unified run 3). Closed-accepted target spectra whose winner
became another peptide, by where that winner landed:

| Winner landed on | Spectra | Winner is a decoy | Closed answer outside spiked | Winner outside spiked |
|---|---|---|---|---|
| exact mass, other peptide (closed rank 2-3 or open zero family) | 89,134 | 0.6% | 4.4% | **2.5%** |
| unexplained delta | 43,408 | **27.3%** | 5.0% | 17.0% |
| retained delta-mass peak | 9,500 | **23.2%** | 5.8% | 19.1% |

- Shifted overrides of an exact answer are mostly wrong: random matches at random deltas,
  relatives' peptides carrying a substitution as a shift, and nested sequences.
- Same-mass alternatives are better than closed rank 1.
- On the closed-accepted target spectra, the open winner is a different peptide in 37% of cases
  (213,149 of 580,247).

**Fecal** (IPX0003851004, 9 runs = 3 samples x 3; cascade and joint arms of one run, 2026-10-21).
Replication = the method accepts the peptide in >= 2 runs.

| | Spectra | Answer replicates |
|---|---|---|
| both methods, same peptide | 402,819 | 98.0% |
| both methods, different peptide | 48,589 | cascade 94.7%, **joint 97.4%** |
| only cascade | 25,819 | 89.0% |
| only joint | 20,780 | 89.1% |

- Cascade accepts more at equal quality: +5,000 spectra and +4,400 peptides. Its extra spectra are
  mostly closed-stage answers (22,393).
- Joint chooses better when the searches disagree on a spectrum.
- The cascade wins on what to accept, the joint model on which candidate to pick.

## 3. Theory, and what it does and does not establish

TDC estimates the false discoveries among a reported set under conditions the procedure has to
respect. They are necessary, not sufficient:

1. **(T1) Label-blind winner.** Each spectrum has one winner, chosen without its label, by a model
   that has not seen that spectrum's label (cross-validation by spectrum, RESET for training,
   Freestone et al., JPR 2025).
   - **Not established by this:** that an incorrect target winner and a decoy winner are
     exchangeable. Label-blind scoring is necessary for that, not sufficient. A feature that a
     decoy formally "can have" can still distribute differently over wrong targets and decoys.
     Recurrence is the measured example: a wrong match to a real, abundant peptide inherits a high
     count; a decoy rarely does (R1).
2. **(T2) The estimate refers to the reported PSMs.** The decoys among the accepted winners
   estimate the false targets among the same winners. Changing a reported PSM's peptide after
   acceptance breaks this. "Accept with the cascade, then swap in the joint model's peptide" is
   therefore excluded: the choice happens inside a stage, before its TDC.
3. **(T3) Sequential stages.** The cascaded-search argument (Kertesz-Farkas, Keich and Noble, JPR
   2015, https://pmc.ncbi.nlm.nih.gov/articles/PMC4533645/) supports running a later stage only on
   spectra an earlier stage did not accept, with its own TDC.
   - **Not established by this:** that the **union** of stages each at expected FDR alpha has
     expected FDR at most alpha. The weighted-average argument bounds a union's realised FDP by the
     largest stage FDP. It does not carry over to expected FDR when the discovery counts are random.
   - **Not established either:** that the cascade result holds for *this* combination of learned
     candidate selection, RESET and Group-walk. The paper does not cover it.
   - The union is therefore reported with its own decoy estimate and validated by entrapment (7.1).
     No guarantee is claimed.
4. **(T4) Groups within a stage** are combined by Group-walk (Freestone et al., Bioinformatics 2022,
   https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii82/6701992). Its guarantee
   rests on its own assumptions (the TDC exchangeability above, within each group), which are no
   stronger here than in T1.

Consequences:
- **Stage membership is decided from mass alone, before any label is read.** Exact vs shifted is
  such a property.
- **Exact before shifted is a prior, not a theorem.** It is justified by section 2 and by decoys
  being unable to model target-target errors (a relative's peptide). The alternative is CONGA's
  separate narrow group, applied here sequentially.
- **What the procedure can claim:** decoy-estimated FDR per stage and for the union, under T1-T4,
  checked by entrapment. Not rigorous FDR control of the complete procedure.

## 4. Design

```
spectra, with candidates from the closed search (top 3) and the open search (top 3)
  |
  Stage 1 - exact-mass candidates: closed top 3 + open zero family (no shift)
  |   winner per spectrum by the stage-1 model (both searches' evidence)
  |   TDC at alpha over the stage-1 winners  ->  accepted A1
  |
  Stage 2 - spectra not in A1; candidates on retained delta-mass peaks
  |   winner per spectrum by the stage-2 model
  |   Group-walk over the peak groups at alpha  ->  accepted A2
  |
  (Stage 3 - unexplained deltas: reported separately, not identifications by default)
  |
  Peptides: each peptide in the first stage where it is accepted, walked within that stage
```

### 4.1 Candidates and stages (mass only)
- **Exact (stage 1):**
  - the closed search's candidates, with their variable modifications and isotope offset 0-2;
  - open candidates whose delta is in the zero family (isotope errors folded in).
- **Merging (R5).** A hypothesis is a **peptidoform**: sequence plus mass explanation.
  - Each closed candidate is its own hypothesis, keyed by spectrum, sequence and its variable
    modifications. The same backbone with different modification placements stays distinct.
  - An open candidate is merged into a closed one, and flagged "found by both", only when they
    explain the same mass: same spectrum and sequence, and open delta = the closed candidate's
    modification mass + isotope offset x 1.00336, within 0.02 Da.
  - Otherwise the open candidate stays a hypothesis of its own, for example the same backbone with
    a different shift.
- **Shifted (stage 2):** open candidates whose delta is on a retained peak. Peaks are detected and
  filtered label-blind as in the unified design (F1-F6), and mapped to UniMod only afterwards.
- **Unexplained:** open candidates on no retained peak. Stage 3 or not reported.
- Neighbour filter within a stage only: more than 50% shared fragment ions, the lower-ranked
  removed. Ranking uses the label-free annotation score, with no credit for a shift.

### 4.2 Stage models
- **One model per stage,** trained only on that stage's candidates. The stage-2 model is trained
  only on the spectra stage 1 left, so it learns that population.
- **Listwise:** one query per spectrum, holding all of the stage's candidates plus a "none" option
  (the unified listwise trees).
- **The winner** is the highest-utility candidate. **The score** is its log-odds against the stage's
  other candidates and "none".
- **Training:** semi-supervised, cross-validated by spectrum, RESET decoy split (T1).

**Features**, all available to decoys too:
- MetaPilot's own fragment annotation;
- the DIA-NN backbone prediction (MS2 similarity, retention-time error);
- both-search agreement (found by closed, by open, by both).
- **Shifted-rival context (stage 1, explicit, R5):** what the open search offers against the exact
  candidate on the same spectrum. Three features:
  - the annotation-score margin over the best shifted candidate;
  - that rival's peak quality (enrichment, fragment-support rate);
  - whether the rival was found by both searches.

  A strong, well-explained shifted rival lowers an exact winner's score; it cannot replace it.
- **Stage 2 adds:** peak quality (enrichment, fragment-support rate, width), delta-placement
  features, and the spectrum's stage-1 outcome (whether it had an exact candidate, and the stage-1
  winner's score).
- **Recurrence (R1): not in the PSM score by default.** On 12mix, with the counts in the PSM score,
  the RT consistency check estimated 6.2% of accepted PSMs RT-inconsistent, against 3.2% without, at
  the same decoy estimate. A wrong assignment to a real, abundant peptide inherits a strong count; a
  decoy usually cannot.
  - The counts (PSMs with the same peptide, PSMs and peptides from the same protein, run share) stay
    available as an **ablation arm** (`--psm-recurrence`).
  - That arm is enabled only after independent validation (entrapment).
  - Recurrence is used where it is justified: the **peptide-presence** decision (peptide model, 4.3).
- **Not used:** search-engine scores (expect, hyperscore, ranks), as in the unified design.

### 4.3 FDR
**RESET roles (R3).** Each decoy *peptide* (base sequence, I = L) is assigned once, by a hash of its
sequence, to *training* or *estimating*. Every PSM of that peptide, in every stage, shares the role.
- **Training decoys** train the stage models and the peptide model and never enter a walk.
- **Estimating decoys** never train anything and are the only decoys in any walk.
- Per-PSM roles (the earlier implementation) split one decoy peptide's PSMs between the two roles,
  and then trained the peptide model on decoy units that were later also estimating. That is not
  RESET's separation.

**PSM level:**
- stage 1 is one TDC over its winners;
- stage 2 is Group-walk over the peak groups, merged by KS test as in CONGA (minimum 2K), with
  pseudo-random tie order;
- estimate (2 x estimating decoys + 1) / targets.

**Stage-1 blocking (what stage 2 sees).**
- tau1 is the lowest stage-1 winner score with PSM q <= alpha.
- Every spectrum whose stage-1 winner scores >= tau1 is removed from stage 2: targets, estimating
  decoys and training decoys alike. The rule depends on the score only, so it reads no label beyond
  what defined tau1.
- Blocking is decided at the PSM level only. A spectrum accepted in stage 1 stays blocked even if
  its peptide later fails the peptide-level filter. It is then not reported, and stage 2 does not
  get it.
- Letting the peptide-level decision unblock spectra would couple the stages through a second,
  label-dependent decision.

**Peptide level:**
- Units are base sequence (I = L) plus label, per stage, built from all of the stage's winners.
- A stage's peptide model is trained on that stage's units:
  - positives are pseudo-target units (targets and estimating decoys) at q <= the training q;
  - negatives are training-decoy units;
  - folds are by peptide.
- Its walk compares target units with estimating-decoy units only (weight 2).
- A peptide accepted in stage 1 is not walked again in stage 2: units with the same key that passed
  stage 1 are excluded from stage 2, for targets and decoys alike.
- Identifications are the PSMs accepted in their stage whose peptide passed in stage 1 or stage 2.

**Union.** Reported with its own decoy estimate (sum over the stages), and validated by entrapment.
No guarantee is claimed for it (T3).

### 4.4 What it changes relative to cascade and joint

| | Cascade | Joint | This design |
|---|---|---|---|
| Who competes in a spectrum | closed rank 1 only, then open on the rest | every candidate of both searches | exact candidates of both searches; shifted ones only where no exact answer was accepted |
| Who chooses the winner | each search's engine rank | one model over everything | one model per stage, seeing both searches |
| Can a shifted hit replace an exact answer | no | yes | no |
| Can a same-mass alternative replace closed rank 1 | no | yes | **yes** |

## 5. Reused components
From `core/unified`, all already tested:
- `SpectrumAnnotator`;
- `DiannBackbone`;
- `DeltaMassBatches` (detection, filters, mapping);
- `ListwiseGbdt`;
- `GroupWalk` (with KS merging and tie order);
- the recurrence counts and the peptide-level model.

New code:
- the stage split;
- per-stage training and walks;
- the stage-2 "rejected exact candidate" feature;
- the stage-wise peptide assignment.

## 6. Risks
- **A true modification whose unmodified form also matches.** A modified peptide that the closed
  search also matches unmodified at exact mass (for example a labile modification lost in source)
  is decided in stage 1. That is correct when it matches the spectrum at exact mass.
- **A wrong exact answer accepted in stage 1** blocks a right shifted answer for that spectrum. The
  stage-1 model sees the shifted-rival context (4.2), so a strong shifted rival lowers the exact
  winner's score without replacing it.
- **Exchangeability (T1)** may fail for features beyond recurrence. The per-feature target/decoy
  AUC report and the entrapment run are where that shows.
- **Small stage-2 groups:** handled by KS merging and the minimum group size.

## 7. Evaluation

### 7.1 Required before promotion: entrapment (R4)
- **Setup:** paired shuffled proteins of the selected genomes, added 1:1 to each refined database
  (`EntrapmentShuffle`, `MGYG8...`) before the closed and open searches.
- **Measured on the exported PSM and peptide sets** (what a user receives), per stage and for the
  union, with the combined and paired estimators (Wen et al., Nature Methods 2025).
- **Reported with uncertainty intervals** (bootstrap over spectra, and over peptides). Pass when
  the interval is consistent with alpha, not when every point estimate is below alpha.
- The genome-level entrapment in `human-gut-shadow-s1` (`MGYG9...`) tests genome selection only.
  No shadow genome was selected there, so no refined database held entrapment proteins.

### 7.2 Supporting evidence, not deciding: replication (R4)
- Fecal (IPX0003851004): against cascade and joint from the same run, overall and on the
  disagreement sets of section 2, with bootstrap intervals.
- **Limits:**
  - replication rewards recurrence, which the peptide model also uses;
  - a systematic wrong assignment can repeat across runs;
  - peptide recurrence does not verify a particular spectrum's assignment.

  It indicates; it does not decide.
- **Held-out runs.** Development uses some runs and datasets; separate runs and datasets are kept
  back for the final validation.

### 7.3 Sanity only
12mix: a defined mock community, so no tuning on it.

### 7.4 Ablations
- no recurrence in the PSM score (default) vs `--psm-recurrence`;
- with and without the DIA-NN block;
- stage-wise vs single competition (the unified arm).

All are judged by 7.1 when it is available.

## 8. Implementation status (2026-10-21)

Implemented in `core/unified/UnifiedIdentification` (`Config.stagewise`, `Config.psmRecurrence`;
CLI `--stagewise`, `--psm-recurrence`). Included:
- peptidoform merging (4.1, R5), tested in `UnifiedIdentificationTest`;
- explicit shifted-rival features (4.2);
- recurrence off in the PSM score by default (R1);
- a RESET role per decoy peptide;
- stage-1 blocking by tau1 whatever the label;
- per-stage peptide walks, with stage-2 units that passed stage 1 not walked again;
- per-stage and union summary lines;
- the arm `stagewise` in `unified_arms.tsv.gz`.

Sanity run, 12mix, closed top-1, no DIA-NN, no tuning:

| | PSMs kept | Peptides |
|---|---|---|
| Stage-wise | 644,682 | 76,628 |
| Unified, same inputs | 612,260 | 73,632 |
| Closed alone | - | 88,236 |
| Cascade | 746,278 | 90,963 |

- Stage 2 adds 1,706 PSMs and 816 peptides. Its peptide model barely separates (AUC 0.53).
- Stage 1 alone keeps 79.7% of the closed-alone peptides. The gap to closed alone is the
  exact-stage model, not the stage structure: the machinery check gave 74,839 peptides with these
  features against 94,046 with the closed engine's features.
- Union estimates: PSM 0.28%, peptide 1.00%.

12mix results are sanity only (7.3); the decision waits for 7.1.

## Sources
- Cascaded search: Kertesz-Farkas, Keich, Noble, J Proteome Res 2015 (cited in `ClosedOpenCascade`).
- Group-walk: https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii82/6701992
- CONGA: https://www.biorxiv.org/content/10.1101/2023.05.02.539167v3.full
- Percolator-RESET: https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00742
- Entrapment / FDRBench: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12240826/

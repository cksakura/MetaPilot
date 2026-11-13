# Fragment-level features against same-mass neighbours (design)

Status: design, not implemented. It follows `identification-principles.md` (P1–P9) and builds on
identification v2 (`identification-v2-design.md`).

## 1. The problem

Under identification v2 the remaining ranking errors are dominated by **same-mass neighbours**. A
neighbour has the same residue composition, and so the same precursor mass, as a peptide that is
present. It wins a weak spectrum of that peptide.

Three mocks changed the learner and left the features alone. None moved yield at matched FDP by
more than 2%:

| Change | Effect on yield at matched FDP |
|---|---|
| Logistic regression → monotone trees | ±2% |
| `scan_has_reassigned` feature | under 1% |
| Within-query ranking objective | under 1% |

The learner is no longer the limit; the information is.

### 1.1 How large the neighbour share is (harness Btm, pseudo-reversed runs)

| Run | Level | Entrapment errors that are composition twins of an accepted real peptide |
|---|---|---|
| fecal | PSM | 881 of 2,373 (37%) |
| fecal | peptide | 374 of 1,069 (35%) |
| 12mix | PSM | 1,499 of 3,867 (39%) |
| 12mix | peptide | 842 of 2,455 (34%) |

Source: `scripts/identification_audit/results/twin_split.txt`.

### 1.2 Is the displaced peptide in the same query?

The table covers exported Btm winners (q ≤ 1%, identifying classes). It shows how often the query
also holds another candidate that is real, or one that is the winner's real composition twin.

| Run | Winner | Other candidates per query (median) | A real candidate in the query | The winner's real twin in the query |
|---|---|---|---|---|
| fecal | entrapment | 4 | 70% | 25% |
| fecal | decoy | 4 | 66% | 15% |
| fecal | real | 4 | 77% | 2.7% |
| 12mix | entrapment | 2 | 63% | 15% |
| 12mix | decoy | 2 | 61% | 13% |
| 12mix | real | 3 | 73% | 1.1% |

Two consequences shape the design:

1. **A rival comparison reaches only part of the errors.** Comparing the winner with a rival in its
   query reaches the 15–25% of neighbour wins where the displaced peptide is a candidate. The rest
   need evidence about the winner's sequence on its own: does the spectrum confirm this residue
   order, or only this composition?
2. **Decoys already make the same kind of error.** Pseudo-reversed decoys win against their real
   twin at a rate close to the entrapment's (13–15% against 15–25%), so they are neighbour errors
   in the same sense. A feature that separates a neighbour from its source should therefore also
   separate a decoy from its source. This is what lets the model learn these features from decoys
   (P4) without the entrapment.

## 2. Theory

### 2.1 Why the current scores cannot separate neighbours

Two peptides with the same composition share the precursor mass. They also share every fragment
that covers the same multiset of residues. That always includes the C-terminal K/R (y1) and the
full-length fragments, and by chance other ions too.

The current scores are dominated by these shared ions:

| Score | Why it does not separate neighbours |
|---|---|
| Hyperscore and expect | Count matched ions without asking which ions discriminate. |
| `ms2_entropy`, `ms2_matched_fraction` | Weight the 24 most intense predicted fragments. These are mostly mid-to-high y ions, which a neighbour shares more often than not, because the y-series of a shuffle keeps the C-terminal block. |
| `hyperscore_margin` | Compares totals, not the ions that differ. |

The information that separates residue orders lies in two places:

- **Distinguishing ions:** fragments whose m/z differs between the two sequences.
- **Ladder continuity:** consecutive backbone cleavages that are each confirmed. A true sequence
  produces runs of adjacent b or y ions. A neighbour matches scattered ions by chance.

### 2.2 Distinguishing ions and their significance

This is the reasoning behind site localisation scores (Ascore, Beausoleil et al. 2006; PhosphoRS,
Taus et al. 2011), applied to residue order rather than a modification site.

For a candidate *i* and a rival *j*:

- Let *D_i* be the theoretical b/y ions of *i* (charge 1 up to min(2, precursor charge)) that are
  not within tolerance of any theoretical ion of *j*. *D_j* is the same with *i* and *j* swapped.
  Ions in both sets carry no information about which sequence is right.
- Let *p* be the probability that one theoretical ion matches a peak by chance. It is the local
  peak density times the tolerance window, counted in a ±50 Th window around each ion and averaged
  over the ions.
- If *k_i* of the *n_i* ions in *D_i* match a peak, then under "*i* is wrong" the count is
  roughly Binomial(*n_i*, *p*). The evidence is
  *S_i* = −log10 P(X ≥ *k_i* | *n_i*, *p*).
- *S_i* − *S_j* is a log-likelihood-type contrast between the two orders, comparable across
  peptide lengths and spectrum densities. That is the property hyperscore lacks.

### 2.3 Evidence without a rival

When the displaced peptide is not a candidate (75–85% of neighbour wins), the winner is judged
against the space of its own rearrangements:

- **Ladder coverage:** a backbone site *s* (between residues *s* and *s* + 1) is confirmed when
  b_s or y_(n−s) matches. A rearrangement of the same composition keeps few confirmed sites in
  their place. The two measures are the fraction of the *n* − 1 sites confirmed, and the longest
  run of consecutive confirmed sites.
- **Whole-sequence ion significance:** *S* over all b/y ions of the candidate, with the same
  density-calibrated null as in 2.2. Unlike hyperscore, this does not reward long peptides in
  dense spectra for chance matches.
- **Unexplained intensity:** the share of intensity in the 20 most intense peaks (excluding ±2 Th
  around the precursor m/z and its isotopes) that no b/y ion of the candidate explains. A
  neighbour that won a present peptide's spectrum leaves that peptide's specific ions unexplained.
  - Peaks explained by the best candidate of the scan's other query count as explained. This keeps
    a co-isolated second precursor from being penalised.
  - This uses the spectrum, not the competition: the queries still compete separately and are
    still FDR groups of their own (P3).

### 2.4 Why calibration can hold

TDC with RESET stays valid when, for incorrect identifications, the features do not separate
estimating decoys from false targets (P4).

- **Built in:** every feature here is computed identically for target, decoy and entrapment
  candidates. The rival is chosen by a label-blind fixed search score, not the model score, which
  also avoids feedback between iterations.
- **Where decoys and entrapment differ:** both keep the first and last residue of every K/R segment
  (`TargetDecoyFasta.pseudoReverse`, `EntrapmentShuffle`). So both share b1 and y1 with their
  source, and the terminal ions are no risk. The interiors differ:
  - The pseudo-reversal reverses the interior, deterministically.
  - The entrapment shuffles it at random.
  - A reversed interior keeps every adjacent residue pair, in mirror order. A shuffle breaks most
    of them. Distinguishing-ion counts and ladder runs against the source could therefore differ
    systematically between a decoy and an entrapment neighbour.
- **Mitigation:** b1, b2, y1 and y2 are left out of every set, which removes the ions that every
  rearrangement shares. Two checks test the result:
  - the exchangeability check (section 4.3, step 2) measures whether the interior difference shows
    up in the features;
  - the calibration gate (P7).

## 3. Features

Per candidate, in its own query (P3). The rival *j* is the best other candidate of the query by
max(closed, open) hyperscore. For the query's best candidate it is the second best. A query with
one candidate has no rival, and those features are missing; the trees learn a default direction
for missing values.

| Name | Definition | Monotone |
|---|---|---|
| `nf_dist_n` | \|*D_i*\|, the number of distinguishing ions against the rival | 0 |
| `nf_dist_frac` | matched share of *D_i* | +1 |
| `nf_rival_dist_frac` | matched share of *D_j* (the rival's own distinguishing ions) | −1 |
| `nf_dist_contrast` | *S_i* − *S_j* (section 2.2) | +1 |
| `nf_dist_pred` | DIA-NN predicted intensity on *D_i* that has a matched peak, over all predicted intensity on *D_i* (missing beyond open rank 3) | +1 |
| `nf_ladder_cov` | confirmed backbone sites / (*n* − 1) | +1 |
| `nf_ladder_run` | longest run of confirmed sites / (*n* − 1) | +1 |
| `nf_ion_sig` | *S* over all b/y ions (b3+, y3+), density-calibrated | +1 |
| `nf_unexplained` | unexplained share of the top-20 intensity (section 2.3) | −1 |

Other specifications:

- **Ions:** b and y, charge 1 up to min(2, z), positions 3 to *n* − 1. No neutral losses in the
  first version.
- **Tolerance:** `joint_fragment_tolerance_ppm` (20 ppm).
- **Modifications:** closed candidates carry their variable modifications. For open candidates,
  the delta sits at the site `PredictedSpectrumFeatures.score` chose (the localised site, or the
  site with the best entropy). That method will return the chosen site so both feature sets use
  the same placement.
- **Monotone constraints:** XGBoost supports −1 directly. The v2 model's constraint string is
  extended with the signs above.

## 4. Evaluation (mock first, from cached outputs)

### 4.1 What the mock needs

- A Java tool, `scripts/identification_audit/FragFeat.java`, reads `cands.parquet` (exported to
  TSV), the run's mzML files (`MzmlScanReader`) and, for `nf_dist_pred`, DIA-NN predictions. It
  writes `frag.tsv.gz` with one row per candidate.
  - If the cached `joint_predictions/` cannot be reread, the mock runs without `nf_dist_pred`.
- `harness.py --frag` joins `frag.tsv.gz` to the feature table.

### 4.2 Arms

All arms are Btm settings (RESET, precursor queries, monotone trees, 10 iterations), on both
pseudo-reversed runs:

| Arm | Features |
|---|---|
| Btm | none of these (reference) |
| BtmF | all features in section 3 |
| BtmFa | rival-free only: `nf_ladder_*`, `nf_ion_sig`, `nf_unexplained` |
| BtmFr | rival only: `nf_dist_*` |

### 4.3 Reading and gates (P6, P7, P9)

1. **Calibration gate.** At every level of the export chain, the entrapment FDP must lie within
   the decoy estimate's 95% interval or below it, on fecal. 12mix must not get worse than Btm.
2. **Exchangeability check.** Compare each feature's distribution between exported decoy winners
   and exported entrapment winners (KS statistic). A feature where they clearly differ is a
   calibration risk, even if the gate passes.
3. **Mechanism check.** The twin share of entrapment errors (`twin_cmp.py`) must fall. So must the
   fraction of entrapment wins whose real twin is in the query (table 1.2).
4. **Yield.** Real identifications at matched entrapment FDP 1%, at the PSM, peptide and protein
   levels.
   - Adopt when fecal peptides gain at least 1% and proteins do not fall.
   - 12mix is a sanity check only: no level may drop by more than 1%, and all 12 mix genomes must
     be accepted.

### 4.4 What to expect

The ceiling is the neighbour share: about 35–39% of entrapment errors. Removing half of them at a
fixed threshold frees about a fifth of the error budget at matched FDP. Spent on the borderline
PSMs, that is roughly 1–3% more PSMs and somewhat less at peptide level. A gain under 1% means the
neighbours are not separable with b/y evidence at this spectrum quality. The next step would then
be better spectrum prediction (Prosit or AlphaPeptDeep), not more features.

## 5. Production, only after the mock passes

- A new class `opensearch/NeighbourFragmentFeatures`, computed in the same mzML pass as
  `PredictedSpectrumFeatures`, with its own unit tests:
  - distinguishing sets for a known pair;
  - binomial significance;
  - the ladder on a synthetic spectrum;
  - label blindness.
- `JointFeatures` gains the nine features under identification v2. `V2_MONOTONE_FEATURES` becomes
  a map from feature to sign.
- The change in cost is one more pass over each candidate's b/y ions. The peaks are already
  loaded, so no new I/O.

## 6. Order

1. `FragFeat.java` and `harness.py --frag`.
2. Run the four arms on both runs, with every level reported.
3. Apply the gates (section 4.3).
4. Only then, the production class and tests.

# Explaining mass shifts by localisation (design)

Status: the BtmX design (section 7) is implemented in production under identification v2, as
`ShiftLocaliser` and the `OPEN_LOCALISED` group of `JointIdentification`, with
`ShiftLocaliserTest`. It has not yet been run on HPC. It follows `identification-principles.md`
(P1–P13) and builds on identification v2.

## 0. Scope (owner decision, 2026-10-25)

The objective is not to explain every shift.

- **Must do:** MS2 analysis and localisation. For each PSM, test whether the fragments carry the
  shift, and where.
- **Assign only what fragments support.** A site and an explanation are assigned only when the
  fragment evidence supports them.
- **Otherwise the shift stays ambiguous.** It is recorded with its mass, and no site or
  explanation is claimed.
- **Out of scope:**
  - precursor-level causes (monoisotopic error, co-isolation, adducts, labile losses);
  - two-site splits;
  - resolving isobaric same-residue pairs beyond reporting the existing secondary evidence.

## 1. Where we start

### 1.1 The pool

Identification v2 exports exact, unmodified and allow-list-shift identifications. Every other
open-search winner is in the **unexplained-shift** class and is never exported. On the fecal
entrapment run (production v2, pseudo-reversed decoys):

| Measure | Value |
|---|---|
| Unexplained-shift PSMs passing their own 1% walk (decoy weight 2) | 46,544 |
| Entrapment FDP of those PSMs | 1.30% |
| Distinct peptides | 17,873 |
| Peptides not otherwise identified | 8,495 (+13% on the 65k exported) |
| 12mix, the same population | 160,724 PSMs, 43,175 peptides |

The peptides are right: the class sits at about 1% FDP. What is missing is an explanation of the
shift.

### 1.2 What mass alone explains (census)

The census is `scripts/identification_audit/census.py`. It used the full Unimod table, a 5 mDa
tolerance, and a chance control: the same tests applied to each shift moved by ±0.02 and ±0.03 Da.

| Explanation | fecal share | fecal chance | 12mix share | 12mix chance |
|---|---|---|---|---|
| Extension or truncation (flanking residues) | 2.3% | 0.01% | 0.2% | 0.01% |
| Substitution, substituted peptide in the database | 2.7% | 0.07% | 0.5% | 0.02% |
| Substitution, substituted peptide not in the database | 19.2% | 5.1% | 9.6% | 3.1% |
| Unimod entry with its site | 5.7% | 4.5% | 8.4% | 4.1% |
| Plus a ±1 isotope error | 8.8% | 16.0% | 8.6% | 10.5% |
| Pair of explanations | 49.0% | 60.4% | 55.4% | 63.4% |
| None | 12.2% | 13.9% | 17.3% | 18.9% |

- **Credible explanations:** the excess over chance is about 20% of the pool on fecal and 12% on
  12mix.
- **Pairs and isotope combinations:** they fit random masses as often as real ones, so they carry no
  information.
- **Substitution against modification:** 98–99% of the substitution matches also fit some Unimod
  entry, so mass alone cannot decide between them.

### 1.3 What fragments add (fecal, MSFragger's `localize_delta_mass`)

| Explanation | PSMs | Shift supported by fragments (score gain > 0) | MSFragger site fits the explanation | Random site would fit | Single site |
|---|---|---|---|---|---|
| Extension | 1,077 | 95% | 58% | 5% | 74% |
| Substitution, in database | 1,275 | 90% | 96% | 42% | 79% |
| Substitution, not in database | 8,944 | 83% | 86% | 31% | 65% |
| Unimod | 2,646 | 40% | 43% | 17% | 57% |
| Plus isotope | 4,096 | 21% | — | — | 51% |
| Pair | 22,813 | 19% | — | — | 44% |
| None | 5,693 | 5% | — | — | 27% |
| **All** | 46,544 | **35%** | | | 43% |

By the size of the shift, the share with fragment support is:

| \|d\| (Da) | ≤ 2 | 2–5 | 5–20 | 20–60 | 60–200 | 200–600 |
|---|---|---|---|---|---|---|
| Fragment support | 0.6% | 2.7% | 57% | 65% | 40% | 20% |

Two facts follow and drive this design:

1. **Two-thirds of the pool carries no shift on its backbone.** The peptide is right, but its
   fragments match the unmodified sequence: placing the shift on any residue does not improve the
   match. The discrepancy is at the precursor. Possible causes are a wrong monoisotopic peak, a
   co-isolated precursor, an adduct, or a labile group lost before fragmentation. This is the
   largest part of the pool, and localisation cannot explain it, by definition.
2. **Where the shift is on the backbone, localisation confirms explanations independently of
   mass.** The site falls on a compatible residue 2–11 times more often than at random.
   Localisation is therefore the test that turns a mass match into an explanation.

### 1.4 What exists already

| Component | What it does | Status |
|---|---|---|
| `DeltaMassPeakFitter`, `DeltaMassKdePeakFinder` | Fit delta-mass peaks | In the workflow |
| `ExpectedPeakSelector` | Allow-list of 8 expected shifts, accepted per peak | In the workflow; the only shifts exported |
| `DeltaMassEvidenceGate` | Peaks whose members lack fragment support (MSFragger score gain) are folded onto the unmodified peptide | In the workflow; **peaks only**, not PSMs outside peaks |
| `SubstitutionReassigner` | Withdraws substitution-shaped peaks; proposes substituted sequences | In the workflow |
| `SequenceVariantCaller` | Variant calls: recurrence, codon reachability, reference observed | Report only |
| `VariantFormRanker`, `VariantFormRanking` | Substitution against modification forms at one site, by DIA-NN predicted intensities | Report only; substitutions only |
| `PangenomeVariantCheck` | Substitution sequences looked up in UHGP-100 against matched random substitutions | Report only |
| `IsobaricModificationTable` | Substitutions with the mass of a modification on the same residue (undecidable by fragments) | Shared table |
| MSFragger localisation | Best site(s) as lowercase residues, a localisation delta score, a score gain | Read into `OpenSearchPsm`; used as features |

**Missing:**
- a site-level localisation score with a probability and an error rate;
- a general explanation step over every unexplained PSM, covering modifications, terminal
  adducts and extensions, not only substitutions;
- explanation-level error control;
- any route by which an explained shift becomes an exported identification.

## 2. Theory

### 2.1 Where a shift can be: two kinds of hypothesis per PSM

For a PSM with peptide *P* (length *n*), charge *z* and shift *d*:

- **H0, not on the backbone:** the fragments are those of the unmodified *P*.
- **Hk, on the backbone at site *k*** (*k* = 1..*n*, or a terminus): every b ion with index ≥ *k*
  and every y ion containing *k* is shifted by *d*. All other ions are unshifted.

Two-site splits are out of scope (section 0). Section 1.2 shows that mass alone cannot support
pairs.

### 2.2 Fragment likelihood with a density null

The test is the same as in `neighbour-fragment-features-design.md`, section 2.2. For a hypothesis
*H*, let *I(H)* be its theoretical b/y ions (charge 1 up to min(2, *z*)) and *p* the chance that
one ion matches a peak, from the local peak density and the tolerance window.

- For two hypotheses *H* and *H′*, the evidence lies only in the ions that differ, *D(H, H′)*.
- If *k* of the *m* ions in *D* match, then under "*H* is wrong" the count is roughly
  Binomial(*m*, *p*). The evidence is *S* = −log10 P(X ≥ *k*).
- The comparison between *H* and *H′* is *S(H)* − *S(H′)*. It does not depend on length or
  spectrum density, unlike the MSFragger hyperscore gain.

### 2.3 Fragment support (H0 against the best Hk)

*S(Hk\*)* − *S(H0)* asks: do the shifted ions of the best site explain the spectrum better than
the unshifted ions they replace? It replaces the MSFragger score gain as the decision between
"not on the backbone" and "on the backbone". It is computed per PSM, not per peak. That generalises
`DeltaMassEvidenceGate` to the two-thirds of the pool outside any accepted peak.

### 2.4 Site probability and false localisation rate

Within "on the backbone", the site posterior is *P(k)* ∝ 10^*S(Hk)*, normalised over the sites
allowed by the shift's candidate explanations plus every residue. Including every residue matters:
it lets a shift land somewhere no explanation allows.

Sites that no candidate explanation of *d* allows act as **decoy sites**. This is the decoy-residue
idea used for phosphosite false localisation rates (FLR), generalised to any shift:

- Every explanation set *E(d)* has a set of allowed sites *A(d)*. Under random localisation, the
  chance of landing in *A(d)* is |*A(d)*| / *n*.
- The rate at which confident localisations fall outside *A(d)* estimates the false localisation
  rate. A localisation inside *A(d)* is false with a probability scaled by |*A(d)*| / (*n* −
  |*A(d)*|), as in the decoy-amino-acid method.
- This gives a site-level FLR per explanation class, with no external truth.

### 2.5 From a site to an explanation

At the localised site *k\**, the candidate explanations are:

- substitutions X→Y with X = *P*[*k\**] and the right mass;
- Unimod entries specific to that residue or terminus;
- an extension or truncation, only at a terminus and only where the protein's flanking residues
  match.

Three cases:

- **Separable** (different residues or termini): decided by localisation.
- **Isobaric on the same residue** (`IsobaricModificationTable`; for example S→T against methyl-S):
  backbone fragments cannot decide. Secondary evidence is reported separately, never folded into
  one score:
  - fragment intensities from predicted spectra (`VariantFormRanker`, generalised to modification
    forms);
  - the catalogue: is the substituted sequence in UHGP-100? (`PangenomeVariantCheck`);
  - recurrence across samples (`SequenceVariantCaller`).
- **No candidate at *k\****: a *localised but chemically unexplained* shift. The mass and site are
  still reported, so it is evidence for an unknown modification at a known residue.

### 2.6 What is assigned, and what stays ambiguous

Each PSM ends in one of three states:

| State | Condition | Assigned |
|---|---|---|
| **Localised, explained** | Fragment support (2.3), a confident site (2.4), and a candidate explanation at that site (2.5) | Mass, site, site probability, explanation; for isobaric cases, the candidate set and secondary evidence |
| **Localised, unexplained** | Fragment support and a confident site, but no candidate explanation there | Mass, site, site probability |
| **Ambiguous** | No fragment support, or no confident site | Mass only |

No precursor-level cause is claimed for ambiguous PSMs. They make up most of the pool on fecal
(section 1.3), and that is acceptable. Whether a state is exported as an identification is decided
by the evaluation (4.3), not by this assignment.

### 2.7 What makes an explanation error-controlled

The census's chance control prices the mass match. The decoy sites (2.4) price the localisation.
Together they give an explanation-level error rate per class. This is separate from the
identification FDR, which target-decoy and entrapment already control (P6).

A shift becomes an exported identification only through the rescoring model (P1). The evidence
from this section enters as features, and each explanation class is its own FDR group (as P3
treats query types).

## 3. Proposed components

| Component | Content |
|---|---|
| `ShiftLocaliser` (core, new) | Per PSM: *S(H0)*, *S(Hk)* for every site, fragment support (2.3), site posterior and best site (2.4). Same mzML pass as `PredictedSpectrumFeatures`. |
| `ShiftExplainer` (core, new) | Per PSM: state (2.6) and candidate explanations at *k\** (2.5), using `SubstitutionMassTable`, `IsobaricModificationTable`, Unimod (`UnimodCatalog`) and protein flanks. |
| `ExplanationFdr` (core, new) | Chance control by mass offsets and decoy-site FLR per explanation class (2.7). |
| `JointFeatures` additions | `shift_fragment_support`, `shift_site_prob`, and state indicators (localised-explained, localised-unexplained, ambiguous). |
| Export | Every unexplained-shift PSM carries its state, mass, and, where assigned, site, site probability and explanation, with the error rate of the explanation. If the evaluation (4.3) supports it, `OPEN_LOCALISED_EXPLAINED` and `OPEN_LOCALISED_UNEXPLAINED` become identification classes, each its own FDR group. Ambiguous PSMs stay in `OPEN_UNEXPLAINED`. |
| Reused, report-only | `VariantFormRanker` (extended to modification forms), `PangenomeVariantCheck`, `SequenceVariantCaller`, as secondary evidence for isobaric cases. |

## 4. Evaluation (mock first, P13)

### 4.1 Localisation accuracy against truth already in the data

1. **Oxidation:** spectra the closed search identified with variable M oxidation on a single
   possible site. Those same spectra in the open search must localise to that M.
2. **Extensions:** the protein context is the truth. The site must be the extended terminus.
3. **Substitutions (12mix):** the true sequence is known when the substituted peptide exists in one
   of the 12 mix genomes.
4. **Decoy-site FLR (2.4):** it must agree with the error rate measured in 1–3.

### 4.2 How the pool divides

Report the share of PSMs and of new peptides in each state (2.6), per run, with the explanations
that dominate the localised-explained state. Ambiguous PSMs are counted and not explained further.

### 4.3 Identification gain

A harness arm adds the new features and FDR groups to Btm. It is read at every level with both
readings (P6) and gated by P7 and P13:
- calibration per new class;
- new peptides at matched entrapment FDP on fecal;
- no loss elsewhere.

The explanation error rate per class is reported next to its counts.

### 4.4 Mock tools

| Tool | Content |
|---|---|
| `LocScore.java` (scratch) | *S(H0)*, *S(Hk)*, support, posterior, best site per census PSM. Extends `FragFeat.java`. |
| `explain.py` (scratch) | States, explanations at the site, chance and decoy-site rates, and the checks in 4.1 and 4.2. |
| `harness.py --shift` | Joins the features and FDR groups for 4.3. |

## 5. Order

1. `LocScore.java` and `explain.py` on both pseudo-reversed runs.
2. Localisation accuracy (4.1) and how the pool divides (4.2).
3. The harness arm (4.3), with every level reported.
4. Only if it passes P13: the production components (section 3) with unit tests.

## 6. Mock results (2026-10-25)

The tools are `scripts/identification_audit/LocScore.java`, `prep_loc.py`, `explain.py` and
`census.py`. They ran on production v2's unexplained-shift PSMs that pass their own 1% walk.
Thresholds: fragment support ≥ 2, site probability ≥ 0.75.

### 6.1 Localisation accuracy

The figures are the share of confident sites that are correct.

| Truth set | fecal n | fecal | 12mix n | 12mix | Decoy-site estimate vs measured error (fecal / 12mix) |
|---|---|---|---|---|---|
| Extension (the protein context gives the terminus) | 1,073 | 84.7% | 370 | 76.6% | 15.6 vs 15.3% / 25.0 vs 23.4% |
| Substitution to a database peptide (its position) | 1,275 | 84.4% | 822 | 84.5% | 9.0 vs 15.6% / 11.6 vs 15.5% |
| Oxidation on a single M | 18 | too few | 6 | too few | — |

- The decoy-site estimate matches the measured error for extensions.
- For substitutions it is lower than the measured error. It cannot see confusion between two
  allowed sites, such as two residues that both fit a CH2 difference.

### 6.2 How the pool divides

| Run | State | PSMs | Share | New peptides | Entrapment FDP |
|---|---|---|---|---|---|
| fecal | Localised, explained | 5,637 | 12.1% | 1,035 | 0.68% |
| fecal | Localised, unexplained | 3,379 | 7.3% | 956 | 4.49% |
| fecal | Ambiguous: no fragment support | 30,824 | 66.2% | | |
| fecal | Ambiguous: site not confident | 6,704 | 14.4% | | |
| 12mix | Localised, explained | 5,414 | 3.4% | 706 | 0.52% |
| 12mix | Localised, unexplained | 3,897 | 2.4% | 967 | 7.66% |
| 12mix | Ambiguous: no fragment support | 130,594 | 81.3% | | |
| 12mix | Ambiguous: site not confident | 20,819 | 13.0% | | |

The shares change little for support thresholds of 1–3 and site-probability thresholds of 0.75–0.9.

### 6.3 What the results show

- **Explained PSMs on fecal are mostly isobaric.** At most fecal localised-explained sites, both a
  substitution and a Unimod entry fit the same residue (4,530 of 5,637), so the explanation stays a
  candidate set.
- **12mix has single-Unimod explanations:** AEBS (765) and Delta:H(4)C(3) (571).
- **A confident site with nothing known there marks doubtful identifications:** FDP 4.5% on fecal
  and 7.7% on 12mix, against 0.5–0.7% when an explanation fits. This should be a rescoring feature.
- **Many Unimod mass matches are coincidences.** For PSMs that only a Unimod mass explains, the
  decoy-site estimate of wrong sites is 66% on fecal and 36% on 12mix, against 16–17% for
  substitutions and extensions.

## 7. Harness arms (2026-10-25)

Every shifted harness candidate was localised by the same label-blind rule: 1.84M on fecal and
3.30M on 12mix. Candidates whose shift is zero within isotope spacings have no localisation
features. The arms build on Btm (RESET, precursor queries, monotone trees):

| Arm | What changes |
|---|---|
| BtmL | `loc_support`, `loc_site_prob` and the state flags added as model features |
| BtmLx | BtmL, plus localised-explained winners exported as class `LOC_EXPLAINED`, its own FDR group |
| BtmX | Btm's own model and winners; localised-explained winners exported as `LOC_EXPLAINED`, its own FDR group |

For the genome test, `LOC_EXPLAINED` counts as non-evidence.

### 7.1 Real identifications at matched entrapment FDP 1%

| Run | Arm | PSMs | Peptides | Proteins after filter |
|---|---|---|---|---|
| fecal | Btm | 525,142 | 65,014 | 15,559 |
| fecal | BtmL | 522,648 (−0.5%) | 64,801 (−0.3%) | 15,301 (−1.7%) |
| fecal | BtmLx | 528,428 (+0.6%) | 65,758 (+1.1%) | 15,467 (−0.6%) |
| fecal | **BtmX** | **530,906 (+1.1%)** | **65,971 (+1.5%)** | **15,705 (+0.9%)** |
| 12mix | Btm | 800,675 | 89,477 | 15,793 |
| 12mix | BtmL | 799,389 | 89,292 | 15,681 |
| 12mix | BtmLx | 805,226 | 89,965 | 15,772 |
| 12mix | **BtmX** | **806,499 (+0.7%)** | **90,155 (+0.8%)** | **15,886 (+0.6%)** |

### 7.2 Calibration (BtmX)

- **Every level of the export chain (fecal):** the entrapment FDP is at or below the decoy FDR.
  PSMs are at 0.97% / 1.00%, final peptides at 1.00% / 0.98%, and proteins after the filter at
  0.33% / 0.26%.
- **12mix:** every level matches Btm, and all 12 mix genomes are accepted.

The new class on its own (decoy-based FDR / entrapment FDP, at its own 1%):

| Run | `LOC_EXPLAINED` | `LOC_UNEXPLAINED` (not exported) |
|---|---|---|
| fecal | 5,497 PSMs, 0.99% / 0.26% | 2,824 PSMs, 0.92% / 2.13% |
| 12mix | 6,012 PSMs, 0.97% / 1.87% | 3,308 PSMs, 0.97% / 1.03% |

### 7.3 What the arms show

- **Localisation helps as an explanation, not as a feature.**
  - As features it costs a little.
  - Exported as an explained class with its own FDR group and the model unchanged (BtmX), it adds
    1.5% peptides and 0.9% proteins on fecal, with calibration intact.
  - This passes the P13 adoption rule.
- **`LOC_UNEXPLAINED` is not calibrated on fecal** (2.1% FDP at 0.9% FDR), so it is not exported.
- **12mix shows one limit.** `LOC_EXPLAINED` is optimistic there (1.87% FDP). Its entrapment hits sit
  mainly on adduct shifts that decoys never take: +40.03 (C3H4, 8 entrapment, 0 decoys) and +26.02
  (C2H2, 12 entrapment, 0 decoys).

## 8. Allow-list shifts must be supported, on an allowed residue (2026-10-26)

The allow-list classes are oxidation, deamidation, pyroGlu, water loss, Na-H, Ca-2H, carbamidomethyl and the
other `ExpectedPeakSelector` entries. They were exported at the search engine's site with no residue check, and
`DeltaMassEvidenceGate` exempts allow-list peaks from its no-support folding. In the 12mix final run, 40% of the
+0.984 "deamidation" PSMs had no N or Q, and none had fragment support. Most were observed at about +1.97 Da.

### 8.1 Rule

Each `OPEN_KNOWN_SHIFT` spectrum winner is checked with `ShiftLocaliser`, using the folded shift (the peak's
mass, not the isotope offset). The identification and its FDR group are unchanged; only what is claimed about
the shift changes:

| Fragments | State | Exported |
|---|---|---|
| No support (support < 2) | `unsupported` | The unmodified peptide; the observed shift stays in `Mass difference` and `normalized_mass_diff` |
| Support, and an allowed residue within a factor of 10 of the best site | `localised_explained` | The modification on the best allowed residue, with its exact mass and explanation |
| Support, and every allowed residue more than 10 times less likely than the best site | `localised_unexplained` | The shift on its fragment-supported site, with no explanation |

**Allowed residue:** a substitution or UNIMOD entry of that mass lists the residue or terminus, as in section 2.5.

**Why a factor of 10 and not the site probability:** sites the fragments cannot tell apart share the probability.
For example, residues 1 and 2 differ only by b1. A correct N-terminal pyroGlu then has a probability of about 0.5
on its residue. The test asks whether the fragments contradict an allowed residue, not whether they single it
out.

### 8.2 Mock

Production v2 known-shift spectrum winners at q ≤ 1%, on the pseudo-reversed entrapment runs (fecal / 12mix):

| State | Share | Entrapment FDP |
|---|---|---|
| Unsupported, exported unmodified | 26% / 35% | 3.9% / 2.0% |
| Allowed residue consistent | 48% / 47% | 1.3% / 0.4% |
| Allowed residues contradicted | 26% / 18% | 0.8% / 0.8% |

- **Unsupported:** mostly +0.98 and +0.97, and on 12mix also +21.98 and +57.02. 81% (fecal) and 91% (12mix) of
  their peptides are identified unmodified elsewhere.
- **Contradicted:** mostly −17.03 on internal residues and +21.98 (Na-H) on V, A, L or I.

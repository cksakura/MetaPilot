# Design draft: confirmed substitutions apart from undetermined forms (2026-11-09)

Status: the 2026-11-09 mock argues against the confirmation search (section 7). Nothing here is implemented.

## 1. Problem

Finding amino acid substitutions against a whole-genome catalogue is a headline feature of the
open-search workflow. Today no substitution reaches the main result, and the folder that should
hold substitutions holds two different things.

**No substitution is ever identified.**
- A PSM whose shift the fragments place on a residue that a substitution could change becomes
  `OPEN_SUBSTITUTION_SITE`, which is never an identification (`JointIdentification.identifies`).
  On 12mix 12-19% of these lay outside the true genomes, against 2-4% for identified groups, and
  the one-residue variant was the true peptide in only 5 of 378 (2026-11-08).
- `sequence_variants.fasta` and `substitution_candidates.fasta` are written "for a follow-up
  search". No step runs one.
- The UHGP-100 pangenome check and the MS2 form ranking are report-only by design.
- A variant peptide reaches the main result only when its own sequence is in the searched
  database. It is then an ordinary peptide, and nothing marks it as a substitution.

**`sequence_variants/` mixes two concepts.**
1. Substitution evidence: variant calls, the UHGP-100 tables, the form ranking.
2. Undetermined forms: the side quantification takes every PSM `SubstitutionReassigner` flagged,
   which needs only that MSFragger's localisation fits some substitution. In the 12mix test
   (2026-11-09):
   - 134 of the quantified "variants" were the same peptide and shift as a main-run form, from
     spectra whose fragments could not place the shift (237 ambiguous, 41 unsupported, 33
     localised but unexplained);
   - 107 were an identified peptide read with one residue changed (`A[+42.0469]ALESTL...` is
     `LALESTL...`).

The user's rule settles where each belongs: a PSM enters the peptide, protein and genome results
only if its sequence is determined.
- A substitution that the data determine is a determined sequence. It belongs in the main result,
  marked as a substitution.
- A spectrum whose forms cannot be told apart belongs in a separate place, with each form listed.

## 2. What "confirmed" must mean

A substitution is confirmed when its variant peptide wins its spectrum in a search where it
competes as an ordinary peptide, at error-controlled FDR. That is the same standard as every
other identification. Each principle it has to meet:

- **P10 (combined, not chosen) and P1 (rescore).** The variant peptide is one more candidate in
  the joint model, competing in its precursor query (P3) with the reference peptide plus a shift
  and with every other closed and open candidate. Within a query no q value from another search
  decides (P2).
- **P11 (decoys make target errors).** Every variant target peptide gets a cleavage-preserving
  pseudo-reversed decoy, as the stage-2 databases do, so a decoy can win a spectrum the way a
  wrong variant can.
- **P12 (label-blind features).** UHGP-100 membership cannot be a model feature: reversed decoys
  are almost never in UHGP-100, so the feature would read the label. It can only gate which
  variant targets are built, before their decoys are made.
- **Group FDR.** Variant winners form their own group (`SUBSTITUTION`), walked separately like
  `EXACT`, `OPEN_UNMODIFIED` and `OPEN_ACCEPTED_SHIFT`, because the prior that a spectrum is a
  variant differs from the others.
- **P13 (adoption).** Mock first on cached outputs. The faecal entrapment run decides; 12mix is a
  sanity check. On 12mix, a confirmed variant should lie in a true strain's proteome: a sequence
  absent from every searched database but present in the truth is a real substitution.

## 3. Proposed flow

### 3.1 Candidate variant peptides

Built after the joint identification, from what the open search proposes:
- every winner whose fragments place the shift on a residue where a substitution fits
  (`OPEN_SUBSTITUTION_SITE`), and the `SubstitutionReassigner` proposals;
- the variant applied in the protein context: tryptic peptide plus flanks, so missed cleavages
  and terminal cleavage are right.

A variant target is kept only if its sequence is absent from the searched database. If present,
the closed search already competed it.

**Open question Q1: gate on the pangenome or not.**
- (a) Every proposal. The search space is small, a few thousand peptides on 12mix.
- (b) Only proposals whose variant peptide is in UHGP-100 (`PangenomeVariantLookup`). This
  uses the catalogue as the prior, which is the feature's point, and reduces the space further.
  The gate is sequence-only and applied before decoys are made, so it is label-blind.
- Recommendation: (b), with (a) measured once in the mock to show what the gate removes.

**Open question Q2: recurrence.** `SequenceVariantCaller` requires three samples for its
`SEQUENCE_VARIANT` class. As a gate it would remove single-sample variants and every variant in
a single-sample run. Recommendation: no recurrence gate. Recurrence stays a reported column.

### 3.2 Confirmation search

- **Database:** variant target proteins (one entry per variant and protein, header naming the
  reference protein, the site and the substitution), their pseudo-reversed decoys, and nothing
  else.
- **Search:** a closed MSFragger search of all spectra, with the stage-2 closed parameters (ppm
  window, isotope_error 0/1/2, the user's variable modifications).
- **Merge:** the results enter the joint candidate set as a new presence source, beside the
  closed and open searches. Searching all spectra, not only the flagged ones, lets a variant beat
  a reference peptide anywhere it fits better, and keeps the competition the same for targets
  and decoys.
- **Cost:** one more closed search, over a database far smaller than the stage-2 one (12mix
  v2final flagged 5,201 PSMs, so at most a few thousand variant peptides). Not yet timed.
- **Reuse:** a separate small database keeps the stage-2 search and its checkpoints untouched. The
  new search folder is self-describing, with its own params next to its pepXML.

### 3.3 Identification and downstream

- A query whose winner is a variant peptide passing the `SUBSTITUTION` group FDR is identified.
  Its PSM goes into `validated_psms.tsv` with `id_group=SUBSTITUTION`, the variant sequence as
  Base Sequence, and the variant protein entry as protein.
- **Quantification:** main FlashLFQ run, like every identification.
- **Proteins and genomes, open question Q3.** The variant protein maps to its reference protein's
  genome. Options:
  - (a) count it as evidence for that genome, since the strain present differs from the MAG at
    that site;
  - (b) protein evidence only, not genome evidence, because a variant can also come from a
    relative absent from the catalogue (the 12mix B3 finding).
  - Recommendation: (b) until the genome test is checked with it.
- **`VariantFormRanking`:** it compares the variant form with the modification form by predicted
  MS2. In the joint search both forms compete within one query, so its question is answered by
  the identification itself. Proposal: retire it once the confirmation search is adopted, or
  keep its score as a reported column only.

## 4. Result layout

```
opensearch_MAG_search/
  combined/                 main result; confirmed substitution PSMs carry id_group=SUBSTITUTION
  final_peptides.tsv ...    unchanged files, now including confirmed substitutions
  substitutions/            the feature: confirmed substitutions and their evidence
  ambiguous_forms/          spectra whose forms the data cannot tell apart
```

**Open question Q4:** folder names. `ambiguous_forms/` follows the user's earlier "peptideform"
idea; `undetermined_forms/` is the alternative.

### 4.1 `substitutions/`

| File | Content | From |
|---|---|---|
| `substitutions.tsv` | one row per confirmed substitution site: protein, site, substitution, variant peptide, PSMs, samples, q, reference peptide observed, UHGP-100 hit and species, codon class | new; joins `sequence_variants.tsv` columns and UHGP-100 tables |
| `substitution_candidates.tsv` | every variant target searched, confirmed or not, with why it was built | replaces `sequence_variants.tsv` as the candidate list |
| `substitution_search/` | variant database with decoys, `fragger_closed.params`, pepXML | new |
| `variant_candidates.tsv`, `variant_sites.tsv`, `variant_pangenome_summary.tsv`, `variant_species_summary.tsv`, `variant_columns.tsv` | UHGP-100 check, now of confirmed substitutions | moved from `sequence_variants/` |
| `substitution_peaks.tsv` | which fitted peaks are substitution-shaped | moved |

`sequence_variants.fasta` and `substitution_candidates.fasta` are replaced by the search
database in `substitution_search/`.

### 4.2 `ambiguous_forms/`

| File | Content | From |
|---|---|---|
| `ambiguous_forms.tsv` | one row per spectrum whose sequence is not determined, and per candidate form: base peptide plus shift, each substitution reading, each modification reading, with the fragment evidence for each | replaces `substitution_forms.tsv` |
| `ambiguous_psms.tsv`, `ambiguous_quantified_peptides.tsv`, `ambiguous_quantified_peaks.tsv` | the side quantification, as now | renamed from `substitution_*` |

Two fixes found in the 2026-11-09 test apply here:
- **Same form as the main run.** A spectrum whose base peptide and shift (within the explanation
  tolerance) match a form in `validated_psms.tsv` is left out of the side input. That form's peak
  is already quantified.
- **Monoisotopic reassignment.** A spectrum whose fragments support the shift one or two ¹³C steps
  below its raw delta has its monoisotopic mass corrected first (separate item, measured before
  adoption).

## 5. Validation plan

1. **Mock (P13).** From the round-5 12mix and faecal outputs, build the variant targets and decoys
   for Q1 (a) and (b). Run only the confirmation search on the HPC, then score with the audit
   harness and the same export chain.
2. **Faecal entrapment decides.** Entrapment FDP of the `SUBSTITUTION` group and overall, peptide
   and PSM level. Calibration gate (P7) first.
3. **12mix sanity.** No level drops by more than 1%, all 12 genomes still accepted, and the share
   of confirmed variant peptides found in a true strain's proteome.
4. **Expectation.** Few confirmations: the variant reading was true in 5 of 378 cases on 12mix.
   A small, clean set is the feature working. A large one would need explaining before adoption.

## 6. Order of work

1. Agree Q1-Q4.
2. Mock the confirmation search (section 5, steps 1-3).
3. If it passes: implement the search step, the `SUBSTITUTION` group and the two folders. Move the
   files, with migration of earlier runs' copies as `inSequenceVariants` does now.
4. Same-form exclusion and monoisotopic reassignment for `ambiguous_forms/`, each measured on its
   own.

## 7. Mock result (2026-11-09): confirmation search rejected

Jobs 24800208 (12mix) and 24800212 (faecal). Variant proteins whole, pseudo-reversed decoys,
twin-matched entrapment variants, no UHGP-100 gate. Searched closed with the stage-2 database.

| | 12mix | faecal |
|---|---|---|
| variant-only target PSMs (peptides) | 10,200 (1,883) | 10,486 (1,649) |
| decoy winners | 12 | 0 |
| entrapment-twin PSMs (peptides) | 100 (12) | 64 (8) |
| overall PSM entrapment FDP, round 5 -> mock | 0.314% -> 0.346% | 0.489% -> 0.509% |
| 12mix off-truth PSMs, round 5 -> mock | 2.52% -> 3.77% | - |
| real variant peptides in a true strain proteome | 10 of 1,871 | - |

- Every variant-only winner checked covers the substituted site; the mapping is right.
- 12mix truth: 8,829 of the 10,100 real variant PSMs (1,600 peptides) are spectra of a true
  reference peptide. The variant reads its modified form as a substitution: G->N +57.021, A->Q
  +57.021, A->S +15.995, S->D +27.995, V->L +14.016.
- 1,895 of these spectra were correctly identified in round 5 (934 as the reference with an
  accepted shift, the others as other true peptides). The variant took them over.
- Decoys (12 and 0) and entrapment twins (about 1-2%) do not see this error. A wrong variant is a
  neighbour of a true peptide, not a random match. So faecal entrapment would have passed it,
  while 12mix truth shows about 0.5% precision.
- A Unimod gate does not rescue it. Variants with no Unimod mass on the substituted residue were
  still 8 of 1,736 peptides in truth. The shift usually sits on a neighbouring residue.
- Cause: the closed variant candidate carries closed-search features (found_by_closed, a narrow
  window), which the open "reference plus shift" reading lacks. The joint model favours the closed
  candidate, so the competition is not like for like.

Recommendation (user decision pending): do not adopt. Substitutions stay reported evidence, as now. A future confirmation needs a
like-for-like test of the variant against the reference with the same shift on every residue,
decided by the site-determining fragments. It also needs an error estimate that sees neighbour
wins, as decoys and twins do not.

### 7.1 Post-hoc UHGP-100 gate (job 24810280, 2026-11-09)

The pipeline's own scan (`PangenomeVariantLookup.scan`, tryptic context required) over
uhgp-100.faa.gz (95.3 M proteins) for all 9,609 variant peptides of both mocks.

| 12mix real variant winners | peptides | in UHGP-100 | PSMs kept |
|---|---|---|---|
| in a true strain proteome | 10 | 10 | 19 of 19 |
| not in a true strain proteome | 1,861 | 839 | 5,151 of 10,081 |

- The gate keeps every true variant but only halves the false ones: peptide precision goes from
  0.5% to 1.2% (10 of 849).
- Searched variants pass at about the same rate as winners (863 of 2,272): the catalogue holds a
  one-residue relative of a common peptide about 40% of the time, so membership says little about
  whether the spectrum is a variant.
- Faecal: 1,228 of 1,641 winners pass. Entrapment twins are shuffled sequences and almost never in
  UHGP-100 (0 of 12 and 2 of 8 winners), so their rate after the gate is no error estimate.

The gate does not rescue the confirmation search; the recommendation in section 7 stands.

### 7.2 Caveat on the 12mix truth (2026-11-10)

The 12mix "truth" is the catalogue's representative genome for each of the 12 species
(`scripts/open_vs_closed/12mix_genomes_mapping.tsv`), not the sequenced strains themselves. A
correct variant that matches the real strain but not the representative genome counts as "not in
truth". The precision figures in 7 and 7.1 (0.5%, 1.2%) are therefore not valid measures for
substitutions. The modification-mass pattern still holds: G->N and A->Q at +57.021, A->S at +15.995.
Those are not plausible strain differences. Methyl-sized changes (V<->L, E<->D) cannot be judged this
way.

The valid test bed is PXD050781: ATCC single strains with UniProt strain proteomes, where the
variant caller's SEQUENCE_VARIANT calls were confirmed for 86-91% (2026-11-02).

# Plan: genome acceptance by expected false counts, reported at the rank the catalogue resolves

Status: version 8. Version 8 applies the review of version 7. Section 8 (identification: closed
search first, open search on the remaining spectra) is now its own document,
`closed-then-open-identification-plan.md`; this one keeps the genome test. Under E: the cause of
the lower `F` in E1 (the peak screen, not a reproducibility fault), the reversed-genome rows with
their own decoys left out of the window, and the E2 verdict. Step 5: the open search's reported
razor proteins follow the step 5 count rule. Interface: `Taxa.tsv` names its open-search column
`kept_unique_peptides`. Section 4 and risk 5 carry the E1 count (14 of 20); step 2 states how `K`
was set.
Versions 5-7 added section 8, from comparing the first build with the closed searches on 12mix
(version 6: what the closed search's answers show about the open search's peak steps; version 7:
the PXD050781 replicate check, the gate's cost there, a per-dataset stop rule, and the offline test
of closed-contradicted PSMs as classifier negatives).
Version 4 adds what the fourth audit asked for before the build: the uniqueness
rule across runs and decoy uniqueness from the digested reversed proteomes (step 2), the written
class choice for `n_g`, `F` and `pi_g` (step 2), the interface to razor assignment and the output files
(section 3, "Interface"), the calibration condition in the B2 rule, and the step 5 calibration note;
and the results of work items B1-B4 (section 2.7), which fix the step 3 neighbourhood at one
substitution.
Version 2 merged the first version (catalogue resolution as an optional report) with the audit's
alternative (expected false count per genome, near-neighbour exclusion, picked protein-group FDR).
Version 3 applies the second audit: step 3 wording, step 4 conditioning, step 5 decoy grouping,
engine scores in step 1, scope of the first build, validation led by 12mix; and it changes how
`pi_g` is estimated after testing the audit's answer to open question 1 (section 2.5).
After version 3, steps 1-3 were prototyped on 12mix (section 2.6): all 12 spiked genomes accepted,
and 15 of 20 non-spiked genomes too. The third audit decided the questions this reopened (section 6:
acceptance stays steps 1-3, a `same_genus_as` label, *Bacillus_A* reported separately); section 7
is the next work plan.
Scope: DDA MAG genome scoring (`DdaMagGenomeScoringService` and its callers). The first build covers
the MetaPilot open search; pFind and FragPipe follow through the same post-analysis service once
they supply what step 1 needs.

## 1. Problem

MetaPilot decides which genomes are present from peptides that are *unique to one genome in the
database the run searched*. Two things make that decision unreliable for any genome except the
dominant ones:

1. **The genome score test cannot fail a genome with one or two good unique peptides.**
   On PXD050781 (21 runs of a pure culture of *Blautia hydrogenotrophica* DSM 10507, target
   genome MGYG000000035) all 39 candidate genomes reach the resampling floor, p = 1/10,001. The
   test is calibrated: the 35 reversed proteomes run through the same statistic give 0 decoy
   genomes at p <= 0.01. So it controls random matches; it does not see a real peptide attributed
   to the wrong genome.
2. **"Unique" is relative to the searched database, and that database is a selection.**
   The open search builds one database per run from its marker step
   (`search_results/<run>/marker/<run>.fasta`); this is the design. Each of the 39 genomes has
   unique PSMs in exactly the runs whose database contained it and in no other run. When a genome
   was not searched, its spectra were still there and were taken by a relative (section 2.2).

The method must work on a **single run**: real datasets rarely have replicates. A run-recurrence
test across replicates was implemented and then removed, because with per-run databases it
measures database inclusion, not reproducibility (section 2.1).

## 2. Evidence (PXD050781 and 12mix, MSFragger-only open search)

Scripts and tables are in the agent store
`C:\Users\daz25hak\AppData\Local\Cursor\AgentStores\cursor_agent_stores\090a6a6e-4a26-4625-a770-6c1635d04c78\files\pxd050781_eval\`
(scripts `26`-`38`, outputs in `out/`). Input: `Z:\PRIDE_dataset\PXD050781_haonan\DDA_10507\opensearch_MAG_search\msfragger_only`.

PXD050781 has no ground truth for low-level non-*Blautia* material: carry-over from gut samples run
on the same instrument is as plausible as misattribution. It shows how the statistics behave; it
cannot say which accepted genomes are wrong. 12mix decides that (section 4).

### 2.1 Run recurrence measured database inclusion (`31_genome_run_profile.py`)

Unique exact-mass PSMs per run are near-constant whenever the genome was searched and zero
otherwise, e.g. MGYG000002455 (*B. cellulosilyticus*): 187, 180, 178, 179 PSMs in runs 1, 4, 11, 21,
the only runs whose database had it. The 11 genomes the recurrence test kept were the 11 searched
in most runs. Restricting the test to runs that searched each genome would pass all 39.

### 2.2 Database swap (`32_database_swap.py`)

For the 31 genomes missing from some runs' databases: their exact-mass unique peptides were looked
up in those runs (same charge, 10 ppm or one isotope above, run-aligned RT within 30 s).

| Measure | Range over 31 genomes |
|---|---|
| Spectrum found | 41-81% of lookups |
| Found by chance (RT shifted 3 min) | 8-18% |
| Spectrum validated as another peptide | 19-72% of lookups |

Validated replacements by relation to the missing genome (pooled): same genus 66%, same family
18%, same order 7%, same class 1%, other phylum 9%. The replacement usually differs by <= 2
residues and is often matched with an open-search mass shift. The spectra are real and present in
every run; which organism they come from is the open question.

### 2.3 Catalogue uniqueness (`34_catalogue_uniqueness.py`)

Query: 22,172 validated, exact-mass (|dm| < 0.05 Da) peptides unique to one genome in the run
database, I = L. Searched in all 95 M UHGP-100 representatives (tryptic context: preceded by K/R,
protein start or Met1; ending K/R or at the protein end), expanded to member genomes through
`uhgp-100.tsv`, mapped to species and GTDB lineage with `genomes-all_metadata.tsv` (4,744 species).
Null: decoy peptides unique to the genome's reversed proteome under the same filters.

| Genome | Observed | Decoys | Species-unique | Genus-unique | Family-unique |
|---|---|---|---|---|---|
| *Blautia_A hydrogenotrophica* (target) | 12,849 | 130 | 5,806 | 6,760 | 9,040 |
| *Escherichia coli_D* | 2,570 | 118 | 16 | 213 | 895 |
| *Phocaeicola dorei* | 1,111 | 156 | 0 | 28 | 161 |
| *Blautia_A wexlerae_A* | 691 | 114 | 0 | 3 | 172 |
| *Bacteroides cellulosilyticus* | 212 | 9 | 0 | 22 | 50 |
| *Lachnoclostridium_A edouardi* | 115 | 57 | 31 | 32 | 49 |

Full table: `out/catalogue_resolution.tsv`. Using decoys as the expected false count and
Poisson p <= 0.01: 1 genome resolves to species (the target), 2 to genus, 10 to family, 26 to none
of these ranks. Median number of species carrying a non-target genome's peptides: 4-156.
These counts are not yet conditioned on the target species (step 4, item 3), so genus and family
counts of *Blautia_A* and *Lachnospiraceae* relatives are upper bounds.

Runtime of the Python prototype (6 processes, uncompressed tar on the network share): 490 s over the
41.6 GB FASTA, 228 s over the 25 GB membership table.

### 2.4 Expected false count with near-neighbour exclusion (`35_expected_false_count_test.py`)

Steps 2 and 3 on the same data: neighbours (identical or one substitution) taken from the full
proteomes of all 39 candidates (as step 3 now specifies, for accepted genomes), decoy side not
reduced, `pi_g` from rank-1 decoys at all scores (the version 2 definition, see 2.5).
`F` = 1,780 accepted decoy peptides unique to one reversed genome.

| Genome | Unique | Kept after exclusion | Expected false | p | Catalogue resolves to |
|---|---|---|---|---|---|
| *B. hydrogenotrophica* (target) | 12,849 | 12,849 | 103 | ~0 | species |
| *E. coli_D* | 2,570 | 2,542 | 130 | ~0 | genus |
| *Phocaeicola dorei* | 1,111 | 1,088 | 141 | ~0 | none |
| *F. prausnitzii_C* | 758 | 746 | 81 | ~0 | none |
| *Gemmiger* sp900540595 | 403 | 331 | 64 | 1e-122 | none |
| *Collinsella* sp900762015 | 154 | 146 | 52 | 1e-26 | none |
| *Bacteroides stercoris* | 85 | 13 | 11 | 0.32 | family |
| *Bacteroides oleiciplenus* | 60 | 5 | 7 | 0.84 | family |

Accepted at p <= 0.01: 24 of 39. Exclusion removes the *Bacteroides* group members that the swap
showed absorbing each other's spectra. The other accepted genomes are typical gut taxa with hundreds
of peptides against an expected false count of tens; whether that is carry-over or misattribution
cannot be settled here (see the note at the top of section 2). Full table:
`out/expected_false_count_test.tsv`.

### 2.5 Where false matches land: the estimate of `pi_g` (`36`-`38`)

Test of open question 1: does the share of decoys unique to `rev_g` at all scores predict the share
among accepted decoys? Accepted set: 1,780 distinct exact-mass decoy peptides unique to one
reversed genome. Chi-square over genomes with expected >= 5.

| `pi_g` from | chi2 / df | p | Spearman(obs/exp, runs searched) |
|---|---|---|---|
| all rank-1 decoys, by PSMs | 161 / 33 | 8e-19 | |
| all rank-1 decoys, distinct peptides | 295 / 35 | 7e-43 | 0.90 |
| per run, distinct peptides, then pooled (`sum_r F_r pi_g,r`, recurrence-corrected) | 155 / 32 | 5e-18 | 0.86 |
| top 1F decoy peptides by score | 273 / 23 | 1e-44 | |
| top 3F | 150 / 22 | 5e-21 | |
| top 10F | 79 / 26 | 3e-7 | |
| top 30F (class_logit >= -9.3) | 34 / 27 | 0.16 | 0.10 |

By runs searched (all-score, distinct peptides): accepted / expected = 0.26 for genomes searched in
1-5 runs, 0.60 for 6-10, 1.09 for 11-15, 1.22 for 16-21. Estimating per run does not remove it, so
it is not a pooling artifact: within a run, high-scoring decoys fall more often on the reversed
proteomes of genomes that are well represented in the data. With all-score `pi_g`, `lambda` is
about 4x too high for marginal genomes (conservative, costs sensitivity) and about 1.2x too low for
the well-represented ones. Two further results:

- Each accepted decoy peptide recurs in 2.15 runs on average (3,834 per-run vs 1,780 distinct), so
  false peptides are reproducible across runs. Distinct-peptide counting on both sides is required,
  and recurrence across runs is not evidence of presence.
- The top-30F window is a fitted constant on one dataset; step 2 therefore reports the agreement on
  every run (calibration rows and the chi-square in the run-info), and 12mix sets the window.

Tables: `out/decoy_share_check.tsv`, `out/decoy_share_by_score.tsv`, `out/decoy_share_per_run.tsv`.

Caveat: `scored_psms.tsv` lists one protein per hit, so these scripts took the all-score decoys'
genome from that one protein without checking that the peptide is unique to it; the accepted decoys
(from `validated_psms.tsv`) were checked against all their proteins. The 12mix prototype (2.6)
decides uniqueness from the digested proteomes on both sides and finds the same pattern.

### 2.6 12mix: steps 1-3 against known composition (`mix12_eval/40_genome_test_prototype.py`)

Data: `Z:\PRIDE_dataset\IPX0003851000\IPX0003851001\DDA\MetaPilot_new\opensearch_MAG_search\msfragger_only`
(search outputs of 2026-10-11; the genome decision is recomputed from `scored_psms.tsv`, so the
genome scoring used in that run does not matter). Truth: `12genomes.txt` (12 spiked genomes).
9 runs, 32 candidate genomes in the union of the per-run databases, all 12 spiked genomes among them.
Outputs in the agent store `...\files\mix12_eval\out\`.

Prototype, as specified in section 3:
- Step 1 on rank-1 hits: 645,251 unmodified, 96,668 known-shift and 511,518 unexplained PSMs;
  accepted peptides at q <= 0.01 (targets / decoys): unmodified 77,013 / 770, known shift
  22,047 / 220, unexplained 25,035 / 250 (not used as evidence).
- Uniqueness from the digested proteomes (strict trypsin, 2 missed cleavages, length 7-50, Met
  clipped, I = L) within each run's database, the same for targets and decoys. Target evidence
  peptides: 43,967 unique to one genome in every run they were seen, 23,683 unique to different
  genomes in different runs (dropped), 11,322 shared, 1,122 not in the digest. `F` = 802.
- Steps 2-3 as specified; acceptance in ascending `log_p` with the Benjamini-Hochberg threshold at
  each position.

| Method | Spiked accepted | Non-spiked accepted | Decoy calibration chi-square p |
|---|---|---|---|
| Current genome score test (`genome_rank_scores.tsv`) | 12 of 12 | 20 of 20 (all at p = 1e-4) | |
| `pi_g` from top 10F decoys | 12 of 12 | 17 of 20 | 0.019 |
| `pi_g` from top 30F decoys (plan default) | 12 of 12 | 15 of 20 | 0.009 |
| `pi_g` from all-score decoys | 11 of 12 | 8 of 20 | 1e-27 |

The all-score row rejects more only because its `lambda` is too high for genomes searched in few
runs (as on PXD050781); that is also how it loses the spiked *Enterococcus_D casseliflavus*.

Non-spiked genomes accepted with `K = 30` (kept peptides after step 3 / expected false count):

| Genome | GTDB species | Kept / lambda | Spiked relative |
|---|---|---|---|
| MGYG000002283 | *Bacillus_A paranthracis* | 485 / 55 | none in the family (possible contaminant) |
| MGYG000000323 | *Enterococcus_D* sp002850555 | 419 / 13 | *E. casseliflavus* (35 unique, 10 kept) |
| MGYG000002511 | *Klebsiella variicola* | 199 / 22 | *K. pneumoniae* |
| MGYG000002350 | *Serratia marcescens_I* | 68 / 3.3 | *Enterobacteriaceae* |
| MGYG000000236 | *Bacteroides fragilis_A* | 41 / 3.8 | *B. fragilis* |
| MGYG000000007 | *Lactobacillus intestinalis* | 21 / 2.9 | *L. acidophilus* |
| MGYG000002323 | *Escherichia* sp000208585 | 20 / 1.7 | *E. coli_D* |
| MGYG000000902 | *Clostridium neonatale* | 19 / 0.7 | *C. butyricum* |
| MGYG000001553 | *Clostridium_P mediterraneense* | 18 / 2.9 | *Clostridiaceae* |
| MGYG000002502 | *Enterobacter chengduensis* | 16 / 3.0 | *E. asburiae_B* |
| MGYG000002536 | *Plesiomonas shigelloides* | 11 / 0.6 | *Enterobacteriaceae* |
| MGYG000002515 | *Escherichia fergusonii* | 9 / 0.2 | *E. coli_D* |
| MGYG000000110 | *Kluyvera* sp902363335 | 6 / 1.5 | *Enterobacteriaceae* |
| MGYG000002288 | *Bacillus_A luti* | 5 / 0.2 | *B. paranthracis* (above) |
| MGYG000002353 | *Enterococcus_B faecium* | 4 / 0.4 | *Enterococcaceae* |

Rejected: *Serratia surfactantfaciens* (26 / 17), *Enterobacter mori*, *E. cloacae*,
*E. cloacae_K*, *Enterococcus_B lactis*.

By rank, the 15 non-spiked genomes accepted with `K = 30`:

| Relation to the spiked strains | Genomes | Error |
|---|---|---|
| Same genus as a spiked strain | 8: *Enterococcus_D* sp002850555, *K. variicola*, *B. fragilis_A*, *Escherichia* sp000208585, *E. fergusonii*, *L. intestinalis*, *C. neonatale*, *E. chengduensis* | species |
| Same family only | 5: *Serratia marcescens_I*, *Plesiomonas*, *Kluyvera*, *Clostridium_P*, *Enterococcus_B faecium* | genus |
| No spiked strain in the family | 2: *Bacillus_A paranthracis*, *B. luti* | possible contaminant and its relative (section 6, Q6) |

Step 3 removes 5-82% of an accepted relative's peptides but leaves 4-485. The *E. casseliflavus*
case is the extreme: the relative *Enterococcus_D* sp002850555 (in 8 of 9 run databases) holds 442
unique peptides, the spiked genome (in 5 of 9) holds 35, and the relative is accepted 11th, the
spiked genome 29th.

**Reading.** The decoy null accepts these genomes correctly as far as random matches go: their
kept peptides exceed `lambda` by factors of 4-40. *Hypothesis, not yet tested:* the evidence is real
spectra from the spiked strains, on peptides that the spiked strain carries but its reference genome
does not (strain variants). It is plausible for the 8 same-genus cases; it is not established for
the 5 family-level ones. Work item B3 measured it (section 2.8): not confirmed in this form.

Reversed genomes accepted: none at q <= 0.01 (section 2.7, B1).

On known composition, steps 1-3 reduce false acceptance from 20 of 20 to 15 of 20 without losing a
spiked genome; 8 of the 15 are species errors within a spiked genus, 5 are genus errors, 2 are the
*Bacillus_A* pair.

### 2.7 12mix: work items B1, B2 and B4 (`mix12_eval/41_b_measurements.py`)

Same data and steps 1-3 as section 2.6. Added: every reversed genome as a calibration row (step 2,
last item; its accepted decoy peptides excluded against the reversed proteomes of the accepted target
genomes other than its own counterpart, `lambda` scaled the same way, Benjamini-Hochberg over the
reversed genomes); the final q-value as Benjamini-Hochberg over the post-exclusion p-values; the
neighbourhood at Hamming distance 1 and 2 (equal length; three-segment pigeonhole join over the 6.7 M
target and 6.5 M decoy digest peptides, verified by Hamming distance); and both window class choices.
"Species error" = non-spiked genome in a spiked genus, "genus error" = in a spiked family only (GTDB
lineage from UHGG `genomes-all_metadata.tsv`). Outputs `out/b_summary.tsv`, `out/b_calibration.tsv`,
`out/b_genome_test_<classes>_<window>_d<1|2>.tsv`.

`log_p` is computed exactly: the regularized incomplete gamma function, and a log-space series where
that underflows. The section 2.6 prototype used `scipy`'s `logsf`, which is `-inf` for the strongest
genomes, so their order came from the tie-break. The counts in the table are the same either way.
Only the *Klebsiella* pair changes order: *K. aerogenes* (log10 p = -805) is now accepted before
*K. pneumoniae* (-619), and the two swap their shared neighbour peptides (806 and 538 kept, against
561 and 666 in section 2.6). The Java build computes `log_p` the same way.

| `pi_g` window | Substitutions | Spiked accepted | Species errors | Genus errors | *Bacillus_A* pair | Reversed genomes at q <= 0.01 | Window chi-square p |
|---|---|---|---|---|---|---|---|
| any class, top 10F | 1 | 12 | 10 | 5 | 2 | 0 | 0.019 |
| any class, top 10F | 2 | 11 | 9 | 4 | 2 | 0 | 0.019 |
| **any class, top 30F** | **1** | **12** | **8** | **5** | **2** | **0** | **0.009** |
| any class, top 30F | 2 | 11 | 8 | 4 | 2 | 0 | 0.009 |
| any class, all scores | 1 | 11 | 6 | 1 | 1 | 1 | 1e-27 |
| any class, all scores | 2 | 11 | 3 | 1 | 1 | 0 | 1e-27 |
| evidence classes (all 5,734) | 1 | 11 | 7 | 3 | 1 | 1 | 4e-24 |
| evidence classes (all 5,734) | 2 | 11 | 5 | 1 | 1 | 1 | 4e-24 |

The genomes where the running threshold of step 3 and the final q-value disagree: none in any row.

**B1, reversed genomes.** With the any-class top-30F window no reversed genome is accepted, under
either neighbourhood; each reversed genome's kept decoy peptides are close to its `lambda` (e.g.
reversed *P. aeruginosa* 167 against 185, *C. butyricum* 59 against 55). The evidence-class window
has fewer decoys than `10 * F`, so its three windows are one all-score window, biased as in section
2.5: reversed *P. aeruginosa* is accepted (168 against 99, q = 7e-9). This decides the class choice
of step 2: the window takes any class.

**B2, two substitutions.** Against one substitution (any class, top 30F): one genus error fewer
(*Kluyvera* sp902363335 rejected), no species error fewer, and the spiked *E. casseliflavus* lost
(10 kept peptides fall to 4, q = 0.29). *K. aerogenes* keeps 647 of 970 instead of 806 and
*K. pneumoniae* 383 of 762 instead of 538; both stay accepted. Reversed genomes accepted: 0 in both.
The adoption rule (section 7, B2) requires every spiked genome to stay accepted, so step 3 keeps
**one substitution**.

**B4, `same_genus_as` (any class, top 30F, one substitution).** The label names the accepted genome
with the most kept peptides in the same GTDB genus, when that is another genome:

| Labelled genome | `same_genus_as` | Spiked? |
|---|---|---|
| *K. variicola* (199 kept) | *K. aerogenes* (806) | no, species error |
| *B. fragilis_A* (41) | *B. fragilis* | no, species error |
| *C. neonatale* (19) | *C. butyricum* | no, species error |
| *Escherichia* sp000208585 (20) | *E. coli_D* | no, species error |
| *E. fergusonii* (9) | *E. coli_D* | no, species error |
| *L. intestinalis* (21) | *L. acidophilus* | no, species error |
| *E. chengduensis* (16) | *E. asburiae_B* | no, species error |
| *B. luti* (5) | *B. paranthracis* | no, *Bacillus_A* pair |
| *K. pneumoniae* (538) | *K. aerogenes* (806) | yes |
| *E. casseliflavus* (10) | *Enterococcus_D* sp002850555 (419) | yes |

It marks 7 of the 8 species errors. The eighth, *Enterococcus_D* sp002850555, outranks the spiked
*E. casseliflavus*, so the label lands on the spiked genome instead. Genus errors are not marked:
they have no accepted genome in their genus.

### 2.8 12mix: work item B3, catalogue carriage (`mix12_eval/42_catalogue_carriage.py`)

For each accepted non-spiked genome X: the fraction of its kept peptides that UHGG genomes of the
spiked species next to it (same genus, else same family) carry, against a control. The control is
tryptic peptides of the same proteins of X that were never a rank-1 target hit, are unique to X among
the candidate genomes, and lie outside the one-substitution neighbourhood of the other accepted
genomes, sampled to the same lengths. UHGP-100 scan with tryptic context and I = L, expanded to
cluster members for representatives of X's family; about 11 minutes. One-sided Fisher test.

| X | Paired spiked species | Kept peptides | Carried | Control carried | Fisher p |
|---|---|---|---|---|---|
| *Enterococcus_D* sp002850555 | *E. casseliflavus* | 419 | 99.3% | 99.0% | 0.5 |
| *K. variicola* | *K. pneumoniae*, *K. aerogenes* | 199 | 74% | 61% | 0.004 |
| *B. fragilis_A* | *B. fragilis* | 41 | 90% | 85% | 0.37 |
| *L. intestinalis* | *L. acidophilus* | 21 | 5% | 0% | 0.5 |
| *Escherichia* sp000208585 | *E. coli_D* | 20 | 90% | 88% (16) | 0.61 |
| *C. neonatale* | *C. butyricum* | 19 | 11% | 5% | 0.5 |
| *E. chengduensis* | *E. asburiae_B* | 16 | 50% | 38% | 0.36 |
| *E. fergusonii* | *E. coli_D* | 9 | 67% | 44% | 0.32 |
| *S. marcescens_I* (genus error) | spiked *Enterobacteriaceae* | 68 | 12% | 0% | 0.003 |
| *Clostridium_P* (genus error) | *C. butyricum* | 18 | 0% | 0% | 1 |
| *Plesiomonas* (genus error) | spiked *Enterobacteriaceae* | 11 | 9% | 0% | 0.5 |
| *Kluyvera* (genus error) | spiked *Enterobacteriaceae* | 6 | 0% | 50% | 1 |
| *Enterococcus_B faecium* (genus error) | *E. faecalis*, *E. casseliflavus* | 4 | 25% | 100% (2) | 1 |

The kept peptides were taken from the section 2.6 order; the relatives' kept counts are identical
under the exact `log_p`.

**Reading.** The hypothesis of section 2.6 is not confirmed in the form stated:

- Same genus: where carriage is high (sp002850555, *B. fragilis_A*, *Escherichia* sp000208585),
  the control is carried as often. UHGG genomes of the spiked species contain nearly all of X's
  peptides, observed or not, so carriage measures how close the two species are in the catalogue,
  not which peptides the sample selected. Only *K. variicola* shows an excess over the control. The
  data cannot tell "the spiked strain's variants" from "a genome nearly identical to the spiked
  species"; for the report both mean the same thing: X cannot be told apart from the spiked
  species at species level. That is what `same_genus_as` labels and what step 4 would resolve to
  genus: with the control carried at 85-99%, X has almost no catalogue-unique peptides at species
  level.
- Family level: carriage is low (0-12%), so the spiked strains in UHGG do not explain these
  genomes' kept peptides (*S. marcescens_I* has a small excess). Their evidence (4-68 kept peptides
  against `lambda` of 0.4-3.3) remains unexplained; it is not random matching (B1 calibrates) and
  not a spiked strain's variants. Candidates: carry-over or contamination, or organisms outside UHGG.
  They are reported as genus errors in E1 and left for step 4 and E2.

## 3. Proposed method

Steps 1-3 decide which genomes are accepted (selection, razor assignment). Step 5 is protein
inference within the accepted genomes. Step 4 decides the rank each accepted genome is reported at;
it never changes acceptance, and it is a separate, later change (see Phasing).

What the current code does, and what happens to it:

| Current | Plan |
|---|---|
| Per-unit weight = min(50, logit) x intensity factor; top-5 sum `S` | drop |
| Resampling null, M-draw top-5, B = 10,000 (floor 1/10,001) | replace by step 2 |
| Global decoy pool, no per-genome share | replace by step 2 |
| Protein p-value by resampling | replace by step 5 |
| Base-peptide unit, I = L | keep |
| Decoy-genome calibration in the run-info | keep, as rows of the genome table |
| Razor restricted to accepted genomes | keep |

### Step 1. Peptide-level FDR per class

Unit: base sequence with I = L; score = best PSM score for that unit. Separate target-decoy walks
per class: unmodified (exact mass), known shift (accepted localized modification), unexplained
shift. Keep q <= 0.01 in each class. Unexplained-shift peptides never count as genome or protein
evidence. Deamidation is a known shift, as everywhere else in MetaPilot.

The genome service reads `MetaPSM.getScore()`, whose meaning is set by each engine:

| Engine | `getScore()` | Decoys at all scores for step 2 |
|---|---|---|
| MetaPilot open search (`OpenSearchMetaPsm`) | `class_logit`; falls back to class probability, then hyperscore, when unset | yes (`scored_psms.tsv`) |
| pFind (`PFindPSM`) | `Final_Score` (direction to verify; the current code treats higher as better) | to verify (`REV_` prefix) |
| FragPipe (`FragPipePSM`) | MSFragger hyperscore | to verify: `psm.tsv` is filtered |
| Sage (`SagePSM`) | hyperscore | to verify |

The first build applies steps 1-3 and 5 to the open search only, and requires `class_logit` on every
PSM (no mixed-scale fallback inside one walk). Its input rows are the spectrum winners of
`scored_psms.tsv` (`is_spectrum_winner`), one per spectrum. In MSFragger-only runs these are the
`hit_rank` 1 rows the prototype used (1,253,437 on 12mix); with MSFragger + Sage, `hit_rank` 1
appears once per engine, and the winner keeps the spectrum counted once. Other engines keep the current genome test until
their score direction and all-score decoys are confirmed.

### Step 2. Expected false count per genome

All counts are **distinct peptides pooled over the runs analysed together**, on both sides.
(Summing per-run decoy counts against distinct target peptides inflates the null; counting per
(peptide, run) on both sides counts the same false peptide once per run, and false peptides recur,
section 2.5.)

**Uniqueness across runs.** Uniqueness is decided per run against that run's database: a peptide is
unique to g in run r when every protein of r's database that contains it (search digest rules:
strict trypsin, up to 2 missed cleavages, length 7-50, protein N-terminal Met clipped, I = L) belongs
to g. Pooled over runs, a peptide counts for g only when every run in which it was accepted calls it
unique to g; a peptide unique to different genomes in different runs, or shared in any run, is
dropped. The decoy side uses the identical rule, with the reversed genomes.

**Decoy uniqueness from the database, not the PSM protein column.** The run's database is the FASTA
stage 2 searched (`marker/<run>.fasta`: the full proteomes of the run's selected genomes plus one
`rev_` decoy per target). Target and decoy uniqueness both come from digesting it; the protein column
of `scored_psms.tsv` names one protein per hit and is never used for uniqueness (section 2.5, caveat).

**Classes: a deliberate difference.** `n_g` and `F` count the evidence classes (unmodified and known
shift). The `pi_g` window takes rank-1 decoy peptides of **any class**, unexplained shift included.
This is chosen on measurement (section 2.6, B1): 12mix has only 5,734 distinct unique rank-1 decoy
peptides in the evidence classes, fewer than `30 * F` = 24,060, so an evidence-class window is the
all-score window of those classes, with its bias (chi-square p = 4e-24; reversed *P. aeruginosa*
accepted at q = 7e-9). The any-class top `30 * F` window gives p = 0.009 and no reversed genome
accepted. The window models where false matches land by score, and the high-scoring decoys of any
class serve that; it is not a count of evidence. E1 reproduces the any-class result.

- `n_g` = accepted target peptides (step 1, classes unmodified and known shift) unique to genome g.
- `F` = accepted decoy peptides (same classes) unique to one reversed genome.
- `pi_g` = share of the top `K * F` distinct rank-1 decoy peptides by score (any class, any
  q-value) that are unique to `rev_g`; each decoy peptide counted once, at its best score. `K = 30`:
  proposed on PXD050781 (section 2.5) and kept on 12mix (B1, section 2.7: any-class top `30 * F`,
  chi-square p = 0.009, no reversed genome accepted); the per-run chi-square watches it (section 6,
  open item 3). All-score decoys are not used: their share is biased
  by up to 4x by how well a genome is represented. A genome absent from a run's database gets no
  decoys from that run, so per-run databases enter through `pi_g`; `runs_searched` is reported.
  The window count gets one added: `pi_g = (w_g + 1) / (W + G)`, with `w_g` the window decoys unique
  to `rev_g`, `W` the window size and `G` the candidate genomes. Without it a genome with no window
  decoy (typically one searched in one or two runs) had `lambda` floored at 1e-9 and was accepted on
  one peptide (E3: five genomes on PXD050781 MSFragger-only, six with Sage). After exclusion
  `w_g` counts the surviving window decoys. The chi-square keeps the raw shares.
- `lambda_g = F * pi_g`.
- `p_g = P(Poisson(lambda_g) >= n_g)`, computed as a log survival probability (`log_p`), because
  most real genomes underflow to p = 0. Benjamini-Hochberg over genomes.
- Every reversed genome goes through the same calculation (`n` = its accepted decoy peptides,
  `lambda = F * pi_g`) and is reported as a calibration row. Expected: decoy genomes accepted at
  q <= 0.01 close to 0, and the chi-square of accepted decoys against `F * pi_g` not significant;
  the run-info reports both.

### Step 3. Sequential near-neighbour exclusion

Accept genomes in ascending `log_p`, ties broken by `(n - lambda) / sqrt(lambda)`. Before testing
genome g, remove from `n_g` each of its unique peptides that is identical to, or within one
substitution of, a tryptic peptide of the **full proteome of an already accepted genome**, whether
or not that genome was in the database of the run where the peptide was found. Per-run databases
create exactly this case: a relative selected in one run absorbs spectra in another.

"One substitution" includes pairs that are equal after deamidation (N/D, Q/E, 0.984 Da) or within
the search's fragment tolerance in mass (e.g. K/Q 0.036 Da; I/L is already merged). Apply the same
exclusion to the decoy side against the reversed proteomes of the accepted genomes, then recompute
`lambda_g` and `p_g`. A genome that fails after exclusion is not accepted, and later genomes are not
excluded against it. When two related genomes are both present, the first accepted takes their
shared neighbour peptides; this is conservative, and `excluded_near_neighbour` shows its cost.

The full proteome of a candidate genome is in the database of every run that selected it (stage 2
writes each selected genome's whole `.faa`), so the union of the runs' databases holds the full
proteome of every genome that can be accepted, and of its reversed counterpart.

The order uses the running Benjamini-Hochberg threshold to decide which genomes exclude later ones.
The reported q-values are Benjamini-Hochberg over all genomes' p-values after exclusion, and
`accepted` means q <= 0.01. The two can differ at the margin; the run-info counts the genomes where
they do.

This step only sees relatives among the candidate genomes; relatives outside them are step 4's job.

### Step 4. Catalogue resolution of accepted genomes (second build)

For each accepted genome X, in acceptance order, using the peptides kept after step 3:

1. **Catalogue lookup.** One streaming pass of all kept peptides (all genomes together) over
   UHGP-100 with `PangenomeVariantLookup.scan` (tryptic context on, species `BitSet` per peptide),
   plus the species of multi-species clusters from `uhgp-100.rep_species.tsv.gz`.
2. **Clade-restricted counts.** For each rank R in {species, genus, family}: `n_R(X)` = kept peptides
   found in the catalogue whose carrying species all lie in X's clade at R (GTDB lineage from
   `genomes-all_metadata.tsv`). Peptides not found in the catalogue are not counted (risk 3).
3. **Conditioning on earlier genomes.** At every rank, drop peptides whose catalogue species set
   includes the species of a genome resolved to species earlier in the order. Step 3 only removes
   neighbours of the accepted genome itself; a peptide absent from MGYG000000035 but present in
   another *B. hydrogenotrophica* genome in UHGG would otherwise count as genus- and family-unique
   for a *Blautia_A* relative, which would then resolve on the target's evidence. The condition uses
   species resolved earlier, not every accepted genome: a genome accepted but resolved only at genus
   or family has not established its species, and dropping by it would let an uncertain call erase
   another genome's evidence. Peptides shared with such a genome's clade still resolve to that
   shared clade, which rolls up into the same taxon.
4. **Null.** `lambda_R(X) = lambda_X * u_R(X)`, with `lambda_X` from step 2 after exclusion and
   `u_R(X)` the fraction of X's own tryptic peptides (search digest rules) that are R-unique in the
   catalogue. `u_R = 1` is the conservative fallback until the background scan exists (risk 1).
   `p_R(X) = P(Poisson(lambda_R) >= n_R)`, Benjamini-Hochberg over accepted genomes per rank.
5. **Resolved rank.** The lowest R with q <= 0.01; `none` when no rank passes. Taxonomy outputs
   report X at that rank: genomes resolved at the same genus or family roll up into one taxon.

Until step 4 is built, and whenever the prepared catalogue is absent, `resolved_rank` is
`not_assessed`, and the run-info says that accepted genomes are accepted against the candidate
genomes only, not resolved against relatives outside them.

### Step 5. Picked protein-group FDR within accepted genomes

Razor peptides go to proteins of accepted genomes only (unchanged). The decoy protein groups are
built by the same grouping and razor procedure applied to the reversed proteins of the accepted
genomes and the accepted decoy peptides; otherwise every decoy peptide counts as unique and the
picked comparison is biased. Protein groups are scored by their best peptide; each target group
competes with its reversed counterpart (picked FDR: Savitski et al. 2015, *Mol Cell Proteomics*
14:2394; The et al. 2016, Percolator 3.0, *J Am Soc Mass Spectrom* 27:1719), and groups are kept at
q <= 0.01.

MetaPilot's DDA reports have no protein groups: each peptide has one razor protein, assigned by
intensity. Decoy peptides from `scored_psms.tsv` are not quantified, so the intensity rule cannot be
applied to them. Step 5 therefore uses one rule on both sides: each step 1 evidence peptide goes to
the protein of an accepted genome (reversed protein on the decoy side) that holds the most evidence
peptides, ties broken by protein ID; a protein is scored by its best peptide; target `P` competes
with `rev_P`.

**The reported razor proteins follow the same rule (open search).** In the first build the protein
tables kept the intensity razor, so `protein_q_value` sat next to a protein whose razor peptides
differ from the ones step 5 scored: on 12mix MSFragger-only, 7,833 of 20,130 razor proteins (39%)
received no step 5 peptide and have empty score cells, holding 18,450 of 106,881 razor peptides
(5,062 of them single-peptide proteins). The open search therefore assigns each reported peptide by
the count rule (most evidence peptides among accepted genomes' proteins, ties by protein ID), and
intensities are summed over those razor peptides. Every reported protein is then scored on its own
razor peptides. Under the closed-then-open plan every reported peptide is in an evidence class, so
the rule covers all of them; until then, a protein whose razor peptides are all unexplained shifts
has empty score cells. pFind, FragPipe and Sage keep the intensity razor.

The protein q-value is not a false-discovery rate where genome acceptance is wrong (E2: all 1,541
non-spiked proteins at q <= 0.01 on 12mix sit on accepted non-spiked genomes, against about 113
false proteins the q-value expects). The tables report it (`protein_q_value`) and nothing is
filtered on it.

### Interface with razor assignment and the output files

- **Razor assignment.** `DdaMagProteinInferenceService.confidentGenomes` reads
  `genomeStats[IDX_Q_VALUE] <= 0.01`. The new test fills the same map, one `double[3]` per target
  genome: `[0]` kept unique peptides, `[1]` p-value, `[2]` Benjamini-Hochberg q-value; the genome
  score map is `1 - q` as now. `confidentGenomes` is unchanged, and razor peptides go to genomes
  with q <= 0.01, the same set as `accepted`.
- **Who writes which file.** The open search writes `genome_scores.tsv` and
  `genome_scores_runinfo.txt`, and no longer writes `genome_rank_scores.tsv` or its run-info. pFind,
  FragPipe and Sage keep the current test and keep writing `genome_rank_scores.tsv`; one run never
  writes both.
- **Readers.** `Genome.tsv` keeps the column `genome_q_value`; for the open search the other two
  stats columns are named `kept_unique_peptides` and `genome_p_value` instead of `genome_score_S`
  and `genome_p_emp`. `Taxa.tsv` and the parquet sum `[0]` over a taxon's genomes; for the open
  search that is the sum of kept unique peptides, and the column is named `kept_unique_peptides`
  (parquet `evidence_score` keeps its name, its meaning documented in the parquet schema); the other
  engines keep `Score`. The only reader of
  the genome files outside the exporter (`Claude outputs/.../01_genomes.py`) reads `genome_q_value`,
  which keeps its meaning. The GUI task writes its own rank file and is not on this path.
- **Protein tables.** For the open search the three score columns of the protein tables are
  `protein_score` (best peptide score), `protein_decoy_score` (score of `rev_P`, empty when it has
  no evidence peptide) and `protein_q_value`.

### Output

`genome_scores.tsv`, replacing `genome_rank_scores.tsv`, one row per genome and per reversed genome:

`genome`, `is_decoy`, `runs_searched`, `unique_peptides`, `excluded_near_neighbour`,
`kept_peptides`, `decoy_unique_peptides`, `expected_false`, `log_p`, `p_value`, `q_value`,
`accepted`, `acceptance_order`, `resolved_rank`, `same_genus_as` (accepted target genomes only;
empty for reversed genomes).

The second build adds `catalogue_species_unique`, `catalogue_genus_unique`,
`catalogue_family_unique`, `q_species`, `q_genus`, `q_family`, `resolved_taxon`.

Run-info: class-wise peptide FDR counts, `F`, `K` and the score at the window edge, decoy genomes
accepted, the chi-square of accepted decoys against `F * pi_g`; in the second build, catalogue path
and version, peptides queried and found, scan time, the `u_R` used.

### Phasing

- **First build: steps 1, 2, 3 and 5** (open search). They replace the current genome score test and
  protein resampling; razor peptides are assigned as now. `resolved_rank` = `not_assessed`.
  Validated on 12mix first (section 4).
- **Second build: step 4.** Needs the one-time catalogue preparation (about 43 minutes) and a 41 GB
  scan per analysis, and is the least validated part. Fills the catalogue columns; taxonomy outputs
  report the resolved rank only after it is validated. Whether genus-only genomes keep their razor
  peptides is decided then.

### Code touch points (expected)

First build:

- `DdaMagGenomeScoringService`: steps 1-3 replace the score and resampling code; keep `peptideKey`
  and `uniqueDecoyGenome`. Needs the rank-1 decoy peptides with scores (not only the validated
  ones) for `pi_g`; the open search has them in `scored_psms.tsv`, so the post-analysis input must
  carry them.
- Near-neighbour exclusion: one digest of the union of the runs' databases (target and `rev_`
  proteins, step 2 digest rules) gives each peptide's genome. The peptides to test (evidence
  peptides, `pi_g` window decoys) are indexed by segments (pigeonhole: `d` substitutions leave one of
  `d + 1` segments exact), and every digested peptide of equal length is checked against them by
  Hamming distance. `TrypticPeptideIndex` (a target-only hash membership index) does not fit.
- `DdaMagProteinInferenceService`: picked protein-group FDR replaces the protein resampling test;
  decoy groups by the same procedure.
- `MagTaxonomyTsvExporter`: `genome_scores.tsv` and its run-info.
- `MagPostAnalysisScoringService` / `DdaMagPostAnalysisService`: route the open search to the new
  test, other engines to the current one.
- `DDAMagOpenSearchParameters` / `DDAMagOpenSearchParameterIO`: `generate_sample_specific_db` is
  ignored by the open-search runner; remove it or document that per-run databases are always built.
- Unit tests: synthetic genomes for steps 2-3 (decoy genomes not accepted; a one-substitution
  relative excluded after its neighbour is accepted, including a relative absent from that run's
  database); picked protein groups with shared decoy peptides.

After the first build (version 8 review):

- `DdaMagProteinInferenceService`: for the open search, razor peptides assigned by the step 5 count
  rule (shared with `GenomeEvidenceTest`'s step 5, so the two cannot drift), protein intensities
  summed over them. Test: every reported protein with an evidence razor peptide has a
  `protein_q_value`, and its scored peptides are its razor peptides.
- `MagTaxonomyTsvExporter` (`Taxa.tsv`): column `kept_unique_peptides` instead of `Score` for the
  open search; the parquet schema documents `evidence_score` for both meanings.
- `GenomeEvidenceTest` calibration rows: each reversed genome's `w_g` and `W` also computed with
  its own evidence decoys left out, reported next to the standard rows (E3).

Second build:

- New class in `core/mag` (e.g. `GenomeCatalogueResolution`) for step 4, with a public entry point
  into the package-private `core/variant` lookup (`PangenomeVariantLookup`, `FastaScanner`,
  `SpeciesMap`).
- `PangenomeResource.locate` and the `uhgp100_path` / `uhgp100_membership_path` settings, shared
  with the pangenome variant check; a one-time `PangenomeResourcePreparer` run on the raw archive
  (`Z:\database\human-gut\v2.0.2\uhgp-100.tar.gz`).
- `MagOpenSearchHeadlessRunner` (message when the catalogue is absent) and a setting such as
  `genome_catalogue_resolution` = `auto` | `on` | `off`.

## 4. Validation plan

First build, on 12mix (known composition, 12 spiked genomes):

1. False genome acceptance among the non-spiked candidate genomes. Prototype: 15 of 20 (2.6);
   first build (E1): 14 of 20 (8 species errors, 4 genus errors, the *Bacillus_A* pair).
2. Sensitivity on the 12 spiked genomes, and `excluded_near_neighbour` for spiked genomes with a
   spiked relative. Prototype: 12 of 12; the *Klebsiella* and *E. coli* genomes lose 13-42%.
3. Reversed-genome calibration rows: decoy genomes accepted at q <= 0.01, and the chi-square of
   accepted decoys against `F * pi_g`; choose `K` here (section 2.5). Prototype: p = 0.019 (10F),
   0.009 (30F), 1e-27 (all scores); `F` = 802 limits the power.
4. False protein rate: accepted proteins from absent genomes, against the picked protein-group q.
   Not prototyped.

Secondary checks:

5. PXD050781: reproduce sections 2.4 and 2.5 in Java; the target accepted first.
6. MSFragger + Sage open search on the same data.

Second build:

7. Reproduce section 2.3 in Java, then with the step 4 conditioning; background scan for `u_R`;
   12mix resolved ranks against the known species; runtime on the HPC with the prepared resource.

## 5. Risks and limitations

1. **Conservative catalogue null.** `u_R = 1` treats every false match as R-unique; it under-calls
   resolution. The background scan replaces it.
2. **Membership is an upper bound.** UHGP-100 clusters with a coverage threshold, so a member can
   lack a peptide its representative contains (`UhgpMembership` Javadoc). Species sets are upper
   bounds and catalogue uniqueness is a lower bound.
3. **Catalogue completeness.** If the organism present is not in UHGG, a peptide it shares only
   with X looks X-unique and X is over-called. Peptides absent from the catalogue (true novel
   variants of a present strain) give no resolution evidence.
4. **Depth, not only identity.** A low-abundance organism is detected mostly through conserved,
   abundant proteins, so it often resolves only to genus or family even when present. The report
   must say "resolved to", not "absent at species level".
5. **Step 2 models random matches only.** Decoys cannot generate a real peptide from a relative;
   step 3 (and later step 4) carries that error class. On 12mix step 3 does not carry it far enough:
   14 of 20 non-spiked genomes are accepted in the first build (E1; prototype 15 of 20, section 2.6).
6. **`pi_g` is window-dependent.** False matches do not land in proportion to the decoy share at
   all scores (section 2.5); the top-`K * F` window fits PXD050781 but `K` is fitted. The per-run
   chi-square makes a bad window visible.
7. **Poisson.** `lambda` is an estimate (its own sampling error is ignored), and false peptides
   recur across runs; distinct-peptide counting removes the recurrence from the counts but the
   accepted false matches are not independent draws.
8. **Per-run databases.** They are the design and the method accounts for them through `pi_g` and
   the full-proteome exclusion in step 3. `generate_sample_specific_db` is ignored by the runner
   (code touch points).

## 6. Decisions from the second and third audits, and what remains open

Decided:

- Q1 `pi_g`: tested (section 2.5). All-score decoys are rejected; the top-`K * F` window with
  `K = 30` fits, distinct peptides as the unit; `K` is set on 12mix.
- Q2 acceptance order: ascending `log_p`; the larger of two present relatives takes their shared
  neighbour peptides; `excluded_near_neighbour` reports the cost.
- Q3: Benjamini-Hochberg per rank, with the step 4 conditioning.
- Q4: razor peptides unchanged in the first build; decided with step 4.
- Q5: the method takes whatever run set it is given; per sample for a cohort, pooled as a summary.

Open:

1. pFind, FragPipe and Sage: score direction and availability of decoys at all scores (step 1
   table).
2. Step 4 conditioning on species resolved earlier versus on every accepted genome (step 4, item 3).
3. Whether `K` transfers between datasets, or needs to be a q-value threshold on the decoy walk
   (the audit's q <= 0.1 admits at least 10F decoys, since targets at q <= 0.1 are at least the
   ~100F accepted at q <= 0.01; the top 10F still disagrees here, chi-square p = 3e-7). On 12mix,
   10F and 30F are both borderline (p = 0.019 and 0.009, `F` = 802).

Reopened by 12mix (section 2.6) and decided by the third audit:

- **Q4, scope: no catalogue conditioning in acceptance.** Conditioning on genomes accepted earlier
  cannot fix the clearest case: *Enterococcus_D* sp002850555 is accepted 11th, the spiked
  *E. casseliflavus* 29th, so conditioning would erase the spiked genome's evidence, not the
  relative's (the same holds for the step 4, item 3 rule, which also depends on order). Acceptance
  stays steps 1-3 and answers only "evidence beyond random matches". The catalogue carriage was
  measured (B3, section 2.8): the same-genus relatives are nearly identical to the spiked species in
  UHGG (control carriage 85-99% where it is high), which step 4's reported rank handles; the
  family-level relatives are not explained by the spiked strains.
- **Q5, a relative holding a present strain's variants.** Choosing a species among close relatives
  is a limit of the database's resolution, not a genome q-value decision: the database genome
  closest to the real strain wins. The report marks it instead: a `same_genus_as` column naming the
  accepted genome with the most kept peptides in the same GTDB genus, when that is another genome.
  It uses the GTDB lineage already in the MAG database, no catalogue and no threshold. On the 12mix
  prototype (B4, section 2.7) it marks 7 of the 8 same-genus relatives (not sp002850555, which
  outranks its spiked neighbour), *B. luti*, and 2 spiked genomes that have a stronger same-genus
  neighbour: *K. pneumoniae* (538 kept, under *K. aerogenes*, 806, also spiked) and
  *E. casseliflavus* (under sp002850555). It is a label, not a rejection.
- **Q6, *Bacillus_A paranthracis*.** Reported separately, not as a false acceptance (485 kept
  against `lambda` = 55, no spiked relative in its family; *B. cereus*-group organisms are common lab
  contaminants). *Bacillus_A luti* counts as its relative, not as an independent error.
- **Open item 3, `K`.** Not decided now. `K = 30` for the first build, with the per-run chi-square
  and the count of accepted reversed genomes in the run-info, so a poorly fitting window shows on
  any dataset.
- **Step 3 neighbourhood: one substitution.** Measured (B2, section 2.7): two substitutions reject
  one genus error and lose the spiked *E. casseliflavus*.

Still open: item 1 (other engines), and the design of step 4's conditioning given the order
problem above.

## 7. Next work plan

Order: A and B are measurements and report fixes on existing 12mix outputs (no Java, nothing written
to `Z:`); C updates this document with their results; D is the first build; E is its validation on
the server. Prototype code stays in the agent store (`mix12_eval`).

Progress at version 4: A1, B1-B4 and C done (sections 2.7 and 2.8; step 2 class choice and step 3
neighbourhood fixed in section 3). D built: `GenomeEvidenceTest` (steps 1, 2, 3, 5),
`OpenSearchGenomeEvidence` (reads `scored_psms.tsv` and the stage-2 databases), routed through
`DdaMagPostAnalysisService` and the open-search exporters; unit tests in `GenomeEvidenceTestTest`
and `OpenSearchGenomeEvidenceTest`. E1 and E3 done on the 2026-10-12 server runs (results under E);
the add-one window count of step 2 came out of E3. A2 is the E1 table.

### A. Report fixes (this document)

- A1. Section 2.6: mechanism stated as a hypothesis; genus-level tally; *Bacillus_A* pair separate.
  Done in this revision.
- A2. Validation table (section 4) scores species errors and genus errors separately, with the
  *Bacillus_A* pair reported on its own line. After B1-B2.

### B. Measurements on 12mix (MSFragger-only open search, 2026-10-11 outputs)

- **B1. Reversed-genome acceptance.** Run every reversed genome through steps 2-3 with the same
  `lambda` and Benjamini-Hochberg thresholds; report the number accepted at q <= 0.01 for
  `K` = 10, 30 and all scores, next to the chi-square. Uses the evidence classes for the `pi_g`
  window (step 2, one class definition), with the any-class window of section 2.6 alongside.
  Extends `40_genome_test_prototype.py`; minutes.
- **B2. Two-substitution neighbourhood in step 3.** Same prototype, step 3 with Hamming distance
  <= 2 at equal length, on both the target and the decoy side. Index: the digested peptides split
  into three segments (pigeonhole: two mismatches leave one segment exact), candidates verified by
  Hamming distance; about 20 M index entries, within memory. Report, against the one-substitution
  run: non-spiked genomes rejected (species errors, genus errors, *Bacillus_A* pair), spiked genomes
  still accepted, and each spiked genome's `excluded_near_neighbour` (the *Klebsiella* and *E. coli*
  genomes lose 13-42% with one substitution). Rule: adopt two substitutions if it rejects relatives,
  every spiked genome stays accepted (watch *K. aerogenes* next to *K. pneumoniae*), and B1's
  reversed-genome acceptance is no worse than with one substitution; otherwise keep one substitution
  and record the result.
- **B3. Catalogue carriage of the relatives' peptides.** UHGP-100 scan (as `34_catalogue_uniqueness.py`,
  about 12 minutes) of the kept peptides of the 15 accepted non-spiked genomes. For each relative X
  and the spiked strain it sits next to (by genus, else family): the fraction of X's kept peptides
  carried in UHGG by the spiked species (other genomes of *E. casseliflavus*, *K. pneumoniae*,
  *E. coli_D* and so on). Control: the same fraction for a random sample of X's tryptic peptides
  that were not observed, so that high carriage is not just similarity between the two species. The
  mechanism holds for X when observed carriage clearly exceeds the control. Reported separately for
  the 8 same-genus and 5 family-level relatives.
- **B4. `same_genus_as` preview.** Computed in the prototype output, to confirm the labels listed
  under Q5 (section 6) before the column is built.

### C. Document update

- C1. Section 2.6 and section 4 with B1-B4: reversed genomes accepted, one vs two substitutions,
  carriage against control, the species/genus split.
- C2. Step 3 text fixed to the neighbourhood B2 chooses; step 4 text records the B3 result and the
  order problem of conditioning.

### D. First build (Java, open search only)

Steps 1, 2, 3 and 5 as in section 3, with:

- D1. `K = 30`; per-run chi-square and the number of reversed genomes accepted at q <= 0.01 in the
  run-info; reversed genomes as calibration rows of `genome_scores.tsv`.
- D2. Step 3 neighbourhood from B2.
- D3. `same_genus_as` column in `genome_scores.tsv` (GTDB genus from the MAG database metadata).
- D4. Step 5 picked protein-group FDR with decoy groups built by the same procedure; reported,
  not used as a filter, until E2.
- D7. Razor interface and output files as in section 3, "Interface". D starts once B2 has chosen
  the neighbourhood.
- D5. Unit tests as listed under the code touch points, plus `same_genus_as` and reversed-genome
  rows.
- D6. Other engines keep the current genome test (open item 1).

### E. Validation of the first build (server runs by the user)

- E1. 12mix MSFragger-only: the numbers of B1-B2 reproduced by the Java output (species errors,
  genus errors, *Bacillus_A* pair, spiked genomes accepted, reversed genomes accepted).
- E2. 12mix false protein rate: accepted protein groups whose proteins all belong to non-spiked
  genomes, at species level and at genus level, against the picked protein-group q. The build is
  not called validated before this is measured.
- E3. PXD050781 and the MSFragger + Sage open search as secondary checks.

Results of the 2026-10-12 runs (before the add-one window count):

- E1, 12mix MSFragger-only: 12 of 12 spiked accepted, 0 reversed genomes, chi-square p = 0.013;
  26 accepted = 12 spiked, 8 species errors, 4 genus errors, 2 *Bacillus_A*. The prototype's fifth
  genus error, *Kluyvera* sp902363335, is rejected (q = 0.027). `same_genus_as` labels as under Q5.
  F = 751 against 802 in the prototype, and unique counts are 1-3% lower: the regenerated
  `scored_psms.tsv` puts 50,034 rank-1 PSMs on accepted peaks against 96,668 in the 2026-10-11
  output, so fewer known-shift peptides enter; the genome test itself gives the same counts on the
  same input. The add-one count changes no 12mix decision.
  Cause (`logging.log` of both runs): peak discovery is identical (7,420 modes, the same 1,279
  non-zero peaks), but the 2026-10-11 run still had the old peak filter ("161/1231 peaks kept at
  1.50% target FDR", 151 accepted delta-mass groups) and the 2026-10-12 run the peak screen of the
  2026-10-11 afternoon edits to `ExpectedPeakSelector`, `ModificationValidator` and
  `DeltaMassEvidenceGate` ("8/1232 peaks accepted", the zero family and 7 allow-list masses; 1,224
  unexplained masses identification-only; 6 accepted groups). About 46,600 rank-1 PSMs on
  non-allow-list peaks moved from known shift to unexplained, which is not evidence; validated PSMs
  barely changed (858,392 against 858,696). An intended change, not a reproducibility fault. The
  prototype numbers of sections 2.6 and 2.7 were computed on the old evidence set and compare with
  E1 only approximately.
- E3, PXD050781: target genome accepted first; 0 reversed genomes in both workflows; MSFragger-only
  35 of 39 accepted (chi-square p = 0.90), MSFragger + Sage 31 of 39 (p = 0.034). Genomes with no
  window decoy were accepted on 1-3 peptides at `lambda` = 1e-9, and the reversed-genome rows cannot
  show it (a reversed genome's own evidence decoys sit in the window). Recomputed from the tables
  with the add-one count: 32 and 29 accepted, the one-peptide genomes dropped (the prototype with
  the add-one count also accepts 32 on MSFragger-only). Not yet rerun.
- Reversed genomes with their own evidence decoys left out of the window
  (`mix12_eval/56_reversed_leave_own_out.py`, MSFragger-only, add-one count for the targets):
  `w_g` and `W` both drop the genome's own evidence decoys, with and without the add-one count.
  Reversed genomes at q <= 0.01: 0 of 32 (12mix) and 0 of 39 (PXD050781) in every variant, and none
  at unadjusted p <= 0.01; summed `lambda` against observed decoys 513 / 498 (12mix) and 188 / 156
  (PXD050781), slightly conservative. The add-one count itself is not tested by this: on PXD050781
  only 9 reversed genomes have `w_g` <= 2 after leaving out their own decoys, holding 1 evidence
  decoy between them against 0.43 expected with the add-one count and 0.13 without, and none on
  12mix. A reversed genome's evidence decoys are its highest-scoring decoys, so one with no window
  decoys rarely has evidence decoys; the regime where the add-one count matters for targets
  (E3 above) has almost no decoy counterpart. The build reports the left-out rows; they check
  calibration on average, not the pseudocount.
- E2, 12mix (`mix12_eval/57_e2_nonspiked_proteins.py`): MSFragger-only accepts 11,346 target
  proteins at q <= 0.01, of which 1,541 (13.6%) are on non-spiked genomes: same genus 1,139, same
  family 197, other 205 (*Bacillus_A*). All 1,541 are on accepted genomes (razor peptides go only
  there), 14 genomes, led by *Klebsiella* MGYG000002511 (537) and *Enterococcus_D* (239); 1,373
  have at least one razor peptide found in no spiked proteome (strain differences and relatives'
  own sequence, or false matches), 168 only peptides that a spiked proteome also holds. The
  q-value expects about 113 false proteins. MSFragger + Sage: 1,504 of 11,082, same pattern;
  FragPipe closed with the resampling genome score: 3,038 of 13,148. Verdict: the non-spiked
  protein share is set by genome acceptance, not by random protein matches, so the protein q is not
  a false-discovery rate wherever a genome is accepted wrongly. It stays reported and unfiltered
  (D4); the build is validated for random matches only, and genome errors are measured by E1.
- 12mix MSFragger + Sage (rerun 2026-10-12 19:45): 12 of 12 spiked accepted, 0 reversed genomes,
  chi-square p = 0.078; 27 accepted = the MSFragger-only set plus *Enterobacter mori* (species error,
  4 kept peptides, `same_genus_as` *E. asburiae_B*). No genome at the 1e-9 floor.
- The comparison with the closed search, and the identification design that follows from it:
  `closed-then-open-identification-plan.md`. E1 and E2 are run again after that build.

Step 4 (catalogue resolution) stays a separate, later change.

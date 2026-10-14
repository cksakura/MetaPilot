package uk.ac.quadram.metapilot.core.variant;

import java.io.IOException;
import java.io.Writer;

/**
 * {@code variant_columns.tsv}: one row per column of the pangenome variant
 * check's tables, saying what it holds. Kept beside the column lists it
 * describes; {@code VariantColumnDictionaryTest} fails when the two drift.
 */
final class VariantColumnDictionary {

	private static final String C = PangenomeVariantCheck.CANDIDATES_FILE;
	private static final String S = PangenomeVariantCheck.SITES_FILE;
	private static final String P = PangenomeVariantCheck.SUMMARY_FILE;
	private static final String Q = PangenomeVariantCheck.SPECIES_FILE;

	static final String[][] ENTRIES = {
			{ "*", "*", "NA: not available or not applicable. A candidate is one PSM x one substitution site. "
					+ "Searched databases: the union of every sample's database in the run, target entries only. "
					+ "Nothing here claims a strain: carrier counts say which catalogued genomes are compatible "
					+ "with a variant, not which genome is in the sample." },

			{ C, "engine", "pfind or open_search" },
			{ C, "file", "Spectrum file of the PSM, as the engine reported it" },
			{ C, "scan", "Scan number of the PSM" },
			{ C, "hit_rank", "Open search: the engine's hit rank; with file and scan, the key back to "
					+ "validated_psms.tsv. NA for pFind" },
			{ C, "charge", "Precursor charge; with file and scan, the key back to pFind's psms.tsv" },
			{ C, "precursor_mz", "Observed precursor m/z" },
			{ C, "psm_q_value", "The PSM's q-value as the workflow reported it: pFind's q-value, or the open-search "
					+ "validator's q-value within the PSM's modification group" },
			{ C, "base_peptide", "Peptide as identified (the database sequence), I and L as reported" },
			{ C, "position", "1-based index of the substituted residue in base_peptide" },
			{ C, "from_aa", "Residue in base_peptide" },
			{ C, "to_aa", "Residue in the variant" },
			{ C, "variant_peptide", "base_peptide with every substitution of the PSM applied (all sites of a pFind "
					+ "multi-substitution PSM, so it can differ from base_peptide at more than this position)" },
			{ C, "n_substitutions_in_variant_peptide", "Residues where variant_peptide differs from base_peptide, I "
					+ "and L counted as one" },
			{ C, "n_candidates_for_psm", "Candidates from the same PSM, this one included: the sites of one pFind "
					+ "multi-substitution peptide (all apply together), or the open-search reassigner's alternative "
					+ "readings of one delta mass (only one can be right)" },
			{ C, "delta_mass", "mass(to_aa) - mass(from_aa), Da, with the engine's cysteine convention" },
			{ C, "assigned_proteins", "Searched-database proteins the PSM was assigned to, as the engine reported "
					+ "them; the assigned_* and in_same_species columns refer to these" },
			{ C, "protein_sites", "protein:position of the substituted residue in each assigned protein whose "
					+ "sequence contains base_peptide in the searched databases (first occurrence); NA when none" },
			{ C, "residue_before", "Residue before base_peptide in the first located protein; '-' at the protein "
					+ "start" },
			{ C, "residue_after", "Residue after base_peptide in the first located protein; '-' at the protein end" },
			{ C, "peak_id", "Open search: the delta-mass peak of the PSM. NA for pFind" },
			{ C, "peak_substitution_shaped", "Open search: whether the substitution reassigner judged the peak "
					+ "substitution-shaped (false: a peak whose centroid matches a substitution mass but that failed "
					+ "its shape or database test; nothing in the workflow acts on these). NA for pFind" },
			{ C, "ms2_decidable", "false when a modification of the same residue has the same mass (within the "
					+ "tolerance), so MS2 cannot tell the two readings apart" },
			{ C, "isobaric_modifications", "Those same-residue modifications, and pFind's own alias if it gave one" },
			{ C, "in_searched_db", "variant_peptide is in a searched database - any sample's, not necessarily the "
					+ "one this spectrum was searched against. Where it was, a direct match was possible and the "
					+ "substitution is not needed to explain the spectrum" },
			{ C, "in_pangenome_all", "variant_peptide occurs in UHGP-100 (any species)" },
			{ C, "in_pangenome_detected_species", "variant_peptide occurs in UHGP-100 in a species of "
					+ "genome_final_selected.tsv" },
			{ C, "in_same_species", "variant_peptide occurs in UHGP-100 in a species of an assigned protein" },
			{ C, "n_species", "UHGP-100 species whose proteins contain variant_peptide" },
			{ C, "assigned_species", "Species representatives of the assigned proteins" },
			{ C, "example_protein", "A UHGP-100 protein containing variant_peptide (the first in the catalogue)" },
			{ C, "n_genomes", "Genomes, any species, with a protein in a UHGP-100 cluster containing "
					+ "variant_peptide; an upper bound, since cluster members are identical to the representative over "
					+ "the clustered length, not necessarily over the peptide" },
			{ C, "n_genomes_same_species", "Of n_genomes, those in assigned_species" },
			{ C, "n_genomes_in_species", "Every catalogued genome of assigned_species, carrying the variant or not" },
			{ C, "carrier_fraction_in_species", "n_genomes_same_species / n_genomes_in_species" },
			{ C, "n_isolates_same_species", "Of n_genomes_same_species, isolate genomes rather than MAGs" },
			{ C, "example_genome_same_species", "A same-species genome carrying the variant: the lowest-numbered "
					+ "isolate, else the lowest-numbered genome" },
			{ C, "ms2_verdict", "Open search MS2 form ranking among the forms at this site: SUBSTITUTION, "
					+ "OTHER_SUBSTITUTION, MODIFICATION, AMBIGUOUS, or NOT_RANKED; NOT_APPLICABLE for pFind, whose "
					+ "calls are taken as made. Secondary evidence: it does not set pangenome_class" },
			{ C, "ms2_score_substitution", "Spectral-entropy similarity of this candidate's variant" },
			{ C, "ms2_score_other_substitution", "Best other same-mass substitution of the same residue" },
			{ C, "ms2_score_modification", "The delta mass as a modification of the same residue" },
			{ C, "ms2_score_elsewhere", "Best form changing any other residue" },
			{ C, "ms2_localises_elsewhere", "A form at another residue beats every same-site form by more than the "
					+ "margin" },
			{ C, "ms2_best_form", "Best-scoring form at any site" },
			{ C, "pangenome_class", "UNDECIDABLE_BY_MS2 (ms2_decidable false), else SUBSTITUTION_IN_SEARCHED_DB "
					+ "(in_searched_db), else SUBSTITUTION_IN_PANGENOME (in_pangenome_all), else SUBSTITUTION_NOVEL" },
			{ C, "assigned_kegg_ko", "eggNOG KEGG KOs of the assigned proteins, union over them" },
			{ C, "assigned_cog", "eggNOG COGs of the assigned proteins, union over them" },
			{ C, "assigned_cog_category", "eggNOG COG categories of the assigned proteins, union over them" },
			{ C, "assigned_preferred_name", "eggNOG preferred names of the assigned proteins, union over them" },
			{ C, "assigned_description", "eggNOG descriptions of the assigned proteins, union over them" },
			{ C, "front_end_call", "What the engine called: pFind's modification entry as written; open search, "
					+ "the reassigner's verdict, or peak_not_substitution_shaped" },
			{ C, "localized_sequence", "Open search: the engine's localized sequence. NA for pFind" },
			{ C, "reassigner_in_db", "Open search: whether the reassigner found the variant in its own database. "
					+ "NA for pFind" },

			{ S, "engine", "pfind or open_search" },
			{ S, "protein", "An assigned protein; a candidate counts under each of its assigned proteins; NA for "
					+ "candidates without one" },
			{ S, "protein_species", "Species representative of protein" },
			{ S, "protein_position", "1-based position of the substituted residue in protein; NA when base_peptide "
					+ "was not found in protein (the row is then one peptide site)" },
			{ S, "from_aa", "Residue in the protein" },
			{ S, "to_aa", "Residue in the variant" },
			{ S, "substitution", "from->to, I and L shown as I/L" },
			{ S, "located", "protein_position is known" },
			{ S, "peptide_sites", "base_peptide@position of every candidate at this site" },
			{ S, "variant_peptides", "Distinct variant peptides at this site" },
			{ S, "n_candidates", "Candidates at this site" },
			{ S, "n_candidates_shared", "Of n_candidates, those assigned to more than one protein (counted under "
					+ "each)" },
			{ S, "n_psms", "Distinct PSMs" },
			{ S, "n_files", "Distinct spectrum files" },
			{ S, "files", "The spectrum files" },
			{ S, "min_psm_q_value", "Lowest psm_q_value of the site's candidates" },
			{ S, "n_candidates_undecidable_by_ms2", "Candidates of class UNDECIDABLE_BY_MS2" },
			{ S, "n_candidates_in_searched_db", "Candidates of class SUBSTITUTION_IN_SEARCHED_DB" },
			{ S, "n_candidates_in_pangenome", "Candidates of class SUBSTITUTION_IN_PANGENOME" },
			{ S, "n_candidates_novel", "Candidates of class SUBSTITUTION_NOVEL" },
			{ S, "in_searched_db", "Any candidate's variant_peptide is in a searched database" },
			{ S, "in_pangenome_all", "Any candidate's variant_peptide is in UHGP-100" },
			{ S, "in_pangenome_detected_species", "Any candidate's variant_peptide is in UHGP-100 in a selected "
					+ "species" },
			{ S, "in_same_species", "Any candidate's variant_peptide is in UHGP-100 in its assigned species" },
			{ S, "n_genomes_same_species", "n_genomes_same_species of the candidate with the most same-species "
					+ "carriers" },
			{ S, "n_genomes_in_species", "n_genomes_in_species of that candidate" },
			{ S, "carrier_fraction_in_species", "carrier_fraction_in_species of that candidate" },
			{ S, "n_ms2_ranked", "Open search: candidates the MS2 form ranking gave a verdict. NA for pFind" },
			{ S, "n_ms2_substitution", "Open search: of n_ms2_ranked, verdict SUBSTITUTION. NA for pFind" },

			{ P, "engine", "pfind or open_search" },
			{ P, "level", "substitution: per substitution type (open search: and per peak_substitution_shaped, "
					+ "pooled over peaks); substitution_peak: open search, per delta-mass peak; all: per "
					+ "ms2_decidable; base: the identified base peptides, a lookup sanity check that should be close to "
					+ "100%" },
			{ P, "substitution", "from->to, I and L shown as I/L; ALL or base for those levels" },
			{ P, "ms2_decidable", "As in variant_candidates.tsv" },
			{ P, "peak_id", "Delta-mass peak (level substitution_peak); ALL when pooled over peaks" },
			{ P, "peak_substitution_shaped", "As in variant_candidates.tsv" },
			{ P, "n_candidates", "Candidates (PSM x site)" },
			{ P, "n_unique_variant_peptides", "Distinct variant peptides; every percentage is over unique "
					+ "peptides, so one abundant peptide counts once" },
			{ P, "n_unique_in_pangenome_all", "Of them, in UHGP-100" },
			{ P, "pct_in_pangenome_all", "100 * n_unique_in_pangenome_all / n_unique_variant_peptides" },
			{ P, "n_unique_in_pangenome_detected", "Of them, in UHGP-100 in a selected species" },
			{ P, "pct_in_pangenome_detected", "100 * n_unique_in_pangenome_detected / n_unique_variant_peptides" },
			{ P, "control_n", "Distinct matched-control peptides: each candidate's substitution applied to a random "
					+ "occurrence of its from-residue in a random identified peptide of the run, with no spectrum "
					+ "behind it - what chance alone gives for that substitution" },
			{ P, "control_n_in_pangenome_all", "Of them, in UHGP-100" },
			{ P, "control_pct_in_pangenome_all", "100 * control_n_in_pangenome_all / control_n" },
			{ P, "control_n_in_pangenome_detected", "Of them, in UHGP-100 in a selected species" },
			{ P, "control_pct_in_pangenome_detected", "100 * control_n_in_pangenome_detected / control_n" },
			{ P, "enrichment_all", "pct_in_pangenome_all / control_pct_in_pangenome_all" },
			{ P, "enrichment_detected", "pct_in_pangenome_detected / control_pct_in_pangenome_detected" },
			{ P, "fisher_p_all", "One-sided Fisher exact p (variants found more often than controls), unique "
					+ "peptides treated as independent" },
			{ P, "fisher_p_detected", "The same for the selected-species scope" },
			{ P, "fisher_q_all", "Benjamini-Hochberg q of fisher_p_all over the rows of the same level" },
			{ P, "fisher_q_detected", "Benjamini-Hochberg q of fisher_p_detected over the rows of the same level" },

			{ Q, "engine", "pfind or open_search" },
			{ Q, "species", "Species representative genome" },
			{ Q, "in_genome_final_selected", "The species has a genome in genome_final_selected.tsv: the "
					+ "post-analysis selected it, which is not proof it is present" },
			{ Q, "n_genomes_in_species", "Catalogued genomes of the species" },
			{ Q, "n_isolates_in_species", "Of them, isolate genomes" },
			{ Q, "n_candidates", "Candidates (PSM x site) with an assigned protein of this species" },
			{ Q, "n_candidates_shared", "Of them, also assigned to another species (counted under each)" },
			{ Q, "n_candidates_undecidable_by_ms2", "Candidates of class UNDECIDABLE_BY_MS2" },
			{ Q, "n_candidates_in_searched_db", "Candidates of class SUBSTITUTION_IN_SEARCHED_DB" },
			{ Q, "n_candidates_in_pangenome", "Candidates of class SUBSTITUTION_IN_PANGENOME" },
			{ Q, "n_candidates_novel", "Candidates of class SUBSTITUTION_NOVEL" },
			{ Q, "n_peptides", "Distinct variant peptides of n_candidates" },
			{ Q, "n_peptides_undecidable_by_ms2", "Distinct variant peptides of class UNDECIDABLE_BY_MS2" },
			{ Q, "n_peptides_in_searched_db", "Distinct variant peptides of class SUBSTITUTION_IN_SEARCHED_DB" },
			{ Q, "n_peptides_in_pangenome", "Distinct variant peptides of class SUBSTITUTION_IN_PANGENOME" },
			{ Q, "n_peptides_in_pangenome_this_species", "Of n_peptides_in_pangenome, found in this species" },
			{ Q, "n_peptides_in_pangenome_other_species_only", "Of n_peptides_in_pangenome, found only in other "
					+ "species" },
			{ Q, "n_peptides_novel", "Distinct variant peptides of class SUBSTITUTION_NOVEL" },
			{ Q, "n_other_genomes_carrying", "Genomes of this species, other than its representative, carrying any "
					+ "of n_peptides_in_pangenome_this_species; an upper bound (see n_genomes)" },
			{ Q, "n_other_isolates_carrying", "Of them, isolate genomes" } };

	private VariantColumnDictionary() {
	}

	static void write(Writer w) throws IOException {
		w.write("file\tcolumn\tmeaning\n");
		for (String[] e : ENTRIES) {
			w.write(e[0]);
			w.write('\t');
			w.write(e[1]);
			w.write('\t');
			w.write(PangenomeVariantCheck.clean(e[2]));
			w.write('\n');
		}
	}
}

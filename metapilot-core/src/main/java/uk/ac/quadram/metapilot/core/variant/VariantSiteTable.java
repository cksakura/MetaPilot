package uk.ac.quadram.metapilot.core.variant;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * {@code variant_sites.tsv}: the candidates rolled up to one row per protein
 * site and substitution, {@code (protein, protein_position, from_aa, to_aa)}.
 * <p>
 * A candidate counts under every protein it was assigned to
 * ({@code n_candidates_shared} says how many of them were assigned to more
 * than one). The site is where the substituted residue sits in the protein,
 * from the base peptide's first occurrence in that protein among the searched
 * databases. When the base peptide was not found in an assigned protein (an
 * accession the databases do not use, or no tryptic occurrence with tryptic
 * context required), the row is per peptide site instead:
 * {@code protein_position} is NA and {@code located} false. Candidates without
 * any assigned protein are grouped under protein NA.
 * <p>
 * Different base peptides (missed cleavages, other substitutions in the same
 * PSM) can give the same site different variant peptides, and so different
 * pangenome answers: the class columns count candidates, the
 * {@code in_*} columns are true when any candidate's variant peptide is, and
 * the genome counts are those of the candidate with the most same-species
 * carriers.
 */
final class VariantSiteTable {

	static final String[] COLUMNS = { "engine", "protein", "protein_species", "protein_position", "from_aa", "to_aa",
			"substitution", "located", "peptide_sites", "variant_peptides", "n_candidates", "n_candidates_shared",
			"n_psms", "n_files", "files", "min_psm_q_value", "n_candidates_undecidable_by_ms2",
			"n_candidates_in_searched_db", "n_candidates_in_pangenome", "n_candidates_novel", "in_searched_db",
			"in_pangenome_all", "in_pangenome_detected_species", "in_same_species", "n_genomes_same_species",
			"n_genomes_in_species", "carrier_fraction_in_species", "n_ms2_ranked", "n_ms2_substitution" };

	private static final String NA = "NA";

	static final class Row {
		final String protein;
		final int proteinPosition;
		final char fromAa;
		final char toAa;
		final String engine;
		final Set<String> peptideSites = new TreeSet<>();
		final Set<String> variantPeptides = new TreeSet<>();
		final Set<String> psms = new TreeSet<>();
		final Set<String> files = new TreeSet<>();
		int candidates;
		int shared;
		int undecidable;
		int inSearched;
		int inPangenome;
		int novel;
		Boolean anyInSearched;
		boolean anyInPangenome;
		Boolean anyInDetected;
		Boolean anyInSameSpecies;
		VariantCandidate mostCarried;
		double minQ = Double.NaN;
		int ranked;
		int ms2Substitution;

		Row(String engine, String protein, int proteinPosition, char fromAa, char toAa) {
			this.engine = engine;
			this.protein = protein;
			this.proteinPosition = proteinPosition;
			this.fromAa = fromAa;
			this.toAa = toAa;
		}
	}

	private VariantSiteTable() {
	}

	static List<Row> summarise(List<VariantCandidate> candidates) {
		Map<String, Row> rows = new HashMap<>();
		for (VariantCandidate c : candidates) {
			List<String> proteins = c.assignedProteins.isEmpty() ? List.of(NA) : new ArrayList<>(c.assignedProteins);
			for (String protein : proteins) {
				int at = NA.equals(protein) ? -1 : c.proteinPosition(protein);
				String key = at >= 0
						? String.join("\u0001", protein, String.valueOf(at), String.valueOf(c.fromAa),
								String.valueOf(c.toAa))
						: String.join("\u0001", protein, "peptide", c.basePeptide, String.valueOf(c.position),
								String.valueOf(c.toAa));
				Row r = rows.computeIfAbsent(key,
						k -> new Row(PangenomeVariantCheck.engineLabel(c.engine), protein, at, c.fromAa, c.toAa));
				add(r, c);
			}
		}
		List<Row> out = new ArrayList<>(rows.values());
		out.sort(Comparator.comparing((Row r) -> NA.equals(r.protein)).thenComparing(r -> r.protein)
				.thenComparing(r -> r.proteinPosition < 0).thenComparingInt(r -> r.proteinPosition)
				.thenComparing(r -> String.join(";", r.peptideSites)).thenComparing(r -> r.toAa));
		return out;
	}

	private static void add(Row r, VariantCandidate c) {
		r.candidates++;
		if (c.assignedProteins.size() > 1) {
			r.shared++;
		}
		r.peptideSites.add(c.basePeptide + "@" + c.position);
		r.variantPeptides.add(c.variantPeptide);
		r.psms.add(c.file + "\u0001" + c.scan + "\u0001" + c.psmId);
		r.files.add(c.file == null ? "" : c.file);
		if (!Double.isNaN(c.qValue) && (Double.isNaN(r.minQ) || c.qValue < r.minQ)) {
			r.minQ = c.qValue;
		}
		switch (c.pangenomeClass) {
		case UNDECIDABLE_BY_MS2:
			r.undecidable++;
			break;
		case SUBSTITUTION_IN_SEARCHED_DB:
			r.inSearched++;
			break;
		case SUBSTITUTION_IN_PANGENOME:
			r.inPangenome++;
			break;
		default:
			r.novel++;
		}
		r.anyInSearched = any(r.anyInSearched, c.inSearchedDb);
		r.anyInPangenome |= c.inPangenomeAll;
		r.anyInDetected = any(r.anyInDetected, c.inPangenomeDetected);
		r.anyInSameSpecies = any(r.anyInSameSpecies, c.inSameSpecies);
		if (c.nGenomesSameSpecies >= 0
				&& (r.mostCarried == null || c.nGenomesSameSpecies > r.mostCarried.nGenomesSameSpecies)) {
			r.mostCarried = c;
		}
		if (c.ms2Verdict != VariantCandidate.Ms2Verdict.NOT_APPLICABLE
				&& c.ms2Verdict != VariantCandidate.Ms2Verdict.NOT_RANKED) {
			r.ranked++;
			if (c.ms2Verdict == VariantCandidate.Ms2Verdict.SUBSTITUTION) {
				r.ms2Substitution++;
			}
		}
	}

	/** True when either is true; null only when both are unknown. */
	private static Boolean any(Boolean a, Boolean b) {
		if (Boolean.TRUE.equals(a) || Boolean.TRUE.equals(b)) {
			return true;
		}
		return a == null && b == null ? null : false;
	}

	static void write(Writer w, List<Row> rows, SpeciesMap species, boolean ranked) throws IOException {
		w.write(String.join("\t", COLUMNS));
		w.write('\n');
		for (Row r : rows) {
			String proteinSpecies = NA;
			if (species != null && !NA.equals(r.protein)) {
				int s = species.speciesOfAccession(r.protein);
				if (s >= 0) {
					proteinSpecies = species.speciesId(s);
				}
			}
			VariantCandidate m = r.mostCarried;
			String[] f = { r.engine, PangenomeVariantCheck.clean(r.protein), proteinSpecies,
					r.proteinPosition < 0 ? NA : String.valueOf(r.proteinPosition), String.valueOf(r.fromAa),
					String.valueOf(r.toAa), VariantCandidate.residueLabel(r.fromAa) + "->"
							+ VariantCandidate.residueLabel(r.toAa),
					String.valueOf(r.proteinPosition >= 0), String.join(";", r.peptideSites),
					String.join(";", r.variantPeptides), String.valueOf(r.candidates), String.valueOf(r.shared),
					String.valueOf(r.psms.size()), String.valueOf(r.files.size()),
					PangenomeVariantCheck.clean(String.join(";", r.files)),
					Double.isNaN(r.minQ) ? NA : String.format(Locale.ROOT, "%.3e", r.minQ),
					String.valueOf(r.undecidable), String.valueOf(r.inSearched), String.valueOf(r.inPangenome),
					String.valueOf(r.novel), PangenomeVariantCheck.bool(r.anyInSearched),
					String.valueOf(r.anyInPangenome), PangenomeVariantCheck.bool(r.anyInDetected),
					PangenomeVariantCheck.bool(r.anyInSameSpecies),
					m == null ? NA : String.valueOf(m.nGenomesSameSpecies),
					m == null || m.nGenomesInSpecies < 0 ? NA : String.valueOf(m.nGenomesInSpecies),
					m == null ? NA : PangenomeVariantCheck.carrierFraction(m), ranked ? String.valueOf(r.ranked) : NA,
					ranked ? String.valueOf(r.ms2Substitution) : NA };
			w.write(String.join("\t", f));
			w.write('\n');
		}
	}
}

package uk.ac.quadram.metapilot.core.variant;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.IntFunction;

/**
 * Pangenome variant check rolled up per species: for each species
 * representative the run selected or assigned a substitution PSM to, how its
 * candidates were classed, and how many of its other catalogued genomes carry
 * the variants the searched database could not explain.
 * <p>
 * A candidate counts under every species its assigned proteins come from
 * ({@code n_candidates_shared} says how many of them are shared).
 * {@code n_candidates_*} count candidates (PSM x site); {@code n_peptides_*}
 * count unique variant peptides. The
 * carrier columns are about the species' genomes other than its
 * representative, and are upper bounds: a UHGP-100 member protein is 100%
 * identical to its representative over the clustered length, not necessarily
 * over the peptide. They say which conspecific genomes are compatible with the
 * variants, not which strain is in the sample.
 */
final class VariantSpeciesSummary {

	static final String[] COLUMNS = { "engine", "species", "in_genome_final_selected", "n_genomes_in_species",
			"n_isolates_in_species", "n_candidates", "n_candidates_shared", "n_candidates_undecidable_by_ms2",
			"n_candidates_in_searched_db", "n_candidates_in_pangenome", "n_candidates_novel", "n_peptides",
			"n_peptides_undecidable_by_ms2", "n_peptides_in_searched_db", "n_peptides_in_pangenome",
			"n_peptides_in_pangenome_this_species", "n_peptides_in_pangenome_other_species_only",
			"n_peptides_novel", "n_other_genomes_carrying", "n_other_isolates_carrying" };

	private static final String NA = "NA";

	static final class Row {
		final int species;
		final String speciesId;
		Boolean selected;
		int candidates;
		int shared;
		int undecidable;
		int inSearched;
		int inPangenome;
		int novel;
		final Set<Integer> peptides = new HashSet<>();
		final Set<Integer> undecidablePeptides = new HashSet<>();
		final Set<Integer> searchedPeptides = new HashSet<>();
		final Set<Integer> novelPeptides = new HashSet<>();
		final Set<Integer> pangenomePeptides = new HashSet<>();
		final Set<Integer> thisSpecies = new HashSet<>();
		final Set<Integer> otherSpeciesOnly = new HashSet<>();
		/** -1 when there is no membership table. */
		int otherGenomesCarrying = -1;
		int otherIsolatesCarrying = -1;

		Row(int species, String speciesId) {
			this.species = species;
			this.speciesId = speciesId;
		}
	}

	private VariantSpeciesSummary() {
	}

	/**
	 * @param selected   species of {@code genome_final_selected.tsv}; null when unavailable
	 * @param carriersOf distinct carrier genome numbers of a query; null without a membership table
	 */
	static List<Row> summarise(List<VariantCandidate> candidates, SpeciesMap species, BitSet selected,
			PangenomeVariantCheck.Evidence ev, IntFunction<int[]> carriersOf) {
		Map<Integer, Row> rows = new HashMap<>();
		if (selected != null) {
			for (int s = selected.nextSetBit(0); s >= 0; s = selected.nextSetBit(s + 1)) {
				row(rows, species, s);
			}
		}
		for (VariantCandidate c : candidates) {
			BitSet own = PangenomeVariantCheck.ownSpecies(c, species);
			int q = ev.index(c.variantPeptide);
			PangenomeVariantLookup.QueryHit hit = ev.uhgp[q];
			for (int s = own.nextSetBit(0); s >= 0; s = own.nextSetBit(s + 1)) {
				Row r = row(rows, species, s);
				r.candidates++;
				if (own.cardinality() > 1) {
					r.shared++;
				}
				r.peptides.add(q);
				switch (c.pangenomeClass) {
				case UNDECIDABLE_BY_MS2:
					r.undecidable++;
					r.undecidablePeptides.add(q);
					break;
				case SUBSTITUTION_IN_SEARCHED_DB:
					r.inSearched++;
					r.searchedPeptides.add(q);
					break;
				case SUBSTITUTION_IN_PANGENOME:
					r.inPangenome++;
					r.pangenomePeptides.add(q);
					if (hit.species != null) {
						(hit.species.get(s) ? r.thisSpecies : r.otherSpeciesOnly).add(q);
					}
					break;
				default:
					r.novel++;
					r.novelPeptides.add(q);
				}
			}
		}
		for (Row r : rows.values()) {
			r.selected = selected == null ? null : selected.get(r.species);
			if (carriersOf != null) {
				countCarriers(r, species, carriersOf);
			}
		}
		List<Row> out = new ArrayList<>(rows.values());
		out.sort(Comparator.comparing((Row r) -> !Boolean.TRUE.equals(r.selected))
				.thenComparing(Comparator.comparingInt((Row r) -> r.candidates).reversed())
				.thenComparing(r -> r.speciesId));
		return out;
	}

	private static Row row(Map<Integer, Row> rows, SpeciesMap species, int s) {
		return rows.computeIfAbsent(s, k -> new Row(k, species.speciesId(k)));
	}

	private static void countCarriers(Row r, SpeciesMap species, IntFunction<int[]> carriersOf) {
		long representative = species.representativeNumber(r.species);
		BitSet genomes = new BitSet();
		for (int q : r.thisSpecies) {
			for (int g : carriersOf.apply(q)) {
				if (g != representative && species.speciesOfGenomeNumber(g) == r.species) {
					genomes.set(g);
				}
			}
		}
		r.otherGenomesCarrying = genomes.cardinality();
		if (species.hasGenomeType()) {
			int isolates = 0;
			for (int g = genomes.nextSetBit(0); g >= 0; g = genomes.nextSetBit(g + 1)) {
				if (species.isIsolate(g)) {
					isolates++;
				}
			}
			r.otherIsolatesCarrying = isolates;
		}
	}

	static void write(Writer w, String engine, List<Row> rows, SpeciesMap species) throws IOException {
		w.write(String.join("\t", COLUMNS));
		w.write('\n');
		int[] isolates = species == null || !species.hasGenomeType() ? null : species.isolatesPerSpecies();
		for (Row r : rows) {
			String[] f = { engine, r.speciesId, PangenomeVariantCheck.bool(r.selected),
					String.valueOf(species.genomesInSpecies(r.species)),
					isolates == null ? NA : String.valueOf(isolates[r.species]), String.valueOf(r.candidates),
					String.valueOf(r.shared), String.valueOf(r.undecidable), String.valueOf(r.inSearched),
					String.valueOf(r.inPangenome), String.valueOf(r.novel), String.valueOf(r.peptides.size()),
					String.valueOf(r.undecidablePeptides.size()), String.valueOf(r.searchedPeptides.size()),
					String.valueOf(r.pangenomePeptides.size()), String.valueOf(r.thisSpecies.size()),
					String.valueOf(r.otherSpeciesOnly.size()), String.valueOf(r.novelPeptides.size()),
					count(r.otherGenomesCarrying), count(r.otherIsolatesCarrying) };
			w.write(String.join("\t", f));
			w.write('\n');
		}
	}

	private static String count(int n) {
		return n < 0 ? NA : String.valueOf(n);
	}
}

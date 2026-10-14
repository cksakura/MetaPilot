package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.io.IOException;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

/**
 * Which query peptides occur in a set of protein FASTA files, and in which
 * species: one streaming pass, never an index of the database.
 * <p>
 * Queries must already be distinct and {@link PeptideAutomaton#normalise
 * normalised} (I and L collapsed). Used twice per run with the same code: once
 * over the databases the run searched (presence only) and once over UHGP-100
 * (presence, protein count, example protein and species).
 */
final class PangenomeVariantLookup {

	/** Where a query sits in one protein: its first occurrence that counts. */
	static final class Location {
		/** 1-based residue index of the query's first residue. */
		final int start;
		/** Residue before the query, {@code '-'} at the protein start. */
		final char before;
		/** Residue after the query, {@code '-'} at the protein end. */
		final char after;

		Location(int start, char before, char after) {
			this.start = start;
			this.before = before;
			this.after = after;
		}
	}

	/** What one query hit. */
	static final class QueryHit {
		int proteins;
		/** The hit protein that came first in the FASTA, so the choice is reproducible. */
		String exampleProtein = "";
		long exampleSerial = Long.MAX_VALUE;
		/** Species of every hit protein's cluster; null when no species map was given. */
		BitSet species;
		/** Hit proteins whose species could not be resolved. */
		int unresolvedProteins;
		/**
		 * {@link SpeciesMap#proteinKey} of every hit protein, for the queries the
		 * scan was asked to record; unsorted, length {@link #proteinKeyCount}.
		 */
		long[] proteinKeys;
		int proteinKeyCount;
		/** Locations in the proteins the scan was asked to locate this query in; null when none was found. */
		Map<String, Location> locations;

		boolean found() {
			return proteins > 0;
		}

		private void addProteinKey(long key) {
			if (proteinKeys == null) {
				proteinKeys = new long[4];
			} else if (proteinKeyCount == proteinKeys.length) {
				proteinKeys = java.util.Arrays.copyOf(proteinKeys, proteinKeyCount * 2);
			}
			proteinKeys[proteinKeyCount++] = key;
		}
	}

	private PangenomeVariantLookup() {
	}

	/**
	 * @param species null for presence only
	 * @return one entry per query, in query order
	 */
	static QueryHit[] scan(List<String> queries, List<File> fastas, SpeciesMap species, boolean requireTrypticContext,
			int threads, String decoyPrefix, Consumer<String> progress, String label) throws IOException {
		return scan(queries, fastas, species, requireTrypticContext, threads, decoyPrefix, progress, label, null);
	}

	/**
	 * @param recordProteins query indices whose hit protein IDs are kept in
	 *                       {@link QueryHit#proteinKeys}; null for none
	 */
	static QueryHit[] scan(List<String> queries, List<File> fastas, SpeciesMap species, boolean requireTrypticContext,
			int threads, String decoyPrefix, Consumer<String> progress, String label, BitSet recordProteins)
			throws IOException {
		return scan(queries, fastas, species, requireTrypticContext, threads, decoyPrefix, progress, label,
				recordProteins, null);
	}

	/**
	 * @param locateIn per query index, the protein IDs whose
	 *                 {@link QueryHit#locations} are wanted; null for none
	 */
	static QueryHit[] scan(List<String> queries, List<File> fastas, SpeciesMap species, boolean requireTrypticContext,
			int threads, String decoyPrefix, Consumer<String> progress, String label, BitSet recordProteins,
			Map<Integer, Set<String>> locateIn) throws IOException {
		QueryHit[] hits = new QueryHit[queries.size()];
		for (int i = 0; i < hits.length; i++) {
			hits[i] = new QueryHit();
			if (species != null) {
				hits[i].species = new BitSet();
			}
		}
		if (queries.isEmpty() || fastas.isEmpty()) {
			return hits;
		}
		PeptideAutomaton automaton = new PeptideAutomaton(queries);
		FastaScanner scanner = new FastaScanner(automaton, queries.size(), requireTrypticContext, threads,
				decoyPrefix);
		long start = System.currentTimeMillis();
		long proteins = scanner.scan(fastas, (q, proteinId, serial, offset, before, after) -> {
			QueryHit h = hits[q];
			Set<String> wanted = locateIn == null ? null : locateIn.get(q);
			synchronized (h) {
				h.proteins++;
				if (serial < h.exampleSerial) {
					h.exampleSerial = serial;
					h.exampleProtein = proteinId;
				}
				if (species != null && !species.addSpeciesOf(proteinId, h.species)) {
					h.unresolvedProteins++;
				}
				if (recordProteins != null && recordProteins.get(q)) {
					long key = SpeciesMap.proteinKey(proteinId);
					if (key >= 0) {
						h.addProteinKey(key);
					}
				}
				if (wanted != null && wanted.contains(proteinId)) {
					if (h.locations == null) {
						h.locations = new HashMap<>();
					}
					h.locations.putIfAbsent(proteinId, new Location(offset + 1, before, after));
				}
			}
		}, progress, label);
		if (progress != null) {
			int found = 0;
			for (QueryHit h : hits) {
				if (h.found()) {
					found++;
				}
			}
			progress.accept(String.format(java.util.Locale.ROOT,
					"Pangenome variant check: %s - %,d proteins scanned in %.0f s; %,d of %,d query peptides found "
							+ "(automaton: %,d states)",
					label, proteins, (System.currentTimeMillis() - start) / 1000.0, found, queries.size(),
					automaton.stateCount()));
		}
		return hits;
	}
}

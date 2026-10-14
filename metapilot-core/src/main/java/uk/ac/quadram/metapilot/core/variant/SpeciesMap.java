package uk.ac.quadram.metapilot.core.variant;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

/**
 * Protein -&gt; genome -&gt; species representative, for the UHGG catalogue.
 * <p>
 * UHGP-100's FASTA holds cluster representatives only, and a representative
 * stands for every member of its 100%-identity cluster - members that can come
 * from other genomes and other species. So the species behind a hit is the
 * union over the cluster, not the species of the representative's own genome.
 * The prepare step writes that union for the 6% of clusters that span more than
 * one species ({@code uhgp-100.rep_species.tsv.gz}); every other representative
 * maps to its own genome's species.
 * <p>
 * Genome -&gt; species comes from {@code genomes-all_metadata.tsv}
 * ({@code Genome} -&gt; {@code Species_rep}). Seventeen genomes are recorded there
 * with a version suffix ({@code MGYG000001406.1}) that their protein IDs drop,
 * so the unversioned ID is registered too.
 * <p>
 * Protein IDs are {@code MGYG<9 digits>_<n>}; they are held as one {@code long}
 * so the 5.7 million multi-species rows fit in about 100 MB.
 */
final class SpeciesMap {

	private final Map<String, Integer> speciesIndex = new HashMap<>();
	private final List<String> speciesIds = new ArrayList<>();
	/** Species index by MGYG genome number, -1 when unknown. */
	private int[] speciesByGenomeNumber = new int[0];
	/** Genomes whose ID is not {@code MGYG<digits>}. */
	private final Map<String, Integer> speciesByOtherGenome = new HashMap<>();
	/** Metadata rows per species index. */
	private int[] genomesPerSpecies = new int[0];
	/** MGYG genome numbers whose {@code Genome_type} is {@code Isolate}. */
	private final BitSet isolates = new BitSet();
	private boolean hasGenomeType;

	private long[] repKeys = new long[0];
	private int[] repOffsets = new int[] { 0 };
	private char[] repSpecies = new char[0];
	private int repCount;

	private SpeciesMap() {
	}

	/**
	 * @param metadataTsv  {@code genomes-all_metadata.tsv}
	 * @param repSpeciesGz {@code uhgp-100.rep_species.tsv.gz}, or null; without it
	 *                     each hit counts only its representative's own species
	 */
	static SpeciesMap load(File metadataTsv, File repSpeciesGz, Consumer<String> progress) throws IOException {
		SpeciesMap map = new SpeciesMap();
		map.loadMetadata(metadataTsv);
		if (repSpeciesGz != null && repSpeciesGz.isFile()) {
			long start = System.currentTimeMillis();
			map.loadRepSpecies(repSpeciesGz);
			if (progress != null) {
				progress.accept(String.format(java.util.Locale.ROOT,
						"Pangenome variant check: %,d multi-species clusters read from %s in %.1f s", map.repCount,
						repSpeciesGz.getName(), (System.currentTimeMillis() - start) / 1000.0));
			}
		}
		return map;
	}

	int speciesCount() {
		return speciesIds.size();
	}

	String speciesId(int index) {
		return speciesIds.get(index);
	}

	/** Species index of a genome, or -1. */
	int speciesOfGenome(String genome) {
		if (genome == null) {
			return -1;
		}
		long n = genomeNumber(genome, genome.length());
		if (n >= 0) {
			return n < speciesByGenomeNumber.length ? speciesByGenomeNumber[(int) n] : -1;
		}
		Integer s = speciesByOtherGenome.get(genome);
		return s == null ? -1 : s;
	}

	/** Species index of an MGYG genome number, or -1. */
	int speciesOfGenomeNumber(long genomeNumber) {
		return genomeNumber >= 0 && genomeNumber < speciesByGenomeNumber.length
				? speciesByGenomeNumber[(int) genomeNumber]
				: -1;
	}

	/** Genomes of a species in the metadata, all types. */
	int genomesInSpecies(int species) {
		return species >= 0 && species < genomesPerSpecies.length ? genomesPerSpecies[species] : 0;
	}

	/** False when the metadata has no {@code Genome_type} column. */
	boolean hasGenomeType() {
		return hasGenomeType;
	}

	boolean isIsolate(long genomeNumber) {
		return genomeNumber >= 0 && genomeNumber <= Integer.MAX_VALUE && isolates.get((int) genomeNumber);
	}

	/** Isolate genomes per species index, from the metadata. */
	int[] isolatesPerSpecies() {
		int[] out = new int[speciesIds.size()];
		for (int g = isolates.nextSetBit(0); g >= 0; g = isolates.nextSetBit(g + 1)) {
			int s = speciesOfGenomeNumber(g);
			if (s >= 0) {
				out[s]++;
			}
		}
		return out;
	}

	/** MGYG number of the species' representative genome, or -1 when its ID is not {@code MGYG<digits>}. */
	long representativeNumber(int species) {
		String id = speciesId(species);
		int dot = id.indexOf('.');
		return genomeNumber(id, dot > 0 ? dot : id.length());
	}

	/** {@code 69512} -&gt; {@code MGYG000069512}. */
	static String genomeIdOf(long genomeNumber) {
		return String.format(java.util.Locale.ROOT, "MGYG%09d", genomeNumber);
	}

	/** Genome number of a {@link #proteinKey}. */
	static long genomeNumberOfKey(long proteinKey) {
		return proteinKey / 8 / 1_000_000L;
	}

	/**
	 * Species of the genome a searched-database accession comes from
	 * ({@code MGYG000002463_01761} -&gt; species of {@code MGYG000002463}), or -1
	 * for a host, contaminant or otherwise unknown accession.
	 */
	int speciesOfAccession(String accession) {
		if (accession == null) {
			return -1;
		}
		String a = accession.trim();
		int space = a.indexOf(' ');
		if (space > 0) {
			a = a.substring(0, space);
		}
		int underscore = a.lastIndexOf('_');
		return speciesOfGenome(underscore > 0 ? a.substring(0, underscore) : a);
	}

	/**
	 * Adds every species the protein's cluster spans to {@code into}.
	 *
	 * @return false when the protein's species could not be resolved at all
	 */
	boolean addSpeciesOf(String proteinId, BitSet into) {
		long key = proteinKey(proteinId);
		if (key >= 0 && repCount > 0) {
			int i = Arrays.binarySearch(repKeys, 0, repCount, key);
			if (i >= 0) {
				for (int k = repOffsets[i]; k < repOffsets[i + 1]; k++) {
					into.set(repSpecies[k]);
				}
				return true;
			}
		}
		int underscore = proteinId.lastIndexOf('_');
		int s = speciesOfGenome(underscore > 0 ? proteinId.substring(0, underscore) : proteinId);
		if (s < 0) {
			return false;
		}
		into.set(s);
		return true;
	}

	private void loadMetadata(File metadataTsv) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(metadataTsv, StandardCharsets.UTF_8), 1 << 20)) {
			String header = br.readLine();
			if (header == null) {
				throw new IOException("Empty genome metadata: " + metadataTsv);
			}
			String[] cols = header.split("\t", -1);
			int gi = -1;
			int si = -1;
			int ti = -1;
			for (int i = 0; i < cols.length; i++) {
				if (cols[i].trim().equals("Genome")) {
					gi = i;
				} else if (cols[i].trim().equals("Species_rep")) {
					si = i;
				} else if (cols[i].trim().equals("Genome_type")) {
					ti = i;
				}
			}
			hasGenomeType = ti >= 0;
			if (gi < 0 || si < 0) {
				throw new IOException("Genome metadata has no Genome / Species_rep columns: " + metadataTsv);
			}
			String line;
			while ((line = br.readLine()) != null) {
				String[] f = line.split("\t", -1);
				if (f.length <= Math.max(gi, si) || f[gi].isEmpty() || f[si].isEmpty()) {
					continue;
				}
				int s = indexOfSpecies(f[si].trim());
				String genome = f[gi].trim();
				putGenome(genome, s);
				if (s >= genomesPerSpecies.length) {
					genomesPerSpecies = Arrays.copyOf(genomesPerSpecies, Math.max(s + 1, genomesPerSpecies.length * 2));
				}
				genomesPerSpecies[s]++;
				long n = genomeNumber(genome, genome.indexOf('.') > 0 ? genome.indexOf('.') : genome.length());
				if (ti >= 0 && ti < f.length && "Isolate".equalsIgnoreCase(f[ti].trim()) && n >= 0
						&& n <= Integer.MAX_VALUE) {
					isolates.set((int) n);
				}
				int dot = genome.indexOf('.');
				if (dot > 0) {
					String unversioned = genome.substring(0, dot);
					if (speciesOfGenome(unversioned) < 0) {
						putGenome(unversioned, s);
					}
				}
			}
		}
	}

	private void putGenome(String genome, int species) {
		long n = genomeNumber(genome, genome.length());
		if (n >= 0 && n < 50_000_000L) {
			if (n >= speciesByGenomeNumber.length) {
				int old = speciesByGenomeNumber.length;
				speciesByGenomeNumber = Arrays.copyOf(speciesByGenomeNumber, (int) Math.max(n + 1, old * 3L / 2));
				Arrays.fill(speciesByGenomeNumber, old, speciesByGenomeNumber.length, -1);
			}
			speciesByGenomeNumber[(int) n] = species;
		} else {
			speciesByOtherGenome.put(genome, species);
		}
	}

	private int indexOfSpecies(String species) {
		Integer s = speciesIndex.get(species);
		if (s == null) {
			if (speciesIds.size() >= Character.MAX_VALUE) {
				throw new IllegalStateException("More than 65535 species representatives");
			}
			s = speciesIds.size();
			speciesIndex.put(species, s);
			speciesIds.add(species);
		}
		return s;
	}

	private void loadRepSpecies(File gz) throws IOException {
		repKeys = new long[1 << 20];
		repOffsets = new int[(1 << 20) + 1];
		repSpecies = new char[1 << 22];
		repCount = 0;
		int used = 0;
		boolean sorted = true;
		try (InputStream in = new GZIPInputStream(new FileInputStream(gz), 1 << 20);
				BufferedReader br = new BufferedReader(new InputStreamReader(in, StandardCharsets.US_ASCII), 1 << 20)) {
			String line;
			while ((line = br.readLine()) != null) {
				int tab = line.indexOf('\t');
				if (tab <= 0) {
					continue;
				}
				long key = proteinKey(line.substring(0, tab));
				if (key < 0) {
					continue;
				}
				if (repCount == repKeys.length) {
					repKeys = Arrays.copyOf(repKeys, repKeys.length * 2);
					repOffsets = Arrays.copyOf(repOffsets, repKeys.length + 1);
				}
				if (repCount > 0 && key <= repKeys[repCount - 1]) {
					sorted = false;
				}
				repKeys[repCount] = key;
				repOffsets[repCount] = used;
				int from = tab + 1;
				while (from < line.length()) {
					int comma = line.indexOf(',', from);
					int to = comma < 0 ? line.length() : comma;
					if (to > from) {
						if (used == repSpecies.length) {
							repSpecies = Arrays.copyOf(repSpecies, repSpecies.length * 2);
						}
						repSpecies[used++] = (char) indexOfSpecies(line.substring(from, to).trim());
					}
					from = to + 1;
				}
				repCount++;
				repOffsets[repCount] = used;
			}
		}
		if (!sorted) {
			sortAndMerge();
		}
	}

	/** Sorts rows by key and unions rows that repeat a representative. */
	void sortAndMerge() {
		Integer[] order = new Integer[repCount];
		for (int i = 0; i < repCount; i++) {
			order[i] = i;
		}
		Arrays.sort(order, (a, b) -> Long.compare(repKeys[a], repKeys[b]));
		long[] keys = new long[repCount];
		int[] offsets = new int[repCount + 1];
		char[] species = new char[repOffsets[repCount]];
		int n = 0;
		int used = 0;
		BitSet merged = new BitSet();
		for (int i = 0; i < repCount;) {
			long key = repKeys[order[i]];
			merged.clear();
			while (i < repCount && repKeys[order[i]] == key) {
				for (int k = repOffsets[order[i]]; k < repOffsets[order[i] + 1]; k++) {
					merged.set(repSpecies[k]);
				}
				i++;
			}
			keys[n] = key;
			offsets[n] = used;
			for (int s = merged.nextSetBit(0); s >= 0; s = merged.nextSetBit(s + 1)) {
				species[used++] = (char) s;
			}
			n++;
			offsets[n] = used;
		}
		repKeys = keys;
		repOffsets = offsets;
		repSpecies = species;
		repCount = n;
	}

	/** A map with genome -&gt; species only, for the prepare step to add cluster rows to. */
	static SpeciesMap metadataOnly(File metadataTsv) throws IOException {
		SpeciesMap map = new SpeciesMap();
		map.loadMetadata(metadataTsv);
		map.repKeys = new long[1 << 20];
		map.repOffsets = new int[(1 << 20) + 1];
		map.repSpecies = new char[1 << 22];
		return map;
	}

	/** Appends one cluster row; call {@link #sortAndMerge()} before reading or writing. */
	void addRepRow(long key, BitSet species) {
		if (repCount == repKeys.length) {
			repKeys = Arrays.copyOf(repKeys, repKeys.length * 2);
			repOffsets = Arrays.copyOf(repOffsets, repKeys.length + 1);
		}
		int used = repOffsets[repCount];
		if (used + species.cardinality() > repSpecies.length) {
			repSpecies = Arrays.copyOf(repSpecies, Math.max(repSpecies.length * 2, used + species.cardinality()));
		}
		repKeys[repCount] = key;
		for (int s = species.nextSetBit(0); s >= 0; s = species.nextSetBit(s + 1)) {
			repSpecies[used++] = (char) s;
		}
		repCount++;
		repOffsets[repCount] = used;
	}

	int repRowCount() {
		return repCount;
	}

	/**
	 * Writes {@code representative \t species,species,...}, sorted by
	 * representative, species sorted.
	 */
	void writeRepSpecies(java.io.OutputStream out) throws IOException {
		java.io.Writer w = new java.io.BufferedWriter(new java.io.OutputStreamWriter(out, StandardCharsets.US_ASCII),
				1 << 20);
		List<String> names = new ArrayList<>();
		for (int i = 0; i < repCount; i++) {
			names.clear();
			for (int k = repOffsets[i]; k < repOffsets[i + 1]; k++) {
				names.add(speciesIds.get(repSpecies[k]));
			}
			java.util.Collections.sort(names);
			w.write(proteinIdOf(repKeys[i]));
			w.write('\t');
			w.write(String.join(",", names));
			w.write('\n');
		}
		w.flush();
	}

	/**
	 * {@code MGYG000004906_00001} -&gt; a {@code long} that round-trips through
	 * {@link #proteinIdOf(long)}: genome number, protein number and the protein
	 * number's digit count. -1 for any other shape, which callers resolve by
	 * genome instead.
	 */
	static long proteinKey(String id) {
		int underscore = id.lastIndexOf('_');
		int width = id.length() - underscore - 1;
		if (underscore != 13 || width < 1 || width > 6) {
			return -1;
		}
		long genome = genomeNumber(id, underscore);
		if (genome < 0) {
			return -1;
		}
		long protein = 0;
		for (int i = underscore + 1; i < id.length(); i++) {
			char c = id.charAt(i);
			if (c < '0' || c > '9') {
				return -1;
			}
			protein = protein * 10 + (c - '0');
		}
		return (genome * 1_000_000L + protein) * 8 + width;
	}

	static String proteinIdOf(long key) {
		int width = (int) (key % 8);
		long base = key / 8;
		return String.format(java.util.Locale.ROOT, "MGYG%09d_%0" + width + "d", base / 1_000_000L,
				base % 1_000_000L);
	}

	/** Numeric part of {@code MGYG<digits>} in {@code s[0, end)}, or -1. */
	private static long genomeNumber(String s, int end) {
		if (end <= 4 || end > 4 + 12 || !s.startsWith("MGYG")) {
			return -1;
		}
		long n = 0;
		for (int i = 4; i < end; i++) {
			char c = s.charAt(i);
			if (c < '0' || c > '9') {
				return -1;
			}
			n = n * 10 + (c - '0');
		}
		return n;
	}
}

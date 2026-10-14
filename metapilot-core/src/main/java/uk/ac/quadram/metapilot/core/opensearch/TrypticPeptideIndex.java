package uk.ac.quadram.metapilot.core.opensearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * Membership test for "is this peptide in the database we searched?".
 *
 * <h2>Why a hash index rather than a set of strings</h2>
 * The refined databases for one 12mix run digest to 6,532,401 distinct tryptic
 * peptides. Held as a {@code HashSet<String>} that is roughly 780 MB - it would
 * dominate the validation step's memory for a question answered a few thousand
 * times. Stored as sorted 64-bit hashes it is 52 MB and answers by binary
 * search.
 * <p>
 * The cost is a false-positive rate: two different peptides can hash alike. At
 * 64 bits with 6.5M entries the expected number of colliding pairs is about
 * 1e-6, so a false "yes" is far rarer than the biology this class is used to
 * reason about. A false positive would, at worst, propose a substitution
 * candidate that does not exist; the re-search that follows would simply fail to
 * identify it. Nothing downstream trusts membership alone.
 *
 * <h2>Digest</h2>
 * Strict trypsin (cut after K/R, no proline rule), the enzyme every search in
 * this workflow uses, with the same missed-cleavage and length bounds as
 * {@code fragger_open.params}. Those are parameters rather than constants
 * precisely so this index cannot silently describe a different peptide space
 * from the one MSFragger searched - if the search settings change, this must be
 * built with the new ones or its answers are about the wrong database.
 */
public final class TrypticPeptideIndex {

	/** Digest settings; must mirror the search that produced the pepXML. */
	public static final class Config {
		public int missedCleavages = 2;
		public int minLength = 7;
		public int maxLength = 50;
		/** Decoy prefix to skip; decoy sequences are not evidence for anything. */
		public String decoyTag = "rev_";
	}

	private final long[] hashes;
	private final int proteinCount;

	private TrypticPeptideIndex(long[] hashes, int proteinCount) {
		this.hashes = hashes;
		this.proteinCount = proteinCount;
	}

	public int size() {
		return hashes.length;
	}

	public int proteinCount() {
		return proteinCount;
	}

	/** True when the peptide occurs in one of the indexed target proteins. */
	public boolean contains(String peptide) {
		if (peptide == null || peptide.isEmpty()) {
			return false;
		}
		return Arrays.binarySearch(hashes, hash(peptide)) >= 0;
	}

	/**
	 * Builds the index from one or more FASTA files, target sequences only.
	 * <p>
	 * Duplicate protein entries across per-sample databases are common - the
	 * samples share most of their selected genomes - so sequences are digested
	 * once per distinct header.
	 */
	public static TrypticPeptideIndex build(List<File> fastas, Config cfg) throws IOException {
		if (cfg == null) {
			cfg = new Config();
		}
		LongBuffer collected = new LongBuffer(1 << 20);
		java.util.HashSet<String> seenProtein = new java.util.HashSet<>();
		int proteins = 0;
		for (File fasta : fastas) {
			if (fasta == null || !fasta.isFile()) {
				continue;
			}
			try (BufferedReader br = new BufferedReader(new FileReader(fasta))) {
				String line;
				StringBuilder seq = new StringBuilder();
				boolean keep = false;
				while ((line = br.readLine()) != null) {
					if (line.startsWith(">")) {
						if (keep && seq.length() > 0) {
							digest(seq.toString(), cfg, collected);
							proteins++;
						}
						seq.setLength(0);
						String name = line.substring(1).trim();
						int sp = name.indexOf(' ');
						if (sp > 0) {
							name = name.substring(0, sp);
						}
						keep = !name.toLowerCase(Locale.ROOT).startsWith(cfg.decoyTag.toLowerCase(Locale.ROOT))
								&& seenProtein.add(name);
					} else if (keep) {
						seq.append(line.trim());
					}
				}
				if (keep && seq.length() > 0) {
					digest(seq.toString(), cfg, collected);
					proteins++;
				}
			}
		}
		long[] arr = collected.toArray();
		Arrays.sort(arr);
		// Deduplicate in place; a peptide shared by many proteins is one entry.
		int n = 0;
		for (int i = 0; i < arr.length; i++) {
			if (i == 0 || arr[i] != arr[i - 1]) {
				arr[n++] = arr[i];
			}
		}
		return new TrypticPeptideIndex(Arrays.copyOf(arr, n), proteins);
	}

	private static void digest(String seq, Config cfg, LongBuffer out) {
		List<Integer> sites = new ArrayList<>();
		sites.add(0);
		for (int i = 0; i < seq.length() - 1; i++) {
			char c = seq.charAt(i);
			if (c == 'K' || c == 'R') {
				sites.add(i + 1);
			}
		}
		sites.add(seq.length());
		int n = sites.size() - 1;
		for (int i = 0; i < n; i++) {
			for (int mc = 0; mc <= cfg.missedCleavages; mc++) {
				int j = i + 1 + mc;
				if (j > n) {
					break;
				}
				int len = sites.get(j) - sites.get(i);
				if (len >= cfg.minLength && len <= cfg.maxLength) {
					out.add(hash(seq, sites.get(i), sites.get(j)));
				}
			}
		}
	}

	/**
	 * Growable {@code long[]}. A {@code List<Long>} boxes every entry, which on a
	 * 6.5M-peptide database is about 150 MB of short-lived objects for a value
	 * that is eight bytes.
	 */
	private static final class LongBuffer {
		private long[] data;
		private int size;

		LongBuffer(int capacity) {
			this.data = new long[Math.max(16, capacity)];
		}

		void add(long v) {
			if (size == data.length) {
				data = Arrays.copyOf(data, data.length + (data.length >> 1));
			}
			data[size++] = v;
		}

		long[] toArray() {
			return Arrays.copyOf(data, size);
		}
	}

	/** FNV-1a over the residue range; stable across runs and JVMs, unlike String.hashCode's 32 bits. */
	private static long hash(String s, int from, int to) {
		long h = 0xcbf29ce484222325L;
		for (int i = from; i < to; i++) {
			h ^= s.charAt(i);
			h *= 0x100000001b3L;
		}
		return h;
	}

	static long hash(String s) {
		return hash(s, 0, s.length());
	}
}

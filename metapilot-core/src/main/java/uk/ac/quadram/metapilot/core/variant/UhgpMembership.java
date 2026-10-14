package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Locale;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

/**
 * Member genomes of selected UHGP-100 clusters, read from MGnify's
 * {@code uhgp-100.tsv} ({@code representative \t member}, no header, about 630
 * million rows and 25 GB) or a gzipped copy of it.
 * <p>
 * One streaming pass that keeps only the requested clusters; nothing is
 * indexed or written. A member ID is {@code MGYG<genome>_<n>}, so the genome
 * comes from the ID itself.
 * <p>
 * MGnify clustered with a coverage threshold, so a member can be shorter than
 * its representative and may lack a peptide the representative contains.
 * Genome counts built from it are an upper bound.
 */
final class UhgpMembership {

	static final String FILE = "uhgp-100.tsv";
	static final String FILE_GZ = "uhgp-100.tsv.gz";

	/** Requested representatives' {@link SpeciesMap#proteinKey}s, sorted and distinct. */
	final long[] representatives;
	/** Member genome numbers per representative, sorted and distinct; empty when not in the file. */
	final int[][] genomes;
	final long rows;
	final int representativesFound;

	private UhgpMembership(long[] representatives, int[][] genomes, long rows, int found) {
		this.representatives = representatives;
		this.genomes = genomes;
		this.rows = rows;
		this.representativesFound = found;
	}

	/**
	 * The membership table: {@code configuredPath} (the file, or a directory
	 * holding it), else beside the prepared FASTA, either directly or in the
	 * {@code uhgp-100/} folder extracting the archive creates, there or one
	 * level up (the MAG version directory). Null when absent.
	 */
	static File locate(String configuredPath, File resourceDir) {
		if (configuredPath != null && !configuredPath.isBlank()) {
			File f = new File(configuredPath.trim());
			if (f.isFile()) {
				return f;
			}
			return f.isDirectory() ? inDirectory(f) : null;
		}
		if (resourceDir == null) {
			return null;
		}
		File parent = resourceDir.getAbsoluteFile().getParentFile();
		for (File dir : new File[] { resourceDir, new File(resourceDir, "uhgp-100"),
				parent == null ? null : new File(parent, "uhgp-100") }) {
			File f = dir == null ? null : inDirectory(dir);
			if (f != null) {
				return f;
			}
		}
		return null;
	}

	private static File inDirectory(File dir) {
		for (String name : new String[] { FILE, FILE_GZ }) {
			File f = new File(dir, name);
			if (f.isFile() && f.length() > 0) {
				return f;
			}
		}
		return null;
	}

	/** Member genomes of {@code representative}, or null when it was not requested. */
	int[] genomesOf(long representative) {
		int i = Arrays.binarySearch(representatives, representative);
		return i < 0 ? null : genomes[i];
	}

	/**
	 * @param wanted representative keys to keep; need not be sorted or distinct
	 */
	static UhgpMembership read(File tsv, long[] wanted, Consumer<String> progress) throws IOException {
		Rows parser = new Rows(Arrays.stream(wanted).sorted().distinct().toArray());
		if (parser.reps.length == 0) {
			return parser.finish();
		}
		long start = System.currentTimeMillis();
		long bytes = 0;
		long nextReport = 4L << 30;
		boolean gz = tsv.getName().toLowerCase(Locale.ROOT).endsWith(".gz");
		byte[] buf = new byte[1 << 22];
		byte[] line = new byte[128];
		int len = 0;
		boolean overflow = false;
		try (InputStream raw = new FileInputStream(tsv);
				InputStream in = gz ? new GZIPInputStream(raw, 1 << 20) : raw) {
			int n;
			while ((n = in.read(buf)) > 0) {
				if (Thread.currentThread().isInterrupted()) {
					throw new IOException("interrupted");
				}
				bytes += n;
				for (int p = 0; p < n; p++) {
					byte b = buf[p];
					if (b != '\n') {
						if (len < line.length) {
							line[len++] = b;
						} else {
							overflow = true;
						}
						continue;
					}
					if (!overflow) {
						parser.row(line, len);
					}
					len = 0;
					overflow = false;
				}
				if (progress != null && bytes >= nextReport) {
					nextReport += 4L << 30;
					progress.accept(String.format(Locale.ROOT,
							"Pangenome variant check: membership %,d GB read%s, %,d rows, %.0f s", bytes >> 30,
							gz ? " (uncompressed)" : " of " + (tsv.length() >> 30), parser.rows,
							(System.currentTimeMillis() - start) / 1000.0));
				}
			}
			if (!overflow) {
				parser.row(line, len);
			}
		}
		UhgpMembership m = parser.finish();
		if (progress != null) {
			progress.accept(String.format(Locale.ROOT,
					"Pangenome variant check: membership %s - %,d rows in %.0f s; %,d of %,d requested clusters found",
					tsv.getName(), m.rows, (System.currentTimeMillis() - start) / 1000.0, m.representativesFound,
					m.representatives.length));
		}
		return m;
	}

	/** Accumulates member genomes row by row. */
	private static final class Rows {
		final long[] reps;
		final int[][] lists;
		final int[] sizes;
		long rows;
		/** Rows of one cluster are contiguous in MGnify's file, so the last lookup is reused. */
		final byte[] previous = new byte[128];
		int previousLen = -1;
		int previousIndex = -1;

		Rows(long[] reps) {
			this.reps = reps;
			this.lists = new int[reps.length][];
			this.sizes = new int[reps.length];
		}

		void row(byte[] line, int len) {
			int end = len > 0 && line[len - 1] == '\r' ? len - 1 : len;
			if (end == 0) {
				return;
			}
			rows++;
			int tab = -1;
			for (int i = 0; i < end; i++) {
				if (line[i] == '\t') {
					tab = i;
					break;
				}
			}
			if (tab <= 0) {
				return;
			}
			int index;
			if (tab == previousLen && Arrays.equals(line, 0, tab, previous, 0, tab)) {
				index = previousIndex;
			} else {
				long key = proteinKey(line, 0, tab);
				index = key < 0 ? -1 : Arrays.binarySearch(reps, key);
				System.arraycopy(line, 0, previous, 0, tab);
				previousLen = tab;
				previousIndex = index;
			}
			if (index < 0) {
				return;
			}
			long genome = genomeNumber(line, tab + 1, end);
			if (genome < 0 || genome > Integer.MAX_VALUE) {
				return;
			}
			if (lists[index] == null) {
				lists[index] = new int[4];
			} else if (sizes[index] == lists[index].length) {
				lists[index] = Arrays.copyOf(lists[index], sizes[index] * 2);
			}
			lists[index][sizes[index]++] = (int) genome;
		}

		UhgpMembership finish() {
			int found = 0;
			int[][] genomes = new int[reps.length][];
			for (int i = 0; i < reps.length; i++) {
				if (lists[i] == null) {
					genomes[i] = new int[0];
					continue;
				}
				found++;
				genomes[i] = Arrays.stream(lists[i], 0, sizes[i]).sorted().distinct().toArray();
			}
			return new UhgpMembership(reps, genomes, rows, found);
		}
	}

	/** {@link SpeciesMap#proteinKey} of {@code b[from, to)}, without making a String. */
	static long proteinKey(byte[] b, int from, int to) {
		int underscore = -1;
		for (int i = to - 1; i >= from; i--) {
			if (b[i] == '_') {
				underscore = i;
				break;
			}
		}
		int width = to - underscore - 1;
		if (underscore - from != 13 || width < 1 || width > 6) {
			return -1;
		}
		long genome = genomeNumber(b, from, to);
		if (genome < 0) {
			return -1;
		}
		long protein = 0;
		for (int i = underscore + 1; i < to; i++) {
			int d = b[i] - '0';
			if (d < 0 || d > 9) {
				return -1;
			}
			protein = protein * 10 + d;
		}
		return (genome * 1_000_000L + protein) * 8 + width;
	}

	/** Genome number of {@code MGYG<digits>_...} in {@code b[from, to)}, or -1. */
	static long genomeNumber(byte[] b, int from, int to) {
		if (to - from < 6 || b[from] != 'M' || b[from + 1] != 'G' || b[from + 2] != 'Y' || b[from + 3] != 'G') {
			return -1;
		}
		long n = 0;
		int digits = 0;
		for (int i = from + 4; i < to && b[i] != '_'; i++) {
			int d = b[i] - '0';
			if (d < 0 || d > 9 || ++digits > 12) {
				return -1;
			}
			n = n * 10 + d;
		}
		return digits == 0 ? -1 : n;
	}
}

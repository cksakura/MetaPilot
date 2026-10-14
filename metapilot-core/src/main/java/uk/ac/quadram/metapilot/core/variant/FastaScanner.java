package uk.ac.quadram.metapilot.core.variant;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;

/**
 * Streams protein FASTA files once and reports which queries of a
 * {@link PeptideAutomaton} occur in which proteins.
 * <p>
 * One reader thread inflates and parses; {@code threads} matcher threads scan
 * batches of whole proteins. Multi-line records are concatenated before
 * matching, so a peptide split across a line break is found. Each
 * (query, protein) pair is reported at most once, however often the peptide
 * occurs in that protein.
 *
 * <h2>Tryptic context</h2>
 * With {@code requireTrypticContext} a match counts only where trypsin would
 * have produced it: the residue before it is K or R or the protein start - or
 * the protein starts with M and the match starts right after it (initiator
 * methionine removal) - and its last residue is K or R or the protein end. No
 * proline rule: every search this is compared against cuts before P.
 */
final class FastaScanner {

	/**
	 * Receives one (query, protein) pair, at the first occurrence that counts;
	 * must be thread-safe.
	 */
	interface HitSink {
		/**
		 * @param offset 0-based start of the match in the protein
		 * @param before residue before the match, {@code '-'} at the protein start
		 * @param after  residue after the match, {@code '-'} at the protein end
		 */
		void hit(int query, String proteinId, long proteinSerial, int offset, char before, char after);
	}

	private static final int READ_BUFFER = 1 << 20;
	private static final int BATCH_RESIDUES = 1 << 20;
	private static final int BATCH_PROTEINS = 8_192;

	private final PeptideAutomaton automaton;
	private final int queryCount;
	private final boolean requireTrypticContext;
	private final int threads;
	private final String decoyPrefix;

	FastaScanner(PeptideAutomaton automaton, int queryCount, boolean requireTrypticContext, int threads,
			String decoyPrefix) {
		this.automaton = automaton;
		this.queryCount = queryCount;
		this.requireTrypticContext = requireTrypticContext;
		this.threads = Math.max(1, threads);
		this.decoyPrefix = decoyPrefix == null ? "" : decoyPrefix.toLowerCase(Locale.ROOT);
	}

	/**
	 * Scans every file in order; {@code .gz} files are inflated on the fly.
	 * Protein serials run on across files, so they are unique for the whole call.
	 *
	 * @return number of proteins scanned (decoys excluded)
	 */
	long scan(List<File> fastas, HitSink sink, Consumer<String> progress, String label) throws IOException {
		BlockingQueue<Batch> queue = new ArrayBlockingQueue<>(threads * 2);
		ExecutorService pool = Executors.newFixedThreadPool(threads, r -> {
			Thread t = new Thread(r, "pangenome-matcher");
			t.setDaemon(true);
			return t;
		});
		// A matcher that stopped taking batches would leave the reader blocked
		// on the full queue for ever, so a failed matcher keeps draining until
		// the poison batch and the reader stops at its next batch.
		AtomicReference<Throwable> failure = new AtomicReference<>();
		List<Future<?>> workers = new ArrayList<>();
		for (int i = 0; i < threads; i++) {
			workers.add(pool.submit(() -> {
				long[] lastSeen = new long[queryCount];
				Arrays.fill(lastSeen, -1L);
				while (true) {
					Batch b = queue.take();
					if (b.count < 0) {
						return null;
					}
					if (failure.get() != null) {
						continue;
					}
					try {
						match(b, lastSeen, sink);
					} catch (Throwable t) {
						failure.compareAndSet(null, t);
					}
				}
			}));
		}
		long proteins = 0;
		long start = System.currentTimeMillis();
		try {
			Parser parser = new Parser(queue, failure, progress, label, start);
			for (File f : fastas) {
				try (InputStream in = open(f)) {
					parser.parse(in);
				}
			}
			parser.flush();
			proteins = parser.serial;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted while scanning " + label, e);
		} finally {
			for (int i = 0; i < threads; i++) {
				try {
					queue.put(Batch.POISON);
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
					break;
				}
			}
			pool.shutdown();
		}
		for (Future<?> w : workers) {
			try {
				w.get();
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				throw new IOException("Interrupted while scanning " + label, e);
			} catch (java.util.concurrent.ExecutionException e) {
				throw new IOException("Matcher thread failed while scanning " + label, e.getCause());
			}
		}
		if (failure.get() != null) {
			throw matcherFailed(label, failure.get());
		}
		return proteins;
	}

	private static IOException matcherFailed(String label, Throwable cause) {
		return new IOException("Matcher thread failed while scanning " + label + ": " + cause, cause);
	}

	private void match(Batch b, long[] lastSeen, HitSink sink) {
		byte[] seq = b.residues;
		for (int p = 0; p < b.count; p++) {
			final int ps = b.start[p];
			final int pe = b.end[p];
			final long serial = b.firstSerial + p;
			final String id = b.ids[p];
			automaton.scan(seq, ps, pe, (q, s, e) -> {
				if (lastSeen[q] == serial) {
					return;
				}
				if (requireTrypticContext && !trypticContext(seq, ps, pe, ps + s, ps + e)) {
					return;
				}
				lastSeen[q] = serial;
				int ms = ps + s;
				int me = ps + e;
				sink.hit(q, id, serial, s, ms == ps ? '-' : upper(seq[ms - 1]), me == pe - 1 ? '-' : upper(seq[me + 1]));
			});
		}
	}

	/**
	 * @param ps protein start, inclusive; {@code pe} protein end, exclusive
	 * @param ms match start, inclusive; {@code me} match end, inclusive
	 */
	static boolean trypticContext(byte[] seq, int ps, int pe, int ms, int me) {
		boolean left = ms == ps || isKR(seq[ms - 1]) || (ms == ps + 1 && (seq[ps] == 'M' || seq[ps] == 'm'));
		boolean right = me == pe - 1 || isKR(seq[me]);
		return left && right;
	}

	private static char upper(byte b) {
		return Character.toUpperCase((char) (b & 0xff));
	}

	private static boolean isKR(byte b) {
		return b == 'K' || b == 'R' || b == 'k' || b == 'r';
	}

	private static InputStream open(File f) throws IOException {
		InputStream in = new FileInputStream(f);
		if (f.getName().toLowerCase(Locale.ROOT).endsWith(".gz")) {
			return new GZIPInputStream(new BufferedInputStream(in, READ_BUFFER), READ_BUFFER);
		}
		return in;
	}

	private static final class Batch {
		static final Batch POISON = new Batch(-1);

		int count;
		long firstSerial;
		String[] ids;
		int[] start;
		int[] end;
		byte[] residues;
		int used;

		Batch(int count) {
			this.count = count;
		}

		static Batch create(long firstSerial) {
			Batch b = new Batch(0);
			b.firstSerial = firstSerial;
			b.ids = new String[BATCH_PROTEINS];
			b.start = new int[BATCH_PROTEINS];
			b.end = new int[BATCH_PROTEINS];
			b.residues = new byte[BATCH_RESIDUES + (64 << 10)];
			return b;
		}

		void append(byte c) {
			if (used == residues.length) {
				residues = Arrays.copyOf(residues, residues.length * 2);
			}
			residues[used++] = c;
		}
	}

	/** Byte-level FASTA parser; avoids a String per line on a 650-million-line file. */
	private final class Parser {
		private final BlockingQueue<Batch> queue;
		private final AtomicReference<Throwable> failure;
		private final Consumer<String> progress;
		private final String label;
		private final long startMillis;
		private Batch batch;
		private long serial;
		private long nextReport = 10_000_000L;
		private boolean atLineStart = true;
		private boolean inHeader;
		private boolean headerIdDone;
		private boolean keep;
		private final StringBuilder header = new StringBuilder();
		private int recordStart = -1;

		Parser(BlockingQueue<Batch> queue, AtomicReference<Throwable> failure, Consumer<String> progress,
				String label, long startMillis) {
			this.queue = queue;
			this.failure = failure;
			this.progress = progress;
			this.label = label;
			this.startMillis = startMillis;
			this.batch = Batch.create(0);
		}

		void parse(InputStream in) throws IOException, InterruptedException {
			byte[] buf = new byte[READ_BUFFER];
			atLineStart = true;
			int n;
			while ((n = in.read(buf)) > 0) {
				for (int i = 0; i < n; i++) {
					byte c = buf[i];
					if (c == '\n' || c == '\r') {
						if (inHeader) {
							endHeader();
						}
						atLineStart = true;
						continue;
					}
					if (atLineStart && c == '>') {
						endRecord();
						inHeader = true;
						headerIdDone = false;
						header.setLength(0);
						atLineStart = false;
						continue;
					}
					atLineStart = false;
					if (inHeader) {
						if (!headerIdDone) {
							if (c == ' ' || c == '\t') {
								headerIdDone = true;
							} else {
								header.append((char) (c & 0xff));
							}
						}
					} else if (keep && c != ' ' && c != '\t') {
						batch.append(c);
					}
				}
			}
			if (inHeader) {
				endHeader();
			}
			endRecord();
		}

		private void endHeader() {
			inHeader = false;
			String id = header.toString();
			keep = !id.isEmpty() && (decoyPrefix.isEmpty() || !id.toLowerCase(Locale.ROOT).startsWith(decoyPrefix));
			if (keep) {
				recordStart = batch.used;
				batch.ids[batch.count] = id;
			}
		}

		private void endRecord() throws IOException, InterruptedException {
			if (keep && recordStart >= 0) {
				batch.start[batch.count] = recordStart;
				batch.end[batch.count] = batch.used;
				batch.count++;
				serial++;
				if (batch.count == BATCH_PROTEINS || batch.used >= BATCH_RESIDUES) {
					put(batch);
					batch = Batch.create(serial);
				}
				if (progress != null && serial >= nextReport) {
					nextReport += 10_000_000L;
					progress.accept(String.format(Locale.ROOT, "Pangenome variant check: %s - %,d proteins scanned "
							+ "(%.0f s)", label, serial, (System.currentTimeMillis() - startMillis) / 1000.0));
				}
			}
			keep = false;
			recordStart = -1;
		}

		void flush() throws IOException, InterruptedException {
			endRecord();
			if (batch.count > 0) {
				put(batch);
			}
			batch = Batch.create(serial);
		}

		private void put(Batch b) throws IOException, InterruptedException {
			if (failure.get() != null) {
				throw matcherFailed(label, failure.get());
			}
			queue.put(b);
		}
	}
}

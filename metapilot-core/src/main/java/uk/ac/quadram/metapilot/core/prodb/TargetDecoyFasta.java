package uk.ac.quadram.metapilot.core.prodb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Builds and validates a 1:1 target/decoy FASTA, reusing
 * {@link FastaManager#addDecoy(String, String)} for the reversal itself.
 * <p>
 * Generation is idempotent, validated and safe under concurrent per-sample
 * tasks: a lock serialises writes, the file is built to a temporary path and
 * moved into place atomically, and an existing file is reused only if it is
 * newer than its target and actually contains one decoy per target. This
 * mirrors the contract of {@code MagFragPipeHeadlessRunner.ensureTargetDecoyFasta},
 * which is private; extracting it here lets the open-search pipeline share the
 * behaviour without touching the FragPipe runner.
 * <p>
 * <b>1:1, reversed, {@code rev_} only.</b> The FragPipe combined stage also
 * emits {@code shuf_} and {@code rev_shuf_} decoys, giving a 1:3 ratio. Those
 * exist solely because Philosopher strips {@code rev_} entries from the
 * quantified peptide list, leaving genome-level FDR without a decoy population
 * ({@code DDAMagFragPipeTask} calls {@code refineGenomes(quanPeps,
 * FastaManager.shuffle)}). The open-search pipeline emits its own PSM and
 * peptide lists with decoys retained throughout, so that workaround is
 * unnecessary here — and a 1:1 ratio keeps the target-decoy FDR estimator
 * {@code decoys/targets} correct as written, which at 1:3 would overstate FDR
 * roughly threefold.
 * <p>
 * <b>Decoy method.</b> {@link DecoyMethod#REVERSE} reverses the whole protein
 * (the default, unchanged). {@link DecoyMethod#PSEUDO_REVERSE} keeps the
 * cleavage sites: the protein is cut after every K/R, and in each segment the
 * first and last residues stay in place while the standard residues between
 * them are reversed - the rule {@code EntrapmentShuffle} uses, with reversal
 * instead of a random permutation. Every tryptic decoy peptide then has the
 * length, composition and precursor mass of its target, so a decoy can win the
 * weak spectra of a present peptide as often as a false same-mass target can.
 * Whole-protein reversal moves the cleavage sites and gives no such partner, so
 * its decoys cannot model those neighbour errors
 * ({@code docs/developer/peptide-fdr-and-precursor-aware-design.md}). The
 * {@code rev_} prefix is the same for both.
 */
public final class TargetDecoyFasta {

	private static final Logger LOGGER = LogManager.getLogger(TargetDecoyFasta.class);

	/** Serialises generation so parallel per-sample tasks don't race. */
	private static final Object LOCK = new Object();

	/** Decoy prefix written by {@link FastaManager#addDecoy}. */
	public static final String DECOY_TAG = "rev_";

	/** Residues {@link DecoyMethod#PSEUDO_REVERSE} moves; anything else (X, U, B, Z, ...) stays in place. */
	private static final String STANDARD = "ACDEFGHILMNPQSTVWY";

	/** How a decoy sequence is made from its target. */
	public enum DecoyMethod {
		/** Whole-protein reversal. */
		REVERSE("reverse"),
		/** Reversal within each cleavage segment, its first and last residues fixed. */
		PSEUDO_REVERSE("pseudo_reverse");

		public final String label;

		DecoyMethod(String label) {
			this.label = label;
		}

		/** The method named {@code label} (case-insensitive); {@link #REVERSE} for null or blank. */
		public static DecoyMethod fromLabel(String label) {
			if (label == null || label.isBlank()) {
				return REVERSE;
			}
			for (DecoyMethod m : values()) {
				if (m.label.equalsIgnoreCase(label.trim())) {
					return m;
				}
			}
			throw new IllegalArgumentException("Unknown decoy method: " + label + " (reverse or pseudo_reverse)");
		}
	}

	private TargetDecoyFasta() {
	}

	/** The decoy sequence of {@code target} under {@code method}. */
	public static String decoySequence(String target, DecoyMethod method) {
		if (method == DecoyMethod.PSEUDO_REVERSE) {
			return pseudoReverse(target);
		}
		return new StringBuilder(target).reverse().toString();
	}

	/**
	 * Within each segment (cut after every K/R, no proline rule) the standard
	 * residues strictly between the first and last positions are reversed. A
	 * segment whose reversal equals itself (I = L) has those residues rotated by
	 * one instead, so a decoy peptide differs from its target wherever the
	 * residues allow it.
	 */
	static String pseudoReverse(String sequence) {
		char[] seq = sequence.toCharArray();
		char[] out = seq.clone();
		int start = 0;
		for (int i = 0; i <= seq.length; i++) {
			boolean end = i == seq.length || ((seq[i] == 'K' || seq[i] == 'R') && i + 1 < seq.length);
			if (!end) {
				continue;
			}
			int stop = i == seq.length ? seq.length : i + 1;
			List<Integer> movable = new ArrayList<>();
			for (int j = start + 1; j < stop - 1; j++) {
				if (STANDARD.indexOf(seq[j]) >= 0) {
					movable.add(j);
				}
			}
			int n = movable.size();
			if (n >= 2) {
				for (int k = 0; k < n; k++) {
					out[movable.get(k)] = seq[movable.get(n - 1 - k)];
				}
				if (sameIgnoringIl(seq, out, start, stop)) {
					for (int k = 0; k < n; k++) {
						out[movable.get(k)] = seq[movable.get((k + 1) % n)];
					}
				}
			}
			start = stop;
			if (i == seq.length) {
				break;
			}
		}
		return new String(out);
	}

	private static boolean sameIgnoringIl(char[] a, char[] b, int from, int to) {
		for (int i = from; i < to; i++) {
			char x = a[i] == 'I' ? 'L' : a[i];
			char y = b[i] == 'I' ? 'L' : b[i];
			if (x != y) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Ensures {@code targetDecoyFasta} exists as a valid 1:1 target/decoy FASTA
	 * built from {@code targetFasta}, regenerating it if absent, stale or
	 * malformed.
	 *
	 * @return true if the file was (re)generated, false if an existing valid one was reused
	 */
	public static boolean ensure(File targetFasta, File targetDecoyFasta) throws IOException {
		return ensure(targetFasta, targetDecoyFasta, DecoyMethod.REVERSE);
	}

	/**
	 * As {@link #ensure(File, File)}, with decoys made by {@code method}. A file built with the other
	 * method is regenerated, and MSFragger's {@code .pepindex} files beside it are deleted with it so
	 * the next search cannot reuse an index of the old decoys.
	 */
	public static boolean ensure(File targetFasta, File targetDecoyFasta, DecoyMethod method) throws IOException {
		if (targetFasta == null || !targetFasta.exists()) {
			throw new IOException("Target FASTA not found: " + targetFasta);
		}
		if (targetDecoyFasta == null) {
			throw new IOException("Target-decoy FASTA path is null");
		}

		synchronized (LOCK) {
			if (isValid(targetFasta, targetDecoyFasta, method)) {
				return false;
			}
			File parent = targetDecoyFasta.getParentFile();
			if (parent != null && !parent.exists()) {
				parent.mkdirs();
			}
			File temp = new File(parent, targetDecoyFasta.getName() + ".tmp." + System.nanoTime());
			try {
				if (method == DecoyMethod.REVERSE) {
					FastaManager.addDecoy(targetFasta.getAbsolutePath(), temp.getAbsolutePath());
				} else {
					writeWithDecoys(targetFasta, temp, method);
				}
				if (!isValid(targetFasta, temp, method)) {
					throw new IOException(
							"Generated target-decoy FASTA failed validation: " + temp.getAbsolutePath());
				}
				moveReplacing(temp, targetDecoyFasta);
			} catch (IOException e) {
				Files.deleteIfExists(temp.toPath());
				throw e;
			}
			deleteStaleIndexes(targetDecoyFasta);
			Counts counts = count(targetDecoyFasta);
			LOGGER.info("Built target-decoy FASTA " + targetDecoyFasta.getAbsolutePath() + " (" + counts.targets
					+ " targets, " + counts.decoys + " decoys, tag=\"" + DECOY_TAG + "\", " + method.label + ")");
			return true;
		}
	}

	/** Each target record, then its decoy, in the layout {@link FastaManager#addDecoy} writes. */
	private static void writeWithDecoys(File targetFasta, File out, DecoyMethod method) throws IOException {
		try (BufferedReader reader = new BufferedReader(new FileReader(targetFasta));
				PrintWriter writer = new PrintWriter(out)) {
			String header = null;
			StringBuilder seq = new StringBuilder();
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					if (header != null) {
						writer.println(">" + DECOY_TAG + header.substring(1));
						writer.println(decoySequence(seq.toString(), method));
					}
					header = line;
					seq.setLength(0);
				} else if (header != null) {
					seq.append(line.trim());
				}
				writer.println(line);
			}
			if (header != null) {
				writer.println(">" + DECOY_TAG + header.substring(1));
				writer.println(decoySequence(seq.toString(), method));
			}
			if (writer.checkError()) {
				throw new IOException("Could not write " + out);
			}
		}
	}

	/** MSFragger's {@code <fasta>.<n>.pepindex} files of a database that was just rewritten. */
	private static void deleteStaleIndexes(File fasta) throws IOException {
		File parent = fasta.getAbsoluteFile().getParentFile();
		File[] indexes = parent == null ? null
				: parent.listFiles((dir, name) -> name.startsWith(fasta.getName() + ".") && name.endsWith(".pepindex"));
		if (indexes == null) {
			return;
		}
		for (File index : indexes) {
			Files.deleteIfExists(index.toPath());
		}
	}

	/**
	 * The method that made the decoys of {@code targetDecoyFasta}, read from the first target record
	 * and the decoy that follows it; null when neither method reproduces it, when both do, or when
	 * the file cannot be read.
	 */
	public static DecoyMethod methodOf(File targetDecoyFasta) {
		String[] pair = firstPair(targetDecoyFasta);
		if (pair == null) {
			return null;
		}
		boolean reverse = decoySequence(pair[0], DecoyMethod.REVERSE).equals(pair[1]);
		boolean pseudo = decoySequence(pair[0], DecoyMethod.PSEUDO_REVERSE).equals(pair[1]);
		return reverse == pseudo ? null : reverse ? DecoyMethod.REVERSE : DecoyMethod.PSEUDO_REVERSE;
	}

	/**
	 * Whether the first decoy is the {@code method} decoy of the first target. A file whose first pair
	 * cannot be read keeps the old rule (reverse only), so files built before the method existed stay
	 * valid.
	 */
	public static boolean decoysMadeBy(File targetDecoyFasta, DecoyMethod method) {
		String[] pair = firstPair(targetDecoyFasta);
		if (pair == null) {
			return method == DecoyMethod.REVERSE;
		}
		return decoySequence(pair[0], method).equals(pair[1]);
	}

	/** {first target sequence, the decoy sequence right after it}, or null. */
	private static String[] firstPair(File targetDecoyFasta) {
		try (BufferedReader reader = new BufferedReader(new FileReader(targetDecoyFasta))) {
			String target = null;
			boolean inDecoy = false;
			StringBuilder seq = new StringBuilder();
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					if (inDecoy) {
						break;
					}
					if (target == null && seq.length() > 0) {
						target = seq.toString();
					}
					if (target != null) {
						if (!line.substring(1).toLowerCase(Locale.ROOT).startsWith(DECOY_TAG)) {
							return null;
						}
						inDecoy = true;
					}
					seq.setLength(0);
				} else {
					// All whitespace: FastaManager.addDecoy reverses the raw lines, spaces included.
					seq.append(line.replaceAll("\\s+", ""));
				}
			}
			return inDecoy ? new String[] { target, seq.toString() } : null;
		} catch (IOException e) {
			return null;
		}
	}

	/**
	 * A target-decoy FASTA is valid when it is non-empty, at least as new as the
	 * target it derives from, and holds exactly one decoy per target. The
	 * one-to-one check is what makes the downstream {@code decoys/targets} FDR
	 * estimator trustworthy, so it is enforced rather than assumed.
	 */
	public static boolean isValid(File targetFasta, File targetDecoyFasta) {
		return isValid(targetFasta, targetDecoyFasta, DecoyMethod.REVERSE);
	}

	/** As {@link #isValid(File, File)}, and its decoys were made by {@code method}. */
	public static boolean isValid(File targetFasta, File targetDecoyFasta, DecoyMethod method) {
		if (targetDecoyFasta == null || !targetDecoyFasta.exists() || targetDecoyFasta.length() == 0) {
			return false;
		}
		if (!decoysMadeBy(targetDecoyFasta, method)) {
			return false;
		}
		if (targetFasta != null && targetDecoyFasta.lastModified() < targetFasta.lastModified()) {
			return false;
		}
		try {
			Counts counts = count(targetDecoyFasta);
			if (counts.targets == 0 || counts.targets != counts.decoys) {
				return false;
			}
			if (targetFasta != null && countHeaders(targetFasta) != counts.targets) {
				return false;
			}
			return true;
		} catch (IOException e) {
			return false;
		}
	}

	public static Counts count(File fasta) throws IOException {
		int targets = 0;
		int decoys = 0;
		try (BufferedReader reader = new BufferedReader(new FileReader(fasta))) {
			String line;
			while ((line = reader.readLine()) != null) {
				if (!line.startsWith(">")) {
					continue;
				}
				String ref = line.substring(1).toLowerCase(Locale.ROOT);
				if (ref.startsWith(DECOY_TAG)) {
					decoys++;
				} else {
					targets++;
				}
			}
		}
		return new Counts(targets, decoys);
	}

	private static int countHeaders(File fasta) throws IOException {
		int count = 0;
		try (BufferedReader reader = new BufferedReader(new FileReader(fasta))) {
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					count++;
				}
			}
		}
		return count;
	}

	private static void moveReplacing(File source, File target) throws IOException {
		try {
			Files.move(source.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING,
					StandardCopyOption.ATOMIC_MOVE);
		} catch (IOException atomicMoveUnsupported) {
			Files.move(source.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}

	public static final class Counts {
		public final int targets;
		public final int decoys;

		Counts(int targets, int decoys) {
			this.targets = targets;
			this.decoys = decoys;
		}

		public boolean isOneToOne() {
			return targets > 0 && targets == decoys;
		}
	}
}

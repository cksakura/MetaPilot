package uk.ac.quadram.metapilot.core.prodb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
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
 */
public final class TargetDecoyFasta {

	private static final Logger LOGGER = LogManager.getLogger(TargetDecoyFasta.class);

	/** Serialises generation so parallel per-sample tasks don't race. */
	private static final Object LOCK = new Object();

	/** Decoy prefix written by {@link FastaManager#addDecoy}. */
	public static final String DECOY_TAG = "rev_";

	private TargetDecoyFasta() {
	}

	/**
	 * Ensures {@code targetDecoyFasta} exists as a valid 1:1 target/decoy FASTA
	 * built from {@code targetFasta}, regenerating it if absent, stale or
	 * malformed.
	 *
	 * @return true if the file was (re)generated, false if an existing valid one was reused
	 */
	public static boolean ensure(File targetFasta, File targetDecoyFasta) throws IOException {
		if (targetFasta == null || !targetFasta.exists()) {
			throw new IOException("Target FASTA not found: " + targetFasta);
		}
		if (targetDecoyFasta == null) {
			throw new IOException("Target-decoy FASTA path is null");
		}

		synchronized (LOCK) {
			if (isValid(targetFasta, targetDecoyFasta)) {
				return false;
			}
			File parent = targetDecoyFasta.getParentFile();
			if (parent != null && !parent.exists()) {
				parent.mkdirs();
			}
			File temp = new File(parent, targetDecoyFasta.getName() + ".tmp." + System.nanoTime());
			try {
				FastaManager.addDecoy(targetFasta.getAbsolutePath(), temp.getAbsolutePath());
				if (!isValid(targetFasta, temp)) {
					throw new IOException(
							"Generated target-decoy FASTA failed validation: " + temp.getAbsolutePath());
				}
				moveReplacing(temp, targetDecoyFasta);
			} catch (IOException e) {
				Files.deleteIfExists(temp.toPath());
				throw e;
			}
			Counts counts = count(targetDecoyFasta);
			LOGGER.info("Built target-decoy FASTA " + targetDecoyFasta.getAbsolutePath() + " (" + counts.targets
					+ " targets, " + counts.decoys + " decoys, tag=\"" + DECOY_TAG + "\")");
			return true;
		}
	}

	/**
	 * A target-decoy FASTA is valid when it is non-empty, at least as new as the
	 * target it derives from, and holds exactly one decoy per target. The
	 * one-to-one check is what makes the downstream {@code decoys/targets} FDR
	 * estimator trustworthy, so it is enforced rather than assumed.
	 */
	public static boolean isValid(File targetFasta, File targetDecoyFasta) {
		if (targetDecoyFasta == null || !targetDecoyFasta.exists() || targetDecoyFasta.length() == 0) {
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

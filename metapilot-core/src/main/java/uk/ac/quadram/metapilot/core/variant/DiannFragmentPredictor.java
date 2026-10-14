package uk.ac.quadram.metapilot.core.variant;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Predicted b/y fragment intensities from DIA-NN 1.8.x, through the same
 * prediction interface MSBooster uses and with the same DIA-NN the open
 * search's stage-1 MSBooster step is pointed at.
 * <p>
 * Input: a TSV with a {@code peptide\tcharge} header, one precursor per row,
 * modifications in DIA-NN's {@code [UniMod:n]} notation. Command:
 * {@code diann --lib <in>.tsv --predict --threads N --strip-unknown-mods
 * --predict-n-frag 100}. Output: {@code <in>.predicted.bin} beside the input,
 * one record per input row in the same order, little-endian:
 * <pre>
 * int32 nFragments, float32 iRT, float32 ion mobility,
 * nFragments x int32: bits 31..16 intensity (uint16), 15..8 ion number,
 *                     bit 2 y (1) / b (0), bits 1..0 charge - 1
 * </pre>
 * The ion number is the residue count of the fragment: b_n holds the first n
 * residues, y_n the last n. Format as read by MSBooster's
 * {@code DiannSpeclibReader} (Apache-2.0).
 */
final class DiannFragmentPredictor {

	static final String INPUT_FILE = "variant_forms.tsv";
	static final String OUTPUT_FILE = "variant_forms.predicted.bin";
	private static final int MAX_FRAGMENTS = 4096;
	private static final int LOG_TAIL = 30;
	/** Rows per DIA-NN call; bounds its memory the way MSBooster's prediction limit does. */
	static final int BATCH_LIMIT = 500_000;
	private static final int EXIT_ACCESS_VIOLATION = -1073741819;
	private static final int EXIT_MISSING_VC_RUNTIME = -1073741515;

	/** Predicted fragments of one precursor. */
	static final class Prediction {
		final int[] ionNumber;
		final boolean[] yIon;
		final int[] charge;
		final float[] intensity;

		Prediction(int[] ionNumber, boolean[] yIon, int[] charge, float[] intensity) {
			this.ionNumber = ionNumber;
			this.yIon = yIon;
			this.charge = charge;
			this.intensity = intensity;
		}

		int size() {
			return intensity.length;
		}
	}

	private DiannFragmentPredictor() {
	}

	static String key(String peptide, int charge) {
		return peptide + "/" + charge;
	}

	/**
	 * Predicts every {@code (peptide, charge)} pair.
	 *
	 * @param peptides        plain sequences, upper-case residues
	 * @param charges         precursor charge of each peptide
	 * @param carbamidomethyl whether C carries the fixed carbamidomethyl
	 * @return predictions keyed by {@link #key}
	 */
	static Map<String, Prediction> predict(File diann, List<String> peptides, List<Integer> charges,
			boolean carbamidomethyl, File workDir, int threads, Consumer<String> log) throws IOException {
		if (!workDir.isDirectory() && !workDir.mkdirs()) {
			throw new IOException("Cannot create " + workDir);
		}
		Map<String, Prediction> out = new HashMap<>();
		int batches = (peptides.size() + BATCH_LIMIT - 1) / BATCH_LIMIT;
		for (int b = 0; b < batches; b++) {
			int from = b * BATCH_LIMIT;
			int to = Math.min(peptides.size(), from + BATCH_LIMIT);
			List<String> pep = peptides.subList(from, to);
			List<Integer> chg = charges.subList(from, to);
			if (batches > 1) {
				log.accept("MS2 form ranking: DIA-NN batch " + (b + 1) + " of " + batches);
			}
			List<Prediction> predictions;
			try {
				predictions = predictBatch(diann, pep, chg, carbamidomethyl, workDir, threads, log);
			} catch (AccessViolation e) {
				log.accept("MS2 form ranking: DIA-NN stopped with an access violation; retrying once");
				predictions = predictBatch(diann, pep, chg, carbamidomethyl, workDir, threads, log);
			}
			for (int i = 0; i < pep.size(); i++) {
				out.put(key(pep.get(i), chg.get(i)), predictions.get(i));
			}
		}
		return out;
	}

	/** DIA-NN's sporadic 0xC0000005 on Windows, which MSBooster also retries once. */
	private static final class AccessViolation extends IOException {
		private static final long serialVersionUID = 1L;

		AccessViolation(String message) {
			super(message);
		}
	}

	private static List<Prediction> predictBatch(File diann, List<String> peptides, List<Integer> charges,
			boolean carbamidomethyl, File workDir, int threads, Consumer<String> log) throws IOException {
		File input = new File(workDir, INPUT_FILE);
		File output = new File(workDir, OUTPUT_FILE);
		Files.deleteIfExists(output.toPath());
		writeInput(input, peptides, charges, carbamidomethyl);

		List<String> cmd = List.of(diann.getAbsolutePath(), "--lib", input.getAbsolutePath(), "--predict",
				"--threads", String.valueOf(Math.max(1, threads)), "--strip-unknown-mods", "--predict-n-frag", "100");
		log.accept("MS2 form ranking: predicting " + peptides.size() + " precursors with " + String.join(" ", cmd));
		ProcessBuilder pb = new ProcessBuilder(cmd);
		pb.directory(workDir);
		pb.redirectErrorStream(true);
		Process p = pb.start();
		Deque<String> tail = new ArrayDeque<>();
		try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (Thread.currentThread().isInterrupted()) {
					p.destroyForcibly();
					throw new IOException("Interrupted while DIA-NN was predicting");
				}
				if (tail.size() == LOG_TAIL) {
					tail.removeFirst();
				}
				tail.addLast(line);
			}
		}
		int exit;
		try {
			exit = p.waitFor();
		} catch (InterruptedException e) {
			p.destroyForcibly();
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted while DIA-NN was predicting", e);
		}
		if (exit == EXIT_ACCESS_VIOLATION) {
			throw new AccessViolation("DIA-NN prediction stopped with an access violation (exit " + exit + ")");
		}
		if (exit == EXIT_MISSING_VC_RUNTIME) {
			throw new IOException("DIA-NN could not start: the Microsoft Visual C++ Redistributable is missing "
					+ "(exit " + exit + ")");
		}
		if (exit != 0 || !output.isFile()) {
			throw new IOException("DIA-NN prediction failed (exit " + exit + (output.isFile() ? "" : ", no " + OUTPUT_FILE)
					+ "). It needs the 1.8.x prediction interface FragPipe bundles; a 2.x build has no --predict. "
					+ "Last output: " + String.join(" | ", tail));
		}
		return readPredictions(output, peptides.size());
	}

	static void writeInput(File input, List<String> peptides, List<Integer> charges, boolean carbamidomethyl)
			throws IOException {
		try (PrintWriter w = new PrintWriter(input, StandardCharsets.UTF_8)) {
			w.print("peptide\tcharge\n");
			for (int i = 0; i < peptides.size(); i++) {
				w.print(diannSequence(peptides.get(i), carbamidomethyl));
				w.print('\t');
				w.print(charges.get(i));
				w.print('\n');
			}
		}
	}

	static String diannSequence(String peptide, boolean carbamidomethyl) {
		if (!carbamidomethyl || peptide.indexOf('C') < 0) {
			return peptide;
		}
		return peptide.replace("C", "C[UniMod:4]");
	}

	/** Reads exactly {@code expected} records; anything else is a mismatch with the input. */
	static List<Prediction> readPredictions(File bin, int expected) throws IOException {
		List<Prediction> out = new ArrayList<>(expected);
		try (InputStream raw = new BufferedInputStream(Files.newInputStream(bin.toPath()), 1 << 16);
				DataInputStream in = new DataInputStream(raw)) {
			for (int r = 0; r < expected; r++) {
				int n;
				try {
					n = Integer.reverseBytes(in.readInt());
				} catch (EOFException e) {
					throw new IOException(bin.getName() + " holds " + r + " predictions for " + expected
							+ " input rows; DIA-NN and the input disagree");
				}
				if (n < 0 || n > MAX_FRAGMENTS) {
					throw new IOException(bin.getName() + ": implausible fragment count " + n + " at record " + r);
				}
				in.readInt();
				in.readInt();
				int[] ion = new int[n];
				boolean[] y = new boolean[n];
				int[] charge = new int[n];
				float[] intensity = new float[n];
				for (int i = 0; i < n; i++) {
					int f = Integer.reverseBytes(in.readInt());
					intensity[i] = (f >>> 16) & 0xFFFF;
					ion[i] = (f >>> 8) & 0xFF;
					y[i] = ((f >>> 2) & 0x1) == 1;
					charge[i] = (f & 0x3) + 1;
				}
				out.add(new Prediction(ion, y, charge, intensity));
			}
			if (in.read() != -1) {
				throw new IOException(bin.getName() + " holds more predictions than the " + expected
						+ " input rows; DIA-NN and the input disagree");
			}
		}
		return out;
	}
}

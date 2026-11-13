package uk.ac.quadram.metapilot.core.variant;

import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

/**
 * Stand-in for {@code diann --lib <in>.tsv --predict --threads N ...} in tests: writes
 * {@code <in>.predicted.bin} with one fragment-less record per input row, whose iRT is the number
 * after "PEP" in the peptide, or -1 when it was not asked to use exactly one thread.
 */
public final class FakeDiann {

	private FakeDiann() {
	}

	public static void main(String[] args) throws IOException {
		File lib = null;
		String threads = null;
		for (int i = 0; i + 1 < args.length; i++) {
			if (args[i].equals("--lib")) {
				lib = new File(args[i + 1]);
			} else if (args[i].equals("--threads")) {
				threads = args[i + 1];
			}
		}
		List<String> rows = Files.readAllLines(lib.toPath(), StandardCharsets.UTF_8);
		File out = new File(lib.getParentFile(), DiannFragmentPredictor.OUTPUT_FILE);
		try (DataOutputStream w = new DataOutputStream(Files.newOutputStream(out.toPath()))) {
			for (String row : rows.subList(1, rows.size())) {
				String peptide = row.substring(0, row.indexOf('\t'));
				float iRT = "1".equals(threads) ? Float.parseFloat(peptide.substring(3)) : -1f;
				w.writeInt(Integer.reverseBytes(0));
				w.writeInt(Integer.reverseBytes(Float.floatToIntBits(iRT)));
				w.writeInt(0);
			}
		}
	}
}

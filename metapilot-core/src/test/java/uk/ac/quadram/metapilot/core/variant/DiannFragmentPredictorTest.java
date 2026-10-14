package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class DiannFragmentPredictorTest {

	@TempDir
	File tmp;

	/** {@code {intensity, ionNumber, y ? 1 : 0, charge}} per fragment. */
	private static void record(ByteBuffer b, int[]... fragments) {
		b.putInt(fragments.length);
		b.putFloat(12.5f);
		b.putFloat(0.9f);
		for (int[] f : fragments) {
			b.putInt(((f[0] & 0xFFFF) << 16) | ((f[1] & 0xFF) << 8) | (f[2] << 2) | (f[3] - 1));
		}
	}

	private File bin(int records) throws IOException {
		ByteBuffer b = ByteBuffer.allocate(1024).order(ByteOrder.LITTLE_ENDIAN);
		record(b, new int[] { 40000, 3, 1, 1 }, new int[] { 1200, 2, 0, 2 });
		if (records > 1) {
			record(b);
		}
		File f = new File(tmp, DiannFragmentPredictor.OUTPUT_FILE);
		byte[] bytes = new byte[b.position()];
		b.flip();
		b.get(bytes);
		Files.write(f.toPath(), bytes);
		return f;
	}

	@Test
	@DisplayName("Input TSV: peptide and charge, carbamidomethyl C in UniMod notation")
	void input() throws IOException {
		File in = new File(tmp, DiannFragmentPredictor.INPUT_FILE);
		DiannFragmentPredictor.writeInput(in, List.of("PEPCK", "AASDK"), List.of(2, 3), true);
		assertEquals(List.of("peptide\tcharge", "PEPC[UniMod:4]K\t2", "AASDK\t3"),
				Files.readAllLines(in.toPath(), StandardCharsets.UTF_8));
		assertEquals("PEPCK", DiannFragmentPredictor.diannSequence("PEPCK", false));
	}

	@Test
	@DisplayName(".predicted.bin: little-endian records decoded in input order")
	void roundTrip() throws IOException {
		List<DiannFragmentPredictor.Prediction> p = DiannFragmentPredictor.readPredictions(bin(2), 2);
		assertEquals(2, p.size());
		DiannFragmentPredictor.Prediction first = p.get(0);
		assertEquals(2, first.size());
		assertEquals(40000f, first.intensity[0], "intensity is unsigned 16-bit");
		assertEquals(3, first.ionNumber[0]);
		assertTrue(first.yIon[0]);
		assertEquals(1, first.charge[0]);
		assertEquals(1200f, first.intensity[1]);
		assertEquals(2, first.ionNumber[1]);
		assertFalse(first.yIon[1]);
		assertEquals(2, first.charge[1]);
		assertEquals(0, p.get(1).size());
	}

	@Test
	@DisplayName("A record count different from the input is an error, not a misalignment")
	void countMismatch() throws IOException {
		File two = bin(2);
		assertThrows(IOException.class, () -> DiannFragmentPredictor.readPredictions(two, 3));
		assertThrows(IOException.class, () -> DiannFragmentPredictor.readPredictions(two, 1));
	}
}

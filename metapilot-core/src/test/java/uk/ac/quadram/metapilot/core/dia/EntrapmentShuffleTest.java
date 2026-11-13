package uk.ac.quadram.metapilot.core.dia;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;

class EntrapmentShuffleTest {

	@TempDir
	File dir;

	private static final String PROTEIN =
			"MLSFLEYSAYEIIATSIYIINFIVILILIFKEKRSIETIVAWTAVLTILPAGGFILYMLVGRGVAKDNMFKIKEQEDEIIKENISKTKDILKPYSSMHDENLK";

	private static List<String> tryptic(String seq) {
		List<String> out = new ArrayList<>();
		for (int[] s : EntrapmentShuffle.segments(seq)) out.add(seq.substring(s[0], s[1]));
		return out;
	}

	@Test
	void shuffledProteinDigestsLikeItsSource() {
		String sh = EntrapmentShuffle.shuffle(PROTEIN, "MGYG800000001_00001", 1);
		assertEquals(PROTEIN.length(), sh.length());
		assertNotEquals(PROTEIN, sh);
		List<String> a = tryptic(PROTEIN), b = tryptic(sh);
		assertEquals(a.size(), b.size());
		for (int i = 0; i < a.size(); i++) {
			String x = a.get(i), y = b.get(i);
			assertEquals(x.length(), y.length());
			assertEquals(x.charAt(0), y.charAt(0));                       // proline rule and N-terminal Met kept
			assertEquals(x.charAt(x.length() - 1), y.charAt(y.length() - 1)); // K/R site kept
			char[] cx = x.toCharArray(), cy = y.toCharArray();
			Arrays.sort(cx);
			Arrays.sort(cy);
			assertEquals(new String(cx), new String(cy));                 // same composition and mass
		}
		assertEquals(sh, EntrapmentShuffle.shuffle(PROTEIN, "MGYG800000001_00001", 1)); // reproducible
		assertNotEquals(sh, EntrapmentShuffle.shuffle(PROTEIN, "MGYG800000001_00001", 2));
	}

	@Test
	void idsKeepSuffixAndVersion() {
		assertEquals("MGYG800000001_00005", EntrapmentShuffle.entrapmentId("MGYG000000001_00005", 8));
		assertEquals("MGYG800000761.1_00014", EntrapmentShuffle.entrapmentId("MGYG000000761.1_00014", 8));
		assertEquals("MGYG800000761.1", EntrapmentShuffle.entrapmentId("MGYG000000761.1", 8));
		assertThrows(IllegalArgumentException.class, () -> EntrapmentShuffle.entrapmentId("MGYG900000001_00005", 8));
	}

	private static DiaNNPrecursor precursor(String... proteins) {
		return new DiaNNPrecursor("run", proteins, "PEPTIDEK", "PEPTIDEK", 2, 0, 0, 0.001, 0, 0.001, 0.001, 0.001, 0, 0);
	}

	@Test
	void targetAndEntrapmentClassification() {
		DiaNNPrecursor target = precursor("MGYG000000001_00001");
		DiaNNPrecursor entr = precursor("MGYG800000001_00001", "REV_MGYG000000002_00003");
		DiaNNPrecursor both = precursor("MGYG800000001_00001", "MGYG000000002_00003");
		assertTrue(EntrapmentShuffle.hasTargetProtein(target, 8));
		assertFalse(EntrapmentShuffle.onlyEntrapmentProteins(target, 8));
		assertFalse(EntrapmentShuffle.hasTargetProtein(entr, 8));
		assertTrue(EntrapmentShuffle.onlyEntrapmentProteins(entr, 8));
		assertTrue(EntrapmentShuffle.hasTargetProtein(both, 8));   // shared with a target: counts as target
		assertFalse(EntrapmentShuffle.onlyEntrapmentProteins(both, 8));
	}

	@Test
	void appendsOneShuffledCopyPerProtein() throws Exception {
		File faa = new File(dir, "MGYG000000001.faa");
		Files.write(faa.toPath(), List.of(">MGYG000000001_00001 Major cardiolipin synthase ClsA",
				PROTEIN.substring(0, 60), PROTEIN.substring(60), ">MGYG000000001_00002 Elongation factor Tu",
				"MAKEKFERTKPHVNVGTIGHVDHGKTTLTAAITK"), StandardCharsets.UTF_8);
		StringWriter sw = new StringWriter();
		try (PrintWriter w = new PrintWriter(sw)) {
			assertEquals(2, EntrapmentShuffle.appendShuffledGenome(faa, w, 8, 1));
		}
		String[] lines = sw.toString().split("\\R");
		assertEquals(">MGYG800000001_00001 Major cardiolipin synthase ClsA", lines[0]);
		assertEquals(EntrapmentShuffle.shuffle(PROTEIN, "MGYG800000001_00001", 1), lines[1]);
		assertEquals(">MGYG800000001_00002 Elongation factor Tu", lines[2]);
	}
}

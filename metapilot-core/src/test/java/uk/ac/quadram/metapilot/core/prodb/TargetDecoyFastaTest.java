package uk.ac.quadram.metapilot.core.prodb;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta.DecoyMethod;

class TargetDecoyFastaTest {
	@TempDir
	Path dir;

	private static final String SEQ = "MACDEFGHILMNPQSTVWYRKACDEFGHILMNPQSTVWYKAAGPLW";

	/** Tryptic peptides (cut after K/R, no proline rule), as the segments the decoy keeps. */
	private static List<String> segments(String seq) {
		List<String> out = new ArrayList<>();
		int start = 0;
		for (int i = 0; i < seq.length(); i++) {
			if ((seq.charAt(i) == 'K' || seq.charAt(i) == 'R') && i + 1 < seq.length()) {
				out.add(seq.substring(start, i + 1));
				start = i + 1;
			}
		}
		out.add(seq.substring(start));
		return out;
	}

	private static String sorted(String s) {
		char[] c = s.toCharArray();
		Arrays.sort(c);
		return new String(c);
	}

	@Test
	void pseudoReverseKeepsEverySegmentsEndsAndComposition() {
		String decoy = TargetDecoyFasta.decoySequence(SEQ, DecoyMethod.PSEUDO_REVERSE);
		List<String> t = segments(SEQ);
		List<String> d = segments(decoy);
		assertEquals(t.size(), d.size(), "cleavage sites unchanged");
		for (int i = 0; i < t.size(); i++) {
			String a = t.get(i), b = d.get(i);
			assertEquals(a.length(), b.length());
			assertEquals(a.charAt(0), b.charAt(0), "first residue fixed");
			assertEquals(a.charAt(a.length() - 1), b.charAt(b.length() - 1), "last residue fixed");
			assertEquals(sorted(a), sorted(b), "same composition, hence same mass");
		}
		assertEquals("MYWVTSQPNMLIHGFEDCAR", d.get(0), "interior reversed");
		assertNotEquals(SEQ, decoy);
	}

	@Test
	void palindromicSegmentIsRotatedInstead() {
		// Interior ACA reads the same reversed; the rotation gives a different peptide.
		String target = "GACAK";
		String decoy = TargetDecoyFasta.decoySequence(target, DecoyMethod.PSEUDO_REVERSE);
		assertNotEquals(target, decoy);
		assertEquals('G', decoy.charAt(0));
		assertEquals('K', decoy.charAt(decoy.length() - 1));
		assertEquals(sorted(target), sorted(decoy));
	}

	@Test
	void reverseIsWholeProteinReversal() {
		assertEquals(new StringBuilder(SEQ).reverse().toString(),
				TargetDecoyFasta.decoySequence(SEQ, DecoyMethod.REVERSE));
	}

	@Test
	void ensureBuildsDetectsAndRebuildsOnMethodChange() throws Exception {
		File target = dir.resolve("t.fasta").toFile();
		Files.writeString(target.toPath(), ">P1 one\n" + SEQ.substring(0, 20) + "\n" + SEQ.substring(20) + "\n>P2\nGACAGKLLR\n");
		File td = dir.resolve("td.fasta").toFile();

		assertTrue(TargetDecoyFasta.ensure(target, td));
		assertEquals(DecoyMethod.REVERSE, TargetDecoyFasta.methodOf(td));
		assertFalse(TargetDecoyFasta.ensure(target, td), "reverse file reused");
		assertFalse(TargetDecoyFasta.isValid(target, td, DecoyMethod.PSEUDO_REVERSE));

		File index = dir.resolve("td.fasta.1.pepindex").toFile();
		Files.writeString(index.toPath(), "stale");
		assertTrue(TargetDecoyFasta.ensure(target, td, DecoyMethod.PSEUDO_REVERSE), "other method rebuilds");
		assertFalse(index.exists(), "index of the old decoys deleted");
		assertEquals(DecoyMethod.PSEUDO_REVERSE, TargetDecoyFasta.methodOf(td));
		assertFalse(TargetDecoyFasta.ensure(target, td, DecoyMethod.PSEUDO_REVERSE), "pseudo file reused");

		List<String> lines = Files.readAllLines(td.toPath());
		assertEquals(">rev_P1 one", lines.get(3));
		assertEquals(TargetDecoyFasta.decoySequence(SEQ, DecoyMethod.PSEUDO_REVERSE), lines.get(4));
		assertEquals(">rev_P2", lines.get(7));
		TargetDecoyFasta.Counts counts = TargetDecoyFasta.count(td);
		assertEquals(2, counts.targets);
		assertEquals(2, counts.decoys);
	}

	@Test
	void methodLabelsParse() {
		assertEquals(DecoyMethod.REVERSE, DecoyMethod.fromLabel(null));
		assertEquals(DecoyMethod.PSEUDO_REVERSE, DecoyMethod.fromLabel("Pseudo_Reverse"));
		assertThrows(IllegalArgumentException.class, () -> DecoyMethod.fromLabel("shuffle"));
	}
}

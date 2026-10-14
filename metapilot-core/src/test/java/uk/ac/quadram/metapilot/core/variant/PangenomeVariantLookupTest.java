package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PangenomeVariantLookupTest {

	@TempDir
	File tmp;

	private File fasta(String name, String content) throws IOException {
		File f = new File(tmp, name);
		if (name.endsWith(".gz")) {
			try (OutputStream out = new GZIPOutputStream(Files.newOutputStream(f.toPath()))) {
				out.write(content.getBytes(StandardCharsets.US_ASCII));
			}
		} else {
			Files.write(f.toPath(), content.getBytes(StandardCharsets.US_ASCII));
		}
		return f;
	}

	private static boolean[] found(List<String> queries, File fasta, boolean tryptic) throws IOException {
		List<String> normalised = new ArrayList<>();
		for (String q : queries) {
			normalised.add(PeptideAutomaton.normalise(q));
		}
		PangenomeVariantLookup.QueryHit[] hits = PangenomeVariantLookup.scan(normalised, List.of(fasta), null,
				tryptic, 2, "rev_", null, "test");
		boolean[] out = new boolean[hits.length];
		for (int i = 0; i < hits.length; i++) {
			out[i] = hits[i].found();
		}
		return out;
	}

	@Test
	@DisplayName("A peptide split across a FASTA line break is found")
	void multiLineRecord() throws IOException {
		File f = fasta("a.fasta", ">P1 desc\nMKAAAGG\nGGKCCCR\n>P2\nWWWW\n");
		assertTrue(found(List.of("AAAGGGGK"), f, true)[0]);
	}

	@Test
	@DisplayName("Tryptic context: preceded by K/R or protein start, ending in K/R or protein end")
	void trypticContext() throws IOException {
		File f = fasta("b.fasta", ">P1\nAAKPEPTIDEKBBBR\n>P2\nGGGTIDEKGG\n");
		boolean[] r = found(List.of("PEPTIDEK", "EPTIDEK", "TIDEK", "AAK", "GGG"), f, true);
		assertTrue(r[0], "after K, ends in K");
		assertFalse(r[1], "starts after P, not after K/R");
		assertFalse(r[2], "TIDEK in P2 is preceded by G");
		assertTrue(r[3], "protein start, ends in K");
		assertFalse(r[4], "protein start, but ends in G away from the protein end");
	}

	@Test
	@DisplayName("Protein N-terminal methionine removal counts as a cleavage site")
	void methionineLoss() throws IOException {
		File f = fasta("c.fasta", ">P1\nMPEPTIDEKAAAA\n>P2\nGMPEPTIDER\n");
		boolean[] r = found(List.of("PEPTIDEK", "PEPTIDER"), f, true);
		assertTrue(r[0]);
		assertFalse(r[1], "M at position 2 is not the initiator");
	}

	@Test
	@DisplayName("I and L are one residue, in queries and proteins")
	void isoleucineLeucine() throws IOException {
		File f = fasta("d.fasta", ">P1\nKPEPTLDEK\n");
		assertTrue(found(List.of("PEPTIDEK"), f, true)[0]);
	}

	@Test
	@DisplayName("Decoy entries are skipped and gzip input is read")
	void decoysAndGzip() throws IOException {
		File f = fasta("e.fasta.gz", ">rev_P1\nKPEPTIDEK\n>P2\nKAAAAK\n");
		boolean[] r = found(List.of("PEPTIDEK", "AAAAK"), f, true);
		assertFalse(r[0]);
		assertTrue(r[1]);
	}

	@Test
	@DisplayName("A query containing X never matches; matching restarts after an X in the protein")
	void nonResidueResets() throws IOException {
		File f = fasta("f.fasta", ">P1\nPEPXTIDEK\n");
		boolean[] r = found(List.of("PEPXTIDEK", "TIDEK"), f, false);
		assertFalse(r[0]);
		assertTrue(r[1]);
	}

	@Test
	@DisplayName("Example protein is the first hit in the file and proteins are counted once each")
	void exampleProtein() throws IOException {
		File f = fasta("g.fasta", ">P1\nGGG\n>P2\nKPEPKPEPK\n>P3\nKPEPK\n");
		PangenomeVariantLookup.QueryHit[] hits = PangenomeVariantLookup.scan(List.of("PEPK"), List.of(f), null,
				true, 3, "", null, "test");
		assertEquals(2, hits[0].proteins);
		assertEquals("P2", hits[0].exampleProtein);
	}

	@Test
	@DisplayName("Without tryptic context the automaton agrees with a brute-force substring search")
	void bruteForceAgreement() throws IOException {
		Random rnd = new Random(7);
		// Small alphabet so short queries hit often; I and L included to exercise the fold.
		String alphabet = "ACILKR";
		StringBuilder content = new StringBuilder();
		List<String> proteins = new ArrayList<>();
		for (int p = 0; p < 300; p++) {
			StringBuilder seq = new StringBuilder();
			int len = 20 + rnd.nextInt(200);
			for (int i = 0; i < len; i++) {
				seq.append(alphabet.charAt(rnd.nextInt(alphabet.length())));
			}
			proteins.add(seq.toString());
			content.append(">prot").append(p).append('\n');
			for (int i = 0; i < seq.length(); i += 60) {
				content.append(seq, i, Math.min(seq.length(), i + 60)).append('\n');
			}
		}
		Set<String> querySet = new LinkedHashSet<>();
		while (querySet.size() < 400) {
			if (rnd.nextBoolean()) {
				String p = proteins.get(rnd.nextInt(proteins.size()));
				int len = 3 + rnd.nextInt(6);
				int start = rnd.nextInt(p.length() - len);
				querySet.add(p.substring(start, start + len));
			} else {
				StringBuilder q = new StringBuilder();
				int len = 3 + rnd.nextInt(5);
				for (int i = 0; i < len; i++) {
					q.append(alphabet.charAt(rnd.nextInt(alphabet.length())));
				}
				querySet.add(q.toString());
			}
		}
		Set<String> normalisedSet = new LinkedHashSet<>();
		for (String q : querySet) {
			normalisedSet.add(PeptideAutomaton.normalise(q));
		}
		List<String> queries = new ArrayList<>(normalisedSet);
		File f = fasta("h.fasta", content.toString());
		PangenomeVariantLookup.QueryHit[] hits = PangenomeVariantLookup.scan(queries, List.of(f), null, false, 4,
				"", null, "test");
		for (int q = 0; q < queries.size(); q++) {
			int expected = 0;
			for (String p : proteins) {
				if (PeptideAutomaton.normalise(p).contains(queries.get(q))) {
					expected++;
				}
			}
			assertEquals(expected, hits[q].proteins, "query " + queries.get(q));
		}
	}

	@Test
	@DisplayName("Tryptic context check on raw offsets")
	void trypticContextOffsets() {
		byte[] seq = "MKPEPTIDERAA".getBytes(StandardCharsets.US_ASCII);
		assertTrue(FastaScanner.trypticContext(seq, 0, seq.length, 2, 9));
		assertFalse(FastaScanner.trypticContext(seq, 0, seq.length, 3, 9));
		assertTrue(FastaScanner.trypticContext(seq, 0, seq.length, 10, 11));
		assertTrue(FastaScanner.trypticContext(seq, 0, seq.length, 1, 1));
	}
}

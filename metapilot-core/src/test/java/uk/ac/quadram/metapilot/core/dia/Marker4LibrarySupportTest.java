package uk.ac.quadram.metapilot.core.dia;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;

class Marker4LibrarySupportTest {

	@TempDir
	File tempDir;

	@Test
	void identityMatchesOnEqualFields() {
		Marker4LibrarySupport.LibraryIdentity a = new Marker4LibrarySupport.LibraryIdentity(
				"Marker4_highpep.fasta", "abc", "mExcision=true;carbamiC=true;oxidationM=false;acetyProN=false;phosphoSTY=false;kGG=false",
				"1.8.1");
		Marker4LibrarySupport.LibraryIdentity b = new Marker4LibrarySupport.LibraryIdentity(
				"Marker4_highpep.fasta", "ABC", "mExcision=true;carbamiC=true;oxidationM=false;acetyProN=false;phosphoSTY=false;kGG=false",
				"1.8.1");
		assertTrue(a.matches(b));
	}

	@Test
	void identityDetectsModificationMismatch() {
		Marker4LibrarySupport.LibraryIdentity cached = new Marker4LibrarySupport.LibraryIdentity(
				"Marker4.fasta", "abc", "mExcision=true;carbamiC=true;oxidationM=false;acetyProN=false;phosphoSTY=false;kGG=false",
				"1.8.1");
		Marker4LibrarySupport.LibraryIdentity expected = new Marker4LibrarySupport.LibraryIdentity(
				"Marker4.fasta", "abc", "mExcision=true;carbamiC=true;oxidationM=true;acetyProN=false;phosphoSTY=false;kGG=false",
				"1.8.1");
		assertFalse(cached.matches(expected));
		assertTrue(cached.describeMismatch(expected).contains("Modifications"));
	}

	@Test
	void writeAndReadInfoRoundTrip() throws Exception {
		File fasta = new File(tempDir, "Marker4_highpep.fasta");
		try (PrintWriter w = new PrintWriter(fasta, StandardCharsets.UTF_8)) {
			w.println(">p1");
			w.println("PEPTIDE");
			w.println(">p2");
			w.println("OTHER");
		}
		String sha = Marker4LibrarySupport.sha256Hex(fasta);
		Marker4LibrarySupport.LibraryIdentity id = new Marker4LibrarySupport.LibraryIdentity(
				fasta.getName(), sha, "mExcision=true;carbamiC=true;oxidationM=false;acetyProN=false;phosphoSTY=false;kGG=false",
				"2.0");
		File info = new File(tempDir, "Marker4_info.txt");
		Marker4LibrarySupport.writeInfoFile(info, 42, id);

		Map<String, String> map = Marker4LibrarySupport.readInfoMap(info);
		assertEquals("42", map.get(Marker4LibrarySupport.KEY_DATA_SET_SIZE));
		Marker4LibrarySupport.LibraryIdentity read = Marker4LibrarySupport.identityFromInfoMap(map);
		assertTrue(id.matches(read));
	}

	@Test
	void formatLibraryInUseMessageIncludesEntryCount() throws Exception {
		File fasta = new File(tempDir, "Marker4_highpep.fasta");
		try (PrintWriter w = new PrintWriter(fasta, StandardCharsets.UTF_8)) {
			w.println(">a");
			w.println("AAA");
			w.println(">b");
			w.println("BBB");
			w.println(">c");
			w.println("CCC");
		}
		String msg = Marker4LibrarySupport.formatLibraryInUseMessage(fasta);
		assertEquals("Marker4 library: Marker4_highpep.fasta, 3 entries", msg);
	}

	@Test
	void highpepRequiredFailsWithoutDeepDetect() throws Exception {
		File current = new File(tempDir, "mag");
		assertTrue(current.mkdirs());
		File marker4 = new File(current, Marker4ResourcePreparer.MARKER4_FASTA_NAME);
		try (PrintWriter w = new PrintWriter(marker4, StandardCharsets.UTF_8)) {
			// Write just over the threshold would be huge; instead call ensureHighpepIfNeeded
			// with a small file after temporarily... we unit-test the error by invoking with
			// a counted file that we stub via writing HIGHPEP_PROTEIN_THRESHOLD+1 headers is too heavy.
			// Cover the predictor-null path with a helper that already has protein count simulated
			// by writing a modest file and invoking ensureHighpepIfNeeded only when count is forced.
			w.println(">only_one");
			w.println("MPEPTIDE");
		}
		// Below threshold: should succeed with null predictor
		Marker4ResourcePreparer.ensureHighpepIfNeeded(marker4, current, null, null);

		// Create a fake highpep requirement by writing many headers is expensive; instead verify
		// the public constants and that missing highpep with null predictor throws when count is high.
		File large = new File(current, "large.fasta");
		try (PrintWriter w = new PrintWriter(large, StandardCharsets.UTF_8)) {
			int n = Marker4ResourcePreparer.HIGHPEP_PROTEIN_THRESHOLD + 1;
			for (int i = 0; i < n; i++) {
				w.println(">" + i);
				w.println("A");
			}
		}
		boolean threw = false;
		try {
			Marker4ResourcePreparer.ensureHighpepIfNeeded(large, current, null, null);
		} catch (Exception e) {
			threw = true;
			assertTrue(e.getMessage().contains("DeepDetect"));
			assertTrue(e.getMessage().contains("Refusing"));
		}
		assertTrue(threw, "expected failure when DeepDetect missing for large catalogue");
		Files.deleteIfExists(large.toPath());
	}

	@Test
	void countFastaProteins() throws Exception {
		File fasta = new File(tempDir, "t.fasta");
		try (PrintWriter w = new PrintWriter(fasta, StandardCharsets.UTF_8)) {
			w.println(">a");
			w.println("AA");
			w.println(">b desc");
			w.println("BB");
		}
		assertEquals(2, Marker4ResourcePreparer.countFastaProteins(fasta));
	}
}

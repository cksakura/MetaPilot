package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTimeoutPreemptively;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.time.Duration;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class FastaScannerTest {

	private static final int PROTEINS = 60_000;

	@TempDir
	Path tmp;

	@Test
	void everyProteinIsScannedWhenTheSinkSucceeds() throws IOException {
		File fasta = manyProteins();
		AtomicLong hits = new AtomicLong();
		long scanned = scanner(2).scan(List.of(fasta), (q, id, serial, offset, before, after) -> {
			assertEquals(3, offset, "MAA|PEPTIDEK");
			assertEquals('A', before);
			assertEquals('L', after);
			hits.incrementAndGet();
		}, null, "test");
		assertEquals(PROTEINS, scanned);
		assertEquals(PROTEINS, hits.get());
	}

	@Test
	void flanksAreDashesAtTheProteinEnds() throws IOException {
		File fasta = tmp.resolve("ends.fasta").toFile();
		try (PrintWriter w = new PrintWriter(fasta, StandardCharsets.UTF_8)) {
			w.println(">whole");
			w.println("PEPTIDEK");
		}
		StringBuilder seen = new StringBuilder();
		scanner(1).scan(List.of(fasta), (q, id, serial, offset, before, after) -> seen.append(offset).append(before)
				.append(after), null, "test");
		assertEquals("0--", seen.toString());
	}

	@Test
	void aFailingMatcherEndsTheScanInsteadOfBlockingTheReader() throws IOException {
		File fasta = manyProteins();
		for (int threads : new int[] { 1, 3 }) {
			IOException e = assertTimeoutPreemptively(Duration.ofSeconds(30),
					() -> assertThrows(IOException.class, () -> scanner(threads).scan(List.of(fasta),
							(q, id, serial, offset, before, after) -> {
								throw new IllegalStateException("sink failed");
							}, null, "test")));
			assertTrue(e.getCause() instanceof IllegalStateException, String.valueOf(e.getCause()));
		}
	}

	@Test
	void anErrorInAMatcherIsReportedAsAnIoException() throws IOException {
		File fasta = manyProteins();
		IOException e = assertTimeoutPreemptively(Duration.ofSeconds(30),
				() -> assertThrows(IOException.class, () -> scanner(1).scan(List.of(fasta), (q, id, serial, offset, before, after) -> {
					throw new OutOfMemoryError("simulated");
				}, null, "test")));
		assertTrue(e.getCause() instanceof OutOfMemoryError, String.valueOf(e.getCause()));
	}

	private static FastaScanner scanner(int threads) {
		return new FastaScanner(new PeptideAutomaton(Collections.singletonList("PEPTIDEK")), 1, false, threads, "");
	}

	private File manyProteins() throws IOException {
		File fasta = tmp.resolve("many.fasta").toFile();
		try (PrintWriter w = new PrintWriter(fasta, StandardCharsets.UTF_8)) {
			for (int i = 0; i < PROTEINS; i++) {
				w.println(">p" + i);
				w.println("MAAPEPTIDEKLLSGR");
			}
		}
		return fasta;
	}
}

package uk.ac.quadram.metapilot.core.dia;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;

class DiaHostLibrarySupportTest {

	@TempDir
	File dir;

	private static DiaNNPrecursor precursor(String seq, String... proteins) {
		return new DiaNNPrecursor("run", proteins, seq, seq, 2, 0, 0, 0.001, 0, 0.001, 0.001, 0.001, 0, 0);
	}

	private File touch(String name) throws IOException {
		File f = new File(dir, name);
		Files.write(f.toPath(), new byte[] { 1 });
		return f;
	}

	@Test
	void hostTakesSharedPeptides() {
		DiaNNPrecursor microOnly = precursor("AAAAK", "MGYG000000001_00001");
		DiaNNPrecursor shared = precursor("SHAREDK", "MGYG000000001_00002");
		DiaNNPrecursor hostOnly = precursor("HUMANK", "P12345");
		DiaHostLibrarySupport.Split split = DiaHostLibrarySupport.split(
				new DiaNNPrecursor[] { microOnly, shared, hostOnly }, Set.of("SHAREDK", "HUMANK"));
		assertArrayEquals(new DiaNNPrecursor[] { microOnly }, split.microbial);
		assertArrayEquals(new DiaNNPrecursor[] { shared, hostOnly }, split.host);
	}

	@Test
	void withoutHostTheLibraryIsUnchangedAndNoRecordIsLeft() throws Exception {
		File combinedDir = new File(dir, "combined");
		File speclib = touch("combined.speclib");
		DiaHostLibrarySupport.writeHostRecord(combinedDir, List.of(touch("stale_host.parquet")));
		String[] libs = DiaHostLibrarySupport.finalSearchLibraries(speclib, (File) null, combinedDir, null, null, null,
				null);
		assertArrayEquals(new String[] { speclib.getAbsolutePath() }, libs);
		assertFalse(new File(combinedDir, DiaHostLibrarySupport.HOST_RECORD).exists());
		assertTrue(DiaHostLibrarySupport.readHostRecord(combinedDir).isEmpty());
	}

	@Test
	void withHostBothParquetsAreSearchedAndRecorded() throws Exception {
		File combinedDir = new File(dir, "combined");
		File micro = touch("seed.parquet");
		File host1 = touch("host_human.parquet");
		File host2 = touch("crap.parquet");
		String[] libs = DiaHostLibrarySupport.finalSearchLibraries(micro, List.of(host1, host2), combinedDir, null,
				null, null, null);
		assertArrayEquals(new String[] { micro.getAbsolutePath(), host1.getAbsolutePath(), host2.getAbsolutePath() },
				libs);
		assertEquals(List.of(host1.getAbsoluteFile(), host2.getAbsoluteFile()),
				DiaHostLibrarySupport.readHostRecord(combinedDir).stream().map(File::getAbsoluteFile).toList());
	}

	@Test
	void hostMustBeAParquetLibrary() throws Exception {
		File micro = touch("seed.parquet");
		File host = touch("host_human.speclib");
		assertThrows(IOException.class, () -> DiaHostLibrarySupport.finalSearchLibraries(micro, host,
				new File(dir, "combined"), null, null, null, null));
	}

	@Test
	void libraryFolderParquetIsPreferred() throws Exception {
		File speclib = touch("Elo_faecal.speclib");
		assertEquals(speclib, DiaHostLibrarySupport.parquetInLibraryFolder(speclib));
		File parquet = touch("Elo_faecal.parquet");
		assertEquals(parquet.getAbsoluteFile(), DiaHostLibrarySupport.parquetInLibraryFolder(speclib).getAbsoluteFile());
		assertEquals(parquet, DiaHostLibrarySupport.parquetInLibraryFolder(parquet));
	}

	@Test
	void hostReportsListHostPeptides() throws Exception {
		DiaHostLibrarySupport.Split split = DiaHostLibrarySupport.split(new DiaNNPrecursor[] {
				precursor("AAAAK", "MGYG000000001_00001"), precursor("HUMANK", "P12345"),
				precursor("HUMANK", "P12345", "Q99999") }, Set.of("HUMANK"));
		File peptides = new File(dir, DiaHostLibrarySupport.HOST_PEPTIDES_TSV);
		File summary = new File(dir, DiaHostLibrarySupport.HOST_SUMMARY_TSV);
		DiaHostLibrarySupport.writeHostReports(split, new File(dir, "missing.pr_matrix.tsv"), peptides, summary);
		List<String> rows = Files.readAllLines(peptides.toPath(), StandardCharsets.UTF_8);
		assertEquals(List.of("Sequence\tProteins", "HUMANK\tP12345;Q99999"), rows);
		assertTrue(Files.readAllLines(summary.toPath(), StandardCharsets.UTF_8)
				.contains("# host precursors (all runs):\t2"));
	}
}

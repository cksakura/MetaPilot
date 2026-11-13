package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;

class OpenSearchHostDatabaseTest {

	@TempDir
	Path dir;

	private File write(String name, String... lines) throws IOException {
		File f = dir.resolve(name).toFile();
		Files.write(f.toPath(), List.of(lines), StandardCharsets.UTF_8);
		return f;
	}

	@Test
	@DisplayName("No host: the genome targets are searched as they are")
	void withoutHostTheGenomeTargetsAreSearched() throws IOException {
		File genomes = write("s_refined_targets.fasta", ">MGYG000000001_00001", "PEPTIDEK");
		assertSame(genomes, OpenSearchHostDatabase.searchTargets(genomes, null));
	}

	@Test
	@DisplayName("Host targets follow the genome targets; host decoys are left out, and the result is reused")
	void appendsHostTargetsAndReusesTheResult() throws IOException {
		File genomes = write("s_refined_targets.fasta", ">MGYG000000001_00001 a", "PEPTIDEK", ">MGYG000000001_00002",
				"ANOTHERK");
		File host = write("human.fasta", ">sp|P02768|ALBU_HUMAN Albumin", "DAHKSEVAHR", "FK", "",
				">rev_sp|P02768|ALBU_HUMAN", "KFRHAVESKHAD", ">sp|P68871|HBB_HUMAN", "VHLTPEEK");
		long past = System.currentTimeMillis() - 60_000;
		genomes.setLastModified(past);
		host.setLastModified(past);
		File searched = OpenSearchHostDatabase.searchTargets(genomes, host);

		assertEquals("s_refined_targets_with_host.fasta", searched.getName());
		List<String> lines = Files.readAllLines(searched.toPath(), StandardCharsets.UTF_8);
		assertEquals(List.of(">MGYG000000001_00001 a", "PEPTIDEK", ">MGYG000000001_00002", "ANOTHERK",
				">sp|P02768|ALBU_HUMAN Albumin", "DAHKSEVAHR", "FK", ">sp|P68871|HBB_HUMAN", "VHLTPEEK"), lines);

		// Unchanged sources: the same file, not rewritten (a rewrite would reset its time).
		long marked = past + 30_000;
		searched.setLastModified(marked);
		File again = OpenSearchHostDatabase.searchTargets(genomes, host);
		assertEquals(searched, again);
		assertEquals(marked, again.lastModified());

		// One decoy per target for genomes and host alike.
		File td = dir.resolve("s.fasta").toFile();
		TargetDecoyFasta.ensure(searched, td, TargetDecoyFasta.DecoyMethod.REVERSE);
		TargetDecoyFasta.Counts counts = TargetDecoyFasta.count(td);
		assertEquals(4, counts.targets);
		assertEquals(4, counts.decoys);
	}

	@Test
	@DisplayName("Another host file rebuilds the combined targets")
	void anotherHostRebuilds() throws IOException {
		File genomes = write("s_refined_targets.fasta", ">MGYG000000001_00001", "PEPTIDEK");
		File first = write("first.fasta", ">sp|P1|A_HUMAN", "AAAAK");
		File second = write("second.fasta", ">sp|P2|B_HUMAN", "CCCCK", ">sp|P3|C_HUMAN", "DDDDK");
		OpenSearchHostDatabase.searchTargets(genomes, first);
		File searched = OpenSearchHostDatabase.searchTargets(genomes, second);
		assertTrue(Files.readString(searched.toPath()).contains("sp|P3|C_HUMAN"));
		assertFalse(Files.readString(searched.toPath()).contains("sp|P1|A_HUMAN"));
	}

	@Test
	@DisplayName("A run that asks for the host fails without it instead of searching the genomes alone")
	void missingHostFails() throws IOException {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		assertNull(OpenSearchHostDatabase.resolve(params));

		params.setAppendHostDb(true);
		params.setHostDB("");
		assertThrows(IOException.class, () -> OpenSearchHostDatabase.resolve(params));
		params.setHostDB(dir.resolve("absent.fasta").toString());
		assertThrows(IOException.class, () -> OpenSearchHostDatabase.resolve(params));

		File host = write("human.fasta", ">sp|P1|A_HUMAN", "AAAAK");
		params.setHostDB(host.getAbsolutePath());
		assertEquals(host.getAbsoluteFile(), OpenSearchHostDatabase.resolve(params).getAbsoluteFile());
	}

	@Test
	@DisplayName("The genome test's candidates leave the host out when given the MAG prefix")
	void genomesOfLeavesTheHostOut() throws IOException {
		File db = write("db.fasta", ">MGYG000000001_00001", "PEPTIDEK", ">rev_MGYG000000001_00001", "KEDITPEP",
				">MGYG000000002_00001", "AAAAK", ">sp|P1|A_HUMAN", "CCCCK", ">rev_sp|P1|A_HUMAN", "KCCCC");
		assertEquals(Set.of("MGYG000000001", "MGYG000000002"), OpenSearchGenomeEvidence.genomesOf(db, "rev_", "MGYG"));
		assertTrue(OpenSearchGenomeEvidence.genomesOf(db, "rev_", null).contains("sp|P1|A"));
		assertTrue(OpenSearchHostDatabase.isHostAccession("sp|P1|A_HUMAN", "MGYG"));
		assertFalse(OpenSearchHostDatabase.isHostAccession("MGYG000000001_00001", "MGYG"));
	}
}

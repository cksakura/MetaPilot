package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.BitSet;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class SpeciesMapTest {

	@TempDir
	File tmp;

	private File metadata() throws IOException {
		File f = new File(tmp, "genomes-all_metadata.tsv");
		Files.write(f.toPath(), ("Genome\tGenome_type\tSpecies_rep\n" + "MGYG000000001\tIsolate\tMGYG000000001\n"
				+ "MGYG000000002\tMAG\tMGYG000000001\n" + "MGYG000000003\tMAG\tMGYG000000003\n"
				+ "MGYG000001406.1\tMAG\tMGYG000000005\n").getBytes(StandardCharsets.UTF_8));
		return f;
	}

	@Test
	@DisplayName("Protein keys round-trip, keeping the protein number's width")
	void proteinKeyRoundTrip() {
		for (String id : new String[] { "MGYG000004906_00001", "MGYG000000001_1", "MGYG999999999_999999",
				"MGYG000290000_01234" }) {
			long key = SpeciesMap.proteinKey(id);
			assertTrue(key >= 0, id);
			assertEquals(id, SpeciesMap.proteinIdOf(key));
		}
		assertEquals(-1, SpeciesMap.proteinKey("GUT_GENOME000001_1"));
		assertEquals(-1, SpeciesMap.proteinKey("MGYG000004906_0000001"));
		assertEquals(-1, SpeciesMap.proteinKey("MGYG000004906"));
	}

	@Test
	@DisplayName("Species come from the cluster table when present, otherwise from the own genome")
	void clusterSpecies() throws IOException {
		File rep = new File(tmp, "rep.tsv.gz");
		try (OutputStream out = new GZIPOutputStream(Files.newOutputStream(rep.toPath()))) {
			// unsorted on purpose, with a repeated representative to be unioned
			out.write(("MGYG000000003_00002\tMGYG000000003,MGYG000000001\n"
					+ "MGYG000000001_00007\tMGYG000000001,MGYG000000003\n"
					+ "MGYG000000003_00002\tMGYG000000005\n").getBytes(StandardCharsets.US_ASCII));
		}
		SpeciesMap map = SpeciesMap.load(metadata(), rep, null);
		BitSet s = new BitSet();
		assertTrue(map.addSpeciesOf("MGYG000000003_00002", s));
		assertEquals(3, s.cardinality());

		s.clear();
		assertTrue(map.addSpeciesOf("MGYG000000002_00010", s));
		assertEquals(1, s.cardinality());
		assertEquals("MGYG000000001", map.speciesId(s.nextSetBit(0)));

		s.clear();
		assertTrue(map.addSpeciesOf("MGYG000001406_00001", s), "versioned genome registered unversioned too");
		assertEquals("MGYG000000005", map.speciesId(s.nextSetBit(0)));

		assertFalse(map.addSpeciesOf("MGYG000777777_00001", new BitSet()));
	}

	@Test
	@DisplayName("The prepare-side table merges repeats and writes sorted rows")
	void writeRepSpecies() throws IOException {
		SpeciesMap map = SpeciesMap.metadataOnly(metadata());
		BitSet a = new BitSet();
		a.set(map.speciesOfGenome("MGYG000000003"));
		a.set(map.speciesOfGenome("MGYG000000001"));
		BitSet b = new BitSet();
		b.set(map.speciesOfGenome("MGYG000001406"));
		b.set(map.speciesOfGenome("MGYG000000001"));
		map.addRepRow(SpeciesMap.proteinKey("MGYG000000003_00002"), a);
		map.addRepRow(SpeciesMap.proteinKey("MGYG000000001_00007"), a);
		map.addRepRow(SpeciesMap.proteinKey("MGYG000000003_00002"), b);
		map.sortAndMerge();
		assertEquals(2, map.repRowCount());
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		map.writeRepSpecies(out);
		assertEquals("MGYG000000001_00007\tMGYG000000001,MGYG000000003\n"
				+ "MGYG000000003_00002\tMGYG000000001,MGYG000000003,MGYG000000005\n",
				out.toString(StandardCharsets.US_ASCII));
	}
}

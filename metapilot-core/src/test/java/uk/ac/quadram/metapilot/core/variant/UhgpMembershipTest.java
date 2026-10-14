package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class UhgpMembershipTest {

	@TempDir
	File tmp;

	private static final String TABLE = "MGYG000004906_00001\tMGYG000004906_00001\r\n"
			+ "MGYG000004906_00001\tMGYG000069512_01816\r\n" + "MGYG000004906_00001\tMGYG000069512_00002\r\n"
			+ "MGYG000004906_00033\tMGYG000004906_00033\n" + "MGYG000004906_00033\tMGYG000002168_01195\n"
			+ "MGYG000000007_1\tMGYG000000007_1\n" + "not-a-protein\tMGYG000000008_00001\n"
			+ "MGYG000004906_00001\tMGYG000000010_00009";

	private static long key(String id) {
		return SpeciesMap.proteinKey(id);
	}

	@Test
	@DisplayName("Keeps only requested clusters; genomes distinct and sorted; CRLF, no final newline, repeats")
	void read() throws IOException {
		File plain = new File(tmp, UhgpMembership.FILE);
		Files.write(plain.toPath(), TABLE.getBytes(StandardCharsets.US_ASCII));
		File gz = new File(tmp, UhgpMembership.FILE_GZ);
		try (OutputStream out = new GZIPOutputStream(Files.newOutputStream(gz.toPath()))) {
			out.write(TABLE.getBytes(StandardCharsets.US_ASCII));
		}
		long[] wanted = { key("MGYG000004906_00001"), key("MGYG000000007_1"), key("MGYG000099999_00001"),
				key("MGYG000004906_00001") };
		for (File f : new File[] { plain, gz }) {
			UhgpMembership m = UhgpMembership.read(f, wanted, null);
			assertEquals(8, m.rows, f.getName());
			assertEquals(3, m.representatives.length);
			assertEquals(2, m.representativesFound);
			assertArrayEquals(new int[] { 10, 4906, 69512 }, m.genomesOf(key("MGYG000004906_00001")),
					"later row of the same cluster is merged");
			assertArrayEquals(new int[] { 7 }, m.genomesOf(key("MGYG000000007_1")));
			assertArrayEquals(new int[0], m.genomesOf(key("MGYG000099999_00001")));
			assertNull(m.genomesOf(key("MGYG000004906_00033")), "not requested");
		}
		assertEquals(0, UhgpMembership.read(plain, new long[0], null).rows, "nothing requested, nothing read");
	}

	@Test
	@DisplayName("Byte parsing matches SpeciesMap.proteinKey")
	void proteinKey() {
		for (String id : new String[] { "MGYG000004906_00001", "MGYG000000007_1", "MGYG000289999_123456",
				"MGYG00000490_00001", "MGYG000004906_", "XGYG000004906_00001", "MGYG000004906_0000001" }) {
			byte[] b = id.getBytes(StandardCharsets.US_ASCII);
			assertEquals(SpeciesMap.proteinKey(id), UhgpMembership.proteinKey(b, 0, b.length), id);
		}
		byte[] member = "MGYG000069512_01816".getBytes(StandardCharsets.US_ASCII);
		assertEquals(69512, UhgpMembership.genomeNumber(member, 0, member.length));
		assertEquals(69512, SpeciesMap.genomeNumberOfKey(SpeciesMap.proteinKey("MGYG000069512_01816")));
	}

	@Test
	@DisplayName("Located from the setting, beside the FASTA, or in the extracted uhgp-100/ folder")
	void locate() throws IOException {
		File catalogue = new File(tmp, "protein_catalogue");
		File extracted = new File(catalogue, "uhgp-100");
		extracted.mkdirs();
		assertNull(UhgpMembership.locate("", catalogue));

		File besideCatalogue = new File(tmp, "uhgp-100/" + UhgpMembership.FILE);
		besideCatalogue.getParentFile().mkdirs();
		Files.write(besideCatalogue.toPath(), TABLE.getBytes(StandardCharsets.US_ASCII));
		assertEquals(besideCatalogue, UhgpMembership.locate("", catalogue), "uhgp-100/ in the MAG version dir");

		File inFolder = new File(extracted, UhgpMembership.FILE);
		Files.write(inFolder.toPath(), TABLE.getBytes(StandardCharsets.US_ASCII));
		assertEquals(inFolder, UhgpMembership.locate("", catalogue));

		File beside = new File(catalogue, UhgpMembership.FILE_GZ);
		Files.write(beside.toPath(), new byte[] { 1 });
		assertEquals(beside, UhgpMembership.locate(null, catalogue), "beside the FASTA wins");

		assertEquals(inFolder, UhgpMembership.locate(inFolder.getPath(), catalogue));
		assertEquals(inFolder, UhgpMembership.locate(extracted.getPath(), catalogue));
		assertNull(UhgpMembership.locate(new File(tmp, "missing.tsv").getPath(), catalogue),
				"a configured path that does not exist is not replaced by the default");
	}
}

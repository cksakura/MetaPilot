package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PangenomeResourcePreparerTest {

	@TempDir
	File tmp;

	static final String TSV = "MGYG000000003_00002\tMGYG000000003_00002\n" + "MGYG000000003_00002\tMGYG000000001_00004\n"
			+ "MGYG000000001_00007\tMGYG000000001_00007\n" + "MGYG000000001_00007\tMGYG000000002_00001\n"
			+ "MGYG000000003_00002\tMGYG000001406_00003\n";
	static final String FAA = ">MGYG000000003_00002 desc\nKPEPTIDEK\n>MGYG000000001_00007\nMAAAAK\nGGGR\n";

	/** One 512-byte tar header; TarStream ignores the checksum, so it is left blank. */
	static byte[] header(String name, long size, char type, boolean base256) {
		byte[] h = new byte[512];
		byte[] n = name.getBytes(StandardCharsets.US_ASCII);
		System.arraycopy(n, 0, h, 0, Math.min(100, n.length));
		if (base256) {
			h[124] = (byte) 0x80;
			for (int i = 135; i > 124; i--) {
				h[i] = (byte) (size & 0xff);
				size >>>= 8;
			}
		} else {
			byte[] s = String.format("%011o", size).getBytes(StandardCharsets.US_ASCII);
			System.arraycopy(s, 0, h, 124, s.length);
		}
		h[156] = (byte) type;
		byte[] magic = "ustar\u000000".getBytes(StandardCharsets.US_ASCII);
		System.arraycopy(magic, 0, h, 257, magic.length);
		return h;
	}

	static void member(OutputStream out, byte[] header, byte[] data) throws IOException {
		out.write(header);
		out.write(data);
		int pad = (512 - data.length % 512) % 512;
		out.write(new byte[pad]);
	}

	static String pax(String key, String value) {
		String body = " " + key + "=" + value + "\n";
		int len = body.length() + 1;
		while (String.valueOf(len).length() + body.length() != len) {
			len = String.valueOf(len).length() + body.length();
		}
		return len + body;
	}

	/** The archive layout MGnify ships: the TSV first, then the FASTA (here via PAX path and a GNU long name). */
	File archive(File dir) throws IOException {
		File f = new File(dir, PangenomeResource.ARCHIVE);
		ByteArrayOutputStream tar = new ByteArrayOutputStream();
		member(tar, header("uhgp-100/", 0, '5', false), new byte[0]);
		byte[] longName = ("uhgp-100/" + "x".repeat(120) + "/../uhgp-100.tsv").getBytes(StandardCharsets.US_ASCII);
		member(tar, header("././@LongLink", longName.length, 'L', false), longName);
		byte[] tsv = TSV.getBytes(StandardCharsets.US_ASCII);
		member(tar, header("truncated", tsv.length, '0', true), tsv);
		byte[] paxBody = (pax("path", "uhgp-100/uhgp-100.faa")).getBytes(StandardCharsets.US_ASCII);
		member(tar, header("PaxHeaders/uhgp-100.faa", paxBody.length, 'x', false), paxBody);
		byte[] faa = FAA.getBytes(StandardCharsets.US_ASCII);
		member(tar, header("short.faa", faa.length, '0', false), faa);
		tar.write(new byte[1024]);
		try (OutputStream out = new GZIPOutputStream(Files.newOutputStream(f.toPath()))) {
			out.write(tar.toByteArray());
		}
		return f;
	}

	File metadata(File dir) throws IOException {
		File f = new File(dir, PangenomeResource.METADATA);
		Files.write(f.toPath(), ("Genome\tSpecies_rep\n" + "MGYG000000001\tMGYG000000001\n"
				+ "MGYG000000002\tMGYG000000001\n" + "MGYG000000003\tMGYG000000003\n"
				+ "MGYG000001406.1\tMGYG000000005\n").getBytes(StandardCharsets.UTF_8));
		return f;
	}

	static String gunzip(File f) throws IOException {
		try (InputStream in = new GZIPInputStream(Files.newInputStream(f.toPath()))) {
			return new String(in.readAllBytes(), StandardCharsets.US_ASCII);
		}
	}

	@Test
	@DisplayName("Tar reader: long names, PAX paths and base-256 sizes")
	void tarStream() throws IOException {
		File a = archive(tmp);
		try (InputStream in = new GZIPInputStream(Files.newInputStream(a.toPath()))) {
			TarStream tar = new TarStream(in);
			TarStream.Entry dir = tar.next();
			assertFalse(dir.isFile());
			TarStream.Entry tsv = tar.next();
			assertTrue(tsv.name.endsWith("/uhgp-100.tsv"), tsv.name);
			assertEquals(TSV.length(), tsv.size);
			TarStream.Entry faa = tar.next();
			assertEquals("uhgp-100/uhgp-100.faa", faa.name);
			assertEquals(FAA, new String(tar.content().readAllBytes(), StandardCharsets.US_ASCII));
			assertNull(tar.next());
		}
		byte[] h = header("x", 9_000_000_000L, '0', true);
		assertEquals(9_000_000_000L, TarStream.parseSize(h));
	}

	@Test
	@DisplayName("Prepare writes the FASTA unchanged, the merged multi-species table and a manifest")
	void prepare() throws IOException {
		File dir = new File(tmp, "protein_catalogue");
		dir.mkdirs();
		File archive = archive(dir);
		File metadata = metadata(tmp);

		PangenomeResource before = PangenomeResource.locate(archive.getPath(), null);
		assertNotNull(before);
		assertTrue(before.needsPrepare());

		JSONObject manifest = PangenomeResourcePreparer.prepare(archive, metadata, dir, null);
		assertEquals(2, manifest.getLong("faa_records"));
		assertEquals(1, manifest.getInt("rep_species_rows"));

		assertEquals(FAA, gunzip(new File(dir, PangenomeResource.FAA)));
		assertEquals("MGYG000000003_00002\tMGYG000000001,MGYG000000003,MGYG000000005\n",
				gunzip(new File(dir, PangenomeResource.REP_SPECIES)));
		assertFalse(new File(dir, PangenomeResource.FAA + ".part").exists());

		PangenomeResource after = PangenomeResource.locate("", tmp);
		assertNotNull(after);
		assertFalse(after.needsPrepare());
		assertNotNull(after.repSpeciesGz);
		assertEquals(metadata.getAbsoluteFile(), after.metadataTsv.getAbsoluteFile());
	}

	@Test
	@DisplayName("A FASTA whose size disagrees with the manifest is not used")
	void incompletePreparation() throws IOException {
		File dir = new File(tmp, "protein_catalogue");
		dir.mkdirs();
		Files.write(new File(dir, PangenomeResource.FAA).toPath(), new byte[10]);
		Files.write(new File(dir, PangenomeResource.MANIFEST).toPath(),
				"{\"out_faa_gz_bytes\": 99}".getBytes(StandardCharsets.UTF_8));
		assertNull(PangenomeResource.locate("", tmp));
		assertNull(PangenomeResource.locate(new File(tmp, "nowhere").getPath(), tmp));
	}
}

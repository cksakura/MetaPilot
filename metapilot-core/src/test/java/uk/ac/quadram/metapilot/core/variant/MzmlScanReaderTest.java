package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Base64;
import java.util.Set;
import java.util.zip.Deflater;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MzmlScanReaderTest {

	@TempDir
	File tmp;

	private static String encode(double[] values, boolean doublePrecision, boolean zlib) {
		ByteBuffer b = ByteBuffer.allocate(values.length * (doublePrecision ? 8 : 4)).order(ByteOrder.LITTLE_ENDIAN);
		for (double v : values) {
			if (doublePrecision) {
				b.putDouble(v);
			} else {
				b.putFloat((float) v);
			}
		}
		byte[] bytes = b.array();
		if (zlib) {
			Deflater d = new Deflater();
			d.setInput(bytes);
			d.finish();
			ByteArrayOutputStream out = new ByteArrayOutputStream();
			byte[] buf = new byte[256];
			while (!d.finished()) {
				out.write(buf, 0, d.deflate(buf));
			}
			d.end();
			bytes = out.toByteArray();
		}
		return Base64.getEncoder().encodeToString(bytes);
	}

	private static String array(String kind, double[] values, boolean doublePrecision, boolean zlib) {
		return "<binaryDataArray encodedLength=\"0\">"
				+ "<cvParam cvRef=\"MS\" accession=\"" + (doublePrecision ? "MS:1000523" : "MS:1000521") + "\"/>"
				+ "<cvParam cvRef=\"MS\" accession=\"" + (zlib ? "MS:1000574" : "MS:1000576") + "\"/>"
				+ "<cvParam cvRef=\"MS\" accession=\"" + kind + "\"/>"
				+ "<binary>" + encode(values, doublePrecision, zlib) + "</binary></binaryDataArray>";
	}

	private static String spectrum(int index, String id, int msLevel, double[] mz, double[] in, boolean dbl,
			boolean zlib) {
		return "<spectrum index=\"" + index + "\" id=\"" + id + "\" defaultArrayLength=\"" + mz.length + "\">"
				+ "<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"" + msLevel + "\"/>"
				+ "<precursorList count=\"1\"><precursor><selectedIonList count=\"1\"><selectedIon>"
				+ "<cvParam cvRef=\"MS\" accession=\"MS:1000744\" value=\"500.0\"/>"
				+ "</selectedIon></selectedIonList></precursor></precursorList>"
				+ "<binaryDataArrayList count=\"2\">" + array("MS:1000514", mz, dbl, zlib)
				+ array("MS:1000515", in, dbl, zlib) + "</binaryDataArrayList></spectrum>";
	}

	@Test
	@DisplayName("Only wanted MS2 scans; 64-bit plain and 32-bit zlib; peaks sorted; numpress skipped")
	void read() throws IOException {
		double[] mz = { 300.5, 200.25, 400.125 };
		double[] in = { 30, 20, 40 };
		String numpress = "<spectrum index=\"5\" id=\"scan=10\" defaultArrayLength=\"1\">"
				+ "<cvParam cvRef=\"MS\" accession=\"MS:1000511\" value=\"2\"/>"
				+ "<binaryDataArrayList count=\"1\"><binaryDataArray encodedLength=\"0\">"
				+ "<cvParam cvRef=\"MS\" accession=\"MS:1002312\"/><cvParam cvRef=\"MS\" accession=\"MS:1000514\"/>"
				+ "<binary>AAAA</binary></binaryDataArray></binaryDataArrayList></spectrum>";
		String xml = "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
				+ "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\"><mzML><run id=\"r\"><spectrumList count=\"6\">"
				+ spectrum(0, "controllerType=0 controllerNumber=1 scan=5", 1, mz, in, true, false)
				+ spectrum(1, "controllerType=0 controllerNumber=1 scan=6", 2, mz, in, true, false)
				+ spectrum(2, "controllerType=0 controllerNumber=1 scan=7", 2, mz, in, false, true)
				+ spectrum(3, "controllerType=0 controllerNumber=1 scan=8", 2, mz, in, true, false)
				+ spectrum(8, "sample=1 period=1 cycle=9 experiment=1", 2, mz, in, true, true)
				+ numpress
				+ "</spectrumList></run></mzML></indexedmzML>";
		File f = new File(tmp, "s.mzML");
		Files.write(f.toPath(), xml.getBytes(StandardCharsets.UTF_8));

		MzmlScanReader.Result r = MzmlScanReader.read(f, Set.of(5, 6, 7, 9, 10));
		assertFalse(r.scans.containsKey(5), "MS1 is never returned");
		assertFalse(r.scans.containsKey(8), "not asked for");
		assertEquals(Set.of(6, 7, 9), r.scans.keySet(), "an id without scan= falls back to index + 1");
		assertArrayEquals(new double[] { 200.25, 300.5, 400.125 }, r.scans.get(9).mz, 1e-9);
		assertArrayEquals(new double[] { 200.25, 300.5, 400.125 }, r.scans.get(6).mz, 1e-9);
		assertArrayEquals(new float[] { 20, 30, 40 }, r.scans.get(6).intensity);
		assertArrayEquals(new double[] { 200.25, 300.5, 400.125 }, r.scans.get(7).mz, 1e-3);
		assertEquals(1, r.unsupportedEncoding);
		assertEquals(30f, r.scans.get(6).maxIntensity(300.4, 300.6));
		assertEquals(0f, r.scans.get(6).maxIntensity(250, 260));
	}

	@Test
	@DisplayName("Scan number from scan=N, else index + 1")
	void scanOf() {
		assertEquals(12, MzmlScanReader.scanOf("controllerType=0 controllerNumber=1 scan=12", "0"));
		assertEquals(5, MzmlScanReader.scanOf("index=4", "4"));
		assertEquals(-1, MzmlScanReader.scanOf(null, null));
	}
}

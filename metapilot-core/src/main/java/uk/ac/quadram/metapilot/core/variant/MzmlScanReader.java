package uk.ac.quadram.metapilot.core.variant;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Base64;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

/**
 * Reads the peaks of selected MS2 scans from an mzML file in one streaming
 * pass, decoding only the spectra asked for.
 * <p>
 * The scan number is taken from {@code scan=N} in the spectrum id (as
 * msconvert and ThermoRawFileParser write it, and as MSFragger reports it),
 * falling back to {@code index + 1}. Binary arrays may be 32- or 64-bit,
 * uncompressed or zlib; a numpress-encoded spectrum is skipped and counted,
 * not guessed at.
 */
public final class MzmlScanReader {

	private static final Pattern SCAN = Pattern.compile("scan=(\\d+)");

	/** One centroided spectrum, sorted by m/z. */
	public static final class Peaks {
		public final double[] mz;
		public final float[] intensity;

		public Peaks(double[] mz, float[] intensity) {
			this.mz = mz;
			this.intensity = intensity;
		}

		/** Highest intensity within {@code [lo, hi]}, or 0. */
		public float maxIntensity(double lo, double hi) {
			int i = Arrays.binarySearch(mz, lo);
			if (i < 0) {
				i = -i - 1;
			}
			float best = 0f;
			for (; i < mz.length && mz[i] <= hi; i++) {
				best = Math.max(best, intensity[i]);
			}
			return best;
		}
	}

	public static final class Result {
		public final Map<Integer, Peaks> scans = new HashMap<>();
		public int unsupportedEncoding;
	}

	private MzmlScanReader() {
	}

	public static Result read(File mzml, Set<Integer> wanted) throws IOException {
		Result result = new Result();
		XMLInputFactory factory = XMLInputFactory.newInstance();
		factory.setProperty(XMLInputFactory.SUPPORT_DTD, false);
		factory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, false);
		try (InputStream in = Files.newInputStream(mzml.toPath())) {
			XMLStreamReader x = factory.createXMLStreamReader(in);
			try {
				Spectrum s = null;
				Array array = null;
				while (x.hasNext()) {
					int event = x.next();
					if (event == XMLStreamConstants.START_ELEMENT) {
						String name = x.getLocalName();
						if ("spectrum".equals(name)) {
							int scan = scanOf(x.getAttributeValue(null, "id"), x.getAttributeValue(null, "index"));
							s = wanted.contains(scan) ? new Spectrum(scan) : null;
						} else if (s != null && "binaryDataArray".equals(name)) {
							array = new Array();
						} else if (s != null && "cvParam".equals(name)) {
							String acc = x.getAttributeValue(null, "accession");
							if (array != null) {
								array.cv(acc);
							} else if ("MS:1000511".equals(acc)) {
								s.msLevel = parseInt(x.getAttributeValue(null, "value"));
							}
						} else if (s != null && array != null && "binary".equals(name)) {
							array.base64 = x.getElementText().trim();
						}
					} else if (event == XMLStreamConstants.END_ELEMENT) {
						String name = x.getLocalName();
						if (s != null && "binaryDataArray".equals(name) && array != null) {
							s.add(array);
							array = null;
						} else if (s != null && "spectrum".equals(name)) {
							if (s.unsupported) {
								result.unsupportedEncoding++;
							} else if (s.msLevel != 1) {
								Peaks p = s.peaks();
								if (p != null) {
									result.scans.put(s.scan, p);
								}
							}
							s = null;
						}
					}
				}
			} finally {
				x.close();
			}
		} catch (XMLStreamException e) {
			throw new IOException("Cannot read " + mzml + ": " + e.getMessage(), e);
		}
		return result;
	}

	static int scanOf(String id, String index) {
		if (id != null) {
			Matcher m = SCAN.matcher(id);
			if (m.find()) {
				return Integer.parseInt(m.group(1));
			}
		}
		int i = parseInt(index);
		return i < 0 ? -1 : i + 1;
	}

	private static int parseInt(String s) {
		try {
			return s == null ? -1 : Integer.parseInt(s.trim());
		} catch (NumberFormatException e) {
			return -1;
		}
	}

	private static final class Array {
		boolean mz;
		boolean intensity;
		boolean doublePrecision;
		boolean zlib;
		boolean unsupported;
		String base64 = "";

		void cv(String acc) {
			if (acc == null) {
				return;
			}
			switch (acc) {
			case "MS:1000514":
				mz = true;
				break;
			case "MS:1000515":
				intensity = true;
				break;
			case "MS:1000523":
				doublePrecision = true;
				break;
			case "MS:1000521":
				doublePrecision = false;
				break;
			case "MS:1000574":
				zlib = true;
				break;
			case "MS:1000576":
				break;
			default:
				// MS-Numpress and other encodings
				if (acc.startsWith("MS:100231") || acc.startsWith("MS:100274")) {
					unsupported = true;
				}
			}
		}

		double[] decode() throws IOException {
			byte[] bytes = Base64.getDecoder().decode(base64);
			if (zlib) {
				bytes = inflate(bytes);
			}
			ByteBuffer b = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
			int n = bytes.length / (doublePrecision ? 8 : 4);
			double[] out = new double[n];
			for (int i = 0; i < n; i++) {
				out[i] = doublePrecision ? b.getDouble() : b.getFloat();
			}
			return out;
		}
	}

	private static byte[] inflate(byte[] in) throws IOException {
		Inflater inflater = new Inflater();
		inflater.setInput(in);
		ByteArrayOutputStream out = new ByteArrayOutputStream(in.length * 4);
		byte[] buf = new byte[1 << 14];
		try {
			while (!inflater.finished()) {
				int n = inflater.inflate(buf);
				if (n == 0 && (inflater.needsInput() || inflater.needsDictionary())) {
					break;
				}
				out.write(buf, 0, n);
			}
		} catch (DataFormatException e) {
			throw new IOException("Corrupt zlib array in mzML", e);
		} finally {
			inflater.end();
		}
		return out.toByteArray();
	}

	private static final class Spectrum {
		final int scan;
		int msLevel = -1;
		boolean unsupported;
		double[] mz;
		double[] intensity;

		Spectrum(int scan) {
			this.scan = scan;
		}

		void add(Array a) throws IOException {
			if (a.unsupported) {
				unsupported = true;
				return;
			}
			if (a.mz) {
				mz = a.decode();
			} else if (a.intensity) {
				intensity = a.decode();
			}
		}

		Peaks peaks() {
			if (mz == null || intensity == null || mz.length != intensity.length) {
				return null;
			}
			Integer[] order = new Integer[mz.length];
			boolean sorted = true;
			for (int i = 0; i < mz.length; i++) {
				order[i] = i;
				if (i > 0 && mz[i] < mz[i - 1]) {
					sorted = false;
				}
			}
			if (!sorted) {
				Arrays.sort(order, (a, b) -> Double.compare(mz[a], mz[b]));
			}
			double[] m = new double[mz.length];
			float[] in = new float[mz.length];
			for (int i = 0; i < mz.length; i++) {
				m[i] = mz[order[i]];
				in[i] = (float) intensity[order[i]];
			}
			return new Peaks(m, in);
		}
	}
}

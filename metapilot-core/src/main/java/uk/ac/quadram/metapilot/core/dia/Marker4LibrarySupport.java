package uk.ac.quadram.metapilot.core.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.Marker4ResourcePreparer;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Shared helper for maintaining Marker4 library metadata under {@code libraries/Marker4/}.
 *
 * <p>{@code Marker4_info.txt} records peptide count plus build identity (source FASTA name and
 * checksum, modification flags, DIA-NN version) so stale speclibs can be detected and rebuilt.
 */
public final class Marker4LibrarySupport {

	private static final String ENZYME = "Trypsin/P";
	private static final int MISS = 1;
	private static final boolean MATCH_TO_MAG = true;

	public static final String KEY_DATA_SET_SIZE = "Data set size";
	public static final String KEY_ENZYME = "Enzyme";
	public static final String KEY_MISS = "Miss cleavages";
	public static final String KEY_MATCH_TO_MAG = "MatchToMag";
	public static final String KEY_SOURCE_FASTA = "Source FASTA";
	public static final String KEY_SOURCE_FASTA_SHA256 = "Source FASTA SHA-256";
	public static final String KEY_MODIFICATIONS = "Modifications";
	public static final String KEY_DIANN_VERSION = "DIA-NN version";

	@FunctionalInterface
	public interface SpeclibToParquetConverter {
		boolean convert(DiannParameter diannPar, File speclib, File parquet, File libFolder,
				ProgressCallback callback, Logger logger) throws Exception;
	}

	private Marker4LibrarySupport() {
	}

	/**
	 * Identity of the FASTA + modifications + DIA-NN used to build Marker4.speclib.
	 */
	public static final class LibraryIdentity {
		public final String sourceFastaName;
		public final String sourceFastaSha256;
		public final String modifications;
		public final String diannVersion;

		public LibraryIdentity(String sourceFastaName, String sourceFastaSha256, String modifications,
				String diannVersion) {
			this.sourceFastaName = sourceFastaName == null ? "" : sourceFastaName.trim();
			this.sourceFastaSha256 = sourceFastaSha256 == null ? "" : sourceFastaSha256.trim().toLowerCase(Locale.ROOT);
			this.modifications = modifications == null ? "" : modifications.trim();
			this.diannVersion = diannVersion == null ? "" : diannVersion.trim();
		}

		public boolean matches(LibraryIdentity other) {
			if (other == null) {
				return false;
			}
			return sourceFastaName.equalsIgnoreCase(other.sourceFastaName)
					&& sourceFastaSha256.equals(other.sourceFastaSha256)
					&& modifications.equals(other.modifications)
					&& diannVersion.equals(other.diannVersion);
		}

		public String describeMismatch(LibraryIdentity expected) {
			if (expected == null) {
				return "expected identity is null";
			}
			StringBuilder sb = new StringBuilder();
			if (!sourceFastaName.equalsIgnoreCase(expected.sourceFastaName)) {
				sb.append("Source FASTA: cached=").append(sourceFastaName).append(" current=")
						.append(expected.sourceFastaName).append("; ");
			}
			if (!sourceFastaSha256.equals(expected.sourceFastaSha256)) {
				sb.append("Source FASTA SHA-256 differs; ");
			}
			if (!modifications.equals(expected.modifications)) {
				sb.append("Modifications: cached=").append(modifications).append(" current=")
						.append(expected.modifications).append("; ");
			}
			if (!diannVersion.equals(expected.diannVersion)) {
				sb.append("DIA-NN version: cached=").append(diannVersion).append(" current=")
						.append(expected.diannVersion).append("; ");
			}
			return sb.length() == 0 ? "identity mismatch" : sb.toString().trim();
		}
	}

	private static DiannParameter buildDiannParForLibMetadata(MetaParameterDia diaPar, boolean speclibOnly) {
		DiannParameter diannPar = new DiannParameter();
		diannPar.setThreads(diaPar.getThreadCount());
		diannPar.setEnzyme(ENZYME);
		diannPar.setMissed_cleavages(MISS);
		diannPar.setM_excision(diaPar.isMExcision());
		diannPar.setCarbami_C(diaPar.isCarbamiC());
		diannPar.setOxidation_M(diaPar.isOxidationM());
		diannPar.setAcety_proN(diaPar.isAcetyProN());
		diannPar.setPhospho_STY(diaPar.isPhosphoSTY());
		diannPar.setK_GG(diaPar.isKGG());
		if (speclibOnly) {
			diannPar.setMatrices(false);
			diannPar.setReanalyse(true);
		}
		return diannPar;
	}

	public static String formatModifications(MetaParameterDia diaPar) {
		if (diaPar == null) {
			return "";
		}
		return "mExcision=" + diaPar.isMExcision()
				+ ";carbamiC=" + diaPar.isCarbamiC()
				+ ";oxidationM=" + diaPar.isOxidationM()
				+ ";acetyProN=" + diaPar.isAcetyProN()
				+ ";phosphoSTY=" + diaPar.isPhosphoSTY()
				+ ";kGG=" + diaPar.isKGG();
	}

	public static String sha256Hex(File file) throws IOException {
		if (file == null || !file.exists()) {
			throw new IOException("Cannot checksum missing file: " + (file == null ? "null" : file.getAbsolutePath()));
		}
		MessageDigest md;
		try {
			md = MessageDigest.getInstance("SHA-256");
		} catch (NoSuchAlgorithmException e) {
			throw new IOException("SHA-256 not available", e);
		}
		byte[] buf = new byte[1 << 20];
		try (FileInputStream in = new FileInputStream(file)) {
			int n;
			while ((n = in.read(buf)) > 0) {
				md.update(buf, 0, n);
			}
		}
		StringBuilder hex = new StringBuilder(64);
		for (byte b : md.digest()) {
			hex.append(String.format(Locale.ROOT, "%02x", b));
		}
		return hex.toString();
	}

	/**
	 * Best-effort DIA-NN version label: parent folder if version-like, else {@code --version}
	 * / help probe, else executable name.
	 */
	public static String resolveDiannVersion(File diannExe) {
		if (diannExe == null) {
			return "unknown";
		}
		File parent = diannExe.getParentFile();
		if (parent != null) {
			String parentName = parent.getName().trim();
			if (parentName.matches("\\d+(\\.\\d+)+([.-]\\w+)?")) {
				return parentName;
			}
		}
		String probed = probeDiannVersion(diannExe);
		if (probed != null && !probed.isBlank()) {
			return probed;
		}
		return diannExe.getName();
	}

	private static String probeDiannVersion(File diannExe) {
		try {
			ProcessBuilder pb = new ProcessBuilder(diannExe.getAbsolutePath(), "--version");
			pb.redirectErrorStream(true);
			Process p = pb.start();
			StringBuilder out = new StringBuilder();
			try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
				String line;
				while ((line = br.readLine()) != null) {
					out.append(line).append('\n');
					if (out.length() > 4000) {
						break;
					}
				}
			}
			if (!p.waitFor(8, TimeUnit.SECONDS)) {
				p.destroyForcibly();
				return null;
			}
			String text = out.toString();
			for (String line : text.split("\\R")) {
				String t = line.trim();
				if (t.isEmpty()) {
					continue;
				}
				// e.g. "DIA-NN 1.8.1" or "Version: 2.0"
				java.util.regex.Matcher m = java.util.regex.Pattern
						.compile("(?i)(?:dia-?nn\\s*)?(?:version[:\\s]*)?(\\d+(?:\\.\\d+)+)")
						.matcher(t);
				if (m.find()) {
					return m.group(1);
				}
			}
		} catch (Exception ignored) {
			// best-effort only
		}
		return null;
	}

	public static LibraryIdentity buildExpectedIdentity(File sourceFasta, MetaParameterDia diaPar, File diannExe)
			throws IOException {
		Objects.requireNonNull(sourceFasta, "sourceFasta");
		return new LibraryIdentity(
				sourceFasta.getName(),
				sha256Hex(sourceFasta),
				formatModifications(diaPar),
				resolveDiannVersion(diannExe));
	}

	public static Map<String, String> readInfoMap(File infoFile) throws IOException {
		Map<String, String> map = new LinkedHashMap<>();
		if (infoFile == null || !infoFile.exists()) {
			return map;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(infoFile))) {
			String line;
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t", 2);
				if (cs.length < 2) {
					continue;
				}
				String key = cs[0].trim().replaceAll(":$", "");
				map.put(key, cs[1].trim());
			}
		}
		return map;
	}

	public static LibraryIdentity identityFromInfoMap(Map<String, String> map) {
		if (map == null || map.isEmpty()) {
			return null;
		}
		String fasta = first(map, KEY_SOURCE_FASTA);
		String sha = first(map, KEY_SOURCE_FASTA_SHA256);
		String mods = first(map, KEY_MODIFICATIONS);
		String ver = first(map, KEY_DIANN_VERSION);
		if ((fasta == null || fasta.isBlank()) && (sha == null || sha.isBlank())
				&& (mods == null || mods.isBlank()) && (ver == null || ver.isBlank())) {
			return null; // legacy info file without identity fields
		}
		return new LibraryIdentity(fasta, sha, mods, ver);
	}

	private static String first(Map<String, String> map, String key) {
		if (map.containsKey(key)) {
			return map.get(key);
		}
		for (Map.Entry<String, String> e : map.entrySet()) {
			if (e.getKey() != null && e.getKey().equalsIgnoreCase(key)) {
				return e.getValue();
			}
		}
		return null;
	}

	/**
	 * If a cached speclib exists but its {@code Marker4_info.txt} identity does not match the
	 * current job, delete speclib/parquet/info so they are rebuilt. Legacy info files without
	 * identity fields trigger a warning and rebuild.
	 *
	 * @return true if the speclib was invalidated and must be regenerated
	 */
	public static boolean invalidateSpeclibIfIdentityMismatch(MagDbItem magDb, File sourceFasta,
			MetaParameterDia diaPar, File diannExe, ProgressCallback callback, Logger logger) throws IOException {
		ProgressCallback cb = callback == null ? new ProgressCallback.NoOp() : callback;
		if (magDb == null) {
			return false;
		}
		File speclib = magDb.getMarker4LibFile();
		if (speclib == null || !speclib.exists()) {
			return false;
		}
		File libFolder = magDb.getMarker4LibraryFolder();
		File infoFile = libFolder == null ? null : new File(libFolder, MagDbItem.MARKER4 + "_info.txt");
		LibraryIdentity expected = buildExpectedIdentity(sourceFasta, diaPar, diannExe);

		LibraryIdentity cached = null;
		if (infoFile != null && infoFile.exists()) {
			cached = identityFromInfoMap(readInfoMap(infoFile));
		}

		if (cached != null && cached.matches(expected)) {
			return false;
		}

		String reason;
		if (cached == null) {
			reason = "Marker4_info.txt missing identity fields (Source FASTA / SHA-256 / Modifications / DIA-NN version)";
		} else {
			reason = cached.describeMismatch(expected);
		}
		String msg = "Marker4.speclib cache identity mismatch — rebuilding. " + reason;
		cb.onMessage("WARNING: " + msg);
		if (logger != null) {
			logger.warn(msg);
		}

		deleteQuietly(speclib, logger);
		if (libFolder != null) {
			deleteQuietly(new File(libFolder, MagDbItem.MARKER4 + ".parquet"), logger);
			deleteQuietly(new File(libFolder, MagDbItem.MARKER4 + "_info.txt"), logger);
			deleteQuietly(new File(libFolder, MagDbItem.MARKER4 + ".predicted.speclib"), logger);
		}
		return true;
	}

	private static void deleteQuietly(File file, Logger logger) {
		if (file == null || !file.exists()) {
			return;
		}
		try {
			Files.deleteIfExists(file.toPath());
		} catch (IOException e) {
			if (logger != null) {
				logger.warn("Could not delete stale Marker4 file {}: {}", file.getAbsolutePath(), e.getMessage());
			}
		}
	}

	public static void writeInfoFile(File infoFile, int peptideCount, LibraryIdentity identity) throws IOException {
		if (infoFile == null) {
			throw new IOException("Marker4 info file path is null");
		}
		File parent = infoFile.getParentFile();
		if (parent != null && !parent.exists() && !parent.mkdirs()) {
			throw new IOException("Could not create Marker4 library folder: " + parent.getAbsolutePath());
		}
		try (PrintWriter writer = new PrintWriter(infoFile, StandardCharsets.UTF_8)) {
			writer.println(KEY_DATA_SET_SIZE + ":\t" + peptideCount);
			writer.println(KEY_ENZYME + ":\t" + ENZYME);
			writer.println(KEY_MISS + ":\t" + MISS);
			writer.println(KEY_MATCH_TO_MAG + ":\t" + MATCH_TO_MAG);
			if (identity != null) {
				writer.println(KEY_SOURCE_FASTA + ":\t" + identity.sourceFastaName);
				writer.println(KEY_SOURCE_FASTA_SHA256 + ":\t" + identity.sourceFastaSha256);
				writer.println(KEY_MODIFICATIONS + ":\t" + identity.modifications);
				writer.println(KEY_DIANN_VERSION + ":\t" + identity.diannVersion);
			}
		}
	}

	/**
	 * Log the Marker4 library FASTA in use (name + entry count) at search / build start.
	 */
	public static String formatLibraryInUseMessage(File fastaFile) {
		if (fastaFile == null) {
			return "Marker4 library: (null FASTA)";
		}
		int entries = Marker4ResourcePreparer.countFastaProteins(fastaFile);
		String name = fastaFile.getName();
		if (entries < 0) {
			return "Marker4 library: " + name + " (entry count unavailable)";
		}
		return String.format(Locale.ROOT, "Marker4 library: %s, %,d entries", name, entries);
	}

	public static void logLibraryInUse(File fastaFile, ProgressCallback callback, Logger logger) {
		String msg = formatLibraryInUseMessage(fastaFile);
		if (callback != null) {
			callback.onMessage(msg);
		}
		if (logger != null) {
			logger.info(msg);
		}
	}

	/**
	 * Ensures {@code Marker4.parquet}, {@code Marker4_info.txt}, and in-memory registration exist for the canonical
	 * speclib (parquet + info + {@link MagDbItem#addPepLib}).
	 */
	public static boolean ensureMarker4LibraryMetadata(MagDbItem magDb, MetaParameterDia diaPar,
			ProgressCallback callback, Logger logger, boolean speclibOnly,
			SpeclibToParquetConverter converter) throws Exception {
		return ensureMarker4LibraryMetadata(magDb, diaPar, null, null, callback, logger, speclibOnly, converter);
	}

	/**
	 * Same as {@link #ensureMarker4LibraryMetadata(MagDbItem, MetaParameterDia, ProgressCallback, Logger, boolean, SpeclibToParquetConverter)}
	 * but records source FASTA / DIA-NN identity in {@code Marker4_info.txt}.
	 */
	public static boolean ensureMarker4LibraryMetadata(MagDbItem magDb, MetaParameterDia diaPar,
			File sourceFasta, File diannExe, ProgressCallback callback, Logger logger, boolean speclibOnly,
			SpeclibToParquetConverter converter) throws Exception {
		ProgressCallback cb = callback == null ? new ProgressCallback.NoOp() : callback;
		if (magDb == null) {
			if (logger != null) logger.error("Marker4 metadata: magDb is null");
			return false;
		}
		if (diaPar == null) {
			if (logger != null) logger.error("Marker4 metadata: DIA parameters are null");
			return false;
		}
		if (converter == null) {
			if (logger != null) logger.error("Marker4 metadata: speclib converter is null");
			return false;
		}

		File speclib = magDb.getMarker4LibFile();
		if (speclib == null || !speclib.exists()) {
			if (logger != null) logger.error("Marker4 metadata: canonical speclib missing");
			return false;
		}
		File libFolder = magDb.getMarker4LibraryFolder();
		if (libFolder == null) {
			if (logger != null) logger.error("Marker4 metadata: library folder null");
			return false;
		}
		if (!magDb.ensureMarker4LibraryFolderExists()) {
			if (logger != null) logger.error("Marker4 metadata: could not create library folder {}", libFolder.getAbsolutePath());
			return false;
		}

		String name = MagDbItem.MARKER4;
		File parquet = new File(libFolder, name + ".parquet");
		File infoFile = new File(libFolder, name + "_info.txt");

		DiannParameter diannPar = buildDiannParForLibMetadata(diaPar, speclibOnly);
		if (!parquet.exists()) {
			cb.onMessage("Marker4: converting speclib to parquet (" + parquet.getName() + ")");
			converter.convert(diannPar, speclib, parquet, libFolder, cb, logger);
		}
		if (!parquet.exists()) {
			if (logger != null) logger.error("Marker4 metadata: parquet conversion failed, expected {}", parquet.getAbsolutePath());
			return false;
		}

		HashSet<String> pepSet = DiaNNParquetReader.getPeptideSet(parquet.getAbsolutePath());
		int n = pepSet != null ? pepSet.size() : 0;

		File fasta = sourceFasta != null ? sourceFasta : magDb.getMarker4FastaForLibrary();
		LibraryIdentity identity = null;
		if (fasta != null && fasta.exists()) {
			identity = buildExpectedIdentity(fasta, diaPar, diannExe);
		}

		writeInfoFile(infoFile, n, identity);
		magDb.addPepLib(name, n, ENZYME, MISS, MATCH_TO_MAG);
		return true;
	}

	public static boolean ensureMarker4LibraryMetadata(File diannExe, File condaPrefix, MagDbItem magDb,
			MetaParameterDia diaPar, ProgressCallback callback, Logger logger) throws Exception {
		File sourceFasta = magDb != null ? magDb.getMarker4FastaForLibrary() : null;
		return ensureMarker4LibraryMetadata(magDb, diaPar, sourceFasta, diannExe, callback, logger, true,
				(diannPar, speclib, parquet, libFolder, cb, log) -> {
					String cmd = diannPar.generateSpeclib2Parquet(speclib.getAbsolutePath(),
							parquet.getAbsolutePath());
					int exit = DiaNNHeadlessRunner.run(diannExe, cmd, libFolder, cb, condaPrefix);
					if (exit != 0 && log != null) {
						log.error("Marker4 metadata: DIA-NN speclib->parquet exit code {}", exit);
					}
					return exit == 0;
				});
	}
}

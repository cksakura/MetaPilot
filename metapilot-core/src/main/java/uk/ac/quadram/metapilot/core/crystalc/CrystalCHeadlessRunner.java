package uk.ac.quadram.metapilot.core.crystalc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

/**
 * OS-agnostic runner for Crystal-C, the open-search precursor-refinement step.
 * <p>
 * Mirrors {@code MSBoosterHeadlessRunner}: one process invocation returning an
 * exit code, living in {@code metapilot-core} so the GUI, CLI and HPC entry
 * points share it. Crystal-C is Apache-2.0, the same licence as MetaPilot;
 * cite Chang HY et al., "Crystal-C: A Computational Tool for Refinement of Open
 * Search Results", J Proteome Res 2020.
 *
 * <h2>Why this step exists</h2>
 * An open search reports {@code massdiff = precursor_neutral_mass -
 * calc_neutral_pep_mass}. When the instrument selected the wrong peak of the
 * precursor's isotope cluster, or a co-eluting neighbour inside the isolation
 * window, that difference is not a modification - it is a precursor-assignment
 * error, and it lands a few Da away from zero. Crystal-C re-reads the MS1 scan,
 * re-picks the monoisotopic peak and rewrites {@code precursor_neutral_mass}.
 * <p>
 * Measured on 12mix (2026-09-14), against FragPipe's own open workflow, which
 * runs this step and which MetaPilot's open search did not:
 * <ul>
 * <li>Both searches find the same number of rank-1 unmodified target PSMs in
 * the raw pepXML - 352,484 (MetaPilot) against 352,360 (FragPipe), pooled over
 * nine files. The searches are not the difference.</li>
 * <li>Crystal-C then takes FragPipe's zero-delta rank-1 target count from
 * 317,000 to 591,641. In {@code YS-194-1-1} alone it moved 29,485 of 106,340
 * rank-1 target scans to |delta| &lt; 0.02 Da.</li>
 * <li>92.8% of those had a pre-correction delta within +/-5 Da, piled on
 * integer Da values. Only 13% are within 0.02 Da of a clean n x 1.00336 13C
 * step, which is why MetaPilot's own isotope-family normalisation inside
 * {@code DeltaMassPeakFitter} recovers only ~45,000 of them: the rest are
 * co-isolation artifacts that arithmetic isotope subtraction cannot reach.</li>
 * <li>The cost of not running it: 316,120 of 785,180 validated target PSMs
 * (40.3%) sat at 0.05-5.5 Da, and 335 of the 500 fitted delta-mass peaks lay
 * inside +/-5.5 Da - artifacts being reported as modifications.</li>
 * </ul>
 *
 * <h2>What it does and does not change</h2>
 * Verified on a real before/after pair: 138,543 {@code spectrum_query} in,
 * 138,543 out, same {@code search_hit} count, {@code alternative_protein}
 * entries preserved. It rewrites {@code precursor_neutral_mass} and
 * {@code massdiff}, keeps {@code uncalibrated_precursor_neutral_mass} so the
 * original is recoverable, and rewrites the peptide only for its
 * missed-cleavage / semi-tryptic artifact class. Nothing is deleted, so this is
 * a refinement of the candidate list, not a filter on it - which is what lets
 * MetaPilot keep deciding acceptance itself.
 *
 * <h2>Ranks above 1</h2>
 * FragPipe feeds Crystal-C a {@code output_report_topN = 1} pepXML; MetaPilot's
 * round two reports up to five. Nothing here depends on which ranks Crystal-C
 * chose to rewrite, because the correction is a property of the <em>spectrum</em>:
 * {@code OpenSearchPepxmlReader} recomputes every hit's delta as
 * {@code precursor_neutral_mass - calc_neutral_pep_mass} rather than trusting
 * the {@code massdiff} attribute. That identity was verified to hold on both an
 * uncorrected MetaPilot pepXML (659,521 hits, max deviation 4.7e-4 Da, a float
 * rounding artefact) and a Crystal-C output (106,340 hits, max deviation
 * 1e-5 Da), so the recomputation is a no-op when this step is off.
 *
 * <h2>Parameters</h2>
 * Key names, value syntax and defaults are taken from the shipped tool, not
 * guessed: {@code crystalc.p_ReadData.ParseParameters} recognises exactly the
 * ten keys written by {@link #writeParams}, and the defaults are FragPipe's own
 * bundled {@code crystalc.params}. FragPipe overrides only the four paths; the
 * numeric settings below are its defaults, so a MetaPilot run and a FragPipe
 * run refine the same way.
 */
public final class CrystalCHeadlessRunner {

	private static final Logger LOGGER = LogManager.getLogger(CrystalCHeadlessRunner.class);

	/**
	 * Crystal-C's entry point. Needed explicitly because the classpath form of
	 * {@code java} does not read the jar's manifest - and this jar declares no
	 * {@code Main-Class} anyway.
	 */
	private static final String MAIN_CLASS = "crystalc.Run";

	/** Suffix Crystal-C appends to the input pepXML's stem for its output. */
	public static final String CORRECTED_SUFFIX = "_c.pepXML";

	/** Suffix of its per-PSM report, written beside the corrected pepXML. */
	public static final String REPORT_SUFFIX = "_c.txt";

	private CrystalCHeadlessRunner() {
	}

	/**
	 * Everything Crystal-C's parameter file can carry.
	 * <p>
	 * The four path fields are the ones FragPipe sets per run; the rest are
	 * FragPipe's bundled defaults and should be left alone without a measured
	 * reason to change them.
	 */
	public static final class Settings {
		/** Folder holding the mzML the pepXML was searched from. Required. */
		public File rawFileLocation;
		/** Extension of the spectrum files in that folder, without a dot. */
		public String rawFileExtension = "mzML";
		/** Where the {@code _c.pepXML} is written. Required. */
		public File outputLocation;
		/**
		 * The FASTA the search used. Required - Crystal-C refuses an empty value,
		 * and it needs the sequences for its missed-cleavage artifact class.
		 */
		public File fasta;
		/** {@code -1} means "cores - 1", which is Crystal-C's own reading of it. */
		public int threads = -1;
		/** Charge range searched for chimeric precursors. */
		public int minCharge = 1;
		public int maxCharge = 6;
		/** Theoretical isotope peaks considered per cluster. */
		public int isotopeNumber = 3;
		/** Precursor mass tolerance, in ppm. */
		public double precursorMassPpm = 20;
		/** Precursor isolation window, in Da. */
		public double precursorIsolationWindow = 0.7;
		/**
		 * Present for completeness and written so the file is self-describing.
		 * <p>
		 * It is dead in Crystal-C 1.5.10: {@code ds_Parameters.CorrectIsotopeError}
		 * is parsed and stored but read by no processing class, and the isotope
		 * correction happens regardless. FragPipe ships it {@code false} and still
		 * gets the corrections measured above. Do not read this field as a switch.
		 */
		public boolean correctIsotopeError = false;
	}

	/**
	 * Writes a Crystal-C parameter file for one pepXML.
	 * <p>
	 * Named after the pepXML rather than fixed, because a sample folder may hold
	 * more than one run's parameters and because FragPipe's own convention -
	 * {@code crystalc-<n>-<pepxml name>.params} - makes a failed run traceable to
	 * its input. Written before the run, not after: a failed run is exactly when
	 * the parameters need inspecting. FragPipe deletes these; MetaPilot keeps
	 * them, for the same provenance reason its MSFragger parameters are kept.
	 *
	 * <h3>Format</h3>
	 * {@code key = value}, one per line; {@code #} starts a comment anywhere in
	 * the line and blank lines are skipped. {@code precursor_charge} is the one
	 * value with internal structure: Crystal-C splits it on a single space and
	 * parses both halves as integers, so it must be exactly
	 * {@code <min> <max>} - two numbers, one space, no padding.
	 */
	public static File writeParams(File workDir, File pepXml, Settings settings) throws IOException {
		if (workDir == null) {
			throw new IOException("workDir is null");
		}
		if (pepXml == null) {
			throw new IOException("pepXML is null");
		}
		validate(settings);
		File spectra = expectedSpectrumFile(settings, pepXml);
		if (!spectra.isFile()) {
			throw new IOException("Crystal-C would look for " + spectra.getAbsolutePath() + " and it is not there. "
					+ "It derives the spectrum file from the pepXML's own name, truncated at its FIRST dot, so "
					+ "the pepXML must be named <spectrum file stem>.pepXML and the stem itself must contain no "
					+ "dot. Rename the sample, or point raw_file_location at the folder that holds its "
					+ settings.rawFileExtension + ".");
		}
		workDir.mkdirs();
		File out = new File(workDir, "crystalc-" + pepXml.getName() + ".params");
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Written by MetaPilot for Crystal-C 1.5.x. Key names and value syntax follow");
			w.println("# crystalc.p_ReadData.ParseParameters; the numeric defaults are FragPipe's own");
			w.println("# bundled crystalc.params, so this refines the same way FragPipe's open workflow does.");
			w.println();
			w.println("thread = " + settings.threads);
			w.println("fasta = " + settings.fasta.getAbsolutePath());
			w.println("raw_file_location = " + settings.rawFileLocation.getAbsolutePath());
			w.println("raw_file_extension = " + settings.rawFileExtension);
			w.println("output_location = " + settings.outputLocation.getAbsolutePath());
			w.println();
			// One space, two integers: Crystal-C does value.split(" ") and parses
			// both halves. Any other spacing throws NumberFormatException inside
			// its parser, which exits the JVM with status 1 and no useful message.
			w.println("precursor_charge = " + settings.minCharge + " " + settings.maxCharge);
			w.println("isotope_number = " + settings.isotopeNumber);
			w.println("precursor_mass = " + trim(settings.precursorMassPpm) + "   # ppm");
			w.println("precursor_isolation_window = " + trim(settings.precursorIsolationWindow));
			w.println("correct_isotope_error = " + settings.correctIsotopeError
					+ "   # parsed but unused in 1.5.10; correction happens regardless");
		}
		return out;
	}

	/**
	 * The spectrum file Crystal-C will open for {@code pepXml}, by its own rule.
	 * <p>
	 * {@code p_ReadData.LoadRawFile} builds
	 * {@code raw_file_location + separator + stem + "." + raw_file_extension},
	 * where {@code stem} is the pepXML's file name truncated at its
	 * <em>first</em> dot - not its last. A sample called {@code Run.2} therefore
	 * resolves to {@code Run.mzML} and the run dies on a missing file, so this is
	 * checked before launching rather than discovered afterwards.
	 */
	public static File expectedSpectrumFile(Settings settings, File pepXml) {
		String name = pepXml.getName();
		int dot = name.indexOf('.');
		String stem = dot >= 0 ? name.substring(0, dot) : name;
		return new File(settings.rawFileLocation, stem + "." + settings.rawFileExtension);
	}

	private static void validate(Settings settings) throws IOException {
		if (settings == null) {
			throw new IOException("Crystal-C settings are null");
		}
		if (settings.fasta == null || !settings.fasta.isFile()) {
			throw new IOException("Crystal-C needs the FASTA the search used; it refuses an empty fasta= value. Got: "
					+ settings.fasta);
		}
		if (settings.rawFileLocation == null || !settings.rawFileLocation.isDirectory()) {
			throw new IOException("Crystal-C needs the folder holding the spectrum files, and it must exist "
					+ "because it re-reads the MS1 scans from there. Got: " + settings.rawFileLocation);
		}
		if (settings.outputLocation == null) {
			throw new IOException("Crystal-C output location is null");
		}
		if (settings.rawFileExtension == null || settings.rawFileExtension.isBlank()) {
			throw new IOException("Crystal-C raw file extension is blank");
		}
		// LoadRawFile switches on exactly these three and returns null for anything
		// else, which surfaces later as a NullPointerException rather than as the
		// unsupported format it is.
		String ext = settings.rawFileExtension.trim();
		if (!ext.equalsIgnoreCase("mzML") && !ext.equalsIgnoreCase("mzXML") && !ext.equalsIgnoreCase("raw")) {
			throw new IOException("Crystal-C reads mzML, mzXML or Thermo raw only; got '" + ext
					+ "'. Convert the spectra, or switch this step off.");
		}
		if (settings.minCharge < 1 || settings.maxCharge < settings.minCharge) {
			throw new IOException("Crystal-C precursor charge range is invalid: " + settings.minCharge + " to "
					+ settings.maxCharge);
		}
	}

	/** Where Crystal-C writes the refined pepXML for {@code pepXml}. */
	public static File correctedPepXml(File outputLocation, File pepXml) {
		String name = pepXml.getName();
		String stem = name.toLowerCase(Locale.ROOT).endsWith(".pepxml")
				? name.substring(0, name.length() - ".pepXML".length())
				: name;
		return new File(outputLocation, stem + CORRECTED_SUFFIX);
	}

	/** Crystal-C's per-PSM report, beside the corrected pepXML. */
	public static File report(File outputLocation, File pepXml) {
		String name = pepXml.getName();
		String stem = name.toLowerCase(Locale.ROOT).endsWith(".pepxml")
				? name.substring(0, name.length() - ".pepXML".length())
				: name;
		return new File(outputLocation, stem + REPORT_SUFFIX);
	}

	/**
	 * Runs Crystal-C over one pepXML.
	 *
	 * @return the process exit code; 0 on success.
	 */
	public static int run(File crystalcJar, String javaMaxHeap, File paramsFile, File pepXml, ProgressCallback cb)
			throws IOException, InterruptedException {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (crystalcJar == null || !crystalcJar.isFile()) {
			throw new IOException("Crystal-C jar not found: "
					+ (crystalcJar == null ? "null" : crystalcJar.getAbsolutePath()));
		}
		if (paramsFile == null || !paramsFile.isFile()) {
			throw new IOException("Crystal-C params file not found: "
					+ (paramsFile == null ? "null" : paramsFile.getAbsolutePath()));
		}
		if (pepXml == null || !pepXml.isFile()) {
			throw new IOException("Crystal-C input pepXML not found: " + (pepXml == null ? "null" : pepXml));
		}
		if (!pepXml.getName().toLowerCase(Locale.ROOT).endsWith(".pepxml")) {
			throw new IOException("Crystal-C only accepts .pepXML input; got " + pepXml.getName());
		}

		// The same absolute-java rule MSBooster's runner documents: MetaPilot
		// launched from an IDE or a javaw shortcut hands its children a PATH with
		// no java on it, and on Windows that surfaces as exit 9009 - which looks
		// exactly like Crystal-C producing nothing.
		String javaCommand = MSFraggerHeadlessRunner.resolveJavaCommand();
		if ("java".equals(javaCommand) || !new File(javaCommand).isFile()) {
			throw new IOException("Could not resolve an absolute java executable from java.home ("
					+ System.getProperty("java.home") + "); the launch would fall through to the PATH, which on "
					+ "Windows surfaces as exit 9009 (command not found) and looks exactly like a Crystal-C "
					+ "failure. Refusing to start Crystal-C without a java it can actually run.");
		}

		List<String> cmd = new ArrayList<>();
		cmd.add(javaCommand);
		if (javaMaxHeap != null && !javaMaxHeap.isBlank()) {
			cmd.add("-Xmx" + javaMaxHeap);
		}
		cmd.add("-cp");
		String classpath = buildClasspath(crystalcJar);
		cmd.add(classpath);
		cmd.add(MAIN_CLASS);
		// Argument order is fixed by crystalc.Run.main: args[0] is the parameter
		// file and args[1..] are pepXML files. Passing one file per invocation
		// rather than the whole batch keeps a failure attributable to a sample.
		cmd.add(paramsFile.getAbsolutePath());
		cmd.add(pepXml.getAbsolutePath());

		cb.onMessage("Running Crystal-C: " + String.join(" ", cmd));

		ProcessBuilder pb = new ProcessBuilder(cmd);
		File workDir = paramsFile.getParentFile();
		if (workDir != null && workDir.exists()) {
			pb.directory(workDir);
		}
		pb.redirectErrorStream(true);

		Process p = pb.start();
		try (BufferedReader br = new BufferedReader(
				new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (Thread.currentThread().isInterrupted()) {
					p.destroyForcibly();
					throw new InterruptedException("Interrupted while running Crystal-C");
				}
				String trimmed = line.trim();
				if (!trimmed.isEmpty()) {
					cb.onMessage(trimmed);
				}
			}
		}
		int exit = p.waitFor();
		if (exit == 0) {
			LOGGER.info("Crystal-C refined {} -> {}", pepXml.getName(),
					correctedPepXml(new File(readOutputLocation(paramsFile)), pepXml).getName());
		}
		if (exit != 0) {
			// Crystal-C calls System.exit(1) from its own parameter parser, so a 1
			// here is as likely to be a malformed params file as a data problem -
			// which is why the params file is kept rather than deleted.
			cb.onMessage("Crystal-C exit code: " + exit + " (java: " + javaCommand + ")"
					+ (exit == 9009 ? " - 9009 is Windows command-not-found: that java executable does not exist"
							: "")
					+ ". Parameters: " + paramsFile.getAbsolutePath() + ". Classpath: " + classpath);
		}
		return exit;
	}

	/**
	 * The three jars FragPipe itself puts on Crystal-C's classpath, in its order.
	 * <p>
	 * Taken from a working FragPipe run's log and confirmed against
	 * {@code CmdCrystalc}, which names {@code batmass-io-1.36.5.jar} and
	 * {@code grppr-0.3.23.jar} literally:
	 * <pre>
	 * java.exe -Xmx108G
	 *     -cp &lt;tools&gt;\original-crystalc-1.5.10.jar;&lt;tools&gt;\batmass-io-1.36.5.jar;&lt;tools&gt;\grppr-0.3.23.jar
	 *     crystalc.Run &lt;params&gt; &lt;pepXML&gt;
	 * </pre>
	 * Neither is optional. Crystal-C reads the MS1 scans through batmass-io and
	 * builds its averagine isotope templates through grppr
	 * ({@code umich.ms.ims.isodistr}), so a missing one dies at
	 * NoClassDefFoundError with nothing written - indistinguishable from
	 * Crystal-C simply doing nothing. Resolved by prefix rather than by the
	 * literal versioned names so an upgraded installation keeps working, and
	 * named rather than a {@code tools/*} wildcard so a version collision cannot
	 * creep in.
	 */
	static String buildClasspath(File crystalcJar) throws IOException {
		List<String> entries = new ArrayList<>();
		entries.add(crystalcJar.getAbsolutePath());
		File toolsDir = crystalcJar.getParentFile();
		for (String prefix : new String[] { "batmass-io", "grppr" }) {
			File jar = toolsDir != null && toolsDir.isDirectory() ? highestVersioned(toolsDir, prefix) : null;
			if (jar == null) {
				throw new IOException("No " + prefix + " jar beside " + crystalcJar.getAbsolutePath()
						+ ". Crystal-C reads the MS1 scans through batmass-io and builds its isotope "
						+ "templates through grppr, so without it the run dies at NoClassDefFoundError "
						+ "having written nothing - which looks like Crystal-C doing nothing. Point "
						+ "Resources > FragPipe at a complete installation.");
			}
			entries.add(jar.getAbsolutePath());
		}
		return String.join(File.pathSeparator, entries);
	}

	/**
	 * Reads {@code output_location} back out of the parameter file, so the log
	 * line names the file Crystal-C actually wrote rather than one this class
	 * assumed. Falls back to the parameter file's own folder, which is where
	 * {@code writeParams} puts them.
	 */
	private static String readOutputLocation(File paramsFile) {
		try (BufferedReader br = new BufferedReader(new java.io.FileReader(paramsFile, StandardCharsets.UTF_8))) {
			String line;
			while ((line = br.readLine()) != null) {
				int eq = line.indexOf('=');
				if (eq > 0 && "output_location".equals(line.substring(0, eq).trim())) {
					String value = line.substring(eq + 1);
					int hash = value.indexOf('#');
					return (hash >= 0 ? value.substring(0, hash) : value).trim();
				}
			}
		} catch (IOException e) {
			LOGGER.debug("Could not read output_location from {}", paramsFile, e);
		}
		return paramsFile.getParent();
	}

	/** The highest-versioned {@code <prefix>-<version>.jar} in a folder. */
	private static File highestVersioned(File dir, String prefix) {
		File[] jars = dir.listFiles(f -> {
			String n = f.getName().toLowerCase(Locale.ROOT);
			return f.isFile() && n.startsWith(prefix.toLowerCase(Locale.ROOT)) && n.endsWith(".jar");
		});
		if (jars == null || jars.length == 0) {
			return null;
		}
		File best = null;
		String bestKey = null;
		for (File jar : jars) {
			java.util.regex.Matcher m = java.util.regex.Pattern.compile("(\\d+(?:\\.\\d+)*)").matcher(jar.getName());
			StringBuilder key = new StringBuilder();
			if (m.find()) {
				for (String part : m.group(1).split("\\.")) {
					key.append(String.format(Locale.ROOT, "%06d", Integer.parseInt(part)));
				}
			}
			String k = key.toString();
			if (bestKey == null || k.compareTo(bestKey) > 0) {
				bestKey = k;
				best = jar;
			}
		}
		return best;
	}

	/** {@code 20.0} reads as a tolerance someone chose; {@code 20} reads as the default it is. */
	private static String trim(double value) {
		if (value == Math.rint(value) && !Double.isInfinite(value)) {
			return Long.toString((long) value);
		}
		return String.format(Locale.ROOT, "%s", value);
	}
}

package uk.ac.quadram.metapilot.core.msfragger;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.msFragger.MSFraggerParameter;

/**
 * OS-agnostic runner for MSFragger in standalone mode, for the open-search
 * pipeline (see {@code docs/developer/open-search-standalone-design.md}).
 * Deliberately mirrors {@link uk.ac.quadram.metapilot.core.fragpipe.FragPipeHeadlessRunner}:
 * a single process invocation returning an exit code, living in
 * {@code metapilot-core} so the GUI, CLI ({@code MetaPilotCliApp}) and HPC
 * ({@code HpcStudyRunner}) entry points can all share it. An earlier draft put
 * this in {@code metapilot-windows-gui} modelled on {@code SageTask}, which is
 * the older architectural generation and unreachable from CLI/HPC.
 * <p>
 * Command shape:
	 * {@code <current JVM>/bin/java -Xmx<heap> -jar MSFragger.jar <params file> <mzML files...>}
 * <p>
 * No FragPipe, no Philosopher. Output is pepXML written next to each input
 * spectra file; MSFragger has no output-directory option.
 */
public final class MSFraggerHeadlessRunner {

	private static final String OPEN_TEMPLATE_RESOURCE = "/msfragger/open_fragger_standalone.params";

	/** Matches the tag written by {@code FastaManager.addDecoy}. */
	public static final String DECOY_TAG = "rev_";

	private MSFraggerHeadlessRunner() {
	}

	/**
	 * Maps a user-facing enzyme name to MSFragger's
	 * {@code {name, cutAfter, notCutBefore}}. Same mapping the FragPipe MAG
	 * panel uses ({@code getEnzymeInfo}), so an enzyme selected in either
	 * workflow digests identically. Without this the enzyme name would be
	 * written while the cleavage rules silently stayed tryptic.
	 */
	/**
	 * Enzyme choices offered by the MSFragger-family DDA MAG panels, in display
	 * order.
	 * <p>
	 * The single source of truth for both the UI list and {@link #enzymeInfo}, so
	 * the two cannot drift: every entry here has a case in that switch, and
	 * nothing outside this list can reach it except through the documented
	 * fallback.
	 * <p>
	 * Deliberately NOT read from {@code enzyme_default.ini}. That file is pFind's
	 * configuration vocabulary and defines {@code EN1=Trypsin KR _ C} - trypsin
	 * with <em>no</em> proline restriction - and {@code EN2=Trypsin_P KR P C} with
	 * one, which is the opposite of what those two names mean everywhere else,
	 * including MetaPilot's own {@code core.enzyme.Enzyme} enum and its DIA-NN
	 * mapping. Reading the list from there is what propagated that inversion into
	 * the MSFragger workflows; see {@link #enzymeInfo}.
	 */
	public static final java.util.List<String> ENZYME_CHOICES = java.util.List.of(
			"Trypsin/P", "Trypsin", "LysC/P", "LysC", "ArgC", "AspN", "GluC", "GluC-bicarb",
			"Chymotrypsin", "Elastase", "Proteinase K", "Thermolysin", "Nonspecific");

	/**
	 * The workflow default, first in {@link #ENZYME_CHOICES} because the panels
	 * preselect element 0.
	 * <p>
	 * Trypsin/P - MSFragger's {@code stricttrypsin}, cleaving after K/R including
	 * before proline - for three reasons. It is close to a superset: where the
	 * proline rule leaves a K-P-spanning peptide as the only candidate, strict
	 * digestion offers that peptide and both split forms. It matches the
	 * prevailing default elsewhere (MaxQuant ships Trypsin/P; FragPipe's own
	 * workflows use {@code stricttrypsin}), so MetaPilot and a hand-run FragPipe
	 * agree. And it is the digest every MetaPilot DDA MAG result to date was
	 * actually produced with, back when the mapping was inverted, so keeping it
	 * means the naming correction does not silently invalidate existing studies.
	 * <p>
	 * The trade is not entirely one-way: under strict digestion a K-P site
	 * consumes one of the allowed missed cleavages, where under the proline rule
	 * it consumes none, so a peptide spanning two K-P sites plus a genuine missed
	 * cleavage fits the budget under "Trypsin" and not under "Trypsin/P". Users
	 * who want the conventional Keil-rule digest should select "Trypsin", which
	 * now genuinely performs one.
	 */
	public static final String DEFAULT_ENZYME = "Trypsin/P";

	/**
	 * Maps a user-facing enzyme name to MSFragger's
	 * {@code {name, cutAfter, notCutBefore}}.
	 * <p>
	 * One definition for every MSFragger-family workflow: the FragPipe MAG runner
	 * and the open-search runner both call this, and the open-search runner also
	 * feeds Sage's {@code cleave_at}/{@code restrict} from it, so an enzyme chosen
	 * in any of them digests identically. It previously existed as two copies with
	 * different normalisation rules that happened to agree.
	 *
	 * <h2>The proline rule, and a correction</h2>
	 * MSFragger's own vocabulary is that {@code trypsin} cleaves after K/R but not
	 * before P, and {@code stricttrypsin} cleaves after K/R regardless. Until this
	 * was corrected, both copies of this mapping sent {@code "Trypsin"} to
	 * {@code stricttrypsin} with an empty restriction and {@code "Trypsin/P"} to
	 * {@code trypsin} with {@code P} - exactly backwards, and backwards from
	 * {@code core.enzyme.Enzyme} ({@code Trypsin("Trypsin", true, "KR", "P")}) and
	 * from the DIA-NN mapping ({@code Trypsin -> K*,R*,!*P}) in this same codebase.
	 * The inversion came from {@code enzyme_default.ini}, which the panels read to
	 * populate the dropdown; see {@link #ENZYME_CHOICES}.
	 * <p>
	 * Selecting "Trypsin" therefore used to perform a Trypsin/P digest. It now
	 * performs a trypsin digest, so results from before and after this change are
	 * not directly comparable - re-run rather than mixing them. "Trypsin/P"
	 * remains available and is what the old "Trypsin" actually did.
	 * <p>
	 * The LysC pair carried the same inversion and is corrected with it. The
	 * remaining enzymes have no proline ambiguity in MSFragger's defaults and are
	 * unchanged.
	 */
	public static String[] enzymeInfo(String enzyme) {
		if (enzyme == null) {
			return new String[] { "trypsin", "KR", "P" };
		}
		switch (enzyme.toLowerCase(java.util.Locale.ROOT).replace("-", "").replace("_", "").replace(" ", "")
				.replace("/", "")) {
		case "trypsin":
			return new String[] { "trypsin", "KR", "P" };
		case "trypsinp":
		case "stricttrypsin":
			return new String[] { "stricttrypsin", "KR", "" };
		case "chymotrypsin":
			return new String[] { "chymotrypsin", "FWY", "" };
		case "gluc":
			return new String[] { "gluc", "DE", "" };
		case "glucbicarb":
			return new String[] { "gluc-bicarb", "DE", "" };
		case "lysc":
			return new String[] { "lysc", "K", "P" };
		case "lyscp":
			return new String[] { "lysc-p", "K", "" };
		case "argc":
			return new String[] { "argc", "R", "" };
		case "aspn":
			return new String[] { "aspn", "D", "" };
		case "elastase":
			return new String[] { "elastase", "AGILV", "" };
		case "proteinasek":
			return new String[] { "proteinasek", "AEFILMVWY", "" };
		case "thermolysin":
			return new String[] { "thermolysin", "AEFILMV", "" };
		case "nonspecific":
		case "noenzyme":
			return new String[] { "nonspecific", "", "" };
		default:
			// Trypsin, not stricttrypsin: an unrecognised name should fall back to
			// the conventional digest rather than the permissive one.
			return new String[] { "trypsin", "KR", "P" };
		}
	}

	/** Settings that vary per run; everything else comes from the bundled template. */
	public static final class OpenSearchSettings {
		public File targetDecoyFasta;
		public int threads = 0; // 0 = MSFragger auto-detect
		public String enzymeName = "trypsin";
		public String cutAfter = "KR";
		public String notCutBefore = "P";
		public int enzymeTermini = 2;
		public int missedCleavages = 2;
		public double precursorMassLower = -150.0;
		public double precursorMassUpper = 500.0;
		/** Units of the precursor window: 0 = Da (open search), 1 = ppm (closed search). */
		public int precursorMassUnits = 0;
		/** MSFragger {@code isotope_error}; "0" for open searches, see {@code UnmodifiedPsmSelector}. */
		public String isotopeError = "0";
		/** {@code variable_mod_01}, {@code _02}, ... in order, e.g. "15.9949 M 3"; none by default. */
		public java.util.List<String> variableMods = new java.util.ArrayList<>();
		/**
		 * PSM candidates reported per spectrum. MSFragger defaults to 1
		 * (top-hit only), which would silently defeat this pipeline's
		 * requirement that the classifier sees every candidate.
		 */
		public int outputReportTopN = 5;
		/**
		 * MSFragger's {@code output_format}. {@code pepXML_pin} additionally writes
		 * the Percolator feature table MSBooster consumes, while leaving the pepXML
		 * this pipeline's reader depends on exactly as it was.
		 */
		public String outputFormat = "pepXML";
		public String decoyTag = DECOY_TAG;
	}

	/**
	 * Writes a fragger.params for an open search into {@code workDir}, starting
	 * from the bundled template and overriding the per-run values.
	 */
	public static File writeOpenSearchParams(File workDir, OpenSearchSettings settings) throws IOException {
		if (workDir == null) {
			throw new IOException("workDir is null");
		}
		if (settings == null || settings.targetDecoyFasta == null || !settings.targetDecoyFasta.exists()) {
			throw new IOException("Target-decoy FASTA not found: "
					+ (settings == null || settings.targetDecoyFasta == null ? "null"
							: settings.targetDecoyFasta.getAbsolutePath()));
		}
		workDir.mkdirs();

		// The template is scratch: MSFraggerParameter can only parse from a file,
		// so the bundled resource is staged to a temp file and deleted straight
		// after. Leaving it in workDir would put a second, stale-looking .params
		// beside the real one in a folder whose purpose is to record exactly what
		// was run.
		File template = File.createTempFile("open_fragger_template", ".params");
		MSFraggerParameter par;
		try {
			try (InputStream in = MSFraggerHeadlessRunner.class.getResourceAsStream(OPEN_TEMPLATE_RESOURCE)) {
				if (in == null) {
					throw new IOException(
							"Bundled MSFragger template missing from classpath: " + OPEN_TEMPLATE_RESOURCE);
				}
				Files.copy(in, template.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
			par = MSFraggerParameter.parse(template);
		} finally {
			template.delete();
		}
		if (par == null) {
			throw new IOException("Could not parse bundled MSFragger template: " + OPEN_TEMPLATE_RESOURCE);
		}

		par.setDatabase(settings.targetDecoyFasta.getAbsolutePath());
		par.setThread(settings.threads);

		// Set the open window directly rather than via setPrecursorMassTol(...),
		// which derives lower/upper from a single symmetric tolerance and so
		// cannot express an asymmetric open-search window such as -150/+500.
		par.set("precursor_mass_lower", String.valueOf(settings.precursorMassLower));
		par.set("precursor_mass_upper", String.valueOf(settings.precursorMassUpper));
		par.set("precursor_mass_units", String.valueOf(settings.precursorMassUnits));
		par.set("isotope_error", settings.isotopeError);
		for (int i = 0; i < settings.variableMods.size(); i++) {
			par.setVariableMod(i + 1, settings.variableMods.get(i));
		}
		par.set("decoy_prefix", settings.decoyTag);
		par.set("output_report_topN", String.valueOf(settings.outputReportTopN));
		if (settings.outputFormat != null && !settings.outputFormat.isBlank()) {
			par.set("output_format", settings.outputFormat.trim());
		}

		// Current-generation MSFragger enzyme keys (search_enzyme_*_1). Not using
		// MSFraggerParameter.setEnzyme(...), which targets the pre-dual-enzyme
		// schema (search_enzyme_name, no _1 suffix) and would silently no-op.
		par.set("search_enzyme_name_1", settings.enzymeName);
		par.set("search_enzyme_cut_1", settings.cutAfter);
		par.set("search_enzyme_nocut_1", settings.notCutBefore == null ? "" : settings.notCutBefore);
		par.set("num_enzyme_termini", String.valueOf(settings.enzymeTermini));
		par.set("allowed_missed_cleavage_1", String.valueOf(settings.missedCleavages));

		File paramsFile = new File(workDir, OpenSearchLayout.SEARCH_PARAMS_FILE);
		par.export(paramsFile);
		return paramsFile;
	}

	/**
	 * Renders exactly what {@link #writeOpenSearchParams} would write, without
	 * touching {@code workDir} or leaving anything behind. Used to detect a
	 * configuration change before deciding whether a previous stage-2 search
	 * can be reused instead of re-run - see the staleness check in
	 * {@code MagOpenSearchHeadlessRunner.runSampleDbSearches}.
	 * <p>
	 * Found 2026-09-17: without this, a changed {@code output_report_topN} (or
	 * any other setting derived from the pipeline parameters) left a stale
	 * fragger_open.params - and the pepXML searched under it - on disk
	 * indefinitely, because the only staleness check that existed compared
	 * timestamps against the refined FASTA, which does not change when only
	 * the rank policy does.
	 */
	public static String previewOpenSearchParamsText(OpenSearchSettings settings) throws IOException {
		Path tempDir = Files.createTempDirectory("open_fragger_preview");
		try {
			File written = writeOpenSearchParams(tempDir.toFile(), settings);
			return new String(Files.readAllBytes(written.toPath()), StandardCharsets.UTF_8);
		} finally {
			// Best-effort cleanup of the scratch directory; a leaked preview
			// params file is harmless, just untidy.
			File[] children = tempDir.toFile().listFiles();
			if (children != null) {
				for (File child : children) {
					child.delete();
				}
			}
			tempDir.toFile().delete();
		}
	}

	/**
	 * Runs MSFragger over one or more spectra files.
	 * <p>
	 * <b>Batching matters.</b> MSFragger reuses its fragment index across files
	 * passed to a single invocation, so searching N files in one call is
	 * substantially faster than N calls. That is available at <b>stage 1</b>,
	 * where every sample searches the same marker/HAP database — pass all
	 * samples at once. It is <b>not</b> available at stage 2, where each sample
	 * has its own refined database, so stage 2 must invoke once per sample.
	 *
	 * @param msfraggerJar MSFragger.jar
	 * @param javaMaxHeap  {@code -Xmx} value, e.g. "32g"
	 * @param paramsFile   from {@link #writeOpenSearchParams}
	 * @param spectraFiles mzML/mzXML inputs
	 * @return process exit code (0 = success)
	 */
	public static int run(File msfraggerJar, String javaMaxHeap, File paramsFile, List<File> spectraFiles,
			ProgressCallback cb) throws IOException, InterruptedException {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (msfraggerJar == null || !msfraggerJar.exists()) {
			throw new IOException("MSFragger jar not found: "
					+ (msfraggerJar == null ? "null" : msfraggerJar.getAbsolutePath()));
		}
		if (paramsFile == null || !paramsFile.exists()) {
			throw new IOException("MSFragger params file not found: "
					+ (paramsFile == null ? "null" : paramsFile.getAbsolutePath()));
		}
		if (spectraFiles == null || spectraFiles.isEmpty()) {
			throw new IOException("No spectra files supplied to MSFragger");
		}

		List<String> cmd = new ArrayList<>();
		cmd.add(resolveJavaCommand());
		if (javaMaxHeap != null && !javaMaxHeap.isBlank()) {
			cmd.add("-Xmx" + javaMaxHeap);
		}
		cmd.add("-jar");
		cmd.add(msfraggerJar.getAbsolutePath());
		cmd.add(paramsFile.getAbsolutePath());
		for (File spectra : spectraFiles) {
			if (spectra == null || !spectra.exists()) {
				throw new IOException("Spectra file not found: " + spectra);
			}
			cmd.add(spectra.getAbsolutePath());
		}

		cb.onMessage("Running MSFragger: " + String.join(" ", cmd));

		ProcessBuilder pb = new ProcessBuilder(cmd);
		// MSFragger resolves its ext/ folder (Thermo/Bruker readers) relative to
		// the jar, so run from the jar's directory.
		File parent = msfraggerJar.getParentFile();
		if (parent != null && parent.exists()) {
			pb.directory(parent);
		}
		pb.redirectErrorStream(true);

		Process p = pb.start();
		try (BufferedReader br = new BufferedReader(
				new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (Thread.currentThread().isInterrupted()) {
					p.destroyForcibly();
					throw new InterruptedException("Interrupted while running MSFragger");
				}
				String trimmed = line.trim();
				if (!trimmed.isEmpty()) {
					cb.onMessage(trimmed);
				}
			}
		}

		int exit = p.waitFor();
		if (exit != 0) {
			cb.onMessage("MSFragger exit code: " + exit);
		}
		return exit;
	}

	/**
	 * Uses the Java runtime that is already running MetaPilot. In a packaged
	 * Windows application this is {@code <MetaPilot>/runtime/bin/java.exe}; in
	 * Eclipse it is the configured launch JDK/JRE. Falling back to the PATH keeps
	 * the runner usable in unusual JVM layouts that do not contain a launcher.
	 */
	public static String resolveJavaCommand() {
		return resolveJavaCommand(System.getProperty("java.home"), System.getProperty("os.name"));
	}

	static String resolveJavaCommand(String javaHome, String osName) {
		boolean windows = osName != null && osName.toLowerCase(java.util.Locale.ROOT).contains("win");
		String executableName = windows ? "java.exe" : "java";
		if (javaHome != null && !javaHome.isBlank()) {
			Path executable = Path.of(javaHome, "bin", executableName);
			if (Files.isRegularFile(executable)) {
				return executable.toAbsolutePath().normalize().toString();
			}
		}
		return "java";
	}

	/**
	 * The pepXML MSFragger produces for a given spectra file. MSFragger writes
	 * output next to the input, named after its base name; there is no
	 * output-directory parameter to redirect this.
	 */
	/**
	 * The {@code .pin} MSFragger writes beside a spectra file when
	 * {@code output_format} includes {@code pin}. Same convention as
	 * {@link #expectedPepXml(File)}.
	 */
	public static File expectedPin(File spectraFile) {
		String name = spectraFile.getName();
		int dot = name.lastIndexOf('.');
		String base = dot > 0 ? name.substring(0, dot) : name;
		return new File(spectraFile.getParentFile(), base + ".pin");
	}

	public static File expectedPepXml(File spectraFile) {
		String name = spectraFile.getName();
		int dot = name.lastIndexOf('.');
		String base = dot > 0 ? name.substring(0, dot) : name;
		return new File(spectraFile.getParentFile(), base + ".pepXML");
	}
}

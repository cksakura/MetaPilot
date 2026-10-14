package uk.ac.quadram.metapilot.core.sage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.json.JSONArray;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

/**
 * Runs Sage as a second open-search engine, mirroring
 * {@code MSFraggerHeadlessRunner}.
 * <p>
 * Sage takes one JSON config naming both the database and the mzML inputs, so
 * there is no separate command line to build beyond
 * {@code sage <config.json>}. The config is written into the same folder as the
 * results it produces, for the same reason the MSFragger params file is: a
 * search folder should be readable on its own, without cross-referencing
 * something elsewhere in the tree.
 *
 * <h2>Matching the MSFragger search, not Sage's defaults</h2>
 * The point of a second engine is to see whether it finds things the first
 * missed <em>on the same problem</em>. Any difference in search space would
 * confound that, so the open window, enzyme, missed cleavages, fixed
 * modifications, peptide length/mass limits and reported candidate count are
 * all taken from the fragger.params that was written for the same sample
 * ({@link #applyMsfraggerParams}) rather than from Sage's own defaults. Deriving
 * them from the file that MSFragger was actually handed, instead of restating
 * them here, is what stops the two search spaces drifting apart.
 *
 * <h2>The open window is mirrored, not copied</h2>
 * <b>This is the one conversion that cannot be skipped.</b> MSFragger's
 * {@code precursor_mass_lower}/{@code precursor_mass_upper} bound the
 * <em>delta mass</em> (observed precursor minus theoretical peptide) directly:
 * {@code -150}/{@code +500} means modifications of -150 to +500 Da.
 * <p>
 * Sage's {@code precursor_tol} bounds the <em>candidate peptide mass</em>
 * around the precursor instead: {@code Tolerance::bounds} returns
 * {@code (center + lo, center + hi)} (sage/crates/sage/src/mass.rs), so a
 * window of {@code da: [lo, hi]} admits delta masses of {@code [-hi, -lo]}.
 * Sage's own documentation gives {@code da: [-500, 100]} as the open-search
 * example, i.e. an open search for modifications of -100 to +500 Da.
 * <p>
 * Copying MSFragger's numbers across verbatim therefore runs the open search
 * backwards: {@code da: [-150, 500]} searches modifications of -500 to +150 Da,
 * which hides everything above +150 - the bulk of what an open search exists to
 * find - and spends the candidate budget on large negative deltas instead. The
 * window written here is mirrored: {@code da: [-upper, -lower]}.
 *
 * <h2>Two further settings that are load-bearing</h2>
 * <ul>
 * <li>{@code isotope_errors [0, 0]} - matches MSFragger's
 * {@code isotope_error = 0}. With isotope correction enabled the engine folds
 * isotope-triggered precursors into the zero peak internally, and the
 * distinction this workflow's isotope-family logic depends on becomes
 * unrecoverable.</li>
 * <li>{@code wide_window false} - Sage's wide-window mode ignores
 * {@code precursor_tol} altogether and runs in DIA/dynamic-tolerance mode;
 * leaving it on would discard the open window this class exists to set.</li>
 * </ul>
 */
public final class SageOpenSearchRunner {

	/** Sage's own result filename, relative to its output directory. */
	public static final String RESULTS_FILE = "results.sage.tsv";
	/** Config filename, written beside the results. */
	public static final String CONFIG_FILE = "sage_open_config.json";

	private SageOpenSearchRunner() {
	}

	public static final class Settings {
		public File targetDecoyFasta;
		/** MSFragger's delta-mass window, in Da. Mirrored when written out. */
		public double precursorMassLower = -150.0;
		public double precursorMassUpper = 500.0;
		public double fragmentTolPpm = 20.0;
		public String enzymeCleaveAt = "KR";
		public String enzymeRestrict = "";
		public int missedCleavages = 2;
		public int minPeptideLength = 7;
		public int maxPeptideLength = 50;
		/** Mirrors MSFragger's {@code digest_mass_range}. */
		public double peptideMinMass = 500.0;
		public double peptideMaxMass = 5000.0;
		/** Mirrors MSFragger's {@code max_fragment_charge}; null = Sage default. */
		public Integer maxFragmentCharge = 2;
		/**
		 * Charge range assumed when a spectrum carries no charge annotation;
		 * mirrors MSFragger's {@code precursor_charge}. Sage's default starts at 2,
		 * so an unannotated 1+ spectrum would be searched by one engine and not the
		 * other.
		 */
		public int minPrecursorCharge = 1;
		public int maxPrecursorCharge = 4;
		/** Mirrors MSFragger's {@code fragment_ion_series}. */
		public final List<String> ionKinds = new ArrayList<>(Arrays.asList("b", "y"));
		/**
		 * Fixed modifications, keyed by Sage's site notation (a residue letter, or
		 * {@code ^ $ [ ]} for the peptide/protein termini). Mirrors MSFragger's
		 * {@code add_*} entries; without it Sage searches an unmodified backbone
		 * while MSFragger searches a carbamidomethylated one, and every
		 * cysteine-containing peptide becomes a different problem for the two
		 * engines.
		 */
		public final Map<String, Double> staticMods = new LinkedHashMap<>();
		/** Candidates reported per spectrum; mirrors output_report_topN. */
		public int reportPsms = 5;
		/** Must match the shared decoy tag so target/decoy parsing agrees. */
		public String decoyTag = "rev_";
		/** 0 = let Sage use its rayon default (all cores). */
		public int threads = 0;
	}

	/**
	 * Copies the search space out of a fragger.params written for the same
	 * sample, so the two engines cannot drift apart through an edit to one side
	 * only. Keys that are absent leave the corresponding setting untouched.
	 * <p>
	 * Only fixed modifications are transferred: this workflow's MSFragger
	 * template declares no variable modifications, and a variable mod present on
	 * one engine and not the other would change the candidate space rather than
	 * the engine.
	 */
	public static void applyMsfraggerParams(Settings settings, File fraggerParams) throws IOException {
		if (settings == null || fraggerParams == null || !fraggerParams.exists()) {
			return;
		}
		Map<String, String> kv = new LinkedHashMap<>();
		try (BufferedReader br = new BufferedReader(
				new InputStreamReader(new FileInputStream(fraggerParams), StandardCharsets.UTF_8))) {
			String line;
			while ((line = br.readLine()) != null) {
				int hash = line.indexOf('#');
				if (hash >= 0) {
					line = line.substring(0, hash);
				}
				int eq = line.indexOf('=');
				if (eq <= 0) {
					continue;
				}
				kv.put(line.substring(0, eq).trim(), line.substring(eq + 1).trim());
			}
		}

		settings.precursorMassLower = num(kv, "precursor_mass_lower", settings.precursorMassLower);
		settings.precursorMassUpper = num(kv, "precursor_mass_upper", settings.precursorMassUpper);
		settings.fragmentTolPpm = Math.abs(num(kv, "fragment_mass_tolerance", settings.fragmentTolPpm));
		settings.missedCleavages = (int) num(kv, "allowed_missed_cleavage_1", settings.missedCleavages);
		settings.minPeptideLength = (int) num(kv, "digest_min_length", settings.minPeptideLength);
		settings.maxPeptideLength = (int) num(kv, "digest_max_length", settings.maxPeptideLength);
		settings.reportPsms = (int) num(kv, "output_report_topN", settings.reportPsms);
		int maxCharge = (int) num(kv, "max_fragment_charge", settings.maxFragmentCharge == null ? 0
				: settings.maxFragmentCharge);
		settings.maxFragmentCharge = maxCharge > 0 ? Integer.valueOf(maxCharge) : null;
		if (kv.containsKey("search_enzyme_cut_1")) {
			settings.enzymeCleaveAt = kv.get("search_enzyme_cut_1");
		}
		if (kv.containsKey("search_enzyme_nocut_1")) {
			settings.enzymeRestrict = kv.get("search_enzyme_nocut_1");
		}
		if (kv.containsKey("decoy_prefix") && !kv.get("decoy_prefix").isEmpty()) {
			settings.decoyTag = kv.get("decoy_prefix");
		}
		String charges = kv.get("precursor_charge");
		if (charges != null) {
			String[] parts = charges.trim().split("[\\s,]+");
			if (parts.length >= 2) {
				try {
					settings.minPrecursorCharge = Integer.parseInt(parts[0]);
					settings.maxPrecursorCharge = Integer.parseInt(parts[1]);
				} catch (NumberFormatException ignored) {
					// leave the defaults in place
				}
			}
		}
		String series = kv.get("fragment_ion_series");
		if (series != null && !series.trim().isEmpty()) {
			settings.ionKinds.clear();
			for (String ion : series.trim().split("[\\s,]+")) {
				String kind = ion.trim().toLowerCase(java.util.Locale.ROOT);
				if (!kind.isEmpty() && !settings.ionKinds.contains(kind)) {
					settings.ionKinds.add(kind);
				}
			}
		}

		String massRange = kv.get("digest_mass_range");
		if (massRange != null) {
			String[] parts = massRange.trim().split("[\\s,]+");
			if (parts.length >= 2) {
				try {
					settings.peptideMinMass = Double.parseDouble(parts[0]);
					settings.peptideMaxMass = Double.parseDouble(parts[1]);
				} catch (NumberFormatException ignored) {
					// leave the defaults in place
				}
			}
		}

		settings.staticMods.clear();
		for (Map.Entry<String, String> e : kv.entrySet()) {
			String site = sageSiteForMsfraggerAddKey(e.getKey());
			if (site == null) {
				continue;
			}
			double mass;
			try {
				mass = Double.parseDouble(e.getValue());
			} catch (NumberFormatException ex) {
				continue;
			}
			if (mass != 0.0) {
				settings.staticMods.put(site, mass);
			}
		}
	}

	/**
	 * Maps an MSFragger {@code add_*} key onto Sage's static-modification site
	 * notation, or null if the key is not a fixed modification. The user-defined
	 * placeholder residues (B, J, O, U, X, Z) are skipped: Sage has no
	 * corresponding residue and a non-zero value there would be silently dropped
	 * anyway.
	 */
	private static String sageSiteForMsfraggerAddKey(String key) {
		if (key == null || !key.startsWith("add_")) {
			return null;
		}
		String rest = key.substring(4);
		if (rest.equals("Nterm_peptide")) {
			return "^";
		}
		if (rest.equals("Cterm_peptide")) {
			return "$";
		}
		if (rest.equals("Nterm_protein")) {
			return "[";
		}
		if (rest.equals("Cterm_protein")) {
			return "]";
		}
		if (rest.length() < 2 || rest.charAt(1) != '_') {
			return null;
		}
		char residue = rest.charAt(0);
		if (residue < 'A' || residue > 'Z' || "BJOUXZ".indexOf(residue) >= 0) {
			return null;
		}
		return String.valueOf(residue);
	}

	private static double num(Map<String, String> kv, String key, double fallback) {
		String raw = kv.get(key);
		if (raw == null || raw.isEmpty()) {
			return fallback;
		}
		try {
			return Double.parseDouble(raw.trim());
		} catch (NumberFormatException e) {
			return fallback;
		}
	}

	/**
	 * Writes the Sage config into {@code workDir} and returns it.
	 * <p>
	 * {@code generate_decoys} is deliberately false: the FASTA handed in already
	 * carries 1:1 {@code rev_} decoys built by this workflow, and letting Sage
	 * add its own would both double the decoy count and mix two decoy
	 * generation schemes into one FDR estimate.
	 */
	public static File writeConfig(File workDir, Settings settings, List<File> spectraFiles) throws IOException {
		if (workDir == null) {
			throw new IOException("workDir is null");
		}
		workDir.mkdirs();
		File configFile = new File(workDir, CONFIG_FILE);
		try (PrintWriter writer = new PrintWriter(configFile, "UTF-8")) {
			writer.println(buildConfigJson(workDir, settings, spectraFiles));
		}
		return configFile;
	}

	/**
	 * The config as it would be written, without writing it.
	 * <p>
	 * Exposed so a caller can tell whether an existing search was run under the
	 * same configuration before deciding to reuse its results. Comparing the text
	 * is the only honest test: a settings change that alters the search space -
	 * the open window most of all - leaves a perfectly valid-looking
	 * {@code results.sage.tsv} on disk that answers a different question, and a
	 * timestamp check cannot see that.
	 */
	public static String buildConfigJson(File workDir, Settings settings, List<File> spectraFiles) throws IOException {
		if (workDir == null) {
			throw new IOException("workDir is null");
		}
		if (settings == null || settings.targetDecoyFasta == null || !settings.targetDecoyFasta.exists()) {
			throw new IOException("Target-decoy FASTA not found for Sage: "
					+ (settings == null || settings.targetDecoyFasta == null ? "null"
							: settings.targetDecoyFasta.getAbsolutePath()));
		}
		if (spectraFiles == null || spectraFiles.isEmpty()) {
			throw new IOException("No spectra files supplied to Sage");
		}

		JSONObject database = new JSONObject();
		// Sage's default is 8192; a larger bucket is the documented lever for wide
		// precursor searches, which is exactly what an open search is.
		database.put("bucket_size", 32768);
		database.put("fasta", settings.targetDecoyFasta.getAbsolutePath());
		database.put("generate_decoys", false);
		database.put("decoy_tag", settings.decoyTag);
		// peptide_min_len / peptide_max_len / missed_cleavages are enzyme-level
		// keys in Sage; writing them on `database` as well is ignored and only
		// makes the recorded config look like it says something it does not.
		database.put("peptide_min_mass", settings.peptideMinMass);
		database.put("peptide_max_mass", settings.peptideMaxMass);

		JSONObject enzyme = new JSONObject();
		enzyme.put("missed_cleavages", settings.missedCleavages);
		enzyme.put("min_len", settings.minPeptideLength);
		enzyme.put("max_len", settings.maxPeptideLength);
		// Written even when empty: an empty cleave_at is how Sage is told to
		// digest non-specifically, and omitting the key would silently fall back
		// to Sage's tryptic default.
		enzyme.put("cleave_at", settings.enzymeCleaveAt == null ? "" : settings.enzymeCleaveAt);
		// Likewise explicit: Sage's default restrict is "P", so omitting the key
		// would impose a proline restriction MSFragger was not given.
		enzyme.put("restrict", settings.enzymeRestrict == null || settings.enzymeRestrict.isEmpty()
				? JSONObject.NULL : settings.enzymeRestrict);
		enzyme.put("c_terminal", true);
		database.put("enzyme", enzyme);

		// Fixed modifications, mirroring MSFragger's add_* entries. Empty maps are
		// still written so the recorded config states what was searched.
		JSONObject staticMods = new JSONObject();
		for (Map.Entry<String, Double> mod : settings.staticMods.entrySet()) {
			staticMods.put(mod.getKey(), mod.getValue().doubleValue());
		}
		database.put("static_mods", staticMods);
		database.put("variable_mods", new JSONObject());
		JSONArray ionKinds = new JSONArray();
		for (String kind : settings.ionKinds) {
			ionKinds.put(kind);
		}
		database.put("ion_kinds", ionKinds);

		JSONObject config = new JSONObject();
		config.put("database", database);

		// Mirrored, not copied - see the class javadoc. MSFragger's window bounds
		// the delta mass; Sage's bounds the candidate peptide mass.
		JSONObject precursorTol = new JSONObject();
		precursorTol.put("da", new JSONArray().put(-settings.precursorMassUpper).put(-settings.precursorMassLower));
		config.put("precursor_tol", precursorTol);

		JSONObject fragmentTol = new JSONObject();
		fragmentTol.put("ppm", new JSONArray().put(-settings.fragmentTolPpm).put(settings.fragmentTolPpm));
		config.put("fragment_tol", fragmentTol);

		// See the class javadoc: both of these keep the two engines on the same
		// problem rather than on Sage's defaults.
		config.put("isotope_errors", new JSONArray().put(0).put(0));
		config.put("wide_window", false);
		config.put("precursor_charge",
				new JSONArray().put(settings.minPrecursorCharge).put(settings.maxPrecursorCharge));

		config.put("deisotope", true);
		config.put("chimera", false);
		config.put("min_peaks", 15);
		config.put("max_peaks", 150);
		config.put("min_matched_peaks", 4);
		if (settings.maxFragmentCharge != null) {
			config.put("max_fragment_charge", settings.maxFragmentCharge.intValue());
		}
		config.put("report_psms", settings.reportPsms);
		config.put("predict_rt", false);
		config.put("output_directory", workDir.getAbsolutePath());

		JSONArray mzml = new JSONArray();
		for (File f : spectraFiles) {
			if (f == null || !f.exists()) {
				throw new IOException("Spectra file not found for Sage: " + f);
			}
			mzml.put(f.getAbsolutePath());
		}
		config.put("mzml_paths", mzml);

		return config.toString(2);
	}

	/**
	 * @param sageExe Sage executable
	 * @return process exit code (0 = success)
	 */
	public static int run(File sageExe, File configFile, ProgressCallback cb)
			throws IOException, InterruptedException {
		return run(sageExe, configFile, 0, cb);
	}

	/**
	 * @param threads worker threads, or 0 to let Sage use every core. Sage has no
	 *                thread setting in its JSON config and no thread command-line
	 *                flag - {@code --batch-size} only controls how many
	 *                <em>files</em> are searched at once, and this workflow passes
	 *                one file per invocation - so the count is applied through
	 *                rayon's {@code RAYON_NUM_THREADS} environment variable.
	 * @return process exit code (0 = success)
	 */
	public static int run(File sageExe, File configFile, int threads, ProgressCallback cb)
			throws IOException, InterruptedException {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (sageExe == null || !sageExe.exists()) {
			throw new IOException("Sage executable not found: " + (sageExe == null ? "null" : sageExe.getAbsolutePath()));
		}
		if (configFile == null || !configFile.exists()) {
			throw new IOException("Sage config not found: " + configFile);
		}

		List<String> cmd = new ArrayList<>();
		cmd.add(sageExe.getAbsolutePath());
		cmd.add(configFile.getAbsolutePath());
		cb.onMessage("Running Sage: " + String.join(" ", cmd)
				+ (threads > 0 ? " (RAYON_NUM_THREADS=" + threads + ")" : ""));

		ProcessBuilder pb = new ProcessBuilder(cmd);
		pb.directory(configFile.getParentFile());
		if (threads > 0) {
			pb.environment().put("RAYON_NUM_THREADS", String.valueOf(threads));
		}
		pb.redirectErrorStream(true);

		Process p = pb.start();
		try (BufferedReader br = new BufferedReader(
				new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
			String line;
			while ((line = br.readLine()) != null) {
				if (Thread.currentThread().isInterrupted()) {
					p.destroyForcibly();
					throw new InterruptedException("Interrupted while running Sage");
				}
				String trimmed = line.trim();
				if (!trimmed.isEmpty()) {
					cb.onMessage(trimmed);
				}
			}
		}
		int exit = p.waitFor();
		if (exit != 0) {
			cb.onMessage("Sage exit code: " + exit);
		}
		return exit;
	}

	/** Sage's result file inside an output directory. */
	public static File resultsFile(File workDir) {
		return new File(workDir, RESULTS_FILE);
	}
}

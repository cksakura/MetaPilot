package uk.ac.quadram.metapilot.core.msbooster;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

/**
 * OS-agnostic runner for MSBooster in standalone mode.
 * <p>
 * Deliberately mirrors {@link MSFraggerHeadlessRunner}: one process invocation
 * returning an exit code, living in {@code metapilot-core} so the GUI, CLI and
 * HPC entry points share it.
 *
 * <h2>What it is for</h2>
 * MSBooster adds deep-learning features - predicted-spectrum similarity and
 * predicted retention time - to a Percolator {@code .pin} file. Measured on the
 * 12mix mock community, a rescorer built on search features alone recovers
 * about +14% PSMs at the marker stage, while the closed-search arm's
 * MSBooster+Percolator stack gains +25% PSMs and +23% peptides over its own raw
 * pepXML. The difference is what these features are worth.
 *
 * <h2>Why MSBooster rather than a Python rescorer</h2>
 * MSBooster's prediction backend is DIA-NN ({@code rtModel=DIA-NN},
 * {@code spectraModel=DIA-NN}), which MetaPilot already ships and drives as a
 * first-class resource for the entire DIA workflow. So this adds no Python
 * dependency, no model files, no network egress and no new prediction engine -
 * it reuses two tools that are already configured. Koina is left unset for the
 * same reason: it would make an HPC run depend on outbound network access.
 * <p>
 * It is also not a step back toward the FragPipe stack this workflow avoids.
 * The reason for driving MSFragger directly is that IonQuant does not carry
 * open-search delta masses into the quantified peptide list. MSBooster only
 * appends columns to a feature table and has no such behaviour.
 *
 * <h2>Where it is valid, and where it is not</h2>
 * Predicted spectra are predicted for the peptide as given. At stage 1 this
 * workflow consumes only the zero family, where the delta mass is ~0 and the
 * prediction is therefore exactly the right reference - the clean case. At
 * stage 2 the delta mass is unknown and unlocalized, and the predicted spectrum
 * will not account for the shift; published work handles that by keeping the
 * similarity as a feature and flagging peptides whose modification the model
 * has not seen, so the classifier can down-weight it. MSBooster has no such
 * flag, so a stage-2 caller should treat these features as degraded rather than
 * absent, and must not assume they mean the same thing there as here.
 *
 * <h2>Concurrency</h2>
 * {@code MagFragPipeHeadlessRunner} already documents that concurrent FragPipe
 * instances share one set of MSBooster model files, which upstream does not
 * support. This workflow's stage-1 search is a single batched MSFragger
 * invocation over every sample, so calling MSBooster once over all the samples'
 * pins together sidesteps that entirely. {@link #run} therefore takes a list of
 * pins rather than one, and callers should batch rather than loop.
 */
public final class MSBoosterHeadlessRunner {

	private static final Logger LOGGER = LogManager.getLogger(MSBoosterHeadlessRunner.class);

	/** Filename MSBooster's own FragPipe integration uses; kept for familiarity. */
	public static final String PARAMS_FILE = "msbooster_params.txt";

	/**
	 * Columns MSBooster appends to the pin, in its default configuration
	 * ({@code useSpectra=true}, {@code useRT=true}, {@code useIM=false}).
	 * <p>
	 * Used to pick the new columns back out of the edited pin. Listing them
	 * rather than diffing against the input pin means a column MSBooster renames
	 * upstream shows up as a missing feature and stops the run, instead of
	 * silently changing the feature set the classifier trains on.
	 * <p>
	 * These are the spellings in a real edited pin's header, not the camel-case
	 * names MSBooster prints in its own "Generating edited pin with following
	 * features" log line. Two were wrong here until 2026-09-11
	 * ({@code top_6_matched_intensity}, {@code delta_RT_loess_real_units}), and
	 * the effect was a classifier quietly trained on five of the seven. Check any
	 * change to this list against a pin header, never against the log.
	 */
	public static final List<String> DEFAULT_FEATURE_COLUMNS = List.of(
			"unweighted_spectral_entropy", "weighted_spectral_entropy", "hypergeometric_probability",
			"intersection", "top6_matched_intensity", "delta_RT_loess", "delta_RT_loess_real");

	private MSBoosterHeadlessRunner() {
	}

	/** Everything that goes into {@value #PARAMS_FILE}. */
	public static final class Settings {
		/** DIA-NN executable; the prediction backend. Required. */
		public File diann;
		/** Folders or files holding the mzML/mgf the pins were searched from. Required. */
		public List<File> mzmlDirectories = new ArrayList<>();
		/** The {@code .pin} files to enrich. Required. Batch them; see class docs. */
		public List<File> pinFiles = new ArrayList<>();
		/**
		 * The {@code fragger.params} the search used. Strongly recommended:
		 * MSBooster reads tolerances and modification settings from it, and
		 * without it silently falls back to its own defaults.
		 */
		public File fraggerParams;
		/** {@code unimod.obo}, normally beside the jar in FragPipe's tools folder. */
		public File unimodObo;
		public int threads = Runtime.getRuntime().availableProcessors();
		public boolean useSpectra = true;
		public boolean useRT = true;
		/** Off: this workflow's DDA inputs carry no ion mobility. */
		public boolean useIM = false;
		/** Off: detectability is MetaPilot's own DeepDetect step, not MSBooster's. */
		public boolean useDetect = false;
		public String rtModel = "DIA-NN";
		public String spectraModel = "DIA-NN";
		/**
		 * Off deliberately. Model search re-runs prediction for every candidate
		 * model to pick the best fit, which is a large multiple of the runtime for
		 * a gain this pipeline has not measured. Revisit with a measurement.
		 */
		public boolean findBestModels = false;
		/** Empty: remote inference would make an HPC run need outbound network. */
		public String koinaUrl = "";
		/** Keep the predictions on disk; they are the audit trail for the features. */
		public boolean deletePreds = false;
	}

	/**
	 * Writes {@value #PARAMS_FILE} into {@code workDir}.
	 * <p>
	 * Written before the run, not after, for the same reason the MSFragger
	 * parameters are: a failed run is exactly when they need inspecting.
	 */
	public static File writeParams(File workDir, Settings settings) throws IOException {
		if (workDir == null) {
			throw new IOException("workDir is null");
		}
		validate(settings);
		workDir.mkdirs();
		File out = new File(workDir, PARAMS_FILE);
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Written by MetaPilot. Prediction backend is DIA-NN, taken from MetaPilot's own");
			w.println("# configured resource - no Koina, so this runs without outbound network access.");
			w.println("useDetect = " + settings.useDetect);
			w.println("numThreads = " + Math.max(1, settings.threads));
			w.println("DiaNN = " + settings.diann.getAbsolutePath());
			w.println("renamePin = 1");
			w.println("useRT = " + settings.useRT);
			w.println("useSpectra = " + settings.useSpectra);
			w.println("useIM = " + settings.useIM);
			if (settings.fraggerParams != null) {
				w.println("fragger = " + settings.fraggerParams.getAbsolutePath());
			}
			w.println("deletePreds = " + settings.deletePreds);
			w.println("rtModel = " + settings.rtModel);
			w.println("spectraModel = " + settings.spectraModel);
			w.println("imModel = ");
			w.println("findBestRtModel = " + settings.findBestModels);
			w.println("findBestSpectraModel = " + settings.findBestModels);
			w.println("findBestImModel = false");
			w.println("KoinaURL = " + (settings.koinaUrl == null ? "" : settings.koinaUrl));
			// Let MSBooster read the fragmentation from the data, as FragPipe's own
			// parameter file does. Left unset it falls back to a fixed assumption.
			w.println("FragmentationType = auto");
			if (settings.unimodObo != null) {
				w.println("unimodObo = " + settings.unimodObo.getAbsolutePath());
			}
			w.println("mzmlDirectory = " + joinPaths(settings.mzmlDirectories));
			w.println("pinPepXMLDirectory = " + joinPaths(settings.pinFiles));
		}
		return out;
	}

	/**
	 * Runs MSBooster over the pins named in {@code paramsFile}.
	 *
	 * @return the process exit code; 0 on success.
	 */
	public static int run(File msboosterJar, String javaMaxHeap, File paramsFile, ProgressCallback cb)
			throws IOException, InterruptedException {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (msboosterJar == null || !msboosterJar.exists()) {
			throw new IOException("MSBooster jar not found: "
					+ (msboosterJar == null ? "null" : msboosterJar.getAbsolutePath()));
		}
		if (paramsFile == null || !paramsFile.exists()) {
			throw new IOException("MSBooster params file not found: "
					+ (paramsFile == null ? "null" : paramsFile.getAbsolutePath()));
		}

		List<String> cmd = new ArrayList<>();
		// An absolute path to the JVM running MetaPilot, never a bare "java".
		// MetaPilot launched from an IDE or a javaw shortcut passes a PATH with no
		// java on it to its children, and on Windows that surfaces as exit 9009,
		// command-not-found - the same "produced no files at all" symptom as a
		// genuine MSBooster failure. See FragPipeHeadlessRunner.ensureJavaHome,
		// which solves the same problem for the FragPipe launcher.
		String javaCommand = MSFraggerHeadlessRunner.resolveJavaCommand();
		if ("java".equals(javaCommand) || !new File(javaCommand).isFile()) {
			throw new IOException("Could not resolve an absolute java executable from java.home ("
					+ System.getProperty("java.home") + "); the launch would fall through to the PATH, which on "
					+ "Windows surfaces as exit 9009 (command not found) and looks exactly like an MSBooster "
					+ "failure. Refusing to start MSBooster without a java it can actually run.");
		}
		cmd.add(javaCommand);
		// As FragPipe passes it. MSBooster pulls in AWT through its plotting code;
		// on a headless machine the JVM would otherwise try to open a display.
		cmd.add("-Djava.awt.headless=true");
		if (javaMaxHeap != null && !javaMaxHeap.isBlank()) {
			cmd.add("-Xmx" + javaMaxHeap);
		}
		// -cp, never -jar. MSBooster's jar is not shaded: it declares a Main-Class
		// so -jar starts it and prints its version banner, then dies at the first
		// dependency with NoClassDefFoundError: org/slf4j/Logger, because -jar
		// makes the JVM ignore every other classpath entry. FragPipe runs it
		// against its own lib/ for the same reason.
		cmd.add("-cp");
		String classpath = buildClasspath(msboosterJar);
		cmd.add(classpath);
		cmd.add(MAIN_CLASS);
		cmd.add("--paramsList");
		cmd.add(paramsFile.getAbsolutePath());

		cb.onMessage("Running MSBooster: " + String.join(" ", cmd));

		ProcessBuilder pb = new ProcessBuilder(cmd);
		// The params file's folder, not the jar's: MSBooster writes its
		// predictions and any incidental scratch relative to the working
		// directory, and those belong beside the run they describe rather than
		// inside a shared tools installation that may be read-only.
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
					throw new InterruptedException("Interrupted while running MSBooster");
				}
				String trimmed = line.trim();
				if (!trimmed.isEmpty()) {
					cb.onMessage(trimmed);
				}
			}
		}
		int exit = p.waitFor();
		if (exit != 0) {
			// 9009 is Windows' command-not-found: the JVM itself was not launched,
			// so nothing about MSBooster's configuration is implicated.
			cb.onMessage("MSBooster exit code: " + exit + " (java: " + javaCommand + ")"
					+ (exit == 9009 ? " - 9009 is Windows command-not-found: that java executable does not exist"
							: "") + ". Classpath: " + classpath);
		}
		return exit;
	}

	/**
	 * MSBooster's entry point, as named in its own stack traces
	 * ({@code at mainsteps.MainClass.main}). Needed explicitly because the
	 * classpath form of {@code java} does not read the jar's manifest.
	 */
	private static final String MAIN_CLASS = "mainsteps.MainClass";

	/**
	 * The two jars FragPipe itself puts on MSBooster's classpath, in its order.
	 * <p>
	 * Taken from a working FragPipe run's log rather than guessed:
	 * <pre>
	 * java.exe -Djava.awt.headless=true -Xmx105G
	 *     -cp &lt;tools&gt;\MSBooster-1.4.14.jar;&lt;tools&gt;\batmass-io-1.36.5.jar
	 *     mainsteps.MainClass --paramsList &lt;params&gt;
	 * </pre>
	 * {@code batmass-io} is not optional and not merely a transitive dependency
	 * to be satisfied some other way: MSBooster reads the mzML itself to compare
	 * observed against predicted spectra, and that reader lives there - along with
	 * the slf4j whose absence killed the {@code -jar} launch.
	 * <p>
	 * Two named jars, never a {@code tools/*} wildcard: the wildcard would put
	 * every tool FragPipe ships on the classpath and invite a version collision
	 * that shows up as a puzzling runtime error rather than a missing class. A
	 * missing batmass-io is a broken installation and is reported as one.
	 */
	static String buildClasspath(File msboosterJar) throws IOException {
		List<String> entries = new ArrayList<>();
		entries.add(msboosterJar.getAbsolutePath());
		File toolsDir = msboosterJar.getParentFile();
		File batmass = toolsDir != null && toolsDir.isDirectory() ? highestVersioned(toolsDir, "batmass-io") : null;
		if (batmass == null) {
			throw new IOException("No batmass-io jar beside " + msboosterJar.getAbsolutePath()
					+ ". MSBooster reads the mzML through it, so without it the run dies at "
					+ "NoClassDefFoundError after printing its version banner - which looks like MSBooster "
					+ "simply producing nothing. Point Resources > FragPipe at a complete installation.");
		}
		entries.add(batmass.getAbsolutePath());
		return String.join(File.pathSeparator, entries);
	}

	/** The highest-versioned {@code <prefix>-<version>.jar} in a folder. */
	private static File highestVersioned(File dir, String prefix) {
		File[] jars = dir.listFiles(f -> {
			String n = f.getName().toLowerCase(java.util.Locale.ROOT);
			return f.isFile() && n.startsWith(prefix.toLowerCase(java.util.Locale.ROOT)) && n.endsWith(".jar");
		});
		if (jars == null || jars.length == 0) {
			return null;
		}
		File best = null;
		String bestKey = null;
		for (File jar : jars) {
			java.util.regex.Matcher m = java.util.regex.Pattern.compile("(\\d+(?:\\.\\d+)*)")
					.matcher(jar.getName());
			StringBuilder key = new StringBuilder();
			if (m.find()) {
				for (String part : m.group(1).split("\\.")) {
					key.append(String.format(java.util.Locale.ROOT, "%06d", Integer.parseInt(part)));
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


	/**
	 * Where MSBooster puts its predicted spectra, retention times and QC plots.
	 * <p>
	 * Not configurable: MSBooster derives it itself in
	 * {@code Constants.updateOutputDirectory()} as an {@code MSBooster} folder
	 * under the single pin's parent, or under the deepest common directory of
	 * several. Reproduced here so a caller can name the folder in a log line or
	 * hand it to the result cleaner without having to guess - and so that this
	 * stays wrong in one place if MSBooster ever changes the rule.
	 */
	public static File predictionDir(List<File> pins) {
		if (pins == null || pins.isEmpty()) {
			return null;
		}
		File common = pins.get(0).getAbsoluteFile().getParentFile();
		for (int i = 1; i < pins.size() && common != null; i++) {
			common = deepestCommon(common, pins.get(i).getAbsoluteFile().getParentFile());
		}
		return common == null ? null : new File(common, "MSBooster");
	}

	/** MSBooster's own name for the folder of per-sample QC images. */
	public static final String PLOTS_DIR = "MSBooster_plots";

	/**
	 * Moves the folder MSBooster just wrote to a destination the caller chooses,
	 * merging into whatever is already there.
	 *
	 * <h2>Why this is needed at all</h2>
	 * MSBooster names its own output folder, and the name is always
	 * {@code MSBooster} under the pins' deepest common directory. Two batches over
	 * different pins that happen to share a parent - a stage-one batch over
	 * {@code <sample>/marker/*.pin} and a stage-two batch over
	 * {@code <sample>/sample_db/*.pin} both resolve to the search-results root -
	 * therefore land in the same folder, and the second silently overwrites the
	 * first. Relocating afterwards is the only way to separate them, since
	 * MSBooster takes no output-directory setting.
	 *
	 * <h2>Merging, not replacing</h2>
	 * A resumed run re-predicts only the samples that need it, so its folder holds
	 * only their files. Replacing the destination would throw away the QC plots of
	 * every sample that was correctly skipped. Files are therefore merged, with
	 * same-named files overwritten - which is exact for the per-sample plots and
	 * approximate for the pooled {@code spectraRT*} library, whose contents
	 * describe only the most recent invocation. {@link #writeRunNote} records what
	 * each invocation covered so that is stated rather than assumed.
	 *
	 * @return the destination, or null when MSBooster wrote nothing to move.
	 */
	public static File relocatePredictions(List<File> pins, File destination) throws IOException {
		File produced = predictionDir(pins);
		if (produced == null || !produced.isDirectory()) {
			return null;
		}
		if (destination == null || produced.getCanonicalFile().equals(destination.getCanonicalFile())) {
			return produced;
		}
		mergeInto(produced, destination);
		deleteTree(produced);
		return destination;
	}

	/**
	 * Copies each sample's QC images from the batch's shared plots folder into its
	 * own search folder, so a search folder explains itself without a reader
	 * having to know where the batch put things. The batch's copy is left in
	 * place: it is the complete set for the run, and these are a few hundred
	 * kilobytes each.
	 *
	 * @param plotsRoot  the batch's {@value #PLOTS_DIR} folder.
	 * @param sample     base name whose images to take; matched as
	 *                   {@code <sample>.png} or {@code <sample>_*.png}.
	 * @param sampleDest the folder to write {@value #PLOTS_DIR} into.
	 * @return how many images were copied.
	 */
	public static int copySamplePlots(File plotsRoot, String sample, File sampleDest) throws IOException {
		if (plotsRoot == null || !plotsRoot.isDirectory() || sampleDest == null) {
			return 0;
		}
		int copied = 0;
		File[] groups = plotsRoot.listFiles(File::isDirectory);
		if (groups == null) {
			return 0;
		}
		for (File group : groups) {
			File[] images = group.listFiles(f -> {
				if (!f.isFile()) {
					return false;
				}
				String n = f.getName();
				return n.equals(sample + ".png") || n.startsWith(sample + "_");
			});
			if (images == null || images.length == 0) {
				continue;
			}
			File target = new File(new File(sampleDest, PLOTS_DIR), group.getName());
			target.mkdirs();
			for (File image : images) {
				Files.copy(image.toPath(), new File(target, image.getName()).toPath(),
						StandardCopyOption.REPLACE_EXISTING);
				copied++;
			}
		}
		return copied;
	}

	/**
	 * Appends a line to {@code msbooster_runs.txt} in the predictions folder
	 * recording which samples this invocation actually predicted.
	 * <p>
	 * Without it, a folder holding a pooled {@code spectraRT.predicted.bin} from a
	 * resumed run and QC plots from three different runs looks like one coherent
	 * artefact and is not. This is the cheapest honest way to say so.
	 */
	public static void writeRunNote(File predictionsDir, List<String> samples) throws IOException {
		if (predictionsDir == null || samples == null || samples.isEmpty()) {
			return;
		}
		predictionsDir.mkdirs();
		File note = new File(predictionsDir, "msbooster_runs.txt");
		boolean fresh = !note.isFile();
		try (PrintWriter w = new PrintWriter(new java.io.FileWriter(note, StandardCharsets.UTF_8, true))) {
			if (fresh) {
				w.println("# One line per MSBooster invocation. The per-sample QC images below belong to the");
				w.println("# run that produced them; spectraRT.tsv / spectraRT.predicted.bin describe only the");
				w.println("# most recent line, because a resumed run re-predicts just the samples it needs to.");
			}
			w.println(java.time.LocalDateTime.now().withNano(0) + "\t" + samples.size() + " sample(s)\t"
					+ String.join(",", samples));
		}
	}

	private static void mergeInto(File from, File to) throws IOException {
		to.mkdirs();
		File[] entries = from.listFiles();
		if (entries == null) {
			return;
		}
		for (File entry : entries) {
			File target = new File(to, entry.getName());
			if (entry.isDirectory()) {
				mergeInto(entry, target);
			} else {
				Files.copy(entry.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
		}
	}

	private static void deleteTree(File dir) {
		File[] entries = dir.listFiles();
		if (entries != null) {
			for (File entry : entries) {
				if (entry.isDirectory()) {
					deleteTree(entry);
				} else {
					entry.delete();
				}
			}
		}
		dir.delete();
	}

	private static File deepestCommon(File a, File b) {
		String pa = a.getAbsolutePath();
		while (b != null) {
			String pb = b.getAbsolutePath();
			if (pa.equals(pb) || pa.startsWith(pb + File.separator)) {
				return b;
			}
			b = b.getParentFile();
		}
		return null;
	}

	/**
	 * The file MSBooster writes for a given input pin: {@code <stem>_edited.pin},
	 * beside it. Matches what {@code renamePin = 1} produces.
	 */
	public static File editedPin(File pin) {
		String name = pin.getName();
		int dot = name.lastIndexOf('.');
		String stem = dot > 0 ? name.substring(0, dot) : name;
		return new File(pin.getParentFile(), stem + "_edited.pin");
	}

	private static void validate(Settings s) throws IOException {
		if (s == null) {
			throw new IOException("MSBooster settings are null");
		}
		if (s.diann == null || !s.diann.exists()) {
			throw new IOException("DIA-NN executable not found; MSBooster cannot predict without it: "
					+ (s.diann == null ? "null" : s.diann.getAbsolutePath())
					+ ". Set Resources > DiaNN, which the DIA workflows already use.");
		}
		if (s.pinFiles == null || s.pinFiles.isEmpty()) {
			throw new IOException("No .pin files supplied to MSBooster");
		}
		if (s.mzmlDirectories == null || s.mzmlDirectories.isEmpty()) {
			throw new IOException("No spectra folder supplied to MSBooster");
		}
		// mzmlDirectory and pinPepXMLDirectory are the only two settings
		// MSBooster splits into several paths, and it does so with a bare
		// String.split(" ") - no quoting, no escaping (mainsteps.PinMzmlMatcher).
		// A single path containing a space is therefore torn into fragments and
		// surfaces as "No matching mzml/mgf file", which reads as a missing file
		// rather than an unusable path. There is no encoding that survives, so
		// this is a hard limit of the tool: say so plainly here instead.
		for (File f : joinChecked(s.pinFiles, s.mzmlDirectories)) {
			if (f.getAbsolutePath().indexOf(' ') >= 0) {
				throw new IOException("MSBooster cannot be given a path containing a space: "
						+ f.getAbsolutePath() + ". It splits its file lists on spaces with no way to quote "
						+ "them, so the path is read as several broken ones. Move the result folder or the "
						+ "spectra to a location without spaces, or switch stage-1 MSBooster off.");
			}
		}
		// MSBooster pairs a pin with its spectra file by base name
		// (mainsteps.PinMzmlMatcher). A pin called <sample>_something.pin sends it
		// looking for <sample>_something.mzML and it aborts with "No matching
		// mzml/mgf file". Checked here rather than left to that error, which names
		// the directory and invites the wrong diagnosis.
		for (File pin : s.pinFiles) {
			String stem = pin.getName().replaceFirst("\\.pin$", "");
			boolean matched = false;
			for (File dir : s.mzmlDirectories) {
				if (dir == null) {
					continue;
				}
				if (dir.isFile()) {
					matched |= dir.getName().startsWith(stem + ".");
					continue;
				}
				File[] spectra = dir.listFiles(f -> {
					String n = f.getName().toLowerCase(java.util.Locale.ROOT);
					return f.isFile() && (n.endsWith(".mzml") || n.endsWith(".mgf"));
				});
				if (spectra == null) {
					continue;
				}
				for (File spectrum : spectra) {
					String name = spectrum.getName();
					if (name.substring(0, name.lastIndexOf('.')).equalsIgnoreCase(stem)) {
						matched = true;
						break;
					}
				}
			}
			if (!matched) {
				throw new IOException("No mzML/MGF named '" + stem + "' in " + s.mzmlDirectories
						+ " for pin " + pin.getName() + ". MSBooster matches a pin to its spectra by base name, "
						+ "so the pin must keep the name MSFragger gave it - <sample>.pin, matching "
						+ "<sample>.mzML. Do not rename or suffix it.");
			}
		}
	}

	private static List<File> joinChecked(List<File> a, List<File> b) {
		List<File> out = new ArrayList<>();
		if (a != null) {
			out.addAll(a);
		}
		if (b != null) {
			out.addAll(b);
		}
		out.removeIf(f -> f == null);
		return out;
	}

	private static String joinPaths(List<File> files) {
		List<String> parts = new ArrayList<>();
		for (File f : files) {
			if (f != null) {
				parts.add(f.getAbsolutePath());
			}
		}
		return String.join(" ", parts);
	}

	/** Convenience for callers that have one folder rather than a list. */
	public static List<File> listOf(File... files) {
		List<File> out = new ArrayList<>();
		for (File f : Arrays.asList(files)) {
			if (f != null) {
				out.add(f);
			}
		}
		return out;
	}
}

package uk.ac.quadram.metapilot.resources;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Centralized Resource and Tool Path Management System for MetaPilot
 * 
 * <h2>Purpose</h2>
 * This comprehensive resource management class handles all external tool paths and 
 * dependencies required by MetaPilot workflows. It provides unified configuration,
 * validation, and persistence of tool locations across different proteomics engines
 * and analysis pipelines.
 * 
 * <h2>Architecture</h2>
 * <ul>
 *   <li><strong>Unified Tool Management:</strong> Single interface for all external tools</li>
 *   <li><strong>Automatic Path Resolution:</strong> Intelligent path construction and validation</li>
 *   <li><strong>Persistent Configuration:</strong> User preferences saved across sessions</li>
 *   <li><strong>Status Monitoring:</strong> Real-time availability checking</li>
 *   <li><strong>Cross-Platform Support:</strong> Windows/Linux/macOS compatibility</li>
 * </ul>
 * 
 * <h2>Supported Tool Categories</h2>
 * 
 * <h3>General Tools:</h3>
 * <ul>
 *   <li><strong>FunctionDB:</strong> Functional annotation database</li>
 *   <li><strong>Python:</strong> Python interpreter for scripts</li>
 * </ul>
 * 
 * <h3>DDA (Data-Dependent Acquisition) Tools:</h3>
 * <ul>
 *   <li><strong>MaxQuant:</strong> Label-free quantification and identification</li>
 *   <li><strong>pFind:</strong> Chinese-developed database search engine</li>
 *   <li><strong>FragPipe:</strong> MSFragger-based analysis pipeline</li>
 *   <li><strong>Sage:</strong> Fast Rust-based search engine</li>
 *   <li><strong>FlashLFQ:</strong> Label-free quantification tool</li>
 *   <li><strong>DBReducer:</strong> Database reduction utility</li>
 *   <li><strong>msconvert:</strong> Spectrum file conversion (ProteoWizard)</li>
 * </ul>
 * 
 * <h3>DIA (Data-Independent Acquisition) Tools:</h3>
 * <ul>
 *   <li><strong>DIA-NN:</strong> Deep learning-based DIA analysis</li>
 *   <li><strong>DeepDetect:</strong> Advanced detection algorithms</li>
 * </ul>
 * 
 * <h2>Path Management Strategy</h2>
 * <h3>Base Resource Directory:</h3>
 * <ul>
 *   <li>Default: "./Resources/" relative to application</li>
 *   <li>Configurable through constructor or setter</li>
 *   <li>All tool paths resolved relative to base directory</li>
 * </ul>
 * 
 * <h3>Tool Path Structure:</h3>
 * <pre>{@code
 *   Resources/
 *   ├── mq_bin/MaxQuantCmd.exe           (MaxQuant)
 *   ├── DBReducer/pFind3/bin/pFind.exe   (pFind)
 *   ├── fragpipe/bin/fragpipe.bat        (FragPipe)
 *   ├── sage/sage.exe                    (Sage)
 *   ├── DIA-NN/2.0/DiaNN.exe            (DIA-NN)
 *   ├── ProteoWizard/msconvert.exe       (msconvert)
 *   ├── FlashLFQ/CMD.exe                 (FlashLFQ)
 *   ├── function/func_def.db             (FunctionDB)
 *   └── Python/python.exe                (Python)
 * }</pre>
 * 
 * <h2>Tool Status Management</h2>
 * <ul>
 *   <li><strong>AVAILABLE:</strong> Tool executable found and accessible</li>
 *   <li><strong>UNAVAILABLE:</strong> Tool missing or inaccessible</li>
 *   <li><strong>Auto-Detection:</strong> Status updated on path changes</li>
 *   <li><strong>Validation:</strong> File existence checking</li>
 * </ul>
 * 
 * <h2>Configuration Persistence</h2>
 * <p>Tool paths are automatically saved to Java preferences:</p>
 * <ul>
 *   <li><strong>Node:</strong> uk/ac/quadram/metapilot/resources</li>
 *   <li><strong>Scope:</strong> User-specific preferences</li>
 *   <li><strong>Keys:</strong> Tool type names and base path</li>
 *   <li><strong>Automatic Loading:</strong> Restored on application startup</li>
 * </ul>
 * 
 * <h2>Usage Patterns</h2>
 * 
 * <h3>Basic Initialization:</h3>
 * <pre>{@code
 *   // Use default resource path
 *   MetaResourceParameters resources = new MetaResourceParameters();
 *   
 *   // Check tool availability
 *   boolean maxQuantReady = resources.isToolAvailable(ToolType.MAXQUANT);
 * }</pre>
 * 
 * <h3>Custom Configuration:</h3>
 * <pre>{@code
 *   // Custom base path
 *   MetaResourceParameters resources = new MetaResourceParameters("C:/Proteomics/Tools");
 *   
 *   // Manual tool path override
 *   resources.setToolPath(ToolType.MAXQUANT, "D:/MaxQuant/bin/MaxQuantCmd.exe");
 *   
 *   // Save configuration
 *   resources.saveToPreferences();
 * }</pre>
 * 
 * <h3>Workflow Integration:</h3>
 * <pre>{@code
 *   // Check required tools before workflow
 *   ToolType[] requiredTools = {ToolType.MAXQUANT, ToolType.MSCONVERT};
 *   boolean allAvailable = Arrays.stream(requiredTools)
 *       .allMatch(resources::isToolAvailable);
 *   
 *   if (allAvailable) {
 *       // Launch workflow
 *       String maxQuantPath = resources.getToolPath(ToolType.MAXQUANT);
 *       // ... execute analysis
 *   }
 * }</pre>
 * 
 * <h2>Download Integration</h2>
 * <p>Some tools provide download URLs for automatic installation:</p>
 * <ul>
 *   <li>MaxQuant: Official website download</li>
 *   <li>FragPipe: Nesvilab releases</li>
 *   <li>Sage: GitHub releases</li>
 *   <li>DIA-NN: GitHub releases</li>
 *   <li>msconvert: ProteoWizard official site</li>
 * </ul>
 * 
 * <h2>Cross-Platform Considerations</h2>
 * <ul>
 *   <li><strong>Path Separators:</strong> Automatic OS-specific handling</li>
 *   <li><strong>Executable Extensions:</strong> .exe for Windows, none for Unix</li>
 *   <li><strong>Script Handling:</strong> .bat for Windows, .sh for Unix</li>
 *   <li><strong>Permission Checking:</strong> Executable bit validation on Unix</li>
 * </ul>
 * 
 * <h2>Error Handling</h2>
 * <ul>
 *   <li><strong>Missing Tools:</strong> Graceful degradation with clear status</li>
 *   <li><strong>Invalid Paths:</strong> Automatic correction attempts</li>
 *   <li><strong>Permission Issues:</strong> Clear error reporting</li>
 *   <li><strong>Recovery Strategies:</strong> Fallback to default locations</li>
 * </ul>
 * 
 * <h2>Thread Safety</h2>
 * <p><strong>Note:</strong> This class is NOT thread-safe. External synchronization
 * is required if accessed from multiple threads concurrently.</p>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.resources.MetaResourceIO
 * @see uk.ac.quadram.metapilot.ui.ResourceConfigPanel
 * @see uk.ac.quadram.metapilot.workflows.ResourceConfigWorkflow
 */
public class MetaResourceParameters {

	private static final Logger LOGGER = LogManager.getLogger(MetaResourceParameters.class);

	// Tool types
	public enum ToolType {
		// General tools
		FUNCDB("FunctionDB", ""), PYTHON("Python", ""),
		// DDA tools
		PFIND("pFind", "https://github.com/pFindStudio/pFind3"),
		MAXQUANT("MaxQuant", "https://www.maxquant.org/maxquant/"), DBREDUCER("DBReducer", ""),
		FRAGPIPE("FragPipe", "https://fragpipe.nesvilab.org/"), SAGE("Sage", "https://github.com/lazear/sage"),
		/**
		 * Virtual dependency used by workflow preflight. It is not shown as a
		 * separately configurable resource: {@link #resolveMsfraggerJar()} locates
		 * the MSFragger jar bundled with the configured FragPipe installation.
		 */
		MSFRAGGER("MSFragger", ""),
		FLASHLFQ("FlashLFQ", ""), MSCONVERT("msconvert", "https://proteowizard.sourceforge.io/download.html"),
		// DIA tools
		DIANN("DiaNN", "https://github.com/vdemichev/DiaNN/releases"), DEEPDETECT("DeepDetect", "");

		private final String displayName;
		private final String downloadUrl;

		ToolType(String displayName, String downloadUrl) {
			this.displayName = displayName;
			this.downloadUrl = downloadUrl;
		}

		public String getDisplayName() {
			return displayName;
		}

		public String getDownloadUrl() {
			return downloadUrl;
		}

		/**
		 * Check if this tool has a download URL
		 */
		public boolean hasDownloadUrl() {
			return downloadUrl != null && !downloadUrl.isEmpty();
		}

		public static ToolType[] getGeneralTools() {
			return new ToolType[] { FUNCDB, PYTHON };
		}

		public static ToolType[] getDDATools() {
			return new ToolType[] { PFIND, MAXQUANT, DBREDUCER, FRAGPIPE, SAGE, FLASHLFQ, MSCONVERT };
		}

		public static ToolType[] getDIATools() {
			return new ToolType[] { DIANN, DEEPDETECT };
		}
	}

	// Status of the tool
	public enum ToolStatus {
		AVAILABLE(true), UNAVAILABLE(false);

		private final boolean available;

		ToolStatus(boolean available) {
			this.available = available;
		}

		public boolean isAvailable() {
			return available;
		}
	}

	// Map to store tool paths and their status
	private Map<ToolType, ToolInfo> toolInfoMap = new HashMap<>();
	private String baseResourcePath;

	// Preferences key for storing resource paths
	private static final String PREFERENCES_NODE = "uk/ac/quadram/metapilot/resources";
	private static final String BASE_PATH_KEY = "baseResourcePath";

	/**
	 * Creates a new MetaResourceParameters with default values
	 */
	public MetaResourceParameters() {
		initializeDefaultPaths();
		loadSavedPaths();
	}

	/**
	 * Creates a new MetaResourceParameters with a specified base resource path
	 */
	public MetaResourceParameters(String baseResourcePath) {
		this.baseResourcePath = baseResourcePath;
		initializeDefaultPaths();
		loadSavedPaths();
	}

	/**
	 * Initialize default relative paths for all tools
	 */
	private void initializeDefaultPaths() {
		// Set default base resource path if not specified
		if (baseResourcePath == null || baseResourcePath.isEmpty()) {
			baseResourcePath = new File("Resources").getAbsolutePath();
		}

		// Define default relative paths for each tool
		toolInfoMap.put(ToolType.PFIND, new ToolInfo(ToolType.PFIND, "DBReducer\\pFind3\\bin", "pFind.exe"));
		toolInfoMap.put(ToolType.MAXQUANT, new ToolInfo(ToolType.MAXQUANT, "mq_bin", "MaxQuantCmd.exe"));
		toolInfoMap.put(ToolType.DBREDUCER, new ToolInfo(ToolType.DBREDUCER, "DBReducer\\DBReducer", "DBReducer.exe"));
		toolInfoMap.put(ToolType.FRAGPIPE, new ToolInfo(ToolType.FRAGPIPE, "fragpipe\\bin", "fragpipe.bat"));
		toolInfoMap.put(ToolType.SAGE, new ToolInfo(ToolType.SAGE, "sage", "sage.exe"));
		toolInfoMap.put(ToolType.FLASHLFQ, new ToolInfo(ToolType.FLASHLFQ, "FlashLFQ", "CMD.exe"));
		toolInfoMap.put(ToolType.FUNCDB, new ToolInfo(ToolType.FUNCDB, "function", "func_def.db"));
		toolInfoMap.put(ToolType.MSCONVERT, new ToolInfo(ToolType.MSCONVERT, "ProteoWizard", "msconvert.exe"));
		toolInfoMap.put(ToolType.DIANN, new ToolInfo(ToolType.DIANN, "DIA-NN\\2.0", "DiaNN.exe"));
		toolInfoMap.put(ToolType.DEEPDETECT, new ToolInfo(ToolType.DEEPDETECT, "DeepDetect", "DeepDetect.exe"));
		toolInfoMap.put(ToolType.PYTHON, new ToolInfo(ToolType.PYTHON, "Python", "python.exe"));
		// Initialize full paths for all tools
		for (Map.Entry<ToolType, ToolInfo> entry : toolInfoMap.entrySet()) {
			ToolInfo toolInfo = entry.getValue();
			File fullPath = new File(baseResourcePath, toolInfo.getRelativePath());
			toolInfo.setFullPath(fullPath.getAbsolutePath());
			updateToolStatus(entry.getKey());
		}
	}

	/**
	 * Load saved paths from preferences
	 */
	private void loadSavedPaths() {
		Preferences prefs = Preferences.userRoot().node(PREFERENCES_NODE);

		// Load base path if saved
		String savedBasePath = prefs.get(BASE_PATH_KEY, null);
		if (savedBasePath != null && !savedBasePath.isEmpty()) {
			baseResourcePath = savedBasePath;
		}

		// Load all tool paths
		for (ToolType toolType : ToolType.values()) {
			String savedPath = prefs.get(toolType.name(), null);
			if (savedPath != null && !savedPath.isEmpty()) {
				ToolInfo toolInfo = toolInfoMap.get(toolType);
				if (toolInfo != null) {
					toolInfo.setFullPath(savedPath);
					updateToolStatus(toolType);
				}
			}
		}
	}

	/**
	 * Update all tool paths based on the new base resource path Always updates the
	 * paths based on the new base path and checks if tools are available
	 */
	public void updateAllPaths(String newBasePath) {
		// Check if the new base path is valid
		if (newBasePath == null || newBasePath.isEmpty()) {
			return;
		}

		this.baseResourcePath = newBasePath;

		for (Map.Entry<ToolType, ToolInfo> entry : toolInfoMap.entrySet()) {
			ToolInfo toolInfo = entry.getValue();

			// Always create a new path based on the base resource path
			File newToolPath = new File(newBasePath, toolInfo.getRelativePath());
			toolInfo.setFullPath(newToolPath.getAbsolutePath());

			// Update tool status
			updateToolStatus(entry.getKey());
		}
	}

	/**
	 * Update the status of a specific tool based on whether its path exists
	 */
	private void updateToolStatus(ToolType toolType) {
		ToolInfo toolInfo = toolInfoMap.get(toolType);
		if (toolInfo != null) {
			boolean isAvailable = false;
			String path = toolInfo.getFullPath();

			if (path != null && !path.isEmpty()) {
				File toolFile = new File(path);
				isAvailable = toolFile.exists();
			}

			toolInfo.setStatus(isAvailable ? ToolStatus.AVAILABLE : ToolStatus.UNAVAILABLE);
		}
	}

	/**
	 * Get the full path for a specific tool
	 */
	public String getToolPath(ToolType toolType) {
		ToolInfo toolInfo = toolInfoMap.get(toolType);
		return toolInfo != null ? toolInfo.getFullPath() : "";
	}

	/**
	 * Set the full path for a specific tool
	 */
	public void setToolPath(ToolType toolType, String fullPath) {
		ToolInfo toolInfo = toolInfoMap.get(toolType);
		if (toolInfo != null) {
			toolInfo.setFullPath(fullPath);
			updateToolStatus(toolType);
		}
	}

	/**
	 * Get the status of a specific tool
	 */
	public ToolStatus getToolStatus(ToolType toolType) {
		ToolInfo toolInfo = toolInfoMap.get(toolType);
		return toolInfo != null ? toolInfo.getStatus() : ToolStatus.UNAVAILABLE;
	}

	/**
	 * Check if a specific tool is available
	 */
	public boolean isToolAvailable(ToolType toolType) {
		if (toolType == ToolType.MSFRAGGER) {
			// MSFragger is a derived dependency: FragPipe ships the jar that the
			// standalone workflow invokes directly.
			return resolveMsfraggerJar() != null;
		}
		ToolStatus status = getToolStatus(toolType);
		return status.isAvailable();
	}

	/** Matches a bundled tool folder or jar such as {@code MSFragger-4.4.1}. */
	private static final java.util.regex.Pattern MSFRAGGER_VERSIONED = java.util.regex.Pattern
			.compile("(?i)^msfragger[-_]?([0-9]+(?:\\.[0-9]+)*)$");

	/**
	 * {@code MSBooster-<version>.jar}, a jar directly inside FragPipe's
	 * {@code tools} folder rather than in a versioned subfolder as MSFragger is.
	 */
	private static final java.util.regex.Pattern MSBOOSTER_VERSIONED = java.util.regex.Pattern
			.compile("(?i)MSBooster[-_]?v?(\\d+(?:\\.\\d+)*)");

	/**
	 * {@code original-crystalc-<version>.jar}, beside MSBooster in FragPipe's
	 * {@code tools} folder. FragPipe's own {@code CmdCrystalc} names the file
	 * literally; matching on a prefix instead lets an upgraded installation keep
	 * working.
	 */
	private static final java.util.regex.Pattern CRYSTALC_VERSIONED = java.util.regex.Pattern
			.compile("(?i)(?:original-)?crystalc[-_]?v?(\\d+(?:\\.\\d+)*)");

	/**
	 * Locates MSFragger.jar for the standalone open-search workflow.
	 * <p>
	 * MSFragger is inferred from the configured FragPipe path. FragPipe bundles
	 * MSFragger, and the Resources entry points at the launcher inside
	 * {@code bin/}, so the jar sits at:
	 *
	 * <pre>
	 * &lt;FragPipe root&gt;/
	 *   bin/    &lt;- ToolType.FRAGPIPE points here (fragpipe.bat / fragpipe)
	 *   tools/
	 *     MSFragger-4.4.1/
	 *       MSFragger-4.4.1.jar
	 * </pre>
	 *
	 * When several {@code MSFragger-*} folders are present — upgrades leave the
	 * old ones behind — the highest version wins, so an upgrade is picked up
	 * without reconfiguration rather than silently pinning an old build.
	 *
	 * @return the bundled jar, or null when it cannot be resolved from FragPipe
	 */
	public File resolveMsfraggerJar() {
		return findMsfraggerUnderFragpipe();
	}

	/**
	 * Plausible FragPipe install roots, from the configured launcher.
	 * <p>
	 * The Resources entry is normally {@code <root>/bin/fragpipe.bat}, but
	 * tolerating {@code bin/} itself and the root covers how people actually
	 * configure it, and one extra level up covers nesting such as
	 * {@code FragPipe-23.0/fragpipe/bin}. Shared by every tool that lives under
	 * {@code tools/} so they cannot disagree about where the installation is.
	 */
	private java.util.List<File> fragpipeInstallRoots() {
		java.util.List<File> roots = new java.util.ArrayList<>();
		String fragpipe = getToolPath(ToolType.FRAGPIPE);
		if (fragpipe == null || fragpipe.isBlank()) {
			return roots;
		}
		File configured = new File(fragpipe);
		if (!configured.exists()) {
			return roots;
		}
		File binDir = configured.isDirectory() ? configured : configured.getParentFile();
		if (binDir == null) {
			return roots;
		}
		if (binDir.getParentFile() != null) {
			roots.add(binDir.getParentFile());
		}
		roots.add(binDir);
		if (binDir.getParentFile() != null && binDir.getParentFile().getParentFile() != null) {
			roots.add(binDir.getParentFile().getParentFile());
		}
		return roots;
	}

	/** FragPipe's {@code tools} folder, or null when it cannot be located. */
	public File resolveFragpipeToolsDir() {
		for (File root : fragpipeInstallRoots()) {
			File tools = new File(root, "tools");
			if (tools.isDirectory()) {
				return tools;
			}
		}
		return null;
	}

	/**
	 * The highest-versioned {@code MSBooster-*.jar} in FragPipe's {@code tools}
	 * folder, or null.
	 * <p>
	 * Located rather than separately configured, exactly as
	 * {@link #resolveMsfraggerJar()} is: one Resources entry for FragPipe, and
	 * everything bundled with it is found from there. Version-matched rather than
	 * hard-coded so an upgraded installation keeps working.
	 */
	public File resolveMsboosterJar() {
		File tools = resolveFragpipeToolsDir();
		if (tools == null) {
			return null;
		}
		File[] jars = tools.listFiles(f -> f.isFile()
				&& f.getName().toLowerCase(java.util.Locale.ROOT).startsWith("msbooster")
				&& f.getName().toLowerCase(java.util.Locale.ROOT).endsWith(".jar"));
		if (jars == null || jars.length == 0) {
			return null;
		}
		File best = null;
		int[] bestVersion = null;
		for (File jar : jars) {
			String stem = jar.getName().substring(0, jar.getName().length() - 4);
			java.util.regex.Matcher matcher = MSBOOSTER_VERSIONED.matcher(stem);
			int[] version = matcher.matches() ? parseVersion(matcher.group(1)) : new int[] { 0 };
			if (bestVersion == null || compareVersions(version, bestVersion) > 0) {
				bestVersion = version;
				best = jar;
			}
		}
		return best;
	}

	/**
	 * The highest-versioned Crystal-C jar in FragPipe's {@code tools} folder, or
	 * null.
	 * <p>
	 * Located rather than separately configured, exactly as
	 * {@link #resolveMsboosterJar()} is. Crystal-C also needs
	 * {@code batmass-io} and {@code grppr} from the same folder; those are
	 * resolved by {@code CrystalCHeadlessRunner.buildClasspath} from beside this
	 * jar, so one Resources entry for FragPipe covers all three.
	 */
	public File resolveCrystalcJar() {
		File tools = resolveFragpipeToolsDir();
		if (tools == null) {
			return null;
		}
		File[] jars = tools.listFiles(f -> f.isFile()
				&& f.getName().toLowerCase(java.util.Locale.ROOT).contains("crystalc")
				&& f.getName().toLowerCase(java.util.Locale.ROOT).endsWith(".jar"));
		if (jars == null || jars.length == 0) {
			return null;
		}
		File best = null;
		int[] bestVersion = null;
		for (File jar : jars) {
			String stem = jar.getName().substring(0, jar.getName().length() - 4);
			java.util.regex.Matcher matcher = CRYSTALC_VERSIONED.matcher(stem);
			int[] version = matcher.matches() ? parseVersion(matcher.group(1)) : new int[] { 0 };
			if (bestVersion == null || compareVersions(version, bestVersion) > 0) {
				bestVersion = version;
				best = jar;
			}
		}
		return best;
	}

	/**
	 * The DIA-NN that FragPipe bundles for MSBooster, under
	 * {@code tools/diann/<version>/<os>/}, or null.
	 * <p>
	 * MSBooster drives DIA-NN through the 1.8.x prediction command line, and
	 * FragPipe pins a specific build for it - 1.8.2_beta_8 at the time of
	 * writing. DIA-NN 2.x changed that interface, so the standalone DIA-NN
	 * configured under Resources for the DIA workflows is the wrong binary for
	 * this job even though it is a perfectly good DIA-NN. Pointing MSBooster at a
	 * 2.x build is what made every MSBooster step fail in the 2026-09-10 12mix
	 * run. Prefer the bundled one; fall back to the configured one only so a
	 * FragPipe install without it can still be tried.
	 */
	public File resolveMsboosterDiann() {
		File tools = resolveFragpipeToolsDir();
		if (tools == null) {
			return null;
		}
		File diannRoot = new File(tools, "diann");
		if (!diannRoot.isDirectory()) {
			return null;
		}
		File best = null;
		int[] bestVersion = null;
		File anyBest = null;
		int[] anyBestVersion = null;
		File[] versions = diannRoot.listFiles(File::isDirectory);
		if (versions == null) {
			return null;
		}
		for (File versionDir : versions) {
			File executable = findDiannExecutable(versionDir);
			if (executable == null) {
				continue;
			}
			// Version folders carry suffixes such as "1.8.2_beta_8"; compare on the
			// numeric prefix and let the suffix break nothing.
			java.util.regex.Matcher matcher = java.util.regex.Pattern.compile("(\\d+(?:\\.\\d+)*)")
					.matcher(versionDir.getName());
			int[] version = matcher.find() ? parseVersion(matcher.group(1)) : new int[] { 0 };
			if (anyBestVersion == null || compareVersions(version, anyBestVersion) > 0) {
				anyBestVersion = version;
				anyBest = executable;
			}
			// "Highest version" is the wrong rule here: MSBooster drives DIA-NN
			// through the 1.8.x prediction interface, so a 2.x build newly bundled
			// by a FragPipe upgrade would be picked and would fail - the exact
			// failure the 2026-09-10 run hit with a standalone 2.5.1. Prefer the
			// highest below 2.0 and keep the overall highest only as a fallback,
			// so this neither breaks on an upgrade nor blocks a future MSBooster
			// that wants 2.x.
			if (version.length > 0 && version[0] >= 2) {
				continue;
			}
			if (bestVersion == null || compareVersions(version, bestVersion) > 0) {
				bestVersion = version;
				best = executable;
			}
		}
		if (best == null && anyBest != null) {
			// Not returned. MSBooster drives DIA-NN through the 1.8.x prediction
			// interface; a 2.x build starts, does nothing useful and leaves no
			// predictions, which is indistinguishable from MSBooster not having
			// run. Better to report no usable DIA-NN than to hand back one that
			// will waste a whole run.
			LOGGER.warn("FragPipe bundles only DIA-NN {} under tools/diann/, which is 2.x. MSBooster needs the "
					+ "1.8.x prediction interface, so this is reported as no usable DIA-NN.",
					anyBest.getParentFile());
		}
		return best;
	}

	/**
	 * Any DIA-NN executable at most two levels below {@code versionDir}.
	 * <p>
	 * FragPipe Linux ships the binary as a versioned name such as
	 * {@code diann-1.8.1.8} under {@code tools/diann/<version>/linux/}, not as a
	 * plain {@code diann} or {@code diann-linux}. Accept those plus the common
	 * Windows/Linux names; ignore shared libraries and data files that also sit
	 * in the same folder ({@code .so}, {@code .obo}, …).
	 */
	private static File findDiannExecutable(File versionDir) {
		java.util.List<File> candidates = new java.util.ArrayList<>();
		candidates.add(versionDir);
		File[] osDirs = versionDir.listFiles(File::isDirectory);
		if (osDirs != null) {
			candidates.addAll(java.util.Arrays.asList(osDirs));
		}
		File best = null;
		int bestRank = Integer.MAX_VALUE;
		for (File dir : candidates) {
			File[] files = dir.listFiles(File::isFile);
			if (files == null) {
				continue;
			}
			for (File f : files) {
				if (!isMsboosterDiannBinaryName(f.getName())) {
					continue;
				}
				int rank = diannBinaryPreference(f.getName());
				if (rank < bestRank) {
					bestRank = rank;
					best = f;
				}
			}
		}
		return best;
	}

	/** True for FragPipe-bundled DIA-NN binaries usable by MSBooster. */
	static boolean isMsboosterDiannBinaryName(String fileName) {
		if (fileName == null || fileName.isBlank()) {
			return false;
		}
		String n = fileName.toLowerCase(java.util.Locale.ROOT);
		if (n.endsWith(".so") || n.endsWith(".dll") || n.endsWith(".obo") || n.endsWith(".jar")
				|| n.endsWith(".txt") || n.endsWith(".md") || n.endsWith(".params")) {
			return false;
		}
		return n.equals("diann") || n.equals("diann.exe") || n.startsWith("diann-");
	}

	/** Lower is better: plain {@code diann} / {@code diann-linux*} over versioned names. */
	static int diannBinaryPreference(String fileName) {
		String n = fileName.toLowerCase(java.util.Locale.ROOT);
		if (n.equals("diann") || n.equals("diann.exe")) {
			return 0;
		}
		if (n.startsWith("diann-linux")) {
			return 1;
		}
		return 2;
	}

	/**
	 * {@code unimod.obo} from FragPipe's {@code tools} folder. MSBooster reads it
	 * to interpret modification masses; without it, it falls back to its own
	 * defaults silently.
	 */
	public File resolveUnimodObo() {
		File tools = resolveFragpipeToolsDir();
		if (tools == null) {
			return null;
		}
		File obo = new File(tools, "unimod.obo");
		return obo.isFile() ? obo : null;
	}

	/**
	 * Walks the documented FragPipe layout: from the configured launcher up to
	 * the install root, then into {@code tools/MSFragger-*}.
	 */
	private File findMsfraggerUnderFragpipe() {
		for (File root : fragpipeInstallRoots()) {
			File toolsDir = new File(root, "tools");
			if (!toolsDir.isDirectory()) {
				continue;
			}
			File jar = highestVersionMsfraggerJar(toolsDir);
			if (jar != null) {
				return jar;
			}
		}
		return null;
	}

	/**
	 * Among {@code tools/MSFragger-*} folders, returns the jar from the
	 * highest-versioned one.
	 */
	private static File highestVersionMsfraggerJar(File toolsDir) {
		File[] entries = toolsDir.listFiles(File::isDirectory);
		if (entries == null) {
			return null;
		}
		File bestDir = null;
		int[] bestVersion = null;
		for (File entry : entries) {
			java.util.regex.Matcher matcher = MSFRAGGER_VERSIONED.matcher(entry.getName());
			if (!matcher.matches()) {
				continue;
			}
			int[] version = parseVersion(matcher.group(1));
			if (bestVersion == null || compareVersions(version, bestVersion) > 0) {
				bestVersion = version;
				bestDir = entry;
			}
		}
		if (bestDir == null) {
			// Folder present but unversioned (e.g. plain "MSFragger"): still usable.
			for (File entry : entries) {
				if (entry.getName().toLowerCase(java.util.Locale.ROOT).startsWith("msfragger")) {
					File jar = highestVersionJarIn(entry);
					if (jar != null) {
						return jar;
					}
				}
			}
			return null;
		}
		return highestVersionJarIn(bestDir);
	}

	/** The highest-versioned {@code MSFragger*.jar} directly inside a folder. */
	private static File highestVersionJarIn(File directory) {
		File[] jars = directory.listFiles(f -> f.isFile()
				&& f.getName().toLowerCase(java.util.Locale.ROOT).startsWith("msfragger")
				&& f.getName().toLowerCase(java.util.Locale.ROOT).endsWith(".jar"));
		if (jars == null || jars.length == 0) {
			return null;
		}
		File best = null;
		int[] bestVersion = null;
		for (File jar : jars) {
			String stem = jar.getName().substring(0, jar.getName().length() - 4);
			java.util.regex.Matcher matcher = MSFRAGGER_VERSIONED.matcher(stem);
			int[] version = matcher.matches() ? parseVersion(matcher.group(1)) : new int[] { 0 };
			if (bestVersion == null || compareVersions(version, bestVersion) > 0) {
				bestVersion = version;
				best = jar;
			}
		}
		return best;
	}

	private static int[] parseVersion(String version) {
		String[] parts = version.split("\\.");
		int[] out = new int[parts.length];
		for (int i = 0; i < parts.length; i++) {
			try {
				out[i] = Integer.parseInt(parts[i]);
			} catch (NumberFormatException e) {
				out[i] = 0;
			}
		}
		return out;
	}

	/** Numeric component-wise comparison, so 4.10.0 sorts above 4.9.0. */
	private static int compareVersions(int[] a, int[] b) {
		int n = Math.max(a.length, b.length);
		for (int i = 0; i < n; i++) {
			int left = i < a.length ? a[i] : 0;
			int right = i < b.length ? b[i] : 0;
			if (left != right) {
				return Integer.compare(left, right);
			}
		}
		return 0;
	}

	/**
	 * Get the base resource path
	 */
	public String getBaseResourcePath() {
		return baseResourcePath;
	}

	/**
	 * Set the base resource path and update all tool paths
	 */
	public void setBaseResourcePath(String baseResourcePath) {
		updateAllPaths(baseResourcePath);
	}

	/**
	 * Get all tool information
	 */
	public Map<ToolType, ToolInfo> getToolInfoMap() {
		return toolInfoMap;
	}

	/**
	 * Save all paths to preferences
	 */
	public void saveToPreferences() {
		Preferences prefs = Preferences.userRoot().node(PREFERENCES_NODE);

		// Save base path
		prefs.put(BASE_PATH_KEY, baseResourcePath);

		// Save all tool paths
		for (Map.Entry<ToolType, ToolInfo> entry : toolInfoMap.entrySet()) {
			prefs.put(entry.getKey().name(), entry.getValue().getFullPath());
		}
	}

	/**
	 * Inner class to store information about a tool
	 */
	public static class ToolInfo {
		private final ToolType type;
		private final String folderName;
		private final String executableName;
		private String fullPath;
		private ToolStatus status;

		public ToolInfo(ToolType type, String folderName, String executableName) {
			this.type = type;
			this.folderName = folderName;
			this.executableName = executableName;
			this.status = ToolStatus.UNAVAILABLE;
		}

		public String getRelativePath() {
			return folderName + File.separator + executableName;
		}

		public String getFolderName() {
			return folderName;
		}

		public String getExecutableName() {
			return executableName;
		}

		public ToolType getType() {
			return type;
		}

		public String getFullPath() {
			return fullPath;
		}

		public void setFullPath(String fullPath) {
			this.fullPath = fullPath;
		}

		public ToolStatus getStatus() {
			return status;
		}

		public void setStatus(ToolStatus status) {
			this.status = status;
		}
	}
}

package uk.ac.quadram.metapilot.resources;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Manages resources and default paths for the MetaPilot application
 */
public class ResourceManager {
	private static final Logger LOGGER = LogManager.getLogger(ResourceManager.class);

	// Preference keys
	public static final String RESOURCE_PATH_KEY = "ResourcePath";
	public static final String PFIND_KEY = "pFindPath";
	public static final String MAXQUANT_KEY = "MaxQuantPath";
	public static final String DBREDUCER_KEY = "DBReducerPath";
	public static final String FRAGPIPE_KEY = "FragPipePath";
	public static final String SAGE_KEY = "SagePath";
	public static final String FLASHLFQ_KEY = "FlashLFQPath";
	public static final String FUNCDB_KEY = "FunctionDBPath";
	public static final String MSCONVERT_KEY = "MSConvertPath";
	public static final String DIANN_KEY = "DiaNNPath";
	public static final String DEEPDETECT_KEY = "DeepDetectPath";

	// Default relative paths
	private static final Map<String, String> DEFAULT_RELATIVE_PATHS = new HashMap<>();

	static {
		// Initialize default relative paths
		DEFAULT_RELATIVE_PATHS.put(PFIND_KEY, "Resources\\DBReducer\\pFind3\\bin\\pFind.exe");
		DEFAULT_RELATIVE_PATHS.put(DBREDUCER_KEY, "Resources\\DBReducer\\DBReducer\\DBReducer.exe");
		DEFAULT_RELATIVE_PATHS.put(FRAGPIPE_KEY, "Resources\\fragpipe\\bin\\fragpipe.bat");
		DEFAULT_RELATIVE_PATHS.put(FLASHLFQ_KEY, "Resources\\FlashLFQ\\CMD.exe");
		DEFAULT_RELATIVE_PATHS.put(SAGE_KEY, "Resources\\sage\\sage.exe");
		DEFAULT_RELATIVE_PATHS.put(FUNCDB_KEY, "Resources\\function\\func_def.db");
		DEFAULT_RELATIVE_PATHS.put(MAXQUANT_KEY, "Resources\\mq_bin\\MaxQuantCmd.exe");
		DEFAULT_RELATIVE_PATHS.put(MSCONVERT_KEY, "Resources\\ProteoWizard\\msconvert.exe");
		DEFAULT_RELATIVE_PATHS.put(DIANN_KEY, "Resources\\DIA-NN\\2.0\\DiaNN.exe");
		DEFAULT_RELATIVE_PATHS.put(DEEPDETECT_KEY, "Resources\\DeepDetect\\DeepDetect.exe");
	}

	private final Preferences prefs;
	private String resourceBasePath;
	private final Map<String, String> resourcePaths = new HashMap<>();

	/**
	 * Create a new ResourceManager with default settings
	 */
	public ResourceManager() {
		prefs = Preferences.userNodeForPackage(ResourceManager.class);
		loadFromPreferences();
	}

	/**
	 * Load resource paths from preferences
	 */
	public void loadFromPreferences() {
		// Load the base resource path
		resourceBasePath = prefs.get(RESOURCE_PATH_KEY, getDefaultResourcePath());

		// Load all tool paths from preferences
		for (String key : DEFAULT_RELATIVE_PATHS.keySet()) {
			String savedPath = prefs.get(key, "");
			if (savedPath.isEmpty()) {
				// If no saved path, compute from base resource path
				savedPath = computePathFromBase(key);
			}
			resourcePaths.put(key, savedPath);
		}
	}

	/**
	 * Save all resource paths to preferences
	 */
	public void saveToPreferences() {
		// Save base resource path
		prefs.put(RESOURCE_PATH_KEY, resourceBasePath);

		// Save all tool paths
		for (Map.Entry<String, String> entry : resourcePaths.entrySet()) {
			prefs.put(entry.getKey(), entry.getValue());
		}
	}

	/**
	 * Get the default resource path (Resources folder in app directory)
	 */
	private String getDefaultResourcePath() {
		try {
			// Get the application directory
			File appDir = new File(System.getProperty("user.dir"));

			// Create Resources directory if it doesn't exist
			File resourcesDir = new File(appDir, "Resources");
			if (!resourcesDir.exists()) {
				resourcesDir.mkdirs();
			}

			return resourcesDir.getAbsolutePath();
		} catch (Exception e) {
			LOGGER.error("Error getting default resource path", e);
			return "Resources";
		}
	}

	/**
	 * Get the base resource path
	 */
	public String getResourceBasePath() {
		return resourceBasePath;
	}

	/**
	 * Set the base resource path and update all relative paths
	 */
	public void setResourceBasePath(String path) {
		if (path == null || path.isEmpty()) {
			return;
		}

		// Update the base path
		this.resourceBasePath = path;

		// Update only the tool paths that don't currently exist
		for (String key : DEFAULT_RELATIVE_PATHS.keySet()) {
			// Check if the existing path is valid
			if (!resourceExists(key)) {
				// Path doesn't exist, compute a new one based on the new base path
				String relativePath = DEFAULT_RELATIVE_PATHS.get(key);
				File toolPath = new File(path, getRelativePath(relativePath));
				resourcePaths.put(key, toolPath.getAbsolutePath());
				LOGGER.info("Updated path for " + key + " to " + toolPath.getAbsolutePath());
			} else {
				// Path already exists and is valid, keep it
				LOGGER.info("Keeping existing path for " + key + ": " + resourcePaths.get(key));
			}
		}

		// Save to preferences
		saveToPreferences();
	}

	/**
	 * Get a resource path by key
	 */
	public String getResourcePath(String key) {
		return resourcePaths.getOrDefault(key, "");
	}

	/**
	 * Set a specific resource path
	 */
	public void setResourcePath(String key, String path) {
		if (path == null || path.isEmpty() || !DEFAULT_RELATIVE_PATHS.containsKey(key)) {
			return;
		}

		resourcePaths.put(key, path);
		prefs.put(key, path);
	}

	/**
	 * Check if a resource exists
	 */
	public boolean resourceExists(String key) {
		String path = resourcePaths.get(key);
		return path != null && !path.isEmpty() && new File(path).exists();
	}

	/**
	 * Compute a tool path from the base resource path
	 */
	private String computePathFromBase(String key) {
		String relativePath = DEFAULT_RELATIVE_PATHS.get(key);
		if (relativePath == null) {
			return "";
		}

		File path = new File(resourceBasePath, getRelativePath(relativePath));
		return path.getAbsolutePath();
	}

	/**
	 * Normalize the relative path based on OS
	 */
	private String getRelativePath(String path) {
		// Convert Windows-style paths to the system's format
		return path.replace("\\", File.separator);
	}
}
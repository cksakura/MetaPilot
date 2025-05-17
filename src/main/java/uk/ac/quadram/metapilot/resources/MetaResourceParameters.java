package uk.ac.quadram.metapilot.resources;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;

/**
 * Comprehensive resource parameters class to handle all tool paths for
 * MetaPilot. This combines functionality from MetaSourcesMag and MetaSourcesDia
 * into a single class.
 */
public class MetaResourceParameters {

	// Tool types
	public enum ToolType {
		// DDA tools
		PFIND("pFind", ""),
		MAXQUANT("MaxQuant", "https://www.maxquant.org/maxquant/"),
		DBREDUCER("DBReducer", ""),
		
		// DIA tools
		FRAGPIPE("FragPipe", "https://fragpipe.nesvilab.org/"),
		SAGE("Sage", "https://github.com/lazear/sage"),
		FLASHLFQ("FlashLFQ", ""),
		FUNCDB("FunctionDB", ""),
		MSCONVERT("msconvert", "https://proteowizard.sourceforge.io/download.html"),
		DIANN("DiaNN", "https://github.com/vdemichev/DiaNN/releases"),
		DEEPDETECT("DeepDetect", "");
		
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
		
		public static ToolType[] getDDATools() {
			return new ToolType[] { PFIND, MAXQUANT, DBREDUCER };
		}
		
		public static ToolType[] getDIATools() {
			return new ToolType[] { FRAGPIPE, SAGE, FLASHLFQ, FUNCDB, MSCONVERT, DIANN, DEEPDETECT };
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
	 * Update all tool paths based on the new base resource path
	 * Always updates the paths based on the new base path and checks if tools are available
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
		ToolStatus status = getToolStatus(toolType);
		return status.isAvailable();
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
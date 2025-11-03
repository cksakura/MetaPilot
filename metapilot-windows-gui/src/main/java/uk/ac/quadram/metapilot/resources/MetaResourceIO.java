package uk.ac.quadram.metapilot.resources;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;

/**
 * IO handler for MetaResourceParameters Handles reading and writing parameter
 * files for all tools
 */
public class MetaResourceIO {

	private static final Logger logger = LogManager.getLogger(MetaResourceIO.class);

	// Default location for parameter files
	private static final String DEFAULT_CONFIG_DIRECTORY = "config";
	private static final String RESOURCE_PARAMS_FILENAME = "meta_resources.params";

	private File configFile;

	/**
	 * Create a new MetaResourceIO with the default config file location
	 */
	public MetaResourceIO() {
		File configDir = new File(DEFAULT_CONFIG_DIRECTORY);
		if (!configDir.exists()) {
			configDir.mkdirs();
		}

		this.configFile = new File(configDir, RESOURCE_PARAMS_FILENAME);
	}

	/**
	 * Create a new MetaResourceIO with a custom config file location
	 */
	public MetaResourceIO(File configFile) {
		this.configFile = configFile;

		// Create parent directory if it doesn't exist
		File parentDir = configFile.getParentFile();
		if (parentDir != null && !parentDir.exists()) {
			parentDir.mkdirs();
		}
	}

	/**
	 * Load resource parameters from the config file If the file doesn't exist,
	 * returns default parameters
	 */
	public MetaResourceParameters loadParameters() {
		// If the file doesn't exist, return default parameters
		if (!configFile.exists()) {
			logger.info("Config file doesn't exist, creating default parameters");
			return new MetaResourceParameters();
		}

		// Read parameters from file
		Map<String, String> paramMap = new HashMap<>();

		try (BufferedReader reader = new BufferedReader(new FileReader(configFile))) {
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}

				String[] parts = line.split("=", 2);
				if (parts.length == 2) {
					paramMap.put(parts[0].trim(), parts[1].trim());
				}
			}

			logger.info("Successfully loaded parameters from {}", configFile.getAbsolutePath());
		} catch (IOException e) {
			logger.error("Error reading config file: {}", e.getMessage());
			return new MetaResourceParameters();
		}

		// Create parameters from the map
		MetaResourceParameters params = new MetaResourceParameters();

		// Set base resource path if present
		String basePath = paramMap.get("BASE_PATH");
		if (basePath != null && !basePath.isEmpty()) {
			params.setBaseResourcePath(basePath);
		}

		// Set tool paths
		for (ToolType toolType : ToolType.values()) {
			String toolPath = paramMap.get(toolType.name());
			if (toolPath != null && !toolPath.isEmpty()) {
				params.setToolPath(toolType, toolPath);
			}
		}

		return params;
	}

	/**
	 * Save resource parameters to the config file
	 */
	public boolean saveParameters(MetaResourceParameters params) {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(configFile))) {
			// Write header
			writer.write("# MetaPilot Resource Configuration");
			writer.newLine();
			writer.write("# Generated on " + new java.util.Date());
			writer.newLine();
			writer.newLine();

			// Write base path
			writer.write("BASE_PATH=" + params.getBaseResourcePath());
			writer.newLine();
			writer.newLine();

			// Write tool paths
			writer.write("# Tool Paths");
			writer.newLine();

			for (Map.Entry<ToolType, MetaResourceParameters.ToolInfo> entry : params.getToolInfoMap().entrySet()) {

				ToolType toolType = entry.getKey();
				String path = entry.getValue().getFullPath();

				writer.write(toolType.name() + "=" + path);
				writer.newLine();
			}

			logger.info("Successfully saved parameters to {}", configFile.getAbsolutePath());
			return true;
		} catch (IOException e) {
			logger.error("Error writing config file: {}", e.getMessage());
			return false;
		}
	}

	/**
	 * Get the current config file
	 */
	public File getConfigFile() {
		return configFile;
	}

	/**
	 * Set a new config file
	 */
	public void setConfigFile(File configFile) {
		this.configFile = configFile;
	}
}
package uk.ac.quadram.metapilot.task.par;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 * Base class for parameter IO operations across all workflows
 * 
 * This class provides common functionality for saving and loading parameters
 * to/from JSON files, ensuring consistency across all workflow implementations.
 */
public abstract class BaseParameterIO {

	protected static final Logger LOGGER = LogManager.getLogger(BaseParameterIO.class);

	// Common file naming constants
	protected static final String PREFERENCES_DIR = ".metapilot";
	protected static final String CONFIG_DIR = "config";
	
	// Registry for workflow type to ParameterIO class mapping
	private static final Map<String, Class<? extends BaseParameterIO>> PARAMETER_IO_REGISTRY = new HashMap<>();
	
	/**
	 * Register a ParameterIO class for a specific workflow type
	 * This should be called during application initialization
	 */
	public static void registerParameterIO(String workflowType, Class<? extends BaseParameterIO> parameterIOClass) {
		PARAMETER_IO_REGISTRY.put(workflowType, parameterIOClass);
		LOGGER.debug("Registered ParameterIO class {} for workflow type: {}", parameterIOClass.getSimpleName(), workflowType);
	}
	
	/**
	 * Load parameters from any parameter file without knowing the type beforehand
	 * This is the main method for CLI usage: -par parameter_file_path
	 */
	public static Object loadParametersFromFilePath(String filePath) {
		try {
			System.out.println("DEBUG: loadParametersFromFilePath called with: " + filePath);
			File inputFile = new File(filePath);
			if (!inputFile.exists()) {
				LOGGER.error("Parameter file not found: " + filePath);
				return null;
			}
			
			// Read the file to detect workflow type
			JSONObject json = readJsonFromFileStatic(inputFile);
			String workflowType = json.optString("workflow", null);
			System.out.println("DEBUG: Found workflow type in JSON: '" + workflowType + "'");
			
			if (workflowType == null || workflowType.trim().isEmpty()) {
				LOGGER.error("No workflow type found in parameter file: " + filePath);
				return null;
			}
			
			// Get the appropriate ParameterIO class
			Class<? extends BaseParameterIO> parameterIOClass = PARAMETER_IO_REGISTRY.get(workflowType);
			System.out.println("DEBUG: Looking for ParameterIO class for workflow type: '" + workflowType + "'");
			System.out.println("DEBUG: Found ParameterIO class: " + (parameterIOClass != null ? parameterIOClass.getSimpleName() : "NULL"));
			System.out.println("DEBUG: Available registered workflow types: " + PARAMETER_IO_REGISTRY.keySet());
			
			if (parameterIOClass == null) {
				LOGGER.error("No ParameterIO class registered for workflow type: " + workflowType);
				LOGGER.error("Available workflow types: " + PARAMETER_IO_REGISTRY.keySet());
				return null;
			}
			
			// Create instance and load parameters
			BaseParameterIO parameterIO = parameterIOClass.getDeclaredConstructor().newInstance();
			System.out.println("DEBUG: Created ParameterIO instance: " + parameterIO.getClass().getSimpleName());
			Object result = parameterIO.loadParametersFromFile(inputFile);
			System.out.println("DEBUG: Loaded parameters object: " + (result != null ? result.getClass().getSimpleName() : "NULL"));
			return result;
			
		} catch (Exception e) {
			LOGGER.error("Failed to load parameters from file: " + filePath, e);
			System.out.println("DEBUG: Exception in loadParametersFromFilePath: " + e.getClass().getSimpleName() + " - " + e.getMessage());
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Static version of readJsonFromFile for use in factory method
	 */
	private static JSONObject readJsonFromFileStatic(File inputFile) {
		try (FileReader reader = new FileReader(inputFile)) {
			StringBuilder content = new StringBuilder();
			char[] buffer = new char[1024];
			int bytesRead;
			while ((bytesRead = reader.read(buffer)) != -1) {
				content.append(buffer, 0, bytesRead);
			}
			String jsonString = content.toString();
			if (jsonString.trim().isEmpty()) {
				throw new RuntimeException("Parameter file is empty");
			}
			return new JSONObject(jsonString);
		} catch (Exception e) {
			LOGGER.error("Failed to read parameters from file: " + e.getMessage(), e);
			throw new RuntimeException("Failed to read parameters from file", e);
		}
	}

	/**
	 * Get the workflow-specific analysis parameter filename
	 */
	protected abstract String getAnalysisParameterFileName();

	/**
	 * Get the workflow-specific last used parameters filename
	 */
	protected abstract String getLastUsedParametersFileName();

	/**
	 * Get the workflow-specific preferences filename
	 */
	protected abstract String getPreferencesFileName();

	/**
	 * Get the workflow type string for JSON headers
	 */
	protected abstract String getWorkflowType();

	/**
	 * Save analysis parameters to a specific output directory
	 */
	public boolean saveAnalysisParameters(String outputDirectory, Object parameters) {
		try {
			if (parameters == null || outputDirectory == null || outputDirectory.trim().isEmpty()) {
				LOGGER.warn("Cannot save analysis parameters: parameters or output directory is null/empty");
				return false;
			}

			File outputFile = new File(outputDirectory, getAnalysisParameterFileName());
			return saveParametersToFile(parameters, outputFile, "Analysis-specific parameters");

		} catch (Exception e) {
			LOGGER.error("Failed to save analysis parameters to " + outputDirectory, e);
			return false;
		}
	}

	/**
	 * Save parameters to an explicitly selected file.
	 */
	public boolean exportParametersToFile(Object parameters, File outputFile, String description) {
		try {
			if (parameters == null || outputFile == null) {
				LOGGER.warn("Cannot export parameters: parameters or output file is null");
				return false;
			}
			return saveParametersToFile(parameters, outputFile, description);
		} catch (Exception e) {
			LOGGER.error("Failed to export parameters to " + outputFile, e);
			return false;
		}
	}

	/**
	 * Save parameters as last used parameters
	 */
	public boolean saveLastUsedParameters(Object parameters) {
		try {
			if (parameters == null) {
				LOGGER.warn("Cannot save last used parameters: parameters is null");
				return false;
			}

			File outputFile = new File(CONFIG_DIR, getLastUsedParametersFileName());
			// Create config directory if it doesn't exist
			outputFile.getParentFile().mkdirs();

			return saveParametersToFile(parameters, outputFile, "Last used parameters");

		} catch (Exception e) {
			LOGGER.error("Failed to save last used parameters", e);
			return false;
		}
	}

	/**
	 * Save parameters to user preferences
	 */
	public boolean savePreferences(Object parameters) {
		try {
			if (parameters == null) {
				LOGGER.warn("Cannot save preferences: parameters is null");
				return false;
			}

			// Get user home directory for storing preferences
			String userHome = System.getProperty("user.home");
			File preferencesDir = new File(userHome, PREFERENCES_DIR);
			if (!preferencesDir.exists()) {
				preferencesDir.mkdirs();
			}

			File outputFile = new File(preferencesDir, getPreferencesFileName());
			return saveParametersToFile(parameters, outputFile, "User preferences");

		} catch (Exception e) {
			LOGGER.error("Failed to save preferences", e);
			return false;
		}
	}

	/**
	 * Load last used parameters
	 */
	public Object loadLastUsedParameters() {
		try {
			File inputFile = new File(CONFIG_DIR, getLastUsedParametersFileName());

			if (!inputFile.exists()) {
				LOGGER.info("Last used parameters file not found: " + inputFile.getAbsolutePath());
				return null;
			}

			return loadParametersFromFile(inputFile);

		} catch (Exception e) {
			LOGGER.error("Failed to load last used parameters", e);
			return null;
		}
	}

	/**
	 * Load parameters from user preferences
	 */
	public Object loadPreferences() {
		try {
			String userHome = System.getProperty("user.home");
			File preferencesDir = new File(userHome, PREFERENCES_DIR);
			File inputFile = new File(preferencesDir, getPreferencesFileName());

			if (!inputFile.exists()) {
				LOGGER.info("Preferences file not found: " + inputFile.getAbsolutePath());
				return null;
			}

			return loadParametersFromFile(inputFile);

		} catch (Exception e) {
			LOGGER.error("Failed to load preferences", e);
			return null;
		}
	}

	/**
	 * Load analysis parameters from a specific result directory
	 */
	public Object loadAnalysisParameters(String resultDirectory) {
		try {
			if (resultDirectory == null || resultDirectory.trim().isEmpty()) {
				LOGGER.warn("Cannot load analysis parameters: result directory is null/empty");
				return null;
			}

			File inputFile = new File(resultDirectory, getAnalysisParameterFileName());

			if (!inputFile.exists()) {
				LOGGER.info("Analysis parameters file not found: " + inputFile.getAbsolutePath());
				return null;
			}

			return loadParametersFromFile(inputFile);

		} catch (Exception e) {
			LOGGER.error("Failed to load analysis parameters from " + resultDirectory, e);
			return null;
		}
	}

	/**
	 * Check if analysis parameters exist in the given directory
	 */
	public boolean hasAnalysisParameters(String resultDirectory) {
		File paramFile = new File(resultDirectory, getAnalysisParameterFileName());
		return paramFile.exists() && paramFile.length() > 0;
	}

	/**
	 * Check if last used parameters exist
	 */
	public boolean hasLastUsedParameters() {
		File paramFile = new File(CONFIG_DIR, getLastUsedParametersFileName());
		return paramFile.exists() && paramFile.length() > 0;
	}

	/**
	 * Check if user preferences exist
	 */
	public boolean hasPreferences() {
		String userHome = System.getProperty("user.home");
		File preferencesDir = new File(userHome, PREFERENCES_DIR);
		File paramFile = new File(preferencesDir, getPreferencesFileName());
		return paramFile.exists() && paramFile.length() > 0;
	}

	/**
	 * Save parameters to a specific file - to be implemented by subclasses
	 */
	protected abstract boolean saveParametersToFile(Object parameters, File outputFile, String description);

	/**
	 * Load parameters from a specific file - to be implemented by subclasses
	 * This is the instance method that each ParameterIO subclass must implement
	 * to handle the actual parsing of their specific parameter type.
	 */
	protected abstract Object loadParametersFromFile(File inputFile);

	/**
	 * Create default parameters - to be implemented by subclasses
	 */
	protected abstract Object createDefaultParameters();
	
	/**
	 * Get the workflow type from a parameter object
	 * This is used to identify which ParameterIO class to use
	 */
	public static String getWorkflowTypeFromParameters(Object parameters) {
		System.out.println("DEBUG: getWorkflowTypeFromParameters called with: " + (parameters != null ? parameters.getClass().getSimpleName() : "NULL"));
		
		if (parameters == null) {
			System.out.println("DEBUG: Parameters object is null, returning null");
			return null;
		}
		
		// Try to get workflow type from common parameter interfaces
		try {
			// If it's a MetaBaseParameter or extends it
			if (parameters instanceof uk.ac.quadram.metapilot.task.par.MetaBaseParameter) {
				uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParam = 
					(uk.ac.quadram.metapilot.task.par.MetaBaseParameter) parameters;
				
				// Get the original workflow type string from the parameter object
				String originalWorkflowType = baseParam.getWorkflowTypeString();
				System.out.println("DEBUG: Found original workflow type string: '" + originalWorkflowType + "'");
				return originalWorkflowType;
			}
			
			System.out.println("DEBUG: Not a MetaBaseParameter, trying reflection...");
			// Try reflection to get workflow type
			java.lang.reflect.Method getWorkflowTypeMethod = parameters.getClass().getMethod("getWorkflowType");
			Object workflowType = getWorkflowTypeMethod.invoke(parameters);
			String result = workflowType != null ? workflowType.toString() : null;
			System.out.println("DEBUG: Found workflow type via reflection: '" + result + "'");
			return result;
			
		} catch (Exception e) {
			System.out.println("DEBUG: Exception in getWorkflowTypeFromParameters: " + e.getClass().getSimpleName() + " - " + e.getMessage());
			LOGGER.debug("Could not determine workflow type from parameter object: " + e.getMessage());
			return null;
		}
	}
	
	/**
	 * Utility method for CLI to validate parameter file before loading
	 */
	public static boolean isValidParameterFilePath(String filePath) {
		try {
			File inputFile = new File(filePath);
			if (!inputFile.exists() || inputFile.length() == 0) {
				return false;
			}
			
			JSONObject json = readJsonFromFileStatic(inputFile);
			String workflowType = json.optString("workflow", null);
			
			return workflowType != null && !workflowType.trim().isEmpty() && 
				   PARAMETER_IO_REGISTRY.containsKey(workflowType);
			
		} catch (Exception e) {
			LOGGER.debug("Parameter file validation failed: " + e.getMessage());
			return false;
		}
	}
	
	/**
	 * Get list of supported workflow types for CLI help
	 */
	public static java.util.Set<String> getSupportedWorkflowTypes() {
		return new java.util.HashSet<>(PARAMETER_IO_REGISTRY.keySet());
	}

	/**
	 * Create a basic JSON header with common metadata
	 */
	protected JSONObject createJsonHeader(String description) {
		JSONObject json = new JSONObject();
		json.put("workflow", getWorkflowType());
		json.put("version", "1.0");
		json.put("generated_at", new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
		if (description != null) {
			json.put("description", description);
		}
		return json;
	}

	/**
	 * Write JSON content to file with proper error handling
	 */
	protected boolean writeJsonToFile(JSONObject json, File outputFile) {
		System.out.println("DEBUG: BaseParameterIO.writeJsonToFile() called");
		System.out.println("DEBUG: Output file: " + outputFile.getAbsolutePath());
		System.out.println("DEBUG: JSON object: " + (json != null ? "EXISTS" : "NULL"));

		if (json == null) {
			LOGGER.error("Cannot write null JSON object to file");
			return false;
		}

		try {
			System.out.println("DEBUG: Attempting to write JSON to file...");
			try (FileWriter writer = new FileWriter(outputFile)) {
				System.out.println("DEBUG: FileWriter created successfully");
				String jsonString = json.toString(2); // Pretty print with 2 spaces
				System.out.println("DEBUG: JSON string length: " + jsonString.length());
				writer.write(jsonString);
				System.out.println("DEBUG: JSON written to file successfully");
				LOGGER.info("Successfully saved parameters to: " + outputFile.getAbsolutePath());
				return true;
			}
		} catch (Exception e) {
			LOGGER.error("Failed to write parameters to file: " + e.getMessage(), e);
			System.err.println(
					"DEBUG: Exception in writeJsonToFile: " + e.getClass().getSimpleName() + " - " + e.getMessage());
			e.printStackTrace();
			return false;
		}
	}

	/**
	 * Read JSON content from file with proper error handling
	 */
	protected JSONObject readJsonFromFile(File inputFile) {
		try (FileReader reader = new FileReader(inputFile)) {
			StringBuilder content = new StringBuilder();
			char[] buffer = new char[1024];
			int bytesRead;
			while ((bytesRead = reader.read(buffer)) != -1) {
				content.append(buffer, 0, bytesRead);
			}
			return new JSONObject(content.toString());
		} catch (Exception e) {
			LOGGER.error("Failed to read parameters from file: " + e.getMessage(), e);
			throw new RuntimeException("Failed to read parameters from file", e);
		}
	}
}

package uk.ac.quadram.metapilot.task.dia.par;

import java.io.File;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

/**
 * Parameter I/O operations for DIA Marker workflow
 *
 * This class handles saving and loading of DIA Focus parameters to/from JSON
 * files for persistence and reproducibility.
 */
public class DiaFocusParameterIO extends BaseParameterIO {

	private static final Logger LOGGER = LogManager.getLogger(DiaFocusParameterIO.class);

	// Workflow-specific file names (new marker naming)
	private static final String PARAMETER_FILENAME = "dia_marker_analysis_parameters.json";
	private static final String LAST_PARAMETERS_FILENAME = "last_used_dia_marker_parameters.json";
	private static final String PREFERENCES_FILENAME = "dia_marker_preferences.json";

	// Legacy filenames kept for backward compatibility with existing runs
	private static final String LEGACY_PARAMETER_FILENAME = "dia_focus_analysis_parameters.json";
	private static final String LEGACY_LAST_PARAMETERS_FILENAME = "last_used_dia_focus_parameters.json";
	private static final String LEGACY_PREFERENCES_FILENAME = "dia_focus_preferences.json";

	@Override
	protected String getAnalysisParameterFileName() {
		return PARAMETER_FILENAME;
	}

	@Override
	protected String getLastUsedParametersFileName() {
		return LAST_PARAMETERS_FILENAME;
	}

	@Override
	protected String getPreferencesFileName() {
		return PREFERENCES_FILENAME;
	}

	// ---------- Backward-compatible loading helpers ----------

	@Override
	public Object loadAnalysisParameters(String resultDirectory) {
		try {
			if (resultDirectory == null || resultDirectory.trim().isEmpty()) {
				LOGGER.warn("Cannot load analysis parameters: result directory is null/empty");
				return null;
			}

			File newFile = new File(resultDirectory, PARAMETER_FILENAME);
			File legacyFile = new File(resultDirectory, LEGACY_PARAMETER_FILENAME);

			File chosen = newFile.exists() ? newFile : (legacyFile.exists() ? legacyFile : null);
			if (chosen == null) {
				LOGGER.info("Analysis parameters file not found (new or legacy) under {}", resultDirectory);
				return null;
			}
			return loadParametersFromFile(chosen);

		} catch (Exception e) {
			LOGGER.error("Failed to load analysis parameters from " + resultDirectory, e);
			return null;
		}
	}

	@Override
	public boolean hasAnalysisParameters(String resultDirectory) {
		File newFile = new File(resultDirectory, PARAMETER_FILENAME);
		File legacyFile = new File(resultDirectory, LEGACY_PARAMETER_FILENAME);
		return (newFile.exists() && newFile.length() > 0) || (legacyFile.exists() && legacyFile.length() > 0);
	}

	@Override
	public Object loadLastUsedParameters() {
		try {
			File newFile = new File(CONFIG_DIR, LAST_PARAMETERS_FILENAME);
			File legacyFile = new File(CONFIG_DIR, LEGACY_LAST_PARAMETERS_FILENAME);

			File chosen = newFile.exists() ? newFile : (legacyFile.exists() ? legacyFile : null);
			if (chosen == null) {
				LOGGER.info("Last used parameters file not found (new or legacy) in config directory");
				return null;
			}
			return loadParametersFromFile(chosen);

		} catch (Exception e) {
			LOGGER.error("Failed to load last used parameters", e);
			return null;
		}
	}

	@Override
	public boolean hasLastUsedParameters() {
		File newFile = new File(CONFIG_DIR, LAST_PARAMETERS_FILENAME);
		File legacyFile = new File(CONFIG_DIR, LEGACY_LAST_PARAMETERS_FILENAME);
		return (newFile.exists() && newFile.length() > 0) || (legacyFile.exists() && legacyFile.length() > 0);
	}

	@Override
	public Object loadPreferences() {
		try {
			String userHome = System.getProperty("user.home");
			File preferencesDir = new File(userHome, PREFERENCES_DIR);
			File newFile = new File(preferencesDir, PREFERENCES_FILENAME);
			File legacyFile = new File(preferencesDir, LEGACY_PREFERENCES_FILENAME);

			File chosen = newFile.exists() ? newFile : (legacyFile.exists() ? legacyFile : null);
			if (chosen == null) {
				LOGGER.info("Preferences file not found (new or legacy) under {}", preferencesDir.getAbsolutePath());
				return null;
			}
			return loadParametersFromFile(chosen);

		} catch (Exception e) {
			LOGGER.error("Failed to load preferences", e);
			return null;
		}
	}

	@Override
	public boolean hasPreferences() {
		String userHome = System.getProperty("user.home");
		File preferencesDir = new File(userHome, PREFERENCES_DIR);
		File newFile = new File(preferencesDir, PREFERENCES_FILENAME);
		File legacyFile = new File(preferencesDir, LEGACY_PREFERENCES_FILENAME);
		return (newFile.exists() && newFile.length() > 0) || (legacyFile.exists() && legacyFile.length() > 0);
	}

	@Override
	protected String getWorkflowType() {
		return uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType();
	}

	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof MetaParameterDiaFocus)) {
			LOGGER.error("Invalid parameter type: expected MetaParameterDiaFocus, got "
					+ (parameters != null ? parameters.getClass().getSimpleName() : "null"));
			return false;
		}

		MetaParameterDiaFocus params = (MetaParameterDiaFocus) parameters;

		try {
			// Ensure parent directories exist before writing the file
			File parentDir = outputFile.getParentFile();
			if (parentDir != null && !parentDir.exists()) {
				boolean created = parentDir.mkdirs();
				if (!created) {
					LOGGER.error("Failed to create parent directory: " + parentDir.getAbsolutePath());
					return false;
				}
			}

			JSONObject json = createJsonHeader(description);
			json.put("workflowType", getWorkflowType());
			json.put("timestamp", System.currentTimeMillis());

			// Add base parameters
			JSONObject baseParams = new JSONObject();
			baseParams.put("workflowType", params.getWorkflowTypeString());
			baseParams.put("result", params.getResult());
			baseParams.put("microDb", params.getMicroDb());
			baseParams.put("appendHostDb", params.isAppendHostDb());
			baseParams.put("hostDB", params.getHostDB());
			baseParams.put("ms2ScanMode", params.getMs2ScanMode());
			baseParams.put("threadCount", params.getThreadCount());
			baseParams.put("coreCount", params.getCoreCount());

			// Add metadata
			if (params.getMetadata() != null) {
				JSONObject metadataJson = new JSONObject();
				MetaData metadata = params.getMetadata();
				
				if (metadata.getRawFiles() != null) {
					JSONArray rawFilesArray = new JSONArray();
					for (String file : metadata.getRawFiles()) {
						rawFilesArray.put(file);
					}
					metadataJson.put("rawFiles", rawFilesArray);
				}
				
				if (metadata.getExpNames() != null) {
					JSONArray expNamesArray = new JSONArray();
					for (String name : metadata.getExpNames()) {
						expNamesArray.put(name);
					}
					metadataJson.put("expNames", expNamesArray);
				}
				
				if (metadata.getFractions() != null) {
					JSONArray fractionsArray = new JSONArray();
					for (int fraction : metadata.getFractions()) {
						fractionsArray.put(fraction);
					}
					metadataJson.put("fractions", fractionsArray);
				}
				
				if (metadata.getReplicates() != null) {
					JSONArray replicatesArray = new JSONArray();
					for (int replicate : metadata.getReplicates()) {
						replicatesArray.put(replicate);
					}
					metadataJson.put("replicates", replicatesArray);
				}
				
				baseParams.put("metadata", metadataJson);
			}

			json.put("baseParameters", baseParams);

			// Add MAG database information
			if (params.getMagDbItem() != null) {
				MagDbItem magDb = params.getMagDbItem();
				JSONObject magDbJson = new JSONObject();
				magDbJson.put("catalogueId", magDb.getCatalogueID());
				magDbJson.put("identifier", magDb.getIdentifier());
				
				// Add repository path if available
				if (magDb.getRepositoryFile() != null) {
					magDbJson.put("repositoryPath", magDb.getRepositoryFile().getAbsolutePath());
				}
				
				if (magDb.getUsedVersion() != null) {
					MagDbItem.MagVersion version = magDb.getUsedVersion();
					JSONObject versionJson = new JSONObject();
					versionJson.put("catalogueVersion", version.getCatalogueVersion());
					versionJson.put("speciesCount", version.getSpeciesCount());
					versionJson.put("localAvailable", version.isLocalAvailable());
					magDbJson.put("version", versionJson);
				}
				
				json.put("magDatabase", magDbJson);
			}

			// Add DIA-specific parameters
			json.put("librarySearch", params.isLibrarySearch());
			json.put("selfModelSearch", params.isSelfModelSearch());
			json.put("useModel", params.isUseModel());
			
			if (params.getLibrary() != null) {
				JSONArray libraryArray = new JSONArray();
				for (String lib : params.getLibrary()) {
					libraryArray.put(lib);
				}
				json.put("library", libraryArray);
			}
			
			if (params.getModels() != null) {
				JSONArray modelsArray = new JSONArray();
				for (String model : params.getModels()) {
					modelsArray.put(model);
				}
				json.put("models", modelsArray);
			}

			// Add genome selection strategy
			json.put("genomeSelectionStrategy", params.getGenomeSelectionStrategy());

			// Add modification settings (Focus-specific)
			json.put("mExcision", params.isMExcision());
			json.put("carbamiC", params.isCarbamiC());
			json.put("oxidationM", params.isOxidationM());
			json.put("acetyProN", params.isAcetyProN());
			json.put("phosphoSTY", params.isPhosphoSTY());
			json.put("kGG", params.isKGG());

			return writeJsonToFile(json, outputFile);

		} catch (Exception e) {
			LOGGER.error("Failed to save parameters to " + outputFile.getAbsolutePath(), e);
			return false;
		}
	}

	@Override
	protected Object loadParametersFromFile(File paramFile) {
		LOGGER.info("Loading parameters from {}", paramFile.getAbsolutePath());

		try {
			JSONObject json = readJsonFromFile(paramFile);
			return jsonToParameters(json);
		} catch (Exception e) {
			LOGGER.error("Failed to load parameters from " + paramFile.getAbsolutePath(), e);
			return createDefaultParameters();
		}
	}

	@Override
	protected Object createDefaultParameters() {
		// Create a default MetaParameterDiaFocus with minimal parameters
		uk.ac.quadram.metapilot.task.par.MetaBaseParameter metaBasePar = new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
				uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_MARKER.getAnalysisType(), new MetaData(), "", "",
				false, "", "DIA", 1, Runtime.getRuntime().availableProcessors(), true);

		return new MetaParameterDiaFocus(metaBasePar, null, true, false, new String[0], false, new String[0]);
	}

	private boolean validateParameters(Object parameters) {
		if (!(parameters instanceof MetaParameterDiaFocus)) {
			LOGGER.error("Invalid parameter type: expected MetaParameterDiaFocus, got {}",
					parameters != null ? parameters.getClass().getSimpleName() : "null");
			return false;
		}

		MetaParameterDiaFocus params = (MetaParameterDiaFocus) parameters;

		// Validate required fields
		if (params.getMagDbItem() == null) {
			LOGGER.error("MAG database item is null");
			return false;
		}

		if (params.getResult() == null || params.getResult().trim().isEmpty()) {
			LOGGER.error("Result directory is not set");
			return false;
		}

		return true;
	}

	/**
	 * Lenient validation for auto-save scenarios - only checks essential fields
	 */
	private boolean validateParametersForSave(Object parameters) {
		if (!(parameters instanceof MetaParameterDiaFocus)) {
			LOGGER.error("Invalid parameter type: expected MetaParameterDiaFocus, got {}",
					parameters != null ? parameters.getClass().getSimpleName() : "null");
			return false;
		}

		// Allow null MagDbItem and empty result directory for auto-save
		// These can be set later by the user

		return true;
	}

	private Object jsonToParameters(JSONObject json) {
		try {
			// Validate workflow type
			String workflowType = json.optString("workflowType", "");
			if (!getWorkflowType().equals(workflowType)) {
				LOGGER.warn("Workflow type mismatch: expected {}, got {}", getWorkflowType(), workflowType);
			}

			// Parse base parameters
			JSONObject baseParamsJson = json.getJSONObject("baseParameters");
			MetaData metadata = null;

			if (baseParamsJson.has("metadata")) {
				JSONObject metadataJson = baseParamsJson.getJSONObject("metadata");

				String[] rawFiles = null;
				String[] expNames = null;
				int[] fractions = null;
				int[] replicates = null;

				if (metadataJson.has("rawFiles")) {
					JSONArray rawFilesArray = metadataJson.getJSONArray("rawFiles");
					rawFiles = new String[rawFilesArray.length()];
					for (int i = 0; i < rawFilesArray.length(); i++) {
						rawFiles[i] = rawFilesArray.getString(i);
					}
				}

				if (metadataJson.has("expNames")) {
					JSONArray expNamesArray = metadataJson.getJSONArray("expNames");
					expNames = new String[expNamesArray.length()];
					for (int i = 0; i < expNamesArray.length(); i++) {
						expNames[i] = expNamesArray.getString(i);
					}
				}

				if (metadataJson.has("fractions")) {
					JSONArray fractionsArray = metadataJson.getJSONArray("fractions");
					fractions = new int[fractionsArray.length()];
					for (int i = 0; i < fractionsArray.length(); i++) {
						fractions[i] = fractionsArray.getInt(i);
					}
				}

				if (metadataJson.has("replicates")) {
					JSONArray replicatesArray = metadataJson.getJSONArray("replicates");
					replicates = new int[replicatesArray.length()];
					for (int i = 0; i < replicatesArray.length(); i++) {
						replicates[i] = replicatesArray.getInt(i);
					}
				}

				metadata = new MetaData(rawFiles, expNames, fractions, replicates);
			}

			String result = baseParamsJson.optString("result", "");
			String microDb = baseParamsJson.optString("microDb", "");
			boolean appendHostDb = baseParamsJson.optBoolean("appendHostDb", false);
			String hostDB = baseParamsJson.optString("hostDB", "");
			String ms2ScanMode = baseParamsJson.optString("ms2ScanMode", "DIA");
			int threadCount = baseParamsJson.optInt("threadCount", Runtime.getRuntime().availableProcessors());
			int coreCount = baseParamsJson.optInt("coreCount", 1);

			uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParams = 
				new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
					getWorkflowType(), metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, 
					coreCount, threadCount, true);

			// Parse MAG database
			MagDbItem magDbItem = null;
			if (json.has("magDatabase")) {
				JSONObject magDbJson = json.getJSONObject("magDatabase");
				String catalogueId = magDbJson.getString("catalogueId");
				String identifier = magDbJson.getString("identifier");
				String repositoryPath = magDbJson.optString("repositoryPath", "");
				
				// Create MagDbItem from JSON data
				if (!repositoryPath.isEmpty()) {
					magDbItem = new MagDbItem(catalogueId, identifier, new File(repositoryPath));
					
					// Set version if available
					if (magDbJson.has("version")) {
						JSONObject versionJson = magDbJson.getJSONObject("version");
						// Use addVersion method to create the version properly
						magDbItem.addVersion(
							versionJson.getString("catalogueVersion"), 
							versionJson.getInt("speciesCount"),
							new java.util.Date(), // Default date
							versionJson.getBoolean("localAvailable"));
					}
				}
			}

			// Parse DIA-specific parameters
			boolean librarySearch = json.optBoolean("librarySearch", true);
			boolean selfModelSearch = json.optBoolean("selfModelSearch", false);
			boolean useModel = json.optBoolean("useModel", false);
			
			String[] library = new String[0];
			if (json.has("library")) {
				JSONArray libraryArray = json.getJSONArray("library");
				library = new String[libraryArray.length()];
				for (int i = 0; i < libraryArray.length(); i++) {
					library[i] = libraryArray.getString(i);
				}
			}
			
			String[] models = new String[0];
			if (json.has("models")) {
				JSONArray modelsArray = json.getJSONArray("models");
				models = new String[modelsArray.length()];
				for (int i = 0; i < modelsArray.length(); i++) {
					models[i] = modelsArray.getString(i);
				}
			}

			// Create MetaParameterDiaFocus
			MetaParameterDiaFocus params = new MetaParameterDiaFocus(baseParams, magDbItem, librarySearch, 
					selfModelSearch, library, useModel, models);

			// Set genome selection strategy
			params.setGenomeSelectionStrategy(json.optString("genomeSelectionStrategy",
					MarkerGenomeSelector.KNEE_MODE_FOCUS));

			// Set modification settings (Focus-specific)
			params.setMExcision(json.optBoolean("mExcision", true));
			params.setCarbamiC(json.optBoolean("carbamiC", true));
			params.setOxidationM(json.optBoolean("oxidationM", false));
			params.setAcetyProN(json.optBoolean("acetyProN", false));
			params.setPhosphoSTY(json.optBoolean("phosphoSTY", false));
			params.setKGG(json.optBoolean("kGG", false));

			return params;

		} catch (Exception e) {
			LOGGER.error("Failed to parse parameters from JSON", e);
			return createDefaultParameters();
		}
	}
}

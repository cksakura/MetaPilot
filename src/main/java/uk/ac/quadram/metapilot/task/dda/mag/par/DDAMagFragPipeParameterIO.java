package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.io.File;

import org.json.JSONArray;
import org.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;

/**
 * Parameter I/O operations for DDA MAG FragPipe workflow
 *
 * This class handles saving and loading of DDA MAG FragPipe parameters
 * to/from JSON files for persistence and reproducibility.
 */
public class DDAMagFragPipeParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LoggerFactory.getLogger(DDAMagFragPipeParameterIO.class);

    // Workflow-specific file names
    private static final String PARAMETER_FILENAME = "mag_fragpipe_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_mag_fragpipe_parameters.json";
    private static final String PREFERENCES_FILENAME = "mag_fragpipe_preferences.json";

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

	@Override
	protected String getWorkflowType() {
		return uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType();
	}


	/**
	 * Load DDA MAG FragPipe parameters from JSON file
	 */
	public DDAMagFragPipeParameters loadParameters(File inputFile) {
		try {
			LOGGER.info("Loading parameters from: {}", inputFile.getAbsolutePath());
			JSONObject json = readJsonFromFile(inputFile);
			
			// Create MetaData object from experiments
			MetaData metadata = new MetaData();
			if (json.has("input_output") && json.getJSONObject("input_output").has("experiments")) {
				JSONArray experiments = json.getJSONObject("input_output").getJSONArray("experiments");
				String[] rawFiles = new String[experiments.length()];
				String[] expNames = new String[experiments.length()];
				int[] fractions = new int[experiments.length()];
				int[] replicates = new int[experiments.length()];

				for (int i = 0; i < experiments.length(); i++) {
					JSONObject exp = experiments.getJSONObject(i);
					rawFiles[i] = exp.optString("file_path", "");
					expNames[i] = exp.optString("experiment_name", "");
					fractions[i] = 1; // Default fraction
					replicates[i] = 1; // Default replicate
				}

				metadata = new MetaData(rawFiles, expNames, fractions, replicates);
				LOGGER.info("Loaded {} experiments from parameter file", experiments.length());
			}

			// Load Input/Output settings
			String result = "";
			String microDb = "";
			boolean appendHostDb = false;
			String hostDB = "";
			if (json.has("input_output")) {
				JSONObject ioSettings = json.getJSONObject("input_output");
				result = ioSettings.optString("output_directory", "");
				microDb = ioSettings.optString("mag_database", "");
				appendHostDb = ioSettings.optBoolean("append_host_db", false);
				hostDB = ioSettings.optString("host_database", "");
			}

			// Load additional parameters
			String ms2ScanMode = "HCD-FTMS";
			int threadCount = 4;
			if (json.has("additional_parameters")) {
				JSONObject additionalParams = json.getJSONObject("additional_parameters");
				ms2ScanMode = additionalParams.optString("ms2_scan_mode", "HCD-FTMS");
				threadCount = additionalParams.optInt("threads", 4);
			}

			// Create DDAMagFragPipeParameters
			DDAMagFragPipeParameters params = new DDAMagFragPipeParameters(
				uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType(),
				metadata,
				result,
				microDb,
				appendHostDb,
				hostDB,
				ms2ScanMode,
				threadCount
			);

			// Load enzyme settings
			if (json.has("enzyme_settings")) {
				JSONObject enzymeSettings = json.getJSONObject("enzyme_settings");
				params.setEnzyme(enzymeSettings.optString("enzyme", "Trypsin"));
				params.setMissCleavages(enzymeSettings.optInt("missed_cleavages", 2));
				params.setDigestMode(enzymeSettings.optString("digestion_mode", "Specific"));
			}

			// Load additional parameters
			if (json.has("additional_parameters")) {
				JSONObject additionalParams = json.getJSONObject("additional_parameters");
				params.setReduceDb(additionalParams.optBoolean("generate_sample_specific_db", true));
				params.setSpectraCluster(additionalParams.optBoolean("spectra_clustering", false));
			}

			// Set quantitative mode (default to label-free)
			params.setQuantitativeMode("Label-Free");

			// Load FragPipe workflow (simple flag)
			if (json.optBoolean("fragpipe_workflow_configured", false)) {
				// Workflow will be configured by UI/task
				params.setFragpipeWorkflow(null);
				LOGGER.info("Workflow configuration flag found, will be configured by UI");
			}

			LOGGER.info("Successfully loaded parameters from: {}", inputFile.getAbsolutePath());
			return params;
			
		} catch (Exception e) {
			LOGGER.error("Error loading DDA MAG FragPipe parameters from: {}", inputFile.getAbsolutePath(), e);
			throw new RuntimeException("Failed to load parameters", e);
		}
	}


	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof DDAMagFragPipeParameters)) {
			LOGGER.error("Invalid parameter type: expected DDAMagFragPipeParameters");
			return false;
		}
		
		DDAMagFragPipeParameters params = (DDAMagFragPipeParameters) parameters;
		
		try {
			LOGGER.info("=== Starting parameter save process ===");
			LOGGER.info("Output file: {}", outputFile.getAbsolutePath());
			LOGGER.info("Description: {}", description);
			
			// Ensure parent directories exist before writing the file
			File parentDir = outputFile.getParentFile();
			if (parentDir != null && !parentDir.exists()) {
				LOGGER.info("Creating parent directory: {}", parentDir.getAbsolutePath());
				boolean created = parentDir.mkdirs();
				if (!created) {
					LOGGER.error("Failed to create parent directory: {}", parentDir.getAbsolutePath());
					return false;
				}
			}

			JSONObject json = createJsonHeader(description);
			LOGGER.info("Created JSON header");

			// Input/Output settings
			JSONObject ioSettings = new JSONObject();
			ioSettings.put("output_directory", params.getResult());
			ioSettings.put("mag_database", params.getMicroDb());
			ioSettings.put("append_host_db", params.isAppendHostDb());
			ioSettings.put("host_database", params.getHostDB());
			LOGGER.info("Added I/O settings - output: {}, mag_db: {}", params.getResult(), params.getMicroDb());

			// Raw files and experiments
			JSONArray experiments = new JSONArray();
			MetaData metaData = params.getMetadata();
			if (metaData != null) {
				String[] rawFiles = metaData.getRawFiles();
				String[] expNames = metaData.getExpNames();
				if (rawFiles != null) {
					LOGGER.info("Adding {} experiments", rawFiles.length);
					for (int i = 0; i < rawFiles.length; i++) {
						JSONObject exp = new JSONObject();
						exp.put("file_path", rawFiles[i]);
						exp.put("experiment_name", (expNames != null && i < expNames.length) ? expNames[i] : "");
						experiments.put(exp);
					}
				}
			}
			ioSettings.put("experiments", experiments);
			json.put("input_output", ioSettings);

			// Enzyme settings
			JSONObject enzymeSettings = new JSONObject();
			enzymeSettings.put("enzyme", params.getEnzyme());
			enzymeSettings.put("missed_cleavages", params.getMissCleavages());
			enzymeSettings.put("digestion_mode", params.getDigestModeString());
			json.put("enzyme_settings", enzymeSettings);
			LOGGER.info("Added enzyme settings - enzyme: {}, missed_cleavages: {}", params.getEnzyme(), params.getMissCleavages());

			// Additional parameters
			JSONObject additionalParams = new JSONObject();
			additionalParams.put("threads", params.getThreadCount());
			additionalParams.put("ms2_scan_mode", params.getMs2ScanMode());
			additionalParams.put("generate_sample_specific_db", params.isReduceDb());
			additionalParams.put("spectra_clustering", params.isSpectraCluster());
			json.put("additional_parameters", additionalParams);
			LOGGER.info("Added additional parameters - threads: {}, ms2_scan_mode: {}", params.getThreadCount(), params.getMs2ScanMode());

			// FragPipe workflow (simple flag)
			if (params.getFragpipeWorkflow() != null) {
				json.put("fragpipe_workflow_configured", true);
				LOGGER.info("Workflow is configured, setting flag");
			} else {
				LOGGER.info("No workflow configured");
			}

			// Write formatted JSON
			LOGGER.info("Writing JSON to file...");
			boolean success = writeJsonToFile(json, outputFile);
			if (success) {
				LOGGER.info("=== Parameter save completed successfully ===");
			} else {
				LOGGER.error("=== Parameter save failed - writeJsonToFile returned false ===");
			}
			return success;
			
		} catch (Exception e) {
			LOGGER.error("=== Parameter save failed with exception ===", e);
			return false;
		}
	}

	@Override
	protected Object createDefaultParameters() {
		return new DDAMagFragPipeParameters();
	}

	@Override
	protected Object loadParametersFromFile(File inputFile) {
		return loadParameters(inputFile);
	}

}

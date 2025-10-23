package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;

/**
 * Parameter I/O operations for DDA MAG pFind workflow
 *
 * This class handles saving and loading of DDA MAG pFind parameters
 * to/from JSON files for persistence and reproducibility.
 */
public class DDAMagPfindParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LogManager.getLogger(DDAMagPfindParameterIO.class);

    // Workflow-specific file names
    private static final String PARAMETER_FILENAME = "mag_pfind_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_mag_pfind_parameters.json";
    private static final String PREFERENCES_FILENAME = "mag_pfind_preferences.json";

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
		return uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_MAG_PFIND.getAnalysisType();
	}

	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof DDAMagPfindParameters)) {
			LOGGER.error("Invalid parameter type: expected DDAMagPfindParameters");
			return false;
		}
		
		DDAMagPfindParameters params = (DDAMagPfindParameters) parameters;
		
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

			// Input/Output settings
			JSONObject ioSettings = new JSONObject();
			ioSettings.put("output_directory", params.getResult());
			ioSettings.put("mag_database", params.getMagDbCatalogueId());
			ioSettings.put("mag_database_version", params.getMagDbVersion());
			ioSettings.put("append_host_db", params.isAppendHostDb());
			ioSettings.put("host_database", params.getHostDB());

			// Raw files and experiments
			JSONArray experiments = new JSONArray();
			MetaData metaData = params.getMetadata();
			if (metaData != null) {
				String[] rawFiles = metaData.getRawFiles();
				String[] expNames = metaData.getExpNames();
				if (rawFiles != null) {
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

			// Additional parameters
			JSONObject additionalParams = new JSONObject();
			additionalParams.put("threads", params.getThreadCount());
			additionalParams.put("ms2_scan_mode", params.getMs2ScanMode());
			additionalParams.put("generate_sample_specific_db", params.isReduceDb());
			additionalParams.put("spectra_clustering", params.isSpectraCluster());
			json.put("additional_parameters", additionalParams);

			// Modifications
			JSONObject modifications = new JSONObject();
			JSONArray fixedMods = new JSONArray();
			if (params.getFixMods() != null) {
				for (String mod : params.getFixMods()) {
					fixedMods.put(mod);
				}
			}
			modifications.put("fixed_modifications", fixedMods);

			JSONArray varMods = new JSONArray();
			if (params.getVariMods() != null) {
				for (String mod : params.getVariMods()) {
					varMods.put(mod);
				}
			}
			modifications.put("variable_modifications", varMods);
			json.put("modifications", modifications);

			// Quantitative mode
			JSONObject quantMode = new JSONObject();
			quantMode.put("mode", params.getQuantitativeMode());

			if ("Isotopic".equals(params.getQuantitativeMode())) {
				quantMode.put("light_labels", new JSONArray(params.getLightLabels()));
				quantMode.put("medium_labels", new JSONArray(params.getMediumLabels()));
				quantMode.put("heavy_labels", new JSONArray(params.getHeavyLabels()));
			} else if ("Isobaric".equals(params.getQuantitativeMode())) {
				quantMode.put("isobaric_type", params.getIsobaricType());
				quantMode.put("combine_label", params.isCombineLabel());
			}

			json.put("quantitative_mode", quantMode);

			// Metadata
			JSONArray metadata = new JSONArray();
			if (params.getMetadataTableData() != null) {
				for (MetadataRow row : params.getMetadataTableData()) {
					JSONObject meta = new JSONObject();
					meta.put("sample_name", row.getSampleName());
					meta.put("label_type", row.getLabelType());
					meta.put("age", row.getAge());
					meta.put("gender", row.getGender());
					meta.put("condition", row.getCondition());
					meta.put("treatment", row.getTreatment());
					meta.put("time_point", row.getTimePoint());
					meta.put("custom_field_1", row.getCustomField1());
					meta.put("custom_field_2", row.getCustomField2());
					meta.put("custom_field_3", row.getCustomField3());
					metadata.put(meta);
				}
			}
			json.put("metadata", metadata);

			// Write formatted JSON
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

			// Create MetaData object
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
			}

			// Create DDAMagPfindParameters
			DDAMagPfindParameters params = new DDAMagPfindParameters(
				uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_MAG_PFIND.getAnalysisType(),
				metadata,
				"", // output_directory will be set below
				"", // mag_database will be set below
				false, // append_host_db will be set below
				"", // host_database will be set below
				"FTMS", // ms2_scan_mode will be set below
				Runtime.getRuntime().availableProcessors() // threads will be set below
			);

			// Load Input/Output settings
			if (json.has("input_output")) {
				JSONObject ioSettings = json.getJSONObject("input_output");
				params.setResult(ioSettings.optString("output_directory", ""));
				params.setMagDbCatalogueId(ioSettings.optString("mag_database", ""));
				params.setMagDbVersion(ioSettings.optString("mag_database_version", ""));
				params.setAppendHostDb(ioSettings.optBoolean("append_host_db", false));
				params.setHostDB(ioSettings.optString("host_database", ""));
			}

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
				params.setThreadCount(additionalParams.optInt("threads", Runtime.getRuntime().availableProcessors()));
				params.setMs2ScanMode(additionalParams.optString("ms2_scan_mode", "FTMS"));
				params.setReduceDb(additionalParams.optBoolean("generate_sample_specific_db", true));
				params.setSpectraCluster(additionalParams.optBoolean("spectra_clustering", true));
			}

			// Load modifications
			if (json.has("modifications")) {
				JSONObject modifications = json.getJSONObject("modifications");
				if (modifications.has("fixed_modifications")) {
					JSONArray fixedMods = modifications.getJSONArray("fixed_modifications");
					String[] fixMods = new String[fixedMods.length()];
					for (int i = 0; i < fixedMods.length(); i++) {
						fixMods[i] = fixedMods.getString(i);
					}
					params.setFixMods(fixMods);
				}
				if (modifications.has("variable_modifications")) {
					JSONArray variableMods = modifications.getJSONArray("variable_modifications");
					String[] variMods = new String[variableMods.length()];
					for (int i = 0; i < variableMods.length(); i++) {
						variMods[i] = variableMods.getString(i);
					}
					params.setVariMods(variMods);
				}
			}

			// Load quantitative mode settings
			if (json.has("quantitative_mode")) {
				JSONObject quantMode = json.getJSONObject("quantitative_mode");
				params.setQuantitativeMode(quantMode.optString("mode", "Label-free"));

				// Load isotopic labels
				if (quantMode.has("light_labels")) {
					JSONArray lightLabels = quantMode.getJSONArray("light_labels");
					List<String> lightLabelsList = new ArrayList<>();
					for (int i = 0; i < lightLabels.length(); i++) {
						lightLabelsList.add(lightLabels.getString(i));
					}
					params.setLightLabels(lightLabelsList);
				}
				if (quantMode.has("medium_labels")) {
					JSONArray mediumLabels = quantMode.getJSONArray("medium_labels");
					List<String> mediumLabelsList = new ArrayList<>();
					for (int i = 0; i < mediumLabels.length(); i++) {
						mediumLabelsList.add(mediumLabels.getString(i));
					}
					params.setMediumLabels(mediumLabelsList);
				}
				if (quantMode.has("heavy_labels")) {
					JSONArray heavyLabels = quantMode.getJSONArray("heavy_labels");
					List<String> heavyLabelsList = new ArrayList<>();
					for (int i = 0; i < heavyLabels.length(); i++) {
						heavyLabelsList.add(heavyLabels.getString(i));
					}
					params.setHeavyLabels(heavyLabelsList);
				}

				// Load isobaric settings
				if (quantMode.has("isobaric_type")) {
					params.setIsobaricType(quantMode.optString("isobaric_type", ""));
				}
				if (quantMode.has("combine_label")) {
					params.setCombineLabel(quantMode.optBoolean("combine_label", false));
				}
			}

			// Load metadata table data
			if (json.has("metadata")) {
				JSONArray metadataData = json.getJSONArray("metadata");
				List<MetadataRow> metadataTableData = new ArrayList<>();
				for (int i = 0; i < metadataData.length(); i++) {
					JSONObject row = metadataData.getJSONObject(i);
					MetadataRow metadataRow = new MetadataRow(row.optString("sample_name", ""));
					metadataRow.setLabelType(row.optString("label_type", ""));
					metadataRow.setAge(row.optString("age", ""));
					metadataRow.setGender(row.optString("gender", ""));
					metadataRow.setCondition(row.optString("condition", ""));
					metadataRow.setTreatment(row.optString("treatment", ""));
					metadataRow.setTimePoint(row.optString("time_point", ""));
					metadataRow.setCustomField1(row.optString("custom_field_1", ""));
					metadataRow.setCustomField2(row.optString("custom_field_2", ""));
					metadataRow.setCustomField3(row.optString("custom_field_3", ""));
					metadataTableData.add(metadataRow);
				}
				params.setMetadataTableData(metadataTableData);
			}

			LOGGER.info("Successfully loaded parameters from {}", paramFile.getAbsolutePath());
			return params;

		} catch (Exception e) {
			LOGGER.error("Failed to load parameters from " + paramFile.getAbsolutePath(), e);
			return createDefaultParameters();
		}
	}

	@Override
	protected Object createDefaultParameters() {
		return new DDAMagPfindParameters();
	}
}

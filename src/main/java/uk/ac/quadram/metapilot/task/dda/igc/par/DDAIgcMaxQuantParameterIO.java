package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.prefs.Preferences;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.IsobaricLabelRow;

/**
 * Parameter IO Handler for DDA IGC MaxQuant Workflow
 * 
 * Handles saving and loading of analysis parameters for reproducibility and
 * user convenience. Supports both result-specific and application-level
 * parameter persistence.
 */
public class DDAIgcMaxQuantParameterIO {

	private static final Logger LOGGER = LogManager.getLogger(DDAIgcMaxQuantParameterIO.class);
	private static final String PREFERENCES_NODE = "uk/ac/quadram/metapilot/dda/igc/maxquant";
	private static final String PARAMETER_FILENAME = "analysis_parameters.json";
	private static final String LAST_PARAMETERS_FILENAME = "config/last_dda_igc_maxquant_parameters.json";

	/**
	 * Save parameters to analysis result directory for reproducibility
	 */
	public static boolean saveAnalysisParameters(String resultDirectory, DDAIgcMaxQuantParameters parameters) {
		File paramFile = new File(resultDirectory, PARAMETER_FILENAME);
		return saveParametersToFile(paramFile, parameters, "Analysis-specific parameters");
	}

	/**
	 * Save parameters as last used settings for user convenience
	 */
	public static boolean saveLastUsedParameters(DDAIgcMaxQuantParameters parameters) {
		File paramFile = new File(LAST_PARAMETERS_FILENAME);
		// Create config directory if it doesn't exist
		paramFile.getParentFile().mkdirs();
		return saveParametersToFile(paramFile, parameters, "Last used parameters");
	}

	/**
	 * Load last used parameters for UI initialization
	 */
	public static DDAIgcMaxQuantParameters loadLastUsedParameters() {
		File paramFile = new File(LAST_PARAMETERS_FILENAME);
		if (!paramFile.exists()) {
			LOGGER.info("No last used parameters found, using defaults");
			return createDefaultParameters();
		}

		try {
			return loadParametersFromFile(paramFile);
		} catch (Exception e) {
			LOGGER.warn("Failed to load last used parameters, using defaults: " + e.getMessage());
			return createDefaultParameters();
		}
	}

	/**
	 * Load parameters from analysis result directory
	 */
	public static DDAIgcMaxQuantParameters loadAnalysisParameters(String resultDirectory) {
		File paramFile = new File(resultDirectory, PARAMETER_FILENAME);
		if (!paramFile.exists()) {
			throw new IllegalArgumentException("Parameter file not found: " + paramFile.getAbsolutePath());
		}

		return loadParametersFromFile(paramFile);
	}

	private static boolean saveParametersToFile(File paramFile, DDAIgcMaxQuantParameters parameters,
			String description) {
		try {
			// Ensure parent directories exist before writing the file
			File parentDir = paramFile.getParentFile();
			if (parentDir != null && !parentDir.exists()) {
				boolean created = parentDir.mkdirs();
				if (!created) {
					LOGGER.error("Failed to create parent directory: " + parentDir.getAbsolutePath());
					return false;
				}
			}

			try (FileWriter writer = new FileWriter(paramFile)) {
				JSONObject json = new JSONObject();

				// Header information
				json.put("workflow", "DDA_IGC_MAXQUANT");
				json.put("version", "1.0");
				json.put("generated_at", new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
				json.put("description", description);

				// Input/Output settings
				JSONObject ioSettings = new JSONObject();
				ioSettings.put("output_directory", parameters.getResult());
				ioSettings.put("microbiome_database", parameters.getMicroDb());
				ioSettings.put("append_host_db", parameters.isAppendHostDb());
				ioSettings.put("host_database", parameters.getHostDB());

				// Raw files and experiments
				JSONArray experiments = new JSONArray();
				MetaData metaData = parameters.getMetadata();
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
				enzymeSettings.put("enzyme", parameters.getEnzyme());
				enzymeSettings.put("missed_cleavages", parameters.getMissCleavages());
				enzymeSettings.put("digestion_mode", parameters.getDigestModeString());
				json.put("enzyme_settings", enzymeSettings);

				// Additional parameters
				JSONObject additionalParams = new JSONObject();
				additionalParams.put("threads", parameters.getThreadCount());
				additionalParams.put("ms2_scan_mode", parameters.getMs2ScanMode());
				additionalParams.put("generate_sample_specific_db", parameters.isReduceDb());
				additionalParams.put("spectra_clustering", parameters.isSpectraCluster());
				json.put("additional_parameters", additionalParams);

				// Modifications
				JSONObject modifications = new JSONObject();
				modifications.put("fixed_modifications", new JSONArray(parameters.getFixedModifications()));
				modifications.put("variable_modifications", new JSONArray(parameters.getVariableModifications()));
				json.put("modifications", modifications);

				// Quantitative mode
				JSONObject quantMode = new JSONObject();
				quantMode.put("mode", parameters.getQuantitativeMode());

				if ("Isotopic".equals(parameters.getQuantitativeMode())) {
					quantMode.put("light_labels", new JSONArray(parameters.getLightLabels()));
					quantMode.put("medium_labels", new JSONArray(parameters.getMediumLabels()));
					quantMode.put("heavy_labels", new JSONArray(parameters.getHeavyLabels()));
				} else if ("Isobaric".equals(parameters.getQuantitativeMode())) {
					quantMode.put("isobaric_type", parameters.getIsobaricType());

					JSONArray isobaricLabels = new JSONArray();
					for (IsobaricLabelRow row : parameters.getIsobaricTableData()) {
						JSONObject label = new JSONObject();
						label.put("label_name", row.getLabelName());
						label.put("minus_2", row.getMinus2());
						label.put("minus_1", row.getMinus1());
						label.put("plus_1", row.getPlus1());
						label.put("plus_2", row.getPlus2());
						isobaricLabels.put(label);
					}
					quantMode.put("isobaric_labels", isobaricLabels);
				}

				json.put("quantitative_mode", quantMode);

				// Taxonomic analysis
				JSONObject taxonomic = new JSONObject();
				taxonomic.put("eukaryota", parameters.isEukaryota());
				taxonomic.put("bacteria", parameters.isBacteria());
				taxonomic.put("archaea", parameters.isArchaea());
				taxonomic.put("viruses", parameters.isViruses());
				taxonomic.put("viroids", parameters.isViroids());
				taxonomic.put("unclassified", parameters.isUnclassified());
				taxonomic.put("other", parameters.isOther());
				taxonomic.put("ignore_rank", parameters.getIgnoreRank());
				taxonomic.put("ignore_environmental_samples", parameters.isIgnoreEnvironmentalSamples());
				json.put("taxonomic_analysis", taxonomic);

				// Metadata
				JSONArray metadata = new JSONArray();
				for (MetadataRow row : parameters.getMetadataTableData()) {
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
				json.put("metadata", metadata);

				// Write formatted JSON
				writer.write(json.toString(2));

				LOGGER.info("Successfully saved {} to {}", description, paramFile.getAbsolutePath());
				return true;
			}
		} catch (IOException e) {
			LOGGER.error("Failed to save parameters to " + paramFile.getAbsolutePath(), e);
			return false;
		}
	}

	private static DDAIgcMaxQuantParameters loadParametersFromFile(File paramFile) {
		LOGGER.info("Loading parameters from {}", paramFile.getAbsolutePath());

		try {
			// Read the JSON file
			String jsonContent = new String(java.nio.file.Files.readAllBytes(paramFile.toPath()));
			JSONObject json = new JSONObject(jsonContent);

			DDAIgcMaxQuantParameters parameters = new DDAIgcMaxQuantParameters();

			// Load Input/Output settings
			if (json.has("input_output")) {
				JSONObject ioSettings = json.getJSONObject("input_output");
				parameters.setResult(ioSettings.optString("output_directory", ""));
				parameters.setMicroDb(ioSettings.optString("microbiome_database", ""));
				parameters.setAppendHostDb(ioSettings.optBoolean("append_host_db", false));
				parameters.setHostDB(ioSettings.optString("host_database", ""));

				// Load experiment data
				if (ioSettings.has("experiments")) {
					JSONArray experiments = ioSettings.getJSONArray("experiments");
					String[] rawFiles = new String[experiments.length()];
					String[] expNames = new String[experiments.length()];
					for (int i = 0; i < experiments.length(); i++) {
						JSONObject exp = experiments.getJSONObject(i);
						rawFiles[i] = exp.optString("file_path", "");
						expNames[i] = exp.optString("experiment_name", "");
					}
					MetaData metaData = new MetaData(rawFiles, expNames);
					parameters.setMetadata(metaData);
				}
			}

			// Load enzyme settings
			if (json.has("enzyme_settings")) {
				JSONObject enzymeSettings = json.getJSONObject("enzyme_settings");
				parameters.setEnzyme(enzymeSettings.optString("enzyme", "Trypsin"));
				parameters.setMissCleavages(enzymeSettings.optInt("missed_cleavages", 2));
				parameters.setDigestMode(enzymeSettings.optString("digestion_mode", "Specific"));
			}

			// Load additional parameters
			if (json.has("additional_parameters")) {
				JSONObject additionalParams = json.getJSONObject("additional_parameters");
				parameters
						.setThreadCount(additionalParams.optInt("threads", Runtime.getRuntime().availableProcessors()));
				parameters.setMs2ScanMode(additionalParams.optString("ms2_scan_mode", "FTMS"));
				parameters.setReduceDb(additionalParams.optBoolean("generate_sample_specific_db", true));
				parameters.setSpectraCluster(additionalParams.optBoolean("spectra_clustering", true));
			}

			// Load modifications
			if (json.has("modifications")) {
				JSONObject modifications = json.getJSONObject("modifications");
				if (modifications.has("fixed_modifications")) {
					JSONArray fixedMods = modifications.getJSONArray("fixed_modifications");
					List<String> fixedModifications = new ArrayList<>();
					for (int i = 0; i < fixedMods.length(); i++) {
						fixedModifications.add(fixedMods.getString(i));
					}
					parameters.setFixedModifications(fixedModifications);
				}
				if (modifications.has("variable_modifications")) {
					JSONArray variableMods = modifications.getJSONArray("variable_modifications");
					List<String> variableModifications = new ArrayList<>();
					for (int i = 0; i < variableMods.length(); i++) {
						variableModifications.add(variableMods.getString(i));
					}
					parameters.setVariableModifications(variableModifications);
				}
			}

			// Load quantitative mode settings
			if (json.has("quantitative_mode")) {
				JSONObject quantMode = json.getJSONObject("quantitative_mode");
				parameters.setQuantitativeMode(quantMode.optString("mode", "Label-free"));

				// Load isotopic labels
				if (quantMode.has("light_labels")) {
					JSONArray lightLabels = quantMode.getJSONArray("light_labels");
					List<String> lightLabelsList = new ArrayList<>();
					for (int i = 0; i < lightLabels.length(); i++) {
						lightLabelsList.add(lightLabels.getString(i));
					}
					parameters.setLightLabels(lightLabelsList);
				}
				if (quantMode.has("medium_labels")) {
					JSONArray mediumLabels = quantMode.getJSONArray("medium_labels");
					List<String> mediumLabelsList = new ArrayList<>();
					for (int i = 0; i < mediumLabels.length(); i++) {
						mediumLabelsList.add(mediumLabels.getString(i));
					}
					parameters.setMediumLabels(mediumLabelsList);
				}
				if (quantMode.has("heavy_labels")) {
					JSONArray heavyLabels = quantMode.getJSONArray("heavy_labels");
					List<String> heavyLabelsList = new ArrayList<>();
					for (int i = 0; i < heavyLabels.length(); i++) {
						heavyLabelsList.add(heavyLabels.getString(i));
					}
					parameters.setHeavyLabels(heavyLabelsList);
				}

				// Load isobaric settings
				if (quantMode.has("isobaric_type")) {
					parameters.setIsobaricType(quantMode.optString("isobaric_type", ""));
				}
				if (quantMode.has("isobaric_table_data")) {
					JSONArray isobaricData = quantMode.getJSONArray("isobaric_table_data");
					List<MetaPilotDdaIgcMaxQuantPanel.IsobaricLabelRow> isobaricTableData = new ArrayList<>();
					for (int i = 0; i < isobaricData.length(); i++) {
						JSONObject row = isobaricData.getJSONObject(i);
						MetaPilotDdaIgcMaxQuantPanel.IsobaricLabelRow labelRow = new MetaPilotDdaIgcMaxQuantPanel.IsobaricLabelRow(
								row.optString("label_name", ""));
						labelRow.setMinus2(row.optString("minus2", "0.0"));
						labelRow.setMinus1(row.optString("minus1", "0.0"));
						labelRow.setPlus1(row.optString("plus1", "0.0"));
						labelRow.setPlus2(row.optString("plus2", "0.0"));
						isobaricTableData.add(labelRow);
					}
					parameters.setIsobaricTableData(isobaricTableData);
				}
			}

			// Load taxonomic analysis settings
			if (json.has("taxonomic_analysis")) {
				JSONObject taxonomic = json.getJSONObject("taxonomic_analysis");
				parameters.setEukaryota(taxonomic.optBoolean("eukaryota", false));
				parameters.setBacteria(taxonomic.optBoolean("bacteria", false));
				parameters.setArchaea(taxonomic.optBoolean("archaea", false));
				parameters.setViruses(taxonomic.optBoolean("viruses", false));
				parameters.setViroids(taxonomic.optBoolean("viroids", false));
				parameters.setUnclassified(taxonomic.optBoolean("unclassified", false));
				parameters.setOther(taxonomic.optBoolean("other", false));
				parameters.setIgnoreRank(taxonomic.optString("ignore_rank", "Species"));
				parameters.setIgnoreEnvironmentalSamples(taxonomic.optBoolean("ignore_environmental_samples", false));
			}

			// Load metadata table data
			if (json.has("metadata_table_data")) {
				JSONArray metadataData = json.getJSONArray("metadata_table_data");
				List<MetaPilotDdaIgcMaxQuantPanel.MetadataRow> metadataTableData = new ArrayList<>();
				for (int i = 0; i < metadataData.length(); i++) {
					JSONObject row = metadataData.getJSONObject(i);
					MetaPilotDdaIgcMaxQuantPanel.MetadataRow metadataRow = new MetaPilotDdaIgcMaxQuantPanel.MetadataRow(
							row.optString("sample_name", ""));
					metadataRow.setLabelType(row.optString("label_type", ""));
					metadataRow.setAge(row.optString("age", ""));
					metadataRow.setGender(row.optString("gender", ""));
					metadataRow.setCondition(row.optString("condition", ""));
					metadataRow.setTreatment(row.optString("treatment", ""));
					metadataRow.setTimePoint(row.optString("time_point", ""));
					metadataRow.setCustomField1(row.optString("custom_field1", ""));
					metadataRow.setCustomField2(row.optString("custom_field2", ""));
					metadataRow.setCustomField3(row.optString("custom_field3", ""));
					metadataTableData.add(metadataRow);
				}
				parameters.setMetadataTableData(metadataTableData);
			}

			LOGGER.info("Successfully loaded parameters from {}", paramFile.getAbsolutePath());
			return parameters;

		} catch (Exception e) {
			LOGGER.error("Failed to load parameters from " + paramFile.getAbsolutePath(), e);
			return createDefaultParameters();
		}
	}

	private static DDAIgcMaxQuantParameters createDefaultParameters() {
		return new DDAIgcMaxQuantParameters();
	}

	/**
	 * Save basic preferences (non-file specific settings) to Java preferences
	 */
	public static void savePreferences(DDAIgcMaxQuantParameters parameters) {
		try {
			Preferences prefs = Preferences.userRoot().node(PREFERENCES_NODE);

			// Save commonly reused settings with null checks
			String enzyme = parameters.getEnzyme();
			prefs.put("default_enzyme", enzyme != null ? enzyme : "Trypsin");
			
			prefs.putInt("default_missed_cleavages", parameters.getMissCleavages());
			
			String digestMode = parameters.getDigestModeString();
			prefs.put("default_digestion_mode", digestMode != null ? digestMode : "Specific");
			
			prefs.putInt("default_threads", parameters.getThreadCount());
			
			String ms2ScanMode = parameters.getMs2ScanMode();
			prefs.put("default_ms2_scan_mode", ms2ScanMode != null ? ms2ScanMode : "FTMS");
			
			prefs.putBoolean("default_generate_sample_specific_db", parameters.isReduceDb());
			prefs.putBoolean("default_spectra_clustering", parameters.isSpectraCluster());

			LOGGER.info("Saved user preferences for DDA IGC MaxQuant workflow");
		} catch (Exception e) {
			LOGGER.warn("Failed to save user preferences: " + e.getMessage());
			// Don't throw the exception - preferences saving is not critical
		}
	}

	/**
	 * Load basic preferences from Java preferences
	 */
	public static DDAIgcMaxQuantParameters loadPreferences() {
		Preferences prefs = Preferences.userRoot().node(PREFERENCES_NODE);
		DDAIgcMaxQuantParameters parameters = new DDAIgcMaxQuantParameters();

		// Load commonly reused settings with defaults
		parameters.setEnzyme(prefs.get("default_enzyme", "Trypsin"));
		parameters.setMissCleavages(prefs.getInt("default_missed_cleavages", 2));
		parameters.setDigestMode(prefs.get("default_digestion_mode", "Specific"));
		parameters.setThreadCount(prefs.getInt("default_threads", Runtime.getRuntime().availableProcessors()));
		parameters.setMs2ScanMode(prefs.get("default_ms2_scan_mode", "FTMS"));
		parameters.setReduceDb(prefs.getBoolean("default_generate_sample_specific_db", true));
		parameters.setSpectraCluster(prefs.getBoolean("default_spectra_clustering", true));

		LOGGER.info("Loaded user preferences for DDA IGC MaxQuant workflow");
		return parameters;
	}
}

package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.SampleMetadata;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;

/**
 * Parameter IO Handler for DDA IGC MaxQuant Workflow
 * 
 * Handles saving and loading of analysis parameters for reproducibility and
 * user convenience. Supports both result-specific and application-level
 * parameter persistence.
 */
public class DDAIgcMaxQuantParameterIO extends BaseParameterIO {

	private static final Logger LOGGER = LogManager.getLogger(DDAIgcMaxQuantParameterIO.class);
	
	// Workflow-specific file names
	private static final String PARAMETER_FILENAME = "analysis_parameters.json";
	private static final String LAST_PARAMETERS_FILENAME = "last_dda_igc_maxquant_parameters.json";
	private static final String PREFERENCES_FILENAME = "maxquant_preferences.json";

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
		return uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_IGC_MAXQUANT.getAnalysisType();
	}

	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof DDAIgcMaxQuantParameters)) {
			LOGGER.error("Invalid parameter type: expected DDAIgcMaxQuantParameters");
			return false;
		}
		
		DDAIgcMaxQuantParameters params = (DDAIgcMaxQuantParameters) parameters;
		
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
			ioSettings.put("microbiome_database", params.getMicroDb());
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
			modifications.put("fixed_modifications", new JSONArray(params.getFixedModifications()));
			modifications.put("variable_modifications", new JSONArray(params.getVariableModifications()));
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

				JSONArray isobaricLabels = new JSONArray();
				for (String[] row : params.getIsobaricTableData()) {
					JSONObject label = new JSONObject();
					label.put("label_name", row[0]);
					label.put("minus_2", row[1]);
					label.put("minus_1", row[2]);
					label.put("plus_1", row[3]);
					label.put("plus_2", row[4]);
					isobaricLabels.put(label);
				}
				quantMode.put("isobaric_labels", isobaricLabels);
			}

			json.put("quantitative_mode", quantMode);

			// Taxonomic analysis
			JSONObject taxonomic = new JSONObject();
			taxonomic.put("eukaryota", params.isEukaryota());
			taxonomic.put("bacteria", params.isBacteria());
			taxonomic.put("archaea", params.isArchaea());
			taxonomic.put("viruses", params.isViruses());
			taxonomic.put("viroids", params.isViroids());
			taxonomic.put("unclassified", params.isUnclassified());
			taxonomic.put("other", params.isOther());
			taxonomic.put("ignore_rank", params.getIgnoreRank());
			taxonomic.put("ignore_environmental_samples", params.isIgnoreEnvironmentalSamples());
			json.put("taxonomic_analysis", taxonomic);

			// Metadata - Handle String[] arrays
			JSONArray metadata = new JSONArray();
			for (String[] row : params.getMetadataTableData()) {
				JSONObject meta = new JSONObject();
				meta.put("sample_name", row[0]);
				meta.put("label_type", row[1]);
				meta.put("age", row[2]);
				meta.put("gender", row[3]);
				meta.put("condition", row[4]);
				meta.put("treatment", row[5]);
				meta.put("time_point", row[6]);
				meta.put("custom_field_1", row[7]);
				meta.put("custom_field_2", row[8]);
				meta.put("custom_field_3", row[9]);
				metadata.put(meta);
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
					List<String[]> isobaricTableData = new ArrayList<>();
					for (int i = 0; i < isobaricData.length(); i++) {
						JSONObject row = isobaricData.getJSONObject(i);
						String[] labelRow = new String[5];
						labelRow[0] = row.optString("label_name", "");
						labelRow[1] = row.optString("minus2", "0.0");
						labelRow[2] = row.optString("minus1", "0.0");
						labelRow[3] = row.optString("plus1", "0.0");
						labelRow[4] = row.optString("plus2", "0.0");
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

			// Load metadata table data - Handle String[] arrays
			if (json.has("metadata_table_data")) {
				JSONArray metadataData = json.getJSONArray("metadata_table_data");
				List<String[]> metadataTableData = new ArrayList<>();
				for (int i = 0; i < metadataData.length(); i++) {
					JSONObject row = metadataData.getJSONObject(i);
					String[] metadataRow = new String[10];
					metadataRow[0] = row.optString("sample_name", "");
					metadataRow[1] = row.optString("label_type", "");
					metadataRow[2] = row.optString("age", "");
					metadataRow[3] = row.optString("gender", "");
					metadataRow[4] = row.optString("condition", "");
					metadataRow[5] = row.optString("treatment", "");
					metadataRow[6] = row.optString("time_point", "");
					metadataRow[7] = row.optString("custom_field_1", "");
					metadataRow[8] = row.optString("custom_field_2", "");
					metadataRow[9] = row.optString("custom_field_3", "");
					metadataTableData.add(metadataRow);
				}
				parameters.setMetadataTableData(metadataTableData);
				SampleMetadata.apply(parameters.getMetadata(), metadataTableData);
			}

			LOGGER.info("Successfully loaded parameters from {}", paramFile.getAbsolutePath());
			return parameters;

		} catch (Exception e) {
			LOGGER.error("Failed to load parameters from " + paramFile.getAbsolutePath(), e);
			return createDefaultParameters();
		}
	}

	@Override
	protected Object createDefaultParameters() {
		return new DDAIgcMaxQuantParameters();
	}
}

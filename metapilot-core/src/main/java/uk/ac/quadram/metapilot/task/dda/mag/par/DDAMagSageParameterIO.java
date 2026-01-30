package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;

/**
 * Parameter I/O operations for DDA MAG Sage workflow
 *
 * This class handles saving and loading of DDA MAG Sage parameters
 * to/from JSON files for persistence and reproducibility.
 * Only includes important parameters that users need to modify.
 */
public class DDAMagSageParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LogManager.getLogger(DDAMagSageParameterIO.class);

    // Workflow-specific file names
    private static final String PARAMETER_FILENAME = "mag_sage_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_mag_sage_parameters.json";
    private static final String PREFERENCES_FILENAME = "mag_sage_preferences.json";

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
		return "DDA_MAG_SAGE";
	}

	/**
	 * Load DDA MAG Sage parameters from JSON file
	 */
	public DDAMagSageParameters loadParameters(File inputFile) {
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
			int threadCount = 4;
			if (json.has("additional_parameters")) {
				JSONObject additionalParams = json.getJSONObject("additional_parameters");
				threadCount = additionalParams.optInt("threads", 4);
			}

			// Create DDAMagSageParameters
			DDAMagSageParameters params = new DDAMagSageParameters(
				"DDA_MAG_SAGE",
				metadata,
				result,
				microDb,
				appendHostDb,
				hostDB,
				"HCD-FTMS", // Default MS2 scan mode for Sage
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
			params.setRepresentativeLevel(additionalParams.optString("representative_level", "GENUS"));
			}

			// Load Sage-specific important parameters
			if (json.has("sage_parameters")) {
				JSONObject sageParams = json.getJSONObject("sage_parameters");
				
				// Load search type first to determine how to load tolerances
				String searchType = sageParams.optString("search_type", "Closed Search");
				
				if (searchType.equals("Open Search")) {
					// For open search, load the saved min/max values directly
					params.setPrecursorToleranceMin(sageParams.optDouble("precursor_tolerance_min", -500.0));
					params.setPrecursorToleranceMax(sageParams.optDouble("precursor_tolerance_max", 150.0));
					params.setPrecursorToleranceUnit("Da");
				} else {
					// For closed search, load symmetric tolerance or fallback to defaults
					double tolerance = sageParams.optDouble("precursor_tolerance", 20.0);
					params.setPrecursorToleranceMin(-tolerance);
					params.setPrecursorToleranceMax(tolerance);
					params.setPrecursorToleranceUnit("ppm");
				}
				
				// Keep legacy single tolerance for backward compatibility
				params.setFragmentTolerance(sageParams.optDouble("fragment_tolerance", 20.0));
				params.setFragmentToleranceUnit(sageParams.optString("fragment_tolerance_unit", "ppm"));
				params.setMinPeptideLength(sageParams.optInt("min_peptide_length", 7));
				params.setMaxPeptideLength(sageParams.optInt("max_peptide_length", 50));
				params.setMaxVariableModifications(sageParams.optInt("max_variable_modifications", 2));
				params.setDecoysEnabled(sageParams.optBoolean("decoys_enabled", true));
				params.setGenerateDecoys(sageParams.optBoolean("generate_decoys", true));
				params.setDecoysTag(sageParams.optString("decoys_tag", "rev_"));
				params.setMaxQ(sageParams.optDouble("max_q", 0.01));
			}

			// Set quantitative mode (default to label-free)
			params.setQuantitativeMode("Label-Free");

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

			// Load metadata table data as String arrays
			if (json.has("metadata")) {
				JSONArray metadataData = json.getJSONArray("metadata");
				List<String[]> metadataTableData = new ArrayList<>();
				for (int i = 0; i < metadataData.length(); i++) {
					JSONObject row = metadataData.getJSONObject(i);
					// MetadataData format: [sampleName, labelType, age, gender, condition, treatment, timePoint, customField1, customField2, customField3]
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
				params.setMetadataTableData(metadataTableData);
			}

			LOGGER.info("Successfully loaded parameters from: {}", inputFile.getAbsolutePath());
			return params;
			
		} catch (Exception e) {
			LOGGER.error("Error loading DDA MAG Sage parameters from: {}", inputFile.getAbsolutePath(), e);
			throw new RuntimeException("Failed to load parameters", e);
		}
	}

	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof DDAMagSageParameters)) {
			LOGGER.error("Invalid parameter type: expected DDAMagSageParameters");
			return false;
		}
		
		DDAMagSageParameters params = (DDAMagSageParameters) parameters;
		
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
			additionalParams.put("generate_sample_specific_db", params.isReduceDb());
			additionalParams.put("spectra_clustering", params.isSpectraCluster());
		additionalParams.put("representative_level", params.getRepresentativeLevel());
			json.put("additional_parameters", additionalParams);
			LOGGER.info("Added additional parameters - threads: {}", params.getThreadCount());

			// Sage-specific important parameters (only the ones users need to modify)
			JSONObject sageParams = new JSONObject();
			// Save separate min/max precursor tolerance values
			sageParams.put("precursor_tolerance_min", params.getPrecursorToleranceMin());
			sageParams.put("precursor_tolerance_max", params.getPrecursorToleranceMax());
			// Save search type to distinguish between Open/Closed Search
			String searchType = params.getPrecursorToleranceUnit().equals("Da") ? "Open Search" : "Closed Search";
			sageParams.put("search_type", searchType);
			// Don't save the old precursor_tolerance field - only keep for loading old files
			sageParams.put("fragment_tolerance", params.getFragmentTolerance());
			sageParams.put("precursor_tolerance_unit", params.getPrecursorToleranceUnit());
			sageParams.put("fragment_tolerance_unit", params.getFragmentToleranceUnit());
			sageParams.put("min_peptide_length", params.getMinPeptideLength());
			sageParams.put("max_peptide_length", params.getMaxPeptideLength());
			sageParams.put("max_variable_modifications", params.getMaxVariableModifications());
			sageParams.put("decoys_enabled", params.isDecoysEnabled());
			sageParams.put("generate_decoys", params.isGenerateDecoys());
			sageParams.put("decoys_tag", params.getDecoysTag());
			sageParams.put("max_q", params.getMaxQ());
			json.put("sage_parameters", sageParams);
			LOGGER.info("Added Sage-specific parameters");

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

			// Metadata - convert String arrays to JSON
			JSONArray metadata = new JSONArray();
			if (params.getMetadataTableData() != null) {
				for (String[] row : params.getMetadataTableData()) {
					// MetadataData format: [sampleName, labelType, age, gender, condition, treatment, timePoint, customField1, customField2, customField3]
					JSONObject meta = new JSONObject();
					meta.put("sample_name", row.length > 0 ? row[0] : "");
					meta.put("label_type", row.length > 1 ? row[1] : "");
					meta.put("age", row.length > 2 ? row[2] : "");
					meta.put("gender", row.length > 3 ? row[3] : "");
					meta.put("condition", row.length > 4 ? row[4] : "");
					meta.put("treatment", row.length > 5 ? row[5] : "");
					meta.put("time_point", row.length > 6 ? row[6] : "");
					meta.put("custom_field_1", row.length > 7 ? row[7] : "");
					meta.put("custom_field_2", row.length > 8 ? row[8] : "");
					meta.put("custom_field_3", row.length > 9 ? row[9] : "");
					metadata.put(meta);
				}
			}
			json.put("metadata", metadata);

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
		return new DDAMagSageParameters();
	}

	@Override
	protected Object loadParametersFromFile(File inputFile) {
		return loadParameters(inputFile);
	}
}

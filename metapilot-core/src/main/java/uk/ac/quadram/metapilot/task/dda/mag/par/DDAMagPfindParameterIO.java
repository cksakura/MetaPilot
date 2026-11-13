package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.task.par.IsobaricLabelSettings;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.SampleMetadata;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

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
		return "DDA_MAG_PFIND";
	}

	/**
	 * Public convenience loader for CLI/headless usage, mirroring
	 * {@code DiaMarkerParameterIO.loadParameters(File)} /
	 * {@code DDAMagFragPipeParameterIO.loadParameters(File)}. Delegates to the
	 * existing protected {@link #loadParametersFromFile(File)} implementation.
	 */
	public DDAMagPfindParameters loadParameters(File inputFile) {
		Object o = loadParametersFromFile(inputFile);
		return (o instanceof DDAMagPfindParameters) ? (DDAMagPfindParameters) o : null;
	}

	@Override
	public Object loadLastUsedParameters() {
		File inputFile = new File(CONFIG_DIR, getLastUsedParametersFileName());
		if (!inputFile.exists()) {
			LOGGER.info("Last used parameters file not found: " + inputFile.getAbsolutePath());
			return null;
		}
		return loadParametersFromFile(inputFile);
	}

	@Override
	public Object loadAnalysisParameters(String resultDirectory) {
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
	}

	private static double requiredDouble(JSONObject object, String key) {
		Object value = object.get(key);
		try {
			double parsed = Double.parseDouble(String.valueOf(value));
			if (Double.isNaN(parsed) || Double.isInfinite(parsed)) {
				throw new NumberFormatException(String.valueOf(value));
			}
			return parsed;
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException(key + " must be a number (got " + value + ")");
		}
	}

	private static int requiredInt(JSONObject object, String key) {
		Object value = object.get(key);
		try {
			double parsed = Double.parseDouble(String.valueOf(value));
			if (parsed != Math.rint(parsed) || parsed > Integer.MAX_VALUE || parsed < Integer.MIN_VALUE) {
				throw new NumberFormatException(String.valueOf(value));
			}
			return (int) parsed;
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException(key + " must be an integer (got " + value + ")");
		}
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
			additionalParams.put("open_search", params.getOsp() != null && params.getOsp().isOpenSearch());
			additionalParams.put("genome_selector_mode", params.getGenomeSelectorMode());
			additionalParams.put("genome_selection_strategy", params.getGenomeSelectionStrategy());
			additionalParams.put("genome_selection_psm_coverage", params.getGenomeSelectionPsmCoverage());
			additionalParams.put("genome_selection_step_guard", params.getGenomeSelectionStepGuard());
			additionalParams.put("genome_selection_max_genomes", params.getGenomeSelectionMaxGenomes());
			additionalParams.put("generate_sample_specific_db", params.isReduceDb());
			additionalParams.put("spectra_clustering", params.isSpectraCluster());
			params.getPangenomeVariantSettings().save(additionalParams);
			additionalParams.put("undecidable_tolerance_da", params.getUndecidableToleranceDa());
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

			String mode = params.getQuantitativeMode();
			if (mode != null && mode.toLowerCase().contains("isotopic")) {
				quantMode.put("light_labels", new JSONArray(params.getLightLabels()));
				quantMode.put("medium_labels", new JSONArray(params.getMediumLabels()));
				quantMode.put("heavy_labels", new JSONArray(params.getHeavyLabels()));
			} else if (IsobaricLabelSettings.isIsobaricMode(mode)) {
				// Save the resolved tag name (e.g. "TMT10plex"), not the mode text, so the
				// tag survives a reload and the GUI can reselect it.
				IsobaricTag tag = params.getIsobaricTag();
				quantMode.put("isobaric_type", tag != null ? tag.getName() : params.getIsobaricType());
				quantMode.put("combine_label", params.isCombineLabel());

				JSONArray isobaricTable = new JSONArray();
				if (params.getIsobaricTableData() != null) {
					for (String[] row : params.getIsobaricTableData()) {
						if (row == null || row.length == 0) {
							continue;
						}
						JSONObject label = new JSONObject();
						label.put("label_name", row[0]);
						label.put("minus2", row.length > 1 ? row[1] : "0");
						label.put("minus1", row.length > 2 ? row[2] : "0");
						label.put("plus1", row.length > 3 ? row[3] : "0");
						label.put("plus2", row.length > 4 ? row[4] : "0");
						isobaricTable.put(label);
					}
				}
				quantMode.put("isobaric_table", isobaricTable);
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
				"DDA_MAG_PFIND",
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
				params.setDigestMode(enzymeSettings.optString("digestion_mode", "Fully-Specific"));
			}

			// Load additional parameters
			if (json.has("additional_parameters")) {
				JSONObject additionalParams = json.getJSONObject("additional_parameters");
				params.setThreadCount(additionalParams.optInt("threads", Runtime.getRuntime().availableProcessors()));
				params.setMs2ScanMode(additionalParams.optString("ms2_scan_mode", "FTMS"));
				if (params.getOsp() != null) {
					params.getOsp().setOpenSearch(additionalParams.optBoolean("open_search", false));
				}
				params.setGenomeSelectorMode(additionalParams.optString("genome_selector_mode", MagDbItem.MARKER2));
				String storedStrategy = additionalParams.optString("genome_selection_strategy", "Focus");
				String strategyWarning = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut
						.strategyReplacementWarning(storedStrategy);
				if (strategyWarning != null) {
					LOGGER.warn(strategyWarning);
				}
				params.setGenomeSelectionStrategy(
						uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.strategyOrDefault(storedStrategy));
				boolean coverageKey = additionalParams.has("genome_selection_psm_coverage");
				boolean guardKey = additionalParams.has("genome_selection_step_guard");
				boolean maxKey = additionalParams.has("genome_selection_max_genomes");
				uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.LoadedCoverage loaded =
						uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.readAdditionalKeys(additionalParams,
								params.getGenomeSelectionPsmCoverage(), params.getGenomeSelectionStepGuard(),
								params.getGenomeSelectionMaxGenomes());
				params.setGenomeSelectionPsmCoverage(loaded.coverage);
				params.setGenomeSelectionStepGuard(loaded.stepGuard);
				params.setGenomeSelectionMaxGenomes(loaded.maxGenomes);
				if (loaded.replacementWarning != null) {
					LOGGER.warn(loaded.replacementWarning);
				}
				String coverageWarning = uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.unusedCoverageKeysWarning(
						params.getGenomeSelectionStrategy(), params.getGenomeSelectionPsmCoverage(),
						params.getGenomeSelectionStepGuard(), params.getGenomeSelectionMaxGenomes());
				if (coverageWarning != null && (coverageKey || guardKey || maxKey)) {
					LOGGER.warn(coverageWarning);
				}
				params.setReduceDb(additionalParams.optBoolean("generate_sample_specific_db", true));
				params.setSpectraCluster(additionalParams.optBoolean("spectra_clustering", true));
				params.getPangenomeVariantSettings().load(additionalParams);
				params.setUndecidableToleranceDa(
						additionalParams.optDouble("undecidable_tolerance_da", params.getUndecidableToleranceDa()));
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
				if (quantMode.has("isobaric_table")) {
					JSONArray isobaricTable = quantMode.getJSONArray("isobaric_table");
					List<String[]> isobaricTableData = new ArrayList<>();
					for (int i = 0; i < isobaricTable.length(); i++) {
						JSONObject row = isobaricTable.getJSONObject(i);
						isobaricTableData.add(new String[] { row.optString("label_name", ""),
								row.optString("minus2", "0"), row.optString("minus1", "0"),
								row.optString("plus1", "0"), row.optString("plus2", "0") });
					}
					params.setIsobaricTableData(isobaricTableData);
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
				SampleMetadata.apply(params.getMetadata(), metadataTableData);
			}

			LOGGER.info("Successfully loaded parameters from {}", paramFile.getAbsolutePath());
			return params;

		} catch (IllegalArgumentException e) {
			LOGGER.error("Rejected parameters in " + paramFile.getAbsolutePath() + ": " + e.getMessage());
			throw e;
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

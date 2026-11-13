package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.SampleMetadata;

/**
 * Parameter I/O operations for DDA MAG FragPipe workflow
 *
 * This class handles saving and loading of DDA MAG FragPipe parameters
 * to/from JSON files for persistence and reproducibility.
 */
public class DDAMagFragPipeParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LogManager.getLogger(DDAMagFragPipeParameterIO.class);

    // Workflow-specific file names
    private static final String PARAMETER_FILENAME = "mag_fragpipe_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_mag_fragpipe_parameters.json";
    private static final String PREFERENCES_FILENAME = "mag_fragpipe_preferences.json";

    /** Marks files whose modification lists are written into the FragPipe workflow. */
    static final String MODIFICATIONS_APPLIED_KEY = "modifications_applied";

    private static List<String> stringList(JSONArray array) {
        List<String> list = new ArrayList<>();
        if (array != null) {
            for (int i = 0; i < array.length(); i++) {
                list.add(array.getString(i));
            }
        }
        return list;
    }

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
		return "DDA_MAG_FRAGPIPE";
	}


	/**
	 * Load DDA MAG FragPipe parameters from JSON file
	 */
	public DDAMagFragPipeParameters loadParameters(File inputFile) {
		try {
			LOGGER.info("Loading parameters from: {}", inputFile.getAbsolutePath());
			JSONObject json = readJsonFromFile(inputFile);
			
			MetaData metadata = new MetaData();
			if (json.has("input_output")) {
				JSONObject ioSettings = json.getJSONObject("input_output");
				if (ioSettings.has("experiments")) {
					JSONArray experiments = ioSettings.getJSONArray("experiments");
					String[] rawFiles = new String[experiments.length()];
					String[] expNames = new String[experiments.length()];
					int[] fractions = new int[experiments.length()];
					int[] replicates = new int[experiments.length()];

					for (int i = 0; i < experiments.length(); i++) {
						JSONObject exp = experiments.getJSONObject(i);
						rawFiles[i] = exp.optString("file_path", "");
						expNames[i] = exp.optString("experiment_name", "");
						fractions[i] = 1;
						replicates[i] = 1;
					}

					metadata = new MetaData(rawFiles, expNames, fractions, replicates);
					LOGGER.info("Loaded {} experiments from parameter file", experiments.length());
				} else if (ioSettings.has("input_dir")) {
					metadata = uk.ac.quadram.metapilot.task.par.HpcInputDiscovery.discoverMetaData(
							ioSettings.optString("input_dir", ""));
					LOGGER.info("Discovered {} raw files from input_dir", metadata.getRawFiles().length);
				}
			}

			// Load Input/Output settings
			String result = "";
			String microDb = "";
			String magDbIdentifier = "MGYG";
			String magDbVersion = "";
			boolean appendHostDb = false;
			String hostDB = "";
			if (json.has("input_output")) {
				JSONObject ioSettings = json.getJSONObject("input_output");
				result = ioSettings.optString("output_directory", "");
				microDb = ioSettings.optString("mag_database", "");
				magDbIdentifier = ioSettings.optString("mag_db_identifier", magDbIdentifier);
				magDbVersion = ioSettings.optString("mag_db_version", "");
				appendHostDb = ioSettings.optBoolean("append_host_db", false);
				hostDB = ioSettings.optString("host_database", "");
			}
			if (json.has("mag_database")) {
				JSONObject magDbSettings = json.getJSONObject("mag_database");
				microDb = magDbSettings.optString("catalogue_id", microDb);
				magDbIdentifier = magDbSettings.optString("identifier", magDbIdentifier);
				magDbVersion = magDbSettings.optString("catalogue_version", magDbVersion);
			}

			// Load additional parameters
			String ms2ScanMode = "HCD-FTMS";
			int threadCount = 4;
			int fragpipeParallelTasks = 0;
			String genomeSelectorMode = uk.ac.quadram.metapilot.task.mag.MagDbItem.MARKER2;
			String genomeSelectionStrategy = "Focus";
			String searchMode = "Closed";
			boolean openSearchAtMarkerStage = false;
			String customBasicWorkflowPath = "";
			String customCombinedWorkflowPath = "";
			if (json.has("additional_parameters")) {
				JSONObject additionalParams = json.getJSONObject("additional_parameters");
				ms2ScanMode = additionalParams.optString("ms2_scan_mode", "HCD-FTMS");
				threadCount = additionalParams.optInt("threads", 4);
				fragpipeParallelTasks = additionalParams.optInt("fragpipe_parallel_tasks", 0);
				genomeSelectorMode = additionalParams.optString("genome_selector_mode", genomeSelectorMode);
				genomeSelectionStrategy = additionalParams.optString("genome_selection_strategy", genomeSelectionStrategy);
				// Missing in older parameter files -> defaults reproduce the historical
				// closed-search-only behavior exactly (no migration needed).
				searchMode = additionalParams.optString("search_mode", "Closed");
				openSearchAtMarkerStage = additionalParams.optBoolean("open_search_marker_stage", false);
				customBasicWorkflowPath = additionalParams.optString("custom_basic_workflow_path", "");
				customCombinedWorkflowPath = additionalParams.optString("custom_combined_workflow_path", "");
			}

			// Create DDAMagFragPipeParameters
			DDAMagFragPipeParameters params = new DDAMagFragPipeParameters(
				"DDA_MAG_FRAGPIPE",
				metadata,
				result,
				microDb,
				appendHostDb,
				hostDB,
				ms2ScanMode,
				threadCount
			);

			params.setGenomeSelectorMode(genomeSelectorMode);
			String strategyWarning = GenomeSelectionCut.strategyReplacementWarning(genomeSelectionStrategy);
			if (strategyWarning != null) {
				LOGGER.warn(strategyWarning);
			}
			params.setGenomeSelectionStrategy(GenomeSelectionCut.strategyOrDefault(genomeSelectionStrategy));
			if (json.has("additional_parameters")) {
				JSONObject coverageKeys = json.getJSONObject("additional_parameters");
				GenomeSelectionCut.LoadedCoverage loaded = GenomeSelectionCut.readAdditionalKeys(coverageKeys,
						params.getGenomeSelectionPsmCoverage(), params.getGenomeSelectionStepGuard(),
						params.getGenomeSelectionMaxGenomes());
				params.setGenomeSelectionPsmCoverage(loaded.coverage);
				params.setGenomeSelectionStepGuard(loaded.stepGuard);
				params.setGenomeSelectionMaxGenomes(loaded.maxGenomes);
				if (loaded.replacementWarning != null) {
					LOGGER.warn(loaded.replacementWarning);
				}
				if (loaded.anyKey) {
					String warning = GenomeSelectionCut.unusedCoverageKeysWarning(params.getGenomeSelectionStrategy(),
							loaded.coverage, loaded.stepGuard, loaded.maxGenomes);
					if (warning != null) {
						LOGGER.warn(warning);
					}
				}
			}
			params.setFragpipeParallelTasks(fragpipeParallelTasks);
			params.setSearchMode(searchMode);
			params.setOpenSearchAtMarkerStage(openSearchAtMarkerStage);
			params.setCustomBasicWorkflowPath(customBasicWorkflowPath);
			params.setCustomCombinedWorkflowPath(customCombinedWorkflowPath);
			// Keep MAG DB info in MAG-specific fields for UI restore / HPC resolve
			params.setMagDbCatalogueId(microDb);
			params.setMagDbIdentifier(magDbIdentifier);
			params.setMagDbVersion(magDbVersion);
			if (microDb != null && !microDb.isBlank()) {
				try {
					uk.ac.quadram.metapilot.task.mag.MagDbItem magDbItem =
							uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.resolve(
									microDb, magDbIdentifier, magDbVersion, null);
					if (magDbItem != null) {
						params.setMagDbItem(magDbItem);
					}
				} catch (Exception e) {
					LOGGER.warn("Failed to resolve MAG database from config: {}", e.getMessage());
				}
			}

			// Load enzyme settings
			if (json.has("enzyme_settings")) {
				JSONObject enzymeSettings = json.getJSONObject("enzyme_settings");
				params.setEnzyme(migrateEnzymeName(enzymeSettings,
						enzymeSettings.optString("enzyme", MSFraggerHeadlessRunner.DEFAULT_ENZYME)));
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

			if (json.optBoolean(MODIFICATIONS_APPLIED_KEY, false)) {
				params.setFixedModifications(stringList(json.optJSONArray("fixed_modifications")));
				params.setVariableModifications(stringList(json.optJSONArray("variable_modifications")));
				LOGGER.info("Loaded {} fixed and {} variable modifications", params.getFixedModifications().size(),
						params.getVariableModifications().size());
			} else {
				// Saved before the modifications reached FragPipe: those runs searched
				// the workflow template's own modifications, so this file still does.
				params.setApplyModifications(false);
				LOGGER.warn("{} predates applying modifications to FragPipe; the workflow template's "
						+ "modifications are searched. Save the parameters again to use the listed ones.",
						inputFile.getName());
			}

			List<String[]> metadataRows = SampleMetadata.fromJson(json);
			if (metadataRows != null) {
				params.setMetadataTableData(metadataRows);
				SampleMetadata.apply(params.getMetadata(), metadataRows);
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
			ioSettings.put("mag_database", params.getMagDbCatalogueId() != null && !params.getMagDbCatalogueId().isBlank()
					? params.getMagDbCatalogueId() : params.getMicroDb());
			ioSettings.put("mag_db_identifier", params.getMagDbIdentifier());
			ioSettings.put("mag_db_version", params.getMagDbVersion());
			ioSettings.put("append_host_db", params.isAppendHostDb());
			ioSettings.put("host_database", params.getHostDB());
			LOGGER.info("Added I/O settings - output: {}, mag_db: {}, identifier: {}",
					params.getResult(), params.getMagDbCatalogueId(), params.getMagDbIdentifier());

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

			// Dedicated MAG database block (aligned with DIA Marker)
			JSONObject magDbSettings = new JSONObject();
			String catalogueId = params.getMagDbCatalogueId();
			if ((catalogueId == null || catalogueId.isBlank()) && params.getMicroDb() != null) {
				catalogueId = params.getMicroDb();
			}
			magDbSettings.put("catalogue_id", catalogueId != null ? catalogueId : "");
			magDbSettings.put("identifier", params.getMagDbIdentifier());
			magDbSettings.put("catalogue_version", params.getMagDbVersion() != null ? params.getMagDbVersion() : "");
			json.put("mag_database", magDbSettings);

			// Enzyme settings
			JSONObject enzymeSettings = new JSONObject();
			enzymeSettings.put("enzyme", params.getEnzyme());
			enzymeSettings.put(ENZYME_NAMING_KEY, ENZYME_NAMING_STANDARD);
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
			additionalParams.put("genome_selector_mode", params.getGenomeSelectorMode());
			additionalParams.put("genome_selection_strategy", params.getGenomeSelectionStrategy());
			GenomeSelectionCut.writeAdditionalKeys(additionalParams, params.getGenomeSelectionPsmCoverage(),
					params.getGenomeSelectionStepGuard(), params.getGenomeSelectionMaxGenomes());
			additionalParams.put("fragpipe_parallel_tasks", params.getFragpipeParallelTasks());
			additionalParams.put("search_mode", params.getSearchMode());
			additionalParams.put("open_search_marker_stage", params.isOpenSearchAtMarkerStage());
			additionalParams.put("custom_basic_workflow_path", params.getCustomBasicWorkflowPath());
			additionalParams.put("custom_combined_workflow_path", params.getCustomCombinedWorkflowPath());
			json.put("additional_parameters", additionalParams);
			LOGGER.info("Added additional parameters - threads: {}, ms2_scan_mode: {}", params.getThreadCount(), params.getMs2ScanMode());

			// Fixed modifications
			JSONArray fixedMods = new JSONArray();
			List<String> fixedModifications = params.getFixedModifications();
			if (fixedModifications != null) {
				for (String mod : fixedModifications) {
					fixedMods.put(mod);
				}
				LOGGER.info("Added {} fixed modifications", fixedModifications.size());
			}
			json.put("fixed_modifications", fixedMods);

			// Variable modifications
			JSONArray variableMods = new JSONArray();
			List<String> variableModifications = params.getVariableModifications();
			if (variableModifications != null) {
				for (String mod : variableModifications) {
					variableMods.put(mod);
				}
				LOGGER.info("Added {} variable modifications", variableModifications.size());
			}
			json.put("variable_modifications", variableMods);
			json.put(MODIFICATIONS_APPLIED_KEY, params.isApplyModifications());
			json.put(SampleMetadata.JSON_KEY, SampleMetadata.toJson(params.getMetadataTableData()));

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

	/**
	 * Marks a parameter file as using the corrected enzyme naming.
	 * <p>
	 * Before the correction, {@code "Trypsin"} in these files meant MSFragger's
	 * {@code stricttrypsin} - cleavage after K/R including before proline - and
	 * {@code "Trypsin/P"} meant the conventional proline-restricted digest. The
	 * names now mean what they say. A file written before the change therefore
	 * asks for a different digest than its author intended, and reloading one
	 * would silently re-run a study with a different enzyme.
	 * <p>
	 * The shared header's {@code version} field is hard-coded and never read
	 * back, so it cannot date a file. The presence of this key can: it is written
	 * from now on and absent from everything older.
	 */
	private static final String ENZYME_NAMING_KEY = "enzyme_naming";
	private static final String ENZYME_NAMING_STANDARD = "standard";

	/**
	 * Translates a pre-correction enzyme name to the one that means the same
	 * digest today. Only the two inverted pairs move; everything else is
	 * unchanged, and a file already carrying {@link #ENZYME_NAMING_KEY} is left
	 * exactly as written.
	 */
	private static String migrateEnzymeName(JSONObject enzymeSettings, String enzyme) {
		if (enzymeSettings.has(ENZYME_NAMING_KEY) || enzyme == null) {
			return enzyme;
		}
		String migrated;
		switch (enzyme.trim().toLowerCase(java.util.Locale.ROOT)) {
		case "trypsin":
			migrated = "Trypsin/P";
			break;
		case "trypsin/p":
		case "trypsin_p":
			migrated = "Trypsin";
			break;
		case "lysc":
			migrated = "LysC/P";
			break;
		case "lysc/p":
		case "lysc_p":
			migrated = "LysC";
			break;
		default:
			return enzyme;
		}
		LOGGER.warn("This parameter file predates the enzyme-naming correction, where \"{}\" meant the digest now "
				+ "called \"{}\". Loading it as \"{}\" so the study runs the enzyme it was written with. Re-save the "
				+ "file to make this permanent.", enzyme, migrated, migrated);
		return migrated;
	}

}

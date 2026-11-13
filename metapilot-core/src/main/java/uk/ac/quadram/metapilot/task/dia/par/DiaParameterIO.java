package uk.ac.quadram.metapilot.task.dia.par;

import java.io.File;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.SampleMetadata;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;

/**
 * Parameter I/O operations for DIA Seed workflow
 *
 * This class handles saving and loading of DIA parameters
 * to/from JSON files for persistence and reproducibility.
 */
public class DiaParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LogManager.getLogger(DiaParameterIO.class);

    // Workflow-specific file names (new seed naming)
    private static final String PARAMETER_FILENAME = "dia_seed_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_dia_seed_parameters.json";
    private static final String PREFERENCES_FILENAME = "dia_seed_preferences.json";

    // Legacy filenames kept for backward compatibility with existing runs
    private static final String LEGACY_PARAMETER_FILENAME = "dia_analysis_parameters.json";
    private static final String LEGACY_LAST_PARAMETERS_FILENAME = "last_used_dia_parameters.json";
    private static final String LEGACY_PREFERENCES_FILENAME = "dia_preferences.json";

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
		return uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType();
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

    public MetaParameterDia loadParameters(File inputFile) {
        Object o = loadParametersFromFile(inputFile);
        return (o instanceof MetaParameterDia) ? (MetaParameterDia) o : null;
    }

    @Override
    protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
        if (!(parameters instanceof MetaParameterDia)) {
            LOGGER.error("Invalid parameter type: expected MetaParameterDia");
            return false;
        }
        
        MetaParameterDia params = (MetaParameterDia) parameters;
        
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
            ioSettings.put("micro_db", params.getMicroDb());
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

            // DIA-specific settings (simplified: no self_model_search / use_model / models for now)
            JSONObject diaSettings = new JSONObject();
            diaSettings.put("library_search", params.isLibrarySearch());

            // Library information
            JSONArray libraries = new JSONArray();
            if (params.getLibrary() != null) {
                for (String lib : params.getLibrary()) {
                    libraries.put(lib);
                }
            }
            diaSettings.put("libraries", libraries);

            // Modification settings (for DIA-NN / library generation)
            diaSettings.put("mExcision", params.isMExcision());
            diaSettings.put("carbamiC", params.isCarbamiC());
            diaSettings.put("oxidationM", params.isOxidationM());
            diaSettings.put("acetyProN", params.isAcetyProN());
            diaSettings.put("phosphoSTY", params.isPhosphoSTY());
            diaSettings.put("kGG", params.isKGG());

            // Seed mode (DIRECT, FOCUS, COMPREHENSIVE, DEEP)
            String seedMode = params.getSeedMode();
            if (seedMode == null || seedMode.isEmpty()) {
                seedMode = "DIRECT";
            }
            diaSettings.put("seed_mode", seedMode);
            diaSettings.put("use_peptide_based_round", params.isUsePeptideBasedRound());
            diaSettings.put("peptide_detectability_threshold", params.getPeptideDetectabilityThreshold());

            json.put("dia_settings", diaSettings);

            // MAG database information
            JSONObject magDbSettings = new JSONObject();
            if (params.getMagDbItem() != null) {
                MagDbItem magDbItem = params.getMagDbItem();
                magDbSettings.put("catalogue_id", magDbItem.getCatalogueID());
                magDbSettings.put("identifier", magDbItem.getIdentifier());
                if (magDbItem.getUsedVersion() != null) {
                    magDbSettings.put("catalogue_version", magDbItem.getUsedVersion().getCatalogueVersion());
                }
            }
            json.put("mag_database", magDbSettings);

            // Basic parameters
            JSONObject basicParams = new JSONObject();
            basicParams.put("threads", params.getThreadCount());
            basicParams.put("core_count", params.getCoreCount());
            basicParams.put("ms2_scan_mode", params.getMs2ScanMode());
            basicParams.put("meta_workflow", params.isMetaWorkflow());
            json.put("basic_parameters", basicParams);

            json.put(SampleMetadata.JSON_KEY, SampleMetadata.toJson(metaData));

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
            SampleMetadata.apply(metadata, SampleMetadata.fromJson(json));

            // Create MetaBaseParameter first
            uk.ac.quadram.metapilot.task.par.MetaBaseParameter metaBasePar = new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
                "DiaNNMAG",
                metadata,
                "", // output_directory will be set below
                "", // micro_db will be set below
                false, // append_host_db will be set below
                "", // host_database will be set below
                "DDA", // ms2_scan_mode will be set below
                1, // core_count will be set below
                Runtime.getRuntime().availableProcessors(), // threads will be set below
                true // meta_workflow
            );

            // Load Input/Output settings
            if (json.has("input_output")) {
                JSONObject ioSettings = json.getJSONObject("input_output");
                metaBasePar.setResult(ioSettings.optString("output_directory", ""));
                metaBasePar.setMicroDb(ioSettings.optString("micro_db", ""));
                metaBasePar.setAppendHostDb(ioSettings.optBoolean("append_host_db", false));
                metaBasePar.setHostDB(ioSettings.optString("host_database", ""));
            }

            // Load basic parameters
            if (json.has("basic_parameters")) {
                JSONObject basicParams = json.getJSONObject("basic_parameters");
                metaBasePar.setThreadCount(basicParams.optInt("threads", Runtime.getRuntime().availableProcessors()));
                metaBasePar.setCoreCount(basicParams.optInt("core_count", 1));
                metaBasePar.setMs2ScanMode(basicParams.optString("ms2_scan_mode", "DDA"));
                metaBasePar.setMetaWorkflow(basicParams.optBoolean("meta_workflow", true));
            }

            // Load MAG database information and find the MagDbItem
            MagDbItem magDbItem = null;
            if (json.has("mag_database")) {
                JSONObject magDbSettings = json.getJSONObject("mag_database");
                String catalogueId = magDbSettings.optString("catalogue_id", "");
                String identifier = magDbSettings.optString("identifier", "");
                String catalogueVersion = magDbSettings.optString("catalogue_version", "");
                
                if (!catalogueId.isEmpty()) {
                    // Try to load MAG databases from config and find the matching one
                    try {
                        File configFile = new File("config/MAGDB_1_0.json");
                        if (!configFile.exists()) {
                            configFile = new File("MAGDB_1_0.json");
                        }
                        if (configFile.exists()) {
                            MagDbItem[][] magDbItems = MagDbConfigIO.parse(configFile);
                            if (magDbItems != null) {
                                // Search through all MAG databases
                                for (MagDbItem[] dbArray : magDbItems) {
                                    if (dbArray != null) {
                                        for (MagDbItem item : dbArray) {
                                            if (item.getCatalogueID().equals(catalogueId) && 
                                                item.getIdentifier().equals(identifier)) {
                                                magDbItem = item;
                                                
                                                // Set the used version if specified
                                                if (!catalogueVersion.isEmpty()) {
                                                    for (MagDbItem.MagVersion version : item.getVersionList()) {
                                                        if (version.getCatalogueVersion().equals(catalogueVersion)) {
                                                            item.setUsedVersion(version);
                                                            break;
                                                        }
                                                    }
                                                }
                                                break;
                                            }
                                        }
                                    }
                                    if (magDbItem != null) break;
                                }
                            }
                        }
                    } catch (Exception e) {
                        LOGGER.warn("Failed to load MAG database from config: " + e.getMessage());
                    }
                    if (magDbItem == null) {
                        // Site MAG registry (CLI / HPC), as DiaMarkerParameterIO does.
                        try {
                            magDbItem = uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.resolve(
                                    catalogueId, identifier, catalogueVersion, null);
                        } catch (Exception e) {
                            LOGGER.warn("Failed to resolve MAG database from registry: " + e.getMessage());
                        }
                    }
                }
            }

            // Load DIA-specific settings
            boolean librarySearch = false;
            String[] library = new String[0];

            // Modification defaults (match DIA Marker defaults)
            boolean mExcision = true;
            boolean carbamiC = true;
            boolean oxidationM = false;
            boolean acetyProN = false;
            boolean phosphoSTY = false;
            boolean kGG = false;

            // Seed mode default
            String seedMode = "DIRECT";
            boolean usePeptideBasedRound = false;
            double peptideDetectabilityThreshold = 0.5;

            if (json.has("dia_settings")) {
                JSONObject diaSettings = json.getJSONObject("dia_settings");
                librarySearch = diaSettings.optBoolean("library_search", false);

                // Load libraries
                if (diaSettings.has("libraries")) {
                    JSONArray libraries = diaSettings.getJSONArray("libraries");
                    library = new String[libraries.length()];
                    for (int i = 0; i < libraries.length(); i++) {
                        library[i] = libraries.getString(i);
                    }
                }

                // Load modification settings if present
                mExcision   = diaSettings.optBoolean("mExcision",   mExcision);
                carbamiC    = diaSettings.optBoolean("carbamiC",    carbamiC);
                oxidationM  = diaSettings.optBoolean("oxidationM",  oxidationM);
                acetyProN   = diaSettings.optBoolean("acetyProN",   acetyProN);
                phosphoSTY  = diaSettings.optBoolean("phosphoSTY",  phosphoSTY);
                kGG         = diaSettings.optBoolean("kGG",         kGG);

                // Load seed mode if present
                seedMode = diaSettings.optString("seed_mode", seedMode);
                usePeptideBasedRound = diaSettings.optBoolean("use_peptide_based_round", usePeptideBasedRound);
                peptideDetectabilityThreshold = diaSettings.optDouble("peptide_detectability_threshold",
                        peptideDetectabilityThreshold);
            }

            // Create MetaParameterDia (no self model / external models for now)
            MetaParameterDia params = new MetaParameterDia(
                metaBasePar,
                magDbItem,
                librarySearch,
                false,      // selfModelSearch (unused)
                library,
                false,      // useModel (unused)
                new String[0]  // models (unused)
            );

            // Apply modification settings to the parameters
            params.setMExcision(mExcision);
            params.setCarbamiC(carbamiC);
            params.setOxidationM(oxidationM);
            params.setAcetyProN(acetyProN);
            params.setPhosphoSTY(phosphoSTY);
            params.setKGG(kGG);

            // Apply seed mode
            params.setSeedMode(seedMode);
            params.setUsePeptideBasedRound(usePeptideBasedRound);
            params.setPeptideDetectabilityThreshold(peptideDetectabilityThreshold);

            LOGGER.info("Successfully loaded parameters from {}", paramFile.getAbsolutePath());
            return params;

        } catch (Exception e) {
            LOGGER.error("Failed to load parameters from " + paramFile.getAbsolutePath(), e);
            return createDefaultParameters();
        }
    }

    @Override
	protected Object createDefaultParameters() {
		// Create a default MetaParameterDia with minimal parameters
		uk.ac.quadram.metapilot.task.par.MetaBaseParameter metaBasePar = new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
				uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType(), new MetaData(), "", "",
				false, "", "DIA", 1, Runtime.getRuntime().availableProcessors(), true);

		return new MetaParameterDia(metaBasePar, null, // No MAG database by default
				false, // librarySearch
				false, // selfModelSearch
				new String[0], // library
				false, // useModel
				new String[0] // models
		);
	}
}

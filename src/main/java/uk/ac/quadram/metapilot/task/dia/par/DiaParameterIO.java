package uk.ac.quadram.metapilot.task.dia.par;

import java.io.File;

import org.json.JSONArray;
import org.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;

/**
 * Parameter I/O operations for DIA workflow
 *
 * This class handles saving and loading of DIA parameters
 * to/from JSON files for persistence and reproducibility.
 */
public class DiaParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LoggerFactory.getLogger(DiaParameterIO.class);

    // Workflow-specific file names
    private static final String PARAMETER_FILENAME = "dia_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_dia_parameters.json";
    private static final String PREFERENCES_FILENAME = "dia_preferences.json";

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
        return uk.ac.quadram.metapilot.task.MetaLabWorkflowType.DiaNNMAG.name();
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

            // DIA-specific settings
            JSONObject diaSettings = new JSONObject();
            diaSettings.put("library_search", params.isLibrarySearch());
            diaSettings.put("self_model_search", params.isSelfModelSearch());
            diaSettings.put("use_model", params.isUseModel());
            
            // Library information
            JSONArray libraries = new JSONArray();
            if (params.getLibrary() != null) {
                for (String lib : params.getLibrary()) {
                    libraries.put(lib);
                }
            }
            diaSettings.put("libraries", libraries);
            
            // Model information
            JSONArray models = new JSONArray();
            if (params.getModels() != null) {
                for (String model : params.getModels()) {
                    models.put(model);
                }
            }
            diaSettings.put("models", models);
            
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

            // Metadata
            JSONArray metadata = new JSONArray();
            if (metaData != null && metaData.getMetaInfo() != null) {
                String[][] metaInfo = metaData.getMetaInfo();
                for (int i = 0; i < metaInfo.length; i++) {
                    JSONObject meta = new JSONObject();
                    meta.put("sample_index", i);
                    if (metaInfo[i].length > 0) meta.put("age", metaInfo[i][0]);
                    if (metaInfo[i].length > 1) meta.put("gender", metaInfo[i][1]);
                    if (metaInfo[i].length > 2) meta.put("condition", metaInfo[i][2]);
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
                }
            }

            // Load DIA-specific settings
            boolean librarySearch = false;
            boolean selfModelSearch = false;
            boolean useModel = false;
            String[] library = new String[0];
            String[] models = new String[0];
            
            if (json.has("dia_settings")) {
                JSONObject diaSettings = json.getJSONObject("dia_settings");
                librarySearch = diaSettings.optBoolean("library_search", false);
                selfModelSearch = diaSettings.optBoolean("self_model_search", false);
                useModel = diaSettings.optBoolean("use_model", false);
                
                // Load libraries
                if (diaSettings.has("libraries")) {
                    JSONArray libraries = diaSettings.getJSONArray("libraries");
                    library = new String[libraries.length()];
                    for (int i = 0; i < libraries.length(); i++) {
                        library[i] = libraries.getString(i);
                    }
                }
                
                // Load models
                if (diaSettings.has("models")) {
                    JSONArray modelsArray = diaSettings.getJSONArray("models");
                    models = new String[modelsArray.length()];
                    for (int i = 0; i < modelsArray.length(); i++) {
                        models[i] = modelsArray.getString(i);
                    }
                }
            }

            // Create MetaParameterDia
            MetaParameterDia params = new MetaParameterDia(
                metaBasePar,
                magDbItem,
                librarySearch,
                selfModelSearch,
                library,
                useModel,
                models
            );

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
            uk.ac.quadram.metapilot.task.MetaLabWorkflowType.DiaNNMAG.name(),
            new MetaData(),
            "",
            "",
            false,
            "",
            "DDA",
            1,
            Runtime.getRuntime().availableProcessors(),
            true
        );
        
        return new MetaParameterDia(
            metaBasePar,
            null, // No MAG database by default
            false, // librarySearch
            false, // selfModelSearch
            new String[0], // library
            false, // useModel
            new String[0] // models
        );
    }
}

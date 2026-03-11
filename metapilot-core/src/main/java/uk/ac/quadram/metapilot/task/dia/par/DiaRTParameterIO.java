package uk.ac.quadram.metapilot.task.dia.par;

import java.io.File;

import org.json.JSONArray;
import org.json.JSONObject;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Parameter I/O operations for DIA Real-Time workflow
 *
 * This class handles saving and loading of DIA Real-Time parameters
 * to/from JSON files for persistence and reproducibility.
 */
public class DiaRTParameterIO extends BaseParameterIO {

    private static final Logger LOGGER = LogManager.getLogger(DiaRTParameterIO.class);

    // Workflow-specific file names
    private static final String PARAMETER_FILENAME = "dia_rt_analysis_parameters.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_used_dia_rt_parameters.json";
    private static final String PREFERENCES_FILENAME = "dia_rt_preferences.json";

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
        return uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_REALTIME.getAnalysisType();
    }

    @Override
    protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
        System.out.println("DEBUG: DiaRTParameterIO.saveParametersToFile() called");
        System.out.println("DEBUG: Parameters type: " + (parameters != null ? parameters.getClass().getSimpleName() : "NULL"));
        System.out.println("DEBUG: Output file: " + outputFile.getAbsolutePath());
        System.out.println("DEBUG: Description: " + description);
        
        if (!(parameters instanceof MetaParameterDiaRT)) {
            LOGGER.error("Invalid parameter type: expected MetaParameterDiaRT, got " + 
                        (parameters != null ? parameters.getClass().getSimpleName() : "null"));
            System.err.println("DEBUG: Parameter type validation FAILED");
            return false;
        }
        
        System.out.println("DEBUG: Parameter type validation PASSED");
        MetaParameterDiaRT params = (MetaParameterDiaRT) parameters;
        
        try {
            System.out.println("DEBUG: Starting JSON creation process...");
            
            // Ensure parent directories exist before writing the file
            File parentDir = outputFile.getParentFile();
            System.out.println("DEBUG: Parent directory: " + (parentDir != null ? parentDir.getAbsolutePath() : "NULL"));
            if (parentDir != null && !parentDir.exists()) {
                System.out.println("DEBUG: Creating parent directory...");
                boolean created = parentDir.mkdirs();
                if (!created) {
                    LOGGER.error("Failed to create parent directory: " + parentDir.getAbsolutePath());
                    System.err.println("DEBUG: Failed to create parent directory");
                    return false;
                }
                System.out.println("DEBUG: Parent directory created successfully");
            } else {
                System.out.println("DEBUG: Parent directory already exists");
            }

            System.out.println("DEBUG: Creating JSON header...");
            JSONObject json = createJsonHeader(description);
            json.put("workflowType", getWorkflowType());
            json.put("timestamp", System.currentTimeMillis());
            System.out.println("DEBUG: JSON header created successfully");
            
            // Add base parameters
            System.out.println("DEBUG: Processing base parameters...");
            System.out.println("DEBUG: Base parameters exist, creating JSON...");
            JSONObject baseParams = new JSONObject();
            baseParams.put("workflowType", params.getWorkflowTypeString());
            baseParams.put("result", params.getResult());
            baseParams.put("microDb", params.getMicroDb());
            baseParams.put("appendHostDb", params.isAppendHostDb());
            baseParams.put("hostDB", params.getHostDB());
            baseParams.put("ms2ScanMode", params.getMs2ScanMode());
            baseParams.put("coreCount", params.getCoreCount());
            baseParams.put("threadCount", params.getThreadCount());
            baseParams.put("metaWorkflow", params.isMetaWorkflow());
            System.out.println("DEBUG: Base parameters JSON created successfully");
            
            // Add metadata if available
            if (params.getMetadata() != null) {
                MetaData metadata = params.getMetadata();
                    JSONObject metadataJson = new JSONObject();
                    
                    if (metadata.getRawFiles() != null) {
                        JSONArray rawFiles = new JSONArray();
                        for (String file : metadata.getRawFiles()) {
                            rawFiles.put(file);
                        }
                        metadataJson.put("rawFiles", rawFiles);
                    }
                    
                    if (metadata.getExpNames() != null) {
                        JSONArray expNames = new JSONArray();
                        for (String name : metadata.getExpNames()) {
                            expNames.put(name);
                        }
                        metadataJson.put("expNames", expNames);
                    }
                    
                    if (metadata.getFractions() != null) {
                        JSONArray fractions = new JSONArray();
                        for (int fraction : metadata.getFractions()) {
                            fractions.put(fraction);
                        }
                        metadataJson.put("fractions", fractions);
                    }
                    
                    if (metadata.getReplicates() != null) {
                        JSONArray replicates = new JSONArray();
                        for (int replicate : metadata.getReplicates()) {
                            replicates.put(replicate);
                        }
                        metadataJson.put("replicates", replicates);
                    }
                    
                    if (metadata.getMetaInfo() != null) {
                        JSONArray metaInfo = new JSONArray();
                        for (String[] row : metadata.getMetaInfo()) {
                            JSONArray metaRow = new JSONArray();
                            for (String cell : row) {
                                metaRow.put(cell);
                            }
                            metaInfo.put(metaRow);
                        }
                        metadataJson.put("metaInfo", metaInfo);
                    }
                    
                    baseParams.put("metadata", metadataJson);
                }
                
            json.put("baseParameters", baseParams);
            
            // Add MAG database information
            if (params.getMagDbItem() != null) {
                JSONObject magDbJson = new JSONObject();
                magDbJson.put("catalogueID", params.getMagDbItem().getCatalogueID());
                magDbJson.put("identifier", params.getMagDbItem().getIdentifier());
                magDbJson.put("repositoryPath", params.getMagDbItem().getRepositoryFile() != null ? 
                             params.getMagDbItem().getRepositoryFile().getAbsolutePath() : "");
                
                // Add version information
                if (params.getMagDbItem().getUsedVersion() != null) {
                    MagDbItem.MagVersion version = params.getMagDbItem().getUsedVersion();
                    JSONObject versionJson = new JSONObject();
                    versionJson.put("catalogueVersion", version.getCatalogueVersion());
                    versionJson.put("speciesCount", version.getSpeciesCount());
                    versionJson.put("localAvailable", version.isLocalAvailable());
                    magDbJson.put("usedVersion", versionJson);
                }
                
                json.put("magDatabase", magDbJson);
            }
            
            // Add DIA-specific parameters (simplified: librarySearch + libraries + modifications)
            System.out.println("DEBUG: Adding DIA-specific parameters...");
            json.put("librarySearch", params.isLibrarySearch());
            
            // Add library array
            if (params.getLibrary() != null) {
                JSONArray libraryArray = new JSONArray();
                for (String lib : params.getLibrary()) {
                    libraryArray.put(lib);
                }
                json.put("library", libraryArray);
            }

            // DIA-NN modification settings (unified across DIA workflows)
            json.put("mExcision", params.isMExcision());
            json.put("carbamiC", params.isCarbamiC());
            json.put("oxidationM", params.isOxidationM());
            json.put("acetyProN", params.isAcetyProN());
            json.put("phosphoSTY", params.isPhosphoSTY());
            json.put("kGG", params.isKGG());
            System.out.println("DEBUG: DIA-specific parameters added successfully");
            
            // Add DIA Real-Time specific parameters
            System.out.println("DEBUG: Adding DIA Real-Time specific parameters...");
            json.put("combineResult", params.isCombineResult());
            json.put("stopNow", params.isStopNow());
            json.put("stopComplete", params.isStopComplete());
            json.put("taskCount", params.getTaskCount());
            // Real-time workflow mode (marker / seed variants)
            if (params.getRtWorkflowMode() != null) {
                json.put("rtWorkflowMode", params.getRtWorkflowMode());
            }
            
            // Add monitor and result folders
            if (params.getMonitorFolder() != null) {
                json.put("monitorFolder", params.getMonitorFolder().getAbsolutePath());
            }
            if (params.getResultFolder() != null) {
                json.put("resultFolder", params.getResultFolder().getAbsolutePath());
            }
            System.out.println("DEBUG: DIA Real-Time specific parameters added successfully");
            
            System.out.println("DEBUG: All JSON data prepared, writing to file...");
            boolean writeResult = writeJsonToFile(json, outputFile);
            System.out.println("DEBUG: writeJsonToFile() returned: " + writeResult);
            
            if (writeResult) {
                System.out.println("DEBUG: Parameters saved successfully to: " + outputFile.getAbsolutePath());
            } else {
                System.err.println("DEBUG: Failed to write JSON to file");
            }
            
            return writeResult;
            
        } catch (Exception e) {
            LOGGER.error("Failed to save parameters to " + outputFile.getAbsolutePath(), e);
            System.err.println("DEBUG: Exception in saveParametersToFile: " + e.getClass().getSimpleName() + " - " + e.getMessage());
            e.printStackTrace();
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
		// Create a default MetaParameterDiaRT with minimal parameters
		uk.ac.quadram.metapilot.task.par.MetaBaseParameter metaBasePar = new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
				uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_REALTIME.getAnalysisType(), new MetaData(), "",
				"", false, "", "DIA", 1, Runtime.getRuntime().availableProcessors(), true);

		return new MetaParameterDiaRT(metaBasePar, null, true, false, new String[0], false, new String[0]);
	}

    private boolean validateParameters(Object parameters) {
        if (!(parameters instanceof MetaParameterDiaRT)) {
            LOGGER.error("Invalid parameter type: expected MetaParameterDiaRT, got {}", 
                        parameters != null ? parameters.getClass().getSimpleName() : "null");
            return false;
        }
        
        MetaParameterDiaRT params = (MetaParameterDiaRT) parameters;
        
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
        if (!(parameters instanceof MetaParameterDiaRT)) {
            LOGGER.error("Invalid parameter type: expected MetaParameterDiaRT, got {}", 
                        parameters != null ? parameters.getClass().getSimpleName() : "null");
            return false;
        }
        
        MetaParameterDiaRT params = (MetaParameterDiaRT) parameters;
        
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
                String[][] metaInfo = null;
                
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
                
                if (metadataJson.has("metaInfo")) {
                    JSONArray metaInfoArray = metadataJson.getJSONArray("metaInfo");
                    metaInfo = new String[metaInfoArray.length()][];
                    for (int i = 0; i < metaInfoArray.length(); i++) {
                        JSONArray row = metaInfoArray.getJSONArray(i);
                        metaInfo[i] = new String[row.length()];
                        for (int j = 0; j < row.length(); j++) {
                            metaInfo[i][j] = row.getString(j);
                        }
                    }
                }
                
                metadata = new MetaData(rawFiles, expNames, fractions, replicates);
                if (metaInfo != null) {
                    metadata.setMetaInfo(metaInfo);
                }
            }
            
            uk.ac.quadram.metapilot.task.par.MetaBaseParameter baseParams = new uk.ac.quadram.metapilot.task.par.MetaBaseParameter(
                baseParamsJson.getString("workflowType"),
                metadata,
                baseParamsJson.getString("result"),
                baseParamsJson.getString("microDb"),
                baseParamsJson.getBoolean("appendHostDb"),
                baseParamsJson.getString("hostDB"),
                baseParamsJson.getString("ms2ScanMode"),
                baseParamsJson.getInt("coreCount"),
                baseParamsJson.getInt("threadCount"),
                baseParamsJson.getBoolean("metaWorkflow")
            );
            
            // Parse MAG database
            MagDbItem magDbItem = null;
            if (json.has("magDatabase")) {
                JSONObject magDbJson = json.getJSONObject("magDatabase");
                magDbItem = new MagDbItem(
                    magDbJson.getString("catalogueID"),
                    magDbJson.getString("identifier"),
                    new File(magDbJson.getString("repositoryPath"))
                );
                
                // Set version if available
                if (magDbJson.has("usedVersion")) {
                    JSONObject versionJson = magDbJson.getJSONObject("usedVersion");
                    // Use addVersion method to create the version properly
                    magDbItem.addVersion(
                        versionJson.getString("catalogueVersion"),
                        versionJson.getInt("speciesCount"),
                        new java.util.Date(), // Default date
                        versionJson.getBoolean("localAvailable")
                    );
                }
            }
            
            // Parse DIA-specific parameters
            boolean librarySearch = json.optBoolean("librarySearch", true);
            
            String[] library = new String[0];
            if (json.has("library")) {
                JSONArray libraryArray = json.getJSONArray("library");
                library = new String[libraryArray.length()];
                for (int i = 0; i < libraryArray.length(); i++) {
                    library[i] = libraryArray.getString(i);
                }
            }
            
            // Create MetaParameterDiaRT (no self model / external models for now)
            MetaParameterDiaRT params = new MetaParameterDiaRT(baseParams, magDbItem, librarySearch, false, library, false, new String[0]);
            
            // DIA-NN modification settings (default values match Seed workflow)
            boolean mExcision = json.optBoolean("mExcision", true);
            boolean carbamiC = json.optBoolean("carbamiC", true);
            boolean oxidationM = json.optBoolean("oxidationM", false);
            boolean acetyProN = json.optBoolean("acetyProN", false);
            boolean phosphoSTY = json.optBoolean("phosphoSTY", false);
            boolean kGG = json.optBoolean("kGG", false);

            params.setMExcision(mExcision);
            params.setCarbamiC(carbamiC);
            params.setOxidationM(oxidationM);
            params.setAcetyProN(acetyProN);
            params.setPhosphoSTY(phosphoSTY);
            params.setKGG(kGG);

            // Set DIA Real-Time specific parameters
            params.setCombineResult(json.optBoolean("combineResult", false));
            params.setStopNow(json.optBoolean("stopNow", false));
            params.setStopComplete(json.optBoolean("stopComplete", false));
            params.setTaskCount(json.optInt("taskCount", 0));
            // Real-time workflow mode; default to SEED_DIRECT if missing
            String rtMode = json.optString("rtWorkflowMode", null);
            if (rtMode == null || rtMode.isEmpty()) {
                rtMode = "SEED_DIRECT";
            }
            params.setRtWorkflowMode(rtMode);
            
            // Set monitor and result folders
            if (json.has("monitorFolder")) {
                params.setMonitorFolder(new File(json.getString("monitorFolder")));
            }
            if (json.has("resultFolder")) {
                params.setResultFolder(new File(json.getString("resultFolder")));
            }
            
            return params;
            
        } catch (Exception e) {
            LOGGER.error("Failed to parse DIA Real-Time parameters from JSON", e);
            throw new RuntimeException("Failed to parse DIA Real-Time parameters", e);
        }
    }
}

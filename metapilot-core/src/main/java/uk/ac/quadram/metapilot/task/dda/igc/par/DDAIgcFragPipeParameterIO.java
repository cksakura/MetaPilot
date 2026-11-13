package uk.ac.quadram.metapilot.task.dda.igc.par;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.SampleMetadata;

/**
 * Parameter IO Handler for DDA IGC FragPipe Workflow
 * 
 * Handles saving and loading of analysis parameters for reproducibility
 * and user convenience. Uses workflow-specific file naming to avoid
 * conflicts with other workflows.
 * 
 * File naming convention:
 * - Analysis-specific: "analysis_parameters_fragpipe.json"
 * - Last used: "config/last_dda_igc_fragpipe_parameters.json"
 * - Preferences: "~/.metapilot/fragpipe_preferences.json"
 */
public class DDAIgcFragPipeParameterIO extends BaseParameterIO {
    
    private static final Logger LOGGER = LogManager.getLogger(DDAIgcFragPipeParameterIO.class);
    
    // Workflow-specific file names to avoid conflicts
    private static final String PARAMETER_FILENAME = "analysis_parameters_fragpipe.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_dda_igc_fragpipe_parameters.json";
    private static final String PREFERENCES_FILENAME = "fragpipe_preferences.json";

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
        return uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_IGC_FRAGPIPE.getAnalysisType();
    }

    @Override
    protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
        if (!(parameters instanceof DDAIgcFragPipeParameters)) {
            LOGGER.error("Invalid parameter type: expected DDAIgcFragPipeParameters");
            return false;
        }
        
        DDAIgcFragPipeParameters params = (DDAIgcFragPipeParameters) parameters;
        
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
            List<String[]> experimentData = params.getExperimentData();
            for (String[] row : experimentData) {
                JSONObject exp = new JSONObject();
                exp.put("file_path", row[0]);
                exp.put("experiment_name", row[1]);
                experiments.put(exp);
            }
            ioSettings.put("experiments", experiments);
            
            json.put("input_output", ioSettings);
            
            // FragPipe-specific parameters
            JSONObject fragpipeParams = new JSONObject();
            fragpipeParams.put("enzyme", params.getEnzyme());
            fragpipeParams.put("missed_cleavages", params.getMissCleavages());
            fragpipeParams.put("digestion_mode", params.getDigestModeString());
            fragpipeParams.put("ms2_scan_mode", params.getMs2ScanMode());
            fragpipeParams.put("thread_count", params.getThreadCount());
            fragpipeParams.put("quantitative_mode", params.getQuantitativeMode());
            fragpipeParams.put("reduce_database", params.isReduceDb());
            fragpipeParams.put("spectra_clustering", params.isSpectraCluster());
            
            // Modifications
            JSONArray fixedMods = new JSONArray();
            String[] fixMods = params.getFixMods();
            if (fixMods != null) {
                for (String mod : fixMods) {
                    fixedMods.put(mod);
                }
            }
            fragpipeParams.put("fixed_modifications", fixedMods);
            
            JSONArray variableMods = new JSONArray();
            String[] variMods = params.getVariMods();
            if (variMods != null) {
                for (String mod : variMods) {
                    variableMods.put(mod);
                }
            }
            fragpipeParams.put("variable_modifications", variableMods);
            
            // Isobaric labeling
            fragpipeParams.put("isobaric_type", params.getIsobaricType());
            fragpipeParams.put("combine_label", params.isCombineLabel());
            
            JSONArray isobaricTable = new JSONArray();
            List<String[]> isobaricData = params.getIsobaricTableData();
            for (String[] row : isobaricData) {
                JSONObject isoRow = new JSONObject();
                isoRow.put("label", row[0]);
                isoRow.put("experiment", row[0]);
                isobaricTable.put(isoRow);
            }
            fragpipeParams.put("isobaric_table", isobaricTable);
            
            // Isotopic labeling
            JSONArray lightLabels = new JSONArray();
            for (String label : params.getLightLabels()) {
                lightLabels.put(label);
            }
            fragpipeParams.put("light_labels", lightLabels);
            
            JSONArray mediumLabels = new JSONArray();
            for (String label : params.getMediumLabels()) {
                mediumLabels.put(label);
            }
            fragpipeParams.put("medium_labels", mediumLabels);
            
            JSONArray heavyLabels = new JSONArray();
            for (String label : params.getHeavyLabels()) {
                heavyLabels.put(label);
            }
            fragpipeParams.put("heavy_labels", heavyLabels);
            
            json.put("fragpipe_parameters", fragpipeParams);
            
            // Metadata
            JSONArray metadata = new JSONArray();
            List<String[]> metadataData = params.getMetadataTableData();
            for (String[] row : metadataData) {
                JSONObject metaRow = new JSONObject();
                metaRow.put("sample_name", row[0]);
                metaRow.put("label_type", row[1]);
                metaRow.put("age", row[2]);
                metaRow.put("gender", row[3]);
                metaRow.put("condition", row[4]);
                metaRow.put("treatment", row[5]);
                metaRow.put("time_point", row[6]);
                metaRow.put("custom_field_1", row[7]);
                metaRow.put("custom_field_2", row[8]);
                metaRow.put("custom_field_3", row[9]);
                metadata.put(metaRow);
            }
            json.put("metadata", metadata);
            
            // Write to file
            try (FileWriter writer = new FileWriter(outputFile)) {
                writer.write(json.toString(2)); // Pretty print with 2-space indentation
                LOGGER.info("Successfully saved FragPipe parameters to: " + outputFile.getAbsolutePath());
                return true;
            }
            
        } catch (IOException e) {
            LOGGER.error("Failed to save FragPipe parameters to: " + outputFile.getAbsolutePath(), e);
            return false;
        }
    }

    @Override
    protected Object loadParametersFromFile(File paramFile) {
        LOGGER.info("Loading FragPipe parameters from: " + paramFile.getAbsolutePath());
        
        try {
            JSONObject json = readJsonFromFile(paramFile);
            return jsonToParameters(json);
        } catch (Exception e) {
            LOGGER.error("Failed to load FragPipe parameters from: " + paramFile.getAbsolutePath(), e);
            return createDefaultParameters();
        }
    }

    @Override
    protected Object createDefaultParameters() {
        return new DDAIgcFragPipeParameters();
    }

    /**
     * Convert JSON object to DDAIgcFragPipeParameters
     */
    private DDAIgcFragPipeParameters jsonToParameters(JSONObject json) {
        DDAIgcFragPipeParameters params = new DDAIgcFragPipeParameters();
        
        try {
            // Load input/output settings
            if (json.has("input_output")) {
                JSONObject io = json.getJSONObject("input_output");
                params.setResult(io.optString("output_directory", ""));
                params.setMicroDb(io.optString("microbiome_database", ""));
                params.setAppendHostDb(io.optBoolean("append_host_db", false));
                params.setHostDB(io.optString("host_database", ""));
                
                // Load experiments
                if (io.has("experiments")) {
                    JSONArray experiments = io.getJSONArray("experiments");
                    List<String[]> experimentData = new ArrayList<>();
                    for (int i = 0; i < experiments.length(); i++) {
                        JSONObject exp = experiments.getJSONObject(i);
                        String[] row = new String[]{
                            exp.optString("file_path", ""),
                            exp.optString("experiment_name", ""),
                            "1"
                        };
                        experimentData.add(row);
                    }
                    params.setExperimentData(experimentData);
                }
            }
            
            // Load FragPipe-specific parameters
            if (json.has("fragpipe_parameters")) {
                JSONObject fragpipe = json.getJSONObject("fragpipe_parameters");
                params.setEnzyme(fragpipe.optString("enzyme", "Trypsin"));
                params.setMissCleavages(fragpipe.optInt("missed_cleavages", 2));
                params.setDigestMode(fragpipe.optString("digestion_mode", "Specific"));
                params.setMs2ScanMode(fragpipe.optString("ms2_scan_mode", "HCD-FTMS"));
                params.setThreadCount(fragpipe.optInt("thread_count", Runtime.getRuntime().availableProcessors()));
                params.setQuantitativeMode(fragpipe.optString("quantitative_mode", "label-free"));
                params.setReduceDb(fragpipe.optBoolean("reduce_database", true));
                params.setSpectraCluster(fragpipe.optBoolean("spectra_clustering", false));
                
                // Load modifications
                if (fragpipe.has("fixed_modifications")) {
                    JSONArray fixedMods = fragpipe.getJSONArray("fixed_modifications");
                    String[] fixMods = new String[fixedMods.length()];
                    for (int i = 0; i < fixedMods.length(); i++) {
                        fixMods[i] = fixedMods.getString(i);
                    }
                    params.setFixMods(fixMods);
                }
                
                if (fragpipe.has("variable_modifications")) {
                    JSONArray variableMods = fragpipe.getJSONArray("variable_modifications");
                    String[] variMods = new String[variableMods.length()];
                    for (int i = 0; i < variableMods.length(); i++) {
                        variMods[i] = variableMods.getString(i);
                    }
                    params.setVariMods(variMods);
                }
                
                // Load isobaric labeling
                params.setIsobaricType(fragpipe.optString("isobaric_type", ""));
                params.setCombineLabel(fragpipe.optBoolean("combine_label", false));
                
                if (fragpipe.has("isobaric_table")) {
                    JSONArray isobaricTable = fragpipe.getJSONArray("isobaric_table");
                    List<String[]> isobaricData = new ArrayList<>();
                    for (int i = 0; i < isobaricTable.length(); i++) {
                        JSONObject isoRow = isobaricTable.getJSONObject(i);
                        String[] row = new String[5];
                        row[0] = isoRow.optString("label", "");
                        row[1] = "0"; // Default values for minus2, minus1, plus1, plus2
                        row[2] = "0";
                        row[3] = "0";
                        row[4] = "0";
                        isobaricData.add(row);
                    }
                    params.setIsobaricTableData(isobaricData);
                }
                
                // Load isotopic labeling
                if (fragpipe.has("light_labels")) {
                    JSONArray lightLabels = fragpipe.getJSONArray("light_labels");
                    List<String> lightList = new ArrayList<>();
                    for (int i = 0; i < lightLabels.length(); i++) {
                        lightList.add(lightLabels.getString(i));
                    }
                    params.setLightLabels(lightList);
                }
                
                if (fragpipe.has("medium_labels")) {
                    JSONArray mediumLabels = fragpipe.getJSONArray("medium_labels");
                    List<String> mediumList = new ArrayList<>();
                    for (int i = 0; i < mediumLabels.length(); i++) {
                        mediumList.add(mediumLabels.getString(i));
                    }
                    params.setMediumLabels(mediumList);
                }
                
                if (fragpipe.has("heavy_labels")) {
                    JSONArray heavyLabels = fragpipe.getJSONArray("heavy_labels");
                    List<String> heavyList = new ArrayList<>();
                    for (int i = 0; i < heavyLabels.length(); i++) {
                        heavyList.add(heavyLabels.getString(i));
                    }
                    params.setHeavyLabels(heavyList);
                }
            }
            
            // Load metadata
            if (json.has("metadata")) {
                JSONArray metadata = json.getJSONArray("metadata");
                List<String[]> metadataData = new ArrayList<>();
                for (int i = 0; i < metadata.length(); i++) {
                    JSONObject metaRow = metadata.getJSONObject(i);
                    String[] row = new String[10];
                    row[0] = metaRow.optString("sample_name", "");
                    row[1] = metaRow.optString("label_type", "");
                    row[2] = metaRow.optString("age", "");
                    row[3] = metaRow.optString("gender", "");
                    row[4] = metaRow.optString("condition", "");
                    row[5] = metaRow.optString("treatment", "");
                    row[6] = metaRow.optString("time_point", "");
                    row[7] = metaRow.optString("custom_field_1", "");
                    row[8] = metaRow.optString("custom_field_2", "");
                    row[9] = metaRow.optString("custom_field_3", "");
                    metadataData.add(row);
                }
                params.setMetadataTableData(metadataData);
                SampleMetadata.apply(params.getMetadata(), metadataData);
            }
            
        } catch (Exception e) {
            LOGGER.error("Error parsing FragPipe parameters from JSON", e);
        }
        
        return params;
    }

    /**
     * Create JSON header with metadata
     */
    public JSONObject createJsonHeader(String description) {
        JSONObject json = new JSONObject();
        json.put("workflow", getWorkflowType());
        json.put("version", "1.2.0");
        json.put("description", description);
        json.put("created", new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
        json.put("parameter_type", "DDAIgcFragPipeParameters");
        appendSoftwareVersion(json, description);
        return json;
    }
}

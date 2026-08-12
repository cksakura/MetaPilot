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

/**
 * Parameter IO Handler for DDA IGC Sage Workflow
 * 
 * Handles saving and loading of analysis parameters for reproducibility
 * and user convenience. Uses workflow-specific file naming to avoid
 * conflicts with other workflows.
 * 
 * File naming convention:
 * - Analysis-specific: "analysis_parameters_sage.json"
 * - Last used: "config/last_dda_igc_sage_parameters.json"
 * - Preferences: "~/.metapilot/sage_preferences.json"
 */
public class DDAIgcSageParameterIO extends BaseParameterIO {
    
    private static final Logger LOGGER = LogManager.getLogger(DDAIgcSageParameterIO.class);
    
    // Workflow-specific file names to avoid conflicts
    private static final String PARAMETER_FILENAME = "analysis_parameters_sage.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_dda_igc_sage_parameters.json";
    private static final String PREFERENCES_FILENAME = "sage_preferences.json";

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
        return uk.ac.quadram.metapilot.workflow.WorkflowType.DDA_IGC_SAGE.getAnalysisType();
    }

    @Override
    protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
        if (!(parameters instanceof DDAIgcSageParameters)) {
            LOGGER.error("Invalid parameter type: expected DDAIgcSageParameters");
            return false;
        }
        
        DDAIgcSageParameters params = (DDAIgcSageParameters) parameters;
        
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
            
            // Sage-specific parameters
            JSONObject sageParams = new JSONObject();
            sageParams.put("enzyme", params.getEnzyme());
            sageParams.put("missed_cleavages", params.getMissCleavages());
            sageParams.put("digestion_mode", params.getDigestModeString());
            sageParams.put("ms2_scan_mode", params.getMs2ScanMode());
            sageParams.put("thread_count", params.getThreadCount());
            sageParams.put("quantitative_mode", params.getQuantitativeMode());
            sageParams.put("reduce_database", params.isReduceDb());
            sageParams.put("spectra_clustering", params.isSpectraCluster());
            
            // Sage-specific search parameters
            sageParams.put("precursor_tolerance_min", params.getPrecursorToleranceMin());
            sageParams.put("precursor_tolerance_max", params.getPrecursorToleranceMax());
            sageParams.put("fragment_tolerance", params.getFragmentTolerance());
            sageParams.put("precursor_tolerance_unit", params.getPrecursorToleranceUnit());
            sageParams.put("fragment_tolerance_unit", params.getFragmentToleranceUnit());
            sageParams.put("min_peptide_length", params.getMinPeptideLength());
            sageParams.put("max_peptide_length", params.getMaxPeptideLength());
            sageParams.put("min_peak_intensity", params.getMinPeakIntensity());
            sageParams.put("max_variable_modifications", params.getMaxVariableModifications());
            sageParams.put("decoys_enabled", params.isDecoysEnabled());
            sageParams.put("generate_decoys", params.isGenerateDecoys());
            sageParams.put("decoys_tag", params.getDecoysTag());
            sageParams.put("max_q", params.getMaxQ());
            sageParams.put("write_checkpoint", params.isWriteCheckpoint());
            
            // Modifications
            JSONArray fixedMods = new JSONArray();
            String[] fixMods = params.getFixMods();
            if (fixMods != null) {
                for (String mod : fixMods) {
                    fixedMods.put(mod);
                }
            }
            sageParams.put("fixed_modifications", fixedMods);
            
            JSONArray variableMods = new JSONArray();
            String[] variMods = params.getVariMods();
            if (variMods != null) {
                for (String mod : variMods) {
                    variableMods.put(mod);
                }
            }
            sageParams.put("variable_modifications", variableMods);
            
            // Isobaric labeling
            sageParams.put("isobaric_type", params.getIsobaricType());
            sageParams.put("combine_label", params.isCombineLabel());
            
            JSONArray isobaricTable = new JSONArray();
            List<String[]> isobaricData = params.getIsobaricTableData();
            for (String[] row : isobaricData) {
                JSONObject isoRow = new JSONObject();
                isoRow.put("label", row[0]);
                isoRow.put("experiment", row[0]);
                isobaricTable.put(isoRow);
            }
            sageParams.put("isobaric_table", isobaricTable);
            
            // Isotopic labeling
            JSONArray lightLabels = new JSONArray();
            for (String label : params.getLightLabels()) {
                lightLabels.put(label);
            }
            sageParams.put("light_labels", lightLabels);
            
            JSONArray mediumLabels = new JSONArray();
            for (String label : params.getMediumLabels()) {
                mediumLabels.put(label);
            }
            sageParams.put("medium_labels", mediumLabels);
            
            JSONArray heavyLabels = new JSONArray();
            for (String label : params.getHeavyLabels()) {
                heavyLabels.put(label);
            }
            sageParams.put("heavy_labels", heavyLabels);
            
            json.put("sage_parameters", sageParams);
            
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
                LOGGER.info("Successfully saved Sage parameters to: " + outputFile.getAbsolutePath());
                return true;
            }
            
        } catch (IOException e) {
            LOGGER.error("Failed to save Sage parameters to: " + outputFile.getAbsolutePath(), e);
            return false;
        }
    }

    @Override
    protected Object loadParametersFromFile(File paramFile) {
        LOGGER.info("Loading Sage parameters from: " + paramFile.getAbsolutePath());
        
        try {
            JSONObject json = readJsonFromFile(paramFile);
            return jsonToParameters(json);
        } catch (Exception e) {
            LOGGER.error("Failed to load Sage parameters from: " + paramFile.getAbsolutePath(), e);
            return createDefaultParameters();
        }
    }

    @Override
    protected Object createDefaultParameters() {
        return new DDAIgcSageParameters();
    }

    /**
     * Convert JSON object to DDAIgcSageParameters
     */
    private DDAIgcSageParameters jsonToParameters(JSONObject json) {
        DDAIgcSageParameters params = new DDAIgcSageParameters();
        
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
            
            // Load Sage-specific parameters
            if (json.has("sage_parameters")) {
                JSONObject sage = json.getJSONObject("sage_parameters");
                params.setEnzyme(sage.optString("enzyme", "Trypsin"));
                params.setMissCleavages(sage.optInt("missed_cleavages", 2));
                params.setDigestMode(sage.optString("digestion_mode", "Specific"));
                params.setMs2ScanMode(sage.optString("ms2_scan_mode", "HCD-FTMS"));
                params.setThreadCount(sage.optInt("thread_count", Runtime.getRuntime().availableProcessors()));
                params.setQuantitativeMode(sage.optString("quantitative_mode", "label-free"));
                params.setReduceDb(sage.optBoolean("reduce_database", true));
                params.setSpectraCluster(sage.optBoolean("spectra_clustering", false));
                
                // Load Sage-specific search parameters
                params.setPrecursorToleranceMin(sage.optDouble("precursor_tolerance_min", -500.0));
                params.setPrecursorToleranceMax(sage.optDouble("precursor_tolerance_max", 150.0));
                params.setFragmentTolerance(sage.optDouble("fragment_tolerance", 20.0));
                params.setPrecursorToleranceUnit(sage.optString("precursor_tolerance_unit", "ppm"));
                params.setFragmentToleranceUnit(sage.optString("fragment_tolerance_unit", "ppm"));
                params.setMinPeptideLength(sage.optInt("min_peptide_length", 7));
                params.setMaxPeptideLength(sage.optInt("max_peptide_length", 50));
                params.setMinPeakIntensity(sage.optDouble("min_peak_intensity", 0.0));
                params.setMaxVariableModifications(sage.optInt("max_variable_modifications", 2));
                params.setDecoysEnabled(sage.optBoolean("decoys_enabled", true));
                params.setGenerateDecoys(sage.optBoolean("generate_decoys", true));
                params.setDecoysTag(sage.optString("decoys_tag", "rev_"));
                params.setMaxQ(sage.optDouble("max_q", 0.01));
                params.setWriteCheckpoint(sage.optBoolean("write_checkpoint", false));
                
                // Load modifications
                if (sage.has("fixed_modifications")) {
                    JSONArray fixedMods = sage.getJSONArray("fixed_modifications");
                    String[] fixMods = new String[fixedMods.length()];
                    for (int i = 0; i < fixedMods.length(); i++) {
                        fixMods[i] = fixedMods.getString(i);
                    }
                    params.setFixMods(fixMods);
                }
                
                if (sage.has("variable_modifications")) {
                    JSONArray variableMods = sage.getJSONArray("variable_modifications");
                    String[] variMods = new String[variableMods.length()];
                    for (int i = 0; i < variableMods.length(); i++) {
                        variMods[i] = variableMods.getString(i);
                    }
                    params.setVariMods(variMods);
                }
                
                // Load isobaric labeling
                params.setIsobaricType(sage.optString("isobaric_type", ""));
                params.setCombineLabel(sage.optBoolean("combine_label", false));
                
                if (sage.has("isobaric_table")) {
                    JSONArray isobaricTable = sage.getJSONArray("isobaric_table");
                    List<String[]> isobaricData = new ArrayList<>();
                    for (int i = 0; i < isobaricTable.length(); i++) {
                        JSONObject isoRow = isobaricTable.getJSONObject(i);
                        String[] row = new String[5];
                        row[0] = isoRow.optString("label", "");
                        row[1] = "0";
                        row[2] = "0";
                        row[3] = "0";
                        row[4] = "0";
                        isobaricData.add(row);
                    }
                    params.setIsobaricTableData(isobaricData);
                }
                
                // Load isotopic labeling
                if (sage.has("light_labels")) {
                    JSONArray lightLabels = sage.getJSONArray("light_labels");
                    List<String> lightList = new ArrayList<>();
                    for (int i = 0; i < lightLabels.length(); i++) {
                        lightList.add(lightLabels.getString(i));
                    }
                    params.setLightLabels(lightList);
                }
                
                if (sage.has("medium_labels")) {
                    JSONArray mediumLabels = sage.getJSONArray("medium_labels");
                    List<String> mediumList = new ArrayList<>();
                    for (int i = 0; i < mediumLabels.length(); i++) {
                        mediumList.add(mediumLabels.getString(i));
                    }
                    params.setMediumLabels(mediumList);
                }
                
                if (sage.has("heavy_labels")) {
                    JSONArray heavyLabels = sage.getJSONArray("heavy_labels");
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
            }
            
        } catch (Exception e) {
            LOGGER.error("Error parsing Sage parameters from JSON", e);
        }
        
        return params;
    }

    /**
     * Create JSON header with metadata
     */
    public JSONObject createJsonHeader(String description) {
        JSONObject json = new JSONObject();
        json.put("workflow", getWorkflowType());
        json.put("version", "1.1.1");
        json.put("description", description);
        json.put("created", new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
        json.put("parameter_type", "DDAIgcSageParameters");
        return json;
    }
}

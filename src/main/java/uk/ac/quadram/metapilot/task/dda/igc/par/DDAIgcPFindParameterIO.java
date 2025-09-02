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

import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;

/**
 * Parameter IO Handler for DDA IGC pFind Workflow
 * 
 * Handles saving and loading of analysis parameters for reproducibility
 * and user convenience. Uses workflow-specific file naming to avoid
 * conflicts with other workflows.
 * 
 * File naming convention:
 * - Analysis-specific: "analysis_parameters_pfind.json"
 * - Last used: "config/last_dda_igc_pfind_parameters.json"
 * - Preferences node: "uk/ac/quadram/metapilot/dda/igc/pfind"
 */
public class DDAIgcPFindParameterIO {
    
    private static final Logger LOGGER = LogManager.getLogger(DDAIgcPFindParameterIO.class);
    private static final String PREFERENCES_NODE = "uk/ac/quadram/metapilot/dda/igc/pfind";
    
    // Workflow-specific file names to avoid conflicts
    private static final String PARAMETER_FILENAME = "analysis_parameters_pfind.json";
    private static final String LAST_PARAMETERS_FILENAME = "config/last_dda_igc_pfind_parameters.json";
    
    /**
     * Save parameters to analysis result directory for reproducibility
     */
    public static boolean saveAnalysisParameters(String resultDirectory, DDAIgcPFindParameters parameters) {
        File paramFile = new File(resultDirectory, PARAMETER_FILENAME);
        return saveParametersToFile(paramFile, parameters, "pFind workflow analysis parameters");
    }
    
    /**
     * Save parameters as last used settings for user convenience
     */
    public static boolean saveLastUsedParameters(DDAIgcPFindParameters parameters) {
        File paramFile = new File(LAST_PARAMETERS_FILENAME);
        // Create config directory if it doesn't exist
        paramFile.getParentFile().mkdirs();
        return saveParametersToFile(paramFile, parameters, "Last used pFind workflow parameters");
    }
    
    /**
     * Load last used parameters for UI initialization
     */
    public static DDAIgcPFindParameters loadLastUsedParameters() {
        File paramFile = new File(LAST_PARAMETERS_FILENAME);
        if (!paramFile.exists()) {
            LOGGER.info("No last used pFind parameters found, using defaults");
            return createDefaultParameters();
        }
        
        try {
            return loadParametersFromFile(paramFile);
        } catch (Exception e) {
            LOGGER.warn("Failed to load last used pFind parameters, using defaults: " + e.getMessage());
            return createDefaultParameters();
        }
    }
    
    /**
     * Load parameters from analysis result directory
     */
    public static DDAIgcPFindParameters loadAnalysisParameters(String resultDirectory) {
        File paramFile = new File(resultDirectory, PARAMETER_FILENAME);
        if (!paramFile.exists()) {
            throw new IllegalArgumentException("pFind parameter file not found: " + paramFile.getAbsolutePath());
        }
        
        return loadParametersFromFile(paramFile);
    }
    
    private static boolean saveParametersToFile(File paramFile, DDAIgcPFindParameters parameters, String description) {
        try (FileWriter writer = new FileWriter(paramFile)) {
            JSONObject json = new JSONObject();
            
            // Header information
            json.put("workflow", "DDA_IGC_PFIND");
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
            List<ExperimentRow> experimentData = parameters.getExperimentData();
            for (ExperimentRow row : experimentData) {
                JSONObject exp = new JSONObject();
                exp.put("file_path", row.getFilePath());
                exp.put("experiment_name", row.getExperimentName());
                experiments.put(exp);
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
            modifications.put("fixed", new JSONArray(parameters.getFixedModifications()));
            modifications.put("variable", new JSONArray(parameters.getVariableModifications()));
            json.put("modifications", modifications);
            
            // Quantitative mode
            JSONObject quantMode = new JSONObject();
            quantMode.put("mode", parameters.getQuantitativeMode());
            
            if ("Isobaric".equals(parameters.getQuantitativeMode())) {
                quantMode.put("isobaric_type", parameters.getIsobaricType());
                
                JSONArray isobaricLabels = new JSONArray();
                for (IsobaricLabelRow row : parameters.getIsobaricTableData()) {
                    JSONObject label = new JSONObject();
                    label.put("label_name", row.getLabelName());
                    label.put("minus2", row.getMinus2());
                    label.put("minus1", row.getMinus1());
                    label.put("plus1", row.getPlus1());
                    label.put("plus2", row.getPlus2());
                    isobaricLabels.put(label);
                }
                quantMode.put("isobaric_table", isobaricLabels);
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
            
        } catch (IOException e) {
            LOGGER.error("Failed to save pFind parameters to " + paramFile.getAbsolutePath(), e);
            return false;
        }
    }
    
    private static DDAIgcPFindParameters loadParametersFromFile(File paramFile) {
        LOGGER.info("Loading pFind parameters from {}", paramFile.getAbsolutePath());
        
        try {
            // Read the JSON file
            String jsonContent = new String(java.nio.file.Files.readAllBytes(paramFile.toPath()));
            JSONObject json = new JSONObject(jsonContent);
            
            DDAIgcPFindParameters parameters = new DDAIgcPFindParameters();
            
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
                    List<ExperimentRow> experimentData = new ArrayList<>();
                    for (int i = 0; i < experiments.length(); i++) {
                        JSONObject exp = experiments.getJSONObject(i);
                        String filePath = exp.optString("file_path", "");
                        String experimentName = exp.optString("experiment_name", "");
                        if (!filePath.isEmpty() && !experimentName.isEmpty()) {
                            ExperimentRow row = new ExperimentRow(filePath, experimentName);
                            experimentData.add(row);
                        }
                    }
                    parameters.setExperimentData(experimentData);
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
                parameters.setThreadCount(additionalParams.optInt("threads", Runtime.getRuntime().availableProcessors()));
                parameters.setMs2ScanMode(additionalParams.optString("ms2_scan_mode", "HCD-FTMS"));
                parameters.setReduceDb(additionalParams.optBoolean("generate_sample_specific_db", true));
                parameters.setSpectraCluster(additionalParams.optBoolean("spectra_clustering", false));
            }
            
            // Load modifications
            if (json.has("modifications")) {
                JSONObject modSettings = json.getJSONObject("modifications");
                
                // Fixed modifications
                if (modSettings.has("fixed")) {
                    JSONArray fixedMods = modSettings.getJSONArray("fixed");
                    List<String> fixedModifications = new ArrayList<>();
                    for (int i = 0; i < fixedMods.length(); i++) {
                        fixedModifications.add(fixedMods.getString(i));
                    }
                    parameters.setFixedModifications(fixedModifications);
                }
                
                // Variable modifications
                if (modSettings.has("variable")) {
                    JSONArray variableMods = modSettings.getJSONArray("variable");
                    List<String> variableModifications = new ArrayList<>();
                    for (int i = 0; i < variableMods.length(); i++) {
                        variableModifications.add(variableMods.getString(i));
                    }
                    parameters.setVariableModifications(variableModifications);
                }
            }
            
            // Load quantitative mode
            if (json.has("quantitative_mode")) {
                JSONObject quantSettings = json.getJSONObject("quantitative_mode");
                parameters.setQuantitativeMode(quantSettings.optString("mode", "Label-free"));
                parameters.setIsobaricType(quantSettings.optString("isobaric_type", ""));
                
                // Load isobaric table data
                if (quantSettings.has("isobaric_table")) {
                    JSONArray isobaricTable = quantSettings.getJSONArray("isobaric_table");
                    List<IsobaricLabelRow> isobaricTableData = new ArrayList<>();
                    for (int i = 0; i < isobaricTable.length(); i++) {
                        JSONObject row = isobaricTable.getJSONObject(i);
                        IsobaricLabelRow labelRow = new IsobaricLabelRow(row.optString("label_name", ""));
                        labelRow.setMinus2(row.optString("minus2", "0"));
                        labelRow.setMinus1(row.optString("minus1", "0"));
                        labelRow.setPlus1(row.optString("plus1", "0"));
                        labelRow.setPlus2(row.optString("plus2", "0"));
                        isobaricTableData.add(labelRow);
                    }
                    parameters.setIsobaricTableData(isobaricTableData);
                }
            }
            
            // Load taxonomic analysis
            if (json.has("taxonomic_analysis")) {
                JSONObject taxSettings = json.getJSONObject("taxonomic_analysis");
                parameters.setEukaryota(taxSettings.optBoolean("eukaryota", true));
                parameters.setBacteria(taxSettings.optBoolean("bacteria", true));
                parameters.setArchaea(taxSettings.optBoolean("archaea", true));
                parameters.setViruses(taxSettings.optBoolean("viruses", true));
                parameters.setViroids(taxSettings.optBoolean("viroids", false));
                parameters.setUnclassified(taxSettings.optBoolean("unclassified", false));
                parameters.setOther(taxSettings.optBoolean("other", false));
                parameters.setIgnoreRank(taxSettings.optString("ignore_rank", "Family"));
                parameters.setIgnoreEnvironmentalSamples(taxSettings.optBoolean("ignore_environmental_samples", true));
            }
            
            // Load metadata table data
            if (json.has("metadata")) {
                JSONArray metadataData = json.getJSONArray("metadata");
                List<MetadataRow> metadataTableData = new ArrayList<>();
                for (int i = 0; i < metadataData.length(); i++) {
                    JSONObject row = metadataData.getJSONObject(i);
                    MetadataRow metadataRow = new MetadataRow(row.optString("sample_name", ""));
                    metadataRow.setLabelType(row.optString("label_type", ""));
                    metadataRow.setAge(row.optString("age", ""));
                    metadataRow.setGender(row.optString("gender", ""));
                    metadataRow.setCondition(row.optString("condition", ""));
                    metadataRow.setTreatment(row.optString("treatment", ""));
                    metadataRow.setTimePoint(row.optString("time_point", ""));
                    metadataRow.setCustomField1(row.optString("custom_field_1", ""));
                    metadataRow.setCustomField2(row.optString("custom_field_2", ""));
                    metadataRow.setCustomField3(row.optString("custom_field_3", ""));
                    metadataTableData.add(metadataRow);
                }
                parameters.setMetadataTableData(metadataTableData);
            }
            
            LOGGER.info("Successfully loaded pFind parameters from {}", paramFile.getAbsolutePath());
            return parameters;
            
        } catch (Exception e) {
            LOGGER.error("Failed to load pFind parameters from " + paramFile.getAbsolutePath(), e);
            return createDefaultParameters();
        }
    }
    
    private static DDAIgcPFindParameters createDefaultParameters() {
        return new DDAIgcPFindParameters();
    }
    
    /**
     * Save basic preferences (non-file specific settings) to Java preferences
     */
    public static void savePreferences(DDAIgcPFindParameters parameters) {
        Preferences prefs = Preferences.userRoot().node(PREFERENCES_NODE);
        
        // Save commonly reused settings
        prefs.put("default_enzyme", parameters.getEnzyme());
        prefs.putInt("default_missed_cleavages", parameters.getMissCleavages());
        prefs.put("default_digestion_mode", parameters.getDigestModeString());
        prefs.putInt("default_threads", parameters.getThreadCount());
        prefs.put("default_ms2_scan_mode", parameters.getMs2ScanMode());
        prefs.putBoolean("default_generate_sample_specific_db", parameters.isReduceDb());
        prefs.putBoolean("default_spectra_clustering", parameters.isSpectraCluster());
        
        LOGGER.info("Saved user preferences for DDA IGC pFind workflow");
    }
    
    /**
     * Load basic preferences from Java preferences
     */
    public static DDAIgcPFindParameters loadPreferences() {
        Preferences prefs = Preferences.userRoot().node(PREFERENCES_NODE);
        DDAIgcPFindParameters parameters = new DDAIgcPFindParameters();
        
        // Load commonly reused settings with defaults
        parameters.setEnzyme(prefs.get("default_enzyme", "Trypsin"));
        parameters.setMissCleavages(prefs.getInt("default_missed_cleavages", 2));
        parameters.setDigestMode(prefs.get("default_digestion_mode", "Specific"));
        parameters.setThreadCount(prefs.getInt("default_threads", Runtime.getRuntime().availableProcessors()));
        parameters.setMs2ScanMode(prefs.get("default_ms2_scan_mode", "HCD-FTMS"));
        parameters.setReduceDb(prefs.getBoolean("default_generate_sample_specific_db", true));
        parameters.setSpectraCluster(prefs.getBoolean("default_spectra_clustering", false));
        
        LOGGER.info("Loaded user preferences for DDA IGC pFind workflow");
        return parameters;
    }
    
    /**
     * Get the workflow-specific parameter file name for analysis results
     */
    public static String getAnalysisParameterFileName() {
        return PARAMETER_FILENAME;
    }
    
    /**
     * Get the workflow-specific last used parameter file name
     */
    public static String getLastUsedParameterFileName() {
        return LAST_PARAMETERS_FILENAME;
    }
    
    /**
     * Check if analysis parameters exist in the given directory
     */
    public static boolean hasAnalysisParameters(String resultDirectory) {
        File paramFile = new File(resultDirectory, PARAMETER_FILENAME);
        return paramFile.exists() && paramFile.length() > 0;
    }
    
    /**
     * Check if last used parameters exist
     */
    public static boolean hasLastUsedParameters() {
        File paramFile = new File(LAST_PARAMETERS_FILENAME);
        return paramFile.exists() && paramFile.length() > 0;
    }
}

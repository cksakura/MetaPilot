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

import uk.ac.quadram.metapilot.task.common.gui.MetaPilotIOPanelFX.ExperimentRow;
import uk.ac.quadram.metapilot.task.common.gui.IsobaricLabelingPanel.IsobaricLabelRow;
import uk.ac.quadram.metapilot.task.dda.igc.gui.MetaPilotDdaIgcMaxQuantPanel.MetadataRow;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;

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
 * - Preferences: "~/.metapilot/pfind_preferences.json"
 */
public class DDAIgcPFindParameterIO extends BaseParameterIO {
    
    private static final Logger LOGGER = LogManager.getLogger(DDAIgcPFindParameterIO.class);
    
    // Workflow-specific file names to avoid conflicts
    private static final String PARAMETER_FILENAME = "analysis_parameters_pfind.json";
    private static final String LAST_PARAMETERS_FILENAME = "last_dda_igc_pfind_parameters.json";
    private static final String PREFERENCES_FILENAME = "pfind_preferences.json";

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
		return uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DDA_IGC_PFIND.getAnalysisType();
	}

	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof DDAIgcPFindParameters)) {
			LOGGER.error("Invalid parameter type: expected DDAIgcPFindParameters");
			return false;
		}
		
		DDAIgcPFindParameters params = (DDAIgcPFindParameters) parameters;
		
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
            List<ExperimentRow> experimentData = params.getExperimentData();
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
            modifications.put("fixed", new JSONArray(params.getFixedModifications()));
            modifications.put("variable", new JSONArray(params.getVariableModifications()));
            json.put("modifications", modifications);
            
            // Quantitative mode
            JSONObject quantMode = new JSONObject();
            quantMode.put("mode", params.getQuantitativeMode());
            
            if ("Isobaric".equals(params.getQuantitativeMode())) {
                quantMode.put("isobaric_type", params.getIsobaricType());
                
                JSONArray isobaricLabels = new JSONArray();
                for (IsobaricLabelRow row : params.getIsobaricTableData()) {
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
            
            // Metadata
            JSONArray metadata = new JSONArray();
            for (MetadataRow row : params.getMetadataTableData()) {
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
            return writeJsonToFile(json, outputFile);
            
        } catch (Exception e) {
            LOGGER.error("Failed to save pFind parameters to " + outputFile.getAbsolutePath(), e);
            return false;
        }
	}

	@Override
	protected Object loadParametersFromFile(File paramFile) {
		LOGGER.info("Loading pFind parameters from {}", paramFile.getAbsolutePath());
        
        try {
            JSONObject json = readJsonFromFile(paramFile);
            
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

	@Override
	protected Object createDefaultParameters() {
		return new DDAIgcPFindParameters();
	}
}

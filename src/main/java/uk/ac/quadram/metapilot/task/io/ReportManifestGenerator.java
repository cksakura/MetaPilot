package uk.ac.quadram.metapilot.task.io;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.time.Instant;

import org.json.JSONArray;
import org.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Generates and updates manifest files for the MetaPilot report system.
 * This class handles both JSON and JS manifest formats for HTML reports.
 */
public class ReportManifestGenerator {
    
    private static final Logger LOGGER = LoggerFactory.getLogger(ReportManifestGenerator.class);
    
    private final File reportsDir;
    private final File manifestJsonFile;
    private final File manifestJsFile;
    
    public ReportManifestGenerator(File reportsDir) {
        this.reportsDir = reportsDir;
        this.manifestJsonFile = new File(reportsDir, "manifest.json");
        this.manifestJsFile = new File(reportsDir, "manifest.js");
    }
    
    /**
     * Adds a report entry to the manifest files.
     * 
     * @param title The report title
     * @param fileName The HTML file name (e.g., "Genome.html", "functions.html")
     * @param summary Optional summary description
     * @param tags Optional tags array
     */
    public void addReportEntry(String title, String fileName, String summary, String[] tags) {
        try {
            // Ensure reports directory exists
            if (!reportsDir.exists()) {
                reportsDir.mkdirs();
            }
            
            // Load existing manifest or create new one
            JSONObject manifest = loadExistingManifest();
            
            // Create new report entry
            JSONObject reportEntry = new JSONObject();
            reportEntry.put("title", title);
            reportEntry.put("file", "reports/" + fileName);
            reportEntry.put("created_at", getCurrentTimestamp());
            reportEntry.put("source", fileName);
            
            if (summary != null && !summary.trim().isEmpty()) {
                reportEntry.put("summary", summary);
            }
            
            if (tags != null && tags.length > 0) {
                JSONArray tagsArray = new JSONArray();
                for (String tag : tags) {
                    if (tag != null && !tag.trim().isEmpty()) {
                        tagsArray.put(tag.trim());
                    }
                }
                if (tagsArray.length() > 0) {
                    reportEntry.put("tags", tagsArray);
                }
            }
            
            // Update reports array
            JSONArray reports = manifest.getJSONArray("reports");
            
            // Remove existing entry with same source file
            for (int i = 0; i < reports.length(); i++) {
                JSONObject existingReport = reports.getJSONObject(i);
                if (fileName.equals(existingReport.optString("source", ""))) {
                    reports.remove(i);
                    break;
                }
            }
            
            // Add new entry
            reports.put(reportEntry);
            
            // Update generated timestamp
            manifest.put("generated_at", getCurrentTimestamp());
            
            // Write manifest files
            writeManifestFiles(manifest);
            
            LOGGER.info("Added report entry to manifest: {} -> {}", title, fileName);
            
        } catch (Exception e) {
            LOGGER.error("Error adding report entry to manifest: " + title, e);
        }
    }
    
    /**
     * Loads existing manifest or creates a new one.
     */
    private JSONObject loadExistingManifest() {
        if (manifestJsonFile.exists()) {
            try {
                String content = Files.readString(manifestJsonFile.toPath());
                return new JSONObject(content);
            } catch (IOException e) {
                LOGGER.warn("Could not read existing manifest, creating new one", e);
            }
        }
        
        // Create new manifest
        JSONObject manifest = new JSONObject();
        manifest.put("generated_at", getCurrentTimestamp());
        manifest.put("reports", new JSONArray());
        return manifest;
    }
    
    /**
     * Writes both JSON and JS manifest files.
     */
    private void writeManifestFiles(JSONObject manifest) throws IOException {
        // Write JSON manifest
        try (FileWriter writer = new FileWriter(manifestJsonFile)) {
            writer.write(manifest.toString(2));
        }
        
        // Write JS manifest
        try (FileWriter writer = new FileWriter(manifestJsFile)) {
            writer.write("window.METAPILOT_MANIFEST = ");
            writer.write(manifest.toString(2));
            writer.write(";\n");
        }
        
        LOGGER.debug("Updated manifest files: {} and {}", 
                    manifestJsonFile.getAbsolutePath(), 
                    manifestJsFile.getAbsolutePath());
    }
    
    /**
     * Gets current timestamp in ISO format.
     */
    private String getCurrentTimestamp() {
        return Instant.now().toString();
    }
    
    /**
     * Convenience method to add Genome.html report entry.
     */
    public void addGenomeReport(File genomeHtmlFile) {
        if (genomeHtmlFile != null && genomeHtmlFile.exists()) {
            addReportEntry(
                "Genome Analysis Report",
                "Genome.html",
                "Taxonomic analysis and genome-level insights from MAG workflow.",
                new String[]{"taxonomy", "genome", "MAG"}
            );
        }
    }
    
    /**
     * Convenience method to add functions.html report entry.
     */
    public void addFunctionsReport(File functionsHtmlFile) {
        if (functionsHtmlFile != null && functionsHtmlFile.exists()) {
            addReportEntry(
                "Protein Functional Analysis Report",
                "functions.html",
                "Functional annotation and pathway analysis from MAG workflow.",
                new String[]{"function", "annotation", "MAG"}
            );
        }
    }
    
    /**
     * Adds Python-generated report entries to manifest.
     * This method handles the standard Python reports generated by MetaPythonReportTask.
     */
    public void addPythonReport(String reportType, String fileName, String workflowType) {
        String title;
        String summary;
        String[] tags;
        
        switch (reportType) {
            case "ID_SUMMARY":
                title = "Identification Summary Report";
                summary = "Overall identification statistics and sample comparisons.";
                tags = new String[]{"identification", "summary", workflowType};
                break;
                
            case "PEPTIDE":
                title = "Peptides Summary Report";
                summary = "Peptide-level analysis with length, charge, and mass distributions.";
                tags = new String[]{"peptides", "analysis", workflowType};
                break;
                
            case "PROTEIN":
                title = "Proteins Summary Report";
                summary = "Protein-level analysis with coverage and peptide counts.";
                tags = new String[]{"proteins", "analysis", workflowType};
                break;
                
            case "TAXONOMY":
                title = "Taxonomy Summary Report";
                summary = "Taxonomic classification results and distributions.";
                tags = new String[]{"taxonomy", "classification", workflowType};
                break;
                
            case "FUNCTION":
                title = "Function Summary Report";
                summary = "Functional annotation results and pathway analysis.";
                tags = new String[]{"function", "annotation", workflowType};
                break;
                
            // TEMPORARILY DISABLED: Bubble chart report is hidden due to issues
            // case "BUBBLE_CHART":
            //     title = "Bubble Chart Report";
            //     summary = "Interactive hierarchical bubble chart showing taxonomic abundance.";
            //     tags = new String[]{"bubble-chart", "taxonomy", "hierarchical", workflowType};
            //     break;
                
            default:
                title = "Analysis Report";
                summary = "Generated analysis report from " + workflowType + " workflow.";
                tags = new String[]{workflowType, "analysis"};
                break;
        }
        
        addReportEntry(title, fileName, summary, tags);
    }
    
    /**
     * Adds multiple Python reports to manifest in batch.
     */
    public void addPythonReports(String[] reportTypes, String[] fileNames, String workflowType) {
        if (reportTypes.length != fileNames.length) {
            LOGGER.warn("Report types and file names arrays have different lengths");
            return;
        }
        
        for (int i = 0; i < reportTypes.length; i++) {
            addPythonReport(reportTypes[i], fileNames[i], workflowType);
        }
    }
    
    /**
     * Creates an initial empty manifest file to ensure the index page can load.
     * This should be called at the beginning of report generation.
     */
    public void createInitialManifest() {
        try {
            // Ensure reports directory exists
            if (!reportsDir.exists()) {
                reportsDir.mkdirs();
            }
            
            // Create empty manifest if it doesn't exist
            if (!manifestJsonFile.exists()) {
                JSONObject manifest = new JSONObject();
                manifest.put("generated_at", getCurrentTimestamp());
                manifest.put("reports", new JSONArray());
                
                writeManifestFiles(manifest);
                LOGGER.info("Created initial empty manifest file at: " + manifestJsonFile.getAbsolutePath());
            } else {
                LOGGER.info("Manifest file already exists at: " + manifestJsonFile.getAbsolutePath());
            }
        } catch (Exception e) {
            LOGGER.error("Error creating initial manifest file", e);
        }
    }
    
    /**
     * Creates a test manifest file with sample data for debugging.
     */
    public void createTestManifest() {
        try {
            // Ensure reports directory exists
            if (!reportsDir.exists()) {
                reportsDir.mkdirs();
            }
            
            JSONObject manifest = new JSONObject();
            manifest.put("generated_at", getCurrentTimestamp());
            
            JSONArray reports = new JSONArray();
            
            // Add a test report
            JSONObject testReport = new JSONObject();
            testReport.put("title", "Test Report");
            testReport.put("file", "test_report.html");
            testReport.put("summary", "This is a test report to verify the manifest system is working.");
            testReport.put("created_at", getCurrentTimestamp());
            testReport.put("tags", new JSONArray().put("test").put("debug"));
            
            reports.put(testReport);
            manifest.put("reports", reports);
            
            writeManifestFiles(manifest);
            LOGGER.info("Created test manifest file at: " + manifestJsonFile.getAbsolutePath());
        } catch (Exception e) {
            LOGGER.error("Error creating test manifest file", e);
        }
    }
}

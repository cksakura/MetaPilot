package uk.ac.quadram.metapilot.task;

import javafx.concurrent.Task;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;

/**
 * Python-based report copy and generation task.
 */
public class MetaPythonReportCopyTask extends Task<Boolean> {

    private static final Logger LOGGER = LogManager.getLogger(MetaPythonReportCopyTask.class);
    private static final String TASK_NAME = "Python report copy task";
    private static final SimpleDateFormat FORMAT = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

    private static final String REPORT_TEMPLATE_PATH = "src/main/java/uk/ac/quadram/metapilot/report";
    
    private final File resultDirectory;
    private final File reportDirectory;
    private final File htmlReportDirectory;
    private MetaPythonReportTask reportTask;

    public MetaPythonReportCopyTask(File resultDirectory) {
        this.resultDirectory = resultDirectory;
        this.reportDirectory = new File(resultDirectory, "report");
        this.htmlReportDirectory = new File(reportDirectory, "reports/rmd_html");
        this.reportTask = new MetaPythonReportTask(2);
        
        LOGGER.info("{}: Initialized for {}", TASK_NAME, resultDirectory.getAbsolutePath());
    }

    public boolean addReportData(String reportType, File dataFile) {
        if (dataFile == null || !dataFile.exists()) {
            LOGGER.warn("{}: Data file does not exist: {}", TASK_NAME, dataFile);
            return false;
        }

        String outputFileName = getReportFileName(reportType);
        File outputFile = new File(htmlReportDirectory, outputFileName);

        boolean added = reportTask.addTask(reportType, dataFile, outputFile);
        
        if (added) {
            LOGGER.info("{}: Added {} report task", TASK_NAME, reportType);
        }
        
        return added;
    }

    @Override
    protected Boolean call() throws Exception {
        try {
            updateMessage("Copying report template...");
            updateProgress(10, 100);

            if (!copyReportTemplate()) {
                return false;
            }

            updateMessage("Cleaning up directories...");
            updateProgress(30, 100);
            cleanupUnusedDirectories();

            updateMessage("Creating report directories...");
            updateProgress(40, 100);
            if (!htmlReportDirectory.exists()) {
                htmlReportDirectory.mkdirs();
            }

            updateMessage("Generating reports...");
            updateProgress(50, 100);

            Thread reportThread = new Thread(reportTask);
            reportThread.start();
            reportThread.join();

            Boolean reportSuccess = reportTask.get();
            if (!reportSuccess) {
                updateMessage("Report generation failed");
                return false;
            }

            updateMessage("Completed successfully");
            updateProgress(100, 100);

            LOGGER.info("{}: Successfully completed", TASK_NAME);
            return true;

        } catch (Exception e) {
            LOGGER.error("{}: Error during execution", TASK_NAME, e);
            updateMessage("Error: " + e.getMessage());
            return false;
        }
    }

    private boolean copyReportTemplate() {
        try {
            File sourceReportDir = new File(REPORT_TEMPLATE_PATH);
            
            if (!sourceReportDir.exists()) {
                LOGGER.error("{}: Source template not found: {}", TASK_NAME, sourceReportDir.getAbsolutePath());
                return false;
            }

            FileUtils.copyDirectory(sourceReportDir, reportDirectory);
            LOGGER.info("{}: Template copied successfully", TASK_NAME);
            return true;

        } catch (IOException e) {
            LOGGER.error("{}: Failed to copy template", TASK_NAME, e);
            return false;
        }
    }

    private void cleanupUnusedDirectories() {
        try {
            File metaMapDir = new File(reportDirectory, "reports/metamap");
            File shinyAppsDir = new File(reportDirectory, "reports/shinyapps");

            if (metaMapDir.exists()) {
                FileUtils.deleteDirectory(metaMapDir);
            }
            if (shinyAppsDir.exists()) {
                FileUtils.deleteDirectory(shinyAppsDir);
            }

        } catch (IOException e) {
            LOGGER.warn("{}: Failed to cleanup directories", TASK_NAME, e);
        }
    }

    private String getReportFileName(String reportType) {
        switch (reportType.toLowerCase()) {
            case "id_summary":
            case "summary":
                return "report_ID_summary.html";
            case "peptides":
            case "peptide":
                return "report_peptides_summary.html";
            case "proteins":
            case "protein":
                return "report_proteinGroups_summary.html";
            case "taxonomy":
            case "taxon":
                return "report_taxonomy_summary.html";
            case "function":
            case "functions":
                return "report_function_summary.html";
            default:
                return "report_" + reportType + "_summary.html";
        }
    }

    public String getTaskName() {
        return TASK_NAME;
    }

    public File getResultDirectory() {
        return resultDirectory;
    }

    public File getReportDirectory() {
        return reportDirectory;
    }

    public File getHtmlReportDirectory() {
        return htmlReportDirectory;
    }
} 
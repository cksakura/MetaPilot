package uk.ac.quadram.metapilot.task;

import javafx.concurrent.Task;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;

/**
 * Python-based report template copy task.
 * 
 * This task is responsible for copying report templates and setting up the directory structure
 * for report generation. It does NOT generate reports - that should be handled separately
 * by MetaPythonReportTask.
 */
public class MetaPythonReportCopyTask extends Task<Boolean> {

    private static final String TASK_NAME = "Python report copy task";
    protected Logger logger;
    
    private static final String REPORT_TEMPLATE_PATH = "src/main/java/uk/ac/quadram/metapilot/report";
    
    private final File resultDirectory;
    private final File reportDirectory;
    private final File htmlReportDirectory;

    public MetaPythonReportCopyTask(File resultDirectory) {
        this.resultDirectory = resultDirectory;
        this.reportDirectory = new File(resultDirectory, "report");
        this.htmlReportDirectory = new File(reportDirectory, "reports");
        
        this.logger = LogManager.getLogger(TASK_NAME);
        logger.info("{}: Initialized for {}", TASK_NAME, resultDirectory.getAbsolutePath());
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
            updateProgress(50, 100);
            cleanupUnusedDirectories();

            updateMessage("Creating report directories...");
            updateProgress(80, 100);
            if (!htmlReportDirectory.exists()) {
                htmlReportDirectory.mkdirs();
            }

            updateMessage("Template setup completed successfully");
            updateProgress(100, 100);

            logger.info("{}: Successfully completed template setup", TASK_NAME);
            return true;

        } catch (Exception e) {
        	logger.error("{}: Error during template setup", TASK_NAME, e);
            updateMessage("Error: " + e.getMessage());
            return false;
        }
    }

    private boolean copyReportTemplate() {
        try {
            File sourceReportDir = new File(REPORT_TEMPLATE_PATH);
            
            if (!sourceReportDir.exists()) {
            	logger.error("{}: Source template not found: {}", TASK_NAME, sourceReportDir.getAbsolutePath());
                return false;
            }

            FileUtils.copyDirectory(sourceReportDir, reportDirectory);
            logger.info("{}: Template copied successfully", TASK_NAME);
            return true;

        } catch (IOException e) {
        	logger.error("{}: Failed to copy template", TASK_NAME, e);
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
        	logger.warn("{}: Failed to cleanup directories", TASK_NAME, e);
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
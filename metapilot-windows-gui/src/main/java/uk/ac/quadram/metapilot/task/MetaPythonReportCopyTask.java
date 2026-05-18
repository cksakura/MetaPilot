package uk.ac.quadram.metapilot.task;

import javafx.concurrent.Task;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.report.ReportTemplateInstaller;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.File;

/**
 * JavaFX adapter for shared report template installation.
 */
public class MetaPythonReportCopyTask extends Task<Boolean> {

    private static final String TASK_NAME = "Python report copy task";
    protected final Logger logger;

    private final File resultDirectory;
    private File reportDirectory;
    private File htmlReportDirectory;

    @Deprecated
    public MetaPythonReportCopyTask(File resultDirectory) {
        this(resultDirectory, LogManager.getLogger(MetaPythonReportCopyTask.class));
    }

    public MetaPythonReportCopyTask(File resultDirectory, Logger logger) {
        this.logger = logger != null ? logger : LogManager.getLogger(MetaPythonReportCopyTask.class);
        this.resultDirectory = resultDirectory;
        this.reportDirectory = new File(resultDirectory, "report");
        this.htmlReportDirectory = new File(reportDirectory, "reports");
        this.logger.info("{}: Initialized for {}", TASK_NAME, resultDirectory.getAbsolutePath());
    }

    @Override
    protected Boolean call() {
        try {
            updateMessage("Copying report template...");
            updateProgress(10, 100);
            ReportTemplateInstaller.Result result = ReportTemplateInstaller.install(resultDirectory,
                    javafxProgressCallback());
            this.reportDirectory = result.reportDirectory;
            this.htmlReportDirectory = result.htmlReportDirectory;
            updateMessage("Template setup completed successfully");
            updateProgress(100, 100);
            logger.info("{}: Successfully completed template setup", TASK_NAME);
            return result.success;
        } catch (Exception e) {
            logger.error("{}: Error during template setup", TASK_NAME, e);
            updateMessage("Error: " + e.getMessage());
            return false;
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

    private ProgressCallback javafxProgressCallback() {
        return new ProgressCallback() {
            @Override
            public void onProgress(double workDone, double max) {
                updateProgress(workDone, max);
            }

            @Override
            public void onMessage(String message) {
                if (message != null) {
                    logger.info("{}: {}", TASK_NAME, message);
                    updateMessage(message);
                }
            }

            @Override
            public void onComplete(Object result) {
                if (result != null) {
                    onMessage(String.valueOf(result));
                }
            }

            @Override
            public void onError(Exception error) {
                if (error != null) {
                    logger.error("{}: {}", TASK_NAME, error.getMessage(), error);
                    updateMessage(error.getMessage());
                }
            }
        };
    }
}

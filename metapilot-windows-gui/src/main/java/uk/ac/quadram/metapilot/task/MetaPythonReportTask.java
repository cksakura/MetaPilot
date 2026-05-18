package uk.ac.quadram.metapilot.task;

import javafx.concurrent.Task;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.report.PythonEnvironmentResolver;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * JavaFX adapter for shared Python HTML report generation.
 */
public class MetaPythonReportTask extends Task<Boolean> {

    private static final String TASK_NAME = "Python report generation task";
    private static final SimpleDateFormat FORMAT = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

    public static final String ID_SUMMARY = MetaPilotHtmlReportService.ID_SUMMARY;
    public static final String PEPTIDE = MetaPilotHtmlReportService.PEPTIDE;
    public static final String PROTEIN = MetaPilotHtmlReportService.PROTEIN;
    public static final String TAXONOMY = MetaPilotHtmlReportService.TAXONOMY;
    public static final String FUNCTION = MetaPilotHtmlReportService.FUNCTION;
    public static final String BUBBLE_CHART = MetaPilotHtmlReportService.BUBBLE_CHART;

    private final Logger logger;
    private final int threadCount;
    private final MetaResourceParameters resourceParameters;
    private final MetaPilotHtmlReportService.Request request;
    private String workflowType = "UNKNOWN";
    private File reportsDir;
    private PythonEnvironmentResolver.Result pythonResult;

    @Deprecated
    public MetaPythonReportTask(int threadCount) {
        this(threadCount, new MetaResourceParameters(), LogManager.getLogger(MetaPythonReportTask.class));
    }

    @Deprecated
    public MetaPythonReportTask(int threadCount, MetaResourceParameters resourceParameters) {
        this(threadCount, resourceParameters, LogManager.getLogger(MetaPythonReportTask.class));
    }

    public MetaPythonReportTask(int threadCount, MetaResourceParameters resourceParameters, Logger logger) {
        this.threadCount = Math.max(1, threadCount);
        this.resourceParameters = resourceParameters != null ? resourceParameters : new MetaResourceParameters();
        this.logger = logger != null ? logger : LogManager.getLogger(MetaPythonReportTask.class);
        this.request = new MetaPilotHtmlReportService.Request();
        this.request.threadCount = this.threadCount;
        this.request.resourceParameters = this.resourceParameters;
        this.request.installTemplate = false;
        this.request.warnOnly = false;
        this.request.progressCallback = javafxProgressCallback();
        this.pythonResult = resolvePython();

        this.logger.info("{}: Initialized with {} threads, Python: {}", TASK_NAME, this.threadCount,
                getPythonExecutable());
        System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME
                + ": Initialized with " + this.threadCount + " threads, Python: " + getPythonExecutable());
    }

    public void setWorkflowType(String workflowType) {
        this.workflowType = workflowType != null ? workflowType : "UNKNOWN";
        this.request.workflowType = this.workflowType;
    }

    public void initializeManifestGenerator(File reportsDir) {
        this.reportsDir = reportsDir;
        this.request.reportsDir = reportsDir;
        if (reportsDir != null) {
            this.request.reportDir = reportsDir.getParentFile();
            this.request.resultDir = this.request.reportDir != null ? this.request.reportDir.getParentFile() : null;
        }
    }

    public boolean addTask(String reportType, File inputFile, File outputFile) {
        return addTask(reportType, inputFile, outputFile, null);
    }

    public boolean addTask(String reportType, File inputFile, File outputFile, File metaFile) {
        if (reportType == null || reportType.trim().isEmpty()) {
            logger.error("{}: Report type cannot be null or empty", TASK_NAME);
            return false;
        }
        if (inputFile == null || !inputFile.exists()) {
            logger.error("{}: Input file does not exist: {}", TASK_NAME, inputFile);
            return false;
        }
        if (outputFile == null) {
            logger.error("{}: Output file cannot be null", TASK_NAME);
            return false;
        }
        if (reportsDir == null && outputFile.getParentFile() != null) {
            initializeManifestGenerator(outputFile.getParentFile());
        }
        boolean added = request.addTask(reportType, inputFile, outputFile, metaFile);
        if (added) {
            logger.info("{}: Added {} report task - Input: {}, Output: {}", TASK_NAME, reportType,
                    inputFile.getName(), outputFile.getName());
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Added " + reportType
                    + " report task");
        }
        return added;
    }

    @Override
    protected Boolean call() {
        try {
            updateMessage("Starting Python report generation...");
            updateProgress(0, 100);
            if (request.tasks.isEmpty()) {
                updateMessage("No reports queued");
                updateProgress(100, 100);
                return true;
            }
            if (!isPythonAvailable()) {
                String errorMsg = "Python environment not available: "
                        + (pythonResult != null ? pythonResult.failureMessage : "unknown");
                logger.error("{}: {}", TASK_NAME, errorMsg);
                updateMessage(errorMsg);
                return false;
            }

            MetaPilotHtmlReportService.Result result = MetaPilotHtmlReportService.run(request);
            updateProgress(100, 100);
            updateMessage(result.success ? "All reports generated successfully" : "Report generation failed");
            return result.success;
        } catch (Exception e) {
            logger.error("{}: Error during report generation", TASK_NAME, e);
            updateMessage("Error during report generation: " + e.getMessage());
            return false;
        }
    }

    public String getTaskName() {
        return TASK_NAME;
    }

    public Logger getLogger() {
        return logger;
    }

    public boolean isPythonAvailable() {
        pythonResult = resolvePython();
        return pythonResult != null && pythonResult.valid;
    }

    public String getPythonExecutable() {
        if (pythonResult != null && pythonResult.resolvedCommand != null) {
            return pythonResult.resolvedCommand;
        }
        try {
            String configured = resourceParameters.getToolPath(MetaResourceParameters.ToolType.PYTHON);
            if (configured != null && !configured.isBlank()) {
                return configured;
            }
        } catch (Exception ignored) {
        }
        return "python";
    }

    private PythonEnvironmentResolver.Result resolvePython() {
        PythonEnvironmentResolver.Request pyRequest = new PythonEnvironmentResolver.Request();
        pyRequest.resourceParameters = resourceParameters;
        pyRequest.progressCallback = new ProgressCallback.NoOp();
        return PythonEnvironmentResolver.resolve(pyRequest);
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
                    System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": " + message);
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

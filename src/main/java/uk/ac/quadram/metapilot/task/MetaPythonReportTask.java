package uk.ac.quadram.metapilot.task;

import javafx.concurrent.Task;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;

/**
 * Python-based report generation task that replaces the R-based MetaReportTask.
 * 
 * This task generates HTML reports using Python scripts instead of R, providing
 * the same functionality but with improved performance and reduced dependencies.
 * 
 * <h3>Features:</h3>
 * <ul>
 *   <li>Generates comprehensive HTML reports for metaproteomics analysis results</li>
 *   <li>Uses Python with Plotly for interactive visualizations</li>
 *   <li>Supports multiple report types: ID summary, peptides, proteins, taxonomy, and function</li>
 *   <li>Thread-safe execution with configurable thread pool</li>
 *   <li>Comprehensive error handling and logging</li>
 * </ul>
 * 
 * <h3>Supported Report Types:</h3>
 * <ul>
 *   <li><strong>ID Summary</strong> - Overall identification statistics and sample comparisons</li>
 *   <li><strong>Peptides Summary</strong> - Peptide-level analysis with length, charge, mass distributions</li>
 *   <li><strong>Proteins Summary</strong> - Protein-level analysis with coverage and peptide counts</li>
 *   <li><strong>Taxonomy Summary</strong> - Taxonomic classification results and distributions</li>
 *   <li><strong>Function Summary</strong> - Functional annotation results and pathway analysis</li>
 * </ul>
 * 
 * <h3>Usage Example:</h3>
 * <pre>
 * MetaPythonReportTask reportTask = new MetaPythonReportTask(4); // 4 threads
 * 
 * // Add different report types
 * reportTask.addTask(MetaPythonReportTask.PEPTIDE, 
 *                   new File("final_peptides.tsv"), 
 *                   new File("report_peptides_summary.html"));
 * 
 * reportTask.addTask(MetaPythonReportTask.PROTEIN,
 *                   new File("final_proteins.tsv"),
 *                   new File("report_proteins_summary.html"));
 * 
 * // Execute the task
 * Thread thread = new Thread(reportTask);
 * thread.start();
 * </pre>
 * 
 * <h3>File Format Requirements:</h3>
 * <p>Input files must be tab-separated values (TSV) with specific column headers:</p>
 * <ul>
 *   <li><strong>Peptides file</strong>: Sequence, Length, Charge, Mass, Score, Intensity columns</li>
 *   <li><strong>Proteins file</strong>: Protein IDs, Peptide counts, Score, Intensity columns</li>
 *   <li><strong>Summary file</strong>: Sample, PSMs, Peptides, Proteins columns</li>
 *   <li><strong>Taxa file</strong>: Taxon, Rank, Intensity columns</li>
 *   <li><strong>Function file</strong>: Function, Category, Intensity columns</li>
 * </ul>
 * 
 * @author MetaPilot Development Team
 * @version 1.0
 * @since MetaPilot 3.0
 */
public class MetaPythonReportTask extends Task<Boolean> {

    private final Logger logger;
    private static final String TASK_NAME = "Python report generation task";
    private static final SimpleDateFormat FORMAT = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

    // Report type constants
    public static final String ID_SUMMARY = "id_summary";
    public static final String PEPTIDE = "peptides";
    public static final String PROTEIN = "proteins";
    public static final String TAXONOMY = "taxonomy";
    public static final String FUNCTION = "function";
    public static final String BUBBLE_CHART = "bubble_chart";

    // Python configuration - will be resolved from MetaResourceParameters
    private String pythonExecutable;

    // Task management
    private final ExecutorService executor;
    private final List<ReportTaskInfo> reportTasks;
    private final List<File> outputFiles;
    private final MetaResourceParameters resourceParameters;
    
    // Manifest integration
    private ReportManifestGenerator manifestGenerator;
    private String workflowType = "UNKNOWN";

    /**
     * Creates a new Python report generation task with specified thread count.
     * 
     * @param threadCount Number of threads to use for parallel report generation
     * @deprecated Use constructor with Logger parameter for proper task-specific logging
     */
    @Deprecated
    public MetaPythonReportTask(int threadCount) {
        this(threadCount, new MetaResourceParameters(), LogManager.getLogger(MetaPythonReportTask.class));
    }
    
    /**
     * Creates a new Python report generation task with specified thread count and resource parameters.
     * 
     * @param threadCount Number of threads to use for parallel report generation
     * @param resourceParameters Resource parameters containing Python executable path
     * @deprecated Use constructor with Logger parameter for proper task-specific logging
     */
    @Deprecated
    public MetaPythonReportTask(int threadCount, MetaResourceParameters resourceParameters) {
        this(threadCount, resourceParameters, LogManager.getLogger(MetaPythonReportTask.class));
    }
    
    /**
     * Creates a new Python report generation task with specified thread count, resource parameters, and logger.
     * This is the recommended constructor for proper task-specific logging.
     * 
     * @param threadCount Number of threads to use for parallel report generation
     * @param resourceParameters Resource parameters containing Python executable path
     * @param logger Logger instance from parent task for task-specific logging
     */
    public MetaPythonReportTask(int threadCount, MetaResourceParameters resourceParameters, Logger logger) {
        this.logger = logger;
        this.resourceParameters = resourceParameters;
        this.executor = Executors.newFixedThreadPool(threadCount);
        this.reportTasks = new ArrayList<>();
        this.outputFiles = new ArrayList<>();
        
        // Resolve Python executable path from resource parameters
        this.pythonExecutable = resolvePythonExecutable();
        
        logger.info("{}: Initialized with {} threads, Python: {}", TASK_NAME, threadCount, pythonExecutable);
        logger.info("{}: Python availability check: {}", TASK_NAME, isPythonAvailable());
        System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Initialized with " + threadCount + " threads, Python: " + pythonExecutable);
        System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Python availability: " + isPythonAvailable());
    }

    /**
     * Sets the workflow type for manifest generation.
     * 
     * @param workflowType The workflow type (e.g., "MAG", "DIA", "IGC")
     */
    public void setWorkflowType(String workflowType) {
        this.workflowType = workflowType;
    }
    
    /**
     * Initializes the manifest generator for the given reports directory.
     * 
     * @param reportsDir The reports directory where manifest files will be created
     */
    public void initializeManifestGenerator(File reportsDir) {
        this.manifestGenerator = new ReportManifestGenerator(reportsDir);
    }

    /**
     * Adds a report generation task to the queue.
     * 
     * @param reportType Type of report to generate (use constants: ID_SUMMARY, PEPTIDE, etc.)
     * @param inputFile Input data file (TSV format)
     * @param outputFile Output HTML report file
     * @return true if task was added successfully, false otherwise
     */
    public boolean addTask(String reportType, File inputFile, File outputFile) {
        return addTask(reportType, inputFile, outputFile, null);
    }

    /**
     * Adds a report generation task to the queue with optional metadata file.
     * 
     * @param reportType Type of report to generate
     * @param inputFile Input data file (TSV format)
     * @param outputFile Output HTML report file
     * @param metaFile Optional metadata file for additional context
     * @return true if task was added successfully, false otherwise
     */
    public boolean addTask(String reportType, File inputFile, File outputFile, File metaFile) {
        // Validate input parameters
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

        // Create output directory if it doesn't exist
        File outputDir = outputFile.getParentFile();
        if (outputDir != null && !outputDir.exists()) {
            outputDir.mkdirs();
        }

        // Create task info
        ReportTaskInfo taskInfo = new ReportTaskInfo(reportType, inputFile, outputFile, metaFile);
        reportTasks.add(taskInfo);
        outputFiles.add(outputFile);

        logger.info("{}: Added {} report task - Input: {}, Output: {}", 
                   TASK_NAME, reportType, inputFile.getName(), outputFile.getName());
        System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Added " + reportType + " report task");

        // Submit task to executor
        executor.submit(() -> generateReport(taskInfo));

        return true;
    }

    /**
     * Main task execution method called by JavaFX Task framework.
     * 
     * @return true if all reports were generated successfully, false otherwise
     */
    @Override
    protected Boolean call() throws Exception {
        try {
            updateMessage("Starting Python report generation...");
            updateProgress(0, 100);

            // Validate Python environment
            if (!isPythonAvailable()) {
                String errorMsg = "Python environment not available. Please check resource configuration.";
                logger.error("{}: {}", TASK_NAME, errorMsg);
                updateMessage(errorMsg);
                return false;
            }

            logger.info("{}: Starting report generation for {} tasks using Python: {}", TASK_NAME, reportTasks.size(), pythonExecutable);
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Starting report generation using Python: " + pythonExecutable);

            // Shutdown executor and wait for completion
            executor.shutdown();
            boolean finished = executor.awaitTermination(30, TimeUnit.MINUTES);

            if (!finished) {
                logger.error("{}: Report generation timed out after 30 minutes", TASK_NAME);
                updateMessage("Report generation timed out");
                return false;
            }

            // Verify all output files were created
            updateMessage("Verifying generated reports...");
            updateProgress(90, 100);

            for (File outputFile : outputFiles) {
                if (!outputFile.exists() || outputFile.length() == 0) {
                    logger.error("{}: Report file not generated or empty: {}", TASK_NAME, outputFile.getName());
                    updateMessage("Failed to generate report: " + outputFile.getName());
                    return false;
                }
            }

            updateMessage("All reports generated successfully");
            updateProgress(100, 100);

            logger.info("{}: All reports generated successfully", TASK_NAME);
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": All reports generated successfully");

            return true;

        } catch (Exception e) {
            logger.error("{}: Error during report generation", TASK_NAME, e);
            updateMessage("Error during report generation: " + e.getMessage());
            return false;
        }
    }

    /**
     * Generates a single report using Python script.
     * 
     * @param taskInfo Information about the report to generate
     */
    private void generateReport(ReportTaskInfo taskInfo) {
        try {
            logger.info("{}: Generating {} report from {}", TASK_NAME, taskInfo.reportType, taskInfo.inputFile.getName());
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Generating " + taskInfo.reportType + " report");

            // Build Python command
            List<String> command = buildPythonCommand(taskInfo);

            // Execute Python script
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.redirectErrorStream(true);
            processBuilder.directory(new File(System.getProperty("user.dir"))); // Set working directory

            logger.info("{}: Executing Python command: {}", TASK_NAME, String.join(" ", command));
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Executing Python command: " + String.join(" ", command));

            Process process = processBuilder.start();

            // Read output from Python script
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    logger.debug("{}: Python output: {}", TASK_NAME, line);
                    System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": " + line);
                }
            }

            // Wait for process completion
            int exitCode = process.waitFor();

            if (exitCode == 0) {
                logger.info("{}: Successfully generated {} report", TASK_NAME, taskInfo.reportType);
                System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Successfully generated " + taskInfo.reportType + " report");
                
                // Add to manifest if generator is available
                if (manifestGenerator != null && taskInfo.outputFile.exists()) {
                    try {
                        String fileName = taskInfo.outputFile.getName();
                        manifestGenerator.addPythonReport(taskInfo.reportType, fileName, workflowType);
                        logger.info("{}: Added {} report to manifest", TASK_NAME, taskInfo.reportType);
                    } catch (Exception e) {
                        logger.warn("{}: Failed to add {} report to manifest: {}", TASK_NAME, taskInfo.reportType, e.getMessage());
                    }
                }
            } else {
                logger.error("{}: Python script failed with exit code {} for {} report", TASK_NAME, exitCode, taskInfo.reportType);
                System.err.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Python script failed for " + taskInfo.reportType + " report");
            }

        } catch (Exception e) {
            logger.error("{}: Error generating {} report", TASK_NAME, taskInfo.reportType, e);
            System.err.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Error generating " + taskInfo.reportType + " report: " + e.getMessage());
        }
    }

    /**
     * Escapes a path for use in Python string literals.
     * 
     * @param path The path to escape
     * @return Escaped path safe for Python
     */
    private String escapePathForPython(String path) {
        if (path == null) return "";
        // Replace backslashes with forward slashes for cross-platform compatibility
        // and escape any single quotes that might be in the path
        return path.replace("\\", "/").replace("'", "\\'");
    }

    /**
     * Builds the Python command for report generation.
     * 
     * @param taskInfo Task information
     * @return List of command arguments
     */
    private List<String> buildPythonCommand(ReportTaskInfo taskInfo) {
        List<String> command = new ArrayList<>();
        
        // Load data and metadata - properly escape paths for Python
        String inputPath = escapePathForPython(taskInfo.inputFile.getAbsolutePath());
        String outputPath = escapePathForPython(taskInfo.outputFile.getAbsolutePath());
        String outputDir = escapePathForPython(taskInfo.outputFile.getParent());
        
        logger.info("{}: Processing paths - Input: {}, Output: {}, OutputDir: {}", TASK_NAME, inputPath, outputPath, outputDir);
        
        // Create temporary Python script file to avoid command line length issues
        try {
            File tempScript = File.createTempFile("metapilot_report_", ".py");
            tempScript.deleteOnExit();
            
            // Python script to execute
            StringBuilder pythonScript = new StringBuilder();
            pythonScript.append("import sys, os\n");
            pythonScript.append("sys.path.append(os.path.join(os.getcwd(), 'src', 'main', 'python'))\n");
            pythonScript.append("from report_generation.metapilot_report_generator import MetaPilotReportGenerator\n");
            pythonScript.append("import pandas as pd\n");
            pythonScript.append("import numpy as np\n");
            pythonScript.append("rg = MetaPilotReportGenerator('").append(outputDir).append("')\n");
            
            // Load metadata if available
            if (taskInfo.metaFile != null && taskInfo.metaFile.exists()) {
                String metaPath = escapePathForPython(taskInfo.metaFile.getAbsolutePath());
                pythonScript.append("metadata = rg.load_metadata('").append(metaPath).append("')\n");
            } else {
                pythonScript.append("metadata = None\n");
            }
        
            switch (taskInfo.reportType) {
                case ID_SUMMARY:
                    pythonScript.append("df = pd.read_csv('").append(inputPath).append("', sep='\\t')\n");
                    pythonScript.append("fig = rg.create_id_summary(df, metadata)\n");
                    pythonScript.append("rg.generate_html_report(fig, 'ID Summary Report', '").append(outputPath).append("')\n");
                    break;
                    
                case PEPTIDE:
                    pythonScript.append("df = pd.read_csv('").append(inputPath).append("', sep='\\t')\n");
                    pythonScript.append("fig = rg.create_peptides_summary(df, metadata)\n");
                    pythonScript.append("rg.generate_html_report(fig, 'Peptides Summary Report', '").append(outputPath).append("')\n");
                    break;
                    
                case PROTEIN:
                    pythonScript.append("df = pd.read_csv('").append(inputPath).append("', sep='\\t')\n");
                    pythonScript.append("fig = rg.create_proteins_summary(df, metadata)\n");
                    pythonScript.append("rg.generate_html_report(fig, 'Proteins Summary Report', '").append(outputPath).append("')\n");
                    break;
                    
                case TAXONOMY:
                    pythonScript.append("df = pd.read_csv('").append(inputPath).append("', sep='\\t')\n");
                    pythonScript.append("fig = rg.create_taxonomy_summary(df, metadata)\n");
                    pythonScript.append("rg.generate_html_report(fig, 'Taxonomy Summary Report', '").append(outputPath).append("')\n");
                    break;
                    
                case FUNCTION:
                    // For function reports, inputFile is already the folder path containing all separate files
                    // No need to load main file - all data comes from separate files
                    pythonScript.append("fig = rg.create_function_summary(metadata, input_file_path='").append(inputPath).append("')\n");
                    pythonScript.append("rg.generate_html_report(fig, 'Function Summary Report', '").append(outputPath).append("')\n");
                    break;
                    
                case BUBBLE_CHART:
                    // TEMPORARILY DISABLED: Bubble chart report generation is currently hidden due to issues
                    pythonScript.append("# Bubble chart generation temporarily disabled\n");
                    pythonScript.append("import os\n");
                    pythonScript.append("with open('").append(outputPath).append("', 'w') as f:\n");
                    pythonScript.append("    f.write('<html><body><h1>Bubble Chart Report</h1><p>Bubble chart generation is temporarily disabled.</p></body></html>')\n");
                    break;
                    
                default:
                    throw new IllegalArgumentException("Unsupported report type: " + taskInfo.reportType);
            }
            
            // Write script to temporary file
            try (FileWriter writer = new FileWriter(tempScript)) {
                writer.write(pythonScript.toString());
            }
            
            // Build command to execute the script file
            command.add(pythonExecutable);
            command.add(tempScript.getAbsolutePath());
            
        } catch (IOException e) {
            logger.error("{}: Failed to create temporary Python script: {}", TASK_NAME, e.getMessage());
            throw new RuntimeException("Failed to create temporary Python script", e);
        }
        
        return command;
    }

    /**
     * Internal class to hold report task information.
     */
    private static class ReportTaskInfo {
        final String reportType;
        final File inputFile;
        final File outputFile;
        final File metaFile;

        ReportTaskInfo(String reportType, File inputFile, File outputFile, File metaFile) {
            this.reportType = reportType;
            this.inputFile = inputFile;
            this.outputFile = outputFile;
            this.metaFile = metaFile;
        }
    }

    /**
     * Gets the task name for logging purposes.
     * 
     * @return Task name
     */
    public String getTaskName() {
        return TASK_NAME;
    }

    /**
     * Gets the logger instance for this task.
     * 
     * @return Logger instance
     */
    public Logger getLogger() {
        return logger;
    }
    
    /**
     * Resolves the Python executable path from resource parameters.
     * 
     * @return Full path to Python executable, or "python" as fallback
     */
    private String resolvePythonExecutable() {
        try {
            // Get Python path from resource parameters
            String pythonPath = resourceParameters.getToolPath(MetaResourceParameters.ToolType.PYTHON);
            logger.info("{}: Retrieved Python path from MetaResourceParameters: '{}'", TASK_NAME, pythonPath);
            
            if (pythonPath != null && !pythonPath.trim().isEmpty()) {
                File pythonFile = new File(pythonPath);
                logger.info("{}: Checking Python file existence: {} -> {}", TASK_NAME, pythonPath, pythonFile.exists());
                if (pythonFile.exists()) {
                    logger.info("{}: Using bundled Python executable: {}", TASK_NAME, pythonPath);
                    return pythonPath;
                } else {
                    logger.warn("{}: Bundled Python executable not found at: {}, falling back to system python", TASK_NAME, pythonPath);
                }
            } else {
                logger.warn("{}: No Python path configured in resource parameters, falling back to system python", TASK_NAME);
            }
        } catch (Exception e) {
            logger.warn("{}: Error resolving Python executable from resource parameters: {}, falling back to system python", TASK_NAME, e.getMessage());
        }
        
        // Fallback to system python
        logger.info("{}: Using system Python executable as fallback", TASK_NAME);
        return "python";
    }
    
    /**
     * Checks if the Python environment is properly configured.
     * 
     * @return true if Python is available, false otherwise
     */
    public boolean isPythonAvailable() {
        return resourceParameters.isToolAvailable(MetaResourceParameters.ToolType.PYTHON);
    }
    
    /**
     * Gets the Python executable path being used.
     * 
     * @return Python executable path
     */
    public String getPythonExecutable() {
        return pythonExecutable;
    }
} 
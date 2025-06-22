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

    private static final Logger LOGGER = LogManager.getLogger(MetaPythonReportTask.class);
    private static final String TASK_NAME = "Python report generation task";
    private static final SimpleDateFormat FORMAT = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

    // Report type constants
    public static final String ID_SUMMARY = "id_summary";
    public static final String PEPTIDE = "peptides";
    public static final String PROTEIN = "proteins";
    public static final String TAXONOMY = "taxonomy";
    public static final String FUNCTION = "function";

    // Python configuration
    private static final String PYTHON_EXECUTABLE = "python";

    // Task management
    private final ExecutorService executor;
    private final List<ReportTaskInfo> reportTasks;
    private final List<File> outputFiles;
    private final int threadCount;

    /**
     * Creates a new Python report generation task with specified thread count.
     * 
     * @param threadCount Number of threads to use for parallel report generation
     */
    public MetaPythonReportTask(int threadCount) {
        this.threadCount = threadCount;
        this.executor = Executors.newFixedThreadPool(threadCount);
        this.reportTasks = new ArrayList<>();
        this.outputFiles = new ArrayList<>();
        
        LOGGER.info("{}: Initialized with {} threads", TASK_NAME, threadCount);
        System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Initialized with " + threadCount + " threads");
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
            LOGGER.error("{}: Report type cannot be null or empty", TASK_NAME);
            return false;
        }

        if (inputFile == null || !inputFile.exists()) {
            LOGGER.error("{}: Input file does not exist: {}", TASK_NAME, inputFile);
            return false;
        }

        if (outputFile == null) {
            LOGGER.error("{}: Output file cannot be null", TASK_NAME);
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

        LOGGER.info("{}: Added {} report task - Input: {}, Output: {}", 
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

            LOGGER.info("{}: Starting report generation for {} tasks", TASK_NAME, reportTasks.size());
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Starting report generation");

            // Shutdown executor and wait for completion
            executor.shutdown();
            boolean finished = executor.awaitTermination(30, TimeUnit.MINUTES);

            if (!finished) {
                LOGGER.error("{}: Report generation timed out after 30 minutes", TASK_NAME);
                updateMessage("Report generation timed out");
                return false;
            }

            // Verify all output files were created
            updateMessage("Verifying generated reports...");
            updateProgress(90, 100);

            for (File outputFile : outputFiles) {
                if (!outputFile.exists() || outputFile.length() == 0) {
                    LOGGER.error("{}: Report file not generated or empty: {}", TASK_NAME, outputFile.getName());
                    updateMessage("Failed to generate report: " + outputFile.getName());
                    return false;
                }
            }

            updateMessage("All reports generated successfully");
            updateProgress(100, 100);

            LOGGER.info("{}: All reports generated successfully", TASK_NAME);
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": All reports generated successfully");

            return true;

        } catch (Exception e) {
            LOGGER.error("{}: Error during report generation", TASK_NAME, e);
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
            LOGGER.info("{}: Generating {} report from {}", TASK_NAME, taskInfo.reportType, taskInfo.inputFile.getName());
            System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Generating " + taskInfo.reportType + " report");

            // Build Python command
            List<String> command = buildPythonCommand(taskInfo);

            // Execute Python script
            ProcessBuilder processBuilder = new ProcessBuilder(command);
            processBuilder.redirectErrorStream(true);

            Process process = processBuilder.start();

            // Read output from Python script
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    LOGGER.debug("{}: Python output: {}", TASK_NAME, line);
                    System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": " + line);
                }
            }

            // Wait for process completion
            int exitCode = process.waitFor();

            if (exitCode == 0) {
                LOGGER.info("{}: Successfully generated {} report", TASK_NAME, taskInfo.reportType);
                System.out.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Successfully generated " + taskInfo.reportType + " report");
            } else {
                LOGGER.error("{}: Python script failed with exit code {} for {} report", TASK_NAME, exitCode, taskInfo.reportType);
                System.err.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Python script failed for " + taskInfo.reportType + " report");
            }

        } catch (Exception e) {
            LOGGER.error("{}: Error generating {} report", TASK_NAME, taskInfo.reportType, e);
            System.err.println(FORMAT.format(new Date()) + "\t" + TASK_NAME + ": Error generating " + taskInfo.reportType + " report: " + e.getMessage());
        }
    }

    /**
     * Builds the Python command for report generation.
     * 
     * @param taskInfo Task information
     * @return List of command arguments
     */
    private List<String> buildPythonCommand(ReportTaskInfo taskInfo) {
        List<String> command = new ArrayList<>();
        
        command.add(PYTHON_EXECUTABLE);
        command.add("-c");
        
        // Python script to execute
        StringBuilder pythonScript = new StringBuilder();
        pythonScript.append("import sys; sys.path.append('src/main/python'); ");
        pythonScript.append("from report_generation.metapilot_report_generator import MetaPilotReportGenerator; ");
        pythonScript.append("import pandas as pd; import numpy as np; ");
        pythonScript.append("rg = MetaPilotReportGenerator('").append(taskInfo.outputFile.getParent()).append("'); ");
        
        // Load data and metadata
        String inputPath = taskInfo.inputFile.getAbsolutePath().replace("\\", "/");
        String outputPath = taskInfo.outputFile.getAbsolutePath().replace("\\", "/");
        
        pythonScript.append("df = pd.read_csv('").append(inputPath).append("', sep='\\t'); ");
        
        // Load metadata if available
        if (taskInfo.metaFile != null && taskInfo.metaFile.exists()) {
            String metaPath = taskInfo.metaFile.getAbsolutePath().replace("\\", "/");
            pythonScript.append("metadata = rg.load_metadata('").append(metaPath).append("'); ");
        } else {
            pythonScript.append("metadata = None; ");
        }
        
        switch (taskInfo.reportType) {
            case ID_SUMMARY:
                pythonScript.append("fig = rg.create_id_summary(df, metadata); ");
                pythonScript.append("rg.generate_html_report(fig, 'ID Summary Report', '").append(outputPath).append("')");
                break;
                
            case PEPTIDE:
                pythonScript.append("fig = rg.create_peptides_summary(df, metadata); ");
                pythonScript.append("rg.generate_html_report(fig, 'Peptides Summary Report', '").append(outputPath).append("')");
                break;
                
            case PROTEIN:
                pythonScript.append("fig = rg.create_proteins_summary(df, metadata); ");
                pythonScript.append("rg.generate_html_report(fig, 'Proteins Summary Report', '").append(outputPath).append("')");
                break;
                
            case TAXONOMY:
                pythonScript.append("fig = rg.create_taxonomy_summary(df, metadata); ");
                pythonScript.append("rg.generate_html_report(fig, 'Taxonomy Summary Report', '").append(outputPath).append("')");
                break;
                
            case FUNCTION:
                pythonScript.append("fig = rg.create_function_summary(df, metadata); ");
                pythonScript.append("rg.generate_html_report(fig, 'Function Summary Report', '").append(outputPath).append("')");
                break;
                
            default:
                throw new IllegalArgumentException("Unsupported report type: " + taskInfo.reportType);
        }
        
        command.add(pythonScript.toString());
        
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
        return LOGGER;
    }
} 
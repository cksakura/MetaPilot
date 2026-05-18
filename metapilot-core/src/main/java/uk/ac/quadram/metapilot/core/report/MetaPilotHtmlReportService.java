package uk.ac.quadram.metapilot.core.report;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.io.ReportManifestGenerator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Shared HTML report generator for GUI and headless workflows.
 */
public final class MetaPilotHtmlReportService {

    private static final Logger LOGGER = LogManager.getLogger(MetaPilotHtmlReportService.class);
    private static final String PYTHON_SCRIPTS_RESOURCE_PATH = "/python_scripts";

    public static final String ID_SUMMARY = "id_summary";
    public static final String PEPTIDE = "peptides";
    public static final String PROTEIN = "proteins";
    public static final String TAXONOMY = "taxonomy";
    public static final String FUNCTION = "function";
    public static final String BUBBLE_CHART = "bubble_chart";

    private MetaPilotHtmlReportService() {
    }

    public static final class ReportTaskSpec {
        public final String reportType;
        public final File inputFile;
        public final File outputFile;
        public final File metadataFile;

        public ReportTaskSpec(String reportType, File inputFile, File outputFile, File metadataFile) {
            this.reportType = reportType;
            this.inputFile = inputFile;
            this.outputFile = outputFile;
            this.metadataFile = metadataFile;
        }
    }

    public static final class Request {
        public File resultDir;
        public File reportDir;
        public File reportsDir;
        public int threadCount = 1;
        public String workflowType = "UNKNOWN";
        public String pythonHint;
        public MetaResourceParameters resourceParameters;
        public MetaPilotResourceConfig resourceConfig;
        public File metadataFile;
        public boolean installTemplate = true;
        public boolean warnOnly;
        public ProgressCallback progressCallback;
        public final List<ReportTaskSpec> tasks = new ArrayList<>();

        public boolean addTask(String reportType, File inputFile, File outputFile) {
            return addTask(reportType, inputFile, outputFile, metadataFile);
        }

        public boolean addTask(String reportType, File inputFile, File outputFile, File metaFile) {
            if (reportType == null || reportType.trim().isEmpty()) {
                return false;
            }
            if (inputFile == null || !inputFile.exists()) {
                return false;
            }
            if (outputFile == null) {
                return false;
            }
            File parent = outputFile.getParentFile();
            if (parent != null) {
                parent.mkdirs();
            }
            tasks.add(new ReportTaskSpec(reportType, inputFile, outputFile, metaFile));
            return true;
        }
    }

    public static final class Result {
        public File reportDir;
        public File reportsDir;
        public File manifestFile;
        public PythonEnvironmentResolver.Result python;
        public boolean success;
        public final List<File> generatedHtmlFiles = new ArrayList<>();
        public final List<String> warnings = new ArrayList<>();
    }

    public static Request buildStandardRequest(
            File resultDir,
            String workflowType,
            int threadCount,
            String pythonHint,
            File metadataFile,
            ProgressCallback cb
    ) {
        Request request = new Request();
        request.resultDir = resultDir;
        request.reportDir = new File(resultDir, "report");
        request.reportsDir = new File(request.reportDir, "reports");
        request.workflowType = workflowType != null ? workflowType : "UNKNOWN";
        request.threadCount = threadCount;
        request.pythonHint = pythonHint;
        request.metadataFile = metadataFile;
        request.progressCallback = cb;
        addStandardTasks(request);
        return request;
    }

    public static void addStandardTasks(Request request) {
        if (request == null || request.resultDir == null) {
            return;
        }
        if (request.reportDir == null) {
            request.reportDir = new File(request.resultDir, "report");
        }
        if (request.reportsDir == null) {
            request.reportsDir = new File(request.reportDir, "reports");
        }

        File resultDir = request.resultDir;
        File reportsDir = request.reportsDir;
        File metadata = request.metadataFile;
        addIfPresent(request, ID_SUMMARY, new File(resultDir, "final_summary.tsv"),
                new File(reportsDir, "report_ID_summary.html"), metadata);
        addIfPresent(request, PEPTIDE, new File(resultDir, "final_peptides.tsv"),
                new File(reportsDir, "report_peptides_summary.html"), metadata);
        addIfPresent(request, PROTEIN, new File(resultDir, "razor_proteins.tsv"),
                new File(reportsDir, "report_proteinGroups_summary.html"), metadata);
        addIfPresent(request, TAXONOMY, new File(resultDir, "taxonomy_analysis/Taxa_report.tsv"),
                new File(reportsDir, "report_taxonomy_summary.html"), metadata);

        File functionDir = new File(resultDir, "functional_annotation");
        if (functionDir.isDirectory() && functionReportInputsExist(functionDir)) {
            request.addTask(FUNCTION, functionDir, new File(reportsDir, "report_function_summary.html"), metadata);
        }
    }

    public static Result run(Request request) throws IOException {
        if (request == null) {
            throw new IllegalArgumentException("request is null");
        }
        ProgressCallback cb = request.progressCallback != null ? request.progressCallback : new ProgressCallback.NoOp();
        Result result = new Result();

        try {
            normalizeDirectories(request);
            result.reportDir = request.reportDir;
            result.reportsDir = request.reportsDir;
            result.manifestFile = new File(request.reportsDir, "manifest.json");

            if (isReportDisabled(request.resourceConfig)) {
                result.success = true;
                result.warnings.add("HTML report generation disabled by GENERATE_REPORT=false");
                cb.onMessage("HTML report generation disabled by GENERATE_REPORT=false");
                return result;
            }

            if (request.installTemplate) {
                ReportTemplateInstaller.Result installed = ReportTemplateInstaller.install(request.resultDir, cb);
                result.warnings.addAll(installed.warnings);
                request.reportDir = installed.reportDirectory;
                request.reportsDir = installed.htmlReportDirectory;
                result.reportDir = request.reportDir;
                result.reportsDir = request.reportsDir;
                result.manifestFile = new File(request.reportsDir, "manifest.json");
            } else {
                request.reportDir.mkdirs();
                request.reportsDir.mkdirs();
            }

            ReportManifestGenerator manifest = new ReportManifestGenerator(request.reportsDir);
            manifest.createInitialManifest();

            if (request.tasks.isEmpty()) {
                result.success = true;
                result.warnings.add("No report input TSV files were found");
                return result;
            }

            PythonEnvironmentResolver.Request pyRequest = new PythonEnvironmentResolver.Request();
            pyRequest.explicitPython = request.pythonHint;
            pyRequest.resourceParameters = request.resourceParameters;
            pyRequest.resourceConfig = request.resourceConfig;
            pyRequest.progressCallback = cb;
            result.python = PythonEnvironmentResolver.resolve(pyRequest);
            if (result.python == null || !result.python.valid) {
                String message = "Python report environment is not available: "
                        + (result.python != null ? result.python.failureMessage : "unknown");
                if (request.warnOnly) {
                    result.warnings.add(message);
                    cb.onMessage("WARNING: " + message);
                    return result;
                }
                throw new IOException(message);
            }
            request.pythonHint = result.python.resolvedCommand;

            cb.onMessage("HTML report generation: using Python " + result.python.resolvedCommand);
            File scriptsDir = extractPythonScripts();
            try {
                runTasks(request, result, manifest, scriptsDir, cb);
            } finally {
                FileUtils.deleteQuietly(scriptsDir);
            }

            result.success = result.warnings.isEmpty() || request.warnOnly;
            return result;
        } catch (Exception e) {
            String message = "HTML report generation failed: " + e.getMessage();
            if (request.warnOnly) {
                result.warnings.add(message);
                cb.onMessage("WARNING: " + message);
                result.success = false;
                return result;
            }
            if (e instanceof IOException) {
                throw (IOException) e;
            }
            throw new IOException(message, e);
        }
    }

    public static Result runStandardReportWarningOnly(
            File resultDir,
            String workflowType,
            int threadCount,
            String pythonHint,
            MetaPilotResourceConfig resourceConfig,
            ProgressCallback cb
    ) {
        try {
            Request request = buildStandardRequest(resultDir, workflowType, threadCount, pythonHint, null, cb);
            request.resourceConfig = resourceConfig;
            request.warnOnly = true;
            request.installTemplate = true;
            return run(request);
        } catch (Exception e) {
            Result result = new Result();
            result.reportDir = resultDir != null ? new File(resultDir, "report") : null;
            result.reportsDir = result.reportDir != null ? new File(result.reportDir, "reports") : null;
            result.success = false;
            result.warnings.add("HTML report generation failed: " + e.getMessage());
            if (cb != null) {
                cb.onMessage("WARNING: HTML report generation failed: " + e.getMessage());
            }
            return result;
        }
    }

    private static void addIfPresent(Request request, String reportType, File input, File output, File metadata) {
        if (input.isFile() && input.length() > 0) {
            request.addTask(reportType, input, output, metadata);
        }
    }

    private static boolean functionReportInputsExist(File functionDir) {
        String[] names = {
                "Annotated_proteins.tsv",
                "Gene_Ontology_id_report.tsv",
                "EC_id_report.tsv",
                "KEGG_ko_report.tsv",
                "COG_accession_report.tsv",
                "NOG_accession_report.tsv"
        };
        for (String name : names) {
            File f = new File(functionDir, name);
            if (f.isFile() && f.length() > 0) {
                return true;
            }
        }
        return false;
    }

    private static void normalizeDirectories(Request request) throws IOException {
        if (request.resultDir == null) {
            if (request.reportDir != null) {
                File parent = request.reportDir.getParentFile();
                request.resultDir = parent != null ? parent : request.reportDir;
            } else {
                throw new IOException("resultDir is null");
            }
        }
        if (request.reportDir == null) {
            request.reportDir = new File(request.resultDir, "report");
        }
        if (request.reportsDir == null) {
            request.reportsDir = new File(request.reportDir, "reports");
        }
        request.threadCount = Math.max(1, request.threadCount);
    }

    private static boolean isReportDisabled(MetaPilotResourceConfig config) {
        if (config == null) {
            return false;
        }
        return config.getRaw("GENERATE_REPORT")
                .map(v -> {
                    String t = v.trim().toLowerCase(Locale.ROOT);
                    return t.equals("false") || t.equals("0") || t.equals("no") || t.equals("n");
                })
                .orElse(false);
    }

    private static void runTasks(Request request, Result result, ReportManifestGenerator manifest,
                                 File scriptsDir, ProgressCallback cb) throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(Math.max(1, request.threadCount));
        List<Future<TaskResult>> futures = new ArrayList<>();
        try {
            for (ReportTaskSpec task : request.tasks) {
                futures.add(executor.submit(new ReportCallable(task, request, scriptsDir)));
            }
            for (Future<TaskResult> future : futures) {
                TaskResult taskResult = future.get();
                if (taskResult.success) {
                    result.generatedHtmlFiles.add(taskResult.outputFile);
                    manifest.addPythonReport(manifestType(taskResult.reportType),
                            taskResult.outputFile.getName(), request.workflowType);
                    cb.onMessage("HTML report generated: " + taskResult.outputFile.getName());
                } else {
                    result.warnings.add(taskResult.message);
                    cb.onMessage("WARNING: " + taskResult.message);
                }
            }
        } finally {
            executor.shutdownNow();
        }
    }

    private static String manifestType(String reportType) {
        if (ID_SUMMARY.equals(reportType)) return "ID_SUMMARY";
        if (PEPTIDE.equals(reportType)) return "PEPTIDE";
        if (PROTEIN.equals(reportType)) return "PROTEIN";
        if (TAXONOMY.equals(reportType)) return "TAXONOMY";
        if (FUNCTION.equals(reportType)) return "FUNCTION";
        if (BUBBLE_CHART.equals(reportType)) return "BUBBLE_CHART";
        return reportType == null ? "UNKNOWN" : reportType.toUpperCase(Locale.ROOT);
    }

    private static File extractPythonScripts() throws IOException {
        File tempDir = new File(System.getProperty("java.io.tmpdir"),
                "metapilot_python_" + System.currentTimeMillis());
        if (!tempDir.exists() && !tempDir.mkdirs()) {
            throw new IOException("Could not create temporary Python script directory: " + tempDir.getAbsolutePath());
        }

        URL resourceUrl = MetaPilotHtmlReportService.class.getResource(PYTHON_SCRIPTS_RESOURCE_PATH);
        if (resourceUrl == null) {
            throw new IOException("Python scripts resource not found at: " + PYTHON_SCRIPTS_RESOURCE_PATH);
        }

        String protocol = resourceUrl.getProtocol();
        if ("jar".equals(protocol)) {
            extractPythonScriptsFromJar(resourceUrl, tempDir);
        } else if ("file".equals(protocol)) {
            extractPythonScriptsFromFileSystem(resourceUrl, tempDir);
        } else {
            throw new IOException("Unsupported protocol for Python scripts: " + protocol);
        }
        return tempDir;
    }

    private static void extractPythonScriptsFromJar(URL resourceUrl, File targetDir) throws IOException {
        String jarPath = resourceUrl.getPath();
        int exclamationIdx = jarPath.indexOf("!");
        if (exclamationIdx > 0) {
            jarPath = jarPath.substring(0, exclamationIdx);
            if (jarPath.startsWith("file:")) {
                jarPath = jarPath.substring(5);
            }
        }

        try (JarFile jarFile = new JarFile(jarPath)) {
            Enumeration<JarEntry> entries = jarFile.entries();
            while (entries.hasMoreElements()) {
                JarEntry entry = entries.nextElement();
                String entryName = entry.getName();
                if (!entryName.startsWith("python_scripts/")) {
                    continue;
                }
                String relativePath = entryName.substring("python_scripts/".length());
                if (relativePath.isEmpty() || relativePath.contains("__pycache__")) {
                    continue;
                }
                File targetFile = new File(targetDir, relativePath);
                if (entry.isDirectory()) {
                    targetFile.mkdirs();
                    continue;
                }
                File parent = targetFile.getParentFile();
                if (parent != null) {
                    parent.mkdirs();
                }
                try (InputStream is = jarFile.getInputStream(entry);
                     FileOutputStream fos = new FileOutputStream(targetFile)) {
                    IOUtils.copy(is, fos);
                }
            }
        }
    }

    private static void extractPythonScriptsFromFileSystem(URL resourceUrl, File targetDir) throws IOException {
        try {
            File sourceDir = new File(resourceUrl.toURI());
            if (!sourceDir.exists()) {
                throw new IOException("Source Python scripts not found: " + sourceDir.getAbsolutePath());
            }
            FileUtils.copyDirectory(sourceDir, targetDir);
        } catch (URISyntaxException e) {
            throw new IOException("Invalid Python scripts resource URL", e);
        }
    }

    private static String escapePathForPython(String path) {
        if (path == null) {
            return "";
        }
        return path.replace("\\", "/").replace("'", "\\'");
    }

    private static List<String> buildPythonCommand(ReportTaskSpec task, Request request, File scriptsDir)
            throws IOException {
        String inputPath = escapePathForPython(task.inputFile.getAbsolutePath());
        String outputPath = escapePathForPython(task.outputFile.getAbsolutePath());
        String outputDir = escapePathForPython(task.outputFile.getParent());
        String scriptsPath = escapePathForPython(scriptsDir.getAbsolutePath());

        StringBuilder pythonScript = new StringBuilder();
        pythonScript.append("import sys, os\n");
        pythonScript.append("sys.path.insert(0, r'").append(scriptsPath).append("')\n");
        pythonScript.append("from report_generation.metapilot_report_generator import MetaPilotReportGenerator\n");
        pythonScript.append("import pandas as pd\n");
        pythonScript.append("import numpy as np\n");
        pythonScript.append("rg = MetaPilotReportGenerator('").append(outputDir).append("')\n");

        if (task.metadataFile != null && task.metadataFile.exists()) {
            String metaPath = escapePathForPython(task.metadataFile.getAbsolutePath());
            pythonScript.append("metadata = rg.load_metadata('").append(metaPath).append("')\n");
        } else {
            pythonScript.append("metadata = None\n");
        }

        switch (task.reportType) {
            case ID_SUMMARY:
                appendTableReport(pythonScript, inputPath, outputPath,
                        "create_id_summary", "ID Summary Report");
                break;
            case PEPTIDE:
                appendTableReport(pythonScript, inputPath, outputPath,
                        "create_peptides_summary", "Peptides Summary Report");
                break;
            case PROTEIN:
                appendTableReport(pythonScript, inputPath, outputPath,
                        "create_proteins_summary", "Proteins Summary Report");
                break;
            case TAXONOMY:
                appendTableReport(pythonScript, inputPath, outputPath,
                        "create_taxonomy_summary", "Taxonomy Summary Report");
                break;
            case FUNCTION:
                pythonScript.append("fig = rg.create_function_summary(metadata, input_file_path='")
                        .append(inputPath).append("')\n");
                pythonScript.append("rg.generate_html_report(fig, 'Function Summary Report', '")
                        .append(outputPath).append("')\n");
                break;
            case BUBBLE_CHART:
                pythonScript.append("with open('").append(outputPath).append("', 'w') as f:\n");
                pythonScript.append("    f.write('<html><body><h1>Bubble Chart Report</h1>"
                        + "<p>Bubble chart generation is temporarily disabled.</p></body></html>')\n");
                break;
            default:
                throw new IOException("Unsupported report type: " + task.reportType);
        }

        File tempScript = File.createTempFile("metapilot_report_", ".py");
        tempScript.deleteOnExit();
        try (FileWriter writer = new FileWriter(tempScript)) {
            writer.write(pythonScript.toString());
        }

        return List.of(request.pythonHint, tempScript.getAbsolutePath());
    }

    private static void appendTableReport(StringBuilder pythonScript, String inputPath, String outputPath,
                                          String methodName, String title) {
        pythonScript.append("df = pd.read_csv('").append(inputPath).append("', sep='\\t')\n");
        pythonScript.append("fig = rg.").append(methodName).append("(df, metadata)\n");
        pythonScript.append("rg.generate_html_report(fig, '").append(title).append("', '")
                .append(outputPath).append("')\n");
    }

    private static final class ReportCallable implements Callable<TaskResult> {
        private final ReportTaskSpec task;
        private final Request request;
        private final File scriptsDir;

        ReportCallable(ReportTaskSpec task, Request request, File scriptsDir) {
            this.task = task;
            this.request = request;
            this.scriptsDir = scriptsDir;
        }

        @Override
        public TaskResult call() {
            try {
                if (request.pythonHint == null || request.pythonHint.isBlank()) {
                    throw new IOException("Resolved Python command is empty");
                }
                List<String> command = buildPythonCommand(task, request, scriptsDir);
                ProcessBuilder pb = new ProcessBuilder(command);
                pb.redirectErrorStream(true);
                Process process = pb.start();
                List<String> output = Collections.synchronizedList(new ArrayList<>());
                try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        output.add(line);
                        LOGGER.debug("Python report output: {}", line);
                    }
                }
                int exit = process.waitFor();
                if (exit != 0) {
                    return TaskResult.fail(task.reportType, task.outputFile,
                            "Python report failed for " + task.reportType + " (exit=" + exit
                                    + "): " + String.join(" | ", output));
                }
                if (!task.outputFile.isFile() || task.outputFile.length() == 0) {
                    return TaskResult.fail(task.reportType, task.outputFile,
                            "Report file not generated or empty: " + task.outputFile.getAbsolutePath());
                }
                return TaskResult.ok(task.reportType, task.outputFile);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                return TaskResult.fail(task.reportType, task.outputFile,
                        "Python report interrupted for " + task.reportType);
            } catch (Exception e) {
                return TaskResult.fail(task.reportType, task.outputFile,
                        "Python report failed for " + task.reportType + ": " + e.getMessage());
            }
        }
    }

    private static final class TaskResult {
        final boolean success;
        final String reportType;
        final File outputFile;
        final String message;

        private TaskResult(boolean success, String reportType, File outputFile, String message) {
            this.success = success;
            this.reportType = reportType;
            this.outputFile = outputFile;
            this.message = message;
        }

        static TaskResult ok(String reportType, File outputFile) {
            return new TaskResult(true, reportType, outputFile, null);
        }

        static TaskResult fail(String reportType, File outputFile, String message) {
            return new TaskResult(false, reportType, outputFile, message);
        }
    }
}

package uk.ac.quadram.metapilot.cli;

import org.json.JSONObject;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.resources.ToolKey;
import uk.ac.quadram.metapilot.core.report.PythonEnvironmentResolver;
import uk.ac.quadram.metapilot.core.runner.DiaMarkerHeadlessRunner;
import uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItemResolver;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.ParameterIORegistry;
import uk.ac.quadram.metapilot.workflow.WorkflowType;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

final class HpcStudyRunner {

    private HpcStudyRunner() {
    }

    static int validateStudy(String paramsFile, String resourcesOverride, String magdbOverride) {
        try {
            HpcSiteConfig.ensureSiteConfig(false);
            Path resourcesCfg = HpcSiteConfig.resolveResourcesConfigPath(resourcesOverride);
            Path magdbCfg = HpcSiteConfig.resolveMagDbConfigPath(magdbOverride);
            MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            ParameterIORegistry.registerAllParameterIOs();

            if (!BaseParameterIO.isValidParameterFilePath(paramsFile)) {
                System.err.println("ERROR: unsupported or invalid study parameter file: " + new File(paramsFile).getAbsolutePath());
                ParameterIORegistry.printSupportedWorkflowTypes();
                return 2;
            }

            Object loaded = BaseParameterIO.loadParametersFromFilePath(paramsFile);
            if (loaded == null) {
                System.err.println("ERROR: failed to load study parameters from: " + new File(paramsFile).getAbsolutePath());
                return 2;
            }

            String workflow = readWorkflowType(paramsFile);
            MetaPilotResourceConfig rc = MetaPilotResourceConfig.load(resourcesCfg);
            List<String> issues = rc.validateMinimalLinuxCli();
            if (!issues.isEmpty()) {
                System.err.println("ERROR: resources config is missing required keys: " + resourcesCfg.toAbsolutePath());
                for (String issue : issues) {
                    System.err.println("  - " + issue);
                }
                return 2;
            }

            if (WorkflowType.DIA_MARKER.getAnalysisType().equals(workflow)) {
                return validateDiaMarker((MetaParameterDiaMarker) loaded, rc, resourcesCfg, magdbCfg.toFile());
            }
            if (WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType().equals(workflow)) {
                return validateMagFragpipe((DDAMagFragPipeParameters) loaded, rc, resourcesCfg, magdbCfg.toFile());
            }

            System.err.println("ERROR: workflow is registered but not supported by the headless CLI yet: " + workflow);
            return 2;
        } catch (ClassCastException e) {
            System.err.println("ERROR: parameter object does not match workflow type in study JSON: " + new File(paramsFile).getAbsolutePath());
            return 2;
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return 1;
        }
    }

    static int runStudy(String paramsFile, String resourcesOverride, String magdbOverride) {
        return runStudy(paramsFile, resourcesOverride, magdbOverride, 0);
    }

    static int runStudy(String paramsFile, String resourcesOverride, String magdbOverride,
            int fragpipeParallelTasksOverride) {
        PrintWriter log = null;
        PrintStream originalErr = System.err;
        try {
            HpcSiteConfig.ensureSiteConfig(false);
            Path resourcesCfg = HpcSiteConfig.resolveResourcesConfigPath(resourcesOverride);
            Path magdbCfg = HpcSiteConfig.resolveMagDbConfigPath(magdbOverride);
            MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            ParameterIORegistry.registerAllParameterIOs();

            int validation = validateStudy(paramsFile, resourcesOverride, magdbOverride);
            if (validation != 0) {
                return validation;
            }

            Object loaded = BaseParameterIO.loadParametersFromFilePath(paramsFile);
            if (loaded == null) {
                System.err.println("ERROR: failed to reload study parameters after validation: " + new File(paramsFile).getAbsolutePath());
                return 2;
            }
            String workflow = readWorkflowType(paramsFile);
            MetaPilotResourceConfig rc = MetaPilotResourceConfig.load(resourcesCfg);

            // Open logging.log in the study result directory before the runner starts,
            // then redirect System.err so stack traces also land in the file.
            File resultDir = resolveResultDir(loaded);
            if (resultDir != null) {
                log = openLogWriter(resultDir);
                if (log != null) {
                    System.setErr(teeStream(originalErr, log));
                    logLine(log, "INFO ", "MetaPilot HPC run started  workflow=" + workflow
                            + "  params=" + paramsFile);
                }
            }

            ProgressCallback cb = loggingCallback(log);

            if (WorkflowType.DIA_MARKER.getAnalysisType().equals(workflow)) {
                File diannExe = rc.resolveToolPath(ToolKey.DIANN, null).map(Path::toFile).orElse(null);
                File pythonExe = rc.resolveToolPath(ToolKey.PYTHON, null).map(Path::toFile).orElse(null);
                File funcDb = rc.resolveToolPath(ToolKey.FUNCDB, null).map(Path::toFile).orElse(null);
                DiaMarkerHeadlessRunner.run(
                        (MetaParameterDiaMarker) loaded,
                        diannExe,
                        magdbCfg.toFile(),
                        cb,
                        PythonEnvironmentResolver.inferCondaPrefixFromPython(pythonExe),
                        funcDb,
                        reportGenerationEnabled(rc),
                        null,
                        rc,
                        exportMetaAtlasEnabled(rc));
                logLine(log, "INFO ", "MetaPilot HPC run finished  workflow=" + workflow);
                return 0;
            }

            if (WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType().equals(workflow)) {
                DDAMagFragPipeParameters params = (DDAMagFragPipeParameters) loaded;
                if (fragpipeParallelTasksOverride > 0) {
                    params.setFragpipeParallelTasks(fragpipeParallelTasksOverride);
                }
                File fragpipeExe = rc.resolveToolPath(ToolKey.FRAGPIPE, null).map(Path::toFile).orElse(null);
                File funcDb = rc.resolveToolPath(ToolKey.FUNCDB, null).map(Path::toFile).orElse(null);
                MagFragPipeHeadlessRunner.runWithInstalledDefaults(params, fragpipeExe, magdbCfg.toFile(),
                        funcDb, cb, reportGenerationEnabled(rc), null, rc, exportMetaAtlasEnabled(rc));
                logLine(log, "INFO ", "MetaPilot HPC run finished  workflow=" + workflow);
                return 0;
            }

            System.err.println("ERROR: unsupported workflow: " + workflow);
            return 2;
        } catch (Exception e) {
            logLine(log, "ERROR", "MetaPilot HPC run failed: " + e.getMessage());
            e.printStackTrace(System.err);
            return 1;
        } finally {
            System.setErr(originalErr);
            if (log != null) {
                log.flush();
                log.close();
            }
        }
    }

    private static int validateDiaMarker(MetaParameterDiaMarker params, MetaPilotResourceConfig rc, Path resourcesCfg, File magdbCfg) {
        DiaRefinementMode markerMode = params.getMarkerRefinementMode();
        if (markerMode == DiaRefinementMode.DIRECT) {
            System.err.println("ERROR: DIA_MARKER does not support DIRECT mode. DIRECT is only valid for DIA Seed workflows.");
            return 2;
        }
        File diannExe = rc.resolveToolPath(ToolKey.DIANN, null).map(Path::toFile).orElse(null);
        if (diannExe == null || !diannExe.exists()) {
            System.err.println("ERROR: DIANN is not set or does not exist in resources config: "
                    + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(diannExe) + ")");
            return 2;
        }
        if (params.getMetadata() == null || params.getMetadata().getRawFiles() == null
                || params.getMetadata().getRawFiles().length == 0) {
            System.err.println("ERROR: no raw files found in study parameters.");
            return 2;
        }
        if (params.getResult() == null || params.getResult().isBlank()) {
            System.err.println("ERROR: output_directory is required in study parameters.");
            return 2;
        }
        MagDbItem mag = params.getMagDbItem();
        if (mag == null && params.getMicroDb() != null && !params.getMicroDb().isBlank()) {
            mag = MagDbItemResolver.resolve(
                    params.getMicroDb(),
                    params.getMagDbIdentifier(),
                    params.getMagDbCatalogueVersion(),
                    magdbCfg);
            if (mag != null) {
                params.setMagDbItem(mag);
            }
        }
        if (mag == null) {
            System.err.println("ERROR: MAG database entry could not be resolved from site registry: "
                    + magdbCfg.getAbsolutePath()
                    + " (micro_db=" + String.valueOf(params.getMicroDb())
                    + ", identifier=" + params.getMagDbIdentifier()
                    + ", version=" + params.getMagDbCatalogueVersion() + ")");
            return 2;
        }
        System.out.println("OK: DIA_MARKER study parameters validated. marker_mode="
                + markerMode.getPersistedValue() + " (" + markerMode.getDisplayLabel() + ")");
        if (markerMode == DiaRefinementMode.DEEP_PER_FILE && params.isRepresentativeClusteringEnabled()) {
            System.out.println("OK: DEEP_PER_FILE representative clustering enabled. preset="
                    + params.getRepresentativeClusteringPreset()
                    + ", peptide_jaccard=" + params.getRepresentativePeptideJaccardThreshold()
                    + ", peptide_containment=" + params.getRepresentativePeptideContainmentThreshold()
                    + ", genome_top_n_diagnostic=" + params.getRepresentativeGenomeTopN()
                    + ", genome_overlap_used_for_clustering=false"
                    + ", min_marker_peptides=" + params.getRepresentativeMinMarkerPeptides()
                    + ", grouping_dry_run=" + params.isRepresentativeClusteringDryRun());
        }
        return 0;
    }

    private static int validateMagFragpipe(DDAMagFragPipeParameters params, MetaPilotResourceConfig rc, Path resourcesCfg, File magdbCfg) {
        File fragpipeExe = rc.resolveToolPath(ToolKey.FRAGPIPE, null).map(Path::toFile).orElse(null);
        if (fragpipeExe == null || !fragpipeExe.exists()) {
            System.err.println("ERROR: FRAGPIPE is not set or does not exist in resources config: "
                    + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(fragpipeExe) + ")");
            return 2;
        }
        if (params.getMetadata() == null || params.getMetadata().getRawFiles() == null
                || params.getMetadata().getRawFiles().length == 0) {
            System.err.println("ERROR: no raw files found in study parameters.");
            return 2;
        }
        if (params.getResult() == null || params.getResult().isBlank()) {
            System.err.println("ERROR: output_directory is required in study parameters.");
            return 2;
        }
        try {
            // Search mode (Closed/Open, custom template overrides) is read entirely from
            // the loaded study parameters (search_mode / open_search_marker_stage /
            // custom_basic_workflow_path / custom_combined_workflow_path) - only
            // pre-validate whichever bundled templates those settings will actually need.
            if (params.getCustomBasicWorkflowPath() == null || params.getCustomBasicWorkflowPath().isBlank()) {
                String basicTemplate = params.isOpenSearchAtMarkerStage()
                        ? MagFragPipeHeadlessRunner.OPEN_SEARCH_WORKFLOW
                        : MagFragPipeHeadlessRunner.BASIC_SEARCH_WORKFLOW;
                MagFragPipeHeadlessRunner.requireInstalledWorkflow(fragpipeExe, basicTemplate);
            }
            if (params.getCustomCombinedWorkflowPath() == null || params.getCustomCombinedWorkflowPath().isBlank()) {
                MagFragPipeHeadlessRunner.requireInstalledWorkflow(fragpipeExe, MagFragPipeHeadlessRunner.LFQ_MBR_WORKFLOW);
                if (params.isOpenSearch()) {
                    MagFragPipeHeadlessRunner.requireInstalledWorkflow(fragpipeExe, MagFragPipeHeadlessRunner.OPEN_SEARCH_WORKFLOW);
                }
            }
        } catch (IOException e) {
            System.err.println("ERROR: required FragPipe workflow is not available in the installed FragPipe folder: "
                    + e.getMessage());
            return 2;
        }
        MagDbItem mag = params.getMagDbItem();
        if (mag == null) {
            mag = MagDbItemResolver.resolve(
                    params.getMagDbCatalogueId(),
                    params.getMagDbIdentifier(),
                    params.getMagDbVersion(),
                    magdbCfg);
            if (mag != null) {
                params.setMagDbItem(mag);
            }
        }
        if (mag == null) {
            System.err.println("ERROR: MAG database entry could not be resolved from site registry: "
                    + magdbCfg.getAbsolutePath()
                    + " (catalogue=" + String.valueOf(params.getMagDbCatalogueId())
                    + ", identifier=" + params.getMagDbIdentifier()
                    + ", version=" + String.valueOf(params.getMagDbVersion()) + ")");
            return 2;
        }
        System.out.println("OK: DDA_MAG_FRAGPIPE study parameters validated.");
        return 0;
    }

    private static String readWorkflowType(String paramsFile) throws Exception {
        try (FileReader reader = new FileReader(paramsFile)) {
            StringBuilder content = new StringBuilder();
            char[] buffer = new char[1024];
            int bytesRead;
            while ((bytesRead = reader.read(buffer)) != -1) {
                content.append(buffer, 0, bytesRead);
            }
            JSONObject json = new JSONObject(content.toString());
            return json.optString("workflow", "").trim();
        }
    }

    // ── logging helpers ──────────────────────────────────────────────────────

    private static final SimpleDateFormat LOG_DATE = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    private static File resolveResultDir(Object loaded) {
        if (loaded instanceof MetaParameterDiaMarker) {
            String r = ((MetaParameterDiaMarker) loaded).getResult();
            return (r != null && !r.isBlank()) ? new File(r) : null;
        }
        if (loaded instanceof DDAMagFragPipeParameters) {
            String r = ((DDAMagFragPipeParameters) loaded).getResult();
            return (r != null && !r.isBlank()) ? new File(r) : null;
        }
        return null;
    }

    private static PrintWriter openLogWriter(File resultDir) {
        try {
            resultDir.mkdirs();
            File logFile = new File(resultDir, "logging.log");
            PrintWriter pw = new PrintWriter(new FileWriter(logFile, true), true);
            logLine(pw, "INFO ", "========================================");
            logLine(pw, "INFO ", "MetaPilot HPC logging started: " + logFile.getAbsolutePath());
            logLine(pw, "INFO ", "========================================");
            System.out.println("Logging to: " + logFile.getAbsolutePath());
            return pw;
        } catch (IOException e) {
            System.err.println("WARNING: could not open logging.log in " + resultDir + ": " + e.getMessage());
            return null;
        }
    }

    /** Writes a timestamped line to the log file (no-op when log is null). */
    private static void logLine(PrintWriter log, String level, String message) {
        if (log == null) return;
        synchronized (LOG_DATE) {
            log.println(LOG_DATE.format(new Date()) + " " + level + " - " + message);
        }
    }

    /**
     * Returns a PrintStream that writes every byte to both {@code original} and
     * the given {@code log} PrintWriter (line-buffered to the log).
     */
    private static PrintStream teeStream(PrintStream original, PrintWriter log) {
        OutputStream tee = new OutputStream() {
            private final StringBuilder line = new StringBuilder();

            @Override
            public void write(int b) throws IOException {
                original.write(b);
                char c = (char) b;
                if (c == '\n') {
                    logLine(log, "ERROR", line.toString());
                    line.setLength(0);
                } else if (c != '\r') {
                    line.append(c);
                }
            }

            @Override
            public void flush() throws IOException {
                original.flush();
                if (line.length() > 0) {
                    logLine(log, "ERROR", line.toString());
                    line.setLength(0);
                }
            }
        };
        return new PrintStream(tee, true);
    }

    /** ProgressCallback that mirrors every message to both stdout and the log file. */
    private static ProgressCallback loggingCallback(PrintWriter log) {
        return new ProgressCallback() {
            @Override
            public void onProgress(double workDone, double max) {
            }

            @Override
            public void onMessage(String message) {
                if (message == null || message.isBlank()) return;
                String line;
                synchronized (LOG_DATE) {
                    line = LOG_DATE.format(new Date()) + " INFO  - " + message;
                }
                System.out.println(line);
                if (log != null) log.println(line);
            }

            @Override
            public void onComplete(Object result) {
                onMessage("Done: " + result);
            }

            @Override
            public void onError(Exception error) {
                String msg = "ERROR: " + (error != null ? error.getMessage() : "unknown");
                System.err.println(msg);
                logLine(log, "ERROR", msg);
            }
        };
    }

    // ────────────────────────────────────────────────────────────────────────

    private static boolean reportGenerationEnabled(MetaPilotResourceConfig rc) {
        if (rc == null) {
            return true;
        }
        return rc.getRaw("GENERATE_REPORT")
                .map(v -> {
                    String t = v.trim().toLowerCase(java.util.Locale.ROOT);
                    return !(t.equals("false") || t.equals("0") || t.equals("no") || t.equals("n"));
                })
                .orElse(true);
    }

    /**
     * Opt-in MetaAtlas Parquet export (taxon_abundance.parquet / function_abundance.parquet),
     * controlled by an {@code EXPORT_METAATLAS} key in the resources config - mirrors
     * {@link #reportGenerationEnabled} but defaults to false (off), since this is a new,
     * separate export layer that must be explicitly requested rather than one enabled by
     * default. Not exposed in the Windows GUI.
     */
    private static boolean exportMetaAtlasEnabled(MetaPilotResourceConfig rc) {
        if (rc == null) {
            return false;
        }
        return rc.getRaw("EXPORT_METAATLAS")
                .map(v -> {
                    String t = v.trim().toLowerCase(java.util.Locale.ROOT);
                    return t.equals("true") || t.equals("1") || t.equals("yes") || t.equals("y");
                })
                .orElse(false);
    }

    private static String describeFile(File file) {
        return file == null ? "<not configured>" : file.getAbsolutePath();
    }
}

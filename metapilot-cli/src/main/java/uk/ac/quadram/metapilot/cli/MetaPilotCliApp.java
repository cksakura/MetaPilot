package uk.ac.quadram.metapilot.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Optional;

/**
 * Headless MetaPilot CLI entrypoint (no JavaFX).
 *
 * This module is intended to become the primary runner for Linux/HPC execution.
 * Initial version provides stable CLI scaffolding; workflow execution will be
 * wired up incrementally as JavaFX dependencies are removed from task code.
 */
public final class MetaPilotCliApp {

    private MetaPilotCliApp() {
    }

    public static void main(String[] args) {
        int exitCode = run(args);
        if (exitCode != 0) {
            System.exit(exitCode);
        }
    }

    static int run(String[] args) {
        HpcSiteConfig.bindInstallHome();
        if (args != null && args.length > 0 && "run".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runStudyCommand(subArgs, false);
        }
        if (args != null && args.length > 0 && "init".equalsIgnoreCase(args[0])) {
            return runInitCommand();
        }
        if (args != null && args.length > 0 && "params-validate".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runStudyCommand(subArgs, true);
        }
        if (args != null && args.length > 0 && "mag-register".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return HpcMagRegistryCommands.runMagRegister(subArgs);
        }
        if (args != null && args.length > 0 && "mag-list".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return HpcMagRegistryCommands.runMagList(subArgs);
        }
        if (args != null && args.length > 0 && "mag-download".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runMagDownload(subArgs);
        }
        if (args != null && args.length > 0 && "dia-marker".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDiaMarker(subArgs);
        }
        if (args != null && args.length > 0 && "mag-fragpipe".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runMagFragpipe(subArgs);
        }
        if (args != null && args.length > 0 && "resources-check".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runResourcesCheck(subArgs);
        }
        if (args != null && args.length > 0 && "resources-set".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runResourcesSet(subArgs);
        }
        if (args != null && args.length > 0 && "deepdetect-predict".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDeepDetectPredict(subArgs);
        }
        if (args != null && args.length > 0 && "deepdetect-magdb".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDeepDetectMagDb(subArgs);
        }

        Options options = buildOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printHelp(options);
            return 0;
        }

        if (cmd.hasOption("version")) {
            System.out.println("MetaPilot CLI " + getVersion());
            return 0;
        }

        printHelp(options);
        return 0;
    }

    private static int runInitCommand() {
        try {
            java.util.List<String> created = HpcSiteConfig.ensureSiteConfig(true);
            if (created.isEmpty()) {
                System.out.println("Site config already present under " + HpcSiteConfig.siteConfigDir());
            }
            return 0;
        } catch (Exception e) {
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runStudyCommand(String[] args, boolean validateOnly) {
        Options options = buildStudyOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printStudyHelp(options, validateOnly);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printStudyHelp(options, validateOnly);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printStudyHelp(options, validateOnly);
            return 2;
        }

        if (validateOnly) {
            return HpcStudyRunner.validateStudy(paramsFile, cmd.getOptionValue("resources"), cmd.getOptionValue("magdb"));
        }
        int fragpipeParallelTasks = parsePositiveIntOption(cmd, "fragpipe-parallel-tasks", "fragpipe-concurrency");
        return HpcStudyRunner.runStudy(paramsFile, cmd.getOptionValue("resources"), cmd.getOptionValue("magdb"),
                fragpipeParallelTasks);
    }

    private static int runDiaMarker(String[] args) {
        Options options = buildDiaMarkerOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDiaMarkerHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printDiaMarkerHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printDiaMarkerHelp(options);
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File diannExe = rc.getRaw("DIANN").map(java.io.File::new).orElse(null);
            if (diannExe == null || !diannExe.exists()) {
                System.err.println("ERROR: DIANN is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(diannExe) + ")");
                return 2;
            }

            // Optional: infer conda prefix from PYTHON (for LD_LIBRARY_PATH fix on older systems).
            java.io.File pythonExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.PYTHON, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            java.io.File condaPrefix = uk.ac.quadram.metapilot.core.report.PythonEnvironmentResolver
                    .inferCondaPrefixFromPython(pythonExe);
            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile).orElse(null);

            uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker p =
                    new uk.ac.quadram.metapilot.task.dia.par.DiaMarkerParameterIO()
                            .loadParameters(new java.io.File(paramsFile));

            if (p == null) {
                System.err.println("ERROR: failed to load DIA Marker params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }
            uk.ac.quadram.metapilot.task.dia.DiaRefinementMode markerMode = p.getMarkerRefinementMode();
            if (markerMode == uk.ac.quadram.metapilot.task.dia.DiaRefinementMode.DIRECT) {
                System.err.println("ERROR: dia-marker does not support DIRECT mode. DIRECT is only valid for DIA Seed workflows.");
                return 2;
            }
            System.out.println("DIA Marker mode: " + markerMode.getPersistedValue() + " ("
                    + markerMode.getDisplayLabel() + ")");

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("dia-marker", error); }
            };

            uk.ac.quadram.metapilot.core.runner.DiaMarkerHeadlessRunner.run(
                    p,
                    diannExe,
                    magdbCfg.toFile(),
                    cb,
                    condaPrefix,
                    funcDb,
                    reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                    null,
                    rc
            );
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: dia-marker failed for params file: " + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runDeepDetectPredict(String[] args) {
        Options options = buildDeepDetectPredictOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDeepDetectPredictHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printDeepDetectPredictHelp(options);
            return 0;
        }

        String input = cmd.getOptionValue("input");
        String output = cmd.getOptionValue("output");
        if (input == null || output == null) {
            System.err.println("ERROR: --input and --output are required.");
            printDeepDetectPredictHelp(options);
            return 2;
        }

        Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        String regular = cmd.getOptionValue("regular");

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File pythonExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.PYTHON, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            if (pythonExe == null || !pythonExe.exists()) {
                System.err.println("ERROR: PYTHON is not configured or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(pythonExe) + ")");
                return 2;
            }

            java.io.File mainPy = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.DEEPDETECT, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            if (mainPy == null || !mainPy.exists()) {
                System.err.println("ERROR: DEEPDETECT main.py is not configured or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(mainPy) + ")");
                return 2;
            }

            uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter par =
                    new uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter();
            par.setInput(new java.io.File(input));
            par.setOutput(new java.io.File(output));
            if (cmd.hasOption("protease")) par.setProtease(cmd.getOptionValue("protease"));
            if (cmd.hasOption("missed_cleavages")) par.setMissCleavage(parseIntOr(cmd.getOptionValue("missed_cleavages"), par.getMissCleavage()));
            if (cmd.hasOption("min_len")) par.setMinLength(parseIntOr(cmd.getOptionValue("min_len"), par.getMinLength()));
            if (cmd.hasOption("max_len")) par.setMaxLength(parseIntOr(cmd.getOptionValue("max_len"), par.getMaxLength()));

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("deepdetect-predict", error); }
            };

            int exit = uk.ac.quadram.metapilot.core.deepdetect.DeepDetectPythonRunner.run(pythonExe, mainPy, par, regular, cb);
            return exit == 0 ? 0 : 1;
        } catch (Exception e) {
            System.err.println("ERROR: deepdetect-predict failed for input=" + new java.io.File(input).getAbsolutePath()
                    + " output=" + new java.io.File(output).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runDeepDetectMagDb(String[] args) {
        Options options = buildDeepDetectMagDbOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDeepDetectMagDbHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printDeepDetectMagDbHelp(options);
            return 0;
        }

        String catalogue = cmd.getOptionValue("catalogue");
        String version = cmd.getOptionValue("version");
        String out = cmd.getOptionValue("out");
        if (catalogue == null || version == null) {
            System.err.println("ERROR: --catalogue and --version are required (e.g. --catalogue human-gut --version v2.0.2).");
            printDeepDetectMagDbHelp(options);
            return 2;
        }
        if (out == null || out.isBlank()) {
            out = HpcSiteConfig.installHome().resolve("magdb").toString();
        }

        Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        String regular = cmd.getOptionValue("regular");
        int threads = parseIntOr(cmd.getOptionValue("threads"), 0);

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File pythonExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.PYTHON, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            if (pythonExe == null || !pythonExe.exists()) {
                System.err.println("ERROR: PYTHON is not configured or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(pythonExe) + ")");
                return 2;
            }

            java.io.File mainPy = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.DEEPDETECT, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            if (mainPy == null || !mainPy.exists()) {
                System.err.println("ERROR: DEEPDETECT main.py is not configured or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(mainPy) + ")");
                return 2;
            }

            java.io.File versionDir = new java.io.File(new java.io.File(out, catalogue), version);
            java.io.File originalDbDir = new java.io.File(versionDir, "original_db");
            if (!originalDbDir.exists()) {
                System.err.println("ERROR: MAG DB not found (missing original_db): " + originalDbDir.getAbsolutePath());
                System.err.println("Hint: run mag-download first, using the same --out/--catalogue/--version.");
                return 2;
            }
            java.io.File predictedDir = new java.io.File(versionDir, "predicted");

            uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter template =
                    new uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter();
            if (cmd.hasOption("protease")) template.setProtease(cmd.getOptionValue("protease"));
            if (cmd.hasOption("missed_cleavages")) template.setMissCleavage(parseIntOr(cmd.getOptionValue("missed_cleavages"), template.getMissCleavage()));
            if (cmd.hasOption("min_len")) template.setMinLength(parseIntOr(cmd.getOptionValue("min_len"), template.getMinLength()));
            if (cmd.hasOption("max_len")) template.setMaxLength(parseIntOr(cmd.getOptionValue("max_len"), template.getMaxLength()));

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) {
                    if (max > 0) System.out.printf("Progress: %.1f%%%n", (workDone / max) * 100.0);
                }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("deepdetect-magdb", error); }
            };

            uk.ac.quadram.metapilot.core.deepdetect.DeepDetectMagDbHeadlessRunner.Result r =
                    uk.ac.quadram.metapilot.core.deepdetect.DeepDetectMagDbHeadlessRunner.runAllGenomes(
                            originalDbDir, predictedDir, pythonExe, mainPy, template, regular, threads, cb);
            System.out.println("Predicted folder: " + r.predictedDir.getAbsolutePath());
            System.out.println("Total genomes: " + r.totalGenomes);
            System.out.println("Skipped existing: " + r.skippedExisting);
            System.out.println("Predicted OK: " + r.predictedGenomes);
            System.out.println("Failed: " + r.failedGenomes);
            return r.failedGenomes == 0 ? 0 : 1;
        } catch (Exception e) {
            System.err.println("ERROR: deepdetect-magdb failed for catalogue=" + catalogue
                    + " version=" + version + " out=" + new java.io.File(out).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runMagFragpipe(String[] args) {
        Options options = buildMagFragpipeOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printMagFragpipeHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printMagFragpipeHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        String workflowFile = cmd.getOptionValue("workflow");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printMagFragpipeHelp(options);
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.util.List<String> issues = rc.validateMinimalLinuxCli();
            if (!issues.isEmpty()) {
                System.err.println("ERROR: resources config is missing required keys: " + resourcesCfg.toAbsolutePath());
                for (String issue : issues) System.err.println("  - " + issue);
                return 2;
            }

            java.io.File fragpipeExe = rc.getRaw("FRAGPIPE").map(java.io.File::new).orElse(null);
            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);
            if (fragpipeExe == null || !fragpipeExe.exists()) {
                System.err.println("ERROR: FRAGPIPE is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(fragpipeExe) + ")");
                return 2;
            }

            uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters p =
                    new uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load MAG FragPipe params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }
            int fragpipeParallelTasks = parsePositiveIntOption(cmd, "fragpipe-parallel-tasks", "fragpipe-concurrency");
            if (fragpipeParallelTasks > 0) {
                p.setFragpipeParallelTasks(fragpipeParallelTasks);
            }
            java.io.File workflow = null;
            if (workflowFile != null && !workflowFile.isBlank()) {
                workflow = new java.io.File(workflowFile);
                if (!workflow.exists()) {
                    System.err.println("ERROR: --workflow does not exist: " + workflow.getAbsolutePath());
                    return 2;
                }
                if (!workflow.isFile()) {
                    System.err.println("ERROR: --workflow is not a file: " + workflow.getAbsolutePath());
                    return 2;
                }
            } else {
                try {
                    uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.requireInstalledWorkflow(
                            fragpipeExe, uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.BASIC_SEARCH_WORKFLOW);
                    uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.requireInstalledWorkflow(
                            fragpipeExe, uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.LFQ_MBR_WORKFLOW);
                } catch (java.io.IOException e) {
                    System.err.println("ERROR: required FragPipe workflows are not available in the installed FragPipe folder: "
                            + e.getMessage());
                    return 2;
                }
            }

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("mag-fragpipe", error); }
            };

            if (workflow != null) {
                uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.run(
                        p,
                        fragpipeExe,
                        workflow,
                        magdbCfg.toFile(),
                        funcDb,
                        cb,
                        reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                        null,
                        rc
                );
            } else {
                uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.runWithInstalledDefaults(
                        p,
                        fragpipeExe,
                        magdbCfg.toFile(),
                        funcDb,
                        cb,
                        reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                        null,
                        rc
                );
            }
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: mag-fragpipe failed for params=" + new java.io.File(paramsFile).getAbsolutePath()
                    + " workflow=" + (workflowFile == null ? "<installed defaults>" : new java.io.File(workflowFile).getAbsolutePath()));
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runResourcesCheck(String[] args) {
        Options options = buildResourcesCheckOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printResourcesCheckHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printResourcesCheckHelp(options);
            return 0;
        }

        java.nio.file.Path cfg = resolveResourcesConfigPath(cmd.getOptionValue("file"));
        System.out.println("Using resources config: " + cfg);

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(cfg);
            java.util.List<String> issues = rc.validateMinimalLinuxCli();
            if (issues.isEmpty()) {
                System.out.println("OK: required tools/resources are configured.");
                return 0;
            }
            System.err.println("Missing/invalid resources:");
            for (String issue : issues) {
                System.err.println("  - " + issue);
            }
            return 1;
        } catch (Exception e) {
            System.err.println("ERROR: failed to load/validate resources config: " + cfg.toAbsolutePath()
                    + " (" + e.getMessage() + ")");
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runResourcesSet(String[] args) {
        Options options = buildResourcesSetOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printResourcesSetHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printResourcesSetHelp(options);
            return 0;
        }

        Path cfg = resolveResourcesConfigPath(cmd.getOptionValue("file"));
        String[] sets = cmd.getOptionValues("set");
        String[] unsets = cmd.getOptionValues("unset");
        boolean dryRun = cmd.hasOption("dry-run");
        boolean doPrint = cmd.hasOption("print");
        boolean template = cmd.hasOption("template");
        boolean force = cmd.hasOption("force");

        boolean hasMutations = !((sets == null || sets.length == 0) && (unsets == null || unsets.length == 0));

        try {
            LinkedHashMap<String, String> kv = new LinkedHashMap<>();
            boolean cfgExists = Files.exists(cfg);

            // First-time usability: allow generating a starter template without manual editing.
            if (!hasMutations) {
                if (template) {
                    // Allow previewing a template even if the file exists (dry-run/print use-case).
                    if (cfgExists && !dryRun && !force) {
                        System.err.println("ERROR: resources config already exists: " + cfg);
                        System.err.println("Use --force to overwrite, or use --set/--unset to modify.");
                        return 2;
                    }
                    kv.putAll(defaultLinuxResourceTemplate());
                } else if (!cfgExists) {
                    kv.putAll(defaultLinuxResourceTemplate());
                } else {
                    System.err.println("ERROR: provide at least one --set KEY=VALUE or --unset KEY.");
                    System.err.println("Tip: use --template to generate a starter file.");
                    printResourcesSetHelp(options);
                    return 2;
                }
            } else if (cfgExists) {
                kv.putAll(readKeyValue(cfg));
            }

            if (unsets != null) {
                for (String k : unsets) {
                    if (k == null) continue;
                    String key = k.trim();
                    if (!key.isEmpty()) kv.remove(key);
                }
            }

            if (sets != null) {
                for (String s : sets) {
                    if (s == null) continue;
                    int eq = s.indexOf('=');
                    if (eq < 1) {
                        System.err.println("ERROR: invalid --set value (expected KEY=VALUE): " + s);
                        return 2;
                    }
                    String key = s.substring(0, eq).trim();
                    String value = s.substring(eq + 1).trim();
                    if (key.isEmpty()) {
                        System.err.println("ERROR: invalid --set value (empty KEY): " + s);
                        return 2;
                    }
                    kv.put(key, value);
                }
            }

            List<String> lines = renderParamsFile(kv);

            System.out.println("Using resources config: " + cfg);
            if (dryRun) {
                System.out.println("Dry-run: not writing file.");
            } else {
                Path parent = cfg.getParent();
                if (parent != null) Files.createDirectories(parent);
                try (BufferedWriter bw = Files.newBufferedWriter(cfg, StandardCharsets.UTF_8)) {
                    for (String line : lines) {
                        bw.write(line);
                        bw.newLine();
                    }
                }
                System.out.println((!cfgExists ? "Created" : "Updated") + " resources config.");
            }

            if (doPrint) {
                System.out.println();
                for (String line : lines) System.out.println(line);
            }

            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: failed to update resources config: " + cfg.toAbsolutePath()
                    + " (" + e.getMessage() + ")");
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runMagDownload(String[] args) {
        Options options = buildMagDownloadOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printMagDownloadHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printMagDownloadHelp(options);
            return 0;
        }

        String catalogue = cmd.getOptionValue("catalogue");
        String version = cmd.getOptionValue("version");
        String out = cmd.getOptionValue("out");
        if (catalogue == null || version == null) {
            System.err.println("ERROR: --catalogue and --version are required.");
            printMagDownloadHelp(options);
            return 2;
        }
        if (out == null || out.isBlank()) {
            out = HpcSiteConfig.installHome().resolve("magdb").toString();
        }

        boolean resume = cmd.hasOption("resume");
        int threads = parseIntOr(cmd.getOptionValue("threads"), 0);
        int waitHours = parseIntOr(cmd.getOptionValue("wait-hours"), 6);

        uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
            private long lastPrint = 0;
            @Override
            public void onProgress(double workDone, double max) {
                long now = System.currentTimeMillis();
                if (now - lastPrint < 500) return;
                lastPrint = now;
                if (max > 0) {
                    double pct = (workDone / max) * 100.0;
                    System.out.printf("Progress: %.1f%%%n", pct);
                }
            }

            @Override
            public void onMessage(String message) {
                if (message != null && !message.isBlank()) {
                    System.out.println(message);
                }
            }

            @Override
            public void onComplete(Object result) {
                System.out.println("Done.");
            }

            @Override
            public void onError(Exception error) {
                reportCallbackError("mag-download", error);
            }
        };

        try {
            uk.ac.quadram.metapilot.core.MGYG.MGYGMagDownloader.downloadCatalogue(
                    catalogue,
                    version,
                    new java.io.File(out),
                    resume,
                    threads,
                    waitHours,
                    cb
            );
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: mag-download failed for catalogue=" + catalogue
                    + " version=" + version + " out=" + new java.io.File(out).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static Options buildOptions() {
        Options opts = new Options();

        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());

        opts.addOption(Option.builder()
                .longOpt("version")
                .desc("Print version and exit")
                .build());

        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("Workflow parameter JSON file (future)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (future)")
                .build());

        return opts;
    }

    private static Options buildStudyOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h").longOpt("help").desc("Show help").build());
        opts.addOption(Option.builder().longOpt("params").hasArg().argName("FILE").desc("Study parameter JSON file").required().build());
        opts.addOption(Option.builder().longOpt("resources").hasArg().argName("FILE").desc("Site resources override (default: <install>/config/meta_resources.params)").build());
        opts.addOption(Option.builder().longOpt("magdb").hasArg().argName("FILE").desc("MAG registry override (default: <install>/config/MAGDB_1_0.json)").build());
        opts.addOption(Option.builder().longOpt("fragpipe-parallel-tasks").hasArg().argName("N").desc("Override MAG FragPipe per-sample concurrent task count; warns if above memory-safe estimate").build());
        opts.addOption(Option.builder().longOpt("fragpipe-concurrency").hasArg().argName("N").desc("Alias for --fragpipe-parallel-tasks").build());
        return opts;
    }

    private static void printStudyHelp(Options options, boolean validateOnly) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        String command = validateOnly ? "params-validate" : "run";
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar " + command + " [options]",
                "\n" + (validateOnly ? "Validate" : "Run") + " a headless study from a parameter JSON file.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar " + command
                        + " --params /path/study_params.json\n",
                true
        );
    }

    private static Options buildResourcesCheckOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());

        opts.addOption(Option.builder()
                .longOpt("file")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (default: <install>/config/meta_resources.params)")
                .build());

        return opts;
    }

    private static Options buildResourcesSetOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());

        opts.addOption(Option.builder()
                .longOpt("file")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (default: <install>/config/meta_resources.params)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("template")
                .desc("Generate a starter template (writes file if missing; use --force to overwrite).")
                .build());

        opts.addOption(Option.builder()
                .longOpt("force")
                .desc("Overwrite existing file (only used with --template or when no --set/--unset are provided).")
                .build());

        opts.addOption(Option.builder()
                .longOpt("set")
                .hasArgs()
                .argName("KEY=VALUE")
                .desc("Set or override a resource key (repeatable). Example: --set DIANN=/path/to/diann")
                .build());

        opts.addOption(Option.builder()
                .longOpt("unset")
                .hasArgs()
                .argName("KEY")
                .desc("Remove a key from the config (repeatable).")
                .build());

        opts.addOption(Option.builder()
                .longOpt("dry-run")
                .desc("Do not write the file; just validate/print what would be written.")
                .build());

        opts.addOption(Option.builder()
                .longOpt("print")
                .desc("Print the resulting config content.")
                .build());

        return opts;
    }

    private static Options buildMagFragpipeOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("MAG FragPipe parameters JSON (e.g. last_used_mag_fragpipe_parameters.json)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("workflow")
                .hasArg()
                .argName("FILE")
                .desc("Optional FragPipe workflow template override (.workflow). Default behavior uses installed Basic-Search.workflow for HAP/sample search and LFQ-MBR.workflow for combined search.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (default resolution same as resources-check)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("magdb")
                .hasArg()
                .argName("FILE")
                .desc("MAG DB config JSON (MAGDB_1_0.json). Default: <install>/config/MAGDB_1_0.json")
                .build());
        opts.addOption(Option.builder()
                .longOpt("no-report")
                .desc("Disable final HTML report generation")
                .build());
        opts.addOption(Option.builder()
                .longOpt("fragpipe-parallel-tasks")
                .hasArg()
                .argName("N")
                .desc("Override per-sample concurrent task count; warns if above memory-safe estimate")
                .build());
        opts.addOption(Option.builder()
                .longOpt("fragpipe-concurrency")
                .hasArg()
                .argName("N")
                .desc("Alias for --fragpipe-parallel-tasks")
                .build());
        return opts;
    }

    private static Options buildDiaMarkerOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("DIA Marker parameters JSON (e.g. last_used_dia_marker_parameters.json)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (default resolution same as resources-check)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("magdb")
                .hasArg()
                .argName("FILE")
                .desc("MAG DB config JSON (MAGDB_1_0.json). Default: <install>/config/MAGDB_1_0.json")
                .build());
        opts.addOption(Option.builder()
                .longOpt("no-report")
                .desc("Disable final HTML report generation")
                .build());
        return opts;
    }

    private static Options buildDeepDetectPredictOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (default resolution same as resources-check)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("input")
                .hasArg()
                .argName("FASTA")
                .desc("Input protein FASTA")
                .build());
        opts.addOption(Option.builder()
                .longOpt("output")
                .hasArg()
                .argName("TSV")
                .desc("Output predicted TSV")
                .build());
        opts.addOption(Option.builder()
                .longOpt("regular")
                .hasArg()
                .argName("REGEX")
                .desc("FASTA header regex (optional; DeepDetect default: >(.*?)\\\\s)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("protease")
                .hasArg()
                .argName("NAME")
                .desc("Protease (Trypsin/ArgC/Chymotrypsin/GluC/LysC/AspN/LysN/LysargiNase)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("missed_cleavages")
                .hasArg()
                .argName("N")
                .desc("Missed cleavages (default from DeepDetectParameter)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("min_len")
                .hasArg()
                .argName("N")
                .desc("Minimum peptide length")
                .build());
        opts.addOption(Option.builder()
                .longOpt("max_len")
                .hasArg()
                .argName("N")
                .desc("Maximum peptide length")
                .build());
        return opts;
    }

    private static Options buildDeepDetectMagDbOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file (default resolution same as resources-check)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("catalogue")
                .hasArg()
                .argName("ID")
                .desc("MAG catalogue id (e.g. human-gut)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("version")
                .hasArg()
                .argName("VER")
                .desc("MAG catalogue version (e.g. v2.0.2)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("out")
                .hasArg()
                .argName("DIR")
                .desc("MAG DB root folder (default: <jar-folder>/magdb). Must match mag-download --out")
                .build());
        opts.addOption(Option.builder()
                .longOpt("threads")
                .hasArg()
                .argName("N")
                .desc("Parallel DeepDetect processes (default: max(2, floor(CPU/2)); 1 only on 1-core machines)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("regular")
                .hasArg()
                .argName("REGEX")
                .desc("FASTA header regex (optional; DeepDetect default: >(.*?)\\\\s)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("protease")
                .hasArg()
                .argName("NAME")
                .desc("Protease (Trypsin/ArgC/Chymotrypsin/GluC/LysC/AspN/LysN/LysargiNase)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("missed_cleavages")
                .hasArg()
                .argName("N")
                .desc("Missed cleavages (default from DeepDetectParameter)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("min_len")
                .hasArg()
                .argName("N")
                .desc("Minimum peptide length")
                .build());
        opts.addOption(Option.builder()
                .longOpt("max_len")
                .hasArg()
                .argName("N")
                .desc("Maximum peptide length")
                .build());
        return opts;
    }

    private static void printHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar [options] <command> [command options]",
                "\nMetaPilot CLI (headless runner).",
                options,
                "\nCommands:\n"
                        + "  init              Create <install>/config templates on first use\n"
                        + "  run               Run a study from a parameter JSON file\n"
                        + "  params-validate   Validate a study parameter JSON file without running\n"
                        + "  resources-check   Validate resource configuration (meta_resources.params)\n"
                        + "  resources-set     Update resource configuration (repeatable --set KEY=VALUE)\n"
                        + "  mag-download      Download one MGYG MAG catalogue+version\n"
                        + "  mag-register      Register a local MAG catalogue repository path\n"
                        + "  mag-list          List MAG catalogues in the site registry\n"
                        + "  deepdetect-predict  Run DeepDetect PyTorch prediction (python main.py)\n"
                        + "  deepdetect-magdb    Predict peptides for all genomes in a MAG DB\n"
                        + "\nExamples:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar init\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar run --params /path/study_params.json\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar params-validate --params /path/study_params.json\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar resources-check\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar mag-download --help\n",
                true
        );
    }

    private static void printDeepDetectPredictHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar deepdetect-predict [options]",
                "\nRun DeepDetect prediction using the Python entrypoint (main.py).",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar deepdetect-predict \\\n"
                        + "    --input /path/genome.faa --output /path/genome.predicted.tsv\n",
                true
        );
    }

    private static void printDeepDetectMagDbHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar deepdetect-magdb [options]",
                "\nRun DeepDetect prediction for all genome FASTA files in a downloaded MAG catalogue+version.\n"
                        + "Folder layout must match mag-download: <out>/<catalogue>/<version>/original_db.\n"
                        + "Outputs: <out>/<catalogue>/<version>/predicted/*.predicted.tsv",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar deepdetect-magdb \\\n"
                        + "    --catalogue human-gut --version v2.0.2 --threads 4\n",
                true
        );
    }

    private static void printResourcesCheckHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar resources-check [options]",
                "\nValidate a MetaPilot resource configuration for Linux/HPC CLI execution.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar resources-check --file ~/tools/metapilot/config/meta_resources.params\n",
                true
        );
    }

    private static void printResourcesSetHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar resources-set [options]",
                "\nUpdate a MetaPilot resource configuration file (KEY=VALUE) without manual editing.",
                options,
                "\nExamples:\n"
                        + "  # Create a starter template (first time)\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar resources-set --template --print\n"
                        + "\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar resources-set \\\n"
                        + "    --set SYSTEM=linux --set RESOURCE_AUTO_DETECTED=true\n"
                        + "\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar resources-set \\\n"
                        + "    --set DIANN=/opt/diann/diann-linux --set FRAGPIPE=/opt/fragpipe/bin/fragpipe --print\n"
                        + "\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar resources-set --unset BASE_PATH --dry-run --print\n",
                true
        );
    }

    private static void printMagFragpipeHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar mag-fragpipe [options]",
                "\nRun the MAG+FragPipe workflow in headless mode (HAP search -> genome selection -> sample-specific search).",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar mag-fragpipe \\\n"
                        + "    --params /path/last_used_mag_fragpipe_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n"
                        + "\n"
                        + "Optional override:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar mag-fragpipe \\\n"
                        + "    --params /path/last_used_mag_fragpipe_parameters.json \\\n"
                        + "    --workflow /path/Basic-Search.workflow\n",
                true
        );
    }

    private static void printDiaMarkerHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar dia-marker [options]",
                "\nRun DIA Marker in headless mode from a GUI-exported parameter file. Supported marker modes: Core/Focus, Comprehensive, Deep Together, and Deep Per File. DIA Seed modes are not supported by this command.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar dia-marker \\\n"
                        + "    --params /path/last_used_dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n",
                true
        );
    }

    private static java.nio.file.Path resolveResourcesConfigPath(String explicit) {
        return HpcSiteConfig.resolveResourcesConfigPath(explicit);
    }

    private static LinkedHashMap<String, String> readKeyValue(Path file) throws IOException {
        LinkedHashMap<String, String> out = new LinkedHashMap<>();
        try (BufferedReader br = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                int eq = line.indexOf('=');
                if (eq < 1) continue;
                String key = line.substring(0, eq).trim();
                String value = line.substring(eq + 1).trim();
                if (!key.isEmpty()) out.put(key, value);
            }
        }
        return out;
    }

    private static List<String> renderParamsFile(LinkedHashMap<String, String> kv) {
        ArrayList<String> out = new ArrayList<>();
        out.add("# MetaPilot resource configuration");
        out.add("# Updated by metapilot-cli resources-set at " + LocalDateTime.now());
        out.add("");

        ArrayList<String> keys = new ArrayList<>(kv.keySet());
        keys.sort(Comparator.comparingInt(MetaPilotCliApp::resourceKeyOrder).thenComparing(String::compareToIgnoreCase));
        for (String k : keys) {
            String v = kv.get(k);
            if (v == null) v = "";
            out.add(k + "=" + v);
        }
        return out;
    }

    private static int resourceKeyOrder(String k) {
        if (k == null) return 1000;
        return switch (k) {
            case "SYSTEM" -> 0;
            case "RESOURCE_AUTO_DETECTED" -> 1;
            case "BASE_PATH" -> 2;
            default -> 10;
        };
    }

    private static LinkedHashMap<String, String> defaultLinuxResourceTemplate() {
        LinkedHashMap<String, String> out = new LinkedHashMap<>();
        out.put("SYSTEM", "linux");
        out.put("RESOURCE_AUTO_DETECTED", "true");
        out.put("BASE_PATH", "~/tools");

        // Keep these as simple defaults; users can override with resources-set --set KEY=VALUE.
        out.put("JAVA_MODE", "which");
        out.put("JAVA", "java");
        out.put("PYTHON_MODE", "which");
        out.put("PYTHON", "python3,python");
        out.put("GENERATE_REPORT", "true");

        // Core tools needed by current CLI workflows (DIA Marker + MAG FragPipe)
        out.put("DIANN_MODE", "path");
        out.put("DIANN", "~/tools/diann-2.5.0/diann-linux");

        out.put("FRAGPIPE_MODE", "path");
        out.put("FRAGPIPE", "~/tools/fragpipe-24.0/bin/fragpipe");
        out.put("FRAGPIPE_TOOLS", "~/tools/fragpipe-24.0/tools");

        // Optional keys (kept for completeness / future workflows)
        out.put("DEEPDETECT_MODE", "path");
        out.put("DEEPDETECT", "~/tools/deepdetect/deepdetect");
        out.put("FUNCDB", "~/tools/func_def.db");

        return out;
    }

    private static java.nio.file.Path resolveMagDbConfigPath(String explicit) {
        return HpcSiteConfig.resolveMagDbConfigPath(explicit);
    }

    private static Options buildMagDownloadOptions() {
        Options opts = new Options();

        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());

        opts.addOption(Option.builder()
                .longOpt("catalogue")
                .hasArg()
                .argName("ID")
                .desc("MGYG catalogue ID (e.g. human-gut)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("version")
                .hasArg()
                .argName("VER")
                .desc("Catalogue version (e.g. v2.0.2)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("out")
                .hasArg()
                .argName("DIR")
                .desc("Output directory (default: <jar-folder>/magdb; creates <out>/<catalogue>/<version>/...)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("resume")
                .desc("Resume downloads if partial files exist")
                .build());

        opts.addOption(Option.builder()
                .longOpt("threads")
                .hasArg()
                .argName("N")
                .desc("Download thread count (default: CPU/2)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("wait-hours")
                .hasArg()
                .argName("H")
                .desc("Maximum hours to wait for downloads (default: 6)")
                .build());

        return opts;
    }

    private static void printMagDownloadHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar mag-download [options]",
                "\nDownload one MGYG MAG catalogue+version and generate Marker FASTA files.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar mag-download \\\n"
                        + "    --catalogue human-gut --version v2.0.2 --out /data/magdb --resume\n"
                        + "\n"
                        + "  # Default output (next to jar): <jar-folder>/magdb/<catalogue>/<version>/...\n"
                        + "  java -jar metapilot-cli-1.0.0-shaded.jar mag-download \\\n"
                        + "    --catalogue human-gut --version v2.0.2 --resume\n",
                true
        );
    }


    private static String getVersion() {
        String v = Optional.ofNullable(MetaPilotCliApp.class.getPackage())
                .map(Package::getImplementationVersion)
                .orElse(null);
        return (v != null && !v.isBlank()) ? v : "dev";
    }

    private static int parseIntOr(String s, int fallback) {
        if (s == null) return fallback;
        try {
            return Integer.parseInt(s.trim());
        } catch (Exception e) {
            return fallback;
        }
    }

    private static int parsePositiveIntOption(CommandLine cmd, String primary, String alias) {
        int value = parseIntOr(cmd.getOptionValue(primary), 0);
        if (value <= 0 && alias != null) {
            value = parseIntOr(cmd.getOptionValue(alias), 0);
        }
        return Math.max(0, value);
    }

    private static void reportCallbackError(String context, Exception error) {
        System.err.println("ERROR [" + context + "]: " + (error != null ? error.getMessage() : "unknown"));
    }

    private static String describeFile(java.io.File file) {
        return file == null ? "<not configured>" : file.getAbsolutePath();
    }

    private static boolean reportGenerationEnabled(
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc,
            boolean cliDisabled
    ) {
        if (cliDisabled) {
            return false;
        }
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
}


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
        if (args != null && args.length > 0 && "dia-marker-deep-entrapment".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDiaMarkerDeepEntrapment(subArgs);
        }
        if (args != null && args.length > 0 && "dia-seed".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDiaSeed(subArgs);
        }
        if (args != null && args.length > 0 && "dia-library".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDiaLibrary(subArgs);
        }
        if (args != null && args.length > 0 && "mag-fragpipe".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runMagFragpipe(subArgs);
        }
        if (args != null && args.length > 0 && "mag-opensearch".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runMagOpenSearch(subArgs);
        }
        if (args != null && args.length > 0 && "mag-opensearch-entrapment".equalsIgnoreCase(args[0])) {
            return runMagOpenSearch(java.util.Arrays.copyOfRange(args, 1, args.length), true);
        }
        if (args != null && args.length > 0 && "mag-opensearch-report".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runMagOpenSearchReport(subArgs);
        }
        if (args != null && args.length > 0 && "pfind-mag".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runPfindMag(subArgs);
        }
        if (args != null && args.length > 0 && "pfind-mag-report".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runPfindMagReport(subArgs);
        }
        if (args != null && args.length > 0 && "diann-mag-report".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDiaMagReport(subArgs);
        }
        if (args != null && args.length > 0 && "diann-mag-peptide-map".equalsIgnoreCase(args[0])) {
            String[] subArgs = new String[args.length - 1];
            System.arraycopy(args, 1, subArgs, 0, subArgs.length);
            return runDiaMagPeptideMap(subArgs);
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
        options.addOption(Option.builder()
                .longOpt("continue-with-successful-files")
                .desc("Deep per file: when some raw files have no genome selection after the marker search, "
                        + "continue with the files that have it (as the GUI's continue option) instead of stopping. "
                        + "The skipped files get no per-file rounds but stay in the final combined search.")
                .build());
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
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
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
            if (cmd.hasOption("continue-with-successful-files")) {
                System.setProperty(uk.ac.quadram.metapilot.core.runner.DiaMarkerDeepPerFileHeadlessRunner
                        .CONTINUE_WITH_SUCCESSFUL_PROPERTY, "true");
                System.out.println("Deep per file: --continue-with-successful-files is on");
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
                    rc,
                    cmd.hasOption("export-metaatlas")
            );
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: dia-marker failed for params file: " + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    /**
     * DIA Marker Deep together with paired shuffled entrapment genomes in every round (evaluation mode; see
     * DiaMarkerDeepEntrapmentHeadlessRunner). The other workflows are not affected.
     */
    private static int runDiaMarkerDeepEntrapment(String[] args) {
        Options options = buildDiaMarkerOptions();
        options.addOption(Option.builder()
                .longOpt("entrapment-seed")
                .hasArg()
                .argName("N")
                .desc("Seed of the shuffled entrapment proteins (default 1)")
                .build());
        options.addOption(Option.builder()
                .longOpt("entrapment-id-digit")
                .hasArg()
                .argName("D")
                .desc("Entrapment genome IDs MGYG<D>xxxxxxxx (default 8; must differ from 0 and from any shadow "
                        + "genomes already in the catalogue, e.g. 9 for human-gut-shadow-s1)")
                .build());
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDiaMarkerDeepEntrapmentHelp(options);
            return 2;
        }
        if (cmd.hasOption("help")) {
            printDiaMarkerDeepEntrapmentHelp(options);
            return 0;
        }
        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printDiaMarkerDeepEntrapmentHelp(options);
            return 2;
        }
        long seed;
        int digit;
        try {
            seed = Long.parseLong(cmd.getOptionValue("entrapment-seed", "1"));
            digit = Integer.parseInt(cmd.getOptionValue("entrapment-id-digit",
                    String.valueOf(uk.ac.quadram.metapilot.core.dia.EntrapmentShuffle.DEFAULT_ID_DIGIT)));
        } catch (NumberFormatException e) {
            System.err.println("ERROR: --entrapment-seed and --entrapment-id-digit must be integers");
            return 2;
        }
        if (digit < 1 || digit > 9) {
            System.err.println("ERROR: --entrapment-id-digit must be 1-9");
            return 2;
        }
        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));
        try {
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);
            java.io.File diannExe = rc.getRaw("DIANN").map(java.io.File::new).orElse(null);
            if (diannExe == null || !diannExe.exists()) {
                System.err.println("ERROR: DIANN is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(diannExe) + ")");
                return 2;
            }
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
            if (p.getMarkerRefinementMode() != uk.ac.quadram.metapilot.task.dia.DiaRefinementMode.DEEP_TOGETHER) {
                System.out.println("NOTE: marker_mode in the params file is " + p.getMarkerRefinementMode().getPersistedValue()
                        + "; dia-marker-deep-entrapment always runs Deep together with entrapment");
            }
            System.out.println("DIA Marker Deep entrapment: shuffled copies MGYG" + digit + "..., seed " + seed);
            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("dia-marker-deep-entrapment", error); }
            };
            uk.ac.quadram.metapilot.core.runner.DiaMarkerDeepEntrapmentHeadlessRunner.run(p, diannExe, magdbCfg.toFile(),
                    cb, condaPrefix, funcDb, rc, seed, digit);
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: dia-marker-deep-entrapment failed for params file: "
                    + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static void printDiaMarkerDeepEntrapmentHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar dia-marker-deep-entrapment [options]",
                "\nEvaluation mode: DIA Marker Deep together with entrapment. Marker search and genome ranking as in "
                        + "Deep together; each round also searches a shuffled copy (first and last residue of every "
                        + "K/R segment fixed) of each of its target genomes, and entrapment peptides are carried "
                        + "forward like target peptides. The stop rules count target peptides only. Output: "
                        + "<output_directory>/marker_deep_entrapment/ with entrapment_summary.tsv, "
                        + "entrapment_peptides.tsv and marker/entrapment_rounds.tsv. The Deep together GUI has the "
                        + "same evaluation switch (off by default). Other workflows are unchanged.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-<version>-shaded.jar dia-marker-deep-entrapment \\\n"
                        + "    --params /path/dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n",
                true
        );
    }

    private static int runDiaSeed(String[] args) {
        Options options = buildDiaSeedOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDiaSeedHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printDiaSeedHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printDiaSeedHelp(options);
            return 2;
        }
        if (cmd.hasOption("host-lib") && !cmd.hasOption("seed-lib")) {
            System.err.println("ERROR: --host-lib requires --seed-lib.");
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            String workflow = readWorkflowField(java.nio.file.Paths.get(paramsFile));
            if (!uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType().equals(workflow)) {
                System.err.println("ERROR: dia-seed needs a DIA Seed parameter file (workflow \""
                        + uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType()
                        + "\"), got \"" + workflow + "\": " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }

            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File diannExe = rc.getRaw("DIANN").map(java.io.File::new).orElse(null);
            if (diannExe == null || !diannExe.exists()) {
                System.err.println("ERROR: DIANN is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(diannExe) + ")");
                return 2;
            }

            java.io.File pythonExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.PYTHON, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            java.io.File condaPrefix = uk.ac.quadram.metapilot.core.report.PythonEnvironmentResolver
                    .inferCondaPrefixFromPython(pythonExe);
            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile).orElse(null);

            uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia seedParams =
                    new uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (seedParams == null) {
                System.err.println("ERROR: failed to load DIA Seed params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }
            if (cmd.hasOption("seed-lib")) {
                seedParams.setLibrary(cmd.hasOption("host-lib")
                        ? new String[] { cmd.getOptionValue("seed-lib"), cmd.getOptionValue("host-lib") }
                        : new String[] { cmd.getOptionValue("seed-lib") });
            }
            if (seedParams.getMagDbItem() == null) {
                System.err.println("ERROR: MAG database in the params file (mag_database) could not be resolved from: "
                        + magdbCfg.toAbsolutePath());
                return 2;
            }

            uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker p;
            try {
                p = uk.ac.quadram.metapilot.core.runner.DiaSeedHeadlessRunner.toRunParameters(seedParams);
                uk.ac.quadram.metapilot.core.runner.DiaSeedHeadlessRunner.resolveLibraries(p.getLibrary(), p.getMagDbItem());
            } catch (java.io.IOException e) {
                System.err.println("ERROR: " + e.getMessage());
                return 2;
            }
            uk.ac.quadram.metapilot.task.dia.DiaRefinementMode seedMode = p.getMarkerRefinementMode();
            System.out.println("DIA Seed mode: " + seedMode.getPersistedValue() + " (" + seedMode.getDisplayLabel() + ")");

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("dia-seed", error); }
            };

            uk.ac.quadram.metapilot.core.runner.DiaSeedHeadlessRunner.run(
                    p,
                    diannExe,
                    magdbCfg.toFile(),
                    cb,
                    condaPrefix,
                    funcDb,
                    reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                    null,
                    rc,
                    cmd.hasOption("export-metaatlas")
            );
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: dia-seed failed for params file: " + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static Options buildDiaLibraryOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h").longOpt("help").desc("Show help").build());
        opts.addOption(Option.builder().longOpt("params").hasArg().argName("FILE").required()
                .desc("DIA parameter JSON (Seed or Marker): MAG database, threads and modifications for the library").build());
        opts.addOption(Option.builder().longOpt("name").hasArg().argName("NAME").required()
                .desc("Library name: folder libraries/<NAME> of the MAG database").build());
        opts.addOption(Option.builder().longOpt("peptides").hasArg().argName("FILE")
                .desc("Peptide table (column 'Base sequence' or 'Sequence', e.g. final_peptides.tsv) or peptide FASTA").build());
        opts.addOption(Option.builder().longOpt("fasta").hasArg().argName("FILE")
                .desc("Protein FASTA (for example a host proteome)").build());
        opts.addOption(Option.builder().longOpt("missed-cleavages").hasArg().argName("N")
                .desc("Missed cleavages (default 1; for a peptide list raised to keep its peptides intact, at most 2)").build());
        opts.addOption(Option.builder().longOpt("standalone")
                .desc("Standalone library (not mapped to MAG proteins, e.g. host); default: Microbiome").build());
        opts.addOption(Option.builder().longOpt("overwrite").desc("Rebuild an existing library of this name").build());
        opts.addOption(Option.builder().longOpt("resources").hasArg().argName("FILE")
                .desc("Site resources override (default: <install>/config/meta_resources.params)").build());
        opts.addOption(Option.builder().longOpt("magdb").hasArg().argName("FILE")
                .desc("MAG registry override (default: <install>/config/MAGDB_1_0.json)").build());
        return opts;
    }

    private static void printDiaLibraryHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar dia-library [options]",
                "\nBuild a DIA spectral library into the MAG database (libraries/<name>/: .fasta, .speclib, .parquet, "
                        + "_info.txt), with the same steps as the GUI library builder. From a peptide list, every peptide "
                        + "is one FASTA entry: protein N-terminal acetylation and N-terminal Met excision are not applied, "
                        + "and the missed cleavages keep the peptides intact. Give --peptides or --fasta.",
                options,
                "\nExamples:\n"
                        + "  dia-library --params dia_seed_parameters.json --name Elo_12mix \\\n"
                        + "    --peptides /path/pfind_MAG_search_marker4_focus/final_peptides.tsv\n"
                        + "  dia-library --params dia_marker_parameters.json --name host_human --standalone \\\n"
                        + "    --fasta /path/uniprot_human_reviewed.fasta\n",
                true);
    }

    private static int runDiaLibrary(String[] args) {
        Options options = buildDiaLibraryOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDiaLibraryHelp(options);
            return 2;
        }
        if (cmd.hasOption("help")) {
            printDiaLibraryHelp(options);
            return 0;
        }
        if (cmd.hasOption("peptides") == cmd.hasOption("fasta")) {
            System.err.println("ERROR: give exactly one of --peptides or --fasta.");
            printDiaLibraryHelp(options);
            return 2;
        }
        int missed;
        try {
            missed = cmd.hasOption("missed-cleavages") ? Integer.parseInt(cmd.getOptionValue("missed-cleavages")) : 1;
        } catch (NumberFormatException e) {
            System.err.println("ERROR: --missed-cleavages must be an integer");
            return 2;
        }
        String paramsFile = cmd.getOptionValue("params");
        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));
        try {
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);
            java.io.File diannExe = rc.getRaw("DIANN").map(java.io.File::new).orElse(null);
            if (diannExe == null || !diannExe.exists()) {
                System.err.println("ERROR: DIANN is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(diannExe) + ")");
                return 2;
            }
            java.io.File pythonExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.PYTHON, null)
                    .map(java.nio.file.Path::toFile).orElse(null);
            java.io.File condaPrefix = uk.ac.quadram.metapilot.core.report.PythonEnvironmentResolver
                    .inferCondaPrefixFromPython(pythonExe);

            String workflow = readWorkflowField(java.nio.file.Paths.get(paramsFile));
            uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia params =
                    uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType().equals(workflow)
                            ? new uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO().loadParameters(new java.io.File(paramsFile))
                            : new uk.ac.quadram.metapilot.task.dia.par.DiaMarkerParameterIO().loadParameters(new java.io.File(paramsFile));
            if (params == null) {
                System.err.println("ERROR: failed to load DIA params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }
            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("dia-library", error); }
            };
            boolean fromPeptides = cmd.hasOption("peptides");
            uk.ac.quadram.metapilot.core.runner.DiaLibraryHeadlessRunner.build(
                    params, magdbCfg.toFile(), diannExe, cmd.getOptionValue("name"),
                    new java.io.File(fromPeptides ? cmd.getOptionValue("peptides") : cmd.getOptionValue("fasta")),
                    fromPeptides ? uk.ac.quadram.metapilot.core.runner.DiaLibraryHeadlessRunner.Source.PEPTIDES
                            : uk.ac.quadram.metapilot.core.runner.DiaLibraryHeadlessRunner.Source.PROTEINS,
                    missed, !cmd.hasOption("standalone"), cmd.hasOption("overwrite"), cb, condaPrefix);
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: dia-library failed: " + e.getMessage());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static String readWorkflowField(java.nio.file.Path paramsFile) throws IOException {
        String text = new String(Files.readAllBytes(paramsFile), StandardCharsets.UTF_8);
        return new org.json.JSONObject(text).optString("workflow", "").trim();
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
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printMagFragpipeHelp(options);
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
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
            boolean identificationOnly = cmd.hasOption("identification-only");
            boolean exportMetaAtlas = cmd.hasOption("export-metaatlas");

            // Workflow template selection is entirely parameter-file driven (search_mode,
            // open_search_marker_stage, custom_basic_workflow_path,
            // custom_combined_workflow_path) - there is no --workflow CLI flag, so the
            // parameter file is the single source of truth for both GUI and CLI/HPC runs.
            // Pre-flight-validate whichever bundled templates the loaded params will
            // actually need, so a missing FragPipe workflow fails fast with a clear message
            // instead of partway through Step A.
            try {
                if (p.getCustomBasicWorkflowPath() == null || p.getCustomBasicWorkflowPath().isBlank()) {
                    String basicTemplate = p.isOpenSearchAtMarkerStage()
                            ? uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.OPEN_SEARCH_WORKFLOW
                            : uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.BASIC_SEARCH_WORKFLOW;
                    uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.requireInstalledWorkflow(fragpipeExe, basicTemplate);
                }
                // Identification-only runs never reach the combined/quant search, so the
                // combined-stage template doesn't need to be installed for this mode.
                if (!identificationOnly
                        && (p.getCustomCombinedWorkflowPath() == null || p.getCustomCombinedWorkflowPath().isBlank())) {
                    // Open mode's combined stage is Open.workflow with its quant keys
                    // switched on (see MagFragPipeHeadlessRunner.buildOpenLfqHybridTemplate) -
                    // confirmed against a FragPipe-authored reference template, it does not
                    // reference LFQ-MBR.workflow at all, so only Open.workflow is required here.
                    String combinedTemplate = p.isOpenSearch()
                            ? uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.OPEN_SEARCH_WORKFLOW
                            : uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.LFQ_MBR_WORKFLOW;
                    uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.requireInstalledWorkflow(fragpipeExe, combinedTemplate);
                }
            } catch (java.io.IOException e) {
                System.err.println("ERROR: required FragPipe workflows are not available in the installed FragPipe folder: "
                        + e.getMessage());
                return 2;
            }

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("mag-fragpipe", error); }
            };

            if (identificationOnly) {
                uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.runIdentificationOnlyWithInstalledDefaults(
                        p, fragpipeExe, magdbCfg.toFile(), cb);
            } else {
                uk.ac.quadram.metapilot.core.runner.MagFragPipeHeadlessRunner.runWithInstalledDefaults(
                        p,
                        fragpipeExe,
                        magdbCfg.toFile(),
                        funcDb,
                        cb,
                        reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                        null,
                        rc,
                        exportMetaAtlas
                );
            }
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: mag-fragpipe failed for params=" + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    /**
     * Headless MAG open-search workflow (closed marker search, then closed+open
     * on the sample-specific database via standalone MSFragger from FragPipe).
     * Same runner as the Windows GUI; the GUI menu entry may be hidden in a
     * release build, but this command is available for HPC testing.
     */
    private static int runMagOpenSearch(String[] args) {
        return runMagOpenSearch(args, false);
    }

    private static int runMagOpenSearch(String[] args, boolean entrapmentMode) {
        Options options = buildMagOpenSearchOptions();
        if (entrapmentMode) {
            options.addOption(Option.builder().longOpt("entrapment-seed").hasArg().argName("N")
                    .desc("Reproducible shuffled-protein seed (default 1)").build());
            options.addOption(Option.builder().longOpt("entrapment-id-digit").hasArg().argName("1-9")
                    .desc("MGYG digit for shuffled protein IDs (default 8)").build());
        }
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printMagOpenSearchHelp(options, entrapmentMode);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printMagOpenSearchHelp(options, entrapmentMode);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printMagOpenSearchHelp(options, entrapmentMode);
            return 2;
        }

        uk.ac.quadram.metapilot.core.opensearch.OpenSearchEntrapment.Config entrapment = null;
        if (entrapmentMode) {
            try {
                entrapment = new uk.ac.quadram.metapilot.core.opensearch.OpenSearchEntrapment.Config(
                        Long.parseLong(cmd.getOptionValue("entrapment-seed", "1")),
                        Integer.parseInt(cmd.getOptionValue("entrapment-id-digit", "8")));
            } catch (IllegalArgumentException e) {
                System.err.println("ERROR: invalid entrapment seed or ID digit: " + e.getMessage());
                return 2;
            }
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters p =
                    new uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load MAG open-search params from: "
                        + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }

            uk.ac.quadram.metapilot.task.mag.MagDbItem mag = p.getMagDbItem();
            if (mag == null) {
                mag = uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.resolve(
                        p.getMagDbCatalogueId(), p.getMagDbIdentifier(), p.getMagDbVersion(), magdbCfg.toFile());
                if (mag != null) {
                    p.setMagDbItem(mag);
                }
            }
            if (mag == null) {
                System.err.println("ERROR: MAG database entry could not be resolved from site registry: "
                        + magdbCfg.toAbsolutePath()
                        + " (catalogue=" + p.getMagDbCatalogueId()
                        + ", identifier=" + p.getMagDbIdentifier()
                        + ", version=" + p.getMagDbVersion() + ")");
                return 2;
            }

            uk.ac.quadram.metapilot.resources.MetaResourceParameters mrp =
                    uk.ac.quadram.metapilot.core.runner.MagOpenSearchToolSetup.fromResourceConfig(rc);
            uk.ac.quadram.metapilot.core.runner.MagOpenSearchToolSetup.Tools tools =
                    (entrapmentMode || p.isEntrapmentEnabled())
                            ? uk.ac.quadram.metapilot.core.runner.MagOpenSearchToolSetup.prepareIdentification(p, mag, mrp)
                            : uk.ac.quadram.metapilot.core.runner.MagOpenSearchToolSetup.prepare(p, mag, mrp);

            if (cmd.hasOption("continue-with-successful-files")) {
                System.setProperty(uk.ac.quadram.metapilot.core.runner.MagOpenSearchHeadlessRunner
                        .CONTINUE_WITH_SUCCESSFUL_PROPERTY, "true");
                System.out.println("MAG open search: --continue-with-successful-files is on");
            }

            System.out.println("MAG open search: engines=" + p.getSearchEngines()
                    + ", round1=" + p.getRound1SearchMode()
                    + ", closed_stage=" + p.isClosedStageEnabled()
                    + ", MSFragger=" + tools.msfraggerJar.getAbsolutePath());

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb =
                    new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) {
                    if (message != null && !message.isBlank()) System.out.println(message);
                }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("mag-opensearch", error); }
            };

            if (entrapmentMode) {
                p.setEntrapmentEnabled(true);
                p.setEntrapmentSeed(entrapment.seed);
                p.setEntrapmentIdDigit(entrapment.idDigit);
                System.out.println("Entrapment evaluation: joint plus closed-only, open-only and cascade arms; "
                        + "shuffles added after genome selection; stops after identification.");
                uk.ac.quadram.metapilot.core.runner.MagOpenSearchHeadlessRunner.runEntrapment(
                        p, tools.msfraggerJar, tools.markerFasta, mag, cb, entrapment);
            } else {
                uk.ac.quadram.metapilot.core.runner.MagOpenSearchHeadlessRunner.run(
                    p,
                    tools.msfraggerJar,
                    tools.markerFasta,
                    mag,
                    tools.flashLfqExe,
                    p.getThreadCount(),
                    tools.funcDbFile,
                    cb);
            }
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: mag-opensearch failed for params="
                    + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    /**
     * Headless-only: rebuild open-search final tables from an already-produced
     * {@code validated_psms.tsv} and {@code QuantifiedPeptides.tsv}. No MSFragger,
     * Crystal-C, classification or FlashLFQ. Applies MagResultFinalisationService
     * and the current export columns (shift explanation, spectrum count, mzML
     * identification rate when mzML is readable).
     */
    private static int runMagOpenSearchReport(String[] args) {
        Options options = buildMagOpenSearchReportOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printMagOpenSearchReportHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printMagOpenSearchReportHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printMagOpenSearchReportHelp(options);
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);

            uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters p =
                    new uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load MAG open-search params from: "
                        + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }
            String out = cmd.getOptionValue("out");
            if (out != null && !out.isBlank()) {
                p.setResult(out);
            }

            uk.ac.quadram.metapilot.task.mag.MagDbItem mag = p.getMagDbItem();
            if (mag == null) {
                mag = uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.resolve(
                        p.getMagDbCatalogueId(), p.getMagDbIdentifier(), p.getMagDbVersion(), magdbCfg.toFile());
                if (mag != null) {
                    p.setMagDbItem(mag);
                }
            }
            if (mag == null) {
                System.err.println("ERROR: MAG database entry could not be resolved from site registry: "
                        + magdbCfg.toAbsolutePath()
                        + " (catalogue=" + p.getMagDbCatalogueId()
                        + ", identifier=" + p.getMagDbIdentifier()
                        + ", version=" + p.getMagDbVersion() + ")");
                return 2;
            }

            java.io.File validatedPsms = null;
            String validatedOpt = cmd.getOptionValue("validated-psms");
            if (validatedOpt != null && !validatedOpt.isBlank()) {
                validatedPsms = new java.io.File(validatedOpt);
            }
            java.io.File quantified = null;
            String quantifiedOpt = cmd.getOptionValue("quantified-peptides");
            if (quantifiedOpt != null && !quantifiedOpt.isBlank()) {
                quantified = new java.io.File(quantifiedOpt);
            }

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb =
                    new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) {
                    if (message != null && !message.isBlank()) System.out.println(message);
                }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("mag-opensearch-report", error); }
            };

            uk.ac.quadram.metapilot.core.runner.MagOpenSearchHeadlessRunner.runPostAnalysisOnly(
                    p, mag, funcDb, validatedPsms, quantified, cb);
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: mag-opensearch-report failed for params="
                    + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static void applyPfindGenomeSelectionOverrides(CommandLine cmd,
            uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters params) {
        if (cmd.hasOption("genome-selection-strategy")) {
            params.setGenomeSelectionStrategy(uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut
                    .requireStrategy(cmd.getOptionValue("genome-selection-strategy")));
        }
        if (cmd.hasOption("psm-coverage")) {
            params.setGenomeSelectionPsmCoverage(parseCliDouble("psm-coverage", cmd.getOptionValue("psm-coverage")));
        }
        if (cmd.hasOption("step-guard")) {
            params.setGenomeSelectionStepGuard(parseCliDouble("step-guard", cmd.getOptionValue("step-guard")));
        }
        if (cmd.hasOption("max-genomes")) {
            params.setGenomeSelectionMaxGenomes(parseCliInt("max-genomes", cmd.getOptionValue("max-genomes")));
        }
        boolean coverageFlag = cmd.hasOption("psm-coverage") || cmd.hasOption("step-guard") || cmd.hasOption("max-genomes");
        if (coverageFlag && !uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut
                .isPsmCoverage(params.getGenomeSelectionStrategy())) {
            System.err.println("WARNING: --psm-coverage, --step-guard and --max-genomes are used only when "
                    + "--genome-selection-strategy is \"PSM coverage\". The current strategy is "
                    + params.getGenomeSelectionStrategy() + ", so those values are saved and ignored.");
        }
    }

    private static double parseCliDouble(String name, String text) {
        try {
            return Double.parseDouble(text.trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("--" + name + " must be a number (got " + text + ")");
        }
    }

    private static int parseCliInt(String name, String text) {
        try {
            return Integer.parseInt(text.trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("--" + name + " must be an integer (got " + text + ")");
        }
    }

    /**
     * Windows only: the whole DDA MAG pFind workflow from a parameter file -
     * spectra conversion, marker search, genome selection, sample-specific
     * search, FlashLFQ (and isobaric) quantification, post-analysis, taxonomy,
     * functional annotation and the HTML report. See
     * {@link uk.ac.quadram.metapilot.core.runner.PfindMagSearchHeadlessRunner}.
     */
    private static int runPfindMag(String[] args) {
        Options options = buildPfindMagOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printPfindMagHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printPfindMagHelp(options);
            return 0;
        }

        if (!uk.ac.quadram.metapilot.core.runner.PfindMagSearchHeadlessRunner.isSupportedPlatform()) {
            System.err.println("ERROR: pfind-mag runs on Windows only (pFind and FlashLFQ are Windows executables); "
                    + "this is " + System.getProperty("os.name") + ". To post-process pFind results produced on "
                    + "Windows, use pfind-mag-report.");
            return 2;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printPfindMagHelp(options);
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.task.mag.MagDbItemResolver.setDefaultConfigFile(magdbCfg.toFile());
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File pFindExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.PFIND, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);
            java.io.File flashLfqExe = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FLASHLFQ, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);
            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);
            if (pFindExe == null || !pFindExe.isFile()) {
                System.err.println("ERROR: PFIND is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(pFindExe) + ")");
                return 2;
            }
            if (flashLfqExe == null || !flashLfqExe.exists()) {
                System.err.println("ERROR: FLASHLFQ is not set or does not exist in resources config: "
                        + resourcesCfg.toAbsolutePath() + " (resolved path: " + describeFile(flashLfqExe)
                        + "). Set FLASHLFQ to CMD.exe, CMD.dll, or the FlashLFQ install directory.");
                return 2;
            }
            if (funcDb == null || !funcDb.isFile()) {
                System.out.println("WARNING: FUNCDB is not set or does not exist (" + describeFile(funcDb)
                        + "); functional annotation will be skipped.");
            }

            uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters p =
                    new uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load MAG pFind params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }
            String out = cmd.getOptionValue("out");
            if (out != null && !out.isBlank()) {
                p.setResult(new java.io.File(out).getAbsolutePath());
            }
            applyPfindGenomeSelectionOverrides(cmd, p);

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("pfind-mag", error); }
            };

            uk.ac.quadram.metapilot.core.runner.PfindMagSearchHeadlessRunner.Request request =
                    new uk.ac.quadram.metapilot.core.runner.PfindMagSearchHeadlessRunner.Request();
            request.params = p;
            request.pFindExe = pFindExe;
            request.flashLfqExe = flashLfqExe;
            request.magdbConfigJson = magdbCfg.toFile();
            request.funcDbFile = funcDb;
            request.parallelTasks = parsePositiveIntOption(cmd, "parallel-tasks", null);
            request.generateHtmlReport = reportGenerationEnabled(rc, cmd.hasOption("no-report"));
            request.reportResourceConfig = rc;
            request.exportMetaAtlasParquet = cmd.hasOption("export-metaatlas");
            request.cb = cb;
            uk.ac.quadram.metapilot.core.runner.PfindMagSearchHeadlessRunner.run(request);
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: pfind-mag failed for params=" + new java.io.File(paramsFile).getAbsolutePath()
                    + ": " + e.getMessage());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    /**
     * Headless-only: resume a DDA MAG pFind workflow from an already-produced
     * peptide quantification (QuantifiedPeptides.tsv for label-free,
     * IsobaricQuantifiedPeptides.tsv for isobaric) + psms.tsv and run
     * post-analysis, taxonomy, functional annotation, and the final HTML
     * report. pFind identification and quantification come from the Windows
     * GUI or from pfind-mag - this only continues from their output.
     */
    private static int runPfindMagReport(String[] args) {
        Options options = buildPfindMagReportOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printPfindMagReportHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printPfindMagReportHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printPfindMagReportHelp(options);
            return 2;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);

            uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters p =
                    new uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load MAG pFind params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }

            java.io.File quantifiedPeptides = null;
            String quantifiedPeptidesOpt = cmd.getOptionValue("quantified-peptides");
            if (quantifiedPeptidesOpt != null && !quantifiedPeptidesOpt.isBlank()) {
                quantifiedPeptides = new java.io.File(quantifiedPeptidesOpt);
            }
            java.io.File psms = null;
            String psmsOpt = cmd.getOptionValue("psms");
            if (psmsOpt != null && !psmsOpt.isBlank()) {
                psms = new java.io.File(psmsOpt);
            }

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("pfind-mag-report", error); }
            };

            uk.ac.quadram.metapilot.core.runner.PfindMagHeadlessRunner.run(
                    p,
                    magdbCfg.toFile(),
                    quantifiedPeptides,
                    psms,
                    funcDb,
                    cb,
                    reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                    null,
                    rc,
                    cmd.hasOption("export-metaatlas")
            );
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: pfind-mag-report failed for params=" + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runDiaMagReport(String[] args) {
        Options options = buildDiaMagReportOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDiaMagReportHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printDiaMagReportHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        if (paramsFile == null) {
            System.err.println("ERROR: --params is required.");
            printDiaMagReportHelp(options);
            return 2;
        }

        Double globalPrecursorQValue = null;
        String globalPrecursorQValueOpt = cmd.getOptionValue("global-precursor-qvalue");
        if (globalPrecursorQValueOpt != null) {
            double value;
            try {
                value = Double.parseDouble(globalPrecursorQValueOpt);
            } catch (NumberFormatException e) {
                value = Double.NaN;
            }
            if (!(value > 0.0 && value <= 0.01)) {
                System.err.println("ERROR: --global-precursor-qvalue must be greater than 0 and at most 0.01.");
                return 2;
            }
            globalPrecursorQValue = value;
        }

        java.nio.file.Path resourcesCfg = resolveResourcesConfigPath(cmd.getOptionValue("resources"));
        java.nio.file.Path magdbCfg = resolveMagDbConfigPath(cmd.getOptionValue("magdb"));

        try {
            uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig rc =
                    uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig.load(resourcesCfg);

            java.io.File funcDb = rc.resolveToolPath(uk.ac.quadram.metapilot.core.resources.ToolKey.FUNCDB, null)
                    .map(java.nio.file.Path::toFile)
                    .orElse(null);

            uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker p =
                    new uk.ac.quadram.metapilot.task.dia.par.DiaMarkerParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load DIA Marker params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }

            java.io.File resultDirOverride = null;
            String outOpt = cmd.getOptionValue("out");
            if (outOpt != null && !outOpt.isBlank()) {
                resultDirOverride = new java.io.File(outOpt);
            }
            java.io.File combinedParquet = null;
            String combinedParquetOpt = cmd.getOptionValue("combined-parquet");
            if (combinedParquetOpt != null && !combinedParquetOpt.isBlank()) {
                combinedParquet = new java.io.File(combinedParquetOpt);
            }
            java.io.File prMatrix = null;
            String prMatrixOpt = cmd.getOptionValue("pr-matrix");
            if (prMatrixOpt != null && !prMatrixOpt.isBlank()) {
                prMatrix = new java.io.File(prMatrixOpt);
            }
            java.io.File roundRoot = null;
            String roundRootOpt = cmd.getOptionValue("round-root");
            if (roundRootOpt != null && !roundRootOpt.isBlank()) {
                roundRoot = new java.io.File(roundRootOpt);
            }
            java.io.File peptideProteinMap = null;
            String peptideProteinMapOpt = cmd.getOptionValue("peptide-protein-map");
            if (peptideProteinMapOpt != null && !peptideProteinMapOpt.isBlank()) {
                peptideProteinMap = new java.io.File(peptideProteinMapOpt);
            }

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("diann-mag-report", error); }
            };

            uk.ac.quadram.metapilot.core.runner.DiaMagReportHeadlessRunner.run(
                    p,
                    magdbCfg.toFile(),
                    resultDirOverride,
                    combinedParquet,
                    prMatrix,
                    peptideProteinMap,
                    roundRoot,
                    funcDb,
                    cb,
                    reportGenerationEnabled(rc, cmd.hasOption("no-report")),
                    null,
                    rc,
                    cmd.hasOption("export-metaatlas"),
                    cmd.hasOption("allow-incomplete-deep-protein-mapping"),
                    cmd.getOptionValue("protein-mapping-policy", "SEARCH").toUpperCase(java.util.Locale.ROOT),
                    cmd.hasOption("selected-genomes") ? new java.io.File(cmd.getOptionValue("selected-genomes")) : null,
                    globalPrecursorQValue
            );
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: diann-mag-report failed for params=" + new java.io.File(paramsFile).getAbsolutePath());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    private static int runDiaMagPeptideMap(String[] args) {
        Options options = buildDiaMagPeptideMapOptions();
        CommandLineParser parser = new DefaultParser();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("ERROR: " + e.getMessage());
            printDiaMagPeptideMapHelp(options);
            return 2;
        }

        if (cmd.hasOption("help")) {
            printDiaMagPeptideMapHelp(options);
            return 0;
        }

        String paramsFile = cmd.getOptionValue("params");
        String outOpt = cmd.getOptionValue("out");
        if (paramsFile == null || outOpt == null) {
            System.err.println("ERROR: --params and --out are required.");
            printDiaMagPeptideMapHelp(options);
            return 2;
        }

        try {
            uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker p =
                    new uk.ac.quadram.metapilot.task.dia.par.DiaMarkerParameterIO()
                            .loadParameters(new java.io.File(paramsFile));
            if (p == null) {
                System.err.println("ERROR: failed to load DIA Marker params from: " + new java.io.File(paramsFile).getAbsolutePath());
                return 2;
            }

            java.io.File markerRootOverride = null;
            String markerRootOpt = cmd.getOptionValue("marker-root");
            if (markerRootOpt != null && !markerRootOpt.isBlank()) {
                markerRootOverride = new java.io.File(markerRootOpt);
            }
            java.io.File outTsv = new java.io.File(outOpt);

            uk.ac.quadram.metapilot.core.workflow.ProgressCallback cb = new uk.ac.quadram.metapilot.core.workflow.ProgressCallback() {
                @Override public void onProgress(double workDone, double max) { }
                @Override public void onMessage(String message) { if (message != null && !message.isBlank()) System.out.println(message); }
                @Override public void onComplete(Object result) { System.out.println("Done: " + result); }
                @Override public void onError(Exception error) { reportCallbackError("diann-mag-peptide-map", error); }
            };

            uk.ac.quadram.metapilot.core.runner.DiaMagPeptideProteinMapExportRunner.run(p, markerRootOverride, outTsv, cb);
            return 0;
        } catch (Exception e) {
            System.err.println("ERROR: diann-mag-peptide-map failed for params=" + new java.io.File(paramsFile).getAbsolutePath());
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

        boolean proteinCatalogueOnly = cmd.hasOption("protein-catalogue-only");
        boolean proteinCatalogue = proteinCatalogueOnly || cmd.hasOption("protein-catalogue");
        boolean deleteProteinArchive = !cmd.hasOption("keep-protein-archive");
        if (proteinCatalogue
                && !uk.ac.quadram.metapilot.core.variant.ProteinCatalogueInstaller.supports(catalogue)) {
            System.err.println("ERROR: the UHGP-100 protein catalogue is only available for the "
                    + uk.ac.quadram.metapilot.core.variant.ProteinCatalogueInstaller.SUPPORTED_CATALOGUE
                    + " catalogue.");
            return 2;
        }

        try {
            if (proteinCatalogueOnly) {
                uk.ac.quadram.metapilot.core.MGYG.MGYGMagDownloader.prepareProteinCatalogue(
                        catalogue,
                        version,
                        new java.io.File(out),
                        deleteProteinArchive,
                        cb
                );
                return 0;
            }
            uk.ac.quadram.metapilot.core.MGYG.MGYGMagDownloader.downloadCatalogue(
                    catalogue,
                    version,
                    new java.io.File(out),
                    resume,
                    threads,
                    waitHours,
                    proteinCatalogue,
                    deleteProteinArchive,
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
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar " + command
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
                .longOpt("identification-only")
                .desc("Headless-only: stop after identification (marker search, genome selection, sample-specific "
                        + "search). Skips the combined FASTA build, combined quantification search, and all "
                        + "post-analysis export (final_peptides.tsv, taxonomy, functional annotation, HTML report). "
                        + "Does not affect the Windows GUI.")
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
        opts.addOption(Option.builder()
                .longOpt("export-metaatlas")
                .desc("Opt-in: also write MetaAtlas-aligned taxon_abundance.parquet / function_abundance.parquet "
                        + "alongside the standard taxonomy/functional annotation reports. Off by default; not "
                        + "available in the Windows GUI.")
                .build());
        return opts;
    }

    private static Options buildMagOpenSearchOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("MAG open-search parameters JSON saved by the GUI (workflow DDA_MAG_OPENSEARCH); "
                        + "round-1 mode, closed-stage cascade, Crystal-C, MSBooster and enzyme settings are "
                        + "all read from it")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration with FRAGPIPE (bundled MSFragger / Crystal-C / MSBooster), "
                        + "FLASHLFQ and optional FUNCDB / SAGE / MSCONVERT")
                .build());
        opts.addOption(Option.builder()
                .longOpt("magdb")
                .hasArg()
                .argName("FILE")
                .desc("MAG DB config JSON (MAGDB_1_0.json). Default: <install>/config/MAGDB_1_0.json")
                .build());
        opts.addOption(Option.builder()
                .longOpt("continue-with-successful-files")
                .desc("When some samples have no stage-1 genome selection (empty unmodified/rescored set), "
                        + "continue stage-2 identification with the successful samples instead of stopping. "
                        + "Skipped samples get no sample-specific search but stay in FlashLFQ quantification "
                        + "(match-between-runs), matching DIA Deep per-file.")
                .build());
        return opts;
    }

    private static Options buildMagOpenSearchReportOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("MAG open-search parameters JSON. Identification and quantification are not re-run.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration (FUNCDB optional). Default resolution same as resources-check")
                .build());
        opts.addOption(Option.builder()
                .longOpt("magdb")
                .hasArg()
                .argName("FILE")
                .desc("MAG DB config JSON (MAGDB_1_0.json). Default: <install>/config/MAGDB_1_0.json")
                .build());
        opts.addOption(Option.builder()
                .longOpt("out")
                .hasArg()
                .argName("DIR")
                .desc("User result folder containing opensearch_MAG_search/. Overrides output_directory in --params. "
                        + "Required when the parameter file still points at an HPC path that is not present here.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("validated-psms")
                .hasArg()
                .argName("FILE")
                .desc("Override path to validated_psms.tsv. Default: "
                        + "<result>/opensearch_MAG_search/combined/validated_psms.tsv")
                .build());
        opts.addOption(Option.builder()
                .longOpt("quantified-peptides")
                .hasArg()
                .argName("FILE")
                .desc("Override path to QuantifiedPeptides.tsv. Default: the same combined/ folder")
                .build());
        return opts;
    }

    private static Options buildPfindMagOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("MAG pFind parameters JSON saved by the Windows GUI (e.g. last_used_mag_pfind_parameters.json): "
                        + "raw files, output directory, MAG DB, enzyme, modifications, MS2 scan mode, search mode, "
                        + "genome selection and quantification are all read from it")
                .build());
        opts.addOption(Option.builder()
                .longOpt("resources")
                .hasArg()
                .argName("FILE")
                .desc("Resource configuration file with PFIND, FLASHLFQ and FUNCDB (the Windows GUI's "
                        + "config\\meta_resources.params works as is). Default resolution same as resources-check")
                .build());
        opts.addOption(Option.builder()
                .longOpt("magdb")
                .hasArg()
                .argName("FILE")
                .desc("MAG DB config JSON (MAGDB_1_0.json). Default: <install>/config/MAGDB_1_0.json")
                .build());
        opts.addOption(Option.builder()
                .longOpt("out")
                .hasArg()
                .argName("DIR")
                .desc("Output directory; overrides output_directory in --params. Results go to <DIR>\\pfind_MAG_search "
                        + "and spectra to <DIR>\\spectra, as in the GUI")
                .build());
        opts.addOption(Option.builder()
                .longOpt("parallel-tasks")
                .hasArg()
                .argName("N")
                .desc("Raw files searched at the same time. Default: one per 5 GB of physical memory, as in the GUI")
                .build());
        opts.addOption(Option.builder()
                .longOpt("genome-selection-strategy")
                .hasArg()
                .argName("MODE")
                .desc("Genome-selection cut: Focus, Comprehensive, or \"PSM coverage\". "
                        + "Overrides genome_selection_strategy in --params. Default in the parameter file: Focus")
                .build());
        opts.addOption(Option.builder()
                .longOpt("psm-coverage")
                .hasArg()
                .argName("FRACTION")
                .desc("PSM coverage target, > 0 and <= 1. Used only when the strategy is \"PSM coverage\". Default 0.95")
                .build());
        opts.addOption(Option.builder()
                .longOpt("step-guard")
                .hasArg()
                .argName("FRACTION")
                .desc("Stop when further genomes add at most this share of marker PSMs. "
                        + "Fraction from 0 to 0.002 (0 to 0.2%); 0 switches it off. "
                        + "The GUI shows the same value as a percentage. "
                        + "Used only when the strategy is \"PSM coverage\". Default 0.001 (0.1%)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("max-genomes")
                .hasArg()
                .argName("N")
                .desc("Upper limit on genomes per raw file for \"PSM coverage\". Automatic stays at 300. Default 300")
                .build());
        opts.addOption(Option.builder()
                .longOpt("no-report")
                .desc("Disable final HTML report generation")
                .build());
        opts.addOption(Option.builder()
                .longOpt("export-metaatlas")
                .desc("Opt-in: also write MetaAtlas-aligned taxon_abundance.parquet / function_abundance.parquet "
                        + "alongside the standard taxonomy/functional annotation reports. Off by default; not "
                        + "available in the Windows GUI.")
                .build());
        return opts;
    }

    private static Options buildPfindMagReportOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("MAG pFind parameters JSON (e.g. last_used_mag_pfind_parameters.json). Only used for the "
                        + "MAG DB catalogue/version, raw file list, and output directory - identification/quant "
                        + "parameters are not re-run.")
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
                .longOpt("quantified-peptides")
                .hasArg()
                .argName("FILE")
                .desc("Override path to the peptide quantification. Default: "
                        + "<output_directory>/pfind_MAG_search/combined/QuantifiedPeptides.tsv for label-free, "
                        + "IsobaricQuantifiedPeptides.tsv for isobaric (same location the Windows GUI writes it to)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("psms")
                .hasArg()
                .argName("FILE")
                .desc("Override path to psms.tsv. Default: "
                        + "<output_directory>/pfind_MAG_search/combined/psms.tsv")
                .build());
        opts.addOption(Option.builder()
                .longOpt("no-report")
                .desc("Disable final HTML report generation")
                .build());
        opts.addOption(Option.builder()
                .longOpt("export-metaatlas")
                .desc("Opt-in: also write MetaAtlas-aligned taxon_abundance.parquet / function_abundance.parquet "
                        + "alongside the standard taxonomy/functional annotation reports. Off by default; not "
                        + "available in the Windows GUI.")
                .build());
        return opts;
    }

    private static Options buildDiaMagReportOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("DIA Marker parameters JSON (e.g. last_used_dia_marker_parameters.json). Only used for the "
                        + "MAG DB catalogue/version, marker mode, and output directory - DIA-NN search itself is "
                        + "not re-run and no raw DIA files are read.")
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
                .desc("MAG DB config JSON (MAGDB_1_0.json). Default: <install>/config/MAGDB_1_0.json. This must be "
                        + "a registry that resolves on THIS machine (e.g. produced by 'mag-register --repository "
                        + "<hpc-path-to-mag-db>') - not the Windows GUI's registry file.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("out")
                .hasArg()
                .argName("DIR")
                .desc("Output directory. Strongly recommended: the 'result' path inside --params is exported by "
                        + "the Windows GUI and points at a local path on that machine, which this command never "
                        + "otherwise needs. When --out is set, params.result is ignored entirely and results are "
                        + "written directly to this directory (no marker_focus/marker_comprehensive/marker_deep "
                        + "subfolder is appended).")
                .build());
        opts.addOption(Option.builder()
                .longOpt("combined-parquet")
                .hasArg()
                .argName("FILE")
                .desc("Override path to combined.parquet. Default: <output_directory>/combined/combined.parquet, "
                        + "where <output_directory> is --out if set, else derived from --params (Windows GUI / "
                        + "dia-marker layout: marker_focus, marker_comprehensive, or marker_deep for both Deep "
                        + "Together and Deep Per File)")
                .build());
        opts.addOption(Option.builder()
                .longOpt("pr-matrix")
                .hasArg()
                .argName("FILE")
                .desc("Override path to combined.pr_matrix.tsv. Default: same folder as combined.parquet")
                .build());
        opts.addOption(Option.builder()
                .longOpt("global-precursor-qvalue")
                .hasArg()
                .argName("Q")
                .desc("Additional inclusive DIA-NN Global.Q.Value cutoff for microbial precursor evidence "
                        + "before genome scoring (0 < Q <= 0.01, e.g. 0.001). Uses the existing search and "
                        + "quantification matrix.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("peptide-protein-map")
                .hasArg()
                .argName("FILE")
                .desc("Deep Together / Deep Per File marker modes only, and the RECOMMENDED input for them: a "
                        + "compact peptide-to-proteins TSV produced once, locally, by the diann-mag-peptide-map "
                        + "command (run where the round parquets already exist). Required for those two modes "
                        + "unless --round-root resolves or --allow-incomplete-deep-protein-mapping is set - see "
                        + "diann-mag-report --help description for why. Takes priority over --round-root when both "
                        + "are given. Ignored for Core/Focus/Comprehensive.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("round-root")
                .hasArg()
                .argName("DIR")
                .desc("Deep Together / Deep Per File marker modes only: fallback path to the actual 'marker' "
                        + "round-parquet folder (round1.parquet, round2.parquet, ... or, for Deep Per File, "
                        + "marker/<rawname>/round*/...), for when you'd rather copy that whole folder to this "
                        + "machine instead of using --peptide-protein-map. Only used if --peptide-protein-map is "
                        + "not given. Ignored for Core/Focus/Comprehensive.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("allow-incomplete-deep-protein-mapping")
                .desc("Deprecated: incomplete Deep mappings are rejected. Supply --peptide-protein-map or --round-root.")
                .build());
        opts.addOption(Option.builder().longOpt("protein-mapping-policy").hasArg()
                .desc("Explicit report mapping: SEARCH (default), ROUND_PROTEINS, SEED_DIRECT, or SEED_SELECTED. Recovery policies require --out.").build());
        opts.addOption(Option.builder().longOpt("selected-genomes").hasArg()
                .desc("Recorded genome_selected.tsv used to build a Seed Focus/Comprehensive library; required for SEED_SELECTED.").build());
        opts.addOption(Option.builder()
                .longOpt("no-report")
                .desc("Disable final HTML report generation")
                .build());
        opts.addOption(Option.builder()
                .longOpt("export-metaatlas")
                .desc("Opt-in: also write MetaAtlas-aligned taxon_abundance.parquet / function_abundance.parquet "
                        + "alongside the standard taxonomy/functional annotation reports. Off by default; not "
                        + "available in the Windows GUI.")
                .build());
        return opts;
    }

    private static Options buildDiaMagPeptideMapOptions() {
        Options opts = new Options();
        opts.addOption(Option.builder("h")
                .longOpt("help")
                .desc("Show help")
                .build());
        opts.addOption(Option.builder()
                .longOpt("params")
                .hasArg()
                .argName("FILE")
                .desc("DIA Marker parameters JSON for a completed Deep Together or Deep Per File local run (same "
                        + "file the GUI/dia-marker used). Must be run on the machine where that run's round "
                        + "parquets still exist - only used for marker mode and the output directory "
                        + "(<result>/marker_deep), nothing is re-run.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("marker-root")
                .hasArg()
                .argName("DIR")
                .desc("Override for the marker_deep folder (the one containing 'marker/' and 'combined/'). "
                        + "Default: <params result>/marker_deep.")
                .build());
        opts.addOption(Option.builder()
                .longOpt("out")
                .hasArg()
                .argName("FILE")
                .desc("Where to write the compact peptide-protein map TSV. This is the one file you then transfer "
                        + "to HPC and pass to 'diann-mag-report --peptide-protein-map'.")
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
        opts.addOption(Option.builder()
                .longOpt("export-metaatlas")
                .desc("Opt-in: also write MetaAtlas-aligned taxon_abundance.parquet / function_abundance.parquet "
                        + "alongside the standard taxonomy/functional annotation reports. Off by default; not "
                        + "available in the Windows GUI.")
                .build());
        return opts;
    }

    private static Options buildDiaSeedOptions() {
        Options opts = buildDiaMarkerOptions();
        opts.getOption("params").setDescription(
                "DIA Seed parameters JSON saved by the GUI (e.g. dia_seed_analysis_parameters.json)");
        opts.addOption(Option.builder()
                .longOpt("seed-lib")
                .hasArg()
                .argName("NAME")
                .desc("Seed library name, overriding dia_settings.libraries. It must be a Microbiome library "
                        + "created in the GUI under <MAG DB>/libraries/<NAME>")
                .build());
        opts.addOption(Option.builder()
                .longOpt("host-lib")
                .hasArg()
                .argName("NAME|FILE")
                .desc("Optional host library: a Standalone library name or a .parquet file (needs --seed-lib)")
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
                        + "  dia-marker        Run the DIA Marker workflow from a GUI parameter file\n"
                        + "  dia-marker-deep-entrapment  Evaluation: DIA Marker Deep together with shuffled entrapment genomes in every round\n"
                        + "  dia-seed          Run the DIA Seed workflow (Direct, Core/Focus, Comprehensive, Deep) from a GUI parameter file\n"
                        + "  dia-library       Build a DIA library (seed from a peptide list, or host from a protein FASTA) into the MAG DB\n"
                        + "  mag-fragpipe      Run the DDA MAG FragPipe workflow from a GUI parameter file\n"
                        + "  mag-opensearch    Run the DDA MAG open-search workflow (closed then closed+open) from a GUI parameter file\n"
                        + "  mag-opensearch-entrapment  Evaluate closed/open/cascade/joint with second-search paired shuffles\n"
                        + "  mag-opensearch-report  Rebuild open-search final tables from existing validated_psms.tsv (no search)\n"
                        + "  pfind-mag         Run the DDA MAG pFind workflow end to end (Windows only)\n"
                        + "  pfind-mag-report  Post-process pFind workflow quantification (any OS)\n"
                        + "  deepdetect-predict  Run DeepDetect PyTorch prediction (python main.py)\n"
                        + "  deepdetect-magdb    Predict peptides for all genomes in a MAG DB\n"
                        + "\nExamples:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar init\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar run --params /path/study_params.json\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar params-validate --params /path/study_params.json\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar resources-check\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-download --help\n",
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
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar deepdetect-predict \\\n"
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
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar deepdetect-magdb \\\n"
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
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar resources-check --file ~/tools/metapilot/config/meta_resources.params\n",
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
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar resources-set --template --print\n"
                        + "\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar resources-set \\\n"
                        + "    --set SYSTEM=linux --set RESOURCE_AUTO_DETECTED=true\n"
                        + "\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar resources-set \\\n"
                        + "    --set DIANN=/opt/diann/diann-linux --set FRAGPIPE=/opt/fragpipe/bin/fragpipe --print\n"
                        + "\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar resources-set --unset BASE_PATH --dry-run --print\n",
                true
        );
    }

    private static void printMagOpenSearchHelp(Options options) {
        printMagOpenSearchHelp(options, false);
    }

    private static void printMagOpenSearchHelp(Options options, boolean entrapmentMode) {
        if (entrapmentMode) {
            HelpFormatter hf = new HelpFormatter();
            hf.setWidth(110);
            hf.printHelp("java -jar metapilot-cli-<version>-shaded.jar mag-opensearch-entrapment [options]",
                    "\nEvaluation: normal marker search and genome selection, then 1:1 paired shuffled MGYG8 proteins "
                            + "in each second-search database, with rev_ decoys of real and shuffled proteins. "
                            + "Runs joint identification and all comparison arms; stops before quantification and genome inference.",
                    options,
                    "\nOutputs: <output_directory>/opensearch_entrapment_MAG_search/combined/. "
                            + "Includes identification_arms.tsv.gz, validated_psms.tsv and descriptive entrapment_summary.tsv. "
                            + "Seed and database hashes are recorded in each marker/*.fasta.entrapment.tsv manifest. "
                            + "The MAG open-search GUI has the same evaluation switch (off by default). "
                            + "A parameter file with additional_parameters.entrapment_enabled=true also selects this path from mag-opensearch.\n", true);
            return;
        }
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar mag-opensearch [options]",
                "\nRun the DDA MAG open-search workflow in headless mode (same pipeline as the Windows GUI): "
                        + "closed or open marker search -> genome selection -> sample-specific closed+open cascade "
                        + "(standalone MSFragger from the FragPipe install) -> FlashLFQ -> post-analysis. "
                        + "The GUI menu may be hidden in a release build; this command is for HPC testing.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-opensearch \\\n"
                        + "    --params /path/last_used_mag_opensearch_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n"
                        + "\n"
                        + "Required resources: FRAGPIPE (bin/fragpipe launcher; MSFragger / Crystal-C / MSBooster "
                        + "are found under tools/), FLASHLFQ. Optional: FUNCDB, SAGE, MSCONVERT.\n"
                        + "Round-1 mode, closed-stage cascade and tool toggles are read from the parameter file.\n",
                true
        );
    }

    private static void printMagOpenSearchReportHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar mag-opensearch-report [options]",
                "\nRebuild open-search final tables from an already-produced validated_psms.tsv and "
                        + "QuantifiedPeptides.tsv. No MSFragger, Crystal-C, classification or FlashLFQ. "
                        + "Applies MagResultFinalisationService (genome confidence gate, freeze selected "
                        + "genomes, excluded_evidence) and the current export columns.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-opensearch-report \\\n"
                        + "    --params /path/mag_opensearch_analysis_parameters.json \\\n"
                        + "    --out /path/Qiao_12mix_v2final \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n",
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
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-fragpipe \\\n"
                        + "    --params /path/last_used_mag_fragpipe_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n"
                        + "\n"
                        + "Search mode (Closed/Open, custom workflow overrides, etc.) is read entirely from\n"
                        + "the --params JSON file (search_mode / open_search_marker_stage /\n"
                        + "custom_basic_workflow_path / custom_combined_workflow_path) - there is no\n"
                        + "separate CLI flag for it.\n"
                        + "\n"
                        + "Identification only (stop after ID, skip combined quantification search and all post-analysis export):\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-fragpipe \\\n"
                        + "    --params /path/last_used_mag_fragpipe_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --identification-only\n"
                        + "\n"
                        + "Also emit MetaAtlas long-format taxon_abundance.parquet / function_abundance.parquet:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-fragpipe \\\n"
                        + "    --params /path/last_used_mag_fragpipe_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --export-metaatlas\n",
                true
        );
    }

    private static void printPfindMagHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar pfind-mag [options]",
                "\nRun the whole DDA MAG pFind workflow without the GUI: spectra conversion, marker search, genome "
                        + "selection, sample-specific search, FlashLFQ (and isobaric) quantification, post-analysis, "
                        + "taxonomy, functional annotation and the HTML report. Windows only. Each step is skipped "
                        + "when its output already exists, so a stopped run - or one started in the GUI - resumes "
                        + "where it stopped.",
                options,
                "\nExample (Command Prompt, one line):\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar pfind-mag "
                        + "--params C:\\path\\last_used_mag_pfind_parameters.json "
                        + "--resources C:\\MetaPilot\\config\\meta_resources.params "
                        + "--magdb C:\\MetaPilot\\config\\MAGDB_1_0.json\n"
                        + "\n"
                        + "Same, writing to another folder and searching two raw files at a time:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar pfind-mag "
                        + "--params C:\\path\\last_used_mag_pfind_parameters.json "
                        + "--resources C:\\MetaPilot\\config\\meta_resources.params "
                        + "--magdb C:\\MetaPilot\\config\\MAGDB_1_0.json "
                        + "--out D:\\results\\study1 --parallel-tasks 2\n"
                        + "\n"
                        + "Raw files must be Thermo .raw files in one folder. Exit code 0 = finished, "
                        + "1 = a step failed, 2 = invalid arguments or configuration.\n",
                true
        );
    }

    private static void printPfindMagReportHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar pfind-mag-report [options]",
                "\nResume a DDA MAG pFind workflow headlessly from an already-produced peptide quantification "
                        + "(QuantifiedPeptides.tsv for label-free, IsobaricQuantifiedPeptides.tsv for isobaric) + "
                        + "psms.tsv, and run post-analysis, taxonomy, functional annotation, and the final HTML "
                        + "report. pFind identification and quantification come from the Windows GUI or from "
                        + "pfind-mag (Windows only); only their output is required here, no other intermediate "
                        + "result files are checked.",
                options,
                "\nExample (files at the same location the Windows GUI wrote them to):\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar pfind-mag-report \\\n"
                        + "    --params /path/last_used_mag_pfind_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n"
                        + "\n"
                        + "Explicit override (files copied elsewhere on the HPC machine):\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar pfind-mag-report \\\n"
                        + "    --params /path/last_used_mag_pfind_parameters.json \\\n"
                        + "    --quantified-peptides /data/transfer/combined/QuantifiedPeptides.tsv \\\n"
                        + "    --psms /data/transfer/combined/psms.tsv\n"
                        + "\n"
                        + "Also emit MetaAtlas long-format taxon_abundance.parquet / function_abundance.parquet:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar pfind-mag-report \\\n"
                        + "    --params /path/last_used_mag_pfind_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --export-metaatlas\n",
                true
        );
    }

    private static void printDiaMagReportHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar diann-mag-report [options]",
                "\nResume a DIA MAG (DIA Marker) workflow headlessly from an already-produced DIA-NN search result "
                        + "(combined.parquet + combined.pr_matrix.tsv), and run only the memory-heavy protein/genome "
                        + "inference, taxonomy, functional annotation, and the final HTML report. The DIA-NN search "
                        + "itself still requires the Windows GUI or dia-marker on a machine with DIA-NN installed; "
                        + "no DIA-NN executable and no raw DIA spectra files are needed here. Supports Core/Focus, "
                        + "Comprehensive, Deep Together, and Deep Per File marker modes. For the two Deep modes, "
                        + "combined.parquet + combined.pr_matrix.tsv are NOT enough on their own for correct "
                        + "protein/genome inference - you also need --peptide-protein-map (recommended: a small "
                        + "file produced locally by diann-mag-peptide-map) or --round-root (the full 'marker' "
                        + "folder); this command fails with an explanatory error for Deep modes if neither is given, "
                        + "rather than silently producing wrong genome/protein assignments.",
                options,
                "\nCore/Focus/Comprehensive - HPC-local site config + explicit --out, so nothing depends on any "
                        + "path that only exists on your local computer:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar diann-mag-report \\\n"
                        + "    --params /data/transfer/last_used_dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --out /data/metapilot/studies/study1 \\\n"
                        + "    --combined-parquet /data/transfer/combined.parquet \\\n"
                        + "    --pr-matrix /data/transfer/combined.pr_matrix.tsv\n"
                        + "\n"
                        + "Deep Together / Deep Per File - also pass the compact map from diann-mag-peptide-map:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar diann-mag-report \\\n"
                        + "    --params /data/transfer/last_used_dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --out /data/metapilot/studies/study1 \\\n"
                        + "    --combined-parquet /data/transfer/combined.parquet \\\n"
                        + "    --pr-matrix /data/transfer/combined.pr_matrix.tsv \\\n"
                        + "    --peptide-protein-map /data/transfer/peptide_protein_map.tsv\n"
                        + "\n"
                        + "(--resources/--magdb must already resolve on this machine - e.g. --magdb from "
                        + "'mag-register --repository <hpc-path-to-mag-db>' - not the Windows GUI's copies. Without "
                        + "--out, the output directory falls back to a path derived from --params' 'result' field, "
                        + "which is only safe to omit when --params was exported/edited to already point at an "
                        + "HPC path.)\n"
                        + "\n"
                        + "Also emit MetaAtlas long-format taxon_abundance.parquet / function_abundance.parquet:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar diann-mag-report \\\n"
                        + "    --params /data/transfer/last_used_dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --out /data/metapilot/studies/study1 \\\n"
                        + "    --combined-parquet /data/transfer/combined.parquet \\\n"
                        + "    --pr-matrix /data/transfer/combined.pr_matrix.tsv \\\n"
                        + "    --export-metaatlas\n",
                true
        );
    }

    private static void printDiaMagPeptideMapHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar diann-mag-peptide-map [options]",
                "\nRun this LOCALLY, on the machine that still has the round-parquet folder from a completed Deep "
                        + "Together or Deep Per File DIA Marker run (Core/Focus/Comprehensive don't need this - "
                        + "their combined.parquet is already complete). Reads round1.parquet, round2.parquet, ... "
                        + "and writes one compact peptide-to-proteins TSV. Transfer only that one small file to "
                        + "HPC and pass it to 'diann-mag-report --peptide-protein-map' - no need to copy the whole "
                        + "'marker' folder.",
                options,
                "\nExample (same params file the local Deep run used):\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar diann-mag-peptide-map \\\n"
                        + "    --params /path/last_used_dia_marker_parameters.json \\\n"
                        + "    --out /path/peptide_protein_map.tsv\n",
                true
        );
    }

    private static void printDiaMarkerHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar dia-marker [options]",
                "\nRun DIA Marker in headless mode from a GUI-exported parameter file. Supported marker modes: Core/Focus, Comprehensive, Deep Together, and Deep Per File. For DIA Seed, use the dia-seed command.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar dia-marker \\\n"
                        + "    --params /path/last_used_dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n"
                        + "\n"
                        + "Also emit MetaAtlas long-format taxon_abundance.parquet / function_abundance.parquet:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar dia-marker \\\n"
                        + "    --params /path/last_used_dia_marker_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json \\\n"
                        + "    --export-metaatlas\n",
                true
        );
    }

    private static void printDiaSeedHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar dia-seed [options]",
                "\nRun DIA Seed in headless mode from a GUI-saved parameter file, with the same steps as the GUI. "
                        + "Supported seed modes (dia_settings.seed_mode): DIRECT, FOCUS (Core), COMPREHENSIVE and "
                        + "DEEP (Deep together, 5 genomes per round as in the GUI). Libraries: one seed "
                        + "library (Microbiome, created in the GUI first) and, for Core/Focus and Comprehensive only, "
                        + "optionally one host library (Standalone, .parquet).",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar dia-seed \\\n"
                        + "    --params /path/dia_seed_analysis_parameters.json \\\n"
                        + "    --resources ~/tools/metapilot/config/meta_resources.params \\\n"
                        + "    --magdb ~/tools/metapilot/config/MAGDB_1_0.json\n"
                        + "\n"
                        + "Choose the libraries on the command line:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar dia-seed \\\n"
                        + "    --params /path/dia_seed_analysis_parameters.json \\\n"
                        + "    --seed-lib MySeedLib --host-lib /path/human.parquet\n",
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
                .desc("Skip files that are already fully downloaded; continue partial downloads from where "
                        + "they stopped. Without this flag every file is downloaded again and overwritten")
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

        opts.addOption(Option.builder()
                .longOpt("protein-catalogue")
                .desc("Also download and prepare the UHGP-100 protein catalogue for the pangenome variant check "
                        + "(human-gut only; the 100% .tar.gz name is read from protein_catalogue/; "
                        + "28-32 GB download, about 45 min preparation, about 23 GB on disk)")
                .build());

        opts.addOption(Option.builder()
                .longOpt("protein-catalogue-only")
                .desc("Only download and prepare the UHGP-100 protein catalogue for an already downloaded "
                        + "<out>/<catalogue>/<version>")
                .build());

        opts.addOption(Option.builder()
                .longOpt("keep-protein-archive")
                .desc("Keep uhgp-100.tar.gz after preparation (default: delete it)")
                .build());

        return opts;
    }

    private static void printMagDownloadHelp(Options options) {
        HelpFormatter hf = new HelpFormatter();
        hf.setWidth(110);
        hf.printHelp(
                "java -jar metapilot-cli-<version>-shaded.jar mag-download [options]",
                "\nDownload one MGYG MAG catalogue+version and generate Marker FASTA files.\n"
                        + "Use --resume when re-running into a folder that already has data: complete files are "
                        + "skipped, partial files continue. Without --resume every file is downloaded again.",
                options,
                "\nExample:\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-download \\\n"
                        + "    --catalogue human-gut --version v3.0 --out /data/magdb --resume\n"
                        + "\n"
                        + "  # Default output (next to jar): <jar-folder>/magdb/<catalogue>/<version>/...\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-download \\\n"
                        + "    --catalogue human-gut --version v2.0.2 --resume\n"
                        + "\n"
                        + "  # Add the UHGP-100 protein catalogue to an existing database\n"
                        + "  java -jar metapilot-cli-1.2.0-shaded.jar mag-download \\\n"
                        + "    --catalogue human-gut --version v2.0.2 --out /data/magdb --protein-catalogue-only\n",
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
        if (error == null) {
            System.err.println("ERROR [" + context + "]: unknown");
            return;
        }
        System.err.println("ERROR [" + context + "]: " + error.getMessage());
        Throwable cause = error.getCause();
        if (cause != null && cause.getMessage() != null && !cause.getMessage().equals(error.getMessage())) {
            System.err.println("  cause: " + cause.getClass().getSimpleName() + ": " + cause.getMessage());
        }
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

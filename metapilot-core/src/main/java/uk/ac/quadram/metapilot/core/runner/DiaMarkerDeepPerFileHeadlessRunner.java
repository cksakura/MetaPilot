package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.report.MetaPilotHtmlReportService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DeepComprehensivePaths;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeCommonSupport;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeStopRule;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Headless DEEP_PER_FILE DIA Marker workflow runner (HPC version).
 *
 * <p>Mirrors the behavior of {@code DiaMarkerDeepPerFileTask} when {@code markerMode == DEEP_PER_FILE},
 * but without any JavaFX dependency. For each raw DIA file:
 * <ol>
 *   <li>run a marker no-quant search ({@code marker.parquet}) using only that raw file,</li>
 *   <li>rank/select genomes per file ({@code genome_evidence/ranked/selected.tsv}),</li>
 *   <li>run iterative rounds (one library per round = 5 genomes + previous-round peptides),</li>
 *   <li>stop on the same trend-based plateau / knee+plateau rule as the GUI workflow.</li>
 * </ol>
 *
 * <p>After all per-file workers finish, the runner builds a single global combined
 * library from the peptides identified across <i>all</i> per-file rounds (for low
 * missingness) and runs one final quantifying DIA-NN search across all DIA files.
 *
 * <p>Per-file parallelism follows the same RAM-budget rule as the GUI:
 * {@code poolSize = totalRAM_GB / 32} (sequential when total RAM &lt; 32 GB);
 * each concurrent worker is given {@code totalThreads / poolSize} DIA-NN threads.
 *
 * <p>Outputs (under {@code <resultDir>/marker_deep/}):
 * <pre>
 *   marker/&lt;rawname&gt;/marker.parquet
 *   marker/&lt;rawname&gt;/genome_evidence.tsv, genome_ranked.tsv, genome_selected.tsv
 *   marker/&lt;rawname&gt;/round1/{combined.fasta, combined.speclib, round1.parquet}
 *   marker/&lt;rawname&gt;/roundN/...
 *   combined/{combined.fasta, combined.speclib,
 *             combined.parquet, combined.pr_matrix.tsv, combined.pg_matrix.tsv}
 * </pre>
 *
 * <p>This runner uses the shared final HTML report service after post-analysis.
 */
public final class DiaMarkerDeepPerFileHeadlessRunner {

    private static final Logger LOGGER = LogManager.getLogger(DiaMarkerDeepPerFileHeadlessRunner.class);

    private static final String TASK_NAME = "DIA Deep Marker (per-file, HPC)";
    private static final int DEFAULT_GENOMES_PER_ROUND = 5;
    private static final int MAX_ITERATIVE_ROUNDS = 100;
    private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
    private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;

    private DiaMarkerDeepPerFileHeadlessRunner() {
    }

    public static final class Result {
        public final File resultDir;
        public final File markerDeepDir;
        public final File combinedDir;
        public final File combinedFasta;
        public final File combinedSpeclib;
        public final File combinedParquet;
        public final File prMatrixTsv;
        public final File pgMatrixTsv;

        Result(File resultDir, File markerDeepDir, File combinedDir, File combinedFasta, File combinedSpeclib,
               File combinedParquet, File prMatrixTsv, File pgMatrixTsv) {
            this.resultDir = resultDir;
            this.markerDeepDir = markerDeepDir;
            this.combinedDir = combinedDir;
            this.combinedFasta = combinedFasta;
            this.combinedSpeclib = combinedSpeclib;
            this.combinedParquet = combinedParquet;
            this.prMatrixTsv = prMatrixTsv;
            this.pgMatrixTsv = pgMatrixTsv;
        }
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, null, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                reportPythonHint, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint,
                             MetaPilotResourceConfig reportResourceConfig) throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true,
                generateHtmlReport, reportPythonHint, reportResourceConfig);
    }

    public static Result runWithoutPostAnalysis(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                                ProgressCallback cb, File condaPrefix, File funcDbFile) throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, false,
                false, null, null);
    }

    private static Result runInternal(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                      ProgressCallback cb, File condaPrefix, File funcDbFile,
                                      boolean exportPostAnalysis, boolean generateHtmlReport,
                                      String reportPythonHint, MetaPilotResourceConfig reportResourceConfig)
            throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (params == null) throw new IllegalArgumentException("params is null");
        if (diannExe == null || !diannExe.exists()) {
            throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));
        }

        File resultDir = new File(params.getResult());
        if (!resultDir.exists()) resultDir.mkdirs();
        if (!resultDir.isDirectory()) {
            throw new IOException("Result dir is not a directory: " + resultDir.getAbsolutePath());
        }

        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }

        MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");

        // Ensure Marker4 FASTA, genomes-taxa.tsv and Marker4.speclib exist (shared across all per-file workers).
        DiaMarkerHeadlessRunner.ensureMarker4Files(magDb, diannExe, params, cb);
        cb.onMessage(TASK_NAME + ": preparing Marker4 library");
        DiaMarkerHeadlessRunner.ensureMarker4Library(magDb, diannExe, params, cb, condaPrefix);

        // Per-file folder layout:  <resultDir>/marker_deep/marker/<rawname>/...
        File markerDeepDir = new File(resultDir, "marker_deep");
        if (!markerDeepDir.exists()) markerDeepDir.mkdirs();
        File markerRoot = new File(markerDeepDir, "marker");
        if (!markerRoot.exists()) markerRoot.mkdirs();

        LinkedHashMap<String, File> perFileDirs = new LinkedHashMap<>();
        LinkedHashMap<String, String> perFileNames = new LinkedHashMap<>();
        for (int i = 0; i < diaFiles.length; i++) {
            String raw = diaFiles[i];
            String rawName = DeepComprehensivePaths.runIdFromRawFilePath(raw);
            File perFileDir = new File(markerRoot, rawName);
            perFileDir.mkdirs();
            perFileDirs.put(raw, perFileDir);
            perFileNames.put(raw, rawName);
        }

        int totalThreads = Math.max(1, params.getThreadCount());
        DeepIterativeCommonSupport.PoolDecision decision =
                DeepIterativeCommonSupport.computeDeepPerFilePoolSizeHpc(diaFiles.length, totalThreads);
        int poolSize = decision.poolSize;
        int perTaskThreads = decision.perTaskThreads;
        String decisionMsg = String.format(
                "per-file processing budget -> poolSize=%d, perTaskThreads=%d, totalThreads=%d, files=%d, "
                        + "memorySource=%s (%s), effectiveMemoryGB=%d, ramPerTaskGB=%d, explicitOverride=%s",
                poolSize, perTaskThreads, totalThreads, diaFiles.length,
                decision.memorySource, decision.memorySourceDetail,
                decision.effectiveMemoryGB, decision.ramPerTaskGB, decision.explicitOverride);
        LOGGER.info("{}: {}", TASK_NAME, decisionMsg);
        cb.onMessage(TASK_NAME + ": " + decisionMsg);
        cb.onMessage(TASK_NAME + ": (override pool with METAPILOT_DEEP_PER_FILE_POOL=N or"
                + " ram-per-task with METAPILOT_DEEP_PER_FILE_RAM_PER_TASK_GB=N)");

        ConcurrentHashMap<String, Integer> perFileLastRound = new ConcurrentHashMap<>();

        boolean ok;
        if (poolSize <= 1) {
            ok = runPerFileSequential(diaFiles, perFileDirs, perFileNames, perFileLastRound,
                    diannExe, params, magDb, condaPrefix, perTaskThreads, cb);
        } else {
            ok = runPerFileParallel(diaFiles, perFileDirs, perFileNames, perFileLastRound,
                    diannExe, params, magDb, condaPrefix, poolSize, perTaskThreads, cb);
        }
        if (!ok) {
            throw new IOException(TASK_NAME + ": one or more per-file workers failed");
        }

        // Order perFileLastRound by original DIA file order for deterministic downstream processing.
        LinkedHashMap<String, Integer> orderedLastRound = new LinkedHashMap<>();
        for (String raw : diaFiles) {
            Integer v = perFileLastRound.get(raw);
            orderedLastRound.put(raw, v == null ? 0 : v);
        }

        // Build global combined.fasta / combined.speclib from all per-file round parquets and run final search.
        cb.onMessage(TASK_NAME + ": building global combined search from all per-file rounds");
        File combinedDir = new File(markerDeepDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();
        File combinedFasta = new File(combinedDir, "combined.fasta");
        File combinedSpeclib = new File(combinedDir, "combined.speclib");
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");

        if (!runGlobalCombinedFromPerFileRounds(diaFiles, perFileDirs, orderedLastRound,
                combinedDir, combinedFasta, combinedSpeclib, combinedParquet,
                diannExe, params, condaPrefix, cb)) {
            throw new IOException(TASK_NAME + ": global combined search failed");
        }

        if (!combinedParquet.exists()) throw new IOException("Combined quant did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Combined quant did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Combined quant did not produce: " + pgMatrix.getAbsolutePath());

        DiaMarkerHeadlessRunner.runHostSearchIfRequested(markerDeepDir, diannExe, params,
                params.getMetadata() != null ? params.getMetadata().getRawFiles() : null, cb, condaPrefix);

        // TSV/report exports go under marker_deep/, matching the Windows task layout.
        if (exportPostAnalysis) {
            try {
                File gtdbNcbi = null;
                File maybe = new File("Resources", "gtdb_ncbi");
                if (maybe.exists()) gtdbNcbi = maybe;
                DiaMagPostAnalysisTsvExporter.export(
                        markerDeepDir,
                        magDb,
                        combinedParquet,
                        prMatrix,
                        gtdbNcbi,
                        params.getMarkerMode(),
                        funcDbFile,
                        params.getMicroDb(),
                        params.getQuantitativeMode(),
                        cb,
                        DiaDeepProteinMappingService.Request.deepPerFile(markerDeepDir, diaFiles, "KR")
                );
                runHtmlReportIfRequested(markerDeepDir, params, generateHtmlReport, reportPythonHint,
                        reportResourceConfig, cb);
            } catch (Exception e) {
                cb.onMessage("WARNING: DIA marker post-analysis export failed: " + e.getMessage());
                throw e;
            }
        }

        cb.onComplete("DIA Deep Marker (per-file) complete: " + resultDir.getAbsolutePath());
        return new Result(resultDir, markerDeepDir, combinedDir, combinedFasta, combinedSpeclib,
                combinedParquet, prMatrix, pgMatrix);
    }

    private static void runHtmlReportIfRequested(File resultDir, MetaParameterDiaMarker params,
                                                 boolean generateHtmlReport, String reportPythonHint,
                                                 MetaPilotResourceConfig reportResourceConfig,
                                                 ProgressCallback cb) {
        if (!generateHtmlReport) {
            cb.onMessage("HTML report generation skipped (disabled)");
            return;
        }
        cb.onMessage(TASK_NAME + ": generating shared HTML report");
        MetaPilotHtmlReportService.Result report = MetaPilotHtmlReportService.runStandardReportWarningOnly(
                resultDir,
                "DIA_MAG",
                params != null ? params.getThreadCount() : 1,
                reportPythonHint,
                reportResourceConfig,
                cb
        );
        if (!report.warnings.isEmpty()) {
            cb.onMessage("WARNING: HTML report generation completed with warnings: "
                    + String.join("; ", report.warnings));
        }
    }

    // ----------------------- per-file orchestration -----------------------

    private static boolean runPerFileSequential(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            ConcurrentHashMap<String, Integer> perFileLastRound, File diannExe, MetaParameterDiaMarker params,
            MagDbItem magDb, File condaPrefix, int perTaskThreads, ProgressCallback cb) {
        int total = diaFiles.length;
        for (int i = 0; i < total; i++) {
            String raw = diaFiles[i];
            File perFileDir = perFileDirs.get(raw);
            String rawName = perFileNames.get(raw);
            cb.onMessage(TASK_NAME + ": per-file workflow (" + (i + 1) + "/" + total + "): " + rawName);
            int lastRound = runOnePerFileWorker(raw, perFileDir, diannExe, params, magDb, condaPrefix, perTaskThreads, cb);
            if (lastRound <= 0) {
                LOGGER.error("{}: per-file workflow failed for {}", TASK_NAME, raw);
                return false;
            }
            perFileLastRound.put(raw, lastRound);
        }
        return true;
    }

    private static boolean runPerFileParallel(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, String> perFileNames,
            ConcurrentHashMap<String, Integer> perFileLastRound, File diannExe, MetaParameterDiaMarker params,
            MagDbItem magDb, File condaPrefix, int poolSize, int perTaskThreads, ProgressCallback cb) {
        int total = diaFiles.length;
        LOGGER.info("{}: running {} per-file workflows in parallel (poolSize={}, perTaskThreads={})",
                TASK_NAME, total, poolSize, perTaskThreads);

        ExecutorService executor = Executors.newFixedThreadPool(poolSize, r -> {
            Thread t = new Thread(r, "DeepPerFile-Worker");
            t.setDaemon(true);
            return t;
        });
        AtomicBoolean anyFailed = new AtomicBoolean(false);
        List<Future<?>> futures = new ArrayList<>();

        try {
            for (String raw : diaFiles) {
                final File perFileDir = perFileDirs.get(raw);
                final String rawName = perFileNames.get(raw);
                final String rawFile = raw;
                futures.add(executor.submit(() -> {
                    try {
                        if (anyFailed.get()) return;
                        LOGGER.info("{}: [parallel] starting per-file workflow for {} (perTaskThreads={})",
                                TASK_NAME, rawName, perTaskThreads);
                        int lastRound = runOnePerFileWorker(rawFile, perFileDir, diannExe, params, magDb,
                                condaPrefix, perTaskThreads, cb);
                        if (lastRound <= 0) {
                            LOGGER.error("{}: [parallel] per-file workflow failed for {}", TASK_NAME, rawFile);
                            anyFailed.set(true);
                            return;
                        }
                        perFileLastRound.put(rawFile, lastRound);
                        LOGGER.info("{}: [parallel] completed per-file workflow for {} (lastRound={})",
                                TASK_NAME, rawName, lastRound);
                    } catch (Throwable t) {
                        LOGGER.error("{}: [parallel] exception in per-file worker for {}", TASK_NAME, rawFile, t);
                        anyFailed.set(true);
                    }
                }));
            }
            for (Future<?> f : futures) {
                try {
                    f.get();
                } catch (Exception e) {
                    LOGGER.error("{}: [parallel] worker future failed", TASK_NAME, e);
                    anyFailed.set(true);
                }
            }
        } finally {
            executor.shutdown();
            try {
                if (!executor.awaitTermination(5, TimeUnit.MINUTES)) {
                    executor.shutdownNow();
                }
            } catch (InterruptedException ignored) {
                executor.shutdownNow();
                Thread.currentThread().interrupt();
            }
        }

        return !anyFailed.get();
    }

    /**
     * Run the iterative phase (marker search + iterative rounds) for a single raw file.
     * Returns lastRound (>0) on success, <=0 on failure.
     */
    private static int runOnePerFileWorker(String raw, File perFileDir, File diannExe,
                                           MetaParameterDiaMarker params, MagDbItem magDb, File condaPrefix,
                                           int perTaskThreads, ProgressCallback cb) {
        try {
            // 1) Marker no-quant search for this raw file only -> marker.parquet
            File markerParquet = new File(perFileDir, "marker.parquet");
            File markerLib = magDb.resolveMarker4LibraryFile();
            if (markerLib == null || !markerLib.exists()) {
                LOGGER.error("{}: Marker4 library not found for per-file marker search: {}",
                        TASK_NAME, markerLib == null ? "null" : markerLib.getAbsolutePath());
                return -1;
            }
            if (!markerParquet.exists()) {
                cb.onMessage(TASK_NAME + ": marker no-quant search for " + perFileDir.getName());
                DiaNNNoQuantLibrarySearchService.run(diannExe, params, new String[] { raw }, markerLib,
                        markerParquet, perTaskThreads, true, false, cb, condaPrefix, LOGGER, TASK_NAME);
            }
            if (!markerParquet.exists()) {
                LOGGER.error("{}: marker.parquet not produced for {}", TASK_NAME, perFileDir.getName());
                return -1;
            }

            // 2) Genome ranking + selection for this raw file only
            File evidenceTsv = new File(perFileDir, "genome_evidence.tsv");
            File rankedTsv = new File(perFileDir, "genome_ranked.tsv");
            File selectedTsv = new File(perFileDir, "genome_selected.tsv");
            if (!rankedTsv.exists() || !selectedTsv.exists()) {
                cb.onMessage(TASK_NAME + ": genome selection for " + perFileDir.getName());
                DiaNNPrecursor[] precursors = DiaMarkerHeadlessRunner.readPrecursors(markerParquet, params);
                if (precursors == null || precursors.length == 0) {
                    LOGGER.error("{}: no precursors in marker.parquet for {}", TASK_NAME, perFileDir.getName());
                    return -1;
                }
                File currentFile = magDb.getCurrentFile();
                File markerFasta = new File(currentFile, "Marker4.fasta");
                File genomesTaxaTsv = new File(currentFile, "genomes-taxa.tsv");
                MagGenomeSelectionService.Request selectionRequest = new MagGenomeSelectionService.Request();
                selectionRequest.magDb = magDb;
                selectionRequest.genomeSelectorMode = MagDbItem.MARKER4;
                selectionRequest.genomeSelectionStrategy = MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
                selectionRequest.diaPrecursors = precursors;
                selectionRequest.markerFasta = markerFasta;
                selectionRequest.genomesTaxaTsv = genomesTaxaTsv;
                selectionRequest.evidenceFile = evidenceTsv;
                selectionRequest.rankedFile = rankedTsv;
                selectionRequest.selectedFile = selectedTsv;
                selectionRequest.maximumSelected = 300;
                MagGenomeSelectionService.selectGenomes(selectionRequest, cb);
            }
            if (!rankedTsv.exists()) {
                LOGGER.error("{}: genome_ranked.tsv missing for {}", TASK_NAME, perFileDir.getName());
                return -1;
            }

            // 3) Iterative rounds
            int genomesPerRound = params.getGenomesPerRound() > 0 ? params.getGenomesPerRound() : DEFAULT_GENOMES_PER_ROUND;
            String cut = "KR";
            String[] singleRawArr = new String[] { raw };
            DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
                    rankedTsv, genomesPerRound, MAX_ITERATIVE_ROUNDS, "Round", cut,
                    cb::onMessage,
                    (round, genomes) -> runRoundAppendPreviousPeptides(round, perFileDir, genomes, singleRawArr,
                            magDb, diannExe, params, condaPrefix, perTaskThreads, cb),
                    DiaMarkerDeepPerFileHeadlessRunner::computeStopRule,
                    LOGGER, TASK_NAME);

            if (loop.lastRound <= 0) {
                LOGGER.error("{}: no round completed for {}", TASK_NAME, perFileDir.getName());
                return -1;
            }
            DeepIterativeCommonSupport.writeIterativeRoundsInfo(perFileDir, loop.roundsRecord, loop.peptidesInRoundRecord,
                    loop.incrementsPerRound, loop.cumulativeRecord, loop.stopReason, loop.lastStopRuleResult,
                    CONSECUTIVE_ROUNDS_FOR_PLATEAU, LOGGER, TASK_NAME);
            return loop.lastRound;
        } catch (Throwable t) {
            LOGGER.error("{}: runOnePerFileWorker exception for {}", TASK_NAME, raw, t);
            return -1;
        }
    }

    // ----------------------- DIA-NN per-round and per-file searches -----------------------

    /**
     * One iterative round for a single raw file:
     * <pre>
     *   roundN/combined.fasta = [proteins from N genomes] [+ identified peptides from previous round]
     *   roundN/combined.speclib = DIA-NN library predicted from the FASTA
     *   roundN/roundN.parquet  = DIA-NN no-quant identification against the new library
     * </pre>
     */
    private static int runRoundAppendPreviousPeptides(int round, File perFileDir, List<String> genomes,
                                                      String[] searchFiles, MagDbItem magDb, File diannExe,
                                                      MetaParameterDiaMarker params, File condaPrefix, int threads,
                                                      ProgressCallback cb) {
        return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, perFileDir, genomes, searchFiles, true,
                "append strategy",
                (gs, fastaOut) -> writeFastaFromGenomes(magDb, gs, fastaOut),
                (gs, prevParquet, fastaOut) -> writeFastaFromGenomesAndPreviousPeptides(magDb, gs, prevParquet, fastaOut),
                (combinedFasta, combinedSpeclib) -> {
                    try {
                        DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
                        diannPar.setThreads(threads);
                        diannPar.setExportDecoy(true);
                        String cmd = diannPar.generateFasta(combinedFasta.getAbsolutePath(),
                                combinedSpeclib.getAbsolutePath(), false);
                        int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedFasta.getParentFile(), cb, condaPrefix);
                        if (exit != 0) LOGGER.warn("{}: DIA-NN round speclib build exit code {}", TASK_NAME, exit);
                        File alt = new File(combinedFasta.getParentFile(), "combined.predicted.speclib");
                        if (!combinedSpeclib.exists() && alt.exists()) {
                            try { Files.move(alt.toPath(), combinedSpeclib.toPath()); } catch (Exception ignored) {}
                        }
                    } catch (Exception e) {
                        LOGGER.error("{}: round {} build speclib failed", TASK_NAME, round, e);
                    }
                },
                (roundParquet, combinedSpeclibPath, files, reanalyse) -> {
                    try {
                        DiaNNNoQuantLibrarySearchService.run(diannExe, params, files, new File(combinedSpeclibPath),
                                roundParquet, threads, reanalyse, false, cb, condaPrefix, LOGGER, TASK_NAME);
                    } catch (Exception e) {
                        LOGGER.error("{}: round {} no-quant search failed", TASK_NAME, round, e);
                    }
                },
                (parquetFile) -> DeepIterativeCommonSupport.countIdentifiedPeptides(parquetFile, "KR", LOGGER, TASK_NAME),
                LOGGER, TASK_NAME);
    }

    // ----------------------- global combined search -----------------------

    private static boolean runGlobalCombinedFromPerFileRounds(String[] diaFiles,
            LinkedHashMap<String, File> perFileDirs, LinkedHashMap<String, Integer> perFileLastRound,
            File combinedDir, File combinedFasta, File combinedSpeclib, File combinedParquet,
            File diannExe, MetaParameterDiaMarker params, File condaPrefix, ProgressCallback cb) {

        // 1) Aggregate peptide -> proteins across ALL per-file rounds.
        HashMap<String, HashSet<String>> pepProMap = new HashMap<>();
        for (java.util.Map.Entry<String, File> entry : perFileDirs.entrySet()) {
            String raw = entry.getKey();
            File perFileDir = entry.getValue();
            Integer lastRound = perFileLastRound.get(raw);
            if (perFileDir == null || lastRound == null || lastRound <= 0) continue;
            HashMap<String, HashSet<String>> one = buildPepToProteinsFromRounds(perFileDir, lastRound);
            for (java.util.Map.Entry<String, HashSet<String>> e : one.entrySet()) {
                if (e.getKey() == null || e.getKey().isEmpty() || e.getValue() == null || e.getValue().isEmpty()) continue;
                HashSet<String> set = pepProMap.computeIfAbsent(e.getKey(), k -> new HashSet<>());
                set.addAll(e.getValue());
            }
        }

        if (pepProMap.isEmpty()) {
            LOGGER.error("{}: no peptides collected from per-file round parquets for global combined FASTA", TASK_NAME);
            return false;
        }
        LOGGER.info("{}: aggregated {} unique peptides across all per-file rounds", TASK_NAME, pepProMap.size());

        // 2) Write combined.fasta
        if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, LOGGER, TASK_NAME)) {
            return false;
        }

        // 3) Build combined.speclib
        if (!DeepIterativeCommonSupport.ensureCombinedSpeclib(combinedDir, combinedSpeclib, () -> {
            try {
                DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
                diannPar.setThreads(params.getThreadCount());
                diannPar.setExportDecoy(true);
                String cmd = diannPar.generateFasta(combinedFasta.getAbsolutePath(), combinedSpeclib.getAbsolutePath(), false);
                int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedDir, cb, condaPrefix);
                if (exit != 0) LOGGER.warn("{}: DIA-NN combined library build exit code {}", TASK_NAME, exit);
                File alt = new File(combinedDir, "combined.predicted.speclib");
                if (!combinedSpeclib.exists() && alt.exists()) {
                    try { Files.move(alt.toPath(), combinedSpeclib.toPath()); } catch (Exception ignored) {}
                }
            } catch (Exception e) {
                LOGGER.error("{}: combined speclib build failed", TASK_NAME, e);
            }
        }, LOGGER, TASK_NAME)) {
            return false;
        }

        // 4) Final quantifying search across ALL DIA files
        boolean ok = DeepIterativeCommonSupport.ensureCombinedParquet(combinedParquet, () -> {
            try {
                DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
                diannPar.setThreads(params.getThreadCount());
                diannPar.setMatrices(true);
                diannPar.setRelaxed_prot_inf(true);
                diannPar.setReanalyse(true);
                diannPar.setExportDecoy(true);

                DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
                libSearchPar.setMbr(true);
                libSearchPar.setDiaFiles(diaFiles);
                libSearchPar.setThreadCount(params.getThreadCount());
                libSearchPar.setqValue(0.01);

                String cmd = diannPar.generateLibSearch(libSearchPar, combinedParquet.getAbsolutePath(),
                        combinedSpeclib.getAbsolutePath());
                cmd = cmd + DiaMarkerHeadlessRunner.buildFArgs(diaFiles);
                int exit = DiaNNHeadlessRunner.run(diannExe, cmd, combinedDir, cb, condaPrefix);
                if (exit != 0) LOGGER.warn("{}: DIA-NN combined quant exit code {}", TASK_NAME, exit);
            } catch (Exception e) {
                LOGGER.error("{}: combined quant failed", TASK_NAME, e);
            }
        });
        return ok;
    }

    /** Build peptide -> proteins map by scanning roundN.parquet (q<=0.01, non-decoy) for one per-file folder. */
    private static HashMap<String, HashSet<String>> buildPepToProteinsFromRounds(File perFileDir, int lastRound) {
        if (perFileDir == null || lastRound <= 0) return new HashMap<>();
        return DiaDeepProteinMappingService.buildPeptideProteinMap(
                DiaDeepProteinMappingService.Request.deepTogether(perFileDir, lastRound, null, "KR"),
                LOGGER, TASK_NAME);
    }

    // ----------------------- per-round FASTA helpers -----------------------

    private static boolean writeFastaFromGenomes(MagDbItem magDb, List<String> genomes, File outputFasta) {
        File dbFolder = magDb.getOriginalDbFolder();
        if (dbFolder == null || !dbFolder.isDirectory()) {
            LOGGER.error("{}: MAG original DB folder not available", TASK_NAME);
            return false;
        }
        try (PrintWriter w = new PrintWriter(outputFasta)) {
            for (String genome : genomes) {
                File genomeFaa = resolveGenomeFile(dbFolder, genome);
                if (genomeFaa == null) {
                    LOGGER.warn("{}: genome file not found: {}", TASK_NAME, genome);
                    continue;
                }
                try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
                    String line;
                    while ((line = reader.readLine()) != null) w.println(line);
                }
            }
        } catch (IOException e) {
            LOGGER.error("{}: error writing FASTA from genomes", TASK_NAME, e);
            return false;
        }
        return outputFasta.exists() && outputFasta.length() > 0;
    }

    /**
     * FASTA = [genome proteins] + one entry per unique identified peptide from the previous round
     * (q <= 0.01, non-decoy). Each peptide entry uses a non-decoy protein id so DIA-NN can link
     * peptide -> protein.
     */
    private static boolean writeFastaFromGenomesAndPreviousPeptides(MagDbItem magDb, List<String> genomes,
            File previousRoundParquet, File outputFasta) {
        File dbFolder = magDb.getOriginalDbFolder();
        if (dbFolder == null || !dbFolder.isDirectory()) {
            LOGGER.error("{}: MAG original DB folder not available", TASK_NAME);
            return false;
        }
        DiaNNPrecursor[] precursors = null;
        if (previousRoundParquet.exists()) {
            try {
                DiaNNParquetReader reader = new DiaNNParquetReader(previousRoundParquet.getAbsolutePath(), "KR");
                precursors = reader.getDiaNNPrecursors(0.01);
            } catch (Exception e) {
                LOGGER.warn("{}: failed to read previous round parquet {}",
                        TASK_NAME, previousRoundParquet.getAbsolutePath(), e);
            }
        }
        Set<String> addedPeptideSequences = new HashSet<>();
        try (PrintWriter w = new PrintWriter(outputFasta)) {
            for (String genome : genomes) {
                File genomeFaa = resolveGenomeFile(dbFolder, genome);
                if (genomeFaa == null) {
                    LOGGER.warn("{}: genome file not found: {}", TASK_NAME, genome);
                    continue;
                }
                try (BufferedReader reader = new BufferedReader(new FileReader(genomeFaa))) {
                    String line;
                    while ((line = reader.readLine()) != null) w.println(line);
                }
            }
            if (precursors != null) {
                for (DiaNNPrecursor p : precursors) {
                    if (p == null || p.getIsDecoy()) continue;
                    String seq = p.getSeqString();
                    if (seq == null || seq.isEmpty()) continue;
                    if (addedPeptideSequences.contains(seq)) continue;
                    String[] proGroups = p.getProGroups();
                    String proteinId = null;
                    if (proGroups != null) {
                        for (String pro : proGroups) {
                            if (pro == null || pro.isEmpty()) continue;
                            // Skip common decoy prefixes
                            if (pro.startsWith("REV_") || pro.startsWith("rev_") || pro.startsWith("DECOY_")
                                    || pro.startsWith("decoy_")) continue;
                            proteinId = pro;
                            break;
                        }
                    }
                    if (proteinId == null) proteinId = "pep_" + seq.length();
                    w.println(">" + proteinId);
                    w.println(seq);
                    addedPeptideSequences.add(seq);
                }
            }
        } catch (IOException e) {
            LOGGER.error("{}: error writing FASTA (genomes + previous peptides)", TASK_NAME, e);
            return false;
        }
        return outputFasta.exists() && outputFasta.length() > 0;
    }

    private static File resolveGenomeFile(File dbFolder, String genome) {
        File f = new File(dbFolder, genome + ".faa");
        if (f.exists()) return f;
        f = new File(dbFolder, genome + ".fasta");
        if (f.exists()) return f;
        return null;
    }

    // ----------------------- stop rule (matches DEEP_MARKER profile) -----------------------

    private static DeepIterativeStopRule.StopRuleResult computeStopRule(List<Integer> increments) {
        return DeepIterativeStopRule.computeStopRoundFromIncrements(increments,
                RELATIVE_GAIN_THRESHOLD, CONSECUTIVE_ROUNDS_FOR_PLATEAU,
                DeepIterativeStopRule.AlgorithmProfile.DEEP_MARKER);
    }

    // ----------------------- misc -----------------------

}

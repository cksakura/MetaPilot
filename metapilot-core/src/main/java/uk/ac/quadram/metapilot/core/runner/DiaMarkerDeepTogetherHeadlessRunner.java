package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
import uk.ac.quadram.metapilot.task.dia.DeepIterativeCommonSupport;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeStopRule;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Headless Deep Together DIA Marker workflow runner.
 *
 * <p>The workflow performs one Marker4 search across all DIA files, ranks genomes with
 * comprehensive selection, runs iterative search rounds across all DIA files together,
 * then builds one final combined library for quantification and post-analysis.
 */
public final class DiaMarkerDeepTogetherHeadlessRunner {

    private static final Logger LOGGER = LogManager.getLogger(DiaMarkerDeepTogetherHeadlessRunner.class);

    private static final String TASK_NAME = "DIA Deep Marker (together, HPC)";
    private static final int DEFAULT_GENOMES_PER_ROUND = 5;
    private static final int MAX_ITERATIVE_ROUNDS = 100;
    private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
    private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;

    private DiaMarkerDeepTogetherHeadlessRunner() {
    }

    public static final class Result {
        public final File resultDir;
        public final File markerDeepDir;
        public final File markerDir;
        public final File combinedDir;
        public final File markerParquet;
        public final File rankedTsv;
        public final File selectedTsv;
        public final File combinedFasta;
        public final File combinedSpeclib;
        public final File combinedParquet;
        public final File prMatrixTsv;
        public final File pgMatrixTsv;

        Result(File resultDir, File markerDeepDir, File markerDir, File combinedDir, File markerParquet,
               File rankedTsv, File selectedTsv, File combinedFasta, File combinedSpeclib,
               File combinedParquet, File prMatrixTsv, File pgMatrixTsv) {
            this.resultDir = resultDir;
            this.markerDeepDir = markerDeepDir;
            this.markerDir = markerDir;
            this.combinedDir = combinedDir;
            this.markerParquet = markerParquet;
            this.rankedTsv = rankedTsv;
            this.selectedTsv = selectedTsv;
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
                             ProgressCallback cb, File condaPrefix, File funcDbFile) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true, null, null);
    }

    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint,
                             MetaPilotResourceConfig reportResourceConfig) throws Exception {
        return run(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, generateHtmlReport,
                reportPythonHint, reportResourceConfig, false);
    }

    /**
     * Same as the 9-argument {@link #run}, with an additional opt-in
     * {@code exportMetaAtlasParquet} flag (see {@code --export-metaatlas} in the CLI). This
     * is a brand new overload, not a change to any existing signature, so it cannot affect
     * any existing caller.
     */
    public static Result run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                             ProgressCallback cb, File condaPrefix, File funcDbFile,
                             boolean generateHtmlReport, String reportPythonHint,
                             MetaPilotResourceConfig reportResourceConfig, boolean exportMetaAtlasParquet)
            throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, true,
                generateHtmlReport, reportPythonHint, reportResourceConfig, exportMetaAtlasParquet);
    }

    public static Result runWithoutPostAnalysis(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                                ProgressCallback cb, File condaPrefix, File funcDbFile) throws Exception {
        return runInternal(params, diannExe, magDbConfigJson, cb, condaPrefix, funcDbFile, false,
                false, null, null, false);
    }

    private static Result runInternal(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson,
                                      ProgressCallback cb, File condaPrefix, File funcDbFile,
                                      boolean exportPostAnalysis, boolean generateHtmlReport,
                                      String reportPythonHint, MetaPilotResourceConfig reportResourceConfig,
                                      boolean exportMetaAtlasParquet)
            throws Exception {
        final ProgressCallback callback = cb == null ? new ProgressCallback.NoOp() : cb;
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

        DiaMarkerHeadlessRunner.ensureMarker4Files(magDb, diannExe, params, callback);
        callback.onMessage(TASK_NAME + ": preparing Marker4 library");
        DiaMarkerHeadlessRunner.ensureMarker4Library(magDb, diannExe, params, callback, condaPrefix);

        File markerDeepDir = new File(resultDir, "marker_deep");
        if (!markerDeepDir.exists()) markerDeepDir.mkdirs();
        File markerDir = new File(markerDeepDir, "marker");
        if (!markerDir.exists()) markerDir.mkdirs();

        File markerParquet = new File(markerDir, "marker.parquet");
        if (!markerParquet.exists()) {
            callback.onMessage(TASK_NAME + ": searching Marker4 library across all DIA files");
            DiaNNNoQuantLibrarySearchService.run(diannExe, params, diaFiles, magDb.resolveMarker4LibraryFile(),
                    markerParquet, params.getThreadCount(), true, false, callback, condaPrefix, LOGGER, TASK_NAME);
        }
        if (!markerParquet.exists()) {
            throw new IOException("Marker search did not produce: " + markerParquet.getAbsolutePath());
        }

        File evidenceTsv = new File(markerDir, "genome_evidence.tsv");
        File rankedTsv = new File(markerDir, "genome_ranked.tsv");
        File selectedTsv = new File(markerDir, "genome_selected.tsv");
        if (!rankedTsv.exists() || !selectedTsv.exists()) {
            callback.onMessage(TASK_NAME + ": selecting genomes with comprehensive marker ranking");
            DiaNNPrecursor[] precursors = DiaMarkerHeadlessRunner.readPrecursors(markerParquet, params);
            if (precursors == null || precursors.length == 0) {
                throw new IOException("No precursors found in marker parquet: " + markerParquet.getAbsolutePath());
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
            MagGenomeSelectionService.selectGenomes(selectionRequest, callback);
        }
        if (!rankedTsv.exists()) {
            throw new IOException("Genome ranked TSV missing: " + rankedTsv.getAbsolutePath());
        }

        int genomesPerRound = params.getGenomesPerRound() > 0 ? params.getGenomesPerRound() : DEFAULT_GENOMES_PER_ROUND;
        DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
                rankedTsv, genomesPerRound, MAX_ITERATIVE_ROUNDS, "Round", "KR",
                callback::onMessage,
                (round, genomes) -> runRoundAppendPreviousPeptides(round, markerDir, genomes, diaFiles,
                        magDb, diannExe, params, condaPrefix, params.getThreadCount(), callback),
                DiaMarkerDeepTogetherHeadlessRunner::computeStopRule,
                LOGGER, TASK_NAME);

        if (loop.lastRound <= 0) {
            throw new IOException(TASK_NAME + ": no iterative round completed");
        }
        DeepIterativeCommonSupport.writeIterativeRoundsInfo(markerDir, loop.roundsRecord, loop.peptidesInRoundRecord,
                loop.incrementsPerRound, loop.cumulativeRecord, loop.stopReason, loop.lastStopRuleResult,
                CONSECUTIVE_ROUNDS_FOR_PLATEAU, LOGGER, TASK_NAME);

        File combinedDir = new File(markerDeepDir, "combined");
        if (!combinedDir.exists()) combinedDir.mkdirs();
        File combinedFasta = new File(combinedDir, "combined.fasta");
        File combinedSpeclib = new File(combinedDir, "combined.speclib");
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");

        if (!runGlobalCombinedFromRounds(diaFiles, markerDir, loop.lastRound, combinedDir, combinedFasta,
                combinedSpeclib, combinedParquet, diannExe, params, condaPrefix, callback)) {
            throw new IOException(TASK_NAME + ": global combined search failed");
        }
        if (!combinedParquet.exists()) throw new IOException("Combined quant did not produce: " + combinedParquet.getAbsolutePath());
        if (!prMatrix.exists()) throw new IOException("Combined quant did not produce: " + prMatrix.getAbsolutePath());
        if (!pgMatrix.exists()) throw new IOException("Combined quant did not produce: " + pgMatrix.getAbsolutePath());

        DiaMarkerHeadlessRunner.runHostSearchIfRequested(markerDeepDir, diannExe, params, diaFiles, callback, condaPrefix);

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
                        callback,
                        DiaDeepProteinMappingService.Request.deepTogether(markerDir, loop.lastRound, diaFiles, "KR"),
                        exportMetaAtlasParquet
                );
                runHtmlReportIfRequested(markerDeepDir, params, generateHtmlReport, reportPythonHint,
                        reportResourceConfig, callback);
            } catch (Exception e) {
                callback.onMessage("WARNING: DIA marker post-analysis export failed: " + e.getMessage());
                throw e;
            }
        }

        callback.onComplete("DIA Deep Marker (together) complete: " + markerDeepDir.getAbsolutePath());
        return new Result(resultDir, markerDeepDir, markerDir, combinedDir, markerParquet, rankedTsv, selectedTsv,
                combinedFasta, combinedSpeclib, combinedParquet, prMatrix, pgMatrix);
    }

    private static int runRoundAppendPreviousPeptides(int round, File markerDir, List<String> genomes,
                                                      String[] searchFiles, MagDbItem magDb, File diannExe,
                                                      MetaParameterDiaMarker params, File condaPrefix,
                                                      int threads, ProgressCallback cb) {
        return DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, markerDir, genomes, searchFiles, true,
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

    private static boolean runGlobalCombinedFromRounds(String[] diaFiles, File markerDir, int lastRound,
                                                       File combinedDir, File combinedFasta, File combinedSpeclib,
                                                       File combinedParquet, File diannExe,
                                                       MetaParameterDiaMarker params, File condaPrefix,
                                                       ProgressCallback cb) {
        HashMap<String, HashSet<String>> pepProMap = buildPepToProteinsFromRounds(markerDir, lastRound);
        if (pepProMap.isEmpty()) {
            LOGGER.error("{}: no peptides collected from round parquets for global combined FASTA", TASK_NAME);
            return false;
        }
        LOGGER.info("{}: aggregated {} unique peptides across all rounds", TASK_NAME, pepProMap.size());

        if (!DeepIterativeCommonSupport.writeCombinedFastaFromPeptideProteinMap(combinedFasta, pepProMap, LOGGER, TASK_NAME)) {
            return false;
        }

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

        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        File pgMatrix = new File(combinedDir, "combined.pg_matrix.tsv");
        if (combinedParquet.exists() && prMatrix.exists() && pgMatrix.exists()) {
            return true;
        }
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
            return false;
        }
        return combinedParquet.exists() && prMatrix.exists() && pgMatrix.exists();
    }

    private static HashMap<String, HashSet<String>> buildPepToProteinsFromRounds(File markerDir, int lastRound) {
        if (markerDir == null || lastRound <= 0) return new HashMap<>();
        return DiaDeepProteinMappingService.buildPeptideProteinMap(
                DiaDeepProteinMappingService.Request.deepTogether(markerDir, lastRound, null, "KR"),
                LOGGER, TASK_NAME);
    }

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
                    if (seq == null || seq.isEmpty() || addedPeptideSequences.contains(seq)) continue;
                    String proteinId = firstNonDecoyProtein(p.getProGroups());
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

    private static String firstNonDecoyProtein(String[] proGroups) {
        if (proGroups == null) return null;
        for (String pro : proGroups) {
            if (pro == null || pro.isEmpty()) continue;
            if (pro.startsWith("REV_") || pro.startsWith("rev_") || pro.startsWith("DECOY_")
                    || pro.startsWith("decoy_")) continue;
            return pro;
        }
        return null;
    }

    private static File resolveGenomeFile(File dbFolder, String genome) {
        File f = new File(dbFolder, genome + ".faa");
        if (f.exists()) return f;
        f = new File(dbFolder, genome + ".fasta");
        if (f.exists()) return f;
        return null;
    }

    private static DeepIterativeStopRule.StopRuleResult computeStopRule(List<Integer> increments) {
        return DeepIterativeStopRule.computeStopRoundFromIncrements(increments,
                RELATIVE_GAIN_THRESHOLD, CONSECUTIVE_ROUNDS_FOR_PLATEAU,
                DeepIterativeStopRule.AlgorithmProfile.DEEP_MARKER);
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
}

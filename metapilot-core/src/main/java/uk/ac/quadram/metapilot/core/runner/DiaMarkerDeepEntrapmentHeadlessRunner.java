package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Predicate;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.dia.EntrapmentShuffle;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeCommonSupport;
import uk.ac.quadram.metapilot.task.dia.DeepIterativeStopRule;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Entrapment version of DIA Marker Deep together ({@code dia-marker-deep-entrapment} and the Deep together
 * GUI switch), for evaluating the error of the iterative round searches. It does not change any other workflow.
 * <ol>
 *   <li>Marker search and genome ranking: identical to Deep together.</li>
 *   <li>Each round's database = the round's target genomes (5 by default) + a shuffled copy of each of them
 *   ({@link EntrapmentShuffle}, IDs {@code MGYG8...}) + the peptides identified in the previous round, target and
 *   entrapment alike. So entrapment peptides that pass a round are carried forward exactly like target peptides.</li>
 *   <li>The entrapment never decides anything: the genomes of each round come from the ranked list, and the stop
 *   rules count target peptides only.</li>
 *   <li>The final combined search is built from all round peptides (target and entrapment), as in Deep together,
 *   followed by the standard post-analysis.</li>
 * </ol>
 * Output under {@code marker_deep_entrapment/}, as {@code marker_deep/} plus:
 * <ul>
 *   <li>{@code marker/entrapment_rounds.tsv}: per round, target and entrapment peptides (q &lt; 0.01), cumulative,
 *   and the combined-estimator FDP of the accumulated peptides;</li>
 *   <li>{@code entrapment_summary.tsv}: final search, precursors and peptides at global q thresholds, and the
 *   reported peptides (final_peptides.tsv), with FDP = N_e (1 + 1/r) / (N_t + N_e), r = 1 (one shuffled copy
 *   per target protein);</li>
 *   <li>{@code entrapment_peptides.tsv}: every entrapment peptide of the final search with its paired source
 *   peptide (same span of the source protein).</li>
 * </ul>
 * Shuffled peptides that also occur as target peptides count as target here (any target protein); the exact
 * shared-sequence classification against the full catalogue is left to the analysis scripts.
 */
public final class DiaMarkerDeepEntrapmentHeadlessRunner {

    private static final Logger LOGGER = LogManager.getLogger(DiaMarkerDeepEntrapmentHeadlessRunner.class);
    private static final String TASK_NAME = "DIA Deep Marker entrapment";
    public static final String RESULT_FOLDER = "marker_deep_entrapment";
    private static final int DEFAULT_GENOMES_PER_ROUND = 5;
    private static final int MAX_ITERATIVE_ROUNDS = 100;
    /** Same stop rule as Deep together. */
    private static final double RELATIVE_GAIN_THRESHOLD = 0.015;
    private static final int CONSECUTIVE_ROUNDS_FOR_PLATEAU = 2;
    private static final double[] Q_THRESHOLDS = { 0.001, 0.002, 0.005, 0.01 };
    /** Entrapment-to-target size ratio: one shuffled copy of every target protein. */
    private static final double R = 1.0;

    private DiaMarkerDeepEntrapmentHeadlessRunner() {
    }

    public static File run(MetaParameterDiaMarker params, File diannExe, File magDbConfigJson, ProgressCallback cb,
                           File condaPrefix, File funcDbFile, MetaPilotResourceConfig resourceConfig, long seed,
                           int idDigit) throws Exception {
        final ProgressCallback callback = cb == null ? new ProgressCallback.NoOp() : cb;
        if (params == null) throw new IllegalArgumentException("params is null");
        if (diannExe == null || !diannExe.exists()) {
            throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));
        }
        if (idDigit < 1 || idDigit > 9) throw new IllegalArgumentException("entrapment ID digit must be 1-9");
        File resultDir = new File(params.getResult());
        if (!resultDir.exists()) resultDir.mkdirs();
        String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files found in params metadata (input_output.experiments)");
        }
        MagDbItem magDb = DiaMarkerHeadlessRunner.resolveMagDb(params, magDbConfigJson);
        if (magDb == null) throw new IOException("MAG DB could not be resolved from params + magdb config");

        // 1. marker search and genome ranking (as Deep together)
        DiaMarkerHeadlessRunner.ensureMarker4Files(magDb, diannExe, params, callback, resourceConfig);
        callback.onMessage(TASK_NAME + ": preparing Marker4 library");
        DiaMarkerHeadlessRunner.ensureMarker4Library(magDb, diannExe, params, callback, condaPrefix);
        File deepDir = new File(resultDir, RESULT_FOLDER);
        File markerDir = new File(deepDir, "marker");
        markerDir.mkdirs();
        File markerParquet = new File(markerDir, "marker.parquet");
        if (!markerParquet.exists()) {
            callback.onMessage(TASK_NAME + ": searching Marker4 library across all DIA files");
            DiaNNNoQuantLibrarySearchService.run(diannExe, params, diaFiles, magDb.resolveMarker4LibraryFile(),
                    markerParquet, params.getThreadCount(), true, false, callback, condaPrefix, LOGGER, TASK_NAME);
        }
        if (!markerParquet.exists()) throw new IOException("Marker search did not produce: " + markerParquet);
        File rankedTsv = new File(markerDir, "genome_ranked.tsv");
        File selectedTsv = new File(markerDir, "genome_selected.tsv");
        if (!rankedTsv.exists() || !selectedTsv.exists()) {
            callback.onMessage(TASK_NAME + ": selecting genomes with comprehensive marker ranking");
            DiaNNPrecursor[] precursors = DiaMarkerHeadlessRunner.readPrecursors(markerParquet, params);
            if (precursors == null || precursors.length == 0) {
                throw new IOException("No precursors found in marker parquet: " + markerParquet);
            }
            MagGenomeSelectionService.Request req = new MagGenomeSelectionService.Request();
            req.magDb = magDb;
            req.genomeSelectorMode = MagDbItem.MARKER4;
            req.genomeSelectionStrategy = MagGenomeSelectionService.SELECTION_STRATEGY_COMPREHENSIVE;
            req.diaPrecursors = precursors;
            req.markerFasta = new File(magDb.getCurrentFile(), "Marker4.fasta");
            req.genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
            req.evidenceFile = new File(markerDir, "genome_evidence.tsv");
            req.rankedFile = rankedTsv;
            req.selectedFile = selectedTsv;
            req.maximumSelected = 300;
            MagGenomeSelectionService.selectGenomes(req, callback);
        }
        if (!rankedTsv.exists()) throw new IOException("Genome ranked TSV missing: " + rankedTsv);

        // 2. rounds: target genomes + their shuffled copies + previous round peptides; target-only counting
        int threads = params.getThreadCount();
        Predicate<DiaNNPrecursor> targetOnly = p -> EntrapmentShuffle.hasTargetProtein(p, idDigit);
        int genomesPerRound = params.getGenomesPerRound() > 0 ? params.getGenomesPerRound() : DEFAULT_GENOMES_PER_ROUND;
        callback.onMessage(TASK_NAME + ": rounds of " + genomesPerRound + " target genomes + " + genomesPerRound
                + " shuffled copies (IDs MGYG" + idDigit + "..., seed " + seed + ")");
        DeepIterativeCommonSupport.IterativeRoundLoopResult loop = DeepIterativeCommonSupport.runIterativeRoundLoop(
                rankedTsv, genomesPerRound, MAX_ITERATIVE_ROUNDS, "Round", "KR", false, 0,
                callback::onMessage,
                (round, genomes) -> DeepIterativeCommonSupport.runAppendPreviousPeptidesRound(round, markerDir,
                        genomes, diaFiles, true, "entrapment",
                        (gs, fastaOut) -> DiaMarkerDeepTogetherHeadlessRunner.writeFastaFromGenomes(magDb, gs, fastaOut)
                                && appendShuffledGenomes(magDb, gs, fastaOut, idDigit, seed),
                        (gs, prevParquet, fastaOut) -> DiaMarkerDeepTogetherHeadlessRunner
                                .writeFastaFromGenomesAndPreviousPeptides(magDb, gs, prevParquet, fastaOut)
                                && appendShuffledGenomes(magDb, gs, fastaOut, idDigit, seed),
                        (combinedFasta, combinedSpeclib) -> buildRoundSpeclib(params, threads, combinedFasta,
                                combinedSpeclib, diannExe, callback, condaPrefix),
                        (roundParquet, speclibPath, files, reanalyse) -> {
                            try {
                                DiaNNNoQuantLibrarySearchService.run(diannExe, params, files, new File(speclibPath),
                                        roundParquet, threads, reanalyse, false, callback, condaPrefix, LOGGER, TASK_NAME);
                            } catch (Exception e) {
                                LOGGER.error("{}: round no-quant search failed", TASK_NAME, e);
                            }
                        },
                        parquet -> countPeptides(parquet, targetOnly),
                        LOGGER, TASK_NAME),
                increments -> DeepIterativeStopRule.computeStopRoundFromIncrements(increments, RELATIVE_GAIN_THRESHOLD,
                        CONSECUTIVE_ROUNDS_FOR_PLATEAU, DeepIterativeStopRule.AlgorithmProfile.DEEP_MARKER),
                targetOnly, LOGGER, TASK_NAME);
        if (loop.lastRound <= 0) throw new IOException(TASK_NAME + ": no iterative round completed");
        DeepIterativeCommonSupport.writeIterativeRoundsInfo(markerDir, loop.roundsRecord, loop.peptidesInRoundRecord,
                loop.incrementsPerRound, loop.cumulativeRecord, loop.stopReason, loop.lastStopRuleResult,
                CONSECUTIVE_ROUNDS_FOR_PLATEAU, LOGGER, TASK_NAME);
        writeRoundsTable(markerDir, loop.lastRound, idDigit);

        // 3. final combined search from all round peptides (as Deep together, incl. the host library if set)
        File combinedDir = new File(deepDir, "combined");
        combinedDir.mkdirs();
        File combinedParquet = new File(combinedDir, "combined.parquet");
        File prMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
        List<File> hostParquets = DiaMarkerHeadlessRunner.resolveHostParquets(params, new File(deepDir, "host_library"),
                diannExe, callback, condaPrefix);
        if (!DiaMarkerDeepTogetherHeadlessRunner.runGlobalCombinedFromRounds(diaFiles, markerDir, loop.lastRound,
                combinedDir, new File(combinedDir, "combined.fasta"), new File(combinedDir, "combined.speclib"),
                combinedParquet, hostParquets, diannExe, params, condaPrefix, callback)) {
            throw new IOException(TASK_NAME + ": final combined search failed");
        }

        // 4. standard post-analysis; the entrapment report below does not depend on its later steps
        try {
            File gtdbNcbi = new File("Resources", "gtdb_ncbi");
            DiaMagPostAnalysisTsvExporter.export(deepDir, magDb, combinedParquet, prMatrix,
                    gtdbNcbi.exists() ? gtdbNcbi : null, params.getMarkerMode(), funcDbFile, params.getMicroDb(),
                    params.getQuantitativeMode(), callback,
                    DiaDeepProteinMappingService.Request.deepTogether(markerDir, loop.lastRound, diaFiles, "KR"),
                    false);
        } catch (Exception e) {
            callback.onMessage("WARNING: post-analysis export stopped (" + e.getMessage()
                    + "); the entrapment report uses combined.parquet and final_peptides.tsv if written");
            LOGGER.warn("{}: post-analysis export failed", TASK_NAME, e);
        }

        // 5. entrapment report
        writeFinalReport(deepDir, combinedParquet, magDb, idDigit, seed, callback);
        callback.onComplete(TASK_NAME + " complete: " + deepDir.getAbsolutePath());
        return deepDir;
    }

    static boolean appendShuffledGenomes(MagDbItem magDb, List<String> genomes, File fasta, int digit, long seed) {
        File dbFolder = magDb.getOriginalDbFolder();
        try (PrintWriter w = new PrintWriter(new java.io.FileWriter(fasta, StandardCharsets.UTF_8, true))) {
            int n = 0;
            for (String g : genomes) {
                File faa = genomeFile(dbFolder, g);
                if (faa == null) {
                    LOGGER.warn("{}: genome file not found for its shuffled copy: {}", TASK_NAME, g);
                    continue;
                }
                n += EntrapmentShuffle.appendShuffledGenome(faa, w, digit, seed);
            }
            LOGGER.info("{}: {} shuffled entrapment proteins appended to {}", TASK_NAME, n, fasta.getName());
            return true;
        } catch (IOException e) {
            LOGGER.error("{}: appending shuffled genomes failed", TASK_NAME, e);
            return false;
        }
    }

    private static File genomeFile(File dbFolder, String genome) {
        if (dbFolder == null) return null;
        File f = new File(dbFolder, genome + ".faa");
        if (f.exists()) return f;
        f = new File(dbFolder, genome + ".fasta");
        return f.exists() ? f : null;
    }

    private static void buildRoundSpeclib(MetaParameterDiaMarker params, int threads, File fasta, File speclib,
                                          File diannExe, ProgressCallback cb, File condaPrefix) {
        try {
            DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
            diannPar.setThreads(threads);
            diannPar.setExportDecoy(true);
            int exit = DiaNNHeadlessRunner.run(diannExe,
                    diannPar.generateFasta(fasta.getAbsolutePath(), speclib.getAbsolutePath(), false),
                    fasta.getParentFile(), cb, condaPrefix);
            if (exit != 0) LOGGER.warn("{}: DIA-NN round speclib build exit code {}", TASK_NAME, exit);
            File alt = new File(fasta.getParentFile(), "combined.predicted.speclib");
            if (!speclib.exists() && alt.exists()) Files.move(alt.toPath(), speclib.toPath());
        } catch (Exception e) {
            LOGGER.error("{}: round speclib build failed", TASK_NAME, e);
        }
    }

    private static int countPeptides(File parquet, Predicate<DiaNNPrecursor> filter) {
        Set<String> peps = new HashSet<>();
        for (DiaNNPrecursor p : readPassing(parquet, 0.01)) {
            if (filter.test(p)) peps.add(p.getSeqString());
        }
        return peps.size();
    }

    private static List<DiaNNPrecursor> readPassing(File parquet, double q) {
        List<DiaNNPrecursor> out = new ArrayList<>();
        if (parquet == null || !parquet.exists()) return out;
        try {
            DiaNNPrecursor[] ps = new DiaNNParquetReader(parquet.getAbsolutePath(), "KR").getDiaNNPrecursors(q);
            if (ps != null) {
                for (DiaNNPrecursor p : ps) {
                    if (p != null && !p.getIsDecoy() && p.getSeqString() != null && !p.getSeqString().isEmpty()) {
                        out.add(p);
                    }
                }
            }
        } catch (Exception e) {
            LOGGER.warn("{}: failed to read {}", TASK_NAME, parquet, e);
        }
        return out;
    }

    static double fdp(int target, int entrapment) {
        int n = target + entrapment;
        return n == 0 ? 0.0 : entrapment * (1.0 + 1.0 / R) / n;
    }

    private static void writeRoundsTable(File markerDir, int lastRound, int digit) throws IOException {
        Set<String> cumT = new HashSet<>();
        Set<String> cumE = new HashSet<>();
        try (PrintWriter w = new PrintWriter(new File(markerDir, "entrapment_rounds.tsv"), StandardCharsets.UTF_8)) {
            w.println("Round\tTarget peptides\tEntrapment peptides\tRound FDP (combined, r=1)\t"
                    + "Cumulative target peptides\tCumulative entrapment peptides\tCumulative FDP (combined, r=1)");
            for (int round = 1; round <= lastRound; round++) {
                File parquet = new File(new File(markerDir, "round" + round), "round" + round + ".parquet");
                Set<String> t = new HashSet<>();
                Set<String> e = new HashSet<>();
                for (DiaNNPrecursor p : readPassing(parquet, 0.01)) {
                    if (EntrapmentShuffle.hasTargetProtein(p, digit)) t.add(p.getSeqString());
                    else if (EntrapmentShuffle.onlyEntrapmentProteins(p, digit)) e.add(p.getSeqString());
                }
                cumT.addAll(t);
                cumE.addAll(e);
                w.println(round + "\t" + t.size() + "\t" + e.size() + "\t" + fmt(fdp(t.size(), e.size())) + "\t"
                        + cumT.size() + "\t" + cumE.size() + "\t" + fmt(fdp(cumT.size(), cumE.size())));
            }
        }
    }

    private static void writeFinalReport(File deepDir, File combinedParquet, MagDbItem magDb, int digit, long seed,
                                         ProgressCallback cb) throws IOException {
        // final search at global q thresholds
        List<DiaNNPrecursor> all = readPassing(combinedParquet, Q_THRESHOLDS[Q_THRESHOLDS.length - 1]);
        List<String[]> rows = new ArrayList<>();
        for (double q : Q_THRESHOLDS) {
            Set<String> tPrec = new HashSet<>(), ePrec = new HashSet<>(), tPep = new HashSet<>(), ePep = new HashSet<>();
            for (DiaNNPrecursor p : all) {
                if (!(p.getGlobalQvalue() < q)) continue;
                String prec = p.getModSeqString() + "/" + p.getCharge();
                if (EntrapmentShuffle.hasTargetProtein(p, digit)) {
                    tPrec.add(prec);
                    tPep.add(p.getSeqString());
                } else if (EntrapmentShuffle.onlyEntrapmentProteins(p, digit)) {
                    ePrec.add(prec);
                    ePep.add(p.getSeqString());
                }
            }
            rows.add(row("final search, global q < " + q, "precursors", tPrec.size(), ePrec.size()));
            rows.add(row("final search, global q < " + q, "peptides", tPep.size(), ePep.size()));
        }
        // reported peptides (after MetaPilot's genome/protein scoring)
        Map<String, String> reported = readFinalPeptides(new File(deepDir, "final_peptides.tsv"));
        int rt = 0, re = 0;
        for (Map.Entry<String, String> en : reported.entrySet()) {
            String[] pros = en.getValue().split("[;,\\s]+");
            boolean target = false, entr = false;
            for (String pro : pros) {
                if (pro.isEmpty() || EntrapmentShuffle.isDecoyId(pro)) continue;
                if (EntrapmentShuffle.isEntrapment(pro, digit)) entr = true;
                else target = true;
            }
            if (target) rt++;
            else if (entr) re++;
        }
        if (!reported.isEmpty()) rows.add(row("reported (final_peptides.tsv)", "peptides", rt, re));
        try (PrintWriter w = new PrintWriter(new File(deepDir, "entrapment_summary.tsv"), StandardCharsets.UTF_8)) {
            w.println("Set\tLevel\tTarget\tEntrapment\tFDP (combined estimator, r=1)");
            for (String[] r : rows) w.println(String.join("\t", r));
            w.println("# entrapment: shuffled copies (IDs MGYG" + digit + "..., seed " + seed + ") of every target genome "
                    + "searched in each round; target = any non-entrapment protein; FDP = N_e (1 + 1/r) / (N_t + N_e)");
        }
        // entrapment peptides of the final search with their paired source peptides
        Map<String, String[]> proteinCache = new HashMap<>();
        TreeMap<String, Object[]> entr = new TreeMap<>();
        for (DiaNNPrecursor p : all) {
            if (!EntrapmentShuffle.onlyEntrapmentProteins(p, digit)) continue;
            Object[] cur = entr.get(p.getSeqString());
            if (cur == null || p.getGlobalQvalue() < (double) cur[1]) {
                entr.put(p.getSeqString(), new Object[] { p.getProGroups(), p.getGlobalQvalue() });
            }
        }
        try (PrintWriter w = new PrintWriter(new File(deepDir, "entrapment_peptides.tsv"), StandardCharsets.UTF_8)) {
            w.println("Peptide\tEntrapment proteins\tBest global q\tReported\tSource protein\tSource peptide");
            for (Map.Entry<String, Object[]> en : entr.entrySet()) {
                String pep = en.getKey();
                String[] pros = (String[]) en.getValue()[0];
                String srcProt = "", srcPep = "";
                for (String pro : pros) {
                    if (!EntrapmentShuffle.isEntrapment(pro, digit)) continue;
                    String[] pair = sourcePair(magDb, pro, digit, seed, proteinCache);
                    if (pair == null) continue;
                    int i = pair[1].indexOf(pep);
                    if (i >= 0) {
                        srcProt = pair[0];
                        srcPep = pair[2].substring(i, i + pep.length());
                        break;
                    }
                }
                w.println(pep + "\t" + String.join(";", pros) + "\t" + en.getValue()[1] + "\t"
                        + (reported.containsKey(pep) ? "yes" : "no") + "\t" + srcProt + "\t" + srcPep);
            }
        }
        cb.onMessage(TASK_NAME + ": entrapment report written (" + entr.size() + " entrapment peptides in the final search)");
    }

    private static String[] row(String set, String level, int t, int e) {
        return new String[] { set, level, String.valueOf(t), String.valueOf(e), fmt(fdp(t, e)) };
    }

    private static String fmt(double v) {
        return String.format(Locale.ROOT, "%.5f", v);
    }

    /** {source protein ID, entrapment sequence, source sequence} for an entrapment protein ID, or null. */
    private static String[] sourcePair(MagDbItem magDb, String entrapmentProtein, int digit, long seed,
                                       Map<String, String[]> cache) {
        String[] hit = cache.get(entrapmentProtein);
        if (hit != null || cache.containsKey(entrapmentProtein)) return hit;
        String sourceProtein = "MGYG0" + entrapmentProtein.substring(5);
        String genome = sourceProtein.contains("_") ? sourceProtein.substring(0, sourceProtein.indexOf('_')) : sourceProtein;
        File faa = genomeFile(magDb.getOriginalDbFolder(), genome);
        if (faa != null && !cache.containsKey("#" + genome)) {
            cache.put("#" + genome, new String[0]);
            for (Map.Entry<String, String> pr : readFasta(faa).entrySet()) {
                String eid = EntrapmentShuffle.entrapmentId(pr.getKey(), digit);
                cache.put(eid, new String[] { pr.getKey(), EntrapmentShuffle.shuffle(pr.getValue(), eid, seed),
                        pr.getValue().toUpperCase() });
            }
        }
        if (!cache.containsKey(entrapmentProtein)) cache.put(entrapmentProtein, null);
        return cache.get(entrapmentProtein);
    }

    private static Map<String, String> readFasta(File f) {
        Map<String, String> out = new LinkedHashMap<>();
        try (BufferedReader r = new BufferedReader(new FileReader(f, StandardCharsets.UTF_8))) {
            String id = null;
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (id != null) out.put(id, sb.toString());
                    String h = line.substring(1).trim();
                    id = h.contains(" ") ? h.substring(0, h.indexOf(' ')) : h;
                    sb.setLength(0);
                } else {
                    sb.append(line.trim());
                }
            }
            if (id != null) out.put(id, sb.toString());
        } catch (IOException e) {
            LOGGER.warn("{}: failed to read {}", TASK_NAME, f, e);
        }
        return out;
    }

    /** Base sequence -> proteins of final_peptides.tsv (empty when missing). */
    private static Map<String, String> readFinalPeptides(File tsv) {
        Map<String, String> out = new HashMap<>();
        if (!tsv.isFile()) return out;
        try (BufferedReader r = new BufferedReader(new FileReader(tsv, StandardCharsets.UTF_8))) {
            String header = r.readLine();
            if (header == null) return out;
            String[] h = header.split("\t", -1);
            int seq = -1, pro = -1;
            for (int i = 0; i < h.length; i++) {
                if (h[i].equalsIgnoreCase("Base sequence")) seq = i;
                if (h[i].equalsIgnoreCase("Mapped proteins") || (pro < 0 && h[i].equalsIgnoreCase("Proteins"))) pro = i;
            }
            if (seq < 0 || pro < 0) return out;
            String line;
            while ((line = r.readLine()) != null) {
                String[] f = line.split("\t", -1);
                if (f.length > Math.max(seq, pro)) out.merge(f[seq], f[pro], (a, b) -> a + ";" + b);
            }
        } catch (IOException e) {
            LOGGER.warn("{}: failed to read {}", TASK_NAME, tsv, e);
        }
        return out;
    }
}

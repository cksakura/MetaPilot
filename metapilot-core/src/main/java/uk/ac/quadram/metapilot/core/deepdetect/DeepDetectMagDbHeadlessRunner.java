package uk.ac.quadram.metapilot.core.deepdetect;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Headless DeepDetect runner for predicting peptides for an entire MAG database.
 *
 * <p>Matches the Windows GUI behaviour:
 * - inputs: FASTA files under {@code magdb.getOriginalDbFolder()}
 * - outputs: {@code <db-parent>/predicted/<genome>.predicted.tsv}
 * - skips genomes that already have predicted TSV output
 */
public final class DeepDetectMagDbHeadlessRunner {

    private DeepDetectMagDbHeadlessRunner() {
    }

    /**
     * Default parallelism for external DeepDetect processes.
     *
     * <p>Rule: use about half the machine cores, but prefer at least 2 when possible.
     * Falls back to 1 only on truly single-core machines.
     */
    private static int defaultThreads() {
        int cpu = Runtime.getRuntime().availableProcessors();
        if (cpu <= 1) return 1;
        if (cpu == 2) return 2;
        return Math.max(2, cpu / 2);
    }

    public static final class Result {
        public final File predictedDir;
        public final int totalGenomes;
        public final int predictedGenomes;
        public final int skippedExisting;
        public final int failedGenomes;

        Result(File predictedDir, int totalGenomes, int predictedGenomes, int skippedExisting, int failedGenomes) {
            this.predictedDir = predictedDir;
            this.totalGenomes = totalGenomes;
            this.predictedGenomes = predictedGenomes;
            this.skippedExisting = skippedExisting;
            this.failedGenomes = failedGenomes;
        }
    }

    public static Result runAllGenomes(
            MagDbItem magdb,
            File pythonExe,
            File deepDetectMainPy,
            DeepDetectParameter template,
            String regular,
            int threads,
            ProgressCallback cb
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (magdb == null) throw new IllegalArgumentException("magdb is null");
        if (pythonExe == null || !pythonExe.exists()) throw new IllegalArgumentException("pythonExe not found: " + pythonExe);
        if (deepDetectMainPy == null || !deepDetectMainPy.exists()) throw new IllegalArgumentException("deepDetectMainPy not found: " + deepDetectMainPy);
        if (template == null) template = new DeepDetectParameter();
        if (threads <= 0) threads = defaultThreads();

        // Ensure current version is set (MagDbItem relies on usedVersion for original_db resolution).
        if (magdb.getUsedVersion() != null) {
            magdb.setUsedVersion(magdb.getUsedVersion());
        }

        File dbFolder = magdb.getOriginalDbFolder();
        if (dbFolder == null || !dbFolder.exists()) {
            throw new IllegalArgumentException("MAG DB original_db folder not found for item: " + magdb.getIdentifier());
        }

        File predictedDir = new File(dbFolder.getParentFile(), "predicted");
        return runAllGenomes(dbFolder, predictedDir, pythonExe, deepDetectMainPy, template, regular, threads, cb);
    }

    /**
     * Predict peptides for all genome FASTA files under {@code originalDbDir}, writing
     * {@code <predictedDir>/<genome>.predicted.tsv}.
     *
     * <p>Intended for CLI usage where MAG DB location is resolved from {@code --out/--catalogue/--version}.
     */
    public static Result runAllGenomes(
            File originalDbDir,
            File predictedDir,
            File pythonExe,
            File deepDetectMainPy,
            DeepDetectParameter template,
            String regular,
            int threads,
            ProgressCallback cb
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (originalDbDir == null || !originalDbDir.exists()) {
            throw new IllegalArgumentException("originalDbDir not found: " + originalDbDir);
        }
        if (predictedDir == null) {
            throw new IllegalArgumentException("predictedDir is null");
        }
        if (pythonExe == null || !pythonExe.exists()) throw new IllegalArgumentException("pythonExe not found: " + pythonExe);
        if (deepDetectMainPy == null || !deepDetectMainPy.exists()) throw new IllegalArgumentException("deepDetectMainPy not found: " + deepDetectMainPy);
        if (template == null) template = new DeepDetectParameter();
        if (threads <= 0) threads = defaultThreads();

        predictedDir.mkdirs();

        File[] files = originalDbDir.listFiles();
        if (files == null) files = new File[0];

        List<File> fastaFiles = new ArrayList<>();
        for (File f : files) {
            if (f == null || !f.isFile()) continue;
            String n = f.getName().toLowerCase(Locale.ROOT);
            if (n.endsWith(".faa") || n.endsWith(".fasta")) {
                fastaFiles.add(f);
            }
        }
        fastaFiles.sort((a, b) -> a.getName().compareToIgnoreCase(b.getName()));

        int total = fastaFiles.size();
        cb.onMessage("DeepDetect MAG DB prediction: " + total + " genome FASTA files detected in " + originalDbDir.getAbsolutePath());

        final ProgressCallback cbFinal = cb;
        final String regularFinal = regular;
        final int totalFinal = total;

        ExecutorService ex = Executors.newFixedThreadPool(threads);
        AtomicInteger done = new AtomicInteger(0);
        AtomicInteger ok = new AtomicInteger(0);
        AtomicInteger skipped = new AtomicInteger(0);
        AtomicInteger failed = new AtomicInteger(0);

        for (File fasta : fastaFiles) {
            String base = fasta.getName();
            int dot = base.lastIndexOf('.');
            if (dot > 0) base = base.substring(0, dot);
            File outTsv = new File(predictedDir, base + ".predicted.tsv");
            if (outTsv.exists() && outTsv.length() > 0) {
                skipped.incrementAndGet();
                done.incrementAndGet();
                continue;
            }

            DeepDetectParameter par = new DeepDetectParameter(template);
            par.setInput(fasta);
            par.setOutput(outTsv);

            final File fastaFinal = fasta;
            final File outTsvFinal = outTsv;
            final DeepDetectParameter parFinal = par;
            ex.submit(() -> {
                try {
                    int exit = DeepDetectPythonRunner.run(pythonExe, deepDetectMainPy, parFinal, regularFinal, cbFinal);
                    if (exit == 0 && outTsvFinal.exists() && outTsvFinal.length() > 0) {
                        ok.incrementAndGet();
                    } else {
                        failed.incrementAndGet();
                    }
                } catch (Exception e) {
                    failed.incrementAndGet();
                    cbFinal.onMessage("DeepDetect failed for " + fastaFinal.getName() + ": " + e.getMessage());
                } finally {
                    int d = done.incrementAndGet();
                    cbFinal.onProgress(d, totalFinal);
                }
            });
        }

        ex.shutdown();
        boolean finished = ex.awaitTermination(Math.max(1, total), TimeUnit.HOURS);
        if (!finished) {
            cb.onMessage("WARNING: DeepDetect did not finish within timeout; some tasks may still be running.");
        }

        return new Result(predictedDir, total, ok.get(), skipped.get(), failed.get());
    }
}


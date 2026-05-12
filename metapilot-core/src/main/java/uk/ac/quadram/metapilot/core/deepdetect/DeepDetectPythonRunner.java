package uk.ac.quadram.metapilot.core.deepdetect;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.deepDetect.DeepDetectParameter;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

/**
 * Headless DeepDetect runner that invokes the PyTorch Python entrypoint (main.py) directly.
 *
 * <p>Rationale: on HPC/Linux, "activating" a conda env inside Java is brittle.
 * The robust approach is to point {@code PYTHON} to the env's python executable, e.g.:
 * {@code /home/user/miniconda3/envs/deepdetect_torch/bin/python}.
 */
public final class DeepDetectPythonRunner {

    private DeepDetectPythonRunner() {
    }

    /**
     * Run DeepDetect prediction for a single input FASTA.
     *
     * @param pythonExe python executable (preferably from the conda env that contains torch)
     * @param mainPy    DeepDetect python entrypoint (your deepdetect_pred_pytorch/main.py)
     * @param par       DeepDetect parameters (input/output/protease/missed cleavages/len)
     * @param regular   FASTA header regex (optional; if null/blank DeepDetect default is used)
     * @param cb        progress callback (optional)
     */
    public static int run(File pythonExe, File mainPy, DeepDetectParameter par, String regular, ProgressCallback cb)
            throws IOException, InterruptedException {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (pythonExe == null || !pythonExe.exists()) {
            throw new IOException("Python executable not found: " + (pythonExe == null ? "null" : pythonExe.getAbsolutePath()));
        }
        if (mainPy == null || !mainPy.exists()) {
            throw new IOException("DeepDetect main.py not found: " + (mainPy == null ? "null" : mainPy.getAbsolutePath()));
        }
        if (par == null || par.getInput() == null || par.getOutput() == null) {
            throw new IOException("DeepDetectParameter must have input and output set.");
        }

        File out = par.getOutput();
        File outDir = out.getParentFile();
        if (outDir != null) outDir.mkdirs();

        File originalInput = par.getInput();
        File inputToUse = ensureFastaHeadersHaveProteinName(originalInput, outDir != null ? outDir : originalInput.getParentFile());
        boolean isTempInput = inputToUse != null
                && originalInput != null
                && !inputToUse.getAbsolutePath().equals(originalInput.getAbsolutePath())
                && inputToUse.getName().endsWith("_deepDetect_input.faa");

        List<String> cmd = new ArrayList<>();
        cmd.add(pythonExe.getAbsolutePath());
        cmd.add(mainPy.getAbsolutePath());
        cmd.add("--input=" + inputToUse.getAbsolutePath());
        cmd.add("--output=" + out.getAbsolutePath());
        if (regular != null && !regular.trim().isEmpty()) cmd.add("--regular=" + regular.trim());
        cmd.add("--protease=" + par.getProtease());
        cmd.add("--missed_cleavages=" + par.getMissCleavage());
        cmd.add("--min_len=" + par.getMinLength());
        cmd.add("--max_len=" + par.getMaxLength());

        cb.onMessage("Running DeepDetect (python): " + String.join(" ", cmd));

        ProcessBuilder pb = new ProcessBuilder(cmd);
        // Set working directory to the folder containing main.py so relative imports work.
        File wd = mainPy.getParentFile();
        if (wd != null && wd.exists()) pb.directory(wd);
        pb.redirectErrorStream(true);

        Process p = pb.start();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (Thread.currentThread().isInterrupted()) {
                    p.destroyForcibly();
                    throw new InterruptedException("Interrupted while running DeepDetect");
                }
                String t = line.trim();
                if (!t.isEmpty()) cb.onMessage(t);
            }
        }
        int exit = p.waitFor();
        if (exit != 0) cb.onMessage("DeepDetect exit code: " + exit);

        // If we created a temporary sanitized FASTA, remove it after a successful run
        // to avoid doubling the file count in the predicted folder.
        if (exit == 0 && isTempInput) {
            try {
                if (out.exists() && out.length() > 0) {
                    Files.deleteIfExists(inputToUse.toPath());
                }
            } catch (Exception ignored) {
                // Best-effort cleanup only.
            }
        }
        return exit;
    }

    // ---- FASTA preparation (ported from Windows DeepDetectTask, but headless) ----

    private static final String STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY";

    /**
     * DeepDetect requires FASTA headers in the format: {@code >accession proteinName}.
     * This also removes trailing '*' and filters proteins with non-standard amino acids.
     */
    public static File ensureFastaHeadersHaveProteinName(File inputFile, File outputDir) throws IOException {
        if (inputFile == null || !inputFile.exists()) return inputFile;
        File dir = (outputDir != null && outputDir.exists()) ? outputDir : inputFile.getParentFile();
        if (dir == null) dir = new File(System.getProperty("java.io.tmpdir"));

        String baseName = inputFile.getName();
        int lastDot = baseName.lastIndexOf('.');
        String nameNoExt = lastDot > 0 ? baseName.substring(0, lastDot) : baseName;
        File modifiedFile = new File(dir, nameNoExt + "_deepDetect_input.faa");

        int skipped = 0;
        int total = 0;
        try (BufferedReader reader = Files.newBufferedReader(inputFile.toPath(), StandardCharsets.UTF_8);
             PrintWriter writer = new PrintWriter(modifiedFile, StandardCharsets.UTF_8)) {

            String header = null;
            StringBuilder seq = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (header != null) {
                        total++;
                        String cleaned = cleanSequence(seq.toString());
                        if (cleaned != null) {
                            writer.println(header);
                            writer.println(cleaned);
                        } else {
                            skipped++;
                        }
                    }
                    String h = line;
                    if (!line.substring(1).contains(" ")) {
                        String acc = line.substring(1).trim();
                        h = ">" + acc + " " + acc;
                    }
                    header = h;
                    seq = new StringBuilder();
                } else {
                    seq.append(line.trim());
                }
            }
            if (header != null) {
                total++;
                String cleaned = cleanSequence(seq.toString());
                if (cleaned != null) {
                    writer.println(header);
                    writer.println(cleaned);
                } else {
                    skipped++;
                }
            }
        }

        if (skipped > 0) {
            // Don’t fail; just report.
            // (DeepDetect will otherwise error on invalid sequences.)
        }
        return modifiedFile;
    }

    private static String cleanSequence(String seq) {
        if (seq == null) return null;
        if (seq.endsWith("*")) seq = seq.substring(0, seq.length() - 1);
        if (seq.isEmpty()) return null;
        for (int i = 0; i < seq.length(); i++) {
            if (STANDARD_AA.indexOf(seq.charAt(i)) < 0) return null;
        }
        return seq;
    }
}


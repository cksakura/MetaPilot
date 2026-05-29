package uk.ac.quadram.metapilot.core.diann;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

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
 * Headless DIA-NN process runner (Linux/HPC safe).
 *
 * <p>Always writes the DIA-NN argument string to a {@code <output>.cfg} file next to the output
 * and invokes DIA-NN with {@code --cfg <cfg>} (matching the Windows GUI convention and the
 * canonical HPC invocation pattern):
 * <pre>
 *   export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
 *   ./diann-linux --cfg /path/to/output.parquet.cfg
 * </pre>
 *
 * <p>Why cfg-mode instead of inline args:
 * <ul>
 *   <li>DEEP_PER_FILE search commands include one {@code --f} per raw file plus several long
 *       {@code --lib}/{@code --out} paths; on shared HPC nodes this can approach the kernel
 *       {@code ARG_MAX} limit and is harder to debug than a saved cfg file.</li>
 *   <li>The Windows GUI workflow already uses cfg-mode, so behavior stays consistent across OSes
 *       and reproducibility tooling that inspects {@code .cfg} files works the same way on HPC.</li>
 *   <li>The cfg file is left on disk after the run as a permanent audit/reproducibility artifact.</li>
 * </ul>
 */
public final class DiaNNHeadlessRunner {

    private DiaNNHeadlessRunner() {
    }

    public static int run(File diannExe, String argsString, File workingDir, ProgressCallback cb)
            throws IOException, InterruptedException {
        return run(diannExe, argsString, workingDir, cb, null);
    }

    /**
     * Run DIA-NN with optional conda prefix lib injection.
     *
     * <p>On Linux, some DIA-NN builds require a newer libstdc++ (GLIBCXX) than the host provides.
     * If a conda environment is available, prepending {@code <condaPrefix>/lib} to {@code LD_LIBRARY_PATH}
     * usually fixes it.
     *
     * @param condaPrefix optional conda environment prefix; if null, falls back to {@code $CONDA_PREFIX} when present
     */
    public static int run(File diannExe, String argsString, File workingDir, ProgressCallback cb, File condaPrefix)
            throws IOException, InterruptedException {
        if (cb == null) cb = new ProgressCallback.NoOp();

        if (diannExe == null || !diannExe.exists()) {
            throw new IOException("DIA-NN executable not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));
        }
        if (argsString == null || argsString.isBlank()) {
            throw new IOException("DIA-NN args string is empty");
        }

        // 1) Decide where to write the cfg file: <output>.cfg next to the parquet/speclib if we can
        //    parse it from the args, otherwise drop a timestamped file in the working directory.
        File cfgFile = resolveCfgFile(argsString, workingDir, diannExe);
        File cfgParent = cfgFile.getParentFile();
        if (cfgParent != null && !cfgParent.isDirectory()) {
            cfgParent.mkdirs();
        }
        // Write the args verbatim. DIA-NN accepts the entire argument string on a single line
        // and tolerates trailing whitespace; we write what the GUI's DiaNNTask writes.
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(cfgFile.toPath(), StandardCharsets.UTF_8))) {
            w.print(argsString);
        }
        cb.onMessage("DIA-NN cfg written: " + cfgFile.getAbsolutePath());

        // 2) Build invocation: <diann-linux> --cfg <cfgFile>
        List<String> cmd = new ArrayList<>();
        cmd.add(diannExe.getAbsolutePath());
        cmd.add("--cfg");
        cmd.add(cfgFile.getAbsolutePath());

        cb.onMessage("Running DIA-NN: " + String.join(" ", cmd));

        ProcessBuilder pb = new ProcessBuilder(cmd);
        File parent = diannExe.getParentFile();
        if (parent != null && parent.exists()) {
            pb.directory(parent);
        } else if (workingDir != null) {
            pb.directory(workingDir);
        }
        pb.redirectErrorStream(true);

        // 3) Linux-only: inject conda lib into LD_LIBRARY_PATH so DIA-NN can find a newer libstdc++/libgomp.
        //    Mirrors the canonical HPC recipe: `export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"`.
        String os = System.getProperty("os.name", "").toLowerCase();
        if (os.contains("linux")) {
            File prefix = condaPrefix;
            if (prefix == null) {
                String envPrefix = System.getenv("CONDA_PREFIX");
                if (envPrefix != null && !envPrefix.isBlank()) {
                    prefix = new File(envPrefix.trim());
                }
            }
            if (prefix != null) {
                File lib = new File(prefix, "lib");
                if (lib.isDirectory()) {
                    String existing = pb.environment().getOrDefault("LD_LIBRARY_PATH", "");
                    String libPath = lib.getAbsolutePath();
                    String merged = existing == null || existing.isBlank()
                            ? libPath
                            : (libPath + ":" + existing);
                    pb.environment().put("LD_LIBRARY_PATH", merged);
                    cb.onMessage("DIA-NN: using LD_LIBRARY_PATH prefix from conda: " + libPath);
                }
            }
        }

        Process p = pb.start();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (Thread.currentThread().isInterrupted()) {
                    p.destroyForcibly();
                    throw new InterruptedException("Interrupted while running DIA-NN");
                }
                String t = line.trim();
                if (!t.isEmpty()) cb.onMessage(t);
            }
        }

        int exit = p.waitFor();
        if (exit != 0) {
            cb.onMessage("DIA-NN exit code: " + exit + " (cfg: " + cfgFile.getAbsolutePath() + ")");
        }
        return exit;
    }

    // ----------------------- cfg path resolution -----------------------

    /**
     * Locate the most appropriate cfg file path for this DIA-NN invocation.
     *
     * <p>Priority:
     * <ol>
     *   <li>If args contain {@code --out "<path>"} (parquet output, e.g. search), write {@code <path>.cfg}.</li>
     *   <li>Else if args contain {@code --out-lib "<path>"} (speclib output, e.g. generateFasta), write {@code <path>.cfg}.</li>
     *   <li>Else write a timestamped fallback in the working directory (or DIA-NN binary folder).</li>
     * </ol>
     */
    private static File resolveCfgFile(String args, File workingDir, File diannExe) {
        String outPath = extractArgValue(args, "--out");
        if (outPath == null) {
            outPath = extractArgValue(args, "--out-lib");
        }
        if (outPath != null && !outPath.isBlank()) {
            return new File(outPath + ".cfg");
        }
        File baseDir = (workingDir != null && workingDir.exists()) ? workingDir
                : (diannExe.getParentFile() != null ? diannExe.getParentFile() : new File("."));
        return new File(baseDir, "diann_run_" + System.currentTimeMillis() + ".cfg");
    }

    /**
     * Extract a single argument's value from a DIA-NN argument string.
     * Handles both quoted ({@code --out "/path with spaces/file"}) and unquoted ({@code --out /path/file}).
     * Returns {@code null} if the key is not found.
     */
    static String extractArgValue(String args, String key) {
        if (args == null || key == null) return null;
        String needle = key + " ";
        // Search for the key followed by a space; tolerate it appearing multiple times by taking the first.
        int idx = args.indexOf(needle);
        if (idx < 0) return null;
        // Make sure it's a whole-token match (not a suffix of another flag like "--no-out").
        if (idx > 0) {
            char before = args.charAt(idx - 1);
            if (before != ' ' && before != '\t' && before != '\n') return null;
        }
        int s = idx + needle.length();
        if (s >= args.length()) return null;
        if (args.charAt(s) == '"') {
            int e = args.indexOf('"', s + 1);
            if (e < 0) return null;
            return args.substring(s + 1, e);
        }
        int e = s;
        while (e < args.length()) {
            char c = args.charAt(e);
            if (c == ' ' || c == '\t' || c == '\n') break;
            e++;
        }
        return args.substring(s, e);
    }
}

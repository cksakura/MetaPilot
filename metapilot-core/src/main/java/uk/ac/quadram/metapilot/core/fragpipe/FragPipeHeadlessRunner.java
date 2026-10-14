package uk.ac.quadram.metapilot.core.fragpipe;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * OS-agnostic runner for FragPipe headless mode.
 *
 * Command shape (FragPipe supports this):
 *   fragpipe --headless --workflow <workflow> --manifest <manifest> --workdir <dir> --threads <n>
 */
public final class FragPipeHeadlessRunner {

    private FragPipeHeadlessRunner() {
    }

    public static int run(File fragpipeExe, File workflowFile, File manifestFile, File workDir, int threads, ProgressCallback cb)
            throws IOException, InterruptedException {
        return run(fragpipeExe, workflowFile, manifestFile, workDir, threads, 0, cb);
    }

    /**
     * @param ramGb maximum memory FragPipe may hand to the tools it launches, in
     *              GB, passed through as {@code --ram}. Pass 0 to omit the flag and
     *              let FragPipe decide, which is its documented default. FragPipe's
     *              own choice is the whole machine's RAM, so every concurrent
     *              instance would claim all of it; callers running several searches
     *              at once must divide a budget between them instead.
     */
    public static int run(File fragpipeExe, File workflowFile, File manifestFile, File workDir, int threads, int ramGb,
            ProgressCallback cb) throws IOException, InterruptedException {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (fragpipeExe == null || !fragpipeExe.exists()) {
            throw new IOException("FragPipe executable not found: " + (fragpipeExe == null ? "null" : fragpipeExe.getAbsolutePath()));
        }
        if (workflowFile == null || !workflowFile.exists()) {
            throw new IOException("Workflow file not found: " + (workflowFile == null ? "null" : workflowFile.getAbsolutePath()));
        }
        if (manifestFile == null || !manifestFile.exists()) {
            throw new IOException("Manifest file not found: " + (manifestFile == null ? "null" : manifestFile.getAbsolutePath()));
        }
        if (workDir == null) {
            throw new IOException("workDir is null");
        }
        workDir.mkdirs();

        List<String> fragpipeArgs = new ArrayList<>();
        fragpipeArgs.add("--headless");
        fragpipeArgs.add("--workflow");
        fragpipeArgs.add(workflowFile.getAbsolutePath());
        fragpipeArgs.add("--manifest");
        fragpipeArgs.add(manifestFile.getAbsolutePath());
        fragpipeArgs.add("--workdir");
        fragpipeArgs.add(workDir.getAbsolutePath());
        fragpipeArgs.add("--threads");
        fragpipeArgs.add(String.valueOf(Math.max(1, threads)));
        if (ramGb > 0) {
            fragpipeArgs.add("--ram");
            fragpipeArgs.add(String.valueOf(ramGb));
        }

        List<String> cmd = new ArrayList<>();
        if (isWindowsBatchFile(fragpipeExe)) {
            cmd.add("cmd.exe");
            cmd.add("/c");
            cmd.add("call");
            cmd.add(fragpipeExe.getAbsolutePath());
        } else {
            cmd.add(fragpipeExe.getAbsolutePath());
        }
        cmd.addAll(fragpipeArgs);

        cb.onMessage("Running FragPipe: " + String.join(" ", cmd));

        ProcessBuilder pb = new ProcessBuilder(cmd);
        // Similar to Windows GUI behavior (cd to FragPipe bin folder).
        File parent = fragpipeExe.getParentFile();
        if (parent != null && parent.exists()) {
            pb.directory(parent);
        }
        pb.redirectErrorStream(true);
        ensureJavaHome(pb, fragpipeExe, cb);

        Process p = pb.start();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (Thread.currentThread().isInterrupted()) {
                    p.destroyForcibly();
                    throw new InterruptedException("Interrupted while running FragPipe");
                }
                String t = line.trim();
                if (!t.isEmpty()) cb.onMessage(t);
            }
        }

        int exit = p.waitFor();
        if (exit != 0) {
            cb.onMessage("FragPipe exit code: " + exit);
        }
        return exit;
    }

    /**
     * Guarantees the launcher can find a JVM.
     * <p>
     * FragPipe's {@code bin/fragpipe.bat} (and the {@code bin/fragpipe} shell
     * script) locate Java from {@code JAVA_HOME}, falling back to a bare
     * {@code java} on the PATH. Neither is guaranteed for a child process:
     * MetaPilot started from a desktop shortcut, an IDE or a packaged launcher
     * inherits an environment where {@code java} is not on the PATH and
     * {@code JAVA_HOME} is unset. The launcher then reports "JAVA_HOME is not
     * set and no 'java' command could be found in your PATH" and exits with the
     * Windows command-not-found code 9009, which surfaces as an opaque FragPipe
     * failure.
     * <p>
     * MetaPilot is itself a JVM, so a valid runtime is always available. An
     * inherited {@code JAVA_HOME} is left untouched (the user configured it
     * deliberately); otherwise FragPipe's own bundled {@code jre} wins, and the
     * runtime running MetaPilot is the fallback. FragPipe requires Java 9+ and
     * MetaPilot targets 17, so reusing this runtime is safe.
     */
    static void ensureJavaHome(ProcessBuilder pb, File fragpipeExe, ProgressCallback cb) {
        String inherited = pb.environment().get("JAVA_HOME");
        if (inherited != null && !inherited.isBlank()) {
            return;
        }
        File javaHome = resolveJavaHome(fragpipeExe, System.getProperty("java.home"), System.getProperty("os.name"));
        if (javaHome == null) {
            // Nothing better to offer than the launcher's own PATH lookup.
            return;
        }
        pb.environment().put("JAVA_HOME", javaHome.getAbsolutePath());
        if (cb != null) {
            cb.onMessage("JAVA_HOME not set for FragPipe; using " + javaHome.getAbsolutePath());
        }
    }

    /**
     * Prefers the JRE bundled next to FragPipe ({@code <root>/jre}), then the
     * runtime running MetaPilot. Returns null when neither contains a launcher.
     */
    static File resolveJavaHome(File fragpipeExe, String metapilotJavaHome, String osName) {
        File bundled = bundledFragpipeJre(fragpipeExe);
        if (hasJavaLauncher(bundled, osName)) {
            return bundled;
        }
        if (metapilotJavaHome != null && !metapilotJavaHome.isBlank()) {
            File own = new File(metapilotJavaHome);
            if (hasJavaLauncher(own, osName)) {
                return own;
            }
        }
        return null;
    }

    /**
     * The configured entry is {@code <root>/bin/fragpipe.bat}, so the install
     * root is two levels up and the bundled runtime sits at {@code <root>/jre}.
     */
    private static File bundledFragpipeJre(File fragpipeExe) {
        if (fragpipeExe == null) {
            return null;
        }
        File binDir = fragpipeExe.getParentFile();
        if (binDir == null) {
            return null;
        }
        File root = binDir.getParentFile();
        return root == null ? null : new File(root, "jre");
    }

    private static boolean hasJavaLauncher(File javaHome, String osName) {
        if (javaHome == null) {
            return false;
        }
        boolean windows = osName != null && osName.toLowerCase(Locale.ROOT).contains("win");
        Path launcher = javaHome.toPath().resolve("bin").resolve(windows ? "java.exe" : "java");
        return Files.isRegularFile(launcher);
    }

    private static boolean isWindowsBatchFile(File file) {
        String name = file.getName().toLowerCase();
        return name.endsWith(".bat") || name.endsWith(".cmd");
    }
}


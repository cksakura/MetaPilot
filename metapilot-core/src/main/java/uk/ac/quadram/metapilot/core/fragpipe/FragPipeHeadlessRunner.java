package uk.ac.quadram.metapilot.core.fragpipe;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

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

    private static boolean isWindowsBatchFile(File file) {
        String name = file.getName().toLowerCase();
        return name.endsWith(".bat") || name.endsWith(".cmd");
    }
}


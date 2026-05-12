package uk.ac.quadram.metapilot.core.resources;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

final class ToolLocator {

    private ToolLocator() {
    }

    static Optional<Path> findOnPath(List<String> commandNames, ProgressCallback cb) {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (commandNames == null || commandNames.isEmpty()) return Optional.empty();

        String os = System.getProperty("os.name", "").toLowerCase(Locale.ROOT);
        boolean isWindows = os.contains("win");
        String whichCmd = isWindows ? "where" : "which";

        for (String name : commandNames) {
            if (name == null || name.trim().isEmpty()) continue;
            try {
                Process p = new ProcessBuilder(whichCmd, name).redirectErrorStream(true).start();
                try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
                    String firstLine = br.readLine();
                    int exit = p.waitFor();
                    if (exit == 0 && firstLine != null && !firstLine.trim().isEmpty()) {
                        Path found = Paths.get(firstLine.trim());
                        cb.onMessage("Auto-detected " + name + ": " + found);
                        return Optional.of(found);
                    }
                }
            } catch (IOException | InterruptedException ignored) {
                if (ignored instanceof InterruptedException) {
                    Thread.currentThread().interrupt();
                    return Optional.empty();
                }
            }
        }
        return Optional.empty();
    }
}


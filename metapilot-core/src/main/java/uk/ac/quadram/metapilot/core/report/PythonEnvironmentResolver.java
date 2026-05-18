package uk.ac.quadram.metapilot.core.report;

import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.resources.ToolKey;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Optional;
import java.util.concurrent.TimeUnit;

/**
 * Shared Python resolver for GUI and headless report generation.
 */
public final class PythonEnvironmentResolver {

    private static final long VALIDATION_TIMEOUT_SECONDS = 20;
    private static final String PACKAGE_CHECK_SCRIPT =
            "import pandas, numpy, plotly; print('OK')";
    private static final String VERSION_CHECK_SCRIPT =
            "import sys; print(sys.executable); print(sys.version.split()[0])";

    private PythonEnvironmentResolver() {
    }

    public static final class Request {
        public String explicitPython;
        public MetaResourceParameters resourceParameters;
        public MetaPilotResourceConfig resourceConfig;
        public boolean validateReportPackages = true;
        public ProgressCallback progressCallback;
    }

    public static final class Result {
        public String resolvedCommand;
        public String source;
        public String versionText;
        public boolean available;
        public boolean valid;
        public String failureMessage;

        public File asFileIfPath() {
            if (resolvedCommand == null || resolvedCommand.isBlank()) {
                return null;
            }
            File file = new File(resolvedCommand);
            return file.isAbsolute() || resolvedCommand.contains(File.separator) ? file : null;
        }
    }

    public static Result resolve(String explicitPython) {
        Request request = new Request();
        request.explicitPython = explicitPython;
        return resolve(request);
    }

    public static Result resolve(Request request) {
        if (request == null) {
            request = new Request();
        }
        ProgressCallback cb = request.progressCallback != null ? request.progressCallback : new ProgressCallback.NoOp();
        List<Candidate> candidates = collectCandidates(request, cb);
        List<String> failures = new ArrayList<>();

        for (Candidate candidate : candidates) {
            if (candidate.command == null || candidate.command.isBlank()) {
                continue;
            }
            Result result = validateCandidate(candidate, request.validateReportPackages);
            if (result.valid) {
                return result;
            }
            if (result.failureMessage != null && !result.failureMessage.isBlank()) {
                failures.add(candidate.source + "=" + candidate.command + " (" + result.failureMessage + ")");
            }
        }

        Result result = new Result();
        result.available = false;
        result.valid = false;
        result.failureMessage = failures.isEmpty()
                ? "No Python candidate was configured or found on PATH"
                : String.join("; ", failures);
        return result;
    }

    public static File inferCondaPrefixFromPython(File pythonExe) {
        if (pythonExe == null) {
            return null;
        }
        try {
            File bin = pythonExe.getParentFile();
            if (bin == null) {
                return null;
            }
            File prefix = bin.getParentFile();
            if (prefix == null) {
                return null;
            }
            String binName = bin.getName().toLowerCase(Locale.ROOT);
            if (!"bin".equals(binName) && !"scripts".equals(binName)) {
                return null;
            }
            File lib = new File(prefix, "lib");
            return lib.isDirectory() ? prefix : null;
        } catch (Exception e) {
            return null;
        }
    }

    private static List<Candidate> collectCandidates(Request request, ProgressCallback cb) {
        List<Candidate> candidates = new ArrayList<>();
        addCandidate(candidates, "explicit", request.explicitPython);
        addCandidate(candidates, "METAPILOT_PYTHON", System.getenv("METAPILOT_PYTHON"));

        if (request.resourceConfig != null) {
            Optional<Path> resolved = request.resourceConfig.resolveToolPath(ToolKey.PYTHON, cb);
            if (resolved.isPresent()) {
                addCandidate(candidates, "resource-config:PYTHON", resolved.get().toString());
            }
            request.resourceConfig.getRaw("PYTHON")
                    .ifPresent(value -> addCandidate(candidates, "resource-config:PYTHON(raw)", value));
        }

        if (request.resourceParameters != null) {
            try {
                addCandidate(candidates, "resource-parameters:PYTHON",
                        request.resourceParameters.getToolPath(MetaResourceParameters.ToolType.PYTHON));
            } catch (Exception ignored) {
            }
        }

        addCandidate(candidates, "PATH:python3", "python3");
        addCandidate(candidates, "PATH:python", "python");
        return candidates;
    }

    private static void addCandidate(List<Candidate> candidates, String source, String command) {
        if (command == null || command.trim().isEmpty()) {
            return;
        }
        String trimmed = command.trim();
        for (Candidate existing : candidates) {
            if (existing.command.equals(trimmed)) {
                return;
            }
        }
        candidates.add(new Candidate(source, trimmed));
    }

    private static Result validateCandidate(Candidate candidate, boolean validateReportPackages) {
        Result result = new Result();
        result.resolvedCommand = candidate.command;
        result.source = candidate.source;

        CommandResult version = runPython(candidate.command, VERSION_CHECK_SCRIPT);
        if (!version.success) {
            result.available = false;
            result.valid = false;
            result.failureMessage = version.message;
            return result;
        }

        result.available = true;
        result.versionText = version.output.trim();
        if (!validateReportPackages) {
            result.valid = true;
            return result;
        }

        CommandResult packages = runPython(candidate.command, PACKAGE_CHECK_SCRIPT);
        result.valid = packages.success;
        if (!packages.success) {
            result.failureMessage = packages.message;
        }
        return result;
    }

    private static CommandResult runPython(String command, String script) {
        ProcessBuilder pb = new ProcessBuilder(command, "-c", script);
        pb.redirectErrorStream(true);
        StringBuilder output = new StringBuilder();
        try {
            Process process = pb.start();
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    output.append(line).append(System.lineSeparator());
                }
            }
            boolean finished = process.waitFor(VALIDATION_TIMEOUT_SECONDS, TimeUnit.SECONDS);
            if (!finished) {
                process.destroyForcibly();
                return CommandResult.fail("Python validation timed out");
            }
            int exit = process.exitValue();
            if (exit != 0) {
                return CommandResult.fail("exit=" + exit + ", output=" + output.toString().trim());
            }
            return CommandResult.ok(output.toString());
        } catch (IOException e) {
            return CommandResult.fail(e.getMessage());
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return CommandResult.fail("interrupted");
        }
    }

    private static final class Candidate {
        final String source;
        final String command;

        Candidate(String source, String command) {
            this.source = source;
            this.command = command;
        }
    }

    private static final class CommandResult {
        final boolean success;
        final String output;
        final String message;

        private CommandResult(boolean success, String output, String message) {
            this.success = success;
            this.output = output;
            this.message = message;
        }

        static CommandResult ok(String output) {
            return new CommandResult(true, output, null);
        }

        static CommandResult fail(String message) {
            return new CommandResult(false, "", message);
        }
    }
}

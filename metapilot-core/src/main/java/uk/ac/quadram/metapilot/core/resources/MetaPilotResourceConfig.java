package uk.ac.quadram.metapilot.core.resources;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;

/**
 * Headless resource configuration for MetaPilot CLI.
 *
 * Backwards compatible with the Windows-style file (KEY=VALUE), but adds optional:
 * - SYSTEM=linux/windows
 * - RESOURCE_AUTO_DETECTED=true/false
 * - <TOOL>_MODE=path|which  (currently used for CLI)
 * - <TOOL>_ENABLED=true/false (ignored unless explicitly checked)
 */
public final class MetaPilotResourceConfig {

    public enum ToolMode {
        PATH, WHICH
    }

    private final Path sourceFile;
    private final Map<String, String> raw;
    private final EnumMap<ToolKey, ToolSpec> tools = new EnumMap<>(ToolKey.class);

    private final String system;
    private final boolean autoDetected;
    private final Path basePath;

    private MetaPilotResourceConfig(Path sourceFile, Map<String, String> raw) {
        this.sourceFile = sourceFile;
        this.raw = raw;
        this.system = raw.getOrDefault("SYSTEM", "").trim();
        this.autoDetected = parseBool(raw.get("RESOURCE_AUTO_DETECTED"), false);
        this.basePath = parsePathOrNull(raw.get("BASE_PATH"));

        for (ToolKey k : ToolKey.values()) {
            String v = raw.getOrDefault(k.key(), "").trim();
            ToolMode mode = parseMode(raw.get(k.key() + "_MODE"));
            boolean enabled = parseBool(raw.get(k.key() + "_ENABLED"), true);
            tools.put(k, new ToolSpec(k, mode, enabled, v));
        }
    }

    public static MetaPilotResourceConfig load(Path file) throws IOException {
        Objects.requireNonNull(file, "file");
        Map<String, String> map = ParamsFile.readKeyValue(file);
        return new MetaPilotResourceConfig(file, map);
    }

    public Path getSourceFile() {
        return sourceFile;
    }

    public String getSystem() {
        return system;
    }

    public boolean isAutoDetected() {
        return autoDetected;
    }

    public Optional<Path> getBasePath() {
        return Optional.ofNullable(basePath);
    }

    public Optional<String> getRaw(String key) {
        if (key == null) return Optional.empty();
        String v = raw.get(key);
        if (v == null) return Optional.empty();
        String t = v.trim();
        return t.isEmpty() ? Optional.empty() : Optional.of(t);
    }

    public ToolSpec tool(ToolKey key) {
        return tools.get(key);
    }

    public Optional<Path> resolveToolPath(ToolKey key, ProgressCallback cb) {
        if (cb == null) cb = new ProgressCallback.NoOp();
        ToolSpec spec = tool(key);
        if (spec == null || !spec.enabled) return Optional.empty();

        // WHICH mode: treat value as a command name hint (NOT a filesystem path).
        // If value is empty, use default names.
        if (spec.mode == ToolMode.WHICH) {
            List<String> names = new ArrayList<>();
            if (spec.value != null && !spec.value.trim().isEmpty()) {
                // allow comma/space separated command names (e.g. "python,python3")
                String[] parts = spec.value.trim().split("[,\\s]+");
                for (String p : parts) {
                    if (p != null) {
                        String t = p.trim();
                        if (!t.isEmpty()) names.add(t);
                    }
                }
            }
            if (names.isEmpty()) names = defaultWhichNames(key);
            return ToolLocator.findOnPath(names, cb);
        }

        // PATH mode: explicit value must exist on disk (supporting ~/ expansion)
        if (spec.value != null && !spec.value.trim().isEmpty()) {
            return Optional.of(expandHome(Paths.get(spec.value.trim())));
        }

        // For now, PATH mode requires explicit values (we keep BASE_PATH for future expansion)
        return Optional.empty();
    }

    public List<String> validateMinimalLinuxCli() {
        List<String> issues = new ArrayList<>();

        // Required for everything
        checkResolvableExists(ToolKey.JAVA, issues);
        checkResolvableExists(ToolKey.PYTHON, issues);
        checkResolvableExists(ToolKey.FUNCDB, issues);

        // Required tools you said we need now
        checkResolvableExists(ToolKey.FRAGPIPE, issues);
        checkResolvableExists(ToolKey.FRAGPIPE_TOOLS, issues);
        checkResolvableExists(ToolKey.DIANN, issues);
        checkResolvableExists(ToolKey.DEEPDETECT, issues);

        return issues;
    }

    private void checkResolvableExists(ToolKey key, List<String> issues) {
        ToolSpec spec = tool(key);
        if (spec == null || !spec.enabled) {
            issues.add(key.key() + " is disabled but required");
            return;
        }

        Optional<Path> pOpt = resolveToolPath(key, new ProgressCallback.NoOp());
        if (pOpt.isEmpty()) {
            if (spec.mode == ToolMode.WHICH) {
                issues.add(key.key() + " not found on PATH (mode=which)");
            } else {
                issues.add(key.key() + " is missing");
            }
            return;
        }

        Path p = pOpt.get();
        if (!Files.exists(p)) {
            issues.add(key.key() + " does not exist: " + p);
        }
    }

    private static ToolMode parseMode(String s) {
        if (s == null) return ToolMode.PATH;
        String t = s.trim().toUpperCase(Locale.ROOT);
        if (t.equals("WHICH")) return ToolMode.WHICH;
        return ToolMode.PATH;
    }

    private static boolean parseBool(String s, boolean def) {
        if (s == null) return def;
        String t = s.trim().toLowerCase(Locale.ROOT);
        if (t.isEmpty()) return def;
        if (t.equals("true") || t.equals("1") || t.equals("yes") || t.equals("y")) return true;
        if (t.equals("false") || t.equals("0") || t.equals("no") || t.equals("n")) return false;
        return def;
    }

    private static Path parsePathOrNull(String s) {
        if (s == null) return null;
        String t = s.trim();
        if (t.isEmpty()) return null;
        return expandHome(Paths.get(t));
    }

    private static Path expandHome(Path p) {
        if (p == null) return null;
        String s = p.toString();
        if (s.equals("~")) {
            return Paths.get(System.getProperty("user.home", ""));
        }
        if (s.startsWith("~/") || s.startsWith("~\\")) {
            String home = System.getProperty("user.home", "");
            if (home == null || home.isBlank()) return p;
            String rest = s.substring(2);
            return Paths.get(home).resolve(rest);
        }
        return p;
    }

    private static List<String> defaultWhichNames(ToolKey key) {
        switch (key) {
            case JAVA:
                return Collections.singletonList("java");
            case PYTHON:
                return List.of("python", "python3");
            case DIANN:
                return List.of("diann", "DiaNN", "diann-linux");
            case FRAGPIPE:
                return List.of("fragpipe");
            case DEEPDETECT:
                return List.of("deepdetect", "DeepDetect");
            default:
                return Collections.emptyList();
        }
    }

    public static final class ToolSpec {
        public final ToolKey key;
        public final ToolMode mode;
        public final boolean enabled;
        public final String value;

        ToolSpec(ToolKey key, ToolMode mode, boolean enabled, String value) {
            this.key = Objects.requireNonNull(key, "key");
            this.mode = mode == null ? ToolMode.PATH : mode;
            this.enabled = enabled;
            this.value = value == null ? "" : value;
        }

        public Optional<File> asFile() {
            if (value == null || value.trim().isEmpty()) return Optional.empty();
            return Optional.of(new File(value.trim()));
        }
    }
}


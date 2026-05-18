package uk.ac.quadram.metapilot.core.fragpipe;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Minimal headless reader/writer for FragPipe .workflow files.
 *
 * Format is basically KEY=VALUE with optional # comments. We preserve insertion order.
 */
public final class FragpipeWorkflowFile {

    private final LinkedHashMap<String, String> kv;

    private FragpipeWorkflowFile(LinkedHashMap<String, String> kv) {
        this.kv = kv;
    }

    public static FragpipeWorkflowFile parse(Path workflow) throws IOException {
        LinkedHashMap<String, String> map = new LinkedHashMap<>();
        try (BufferedReader br = Files.newBufferedReader(workflow, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) continue;
                int eq = t.indexOf('=');
                if (eq < 1) continue;
                String key = t.substring(0, eq).trim();
                String value = t.substring(eq + 1).trim();
                if (!key.isEmpty()) {
                    map.put(key, value);
                }
            }
        }
        return new FragpipeWorkflowFile(map);
    }

    public FragpipeWorkflowFile set(String key, String value) {
        if (key == null || key.trim().isEmpty()) return this;
        kv.put(key.trim(), value == null ? "" : value);
        return this;
    }

    public FragpipeWorkflowFile setDatabase(String fastaPath) {
        return set("database.db-path", FragpipeWorkflowPaths.formatDatabasePath(fastaPath));
    }

    public String get(String key) {
        return kv.get(key);
    }

    public void write(Path out) throws IOException {
        Files.createDirectories(out.getParent());
        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(out, StandardCharsets.UTF_8))) {
            for (Map.Entry<String, String> e : kv.entrySet()) {
                String value = e.getValue() == null ? "" : e.getValue();
                if ("database.db-path".equals(e.getKey())) {
                    value = FragpipeWorkflowPaths.formatDatabasePath(value);
                }
                pw.println(e.getKey() + "=" + value);
            }
        }
    }

}


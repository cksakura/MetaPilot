package uk.ac.quadram.metapilot.core.resources;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Minimal .params reader: KEY=VALUE with # comments.
 */
final class ParamsFile {

    private ParamsFile() {
    }

    static Map<String, String> readKeyValue(Path file) throws IOException {
        Map<String, String> out = new LinkedHashMap<>();
        try (BufferedReader br = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                int eq = line.indexOf('=');
                if (eq < 1) continue;
                String key = line.substring(0, eq).trim();
                String value = line.substring(eq + 1).trim();
                if (!key.isEmpty()) {
                    out.put(key, value);
                }
            }
        }
        return out;
    }
}


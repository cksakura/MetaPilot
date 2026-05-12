package uk.ac.quadram.metapilot.cli;

import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;

import java.io.BufferedWriter;
import java.io.IOException;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;

final class HpcSiteConfig {

    private HpcSiteConfig() {
    }

    static void bindInstallHome() {
        System.setProperty("metapilot.install.dir", installHome().toString());
    }

    static Path installHome() {
        String env = System.getenv("METAPILOT_HOME");
        if (env != null && !env.isBlank()) {
            return Paths.get(env.trim()).toAbsolutePath().normalize();
        }
        return jarInstallDir();
    }

    static Path siteConfigDir() {
        return installHome().resolve("config");
    }

    static Path resolveResourcesConfigPath(String explicit) {
        if (explicit != null && !explicit.isBlank()) {
            return Paths.get(explicit.trim());
        }
        String env = System.getenv("METAPILOT_RESOURCES");
        if (env != null && !env.isBlank()) {
            return Paths.get(env.trim());
        }
        return siteConfigDir().resolve("meta_resources.params");
    }

    static Path resolveMagDbConfigPath(String explicit) {
        if (explicit != null && !explicit.isBlank()) {
            return Paths.get(explicit.trim());
        }
        String env = System.getenv("METAPILOT_MAGDB");
        if (env != null && !env.isBlank()) {
            return Paths.get(env.trim());
        }
        return siteConfigDir().resolve(MagDbConfigIO.versionFile);
    }

    static List<String> ensureSiteConfig(boolean verbose) throws IOException {
        Path dir = siteConfigDir();
        Files.createDirectories(dir);
        List<String> created = new ArrayList<>();

        Path resources = dir.resolve("meta_resources.params");
        if (!Files.isRegularFile(resources)) {
            writeResourcesTemplate(resources);
            created.add(resources.toString());
            if (verbose) {
                System.out.println("Created site resources template: " + resources);
            }
        }

        Path magdb = dir.resolve(MagDbConfigIO.versionFile);
        if (!Files.isRegularFile(magdb)) {
            MagDbConfigIO.exportBlank(magdb.toFile());
            created.add(magdb.toString());
            if (verbose) {
                System.out.println("Created blank MAG registry: " + magdb);
            }
        }
        return created;
    }

    static LinkedHashMap<String, String> defaultLinuxResourceTemplate() {
        LinkedHashMap<String, String> out = new LinkedHashMap<>();
        out.put("SYSTEM", "linux");
        out.put("RESOURCE_AUTO_DETECTED", "true");
        out.put("BASE_PATH", "~/tools");
        out.put("JAVA", "java");
        out.put("PYTHON", "python");
        out.put("DIANN_MODE", "path");
        out.put("DIANN", "~/tools/diann-2.5.0/diann-linux");
        out.put("FRAGPIPE_MODE", "path");
        out.put("FRAGPIPE", "~/tools/fragpipe-24.0/bin/fragpipe");
        out.put("FRAGPIPE_TOOLS", "~/tools/fragpipe-24.0/tools");
        out.put("DEEPDETECT_MODE", "path");
        out.put("DEEPDETECT", "~/tools/deepdetect/deepdetect");
        out.put("FUNCDB", "~/tools/func_def.db");
        return out;
    }

    private static Path jarInstallDir() {
        try {
            URI uri = MetaPilotCliApp.class.getProtectionDomain().getCodeSource().getLocation().toURI();
            Path codeLocation = Paths.get(uri).toAbsolutePath().normalize();
            if (Files.isDirectory(codeLocation)) {
                return Paths.get("").toAbsolutePath().normalize();
            }
            Path parent = codeLocation.getParent();
            return parent != null ? parent : Paths.get("").toAbsolutePath().normalize();
        } catch (Exception ignored) {
            return Paths.get("").toAbsolutePath().normalize();
        }
    }

    private static void writeResourcesTemplate(Path file) throws IOException {
        LinkedHashMap<String, String> kv = defaultLinuxResourceTemplate();
        List<String> lines = new ArrayList<>();
        lines.add("# MetaPilot resource configuration");
        lines.add("# Created by metapilot-cli init at " + LocalDateTime.now());
        lines.add("");
        ArrayList<String> keys = new ArrayList<>(kv.keySet());
        keys.sort(Comparator.comparingInt(HpcSiteConfig::resourceKeyOrder).thenComparing(String::compareToIgnoreCase));
        for (String key : keys) {
            lines.add(key + "=" + kv.get(key));
        }
        try (BufferedWriter writer = Files.newBufferedWriter(file, StandardCharsets.UTF_8)) {
            for (String line : lines) {
                writer.write(line);
                writer.newLine();
            }
        }
    }

    private static int resourceKeyOrder(String key) {
        if (key == null) {
            return 1000;
        }
        return switch (key) {
            case "SYSTEM" -> 0;
            case "RESOURCE_AUTO_DETECTED" -> 1;
            case "BASE_PATH" -> 2;
            default -> 10;
        };
    }
}

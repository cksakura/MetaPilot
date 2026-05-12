package uk.ac.quadram.metapilot.task.mag;

import uk.ac.quadram.metapilot.core.resources.MetaPilotInstallPaths;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Resolves {@link MagDbItem} entries from the site MAG registry JSON.
 */
public final class MagDbItemResolver {

    private static volatile File defaultConfigFile;

    private MagDbItemResolver() {
    }

    public static void setDefaultConfigFile(File configFile) {
        defaultConfigFile = configFile;
    }

    public static File resolveConfigFile(File explicit) {
        if (explicit != null && explicit.isFile()) {
            return explicit;
        }
        if (defaultConfigFile != null && defaultConfigFile.isFile()) {
            return defaultConfigFile;
        }
        String env = System.getenv("METAPILOT_MAGDB");
        if (env != null && !env.isBlank()) {
            Path envPath = Path.of(env.trim());
            if (Files.isRegularFile(envPath)) {
                return envPath.toFile();
            }
        }
        Path installCfg = MetaPilotInstallPaths.magRegistryPath();
        if (Files.isRegularFile(installCfg)) {
            return installCfg.toFile();
        }
        File cwd = new File(MagDbConfigIO.versionFile);
        if (cwd.isFile()) {
            return cwd;
        }
        File configDir = new File("config", MagDbConfigIO.versionFile);
        if (configDir.isFile()) {
            return configDir;
        }
        return explicit != null ? explicit : installCfg.toFile();
    }

    public static MagDbItem resolve(String catalogueId, String identifier, String catalogueVersion, File configFile) {
        if (catalogueId == null || catalogueId.isBlank()) {
            return null;
        }
        File cfg = resolveConfigFile(configFile);
        if (cfg == null || !cfg.isFile()) {
            return null;
        }
        MagDbItem[][] parsed = MagDbConfigIO.parse(cfg);
        if (parsed == null) {
            return null;
        }
        String id = identifier == null || identifier.isBlank() ? "MGYG" : identifier.trim();
        for (MagDbItem[] group : parsed) {
            if (group == null) {
                continue;
            }
            for (MagDbItem item : group) {
                if (item == null) {
                    continue;
                }
                if (!catalogueId.equals(item.getCatalogueID()) || !id.equals(item.getIdentifier())) {
                    continue;
                }
                if (catalogueVersion != null && !catalogueVersion.isBlank()) {
                    for (MagDbItem.MagVersion version : item.getVersionList()) {
                        if (catalogueVersion.equals(version.getCatalogueVersion())) {
                            item.setUsedVersion(version);
                            break;
                        }
                    }
                }
                return item;
            }
        }
        return null;
    }
}

package uk.ac.quadram.metapilot.core.resources;

import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Install-relative paths for a MetaPilot deployment (for example {@code ~/tools/metapilot/config}).
 */
public final class MetaPilotInstallPaths {

    private MetaPilotInstallPaths() {
    }

    public static Path installHome() {
        String env = System.getenv("METAPILOT_HOME");
        if (env != null && !env.isBlank()) {
            return Paths.get(env.trim()).toAbsolutePath().normalize();
        }
        String property = System.getProperty("metapilot.install.dir");
        if (property != null && !property.isBlank()) {
            return Paths.get(property.trim()).toAbsolutePath().normalize();
        }
        return Paths.get("").toAbsolutePath().normalize();
    }

    public static Path siteConfigDir() {
        return installHome().resolve("config");
    }

    public static Path resourcesConfigPath() {
        return siteConfigDir().resolve("meta_resources.params");
    }

    public static Path magRegistryPath() {
        return siteConfigDir().resolve(MagDbConfigIO.versionFile);
    }

    public static Path magDownloadDir() {
        return installHome().resolve("magdb");
    }
}

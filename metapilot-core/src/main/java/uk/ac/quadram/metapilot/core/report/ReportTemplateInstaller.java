package uk.ac.quadram.metapilot.core.report;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Installs the static HTML report template from classpath resources.
 */
public final class ReportTemplateInstaller {

    private static final String REPORT_RESOURCE_PATH = "/report";

    private ReportTemplateInstaller() {
    }

    public static final class Result {
        public File resultDirectory;
        public File reportDirectory;
        public File htmlReportDirectory;
        public boolean success;
        public final List<String> warnings = new ArrayList<>();
    }

    public static Result install(File resultDirectory, ProgressCallback cb) throws IOException {
        if (cb == null) {
            cb = new ProgressCallback.NoOp();
        }
        if (resultDirectory == null) {
            throw new IOException("resultDirectory is null");
        }

        Result result = new Result();
        result.resultDirectory = resultDirectory;
        result.reportDirectory = new File(resultDirectory, "report");
        result.htmlReportDirectory = new File(result.reportDirectory, "reports");

        cb.onMessage("Installing report template: " + result.reportDirectory.getAbsolutePath());
        URL resourceUrl = ReportTemplateInstaller.class.getResource(REPORT_RESOURCE_PATH);
        if (resourceUrl == null) {
            throw new IOException("Report resource not found at: " + REPORT_RESOURCE_PATH);
        }

        if (!result.reportDirectory.exists() && !result.reportDirectory.mkdirs()) {
            throw new IOException("Could not create report directory: " + result.reportDirectory.getAbsolutePath());
        }

        String protocol = resourceUrl.getProtocol();
        if ("jar".equals(protocol)) {
            copyFromJar(resourceUrl, result.reportDirectory);
        } else if ("file".equals(protocol)) {
            copyFromFileSystem(resourceUrl, result.reportDirectory);
        } else {
            throw new IOException("Unsupported protocol for report resource: " + protocol);
        }

        cleanupUnusedDirectories(result.reportDirectory, result.warnings);
        if (!result.htmlReportDirectory.exists() && !result.htmlReportDirectory.mkdirs()) {
            throw new IOException("Could not create HTML report directory: "
                    + result.htmlReportDirectory.getAbsolutePath());
        }

        result.success = true;
        return result;
    }

    private static void copyFromJar(URL resourceUrl, File reportDirectory) throws IOException {
        String jarPath = resourceUrl.getPath();
        int exclamationIdx = jarPath.indexOf("!");
        if (exclamationIdx > 0) {
            jarPath = jarPath.substring(0, exclamationIdx);
            if (jarPath.startsWith("file:")) {
                jarPath = jarPath.substring(5);
            }
        }

        try (JarFile jarFile = new JarFile(jarPath)) {
            Enumeration<JarEntry> entries = jarFile.entries();
            while (entries.hasMoreElements()) {
                JarEntry entry = entries.nextElement();
                String entryName = entry.getName();
                if (!entryName.startsWith("report/")) {
                    continue;
                }
                String relativePath = entryName.substring("report/".length());
                if (relativePath.isEmpty()) {
                    continue;
                }

                File targetFile = new File(reportDirectory, relativePath);
                if (entry.isDirectory()) {
                    targetFile.mkdirs();
                    continue;
                }
                File parent = targetFile.getParentFile();
                if (parent != null) {
                    parent.mkdirs();
                }
                try (InputStream is = jarFile.getInputStream(entry);
                     FileOutputStream fos = new FileOutputStream(targetFile)) {
                    IOUtils.copy(is, fos);
                }
            }
        }
    }

    private static void copyFromFileSystem(URL resourceUrl, File reportDirectory)
            throws IOException {
        try {
            File sourceDir = new File(resourceUrl.toURI());
            if (!sourceDir.exists()) {
                throw new IOException("Report template source not found: " + sourceDir.getAbsolutePath());
            }
            FileUtils.copyDirectory(sourceDir, reportDirectory);
        } catch (URISyntaxException e) {
            throw new IOException("Invalid report resource URL", e);
        }
    }

    private static void cleanupUnusedDirectories(File reportDirectory, List<String> warnings) {
        deleteQuietly(new File(reportDirectory, "reports/metamap"), warnings);
        deleteQuietly(new File(reportDirectory, "reports/shinyapps"), warnings);
    }

    private static void deleteQuietly(File directory, List<String> warnings) {
        if (!directory.exists()) {
            return;
        }
        try {
            FileUtils.deleteDirectory(directory);
        } catch (IOException e) {
            warnings.add("Could not delete unused report directory " + directory.getAbsolutePath()
                    + ": " + e.getMessage());
        }
    }
}

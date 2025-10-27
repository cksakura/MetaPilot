package uk.ac.quadram.metapilot.task;

import javafx.concurrent.Task;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Python-based report template copy task.
 * 
 * This task is responsible for copying report templates and setting up the directory structure
 * for report generation. It does NOT generate reports - that should be handled separately
 * by MetaPythonReportTask.
 * 
 * This implementation properly loads resources from the classpath, working correctly both
 * in development (IDE) and production (JAR/EXE).
 */
public class MetaPythonReportCopyTask extends Task<Boolean> {

    private static final String TASK_NAME = "Python report copy task";
    protected Logger logger;
    
    // Classpath resource path (works from inside JAR)
    private static final String REPORT_RESOURCE_PATH = "/report";
    
    private final File resultDirectory;
    private final File reportDirectory;
    private final File htmlReportDirectory;

    public MetaPythonReportCopyTask(File resultDirectory) {
        this.resultDirectory = resultDirectory;
        this.reportDirectory = new File(resultDirectory, "report");
        this.htmlReportDirectory = new File(reportDirectory, "reports");
        
        this.logger = LogManager.getLogger(TASK_NAME);
        logger.info("{}: Initialized for {}", TASK_NAME, resultDirectory.getAbsolutePath());
    }


    @Override
    protected Boolean call() throws Exception {
        try {
            updateMessage("Copying report template...");
            updateProgress(10, 100);

            if (!copyReportTemplate()) {
                return false;
            }

            updateMessage("Cleaning up directories...");
            updateProgress(50, 100);
            cleanupUnusedDirectories();

            updateMessage("Creating report directories...");
            updateProgress(80, 100);
            if (!htmlReportDirectory.exists()) {
                htmlReportDirectory.mkdirs();
            }

            updateMessage("Template setup completed successfully");
            updateProgress(100, 100);

            logger.info("{}: Successfully completed template setup", TASK_NAME);
            return true;

        } catch (Exception e) {
        	logger.error("{}: Error during template setup", TASK_NAME, e);
            updateMessage("Error: " + e.getMessage());
            return false;
        }
    }

    /**
     * Copy report template from classpath resources to filesystem.
     * Works both from IDE (directory) and from JAR (packaged resources).
     */
    private boolean copyReportTemplate() {
        try {
            logger.info("{}: Extracting report template from classpath to {}", TASK_NAME, reportDirectory.getAbsolutePath());
            
            // Get the resource URL
            URL resourceUrl = getClass().getResource(REPORT_RESOURCE_PATH);
            
            if (resourceUrl == null) {
                logger.error("{}: Report resource not found at: {}", TASK_NAME, REPORT_RESOURCE_PATH);
                return false;
            }
            
            // Create report directory if it doesn't exist
            if (!reportDirectory.exists()) {
                reportDirectory.mkdirs();
            }
            
            // Check if we're running from a JAR
            String protocol = resourceUrl.getProtocol();
            
            if ("jar".equals(protocol)) {
                // Running from JAR - extract from JAR
                return copyFromJar(resourceUrl);
            } else if ("file".equals(protocol)) {
                // Running from IDE/file system - direct copy
                return copyFromFileSystem(resourceUrl);
            } else {
                logger.error("{}: Unsupported protocol: {}", TASK_NAME, protocol);
                return false;
            }

        } catch (Exception e) {
            logger.error("{}: Failed to copy template", TASK_NAME, e);
            return false;
        }
    }
    
    /**
     * Copy resources from JAR file
     */
    private boolean copyFromJar(URL resourceUrl) throws IOException {
        // Extract JAR path from jar:file:/path/to/app.jar!/report
        String jarPath = resourceUrl.getPath();
        int exclamationIdx = jarPath.indexOf("!");
        if (exclamationIdx > 0) {
            jarPath = jarPath.substring(0, exclamationIdx);
            if (jarPath.startsWith("file:")) {
                jarPath = jarPath.substring(5);
            }
        }
        
        logger.info("{}: Extracting from JAR: {}", TASK_NAME, jarPath);
        
        try (JarFile jarFile = new JarFile(jarPath)) {
            Enumeration<JarEntry> entries = jarFile.entries();
            
            while (entries.hasMoreElements()) {
                JarEntry entry = entries.nextElement();
                String entryName = entry.getName();
                
                // Only process entries in the report directory
                if (entryName.startsWith("report/")) {
                    // Get relative path (remove "report/" prefix)
                    String relativePath = entryName.substring("report/".length());
                    
                    if (relativePath.isEmpty()) {
                        continue; // Skip the directory entry itself
                    }
                    
                    File targetFile = new File(reportDirectory, relativePath);
                    
                    if (entry.isDirectory()) {
                        targetFile.mkdirs();
                    } else {
                        // Create parent directories if needed
                        targetFile.getParentFile().mkdirs();
                        
                        // Copy file from JAR
                        try (InputStream is = jarFile.getInputStream(entry);
                             FileOutputStream fos = new FileOutputStream(targetFile)) {
                            IOUtils.copy(is, fos);
                        }
                        
                        logger.debug("{}: Extracted: {}", TASK_NAME, relativePath);
                    }
                }
            }
        }
        
        logger.info("{}: Successfully extracted template from JAR", TASK_NAME);
        return true;
    }
    
    /**
     * Copy resources from file system (development mode)
     */
    private boolean copyFromFileSystem(URL resourceUrl) throws IOException, URISyntaxException {
        File sourceDir = new File(resourceUrl.toURI());
        
        if (!sourceDir.exists()) {
            logger.error("{}: Source template not found: {}", TASK_NAME, sourceDir.getAbsolutePath());
            return false;
        }
        
        logger.info("{}: Copying from filesystem: {}", TASK_NAME, sourceDir.getAbsolutePath());
        FileUtils.copyDirectory(sourceDir, reportDirectory);
        logger.info("{}: Template copied successfully from filesystem", TASK_NAME);
        return true;
    }

    private void cleanupUnusedDirectories() {
        try {
            File metaMapDir = new File(reportDirectory, "reports/metamap");
            File shinyAppsDir = new File(reportDirectory, "reports/shinyapps");

            if (metaMapDir.exists()) {
                FileUtils.deleteDirectory(metaMapDir);
            }
            if (shinyAppsDir.exists()) {
                FileUtils.deleteDirectory(shinyAppsDir);
            }

        } catch (IOException e) {
        	logger.warn("{}: Failed to cleanup directories", TASK_NAME, e);
        }
    }


    public String getTaskName() {
        return TASK_NAME;
    }

    public File getResultDirectory() {
        return resultDirectory;
    }

    public File getReportDirectory() {
        return reportDirectory;
    }

    public File getHtmlReportDirectory() {
        return htmlReportDirectory;
    }
} 
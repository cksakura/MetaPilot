package uk.ac.quadram.metapilot.task.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Utility class for testing and debugging the manifest system.
 * This can be used to manually create manifest files for testing.
 */
public class ManifestTestUtil {
    
    private static final Logger LOGGER = LogManager.getLogger(ManifestTestUtil.class);
    
    /**
     * Creates a test manifest file in the specified reports directory.
     * 
     * @param reportsDirPath Path to the reports directory
     * @return true if successful, false otherwise
     */
    public static boolean createTestManifest(String reportsDirPath) {
        try {
            File reportsDir = new File(reportsDirPath);
            if (!reportsDir.exists()) {
                reportsDir.mkdirs();
                LOGGER.info("Created reports directory: " + reportsDir.getAbsolutePath());
            }
            
            ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
            manifestGen.createTestManifest();
            
            LOGGER.info("Test manifest created successfully in: " + reportsDir.getAbsolutePath());
            return true;
            
        } catch (Exception e) {
            LOGGER.error("Error creating test manifest", e);
            return false;
        }
    }
    
    /**
     * Scans a reports directory and creates a manifest file based on existing HTML files.
     * 
     * @param reportsDirPath Path to the reports directory
     * @return true if successful, false otherwise
     */
    public static boolean createManifestFromFiles(String reportsDirPath) {
        try {
            File reportsDir = new File(reportsDirPath);
            if (!reportsDir.exists()) {
                LOGGER.warn("Reports directory does not exist: " + reportsDirPath);
                return false;
            }
            
            ReportManifestGenerator manifestGen = new ReportManifestGenerator(reportsDir);
            
            // Scan for HTML files
            File[] files = reportsDir.listFiles((dir, name) -> name.endsWith(".html"));
            if (files != null && files.length > 0) {
                for (File file : files) {
                    if (!file.getName().equals("index.html")) {
                        String title = file.getName().replaceAll("[_-]", " ");
                        title = title.substring(0, title.lastIndexOf('.'));
                        title = title.substring(0, 1).toUpperCase() + title.substring(1);
                        
                        manifestGen.addReportEntry(
                            title,
                            file.getName(),
                            "Generated report file",
                            new String[]{"generated", "html"}
                        );
                    }
                }
                LOGGER.info("Created manifest from " + files.length + " HTML files");
            } else {
                manifestGen.createInitialManifest();
                LOGGER.info("No HTML files found, created empty manifest");
            }
            
            return true;
            
        } catch (Exception e) {
            LOGGER.error("Error creating manifest from files", e);
            return false;
        }
    }
    
    /**
     * Main method for command-line testing.
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Usage: java ManifestTestUtil <reports_directory_path> [test|scan]");
            System.out.println("  test - Create a test manifest with sample data");
            System.out.println("  scan - Create manifest from existing HTML files");
            return;
        }
        
        String reportsDirPath = args[0];
        String mode = args.length > 1 ? args[1] : "test";
        
        boolean success = false;
        switch (mode.toLowerCase()) {
            case "test":
                success = createTestManifest(reportsDirPath);
                break;
            case "scan":
                success = createManifestFromFiles(reportsDirPath);
                break;
            default:
                System.out.println("Unknown mode: " + mode);
                return;
        }
        
        if (success) {
            System.out.println("Manifest creation completed successfully!");
        } else {
            System.out.println("Manifest creation failed!");
        }
    }
}

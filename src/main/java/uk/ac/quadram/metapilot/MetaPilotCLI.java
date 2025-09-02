package uk.ac.quadram.metapilot;

import java.io.File;
import java.util.Date;
import java.text.SimpleDateFormat;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;

/**
 * MetaPilot CLI - Command Line Interface for HPC environments
 * 
 * This class provides a command-line interface for running MetaPilot workflows
 * without requiring JavaFX runtime components.
 */
public class MetaPilotCLI {
    
    private static final SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
    
    public static void main(String[] args) {
        System.out.println("MetaPilot CLI Mode - DDA IGC MaxQuant Workflow");
        System.out.println("==============================================");
        System.out.println("Timestamp: " + format.format(new Date()));
        System.out.println("Java Version: " + System.getProperty("java.version"));
        System.out.println("Working Directory: " + System.getProperty("user.dir"));
        System.out.println();
        
        try {
            // Parse command line arguments
            System.out.println("Parsing command line arguments...");
            String paramFile = "config/last_dda_igc_maxquant_parameters.json";
            String resourceFile = "config/meta_resources.params";
            
            for (int i = 0; i < args.length; i++) {
                System.out.println("  Processing argument " + i + ": " + args[i]);
                if (args[i].equals("--params") && i + 1 < args.length) {
                    paramFile = args[i + 1];
                    System.out.println("    Set parameter file to: " + paramFile);
                    i++; // Skip next argument
                } else if (args[i].equals("--resources") && i + 1 < args.length) {
                    resourceFile = args[i + 1];
                    System.out.println("    Set resource file to: " + resourceFile);
                    i++; // Skip next argument
                } else if (args[i].equals("--help")) {
                    printCLIHelp();
                    return;
                }
            }
            
            System.out.println("Final configuration:");
            System.out.println("  Parameter file: " + paramFile);
            System.out.println("  Resource file: " + resourceFile);
            System.out.println();
            
            // Validate files exist
            File paramFileObj = new File(paramFile);
            File resourceFileObj = new File(resourceFile);
            
            System.out.println("Validating input files...");
            if (!paramFileObj.exists()) {
                System.err.println("ERROR: Parameter file does not exist: " + paramFileObj.getAbsolutePath());
                System.exit(1);
            }
            System.out.println("  ✓ Parameter file exists: " + paramFileObj.getAbsolutePath());
            
            if (!resourceFileObj.exists()) {
                System.err.println("ERROR: Resource file does not exist: " + resourceFileObj.getAbsolutePath());
                System.exit(1);
            }
            System.out.println("  ✓ Resource file exists: " + resourceFileObj.getAbsolutePath());
            System.out.println();
            
            // Load parameters and resources
            System.out.println("Loading resource parameters...");
            MetaResourceIO resourceIO = new MetaResourceIO();
            resourceIO.setConfigFile(resourceFileObj);
            MetaResourceParameters resourceParams = resourceIO.loadParameters();
            System.out.println("  ✓ Loaded resource parameters from: " + resourceFile);
            
            // Validate critical tool paths
            System.out.println("Validating tool paths...");
            String maxquantPath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.MAXQUANT);
            if (maxquantPath == null || maxquantPath.trim().isEmpty()) {
                System.err.println("ERROR: MaxQuant path is not configured in resource file");
                System.exit(1);
            }
            File maxquantFile = new File(maxquantPath);
            if (!maxquantFile.exists()) {
                System.err.println("ERROR: MaxQuant executable not found at: " + maxquantPath);
                System.exit(1);
            }
            System.out.println("  ✓ MaxQuant found at: " + maxquantPath);
            
            String pfindPath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.PFIND);
            if (pfindPath != null && !pfindPath.trim().isEmpty()) {
                File pfindFile = new File(pfindPath);
                if (pfindFile.exists()) {
                    System.out.println("  ✓ pFind found at: " + pfindPath);
                } else {
                    System.out.println("  ⚠ pFind not found at: " + pfindPath);
                }
            }
            System.out.println();
            
            // Load DDA IGC MaxQuant parameters
            System.out.println("Loading workflow parameters...");
            uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameters parameters = 
                uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameterIO.loadLastUsedParameters();
            System.out.println("  ✓ Loaded workflow parameters from: " + paramFile);
            
            // Validate parameters
            System.out.println("Validating workflow parameters...");
            if (parameters.getMetadata() == null) {
                System.err.println("ERROR: No metadata found in parameters");
                System.exit(1);
            }
            String[] rawFiles = parameters.getMetadata().getRawFiles();
            if (rawFiles == null || rawFiles.length == 0) {
                System.err.println("ERROR: No raw files specified in parameters");
                System.exit(1);
            }
            System.out.println("  ✓ Found " + rawFiles.length + " raw files:");
            for (int i = 0; i < rawFiles.length; i++) {
                File rawFile = new File(rawFiles[i]);
                if (rawFile.exists()) {
                    System.out.println("    ✓ " + (i+1) + ": " + rawFiles[i]);
                } else {
                    System.out.println("    ✗ " + (i+1) + ": " + rawFiles[i] + " (NOT FOUND)");
                }
            }
            
            String resultDir = parameters.getResult();
            if (resultDir == null || resultDir.trim().isEmpty()) {
                System.err.println("ERROR: No result directory specified in parameters");
                System.exit(1);
            }
            System.out.println("  ✓ Result directory: " + resultDir);
            System.out.println();
            
            // Create and execute the task
            System.out.println("Creating DDA IGC MaxQuant task...");
            uk.ac.quadram.metapilot.task.dda.igc.DDAIgcMaxQuantTask task = 
                new uk.ac.quadram.metapilot.task.dda.igc.DDAIgcMaxQuantTask(parameters, resourceParams, null, null);
            System.out.println("  ✓ Task created successfully");
            System.out.println();
            
            System.out.println("Starting DDA IGC MaxQuant workflow...");
            System.out.println("=====================================");
            long startTime = System.currentTimeMillis();
            
            // Execute the task
            task.execute();
            
            // Wait for completion and get result
            boolean success = task.get();
            
            long endTime = System.currentTimeMillis();
            long duration = endTime - startTime;
            
            System.out.println("=====================================");
            if (success) {
                System.out.println("✓ Workflow completed successfully!");
                System.out.println("  Duration: " + (duration / 1000.0) + " seconds");
                System.out.println("  Results saved to: " + parameters.getResult());
            } else {
                System.out.println("✗ Workflow failed!");
                System.out.println("  Duration: " + (duration / 1000.0) + " seconds");
                System.exit(1);
            }
            
        } catch (Exception e) {
            System.err.println("ERROR: Exception occurred during CLI workflow execution:");
            System.err.println("  Message: " + e.getMessage());
            System.err.println("  Type: " + e.getClass().getSimpleName());
            System.err.println("  Stack trace:");
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private static void printCLIHelp() {
        System.out.println("MetaPilot CLI Usage:");
        System.out.println("  java -jar metapilot.jar --cli [options]");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  --params <file>     Parameter file (default: config/last_dda_igc_maxquant_parameters.json)");
        System.out.println("  --resources <file>  Resource file (default: config/meta_resources.params)");
        System.out.println("  --help              Show this help message");
        System.out.println();
        System.out.println("Example:");
        System.out.println("  java -jar metapilot.jar --cli --params my_params.json --resources my_resources.params");
    }
}

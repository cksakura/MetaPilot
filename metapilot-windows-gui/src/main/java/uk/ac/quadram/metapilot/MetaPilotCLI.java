package uk.ac.quadram.metapilot;

import java.io.File;
import java.util.Date;
import java.text.SimpleDateFormat;

import javafx.application.Platform;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.ParameterIORegistry;
import uk.ac.quadram.metapilot.workflow.WorkflowType;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipeWorkflow;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;

// DDA IGC parameter classes
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcFragPipeParameters;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcSageParameters;

// DDA MAG parameter classes
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameters;

// DIA parameter classes
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaTg;

// DDA IGC tasks
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcMaxQuantTask;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcPFindTask;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcFragPipeTask;
import uk.ac.quadram.metapilot.task.dda.igc.DDAIgcSageTask;

// DDA MAG tasks
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagPfindTask;
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagFragPipeTask;
import uk.ac.quadram.metapilot.task.dda.mag.task.DDAMagSageTask;

// DIA tasks
import uk.ac.quadram.metapilot.task.dia.MetaDiaMarkerTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaDeepMarkerTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaSeedDirectTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaSeedIterativeTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaTargetTask;
import uk.ac.quadram.metapilot.task.dia.MetaDiaDeepComprehensiveTask;

/**
 * MetaPilot CLI - Command Line Interface for HPC environments
 * 
 * This class provides a command-line interface for running MetaPilot workflows
 * without requiring JavaFX runtime components.
 */
public class MetaPilotCLI {
    
    private static final SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
    private static volatile boolean fxToolkitInitialized = false;
    
    public static void main(String[] args) {
        System.out.println("MetaPilot CLI Mode - Generic Workflow Runner");
        System.out.println("============================================");
        System.out.println("Timestamp: " + format.format(new Date()));
        System.out.println("Java Version: " + System.getProperty("java.version"));
        System.out.println("Working Directory: " + System.getProperty("user.dir"));
        System.out.println();
        
        try {
            // Ensure JavaFX toolkit is initialized so Task/Platform.runLater can be used safely
            ensureFxToolkit();
            // Parse command line arguments
            System.out.println("Parsing command line arguments...");
            String paramFile = "config/last_dda_igc_maxquant_parameters.json";
            String resourceFile = "config/meta_resources.params";
            String fragpipeWorkflowPath = null;
            String magDbConfigPath = null;
            
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
                } else if (args[i].equals("--fragpipe-workflow") && i + 1 < args.length) {
                    fragpipeWorkflowPath = args[i + 1];
                    System.out.println("    Set FragPipe workflow file to: " + fragpipeWorkflowPath);
                    i++; // Skip next argument
                } else if ((args[i].equals("--magdb-config") || args[i].equals("--magdb")) && i + 1 < args.length) {
                    magDbConfigPath = args[i + 1];
                    System.out.println("    Set MAG database config file to: " + magDbConfigPath);
                    i++; // Skip next argument
                } else if (args[i].equals("--help")) {
                    printCLIHelp();
                    return;
                }
            }
            
            System.out.println("Final configuration:");
            System.out.println("  Parameter file: " + paramFile);
            System.out.println("  Resource file: " + resourceFile);
            if (fragpipeWorkflowPath != null) {
                System.out.println("  FragPipe workflow file: " + fragpipeWorkflowPath);
            }
            if (magDbConfigPath != null) {
            System.out.println("  MAG DB config file: " + magDbConfigPath);
            }
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
            
            // Validate critical tool paths for common tools where possible
            System.out.println("Validating tool paths...");
            String maxquantPath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.MAXQUANT);
            if (maxquantPath != null && !maxquantPath.trim().isEmpty()) {
                File maxquantFile = new File(maxquantPath);
                if (maxquantFile.exists()) {
                    System.out.println("  ✓ MaxQuant found at: " + maxquantPath);
                } else {
                    System.out.println("  ⚠ MaxQuant executable not found at: " + maxquantPath);
                }
            }
            
            String pfindPath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.PFIND);
            if (pfindPath != null && !pfindPath.trim().isEmpty()) {
                File pfindFile = new File(pfindPath);
                if (pfindFile.exists()) {
                    System.out.println("  ✓ pFind found at: " + pfindPath);
                } else {
                    System.out.println("  ⚠ pFind not found at: " + pfindPath);
                }
            }
            
            String fragpipePath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.FRAGPIPE);
            if (fragpipePath != null && !fragpipePath.trim().isEmpty()) {
                File fragpipeFile = new File(fragpipePath);
                if (fragpipeFile.exists()) {
                    System.out.println("  ✓ FragPipe found at: " + fragpipePath);
                } else {
                    System.out.println("  ⚠ FragPipe not found at: " + fragpipePath);
                }
            }
            
            String sagePath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.SAGE);
            if (sagePath != null && !sagePath.trim().isEmpty()) {
                File sageFile = new File(sagePath);
                if (sageFile.exists()) {
                    System.out.println("  ✓ Sage found at: " + sagePath);
                } else {
                    System.out.println("  ⚠ Sage not found at: " + sagePath);
                }
            }
            
            String diannPath = resourceParams.getToolPath(uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.DIANN);
            if (diannPath != null && !diannPath.trim().isEmpty()) {
                File diannFile = new File(diannPath);
                if (diannFile.exists()) {
                    System.out.println("  ✓ DIA-NN found at: " + diannPath);
                } else {
                    System.out.println("  ⚠ DIA-NN not found at: " + diannPath);
                }
            }
            System.out.println();
            
            // Load parameters from the specified parameter file
            System.out.println("Loading workflow parameters...");
            // Ensure all ParameterIO handlers are registered so BaseParameterIO can resolve the workflow type
            ParameterIORegistry.registerAllParameterIOs();
            
            Object parameterObject = BaseParameterIO.loadParametersFromFilePath(paramFile);
            if (parameterObject == null) {
                System.err.println("ERROR: Failed to load parameters from file: " + paramFileObj.getAbsolutePath());
                System.exit(1);
            }
            
            String workflowTypeStr = BaseParameterIO.getWorkflowTypeFromParameters(parameterObject);
            System.out.println("Detected workflow type string: " + workflowTypeStr);
            WorkflowType workflowType = WorkflowType.fromAnalysisType(workflowTypeStr);
            if (workflowType == null) {
                System.err.println("ERROR: Unsupported or unknown workflow type: " + workflowTypeStr);
                System.exit(1);
            }
            System.out.println("Resolved workflow: " + workflowType.name());
            
            // Create the appropriate task based on workflow and parameter type
            Object task = null;
            Object typedParams = null;
            String workflowLabel = workflowType.name();

            switch (workflowType) {
                case DDA_IGC_MAXQUANT:
                    if (!(parameterObject instanceof DDAIgcMaxQuantParameters)) {
                        System.err.println("ERROR: Parameters are not DDAIgcMaxQuantParameters for DDA_IGC_MAXQUANT workflow.");
                        System.exit(1);
                    }
                    typedParams = (DDAIgcMaxQuantParameters) parameterObject;
                    task = new DDAIgcMaxQuantTask((DDAIgcMaxQuantParameters) typedParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    break;

                case DDA_IGC_PFIND:
                    if (!(parameterObject instanceof DDAIgcPFindParameters)) {
                        System.err.println("ERROR: Parameters are not DDAIgcPFindParameters for DDA_IGC_PFIND workflow.");
                        System.exit(1);
                    }
                    typedParams = (DDAIgcPFindParameters) parameterObject;
                    task = new DDAIgcPFindTask((DDAIgcPFindParameters) typedParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    break;

                case DDA_IGC_FRAGPIPE:
                    if (!(parameterObject instanceof DDAIgcFragPipeParameters)) {
                        System.err.println("ERROR: Parameters are not DDAIgcFragPipeParameters for DDA_IGC_FRAGPIPE workflow.");
                        System.exit(1);
                    }
                    if (fragpipeWorkflowPath == null || fragpipeWorkflowPath.trim().isEmpty()) {
                        System.err.println("ERROR: --fragpipe-workflow <file> is required for DDA_IGC_FRAGPIPE workflow.");
                        System.exit(1);
                    }
                    File igcWorkflowFile = new File(fragpipeWorkflowPath);
                    if (!igcWorkflowFile.exists()) {
                        System.err.println("ERROR: FragPipe workflow file does not exist: " + igcWorkflowFile.getAbsolutePath());
                        System.exit(1);
                    }
                    FragpipeWorkflow igcWorkflow = FragpipeWorkflow.parse(igcWorkflowFile);
                    if (igcWorkflow == null) {
                        System.err.println("ERROR: Failed to parse FragPipe workflow file: " + igcWorkflowFile.getAbsolutePath());
                        System.exit(1);
                    }
                    typedParams = (DDAIgcFragPipeParameters) parameterObject;
                    task = new DDAIgcFragPipeTask((DDAIgcFragPipeParameters) typedParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar(), igcWorkflow);
                    break;

                case DDA_IGC_SAGE:
                    if (!(parameterObject instanceof DDAIgcSageParameters)) {
                        System.err.println("ERROR: Parameters are not DDAIgcSageParameters for DDA_IGC_SAGE workflow.");
                        System.exit(1);
                    }
                    typedParams = (DDAIgcSageParameters) parameterObject;
                    task = new DDAIgcSageTask((DDAIgcSageParameters) typedParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    break;

                case DDA_MAG_PFIND:
                    if (!(parameterObject instanceof DDAMagPfindParameters)) {
                        System.err.println("ERROR: Parameters are not DDAMagPfindParameters for DDA_MAG_PFIND workflow.");
                        System.exit(1);
                    }
                    typedParams = (DDAMagPfindParameters) parameterObject;
                    // Resolve MAG database from config using catalogue ID and version
                    DDAMagPfindParameters pfindParams = (DDAMagPfindParameters) typedParams;
                    MagDbItem magDbForPfind = resolveMagDbFromConfig(
                            pfindParams.getMagDbCatalogueId(),
                            pfindParams.getMagDbVersion(),
                            magDbConfigPath);
                    pfindParams.setMagDbItem(magDbForPfind);
                    task = new DDAMagPfindTask(pfindParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar(), null);
                    break;

                case DDA_MAG_FRAGPIPE:
                    if (!(parameterObject instanceof DDAMagFragPipeParameters)) {
                        System.err.println("ERROR: Parameters are not DDAMagFragPipeParameters for DDA_MAG_FRAGPIPE workflow.");
                        System.exit(1);
                    }
                    if (fragpipeWorkflowPath == null || fragpipeWorkflowPath.trim().isEmpty()) {
                        System.err.println("ERROR: --fragpipe-workflow <file> is required for DDA_MAG_FRAGPIPE workflow.");
                        System.exit(1);
                    }
                    File magWorkflowFile = new File(fragpipeWorkflowPath);
                    if (!magWorkflowFile.exists()) {
                        System.err.println("ERROR: FragPipe workflow file does not exist: " + magWorkflowFile.getAbsolutePath());
                        System.exit(1);
                    }
                    FragpipeWorkflow magWorkflow = FragpipeWorkflow.parse(magWorkflowFile);
                    if (magWorkflow == null) {
                        System.err.println("ERROR: Failed to parse FragPipe workflow file: " + magWorkflowFile.getAbsolutePath());
                        System.exit(1);
                    }
                    typedParams = (DDAMagFragPipeParameters) parameterObject;
                    DDAMagFragPipeParameters magFragParams = (DDAMagFragPipeParameters) typedParams;
                    // Resolve MAG database from config using catalogue ID and version
                    MagDbItem magDbForFrag = resolveMagDbFromConfig(
                            magFragParams.getMagDbCatalogueId(),
                            magFragParams.getMagDbVersion(),
                            magDbConfigPath);
                    magFragParams.setMagDbItem(magDbForFrag);
                    // statusLabel is GUI-only; pass null in CLI mode
                    task = new DDAMagFragPipeTask(magFragParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar(), null, magWorkflow);
                    break;

                case DDA_MAG_SAGE:
                    if (!(parameterObject instanceof DDAMagSageParameters)) {
                        System.err.println("ERROR: Parameters are not DDAMagSageParameters for DDA_MAG_SAGE workflow.");
                        System.exit(1);
                    }
                    typedParams = (DDAMagSageParameters) parameterObject;
                    DDAMagSageParameters magSageParams = (DDAMagSageParameters) typedParams;
                    // Resolve MAG database from config using catalogue ID and version
                    MagDbItem magDbForSage = resolveMagDbFromConfig(
                            magSageParams.getMagDbCatalogueId(),
                            magSageParams.getMagDbVersion(),
                            magDbConfigPath);
                    magSageParams.setMagDbItem(magDbForSage);
                    task = new DDAMagSageTask(magSageParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar(), null);
                    break;

                case DIA_MARKER:
                    if (!(parameterObject instanceof MetaParameterDiaMarker)) {
                        System.err.println("ERROR: Parameters are not MetaParameterDiaMarker for DIA_MARKER workflow.");
                        System.exit(1);
                    }
                    typedParams = (MetaParameterDiaMarker) parameterObject;
                    MetaParameterDiaMarker diaMarkerParams = (MetaParameterDiaMarker) typedParams;
                    String markerMode = diaMarkerParams.getMarkerMode();
                    if ("DEEP".equalsIgnoreCase(markerMode)) {
                        workflowLabel = "DIA_MARKER (deep)";
                        task = new MetaDiaDeepMarkerTask(diaMarkerParams, resourceParams,
                                new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    } else {
                        workflowLabel = "DIA_MARKER";
                        task = new MetaDiaMarkerTask(diaMarkerParams, resourceParams,
                                new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    }
                    break;

                case DIA_SEED:
                    // Seed workflows can use MetaParameterDia (saved format) and derive the concrete mode from seed_mode
                    if (parameterObject instanceof MetaParameterDia) {
                        MetaParameterDia seedParams = (MetaParameterDia) parameterObject;
                        String seedMode = seedParams.getSeedMode();
                        if (seedMode == null || seedMode.isEmpty()) {
                            // Backward-compatible default: use librarySearch to decide between DIRECT and FOCUS
                            seedMode = seedParams.isLibrarySearch() ? "DIRECT" : "FOCUS";
                        }
                        String upperMode = seedMode.toUpperCase(java.util.Locale.ROOT);

                        if ("DIRECT".equals(upperMode)) {
                            typedParams = seedParams;
                            workflowLabel = "DIA_SEED (direct)";
                            task = new MetaDiaSeedDirectTask(seedParams, resourceParams,
                                    new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                        } else if ("FOCUS".equals(upperMode) || "COMPREHENSIVE".equals(upperMode)) {
                            // Iterative seed workflows (focus/comprehensive) use MetaParameterDiaMarker + MetaDiaSeedIterativeTask
                            MetaParameterDiaMarker focusParams = new MetaParameterDiaMarker(
                                    seedParams, seedParams.getMagDbItem(),
                                    false, false, seedParams.getLibrary(), false, new String[0]);
                            focusParams.setMarkerMode(upperMode);
                            typedParams = focusParams;
                            workflowLabel = "DIA_SEED (" + upperMode.toLowerCase(java.util.Locale.ROOT) + " iterative)";
                            task = new MetaDiaSeedIterativeTask(focusParams, resourceParams,
                                    new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                        } else if ("DEEP".equals(upperMode)) {
                            // Deep comprehensive seed workflow mirrors GUI behaviour
                            MetaParameterDiaMarker deepParams = new MetaParameterDiaMarker(
                                    seedParams, seedParams.getMagDbItem(),
                                    false, false, seedParams.getLibrary(), false, new String[0]);
                            deepParams.setMarkerMode("DEEP");
                            typedParams = deepParams;
                            workflowLabel = "DIA_SEED (deep comprehensive)";
                            task = new MetaDiaDeepComprehensiveTask(deepParams, resourceParams,
                                    new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                        } else {
                            System.err.println("ERROR: Unsupported DIA_SEED mode in parameters: " + seedMode);
                            System.exit(1);
                        }
                    } else if (parameterObject instanceof MetaParameterDiaMarker) {
                        // Fallback: treat marker-style params as iterative focus search
                        MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) parameterObject;
                        typedParams = markerParams;
                        workflowLabel = "DIA_SEED (iterative)";
                        task = new MetaDiaSeedIterativeTask(markerParams, resourceParams,
                                new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    } else {
                        System.err.println("ERROR: Parameters are not MetaParameterDia/MetaParameterDiaMarker for DIA_SEED workflow.");
                        System.exit(1);
                    }
                    break;

                case DIA_TARGET:
                    if (!(parameterObject instanceof MetaParameterDiaTg)) {
                        System.err.println("ERROR: Parameters are not MetaParameterDiaTg for DIA_TARGET workflow.");
                        System.exit(1);
                    }
                    typedParams = (MetaParameterDiaTg) parameterObject;
                    task = new MetaDiaTargetTask((MetaParameterDiaTg) typedParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    break;

                case DIA_DEEP_COMPREHENSIVE:
                    if (!(parameterObject instanceof MetaParameterDiaMarker)) {
                        System.err.println("ERROR: Parameters are not MetaParameterDiaMarker for DIA_DEEP_COMPREHENSIVE workflow.");
                        System.exit(1);
                    }
                    typedParams = (MetaParameterDiaMarker) parameterObject;
                    task = new MetaDiaDeepComprehensiveTask((MetaParameterDiaMarker) typedParams, resourceParams,
                            new javafx.scene.control.ProgressBar(), new javafx.scene.control.ProgressBar());
                    break;

                case DIA_REALTIME:
                    // Real-time DIA monitoring is tightly coupled to the GUI (tables, console panel, etc.)
                    System.err.println("ERROR: DIA_REALTIME workflow is not supported in MetaPilotCLI (GUI-only workflow).");
                    System.exit(1);
                    break;

                default:
                    System.err.println("ERROR: Workflow " + workflowType + " is not supported in MetaPilotCLI.");
                    System.exit(1);
            }

            // At this point we have a MetaAbstractTask<Boolean>-style object
            System.out.println("  ✓ Loaded workflow parameters from: " + paramFileObj.getAbsolutePath());
            System.out.println("Creating task for workflow: " + workflowLabel + " ...");
            System.out.println("  ✓ Task created successfully");
            System.out.println();
            
            System.out.println("Starting workflow: " + workflowLabel);
            System.out.println("============================================");
            long startTime = System.currentTimeMillis();
            
            // Execute the task (MetaAbstractTask.execute()) and wait for completion
            uk.ac.quadram.metapilot.task.MetaAbstractTask metaTask =
                    (uk.ac.quadram.metapilot.task.MetaAbstractTask) task;
            metaTask.execute();
            
            // Wait for completion and get result
            boolean success = metaTask.get();
            
            long endTime = System.currentTimeMillis();
            long duration = endTime - startTime;
            
            System.out.println("============================================");
            if (success) {
                System.out.println("✓ Workflow completed successfully!");
                System.out.println("  Duration: " + (duration / 1000.0) + " seconds");
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
    
    /**
     * Initialize JavaFX toolkit once so that Platform.runLater and Task callbacks work in CLI mode.
     * Safe to call even if the toolkit is already running.
     */
    private static void ensureFxToolkit() {
        if (fxToolkitInitialized) {
            return;
        }
        synchronized (MetaPilotCLI.class) {
            if (fxToolkitInitialized) {
                return;
            }
            try {
                Platform.startup(() -> { /* no-op */ });
            } catch (IllegalStateException e) {
                // Toolkit is already initialized (e.g. running from a GUI context); ignore
            }
            fxToolkitInitialized = true;
        }
    }
    
    /**
     * Resolve a MagDbItem from MAG DB config using catalogue ID and version.
     * Falls back to default MAGDB_1_0.json in the working directory when magDbConfigPath is null.
     */
    private static MagDbItem resolveMagDbFromConfig(String catalogueId, String version, String magDbConfigPath) {
        if (catalogueId == null || catalogueId.trim().isEmpty()) {
            System.err.println("ERROR: MAG database catalogue ID is not set in parameters.");
            System.exit(1);
        }
        if (version == null || version.trim().isEmpty()) {
            System.err.println("ERROR: MAG database version is not set in parameters for catalogue ID: " + catalogueId);
            System.exit(1);
        }

        File configFile;
        if (magDbConfigPath != null && !magDbConfigPath.trim().isEmpty()) {
            configFile = new File(magDbConfigPath);
        } else {
            configFile = new File(MagDbConfigIO.versionFile); // defaults to MAGDB_1_0.json
        }

        if (!configFile.exists()) {
            System.err.println("ERROR: MAG DB config file does not exist: " + configFile.getAbsolutePath());
            System.exit(1);
        }

        MagDbItem[][] magDbItems = MagDbConfigIO.parse(configFile);
        if (magDbItems == null) {
            System.err.println("ERROR: Failed to parse MAG DB config file: " + configFile.getAbsolutePath());
            System.exit(1);
        }

        MagDbItem matched = null;
        MagVersion matchedVersion = null;

        for (MagDbItem[] group : magDbItems) {
            if (group == null) continue;
            for (MagDbItem item : group) {
                if (item == null) continue;
                if (!catalogueId.equals(item.getCatalogueID())) continue;

                if (item.getVersionList() != null) {
                    for (MagVersion v : item.getVersionList()) {
                        if (v != null && version.equals(v.getCatalogueVersion())) {
                            matched = item;
                            matchedVersion = v;
                            break;
                        }
                    }
                }
                if (matched != null) break;
            }
            if (matched != null) break;
        }

        if (matched == null || matchedVersion == null) {
            System.err.println("ERROR: Could not find MAG database for catalogue ID '" + catalogueId +
                    "' and version '" + version + "' in config file: " + configFile.getAbsolutePath());
            System.exit(1);
        }

        matched.setUsedVersion(matchedVersion);
        if (!matchedVersion.isLocalAvailable()) {
            System.err.println("ERROR: Selected MAG database version is not marked as localAvailable: catalogue ID '" +
                    catalogueId + "', version '" + version + "'");
            System.exit(1);
        }

        return matched;
    }
    
    private static void printCLIHelp() {
        System.out.println("MetaPilot CLI Usage:");
        System.out.println("  java -jar metapilot.jar --cli [options]");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  --params <file>            Parameter file (required, any supported workflow)");
        System.out.println("  --resources <file>         Resource file (default: config/meta_resources.params)");
        System.out.println("  --fragpipe-workflow <file> FragPipe .workflow file (required for DDA_IGC_FRAGPIPE and DDA_MAG_FRAGPIPE)");
        System.out.println("  --magdb-config <file>      MAG DB config JSON (default: MAGDB_1_0.json in working directory)");
        System.out.println("  --help                     Show this help message");
        System.out.println();
        System.out.println("Example:");
        System.out.println("  java -jar metapilot.jar --cli --params my_params.json --resources my_resources.params");
        System.out.println("  java -jar metapilot.jar --cli --params my_fragpipe_params.json --resources my_resources.params \\");
        System.out.println("       --fragpipe-workflow fragpipe/Basic-Search.workflow");
    }
}

package uk.ac.quadram.metapilot.task.dda.igc;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.spectra.io.MSConvertor;
import uk.ac.quadram.metapilot.task.igc.MetaPilotIgcTask;
import uk.ac.quadram.metapilot.task.dda.MetaDbReducerTask;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcSageParameters;
import uk.ac.quadram.metapilot.dbSearch.sage.SageTask;
import uk.ac.quadram.metapilot.dbSearch.sage.SageResultReader;
import uk.ac.quadram.metapilot.dbSearch.sage.SagePSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

/**
 * DDA IGC Sage Task - Implements the complete DDA IGC workflow using Sage
 * 
 * This task follows the same workflow pattern as DDAIgcPFindTask but uses Sage
 * as the search engine instead of pFind. The workflow includes:
 * 1. Spectra conversion (using MSConvertor instead of PFindSpectraConvertTask)
 * 2. Database reduction
 * 3. Database search with Sage
 * 4. Quantification
 * 5. Report generation
 * 6. Taxonomy analysis
 * 7. Functional annotation
 */
public class DDAIgcSageTask extends MetaPilotIgcTask {

    // Logger is inherited from MetaPilotIgcTask -> MetaPilotTask as an instance field
    private static final String taskName = "DDA IGC database search by Sage workflow";

    // Task-specific fields
    protected String[] mzMLFiles;
    protected String[] mgfFiles;
    protected File msConvertExe;

    public DDAIgcSageTask(DDAIgcSageParameters metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
            ProgressBar secondBar) {
        super(metaPar, mrp, progressBar, secondBar);
        
        // Debug logging for metadata
        logger.info("{}: Constructor called with parameters", taskName);
        if (metaPar != null) {
            logger.info("{}: Parameters object is not null", taskName);
            if (metaPar.getMetadata() != null) {
                String[] rawFiles = metaPar.getMetadata().getRawFiles();
                String[] expNames = metaPar.getMetadata().getExpNames();
                logger.info("{}: Metadata contains {} raw files", taskName, rawFiles != null ? rawFiles.length : 0);
                if (rawFiles != null) {
                    for (int i = 0; i < rawFiles.length; i++) {
                        logger.info("{}: Raw file {}: {} -> {}", taskName, i, rawFiles[i], 
                            expNames != null && i < expNames.length ? expNames[i] : "N/A");
                    }
                }
            } else {
                logger.warn("{}: Metadata is null in parameters", taskName);
            }
        } else {
            logger.error("{}: Parameters object is null!", taskName);
            return;
        }
        this.spectraDirFile = metaPar.getSpectraFile();
        // Initialize mzML files array and experiment name mapping
        initializeFileArrays();
        
        // Initialize MSConvert executable
        this.msConvertExe = new File(mrp.getToolPath(ToolType.MSCONVERT));
    }

    // Fields to track running subtasks
    private MetaDbReducerTask dbReducerTask;
    private SageTask sageTask;
    private Thread spectraConvertThread;
    private Thread mgfConvertThread;

    @Override
    public void forceStop() {
        logger.info("{}: Force stop requested", taskName);
        System.out.println(format.format(new Date()) + "\t" + taskName + ": Force stop requested");
        
        // Cancel database reduction task if running
        if (dbReducerTask != null && !dbReducerTask.isDone()) {
            logger.info("{}: Cancelling database reduction task", taskName);
            dbReducerTask.cancel(true);
        }
        
        // Cancel Sage task if running
        if (sageTask != null) {
            logger.info("{}: Stopping Sage identification task", taskName);
            // SageTask doesn't have forceStop method, but we can interrupt any running threads
            // The task will be stopped when this main task is cancelled
        }
        
        // Cancel spectra conversion thread if running
        if (spectraConvertThread != null && spectraConvertThread.isAlive()) {
            logger.info("{}: Interrupting spectra conversion thread", taskName);
            spectraConvertThread.interrupt();
        }
        
        // Cancel MGF conversion thread if running
        if (mgfConvertThread != null && mgfConvertThread.isAlive()) {
            logger.info("{}: Interrupting MGF conversion thread", taskName);
            mgfConvertThread.interrupt();
        }
        
        // Cancel this task itself
        this.cancel(true);
        
        logger.info("{}: Force stop completed", taskName);
        System.out.println(format.format(new Date()) + "\t" + taskName + ": Force stop completed");
    }

    @Override
    protected String getTaskName() {
        return taskName;
    }

    @Override
    protected Logger getLogger() {
        return logger;
    }

    @Override
    protected String getVersion() {
        return "1.0.0";
    }

    /**
     * Initialize file arrays and experiment name mapping
     */
    private void initializeFileArrays() {
        DDAIgcSageParameters sageParams = (DDAIgcSageParameters) metaPar;
        String[] rawFiles = sageParams.getMetadata().getRawFiles();
        this.mzMLFiles = new String[rawFiles.length];
        this.mgfFiles = new String[rawFiles.length];
        this.expNameMap = new HashMap<String, String>();
        this.spectraDirFile = sageParams.getSpectraFile();
        
        for (int i = 0; i < rawFiles.length; i++) {
            File rawi = new File(rawFiles[i]);
            this.rawDir = rawi.getParentFile();
            String fileName = rawi.getName();
            String baseName = fileName;
            
            // Handle different file formats
            if (fileName.toLowerCase().endsWith(".mzml")) {
                // File is already in mzML format, use original path directly
                mzMLFiles[i] = rawFiles[i];
            } else if (fileName.toLowerCase().endsWith(".raw")) {
                // File is in .raw format
                baseName = fileName.substring(0, fileName.length() - ".raw".length());
                File mzMLFile = new File(spectraDirFile, baseName + ".mzML");
                mzMLFiles[i] = mzMLFile.getAbsolutePath();
            } else if (fileName.toLowerCase().endsWith(".d")) {
                // File is in .d format
                baseName = fileName.substring(0, fileName.length() - ".d".length());
                File mzMLFile = new File(spectraDirFile, baseName + ".mzML");
                mzMLFiles[i] = mzMLFile.getAbsolutePath();
            } else {
                // Unknown format, assume it needs conversion and use the original name
                File mzMLFile = new File(spectraDirFile, baseName + ".mzML");
                mzMLFiles[i] = mzMLFile.getAbsolutePath();
            }
            
            // Set corresponding MGF file path
            File mgfFile = new File(spectraDirFile, baseName + ".mgf");
            mgfFiles[i] = mgfFile.getAbsolutePath();
        }
    }

    @Override
    protected Boolean call() throws Exception {
        // Update progress at start
        updateProgress(0, 100);
        updateMessage("Starting DDA IGC Sage workflow...");

        updateMessage(getTaskName());

        logger.info(getTaskName() + ": start");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": start");

        // Update progress before spectra conversion
        updateProgress(5, 100);
        updateMessage("Starting spectra conversion...");

        // Step 1: Spectra Conversion using MSConvertor (similar to DDAMagSageTask)
        boolean spectraConverted = convertSpectraToMzML();
        if (!spectraConverted) {
            logger.error(getTaskName() + ": error in spectra conversion workflow");
            System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in spectra conversion workflow");
            return false;
        }

        // Update progress after spectra conversion
        updateProgress(20, 100);
        updateMessage("Spectra conversion completed, starting MGF conversion...");

        // Step 2: Convert mzML files to MGF format
        boolean mgfConverted = convertMzMLToMgf();
        if (!mgfConverted) {
            logger.error(getTaskName() + ": error in MGF conversion workflow");
            System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in MGF conversion workflow");
            return false;
        }

        // Update progress after MGF conversion
        updateProgress(35, 100);
        updateMessage("MGF conversion completed, starting database reduction...");

        // Step 3: Database reduction using MetaDbReducerTask
        dbReducerTask = new MetaDbReducerTask(metaPar, mrp, secondBar, logger);
        logger.info("{}: MetaDbReducerTask created, executing...", taskName);
        dbReducerTask.execute();
        try {
            boolean dbReduced = dbReducerTask.get(); // Wait for completion and get result
            if (!dbReduced) {
                logger.error(getTaskName() + ": error in sample-specific database generation");
                System.err.println(format.format(new Date()) + "\t" + getTaskName()
                        + ": error in sample-specific database generation");
                return false;
            }
        } catch (InterruptedException | ExecutionException e) {
            logger.error(getTaskName() + ": error in sample-specific database generation", e);
            System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in sample-specific database generation: " + e.getMessage());
            return false;
        }

        // Update progress after database reduction
        updateProgress(50, 100);
        updateMessage("Database reduction completed, starting identification...");

        logger.info(getTaskName() + ": sample specific database generation: finished");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": sample specific database generation: finished");

        // Step 4: Database search with Sage
        boolean identification = performIdentification();
        if (!identification) {
            logger.error(getTaskName() + ": error in database search and identification");
            System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in database search and identification");
            return false;
        }

        // Update progress after identification
        updateProgress(70, 100);
        updateMessage("Identification completed, starting quantification...");

        logger.info(getTaskName() + ": database search and identification: finished");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": database search and identification: finished");

        

        // Read Sage results
        File searchResultDir = metaPar.getDbSearchResultFile();
        File pepResultFile = new File(searchResultDir, "pep.tsv");
        File proResultFile = new File(searchResultDir, "pro.tsv");

        if (!pepResultFile.exists() || !proResultFile.exists()) {
            logger.error("Sage result files not found: {} or {}", pepResultFile.getAbsolutePath(), proResultFile.getAbsolutePath());
            return false;
        }
        
		SageResultReader sageReader = new SageResultReader(searchResultDir);
		super.metaPeptides = sageReader.getMetaPeptides();
		if (metaPeptides == null || metaPeptides.length == 0) {
			logger.error("Sage result files not found: {} or {}", pepResultFile.getAbsolutePath(),
					proResultFile.getAbsolutePath());
			return false;
		}
		
		
        
        // Step 5: Report generation
        boolean reportGenerated = exportReport();
        if (!reportGenerated) {
            logger.error(getTaskName() + ": error in report generation");
            System.err.println(format.format(new Date()) + "\t" + getTaskName() + ": error in report generation");
            return false;
        }

        // Update progress to completion
        updateProgress(100, 100);
        updateMessage("Workflow completed successfully!");

        return true;
    }

    /**
     * Convert spectra files to mzML format for Sage processing using multiple threads
     * This method is similar to the one in DDAMagSageTask
     */
    private boolean convertSpectraToMzML() {
        logger.info("Starting spectra conversion to mzML format using MSConvert");
        
        // Validate MSConvert executable
        if (msConvertExe == null || !msConvertExe.exists()) {
            logger.error("MSConvert executable not found: {}", 
                msConvertExe != null ? msConvertExe.getAbsolutePath() : "null");
            return false;
        }
        
        DDAIgcSageParameters sageParams = (DDAIgcSageParameters) metaPar;
        String[] rawFiles = sageParams.getMetadata().getRawFiles();
        if (rawFiles == null || rawFiles.length == 0) {
            logger.error("No raw files found in metadata");
            return false;
        }
        
        // Ensure spectra directory exists
        File spectraDir = sageParams.getSpectraFile();
        if (!spectraDir.exists()) {
            spectraDir.mkdirs();
            logger.info("Created spectra directory: {}", spectraDir.getAbsolutePath());
        }
        
        // Collect files that need conversion
        List<Integer> filesToConvert = new ArrayList<>();
        for (int i = 0; i < rawFiles.length; i++) {
            String rawFile = rawFiles[i];
            File rawFileObj = new File(rawFile);
            
            if (!rawFileObj.exists()) {
                logger.error("Raw file does not exist: {}", rawFile);
                return false;
            }
            
            // Check if the mzML file already exists (mzMLFiles[i] was set in constructor)
            File mzMLFile = new File(mzMLFiles[i]);
            
            if (mzMLFile.exists() && mzMLFile.length() > 0) {
                logger.info("mzML file already exists, skipping conversion: {}", mzMLFile.getAbsolutePath());
            } else {
                filesToConvert.add(i);
            }
        }
        
        // If no files need conversion, return success
        if (filesToConvert.isEmpty()) {
            logger.info("All mzML files already exist, no conversion needed");
            return true;
        }
        
        // Use multiple threads for conversion
        int threadCount = Math.min(filesToConvert.size(), sageParams.getThreadCount());
        ExecutorService executor = Executors.newFixedThreadPool(threadCount);
        
        logger.info("Processing {} files using {} threads", filesToConvert.size(), threadCount);
        
        // Submit conversion/copy tasks
        for (int fileIndex : filesToConvert) {
            executor.submit(() -> {
                String rawFile = rawFiles[fileIndex];
                File rawFileObj = new File(rawFile);
                String fileName = rawFileObj.getName().toLowerCase();
                
                try {
                    String convertedPath = null;
                    
                    if (fileName.endsWith(".raw")) {
                        logger.info("Converting .raw file to mzML: {}", rawFile);
                        convertedPath = MSConvertor.raw2mzML(msConvertExe.getAbsolutePath(), rawFile);
                    } else if (fileName.endsWith(".d")) {
                        logger.info("Converting .d file to mzML: {}", rawFile);
                        convertedPath = MSConvertor.d2mzML(msConvertExe.getAbsolutePath(), rawFile, spectraDir.getAbsolutePath());
                    } else {
                        logger.error("Unsupported file format: {}. Only .raw and .d files are supported for conversion.", fileName);
                        return;
                    }
                    
                    if (convertedPath != null) {
                        // Verify the converted file exists
                        File convertedFile = new File(convertedPath);
                        if (convertedFile.exists() && convertedFile.length() > 0) {
                            // Move the converted file to the spectra directory
                            File targetFile = new File(mzMLFiles[fileIndex]);
                            try {
                                Files.move(convertedFile.toPath(), targetFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                                logger.info("Successfully moved converted file from {} to {}", convertedPath, mzMLFiles[fileIndex]);
                            } catch (Exception e) {
                                logger.error("Failed to move converted file from {} to {}: {}", convertedPath, mzMLFiles[fileIndex], e.getMessage());
                            }
                        } else {
                            logger.error("Converted mzML file does not exist or is empty: {}", convertedPath);
                        }
                    } else {
                        logger.error("Failed to convert file to mzML: {}", rawFile);
                    }
                } catch (Exception e) {
                    logger.error("Error converting file {}: {}", rawFile, e.getMessage());
                }
            });
        }
        
        // Shutdown executor and wait for completion
        executor.shutdown();
        try {
            // Wait for all tasks to complete
            while (!executor.isTerminated()) {
                Thread.sleep(1000);
            }
        } catch (InterruptedException e) {
            logger.error("Interrupted while waiting for conversion tasks to complete", e);
            return false;
        }
        
        // Verify all conversions were successful
        for (int fileIndex : filesToConvert) {
            File mzMLFile = new File(mzMLFiles[fileIndex]);
            if (!mzMLFile.exists() || mzMLFile.length() == 0) {
                logger.error("Conversion failed for file: {}", mzMLFiles[fileIndex]);
                return false;
            }
        }
        
        logger.info("All spectra files converted successfully to mzML format");
        return true;
    }

    /**
     * Convert mzML files to MGF format for Sage processing using multiple threads
     */
    private boolean convertMzMLToMgf() {
        logger.info("Starting mzML to MGF conversion using MSConvert");
        
        // Validate MSConvert executable
        if (msConvertExe == null || !msConvertExe.exists()) {
            logger.error("MSConvert executable not found: {}", 
                msConvertExe != null ? msConvertExe.getAbsolutePath() : "null");
            return false;
        }
        
        DDAIgcSageParameters sageParams = (DDAIgcSageParameters) metaPar;
        String[] rawFiles = sageParams.getMetadata().getRawFiles();
        if (rawFiles == null || rawFiles.length == 0) {
            logger.error("No raw files found in metadata");
            return false;
        }
        
        // Collect files that need MGF conversion
        List<Integer> filesToConvert = new ArrayList<>();
        for (int i = 0; i < mzMLFiles.length; i++) {
            File mzMLFile = new File(mzMLFiles[i]);
            File mgfFile = new File(mgfFiles[i]);
            
            if (!mzMLFile.exists() || mzMLFile.length() == 0) {
                logger.error("mzML file does not exist or is empty: {}", mzMLFile.getAbsolutePath());
                return false;
            }
            
            if (mgfFile.exists() && mgfFile.length() > 0) {
                logger.info("MGF file already exists, skipping conversion: {}", mgfFile.getAbsolutePath());
            } else {
                filesToConvert.add(i);
            }
        }
        
        // If no files need conversion, return success
        if (filesToConvert.isEmpty()) {
            logger.info("All MGF files already exist, no conversion needed");
            return true;
        }
        
        // Use multiple threads for conversion
        int threadCount = Math.min(filesToConvert.size(), sageParams.getThreadCount());
        ExecutorService executor = Executors.newFixedThreadPool(threadCount);
        
        logger.info("Converting {} mzML files to MGF using {} threads", filesToConvert.size(), threadCount);
        
        // Submit MGF conversion tasks
        for (int fileIndex : filesToConvert) {
            executor.submit(() -> {
                String mzMLFile = mzMLFiles[fileIndex];
                String mgfFile = mgfFiles[fileIndex];
                
                try {
                    logger.info("Converting mzML to MGF: {} -> {}", mzMLFile, mgfFile);
                    // MSConvertor.mzml2mgf returns void, so we need to check if the file was created
                    MSConvertor.mzml2mgf(msConvertExe.getAbsolutePath(), mzMLFile, spectraDirFile.getAbsolutePath());
                    
                    // Check if the MGF file was created in the spectra directory
                    File mzMLFileObj = new File(mzMLFile);
                    String baseName = mzMLFileObj.getName();
                    baseName = baseName.substring(0, baseName.lastIndexOf("."));
                    File expectedMgfFile = new File(spectraDirFile, baseName + ".mgf");
                    
                    if (expectedMgfFile.exists() && expectedMgfFile.length() > 0) {
                        // Move the converted file to the target location if it's different
                        File targetFile = new File(mgfFile);
                        if (!expectedMgfFile.getAbsolutePath().equals(targetFile.getAbsolutePath())) {
                            try {
                                Files.move(expectedMgfFile.toPath(), targetFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                                logger.info("Successfully moved MGF file from {} to {}", expectedMgfFile.getAbsolutePath(), mgfFile);
                            } catch (Exception e) {
                                logger.error("Failed to move MGF file from {} to {}: {}", expectedMgfFile.getAbsolutePath(), mgfFile, e.getMessage());
                            }
                        } else {
                            logger.info("MGF file already in correct location: {}", mgfFile);
                        }
                    } else {
                        logger.error("MGF file was not created: {}", expectedMgfFile.getAbsolutePath());
                    }
                } catch (Exception e) {
                    logger.error("Error converting mzML to MGF {}: {}", mzMLFile, e.getMessage());
                }
            });
        }
        
        // Shutdown executor and wait for completion
        executor.shutdown();
        try {
            // Wait for all tasks to complete
            while (!executor.isTerminated()) {
                Thread.sleep(1000);
            }
        } catch (InterruptedException e) {
            logger.error("Interrupted while waiting for MGF conversion tasks to complete", e);
            return false;
        }
        
        // Verify all conversions were successful
        for (int fileIndex : filesToConvert) {
            File mgfFile = new File(mgfFiles[fileIndex]);
            if (!mgfFile.exists() || mgfFile.length() == 0) {
                logger.error("MGF conversion failed for file: {}", mgfFiles[fileIndex]);
                return false;
            }
        }
        
        logger.info("All mzML files converted successfully to MGF format");
        return true;
    }

    /**
     * Perform database search and identification using Sage
     * This method follows the same pattern as DDAMagSageTask
     */
    private boolean performIdentification() {
        logger.info("Starting Sage database search and identification");
        
        // Get Sage executable path
        String sageExe = mrp.getToolPath(ToolType.SAGE);
        if (sageExe == null || sageExe.isEmpty()) {
            logger.error("Sage executable path not found in MetaResourceParameters");
            return false;
        }
        
        // Create Sage task
        sageTask = new SageTask(sageExe);
        
        // Get current database path (should be the reduced database)
        String currentDb = metaPar.getCurrentDb();
        if (currentDb == null || currentDb.isEmpty()) {
            logger.error("Current database path is not set");
            return false;
        }
        
        // Create search result directory
        File sageSearchFile = new File(metaPar.getResult(), "sage_IGC_search");
		if (!sageSearchFile.exists()) {
			sageSearchFile.mkdir();
			logger.info("  Created Sage search directory: " + sageSearchFile.getAbsolutePath());
		} else {
			logger.info("  Sage search directory already exists: " + sageSearchFile.getAbsolutePath());
		}
        
        // Set the database search result file
        metaPar.setDbSearchResultFile(sageSearchFile);
        
        this.combinedFile = new File(sageSearchFile, "combined");
		if (!combinedFile.exists()) {
			combinedFile.mkdir();
			logger.info("  Created combined search directory: " + combinedFile.getAbsolutePath());
		} else {
			logger.info("  Combined search directory already exists: " + combinedFile.getAbsolutePath());
		}
		
        try {
            // Run Sage search against the database
            DDAIgcSageParameters sageParams = (DDAIgcSageParameters) metaPar;

			sageTask.searchMzML(mzMLFiles, combinedFile.getAbsolutePath(), currentDb, sageParams, true);

            // Run Sage
            logger.info("Starting Sage execution with {} threads", metaPar.getThreadCount());
            sageTask.run(metaPar.getThreadCount());
            
            // Check if results were generated
            File pepResultFile = new File(combinedFile, "pep.tsv");
            File proResultFile = new File(combinedFile, "pro.tsv");
            
            if (pepResultFile.exists() && pepResultFile.length() > 0 && 
                proResultFile.exists() && proResultFile.length() > 0) {
                logger.info("Sage identification completed successfully");
                return true;
            } else {
                logger.error("Sage identification failed - result files not found or empty");
                return false;
            }
            
        } catch (Exception e) {
            logger.error("Error during Sage identification: {}", e.getMessage(), e);
            return false;
        }
    }

    /**
     * Export report including peptide/protein reports, taxonomy analysis, and functional annotation
     * This method follows the same pattern as DDAIgcPFindTask but uses Sage readers
     */
    protected boolean exportReport() {
        logger.info(getTaskName() + ": exporting peptide report started");
        System.out.println(format.format(new Date()) + "\t" + getTaskName() + ": exporting peptide report started");

        // Build file list and mapping
        ArrayList<String> fileList = new ArrayList<String>();
        HashMap<String, Integer> fileIdMap = new HashMap<String, Integer>();
        
        // Get experiment names from raw files
        String[] rawFiles = this.metaPar.getMetadata().getRawFiles();
        String[] quanExpNames = new String[rawFiles.length];
        for (int i = 0; i < rawFiles.length; i++) {
            quanExpNames[i] = rawFiles[i].substring(rawFiles[i].lastIndexOf("\\") + 1, rawFiles[i].lastIndexOf("."));
            fileIdMap.put(quanExpNames[i], fileIdMap.size());
            fileList.add(quanExpNames[i]);
        }

        if (!combinedFile.exists()) {
            logger.error(getTaskName() + ": Sage combined result file not found: " + combinedFile.getAbsolutePath());
            return false;
        }

        SageResultReader sageResultReader = new SageResultReader(combinedFile);

        MetaPeptide[] quanPeps = sageResultReader.getMetaPeptides();
        SagePSM[] psms = sageResultReader.getPSMs();
        String[] intensityTitles = quanExpNames; // Use experiment names as intensity titles

        // Create and execute the post-analysis task directly
        DDAIgcPostAnalysisTask postAnalysisTask = new DDAIgcPostAnalysisTask(
            metaPar, mrp, secondBar, quanPeps, psms, fileIdMap, fileList, MetaConstants.labelFree, intensityTitles, getTaskName(), this);
        
        // Execute the task and wait for completion
        Thread taskThread = new Thread(postAnalysisTask);
        taskThread.setDaemon(true);
        taskThread.start();
        
        try {
            taskThread.join();
            boolean success = postAnalysisTask.get();
            
            if (success) {
                logger.info(getTaskName() + ": Post-analysis completed successfully");
                return true;
                } else {
                logger.error(getTaskName() + ": Post-analysis task failed");
            return false;
        }
        } catch (InterruptedException e) {
            logger.error(getTaskName() + ": Post-analysis task was interrupted", e);
            Thread.currentThread().interrupt();
        return false;
        } catch (ExecutionException e) {
            logger.error(getTaskName() + ": Post-analysis task failed with exception", e);
            return false;
        }
    }
}

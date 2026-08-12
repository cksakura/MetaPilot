package uk.ac.quadram.metapilot.web.service;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;
import uk.ac.quadram.metapilot.web.websocket.ProgressBroadcastService;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;

/**
 * Service for executing MetaPilot workflows from the Windows CLI.
 * 
 * <p>This service wraps the existing MetaPilot CLI to execute workflows
 * from the web API. Since MetaPilot tasks currently use JavaFX dependencies,
 * we leverage the existing CLI implementation rather than duplicating logic.</p>
 * 
 * <p>In the future, this should be refactored to use a pure Java
 * implementation with the ProgressCallback interface to avoid JavaFX
 * dependencies in the web module.</p>
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
@Service
public class WorkflowService {
    
    private static final Logger log = LoggerFactory.getLogger(WorkflowService.class);
    
    private final ExecutorService executorService;
    private final ProgressBroadcastService broadcastService;
    private final String configuredJarPath;
    
    // Track running jobs
    private final Map<String, Future<?>> runningJobs = new ConcurrentHashMap<>();
    
    public WorkflowService(ProgressBroadcastService broadcastService,
                           @Value("${metapilot.cli.jar.path:}") String configuredJarPath) {
        // Create a thread pool for workflow execution
        executorService = Executors.newFixedThreadPool(4); // Max 4 concurrent jobs
        this.broadcastService = broadcastService;
        this.configuredJarPath = configuredJarPath;
    }
    
    /**
     * Submit a workflow job for execution.
     * 
     * @param jobId Unique identifier for this job
     * @param paramFile Path to the parameter JSON file
     * @param resourceFile Path to the resource parameters file
     * @return Job ID for tracking
     */
    public String submitJob(String jobId, String paramFile, String resourceFile) {
        log.info("Submitting job {} with paramFile={}, resourceFile={}", jobId, paramFile, resourceFile);
        
        // Validate files exist
        File param = new File(paramFile);
        File resource = new File(resourceFile);
        
        if (!param.exists()) {
            throw new IllegalArgumentException("Parameter file does not exist: " + paramFile);
        }
        if (!resource.exists()) {
            throw new IllegalArgumentException("Resource file does not exist: " + resourceFile);
        }
        
        // Submit job asynchronously
        Future<?> future = executorService.submit(() -> {
            try {
                executeWorkflow(jobId, paramFile, resourceFile);
            } catch (Exception e) {
                log.error("Job {} failed with exception", jobId, e);
            }
        });
        
        runningJobs.put(jobId, future);
        log.info("Job {} submitted successfully", jobId);
        
        return jobId;
    }
    
    /**
     * Check if a job is still running.
     * 
     * @param jobId Job identifier
     * @return true if job is still running
     */
    public boolean isJobRunning(String jobId) {
        Future<?> future = runningJobs.get(jobId);
        return future != null && !future.isDone();
    }
    
    /**
     * Cancel a running job.
     * 
     * @param jobId Job identifier
     * @return true if job was cancelled
     */
    public boolean cancelJob(String jobId) {
        Future<?> future = runningJobs.get(jobId);
        if (future != null && !future.isDone()) {
            boolean cancelled = future.cancel(true);
            if (cancelled) {
                runningJobs.remove(jobId);
                broadcastService.broadcastCancellation(jobId, "Job cancelled by user");
                log.info("Job {} cancelled", jobId);
            }
            return cancelled;
        }
        return false;
    }
    
    /**
     * Execute a workflow using the MetaPilot CLI.
     * This calls the existing CLI implementation to run workflows as a subprocess.
     * 
     * @param jobId Job identifier
     * @param paramFile Parameter file path
     * @param resourceFile Resource file path
     * @throws IOException if execution fails
     */
    private void executeWorkflow(String jobId, String paramFile, String resourceFile) throws IOException {
        log.info("Executing workflow for job {}", jobId);
        
        Process process = null;
        BufferedReader outputReader = null;
        BufferedReader errorReader = null;
        
        try {
            // Broadcast initial progress
            broadcastService.broadcastProgress(jobId, "RUNNING", 0.0, "Starting workflow execution...");
            
            // Find the MetaPilot Windows GUI JAR
            String jarPath = findMetaPilotJar();
            if (jarPath == null) {
                throw new IOException("MetaPilot JAR not found. Please ensure the Windows GUI JAR is available.");
            }
            
            log.info("Using MetaPilot JAR: {}", jarPath);
            broadcastService.broadcastProgress(jobId, "RUNNING", 0.05, "Initializing MetaPilot CLI...");
            
            // Build command to execute
            ProcessBuilder pb = new ProcessBuilder(
                "java", "-jar", jarPath,
                "--params", paramFile,
                "--resources", resourceFile
            );
            
            // Set working directory to the same as the JAR
            File jarFile = new File(jarPath);
            pb.directory(jarFile.getParentFile());
            
            // Redirect error stream to output stream
            pb.redirectErrorStream(true);
            
            log.info("Starting MetaPilot CLI process for job {}", jobId);
            process = pb.start();
            
            // Monitor output for progress updates
            outputReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
            
            String line;
            int progressStage = 0; // Track current stage
            final int TOTAL_STAGES = 10;
            
            // Read output line by line
            while ((line = outputReader.readLine()) != null) {
                if (Thread.currentThread().isInterrupted()) {
                    log.warn("Job {} was interrupted, destroying process", jobId);
                    process.destroyForcibly();
                    broadcastService.broadcastCancellation(jobId, "Job was cancelled");
                    return;
                }
                
                log.info("CLI output: {}", line);
                
                // Try to extract progress information from output
                // Look for percentage patterns or stage indicators
                String progressMessage = extractProgressMessage(line);
                if (progressMessage != null) {
                    progressStage = Math.min(progressStage + 1, TOTAL_STAGES);
                    double progress = (double) progressStage / TOTAL_STAGES;
                    broadcastService.broadcastProgress(jobId, "RUNNING", progress, progressMessage);
                }
                
                // Check for completion or error indicators
                if (line.contains("Workflow completed successfully") || 
                    line.contains("✓ Workflow completed")) {
                    broadcastService.broadcastProgress(jobId, "RUNNING", 0.95, "Finalizing results...");
                } else if (line.contains("ERROR:") || line.contains("failed") || 
                          line.contains("Exception")) {
                    String errorMsg = "Error detected: " + line;
                    log.error(errorMsg);
                    broadcastService.broadcastCompletion(jobId, false, errorMsg);
                }
            }
            
            // Wait for process to complete
            int exitCode = process.waitFor();
            
            if (Thread.currentThread().isInterrupted()) {
                log.warn("Job {} was interrupted during wait", jobId);
                broadcastService.broadcastCancellation(jobId, "Job was cancelled");
                return;
            }
            
            if (exitCode == 0) {
                broadcastService.broadcastCompletion(jobId, true, "Workflow completed successfully");
                log.info("Job {} execution completed successfully", jobId);
            } else {
                String errorMsg = "Workflow execution failed with exit code: " + exitCode;
                broadcastService.broadcastCompletion(jobId, false, errorMsg);
                log.error("Job {} execution failed with exit code {}", jobId, exitCode);
                throw new IOException(errorMsg);
            }
            
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            if (process != null) {
                process.destroyForcibly();
            }
            broadcastService.broadcastCancellation(jobId, "Job was cancelled");
            log.warn("Job {} was interrupted", jobId);
        } catch (Exception e) {
            if (process != null) {
                process.destroyForcibly();
            }
            String errorMsg = "Workflow execution failed: " + e.getMessage();
            broadcastService.broadcastCompletion(jobId, false, errorMsg);
            log.error("Job {} failed", jobId, e);
            throw new IOException(errorMsg, e);
        } finally {
            // Clean up resources
            if (outputReader != null) {
                try { outputReader.close(); } catch (IOException e) { log.warn("Error closing output reader", e); }
            }
            if (errorReader != null) {
                try { errorReader.close(); } catch (IOException e) { log.warn("Error closing error reader", e); }
            }
            if (process != null && process.isAlive()) {
                process.destroyForcibly();
            }
        }
    }
    
    /**
     * Find the MetaPilot Windows GUI JAR file.
     * First checks configured path, then searches common locations.
     * 
     * @return Path to the MetaPilot JAR, or null if not found
     */
    private String findMetaPilotJar() {
        // If configured path is provided, use it
        if (configuredJarPath != null && !configuredJarPath.trim().isEmpty()) {
            File jarFile = new File(configuredJarPath);
            if (jarFile.exists()) {
                log.info("Found MetaPilot JAR at configured path: {}", jarFile.getAbsolutePath());
                return jarFile.getAbsolutePath();
            } else {
                log.warn("Configured JAR path does not exist: {}", configuredJarPath);
            }
        }
        
        // Otherwise, try common locations
        String[] possibleJars = {
            "../metapilot-windows-gui/target/metapilot-windows-gui-1.1.1-shaded.jar",
            "../../metapilot-windows-gui/target/metapilot-windows-gui-1.1.1-shaded.jar",
            "./metapilot-windows-gui-1.1.1-shaded.jar",
            "../metapilot-windows-gui-1.1.1-shaded.jar"
        };
        
        for (String jarPath : possibleJars) {
            File jarFile = new File(jarPath);
            if (jarFile.exists()) {
                log.info("Found MetaPilot JAR at: {}", jarFile.getAbsolutePath());
                return jarFile.getAbsolutePath();
            }
        }
        
        log.error("Could not find MetaPilot JAR in any expected location");
        return null;
    }
    
    /**
     * Extract progress information from CLI output line.
     * 
     * @param line Output line from CLI
     * @return Progress message if progress detected, null otherwise
     */
    private String extractProgressMessage(String line) {
        // Look for common progress indicators in MetaPilot CLI output
        if (line.contains("Loading") || line.contains("Validating")) {
            return line;
        }
        if (line.contains("Creating") || line.contains("Starting")) {
            return line;
        }
        if (line.contains("Processing") || line.contains("Completed")) {
            return line;
        }
        if (line.contains("%")) {
            return line;
        }
        return null;
    }
    
    /**
     * Get job status information.
     * 
     * @param jobId Job identifier
     * @return Job status details
     */
    public Map<String, Object> getJobStatus(String jobId) {
        Map<String, Object> status = new HashMap<>();
        status.put("jobId", jobId);
        status.put("running", isJobRunning(jobId));
        
        Future<?> future = runningJobs.get(jobId);
        if (future != null) {
            status.put("done", future.isDone());
            if (future.isDone()) {
                try {
                    future.get(); // Retrieve any exception
                    status.put("completed", true);
                    runningJobs.remove(jobId);
                } catch (Exception e) {
                    status.put("completed", false);
                    status.put("error", e.getMessage());
                    runningJobs.remove(jobId);
                }
            }
        } else {
            status.put("notFound", true);
        }
        
        return status;
    }
    
    /**
     * Cleanup resources on shutdown.
     */
    public void shutdown() {
        log.info("Shutting down WorkflowService");
        executorService.shutdown();
        try {
            if (!executorService.awaitTermination(30, TimeUnit.SECONDS)) {
                executorService.shutdownNow();
            }
        } catch (InterruptedException e) {
            executorService.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }
}


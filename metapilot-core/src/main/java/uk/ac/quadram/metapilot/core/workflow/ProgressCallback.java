package uk.ac.quadram.metapilot.core.workflow;

/**
 * Generic progress callback interface for workflow execution.
 * 
 * This interface provides a JavaFX-independent way to report progress
 * from workflow tasks. Implementations can be:
 * - JavaFX Property-based (for GUI)
 * - WebSocket-based (for web UI)
 * - Console-based (for CLI)
 * - Log-based (for batch processing)
 * 
 * @author MetaPilot Team
 * @version 1.0
 * @since 1.0.0
 */
public interface ProgressCallback {
    
    /**
     * Report progress on a task.
     * 
     * @param workDone Amount of work completed (e.g., files processed, bytes downloaded)
     * @param max Total amount of work (e.g., total files, total bytes)
     */
    void onProgress(double workDone, double max);
    
    /**
     * Send a message/status update.
     * 
     * @param message Status message or log entry
     */
    void onMessage(String message);
    
    /**
     * Indicate that the task has completed successfully.
     * 
     * @param result The result of the task
     */
    void onComplete(Object result);
    
    /**
     * Indicate that the task has failed with an error.
     * 
     * @param error The exception that occurred
     */
    void onError(Exception error);
    
    /**
     * Default no-op implementation for tasks that don't need progress reporting.
     */
    class NoOp implements ProgressCallback {
        @Override
        public void onProgress(double workDone, double max) {
            // Do nothing
        }
        
        @Override
        public void onMessage(String message) {
            // Do nothing
        }
        
        @Override
        public void onComplete(Object result) {
            // Do nothing
        }
        
        @Override
        public void onError(Exception error) {
            // Do nothing
        }
    }
}


package uk.ac.quadram.metapilot.web.websocket;

import java.time.LocalDateTime;

/**
 * Progress update message sent via WebSocket.
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
public class ProgressMessage {
    
    private String jobId;
    private String status;  // RUNNING, COMPLETED, FAILED, CANCELLED
    private double progress; // 0.0 to 1.0
    private String message;
    private LocalDateTime timestamp;
    
    public ProgressMessage() {
        this.timestamp = LocalDateTime.now();
    }
    
    public ProgressMessage(String jobId, String status, double progress, String message) {
        this.jobId = jobId;
        this.status = status;
        this.progress = progress;
        this.message = message;
        this.timestamp = LocalDateTime.now();
    }
    
    // Getters and setters
    public String getJobId() {
        return jobId;
    }
    
    public void setJobId(String jobId) {
        this.jobId = jobId;
    }
    
    public String getStatus() {
        return status;
    }
    
    public void setStatus(String status) {
        this.status = status;
    }
    
    public double getProgress() {
        return progress;
    }
    
    public void setProgress(double progress) {
        this.progress = progress;
    }
    
    public String getMessage() {
        return message;
    }
    
    public void setMessage(String message) {
        this.message = message;
    }
    
    public LocalDateTime getTimestamp() {
        return timestamp;
    }
    
    public void setTimestamp(LocalDateTime timestamp) {
        this.timestamp = timestamp;
    }
    
    @Override
    public String toString() {
        return "ProgressMessage{" +
                "jobId='" + jobId + '\'' +
                ", status='" + status + '\'' +
                ", progress=" + progress +
                ", message='" + message + '\'' +
                ", timestamp=" + timestamp +
                '}';
    }
}


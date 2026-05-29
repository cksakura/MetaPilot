package uk.ac.quadram.metapilot.web.websocket;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.messaging.simp.SimpMessagingTemplate;
import org.springframework.stereotype.Service;

/**
 * Service for broadcasting progress updates via WebSocket.
 * 
 * <p>This service uses Spring's SimpMessagingTemplate to send
 * progress messages to subscribed clients. Messages are sent to
 * the /topic/jobs/{jobId}/progress destination.</p>
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
@Service
public class ProgressBroadcastService {
    
    private static final Logger log = LoggerFactory.getLogger(ProgressBroadcastService.class);
    
    private final SimpMessagingTemplate messagingTemplate;
    
    public ProgressBroadcastService(SimpMessagingTemplate messagingTemplate) {
        this.messagingTemplate = messagingTemplate;
    }
    
    /**
     * Broadcast a progress update for a specific job.
     * 
     * @param jobId Job identifier
     * @param status Job status (RUNNING, COMPLETED, FAILED, CANCELLED)
     * @param progress Progress value (0.0 to 1.0)
     * @param message Status message
     */
    public void broadcastProgress(String jobId, String status, double progress, String message) {
        log.debug("Broadcasting progress for job {}: {} - {}% - {}", jobId, status, progress * 100, message);
        
        ProgressMessage progressMessage = new ProgressMessage(jobId, status, progress, message);
        
        // Send to /topic/jobs/{jobId}/progress
        String destination = "/topic/jobs/" + jobId + "/progress";
        messagingTemplate.convertAndSend(destination, progressMessage);
    }
    
    /**
     * Broadcast a message-only update (no progress percentage).
     * 
     * @param jobId Job identifier
     * @param message Status message
     */
    public void broadcastMessage(String jobId, String message) {
        broadcastProgress(jobId, "RUNNING", -1.0, message);
    }
    
    /**
     * Broadcast job completion.
     * 
     * @param jobId Job identifier
     * @param success Whether the job succeeded
     * @param message Result message
     */
    public void broadcastCompletion(String jobId, boolean success, String message) {
        String status = success ? "COMPLETED" : "FAILED";
        double progress = success ? 1.0 : 0.0;
        broadcastProgress(jobId, status, progress, message);
    }
    
    /**
     * Broadcast job cancellation.
     * 
     * @param jobId Job identifier
     * @param message Cancellation message
     */
    public void broadcastCancellation(String jobId, String message) {
        broadcastProgress(jobId, "CANCELLED", 0.0, message);
    }
}


package uk.ac.quadram.metapilot.web.rest;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import uk.ac.quadram.metapilot.web.service.WorkflowService;

import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

/**
 * REST Controller for MetaPilot job management.
 * 
 * This controller provides endpoints for:
 * - Submitting workflow jobs
 * - Monitoring job status
 * - Cancelling jobs
 * - Retrieving results
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
@RestController
@RequestMapping("/api/jobs")
public class JobController {
    
    @Autowired
    private WorkflowService workflowService;
    
    /**
     * Health check endpoint to verify the API is running.
     * 
     * @return Status information
     */
    @GetMapping("/status")
    public ResponseEntity<Map<String, String>> getStatus() {
        Map<String, String> status = new HashMap<>();
        status.put("status", "running");
        status.put("message", "MetaPilot Web API is operational");
        return ResponseEntity.ok(status);
    }
    
    /**
     * Submit a new workflow job.
     * 
     * @param request Job submission request containing parameter and resource file paths
     * @return Job submission response with job ID
     */
    @PostMapping("/submit")
    public ResponseEntity<Map<String, Object>> submitJob(@RequestBody Map<String, String> request) {
        String paramFile = request.get("paramFile");
        String resourceFile = request.get("resourceFile");
        
        if (paramFile == null || resourceFile == null) {
            Map<String, Object> error = new HashMap<>();
            error.put("error", "Missing required fields: paramFile, resourceFile");
            return ResponseEntity.badRequest().body(error);
        }
        
        // Generate unique job ID
        String jobId = UUID.randomUUID().toString();
        
        try {
            workflowService.submitJob(jobId, paramFile, resourceFile);
            
            Map<String, Object> response = new HashMap<>();
            response.put("jobId", jobId);
            response.put("status", "submitted");
            response.put("message", "Job submitted successfully");
            return ResponseEntity.ok(response);
        } catch (Exception e) {
            Map<String, Object> error = new HashMap<>();
            error.put("error", e.getMessage());
            return ResponseEntity.badRequest().body(error);
        }
    }
    
    /**
     * Get the status of a job.
     * 
     * @param jobId Job identifier
     * @return Job status information
     */
    @GetMapping("/{jobId}/status")
    public ResponseEntity<Map<String, Object>> getJobStatus(@PathVariable String jobId) {
        Map<String, Object> status = workflowService.getJobStatus(jobId);
        
        if (status.containsKey("notFound")) {
            return ResponseEntity.notFound().build();
        }
        
        return ResponseEntity.ok(status);
    }
    
    /**
     * Cancel a running job.
     * 
     * @param jobId Job identifier
     * @return Cancellation result
     */
    @PostMapping("/{jobId}/cancel")
    public ResponseEntity<Map<String, Object>> cancelJob(@PathVariable String jobId) {
        boolean cancelled = workflowService.cancelJob(jobId);
        
        Map<String, Object> response = new HashMap<>();
        response.put("jobId", jobId);
        response.put("cancelled", cancelled);
        
        if (cancelled) {
            response.put("message", "Job cancelled successfully");
            return ResponseEntity.ok(response);
        } else {
            response.put("message", "Job not found or already completed");
            return ResponseEntity.badRequest().body(response);
        }
    }
}

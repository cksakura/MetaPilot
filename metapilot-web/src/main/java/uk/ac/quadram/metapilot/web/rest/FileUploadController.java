package uk.ac.quadram.metapilot.web.rest;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.Map;

/**
 * REST Controller for file upload management.
 * 
 * Handles uploading of parameter files and resource files to the server.
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
@RestController
@RequestMapping("/api/files")
public class FileUploadController {
    
    private static final Logger log = LoggerFactory.getLogger(FileUploadController.class);
    
    @Value("${upload.directory:/opt/metapilot/data/uploads}")
    private String uploadDirectory;
    
    /**
     * Upload a parameter file.
     * 
     * @param file The uploaded file
     * @param workflowType Type of workflow (optional)
     * @return Upload result with file path
     */
    @PostMapping("/parameters")
    public ResponseEntity<Map<String, Object>> uploadParameterFile(
            @RequestParam("file") MultipartFile file,
            @RequestParam(value = "workflowType", required = false) String workflowType) {
        
        return uploadFile(file, "parameters", workflowType);
    }
    
    /**
     * Upload a resource configuration file.
     * 
     * @param file The uploaded file
     * @return Upload result with file path
     */
    @PostMapping("/resources")
    public ResponseEntity<Map<String, Object>> uploadResourceFile(
            @RequestParam("file") MultipartFile file) {
        
        return uploadFile(file, "resources", null);
    }
    
    /**
     * Upload a raw data file (mass spectrometry data).
     * 
     * @param file The uploaded file
     * @return Upload result with file path
     */
    @PostMapping("/raw")
    public ResponseEntity<Map<String, Object>> uploadRawFile(
            @RequestParam("file") MultipartFile file) {
        
        return uploadFile(file, "raw", null);
    }
    
    /**
     * Common file upload handler.
     * 
     * @param file The uploaded file
     * @param category File category (parameters, resources, raw)
     * @param subcategory Optional subcategory (workflow type)
     * @return Upload result
     */
    private ResponseEntity<Map<String, Object>> uploadFile(
            MultipartFile file, String category, String subcategory) {
        
        Map<String, Object> response = new HashMap<>();
        
        try {
            // Validate file
            if (file.isEmpty()) {
                response.put("error", "File is empty");
                return ResponseEntity.badRequest().body(response);
            }
            
            // Create upload directory structure
            String baseDir = uploadDirectory + "/" + category;
            if (subcategory != null) {
                baseDir += "/" + subcategory;
            }
            Path uploadPath = Paths.get(baseDir);
            
            if (!Files.exists(uploadPath)) {
                Files.createDirectories(uploadPath);
                log.info("Created directory: {}", uploadPath);
            }
            
            // Generate unique filename
            String originalFilename = file.getOriginalFilename();
            String extension = "";
            if (originalFilename != null && originalFilename.contains(".")) {
                extension = originalFilename.substring(originalFilename.lastIndexOf("."));
            }
            String filename = System.currentTimeMillis() + "_" + 
                            (originalFilename != null ? originalFilename : "file") + 
                            extension;
            
            // Save file
            Path filePath = uploadPath.resolve(filename);
            Files.copy(file.getInputStream(), filePath, StandardCopyOption.REPLACE_EXISTING);
            
            log.info("File uploaded: {} ({} bytes)", filePath, file.getSize());
            
            // Return success response
            response.put("success", true);
            response.put("filename", filename);
            response.put("path", filePath.toString());
            response.put("size", file.getSize());
            response.put("originalName", originalFilename);
            
            return ResponseEntity.ok(response);
            
        } catch (IOException e) {
            log.error("File upload failed", e);
            response.put("error", "Failed to upload file: " + e.getMessage());
            return ResponseEntity.internalServerError().body(response);
        }
    }
    
    /**
     * List uploaded parameter files.
     * 
     * @return List of parameter files
     */
    @GetMapping("/parameters")
    public ResponseEntity<Map<String, Object>> listParameterFiles() {
        return listFiles("parameters");
    }
    
    /**
     * List uploaded resource files.
     * 
     * @return List of resource files
     */
    @GetMapping("/resources")
    public ResponseEntity<Map<String, Object>> listResourceFiles() {
        return listFiles("resources");
    }
    
    /**
     * Common file listing handler.
     */
    private ResponseEntity<Map<String, Object>> listFiles(String category) {
        Map<String, Object> response = new HashMap<>();
        
        try {
            Path dirPath = Paths.get(uploadDirectory + "/" + category);
            
            if (!Files.exists(dirPath)) {
                response.put("files", new String[0]);
                return ResponseEntity.ok(response);
            }
            
            File[] files = dirPath.toFile().listFiles();
            if (files == null) {
                response.put("files", new String[0]);
                return ResponseEntity.ok(response);
            }
            
            // Return file information
            Map<String, Object>[] fileList = new HashMap[files.length];
            for (int i = 0; i < files.length; i++) {
                if (files[i].isFile()) {
                    Map<String, Object> fileInfo = new HashMap<>();
                    fileInfo.put("name", files[i].getName());
                    fileInfo.put("path", files[i].getAbsolutePath());
                    fileInfo.put("size", files[i].length());
                    fileInfo.put("modified", files[i].lastModified());
                    fileList[i] = fileInfo;
                }
            }
            
            response.put("files", fileList);
            return ResponseEntity.ok(response);
            
        } catch (Exception e) {
            log.error("Failed to list files", e);
            response.put("error", "Failed to list files: " + e.getMessage());
            return ResponseEntity.internalServerError().body(response);
        }
    }
}


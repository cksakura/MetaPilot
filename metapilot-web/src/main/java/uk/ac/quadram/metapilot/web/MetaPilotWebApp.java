package uk.ac.quadram.metapilot.web;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

/**
 * MetaPilot Web Application - Spring Boot Main Entry Point
 * 
 * <h2>Purpose</h2>
 * This is the main entry point for the MetaPilot web application designed
 * to run on Linux HPC environments. It provides a web-based UI for:
 * <ul>
 *   <li>Configuring MetaPilot parameters through a browser interface</li>
 *   <li>Submitting and monitoring workflow jobs</li>
 *   <li>Managing MAG databases</li>
 *   <li>Viewing results and reports</li>
 * </ul>
 * 
 * <h2>Technology Stack</h2>
 * <ul>
 *   <li>Spring Boot 3.2 for REST API and web server</li>
 *   <li>WebSocket for real-time progress updates</li>
 *   <li>Java 17 for core workflow execution</li>
 * </ul>
 * 
 * <h2>Usage</h2>
 * <pre>{@code
 *   // From command line:
 *   java -jar metapilot-web-1.1.0.jar
 *   
 *   // Access at: http://localhost:8080
 * }</pre>
 * 
 * @author MetaPilot Team
 * @version 1.0.0
 * @since 1.0.0
 */
@SpringBootApplication
public class MetaPilotWebApp {
    
    /**
     * Main method to launch the Spring Boot application.
     * 
     * @param args Command-line arguments passed to the application
     */
    public static void main(String[] args) {
        SpringApplication.run(MetaPilotWebApp.class, args);
    }
}

# MetaPilot Web Application

Web-based interface for MetaPilot workflows on Linux HPC environments.

## Overview

MetaPilot Web is a Spring Boot application that provides a REST API and web UI for submitting, monitoring, and managing MetaPilot workflow jobs. It is designed to run on Linux servers where users can access it via a web browser.

## Features

- ✅ **REST API** for programmatic job submission and monitoring
- ✅ **Web UI** for easy job configuration and submission
- ✅ **Job Management** with status tracking and cancellation
- ✅ **Real-time Progress** (coming soon via WebSocket)

## Quick Start

### Prerequisites

- Java 17 or higher
- Maven 3.6+
- MetaPilot Windows GUI JAR (for workflow execution)

### Building

```bash
# Build the web application
mvn clean package -pl metapilot-web -am -DskipTests

# The JAR will be in metapilot-web/target/
```

### Running

```bash
# Run the web server
java -jar metapilot-web/target/metapilot-web-1.0.0-beta.2.jar

# Or use Spring Boot Maven plugin
mvn spring-boot:run -pl metapilot-web
```

### Accessing the Application

1. Open a web browser and navigate to: `http://localhost:8080`
2. You should see the MetaPilot web interface
3. The API is available at: `http://localhost:8080/api/jobs`

## API Endpoints

### Health Check
```
GET /api/jobs/status
```

Returns the status of the API server.

### Submit Job
```
POST /api/jobs/submit
Content-Type: application/json

{
  "paramFile": "/path/to/parameters.json",
  "resourceFile": "/path/to/resources.params"
}
```

Submits a new workflow job for execution.

### Get Job Status
```
GET /api/jobs/{jobId}/status
```

Retrieves the status of a specific job.

### Cancel Job
```
POST /api/jobs/{jobId}/cancel
```

Cancels a running job.

## Configuration

Edit `metapilot-web/src/main/resources/application.properties` to customize:

- Server port (default: 8080)
- Server address (default: 0.0.0.0 - all interfaces)
- Logging levels
- File upload limits

## Project Structure

```
metapilot-web/
├── src/main/java/uk/ac/quadram/metapilot/web/
│   ├── MetaPilotWebApp.java          # Main application class
│   ├── rest/
│   │   └── JobController.java        # REST API endpoints
│   └── service/
│       └── WorkflowService.java      # Job execution service
├── src/main/resources/
│   ├── application.properties        # Configuration
│   └── static/
│       └── index.html                # Web UI
└── pom.xml                           # Maven configuration
```

## Current Limitations

⚠️ **Important**: The current implementation includes a placeholder workflow execution. To enable actual workflow execution:

1. **Option A**: Add JavaFX dependency to the web module (not recommended for server deployment)
2. **Option B**: Refactor tasks to use `ProgressCallback` interface (recommended long-term solution)
3. **Option C**: Execute workflows via separate CLI process (workaround)

## Development

### Adding New Endpoints

1. Add methods to `JobController` in the `rest` package
2. Update the service layer in `WorkflowService`
3. Test with curl or Postman

### Frontend Development

The frontend is a simple HTML/CSS/JavaScript application in `src/main/resources/static/index.html`. To enhance it:

1. Add more JavaScript for job monitoring
2. Implement WebSocket client for real-time updates
3. Add parameter file upload functionality

## Troubleshooting

### Port Already in Use

Change the port in `application.properties`:
```
server.port=8081
```

### Cannot Connect to API

Check that the Spring Boot application is running:
```bash
curl http://localhost:8080/api/jobs/status
```

### Workflow Execution Not Working

The current implementation includes placeholder execution logic. See "Current Limitations" above for solutions.

## Future Enhancements

- [ ] WebSocket support for real-time progress updates
- [ ] Parameter file upload via web UI
- [ ] Job history and result download
- [ ] MAG database management UI
- [ ] User authentication and authorization
- [ ] Docker containerization
- [ ] Kubernetes deployment manifests

## License

This project is part of MetaPilot, developed at Quadram Institute Bioscience.

## Contact

For questions or issues, contact the MetaPilot development team.


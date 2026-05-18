# MetaPilot HPC Linux Implementation Status

> Historical archive (2025-10). This file is retained for project history and is not current release guidance. For current CLI, workflow and deployment status, use [../README.md](../README.md).

## Project Structure

MetaPilot has been successfully restructured as a Maven multi-module project:

```
MetaPilot/
├── pom.xml                                    # Parent POM (packaging: pom)
├── metapilot-core/                            # Shared core logic (NO GUI dependencies)
│   └── src/main/java/uk/ac/quadram/metapilot/core/
│       ├── workflow/
│       │   └── ProgressCallback.java         # Generic progress reporting interface
│       └── (more core classes to be moved here)
│
├── metapilot-windows-gui/                     # Windows GUI application
│   └── src/main/java/uk/ac/quadram/metapilot/
│       ├── MetaPilotApp.java                 # Main Windows application
│       ├── gui/                               # JavaFX GUI components
│       ├── task/                              # All workflow tasks (currently here)
│       └── (640 Java files - complete Windows application)
│
└── metapilot-web/                             # Web application for Linux HPC
    └── src/main/java/uk/ac/quadram/metapilot/web/
        ├── MetaPilotWebApp.java              # Spring Boot main class
        └── rest/
            └── JobController.java            # REST API controller
```

## Completed Tasks ✅

### Phase 1: Module Restructuring
- ✅ Created Maven multi-module structure with parent POM
- ✅ Created `metapilot-core` module with shared dependencies
- ✅ Created `metapilot-windows-gui` module with all GUI code
- ✅ Created `metapilot-web` module placeholder
- ✅ Verified Windows GUI still builds and runs after restructuring

### Phase 2: Progress Callback Interface
- ✅ Created `ProgressCallback` interface in `metapilot-core`
  - Generic interface for reporting progress
  - Can be implemented by JavaFX, WebSocket, or CLI
  - Includes `NoOp` default implementation

### Phase 3: Spring Boot Web Application
- ✅ Configured `metapilot-web` with Spring Boot 3.2.0
- ✅ Added dependencies:
  - Spring Boot Web Starter
  - Spring Boot WebSocket (for real-time progress)
  - Spring Boot DevTools (for development)
- ✅ Created main application class `MetaPilotWebApp`
- ✅ Created REST controller `JobController` with health check
- ✅ Verified successful compilation and JAR generation (~128 MB)

## Build Status

### All Modules Build Successfully ✅

```bash
# Build all modules
mvn clean package -DskipTests

# Individual modules
mvn clean package -pl metapilot-windows-gui   # Windows GUI JAR
mvn clean package -pl metapilot-web           # Spring Boot JAR
mvn clean package -pl metapilot-core          # Core library
```

## Pending Tasks 📋

### Phase 2 (Partial)
- ⏳ Refactor `MGYGFtpTaskJavaFX` to non-JavaFX version using `ProgressCallback`
- ⏳ Move task classes to `metapilot-core` (carefully separating GUI from logic)

### Phase 4: Linux CLI Support
- ⏳ Create Maven build profiles for windows-gui and linux-cli
- ⏳ Implement OS detection and platform-aware `ToolPathResolver`
- ⏳ Expand `MetaPilotCLI` with full workflow support
- ⏳ Implement `BatchRunner` for processing job manifests

### Phase 5: Web UI Implementation
- ✅ Create workflow submission API endpoints
- ✅ Build web frontend for job configuration  
- ⏳ Implement WebSocket for real-time progress updates
- ⏳ Implement MAG database management UI

### Phase 6: Linux Tool Integration
- ⏳ Create Linux tool path defaults
- ⏳ Develop installation scripts for external tools
- ⏳ Create deployment shell scripts
- ⏳ Test FragPipe, Sage, DIA-NN workflows on Linux

### Phase 7: Documentation
- ⏳ Write Linux Installation Guide
- ⏳ Write HPC Deployment Guide
- ⏳ Write CLI Reference Manual

## Next Steps

### Immediate Next Step:
The web application is now functional with basic job submission and monitoring. Key remaining items:

1. **Workflow Execution**: Implement actual workflow execution (currently placeholder)
   - Add JavaFX dependency OR
   - Refactor tasks to use ProgressCallback OR  
   - Execute via separate CLI process

2. **WebSocket Progress**: Implement real-time progress updates from running jobs

3. **Parameter File Upload**: Allow users to upload parameter JSON files via web UI

4. **Linux Tool Integration**: Configure and test FragPipe, Sage, DIA-NN on Linux

5. **Testing**: Deploy to Linux VM and test with real workflows

## Architecture Decisions

### Why Multi-Module?
- **Separation of Concerns**: GUI code separate from core logic
- **Platform Flexibility**: Different frontends for Windows vs Linux
- **Shared Logic**: One core implementation for all platforms
- **Build Profiles**: Can build only what's needed for each platform

### Why Spring Boot for Web?
- **Embedded Server**: No separate Tomcat installation needed
- **REST API**: Standard HTTP endpoints for job submission
- **WebSocket Support**: Built-in real-time communication
- **Auto-Configuration**: Minimal boilerplate
- **Production Ready**: Well-tested, widely used in enterprise

### ProgressCallback Design
- **Platform Agnostic**: Works with JavaFX, WebSocket, Console
- **Decoupled**: Core tasks don't depend on UI framework
- **Testable**: Easy to mock for unit testing
- **Extensible**: Can add more callback types as needed

## Compatibility Matrix

| Component | Windows | Linux HPC | Notes |
|-----------|---------|-----------|-------|
| JavaFX GUI | ✅ | ❌ | Desktop-only |
| Spring Boot Web | ❌ | ✅ | Server-based |
| FragPipe | ✅ | ✅ | Native Linux support |
| Sage | ✅ | ✅ | Native Linux support |
| DIA-NN | ✅ | ✅ | Native Linux support |
| MaxQuant | ✅ | ❌ | Windows-only |
| pFind | ✅ | ❌ | Windows-only |

## File Counts

- **metapilot-windows-gui**: 640 Java files (complete application)
- **metapilot-core**: 2 Java files (ProgressCallback + structure)
- **metapilot-web**: 2 Java files (WebApp + JobController)

## Eclipse Integration

⚠️ **Eclipse users**: You need to refresh your Eclipse workspace:
1. Right-click on the `metapilot` project
2. Select **Maven → Update Project**
3. Check "Force Update of Snapshots/Releases"
4. Click OK

This will regenerate Eclipse project files for all three modules.

## Current Status Summary

✅ **Windows GUI**: Fully functional, builds successfully  
✅ **Maven Structure**: Multi-module setup working correctly  
✅ **Spring Boot**: Web application fully implemented with REST API  
✅ **Web UI**: Professional HTML interface created and working  
✅ **Job Management**: API endpoints for submission, monitoring, and cancellation  
⏳ **Workflow Execution**: Placeholder implementation (needs refactoring or JavaFX)  
⏳ **WebSocket Progress**: Real-time updates not yet implemented  
⏳ **Linux Tools**: Needs tool path configuration and testing  

## Contact

For questions or issues, contact the MetaPilot development team.

---
Last Updated: 2024
Version: 1.0.0-beta.3


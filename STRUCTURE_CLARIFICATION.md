# MetaPilot Project Structure Clarification

## Important: This is NOT the Project Structure

❌ **Wrong understanding**: "frontend/src and backend/src/main/java are empty"

✅ **Correct understanding**: The modules are:

```
MetaPilot (Root)
├── pom.xml                           # Parent POM
├── metapilot-core/                   # Module 1: Shared core
│   ├── pom.xml
│   └── src/main/java/uk/ac/quadram/metapilot/core/
│       └── workflow/
│           └── ProgressCallback.java
│
├── metapilot-windows-gui/           # Module 2: Windows Desktop App
│   ├── pom.xml
│   └── src/main/java/uk/ac/quadram/metapilot/
│       ├── MetaPilotApp.java
│       ├── MetaPilotCLI.java
│       └── (640 Java files!)
│
└── metapilot-web/                    # Module 3: Linux Web App
    ├── pom.xml
    └── src/main/java/uk/ac/quadram/metapilot/web/
        ├── MetaPilotWebApp.java
        ├── config/
        │   └── WebSocketConfig.java
        ├── rest/
        │   └── JobController.java
        ├── service/
        │   └── WorkflowService.java
        └── websocket/
            ├── ProgressBroadcastService.java
            └── ProgressMessage.java
    └── src/main/resources/
        ├── application.properties
        └── static/
            └── index.html
```

## What "Empty" Means

The **top-level `src/`** directory is **intentionally empty** now because we're using multi-module structure!

Previously (old structure):
- Top-level `pom.xml` with `<packaging>jar</packaging>`
- Top-level `src/main/java/` with all code

Now (new structure):
- Top-level `pom.xml` with `<packaging>pom</packaging>` (parent only)
- Each module has its own `src/main/java/` directory
- `metapilot-windows-gui/src/main/java/` = 640 Java files
- `metapilot-web/src/main/java/` = 6 Java files (complete web app)

## Building the Project

### Build All Modules (Default)
```bash
mvn clean package -DskipTests
```

This builds:
- metapilot-core (shared library)
- metapilot-windows-gui (Windows desktop app)
- metapilot-web (Linux web app)

### Build Specific Module
```bash
# Windows GUI only
mvn clean package -pl metapilot-windows-gui -am -DskipTests

# Web only
mvn clean package -pl metapilot-web -am -DskipTests

# Core only
mvn clean package -pl metapilot-core -am -DskipTests
```

### Using Build Scripts
```bash
# Windows
build_windows.bat

# Linux
./build_linux.sh
```

## Module Breakdown

### metapilot-core
**Purpose**: Shared core logic (no GUI dependencies)  
**Status**: 2 files (ProgressCallback + structure)  
**Future**: More shared code will move here

### metapilot-windows-gui
**Purpose**: Windows desktop application  
**Status**: 640 Java files (complete application)  
**Entry Point**: `MetaPilotApp.java`  
**Output**: `metapilot-windows-gui-1.0.0-beta.2-shaded.jar`

### metapilot-web
**Purpose**: Linux web server application  
**Status**: 6 Java files (complete web app)  
**Entry Point**: `MetaPilotWebApp.java`  
**Output**: `metapilot-web-1.0.0-beta.2.jar`  
**Components**:
- MetaPilotWebApp: Spring Boot main
- JobController: REST API
- WorkflowService: Job execution
- WebSocketConfig: Real-time updates
- ProgressBroadcastService: Progress messages
- ProgressMessage: Data class
- index.html: Web UI

## Confusion Points

### "Empty Directories"

The **empty directories** you saw (`metapilot-web/backend`, `metapilot-web/frontend`) were **leftovers from initial planning** and have been removed.

### Top-Level src/

The top-level `src/` directory is **legacy from the old single-module structure**. It's not used anymore but was kept for backward compatibility. All new code goes in module-specific `src/` directories.

### "Where is the Code?"

All code is in the **module directories**:

✅ **Windows GUI**: `metapilot-windows-gui/src/main/java/`  
✅ **Web App**: `metapilot-web/src/main/java/`  
⏳ **Core**: `metapilot-core/src/main/java/` (just started)

## Verification

To verify code exists:

```bash
# Windows GUI
ls metapilot-windows-gui/src/main/java/uk/ac/quadram/metapilot/
# Should show 640+ files

# Web App
ls metapilot-web/src/main/java/uk/ac/quadram/metapilot/web/
# Should show 6 files

# Check builds
mvn clean package -DskipTests
# Should create JARs in target/ directories
```

## Summary

**Your concern was valid** - but the structure is actually correct! 

- ❌ **NOT**: "frontend/src and backend/src" (those were just planning folders)
- ✅ **YES**: "metapilot-windows-gui/src" and "metapilot-web/src" (these have the code!)

The confusion came from seeing `metapilot-web/backend` and `metapilot-web/frontend` directories, which were empty planning structures. The real code is in `metapilot-web/src/main/java/` and that's **complete and working**.

Think of it as:
- **Module 1**: Windows Desktop (640 files in its src/)
- **Module 2**: Linux Web (6 files in its src/)
- **Module 3**: Core (1 file so far, will grow)

Each module is independent and builds its own JAR file!


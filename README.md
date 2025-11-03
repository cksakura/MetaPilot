# MetaPilot - Multi-Platform Metaproteomics Analysis

MetaPilot is a comprehensive metaproteomics analysis platform supporting both Windows desktop GUI and Linux web-based deployments.

## Project Structure

```
MetaPilot/
├── metapilot-core/          # Shared core logic (no GUI dependencies)
├── metapilot-windows-gui/   # Windows desktop application (JavaFX)
└── metapilot-web/           # Linux web application (Spring Boot)
```

## Quick Start

### Windows GUI

```bash
# Build
build_windows.bat

# Run
java -jar metapilot-windows-gui/target/metapilot-windows-gui-1.0.0-beta.2-shaded.jar
```

### Linux Web

```bash
# Build
./build_linux.sh

# Run
java -jar metapilot-web/target/metapilot-web-1.0.0-beta.2.jar

# Access
open http://localhost:8080
```

## Building

### All Modules
```bash
mvn clean package -DskipTests
```

### Specific Module
```bash
mvn clean package -pl metapilot-windows-gui -DskipTests  # Windows
mvn clean package -pl metapilot-web -DskipTests          # Web
```

## Features

### Windows GUI
- ✅ Full desktop application with 640 Java files
- ✅ Support for MaxQuant, pFind, FragPipe, Sage
- ✅ DDA and DIA workflows
- ✅ MAG database management

### Linux Web
- ✅ REST API for job submission
- ✅ Real-time progress updates via WebSocket
- ✅ Professional web UI
- ✅ Job management and monitoring
- ✅ Support for FragPipe, Sage, DIA-NN

## Documentation

- [Linux Installation Guide](docs/LINUX_INSTALLATION_GUIDE.md)
- [HPC Deployment Guide](docs/HPC_DEPLOYMENT_GUIDE.md)
- [CLI Reference Manual](docs/CLI_REFERENCE_MANUAL.md)
- [Implementation Status](IMPLEMENTATION_STATUS.md)
- [Project Summary](FINAL_STATUS.md)

## Requirements

- **Java**: JDK 17 or higher
- **Maven**: 3.6+
- **Operating Systems**: Windows 10+, Linux (Ubuntu/CentOS/Debian)

## Contact

MetaPilot Development Team - Quadram Institute Bioscience

## License

Copyright © Quadram Institute Bioscience

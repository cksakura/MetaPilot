# MetaPilot HPC Migration Plan
*Comprehensive Technical Report for Linux Server Deployment*

---

## Table of Contents
1. [Executive Summary](#executive-summary)
2. [Current System Analysis](#current-system-analysis)
3. [HPC Architecture Design](#hpc-architecture-design)
4. [File Format Compatibility](#file-format-compatibility)
5. [Implementation Strategy](#implementation-strategy)
6. [Linux File System Guide](#linux-file-system-guide)
7. [Technical Implementation](#technical-implementation)
8. [Testing & Validation](#testing--validation)
9. [Deployment Guide](#deployment-guide)
10. [Migration Timeline](#migration-timeline)

---

## Executive Summary

### Project Goal
Migrate MetaPilot from Windows-only desktop application to a **dual-mode system** supporting both desktop GUI and HPC Linux server execution.

### Key Strategy: Client-Server Architecture
- **Desktop Client**: Windows JavaFX GUI for parameter configuration and result visualization
- **HPC Server**: Headless Linux JAR for high-performance computation
- **Communication**: Parameter files (JSON) and data transfer via standard protocols

### Benefits
- ✅ **Scalability**: Leverage HPC cluster resources (hundreds of cores, TBs of RAM)
- ✅ **Compatibility**: No platform-specific rewrites needed
- ✅ **User Experience**: Familiar GUI for configuration, powerful backend for computation
- ✅ **Cost Effective**: Gradual migration, maintains existing functionality

---

## Current System Analysis

### Platform Dependencies Identified

#### **Critical Windows Dependencies:**
1. **Hardware ID Detection**: Uses `WMIC.exe` for license verification
2. **JavaFX GUI**: Extensive use of desktop UI components
3. **File Paths**: Windows-specific path handling (`C:\Program Files\...`)
4. **Process Execution**: Windows-specific external tool invocation

#### **File Format Challenges:**
1. **Thermo .raw Files**: Require Windows-only vendor libraries
2. **Vendor Formats**: .wiff (AB Sciex), .d (Bruker) are Windows-dependent
3. **Solution**: Convert to open formats (mzML) on Windows, transfer to Linux

#### **External Tool Dependencies:**
- MaxQuant: Windows-only
- msconvert: Cross-platform alternative available
- Search engines: Most have Linux versions

### Current CLI Capabilities
MetaPilot already has CLI interfaces:
- `MetaLabCli.java` - Legacy v1 workflow
- `MetaPilotTask.java` - Current pFind workflow  
- `PrideClusterCliMain.java` - Spectra clustering

---

## HPC Architecture Design

### Client-Server Workflow
```
┌─────────────────────┐    ┌──────────────────────────┐    ┌─────────────────────┐
│   Windows Desktop   │    │      HPC Linux Server   │    │    Results & Reports │
│  (Parameter Creator)│    │   (Computation Engine)  │    │                     │
├─────────────────────┤    ├──────────────────────────┤    ├─────────────────────┤
│ ✅ JavaFX GUI       │    │ ✅ Headless JAR         │    │ ✅ Web Reports      │
│ ✅ File Dialogs     │    │ ✅ CLI Interface        │    │ ✅ Results Download │
│ ✅ Parameter Config │───►│ ✅ Parameter Processing │───►│ ✅ Local Viewing    │
│ ✅ mzML Conversion  │    │ ✅ Heavy Computation     │    │                     │
│ ✅ Upload Helper    │    │ ✅ Batch Processing     │    │                     │
└─────────────────────┘    └──────────────────────────┘    └─────────────────────┘
```

### Dual-Mode Implementation

#### **Desktop JAR (metapilot-gui.jar):**
- Full JavaFX dependencies
- GUI parameter configuration
- File format conversion tools
- Data upload/download utilities
- Local result visualization

#### **Server JAR (metapilot-headless.jar):**
- Zero GUI dependencies
- CLI-only interface
- Core computational engines
- Multi-threading optimization
- Web-based report generation

---

## File Format Compatibility

### Current Format Support Analysis

#### **Problematic Formats (Windows-only):**
```java
// Current supported formats that are Windows-dependent:
.raw (Thermo)     → Requires Windows vendor libraries
.wiff (AB Sciex)  → Windows-only access
.d (Bruker)       → Windows-specific tools
```

#### **Linux-Compatible Formats:**
```java
// Cross-platform formats:
.mzML            → ✅ Open standard, full Linux support
.mzXML           → ✅ Cross-platform
.mgf             → ✅ Text-based, universal
```

### Conversion Strategy

#### **On Windows Desktop:**
1. Use existing msconvert integration
2. Convert .raw → .mzML before upload
3. Validate conversion quality
4. Transfer converted files to Linux

#### **Linux Processing:**
1. Process mzML files directly
2. Use cross-platform search engines
3. Generate standard output formats
4. Create web-based reports

### Alternative Linux Tools
- **ThermoRawFileParser**: Open-source .raw file reader for Linux
- **pwiz-docker**: Containerized ProteoWizard for Linux
- **OpenMS**: Comprehensive Linux-native tools

---

## Implementation Strategy

### Phase 1: Code Separation (Month 1)
#### **Objective**: Separate GUI and computational logic

#### **Tasks:**
1. **Create Headless Main Class**
   ```java
   // New entry point: MetaPilotHeadless.java
   public class MetaPilotHeadless {
       public static void main(String[] args) {
           // CLI-only execution
           // No JavaFX dependencies
       }
   }
   ```

2. **Dependency Separation**
   - Extract core computation classes
   - Remove JavaFX imports from computational modules
   - Create interface abstraction for progress reporting

3. **Build Configuration**
   ```xml
   <!-- Maven profile for headless build -->
   <profile>
       <id>headless</id>
       <dependencies>
           <!-- Exclude JavaFX dependencies -->
           <exclusions>
               <exclusion>
                   <groupId>org.openjfx</groupId>
                   <artifactId>*</artifactId>
               </exclusion>
           </exclusions>
       </dependencies>
   </profile>
   ```

### Phase 2: Enhanced CLI (Month 2)
#### **Objective**: Robust command-line interface

#### **Enhanced CLI Features:**
```bash
# Basic execution
java -jar metapilot-headless.jar --par config.json --output /results/ --threads 32

# Advanced options
java -jar metapilot-headless.jar \
  --par config.json \
  --output /scratch/results/ \
  --threads 64 \
  --memory 128G \
  --temp-dir /tmp/metapilot/ \
  --log-level INFO \
  --batch-mode \
  --web-report
```

#### **Parameter File Enhancement:**
```json
{
  "version": "HPC-1.0",
  "workflowType": "MaxQuantWorkflow",
  "hpcConfig": {
    "threads": 64,
    "memoryGB": 128,
    "tempDirectory": "/tmp/metapilot/",
    "scratchDirectory": "/scratch/user/",
    "logLevel": "INFO"
  },
  "inputFiles": [
    "/data/converted/sample1.mzML",
    "/data/converted/sample2.mzML"
  ],
  "databases": {
    "primary": "/shared/databases/uniprot_bacteria.fasta",
    "host": "/shared/databases/human.fasta"
  },
  "outputDirectory": "/results/experiment1/",
  "searchParameters": {
    "enzyme": "Trypsin",
    "missedCleavages": 2,
    "fixedMods": ["Carbamidomethyl (C)"],
    "variableMods": ["Oxidation (M)", "Acetyl (Protein N-term)"],
    "tolerance": {
      "precursor": "20 ppm",
      "fragment": "0.02 Da"
    }
  },
  "fdrSettings": {
    "psmFdr": 0.01,
    "proteinFdr": 0.01
  }
}
```

### Phase 3: Linux Compatibility (Month 3)
#### **Objective**: Full Linux server compatibility

#### **Platform Abstraction:**
```java
// Abstract platform-specific operations
public interface PlatformManager {
    String getHardwareId();          // License verification
    File getTempDirectory();         // Temporary file location
    ProcessBuilder createProcess();   // External tool execution
}

// Linux implementation
public class LinuxPlatformManager implements PlatformManager {
    @Override
    public String getHardwareId() {
        // Use /sys/class/dmi/id/ or /proc/cpuinfo
        // Alternative: MAC address, disk serial
    }
}
```

#### **License System Adaptation:**
```java
// Alternative hardware ID methods for Linux
public class LinuxHardwareId {
    public static String getSystemId() {
        // Method 1: DMI system UUID
        String uuid = readFile("/sys/class/dmi/id/product_uuid");
        
        // Method 2: CPU info
        String cpuInfo = readFile("/proc/cpuinfo");
        
        // Method 3: Network interface MAC
        String macAddress = getNetworkMac();
        
        return combineIds(uuid, cpuInfo, macAddress);
    }
}
```

### Phase 4: Web Interface Foundation (Month 4)
#### **Objective**: Prepare for web-based interface

#### **Web Report Generation:**
```java
// Generate self-contained HTML reports
public class WebReportGenerator {
    public void generateReport(Results results, File outputDir) {
        // Create interactive HTML5 reports
        // Embed CSS, JavaScript, and data
        // No server dependency - static files
    }
}
```

#### **REST API Foundation:**
```java
// Optional: Simple REST API for job monitoring
@RestController
public class MetaPilotJobController {
    @GetMapping("/jobs/{id}/status")
    public JobStatus getJobStatus(@PathVariable String id) {
        // Return job progress and status
    }
    
    @PostMapping("/jobs/submit")
    public String submitJob(@RequestBody JobParameters params) {
        // Submit new analysis job
    }
}
```

---

## Linux File System Guide

### Essential Directory Structure
```bash
/opt/metapilot/                    # Main installation
├── bin/
│   ├── metapilot-headless.jar    # Server JAR
│   ├── metapilot-gui.jar         # Desktop JAR
│   └── run_metapilot.sh          # Launch script
├── lib/                          # Java libraries
├── config/
│   ├── default_config.yaml
│   └── resources.params
├── databases/                    # Reference databases
├── tools/                        # External binaries
└── logs/                         # System logs

/home/username/metapilot/          # User workspace
├── parameters/                   # Parameter files from Windows
├── data/                         # Input data (mzML files)
├── projects/                     # Individual projects
├── results/                      # Analysis outputs
└── temp/                         # Temporary files

/scratch/metapilot/               # High-performance storage
├── input/                        # Large input files
├── working/                      # Temporary processing
└── output/                       # Results before transfer
```

### Path Translation Examples
```java
// Windows → Linux path translation
public class PathTranslator {
    public static String translatePath(String windowsPath) {
        // C:\MetaPilot\data\sample.raw → /home/user/metapilot/data/sample.mzML
        // C:\Program Files\MetaPilot\Resources\ → /opt/metapilot/config/
        // %TEMP%\metapilot\ → /tmp/metapilot/
    }
}
```

### Essential Linux Commands
```bash
# Navigation
pwd                              # Show current directory
ls -la                          # List files with details
cd /path/to/directory           # Change directory
find / -name "*.jar" 2>/dev/null # Find files

# File operations
mkdir -p /path/to/directory     # Create directory structure
cp source destination           # Copy files
mv old_name new_name           # Move/rename
rm -rf directory/              # Delete directory (CAREFUL!)

# Permissions
chmod +x script.sh             # Make executable
chmod 755 metapilot.jar       # Set permissions

# Process monitoring
top                            # System monitor
ps aux | grep metapilot       # Find MetaPilot processes
kill -9 PID                   # Force kill process

# File transfer
scp local_file user@server:/remote/path/     # Upload
scp user@server:/remote/file local_path/     # Download
rsync -av local/ user@server:/remote/        # Sync directories
```

---

## Technical Implementation

### Build System Configuration

#### **Maven Configuration:**
```xml
<profiles>
    <!-- Desktop build with JavaFX -->
    <profile>
        <id>desktop</id>
        <activation>
            <activeByDefault>true</activeByDefault>
        </activation>
        <dependencies>
            <dependency>
                <groupId>org.openjfx</groupId>
                <artifactId>javafx-controls</artifactId>
                <version>21.0.7</version>
            </dependency>
        </dependencies>
        <build>
            <finalName>metapilot-gui</finalName>
        </build>
    </profile>
    
    <!-- Headless build for servers -->
    <profile>
        <id>headless</id>
        <dependencies>
            <!-- Exclude JavaFX -->
            <dependency>
                <groupId>org.openjfx</groupId>
                <artifactId>javafx-controls</artifactId>
                <scope>provided</scope>
            </dependency>
        </dependencies>
        <build>
            <finalName>metapilot-headless</finalName>
            <plugins>
                <plugin>
                    <groupId>org.springframework.boot</groupId>
                    <artifactId>spring-boot-maven-plugin</artifactId>
                    <configuration>
                        <mainClass>uk.ac.quadram.metapilot.HeadlessMain</mainClass>
                    </configuration>
                </plugin>
            </plugins>
        </build>
    </profile>
</profiles>
```

#### **Gradle Alternative:**
```groovy
// Desktop configuration
configurations {
    desktop {
        extendsFrom implementation
    }
    headless {
        extendsFrom implementation
        exclude group: 'org.openjfx'
    }
}

task desktopJar(type: Jar) {
    classifier = 'gui'
    from sourceSets.main.output
    configurations = [project.configurations.desktop]
}

task headlessJar(type: Jar) {
    classifier = 'headless'
    from sourceSets.main.output
    configurations = [project.configurations.headless]
    manifest {
        attributes 'Main-Class': 'uk.ac.quadram.metapilot.HeadlessMain'
    }
}
```

### Core Classes Architecture

#### **Platform Abstraction:**
```java
// Core abstraction interface
public interface PlatformProvider {
    ProgressReporter createProgressReporter();
    FileManager createFileManager();
    ProcessExecutor createProcessExecutor();
    LicenseValidator createLicenseValidator();
}

// Desktop implementation
public class DesktopPlatformProvider implements PlatformProvider {
    @Override
    public ProgressReporter createProgressReporter() {
        return new JavaFXProgressReporter();
    }
}

// Server implementation  
public class ServerPlatformProvider implements PlatformProvider {
    @Override
    public ProgressReporter createProgressReporter() {
        return new ConsoleProgressReporter();
    }
}
```

#### **Progress Reporting:**
```java
// Abstract progress interface
public interface ProgressReporter {
    void setProgress(double progress);
    void setMessage(String message);
    void setStatus(TaskStatus status);
}

// Console implementation for headless
public class ConsoleProgressReporter implements ProgressReporter {
    @Override
    public void setProgress(double progress) {
        System.out.printf("Progress: %.1f%%\r", progress * 100);
    }
    
    @Override
    public void setMessage(String message) {
        System.out.println(LocalDateTime.now() + " - " + message);
    }
}

// JavaFX implementation for desktop
public class JavaFXProgressReporter implements ProgressReporter {
    private ProgressBar progressBar;
    private Label messageLabel;
    
    @Override
    public void setProgress(double progress) {
        Platform.runLater(() -> progressBar.setProgress(progress));
    }
}
```

#### **Headless Main Class:**
```java
public class MetaPilotHeadless {
    private static final Logger LOGGER = LogManager.getLogger(MetaPilotHeadless.class);
    
    public static void main(String[] args) {
        // Parse command line arguments
        CommandLineOptions options = parseArguments(args);
        
        if (options.hasHelp()) {
            printUsage();
            return;
        }
        
        // Validate environment
        validateLinuxEnvironment();
        
        // Load parameter file
        MetaBaseParameter parameters = loadParameters(options.getParameterFile());
        
        // Set up platform provider
        PlatformProvider platform = new ServerPlatformProvider();
        
        // Initialize headless execution context
        HeadlessExecutionContext context = new HeadlessExecutionContext(
            parameters, platform, options
        );
        
        // Execute workflow
        WorkflowExecutor executor = new WorkflowExecutor(context);
        executor.execute();
    }
    
    private static void validateLinuxEnvironment() {
        // Check Java version
        String javaVersion = System.getProperty("java.version");
        LOGGER.info("Java version: " + javaVersion);
        
        // Check available memory
        long maxMemory = Runtime.getRuntime().maxMemory();
        LOGGER.info("Max memory: " + (maxMemory / 1024 / 1024) + " MB");
        
        // Check temporary directory
        String tempDir = System.getProperty("java.io.tmpdir");
        File temp = new File(tempDir);
        if (!temp.exists() || !temp.canWrite()) {
            throw new RuntimeException("Cannot access temp directory: " + tempDir);
        }
        
        // Validate required tools
        validateExternalTools();
    }
    
    private static void validateExternalTools() {
        // Check for required external tools
        String[] requiredTools = {"java", "python3"};
        
        for (String tool : requiredTools) {
            if (!isToolAvailable(tool)) {
                LOGGER.warn("Tool not found: " + tool);
            }
        }
    }
}
```

### Parameter File Schema

#### **Enhanced Parameter Structure:**
```json
{
  "$schema": "https://metapilot.org/schemas/hpc-parameters-v1.0.json",
  "metaPilotVersion": "1.0.0",
  "parameterVersion": "HPC-1.0",
  "createdBy": "MetaPilot GUI v1.0.0",
  "createdDate": "2024-01-15T10:30:00Z",
  "platform": {
    "type": "HPC",
    "os": "linux",
    "architecture": "x86_64"
  },
  "execution": {
    "mode": "batch",
    "priority": "normal",
    "wallTime": "24:00:00",
    "nodes": 1,
    "coresPerNode": 32,
    "memoryPerNode": "128GB",
    "scratchSize": "1TB"
  },
  "workflow": {
    "type": "MaxQuantWorkflow",
    "steps": [
      "preprocessing",
      "search",
      "postprocessing",
      "taxonomy",
      "functional"
    ],
    "parallelization": "auto"
  },
  "input": {
    "dataFiles": [
      {
        "path": "/scratch/data/sample1.mzML",
        "type": "mzML",
        "experiment": "Control_1",
        "replicate": 1
      }
    ],
    "databases": {
      "primary": "/shared/databases/uniprot_bacteria_2024.fasta",
      "host": "/shared/databases/human_proteome.fasta",
      "contaminants": "/shared/databases/contaminants.fasta"
    }
  },
  "output": {
    "baseDirectory": "/results/experiment_20240115/",
    "format": "standard",
    "compression": "gzip",
    "reports": {
      "html": true,
      "pdf": false,
      "json": true
    }
  },
  "searchParameters": {
    "enzyme": {
      "name": "Trypsin",
      "cleavageRule": "(?<=[KR])(?!P)",
      "missedCleavages": 2
    },
    "modifications": {
      "fixed": [
        {
          "name": "Carbamidomethyl",
          "composition": "H(3) C(2) N O",
          "residues": ["C"],
          "position": "anywhere"
        }
      ],
      "variable": [
        {
          "name": "Oxidation",
          "composition": "O",
          "residues": ["M"],
          "position": "anywhere"
        }
      ]
    },
    "tolerance": {
      "precursor": {
        "value": 20,
        "unit": "ppm"
      },
      "fragment": {
        "value": 0.02,
        "unit": "Da"
      }
    },
    "fdr": {
      "psm": 0.01,
      "protein": 0.01,
      "site": 0.05
    }
  },
  "quantification": {
    "method": "LFQ",
    "normalization": "median",
    "matchBetweenRuns": true,
    "matchTimeWindow": 0.7,
    "alignmentTimeWindow": 20,
    "minRatioCount": 2
  },
  "taxonomy": {
    "enabled": true,
    "database": "/shared/taxonomy/ncbi_taxonomy.db",
    "ranks": ["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
    "lca": {
      "minPeptides": 2,
      "confidenceThreshold": 0.95
    }
  },
  "functional": {
    "enabled": true,
    "databases": {
      "go": "/shared/functional/go.obo",
      "kegg": "/shared/functional/kegg.db",
      "cog": "/shared/functional/cog.db"
    },
    "enrichment": {
      "method": "fisher",
      "correction": "benjamini"
    }
  },
  "performance": {
    "threading": {
      "search": 16,
      "postprocessing": 8,
      "io": 4
    },
    "memory": {
      "jvmHeap": "64G",
      "searchEngine": "32G",
      "bufferSize": "1G"
    },
    "io": {
      "tempDirectory": "/tmp/metapilot/",
      "bufferSize": "64MB",
      "compressionLevel": 6
    }
  },
  "logging": {
    "level": "INFO",
    "file": "/logs/metapilot_20240115.log",
    "maxSize": "100MB",
    "retention": "30 days"
  }
}
```

---

## Testing & Validation

### Testing Strategy

#### **Phase 1: Unit Testing**
```java
// Test platform abstraction
@Test
public void testPlatformAbstraction() {
    PlatformProvider desktop = new DesktopPlatformProvider();
    PlatformProvider server = new ServerPlatformProvider();
    
    // Both should implement same interface
    assertNotNull(desktop.createProgressReporter());
    assertNotNull(server.createProgressReporter());
}

// Test parameter file parsing
@Test  
public void testParameterParsing() {
    File paramFile = new File("test_parameters.json");
    MetaBaseParameter params = loadParameters(paramFile);
    
    assertNotNull(params);
    assertEquals("HPC-1.0", params.getVersion());
}
```

#### **Phase 2: Integration Testing**
```bash
#!/bin/bash
# Test script for Linux compatibility

# Test 1: Basic execution
echo "Testing basic execution..."
java -jar metapilot-headless.jar --help

# Test 2: Parameter file validation
echo "Testing parameter validation..."
java -jar metapilot-headless.jar --validate test_params.json

# Test 3: Small dataset processing
echo "Testing small dataset..."
java -jar metapilot-headless.jar \
  --par small_test.json \
  --output /tmp/test_output/ \
  --threads 2

# Test 4: Memory stress test
echo "Testing memory usage..."
java -Xmx4G -jar metapilot-headless.jar \
  --par memory_test.json \
  --output /tmp/memory_test/
```

#### **Phase 3: HPC Validation**
```bash
#!/bin/bash
#SBATCH --job-name=metapilot-test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00
#SBATCH --mem=64GB

# Load required modules
module load java/11
module load python/3.8

# Set up environment
export JAVA_OPTS="-Xmx60G -XX:+UseG1GC"
export METAPILOT_HOME=/opt/metapilot

# Run test
cd $SLURM_SUBMIT_DIR
java $JAVA_OPTS -jar $METAPILOT_HOME/bin/metapilot-headless.jar \
  --par hpc_test_parameters.json \
  --output $SCRATCH/metapilot_test/ \
  --threads $SLURM_NTASKS
```

### Validation Criteria

#### **Functional Requirements:**
- ✅ All existing workflows execute correctly
- ✅ Results identical to Windows version
- ✅ Parameter files validate properly
- ✅ Error handling robust
- ✅ Log files comprehensive

#### **Performance Requirements:**
- ✅ Linear scaling with thread count
- ✅ Memory usage within bounds
- ✅ I/O performance acceptable
- ✅ Network transfer stable

#### **Compatibility Requirements:**
- ✅ Java 11+ on Linux
- ✅ Common HPC schedulers (SLURM, PBS, SGE)
- ✅ Standard Linux distributions
- ✅ Container environments (Docker, Singularity)

---

## Deployment Guide

### HPC Installation

#### **System Requirements:**
```bash
# Minimum requirements
Java 11+                     # OpenJDK or Oracle JDK
RAM: 8GB                    # Minimum, 64GB+ recommended
Storage: 100GB              # For software + temporary files
CPU: 4+ cores               # 16+ cores recommended

# Recommended HPC environment
Scheduler: SLURM/PBS/SGE    # Job scheduling
Filesystem: Lustre/GPFS     # High-performance filesystem
Network: InfiniBand         # High-speed interconnect
Containers: Singularity     # For reproducibility
```

#### **Installation Script:**
```bash
#!/bin/bash
# MetaPilot HPC Installation Script

set -e

# Configuration
INSTALL_DIR="/opt/metapilot"
USER_DIR="/home/$USER/metapilot"
SHARED_DB_DIR="/shared/databases"

echo "Installing MetaPilot HPC version..."

# Create directories
sudo mkdir -p $INSTALL_DIR/{bin,lib,config,logs}
mkdir -p $USER_DIR/{parameters,data,results,logs}

# Download and install
cd $INSTALL_DIR
sudo wget https://github.com/your-org/metapilot/releases/download/v1.0.0/metapilot-headless.jar
sudo wget https://github.com/your-org/metapilot/releases/download/v1.0.0/metapilot-gui.jar

# Set permissions
sudo chmod +x bin/*.jar
sudo chown -R root:metapilot $INSTALL_DIR
sudo chmod -R 755 $INSTALL_DIR

# Create launch script
sudo tee $INSTALL_DIR/bin/metapilot << 'EOF'
#!/bin/bash
METAPILOT_HOME="/opt/metapilot"
JAVA_OPTS="-Xmx${METAPILOT_MEMORY:-8G} -XX:+UseG1GC"

cd $METAPILOT_HOME
exec java $JAVA_OPTS -jar bin/metapilot-headless.jar "$@"
EOF

sudo chmod +x $INSTALL_DIR/bin/metapilot

# Add to PATH
echo 'export PATH=$PATH:/opt/metapilot/bin' >> ~/.bashrc

# Download reference databases (optional)
if [ -d "$SHARED_DB_DIR" ]; then
    echo "Setting up shared databases..."
    # Add database download/setup scripts here
fi

echo "MetaPilot HPC installation completed!"
echo "Test with: metapilot --help"
```

#### **SLURM Job Template:**
```bash
#!/bin/bash
#SBATCH --job-name=metapilot
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --output=metapilot_%j.out
#SBATCH --error=metapilot_%j.err

# Load required modules
module load java/11
module load python/3.8

# Set environment variables
export JAVA_OPTS="-Xmx120G -XX:+UseG1GC -XX:ParallelGCThreads=4"
export METAPILOT_HOME=/opt/metapilot
export SCRATCH_DIR=$SCRATCH/metapilot_$SLURM_JOB_ID

# Create scratch directory
mkdir -p $SCRATCH_DIR

# Copy input files to scratch
echo "Copying input files to scratch..."
cp $SLURM_SUBMIT_DIR/parameters.json $SCRATCH_DIR/
cp -r $SLURM_SUBMIT_DIR/data/ $SCRATCH_DIR/

# Run MetaPilot
echo "Starting MetaPilot analysis..."
cd $SCRATCH_DIR

$METAPILOT_HOME/bin/metapilot \
  --par parameters.json \
  --output results/ \
  --threads $SLURM_NTASKS \
  --memory ${SLURM_MEM_PER_NODE}MB \
  --temp-dir $SCRATCH_DIR/temp/

# Copy results back
echo "Copying results back..."
cp -r results/ $SLURM_SUBMIT_DIR/

# Cleanup
rm -rf $SCRATCH_DIR

echo "MetaPilot analysis completed!"
```

### Desktop Configuration

#### **Parameter Export Wizard:**
```java
// Enhanced GUI for HPC parameter generation
public class HpcParameterExportWizard extends Dialog<MetaBaseParameter> {
    
    @FXML private ComboBox<String> hpcSystemType;
    @FXML private TextField nodeCount;
    @FXML private TextField coresPerNode;
    @FXML private TextField memoryPerNode;
    @FXML private TextField wallTime;
    @FXML private TextField scratchDirectory;
    @FXML private TextField outputDirectory;
    
    public void initialize() {
        hpcSystemType.getItems().addAll("SLURM", "PBS", "SGE", "Generic");
        
        // Set reasonable defaults
        nodeCount.setText("1");
        coresPerNode.setText("32");
        memoryPerNode.setText("128GB");
        wallTime.setText("24:00:00");
    }
    
    public MetaBaseParameter generateHpcParameters() {
        MetaBaseParameter params = getCurrentParameters();
        
        // Add HPC-specific configuration
        params.setThreadCount(Integer.parseInt(coresPerNode.getText()));
        params.setOutputDirectory(outputDirectory.getText());
        
        // Validate paths for Linux
        validateLinuxPaths(params);
        
        return params;
    }
}
```

#### **Data Upload Assistant:**
```java
public class HpcUploadAssistant {
    
    public void uploadData(List<File> files, String serverInfo) {
        ProgressDialog progress = new ProgressDialog();
        
        Task<Void> uploadTask = new Task<Void>() {
            @Override
            protected Void call() throws Exception {
                for (int i = 0; i < files.size(); i++) {
                    File file = files.get(i);
                    updateMessage("Uploading " + file.getName());
                    updateProgress(i, files.size());
                    
                    // Use SCP or SFTP to upload
                    uploadFile(file, serverInfo);
                }
                return null;
            }
        };
        
        progress.showAndWait(uploadTask);
    }
    
    private void uploadFile(File localFile, String serverInfo) {
        // Implement SCP/SFTP upload logic
        // Could use JSch library or ProcessBuilder with scp command
    }
}
```

---

## Migration Timeline

### Month 1: Foundation
**Week 1-2: Analysis & Planning**
- ✅ Complete dependency analysis
- ✅ Design platform abstraction layer
- ✅ Create project structure for dual builds

**Week 3-4: Code Separation**
- ✅ Extract GUI-independent computational core
- ✅ Create headless main class
- ✅ Set up build profiles (Maven/Gradle)
- ✅ Basic CLI argument parsing

### Month 2: CLI Enhancement
**Week 1-2: Enhanced Command Line**
- ✅ Robust parameter file validation
- ✅ Comprehensive error handling
- ✅ Progress reporting for console
- ✅ Logging system integration

**Week 3-4: Testing & Validation**
- ✅ Unit tests for core functionality
- ✅ Integration tests on Linux VMs
- ✅ Performance benchmarking
- ✅ Memory usage optimization

### Month 3: Linux Compatibility
**Week 1-2: Platform Adaptation**
- ✅ Linux-specific path handling
- ✅ Alternative hardware ID methods
- ✅ External tool integration (Linux versions)
- ✅ File permission handling

**Week 3-4: HPC Integration**
- ✅ SLURM/PBS job script templates
- ✅ Scratch filesystem handling
- ✅ Module system integration
- ✅ Performance tuning for HPC

### Month 4: Polish & Documentation
**Week 1-2: User Experience**
- ✅ Enhanced parameter export wizard
- ✅ Upload/download assistants
- ✅ Web-based progress monitoring
- ✅ Comprehensive error messages

**Week 3-4: Documentation & Training**
- ✅ User documentation for HPC workflow
- ✅ Administrator installation guides
- ✅ Video tutorials for migration
- ✅ Best practices documentation

### Month 5: Deployment & Support
**Week 1-2: Beta Testing**
- ✅ Deploy on target HPC systems
- ✅ User acceptance testing
- ✅ Performance validation
- ✅ Bug fixes and optimization

**Week 3-4: Production Release**
- ✅ Production deployment
- ✅ User training sessions
- ✅ Support documentation
- ✅ Monitoring and maintenance procedures

---

## Risk Mitigation

### Technical Risks

#### **Risk: JavaFX Dependency Removal**
- **Mitigation**: Incremental separation, maintain desktop version
- **Fallback**: Use headless Java with Xvfb if needed

#### **Risk: External Tool Compatibility**
- **Mitigation**: Identify Linux alternatives early
- **Fallback**: Container-based tool execution

#### **Risk: Performance Degradation**
- **Mitigation**: Extensive benchmarking and optimization
- **Fallback**: Hybrid approach with critical sections optimized

### Operational Risks

#### **Risk: User Adoption Resistance**
- **Mitigation**: Maintain familiar GUI, gradual transition
- **Fallback**: Support both modes indefinitely

#### **Risk: HPC Environment Variations**
- **Mitigation**: Support major schedulers, containerization
- **Fallback**: Generic Linux version with manual setup

#### **Risk: Data Transfer Complexity**
- **Mitigation**: Automated tools, clear documentation
- **Fallback**: Standard file transfer protocols

---

## Success Metrics

### Technical Metrics
- ✅ **Performance**: 90%+ of Windows performance on equivalent hardware
- ✅ **Scalability**: Linear scaling up to 64+ cores
- ✅ **Reliability**: <1% job failure rate
- ✅ **Compatibility**: Works on 95%+ of HPC systems

### User Experience Metrics
- ✅ **Adoption**: 50%+ of users try HPC version within 6 months
- ✅ **Satisfaction**: 80%+ user satisfaction score
- ✅ **Support**: <24 hour response time for issues
- ✅ **Documentation**: Comprehensive guides and tutorials

### Business Metrics
- ✅ **Throughput**: 10x increase in analysis capacity
- ✅ **Cost**: Reduced per-analysis computational cost
- ✅ **Accessibility**: Available to users without Windows servers
- ✅ **Competitiveness**: Feature parity with cloud solutions

---

## Conclusion

This migration plan provides a comprehensive roadmap for transforming MetaPilot from a Windows-only desktop application to a powerful, scalable HPC solution. The dual-mode approach ensures we maintain existing user workflows while opening new possibilities for high-throughput analysis.

The key to success will be:
1. **Incremental implementation** to minimize disruption
2. **Thorough testing** across diverse HPC environments  
3. **Clear documentation** for smooth user transition
4. **Strong support** during the migration period

With this plan, MetaPilot will be positioned as a leading solution for large-scale metaproteomics analysis in both academic and commercial HPC environments.

---

*Document Version: 1.0*  
*Last Updated: 2024-01-15*  
*Author: MetaPilot Development Team* 
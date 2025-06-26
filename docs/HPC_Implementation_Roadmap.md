# MetaPilot HPC Implementation Roadmap
*Quick Start Guide for HPC Migration Project*

## 🎯 Project Overview

**Goal**: Create dual-mode MetaPilot (Desktop GUI + HPC CLI) for Linux server deployment

**Strategy**: Client-Server architecture with JSON parameter files

**Timeline**: 5 months

---

## 🚀 Quick Start Checklist

### Before You Begin
- [ ] Review the comprehensive [Technical Report](MetaPilot_HPC_Migration_Plan.md)
- [ ] Set up Linux development environment (VM or WSL)
- [ ] Identify target HPC system specifications
- [ ] Test current CLI functionality: `java -jar metapilot.jar --help`

---

## 📋 Phase-by-Phase Implementation

### Phase 1: Code Separation (Month 1)
**Priority**: Create headless JAR without JavaFX dependencies

#### Critical Tasks:
```bash
# 1. Create new main class
src/main/java/uk/ac/quadram/metapilot/HeadlessMain.java

# 2. Extract computational core
src/main/java/uk/ac/quadram/metapilot/core/headless/

# 3. Set up build profiles
pom.xml (add headless profile)

# 4. Test basic functionality
mvn clean compile -P headless
java -jar target/metapilot-headless.jar --help
```

#### Success Criteria:
- ✅ Headless JAR builds without JavaFX
- ✅ Basic CLI commands work
- ✅ Core algorithms execute correctly

---

### Phase 2: Enhanced CLI (Month 2)
**Priority**: Robust command-line interface

#### Critical Tasks:
```bash
# Enhanced CLI with comprehensive options
java -jar metapilot-headless.jar \
  --par config.json \
  --output /results/ \
  --threads 32 \
  --memory 64G \
  --log-level INFO
```

#### Key Files to Modify:
- `MetaLabCli.java` - Enhance existing CLI
- `MetaPilotTask.java` - Add headless execution
- Parameter validation classes

#### Success Criteria:
- ✅ Full parameter file support
- ✅ Comprehensive error handling
- ✅ Progress reporting to console

---

### Phase 3: Linux Compatibility (Month 3)
**Priority**: Platform abstraction and Linux-specific adaptations

#### Critical Tasks:
```java
// Platform abstraction interface
public interface PlatformProvider {
    String getHardwareId();
    ProcessBuilder createProcess(String... command);
    File getTempDirectory();
}

// Linux implementation
public class LinuxPlatformProvider implements PlatformProvider {
    // Linux-specific implementations
}
```

#### Success Criteria:
- ✅ License system works on Linux
- ✅ External tools integrate properly
- ✅ File paths handle correctly

---

### Phase 4: HPC Integration (Month 4)
**Priority**: HPC scheduler integration and optimization

#### SLURM Job Template:
```bash
#!/bin/bash
#SBATCH --job-name=metapilot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mem=128GB

module load java/11
java -Xmx120G -jar metapilot-headless.jar \
  --par parameters.json \
  --output results/ \
  --threads $SLURM_NTASKS
```

#### Success Criteria:
- ✅ SLURM/PBS integration
- ✅ Scratch filesystem handling
- ✅ Performance optimization

---

### Phase 5: Desktop Integration (Month 5)
**Priority**: Enhanced GUI for HPC workflow

#### Desktop Features:
- Parameter export wizard for HPC
- File upload/download assistants
- Job monitoring dashboard
- Results visualization

#### Success Criteria:
- ✅ Seamless parameter generation
- ✅ Easy data transfer tools
- ✅ Remote job monitoring

---

## 🛠️ Technical Implementation Notes

### Build System Setup
```xml
<!-- Maven profiles for dual builds -->
<profiles>
  <profile>
    <id>desktop</id>
    <!-- Include JavaFX dependencies -->
  </profile>
  <profile>
    <id>headless</id>
    <!-- Exclude JavaFX, CLI-only -->
  </profile>
</profiles>
```

### Parameter File Structure
```json
{
  "version": "HPC-1.0",
  "workflowType": "MaxQuantWorkflow",
  "hpcConfig": {
    "threads": 32,
    "memoryGB": 128,
    "wallTime": "24:00:00"
  },
  "inputFiles": ["/data/sample1.mzML"],
  "outputDirectory": "/results/",
  "searchParameters": {
    "enzyme": "Trypsin",
    "missedCleavages": 2,
    "fixedMods": ["Carbamidomethyl (C)"],
    "variableMods": ["Oxidation (M)"]
  }
}
```

### Linux Path Handling
```java
// Windows → Linux path translation
public class PathTranslator {
    public static String translatePath(String windowsPath) {
        // C:\data\sample.raw → /home/user/data/sample.mzML
        return windowsPath
            .replace("\\", "/")
            .replace("C:/", "/home/user/")
            .replace(".raw", ".mzML");
    }
}
```

---

## 📊 Testing Strategy

### Unit Testing
```java
@Test
public void testHeadlessExecution() {
    MetaPilotHeadless headless = new MetaPilotHeadless();
    String[] args = {"--par", "test.json", "--output", "/tmp/"};
    
    // Should execute without GUI dependencies
    assertDoesNotThrow(() -> headless.main(args));
}
```

### Integration Testing
```bash
#!/bin/bash
# Test Linux compatibility
echo "Testing headless execution..."
java -jar metapilot-headless.jar --validate test_params.json

echo "Testing small dataset..."
java -jar metapilot-headless.jar \
  --par small_test.json \
  --output /tmp/test/ \
  --threads 2
```

### HPC Testing
```bash
# SLURM test job
sbatch --job-name=test-metapilot \
       --time=01:00:00 \
       --mem=8GB \
       test_metapilot.sh
```

---

## 🔧 File Format Migration

### Current Challenge: Thermo .raw Files
```
Windows: sample.raw → [Vendor Tools] → Analysis
Linux:   sample.raw → [Convert to mzML] → Transfer → Analysis
```

### Solution: Desktop Conversion
```java
// On Windows desktop
public class FileFormatConverter {
    public void convertRawToMzML(File rawFile, File outputDir) {
        // Use existing msconvert integration
        MSConvertor convertor = new MSConvertor();
        convertor.convertToMzML(rawFile, outputDir);
    }
}
```

### Workflow:
1. **Windows**: Convert .raw → .mzML
2. **Transfer**: Upload .mzML to Linux
3. **Linux**: Process .mzML files
4. **Results**: Download analysis results

---

## 🚨 Critical Success Factors

### Technical
- ✅ **Zero JavaFX in headless JAR** - Critical for server deployment
- ✅ **Identical results** - Windows vs Linux output must match
- ✅ **Linear scaling** - Performance must scale with cores
- ✅ **Robust error handling** - Graceful failure in HPC environment

### User Experience
- ✅ **Familiar workflow** - Minimal disruption to existing users
- ✅ **Clear documentation** - Step-by-step migration guides
- ✅ **Support tools** - Upload/download assistants
- ✅ **Fallback options** - Desktop version remains available

### Operational
- ✅ **HPC compatibility** - Works with major schedulers
- ✅ **Container support** - Docker/Singularity ready
- ✅ **Security compliance** - Meets HPC security requirements
- ✅ **Performance monitoring** - Comprehensive logging and metrics

---

## 📖 Key Resources

### Documentation
- [Complete Technical Report](MetaPilot_HPC_Migration_Plan.md) - Comprehensive implementation guide
- [Linux File System Guide](MetaPilot_HPC_Migration_Plan.md#linux-file-system-guide) - Navigation basics
- [Parameter File Schema](MetaPilot_HPC_Migration_Plan.md#parameter-file-schema) - JSON structure

### Code Examples
- Current CLI implementations: `MetaLabCli.java`, `MetaPilotTask.java`
- Parameter handlers: `JsonParaHandler.java`, `MetaParameterIOHandler.java`
- Build configurations: `pom.xml`

### External Tools
- **ProteoWizard**: Cross-platform file conversion
- **OpenMS**: Linux-native proteomics tools
- **ThermoRawFileParser**: Open-source .raw file reader

---

## 🎯 Next Steps

### Immediate Actions (Week 1)
1. **Environment Setup**
   ```bash
   # Set up Linux development environment
   # Install Java 11+, Maven, Git
   git clone <your-repo>
   cd MetaPilot
   mvn clean compile
   ```

2. **Current State Analysis**
   ```bash
   # Test existing CLI
   java -jar metapilot.jar --help
   java -jar metapilot.jar --par sample_parameters.json
   ```

3. **Plan Validation**
   - Review technical requirements with HPC team
   - Identify target HPC systems
   - Confirm file format requirements

### Week 2-4: Foundation
- Create headless main class
- Set up build profiles
- Extract core computational logic
- Basic CLI testing

### Month 2 and Beyond
- Follow the detailed phase-by-phase plan
- Regular testing on Linux environments
- User feedback collection
- Performance optimization

---

**Remember**: This is a gradual migration, not a complete rewrite. The existing codebase provides a solid foundation - we're adapting it for new environments while preserving functionality.

**Success Metric**: Users can configure analyses on Windows desktop and execute them on Linux HPC clusters with identical results and 10x improved performance.

---

*Quick Reference Version: 1.0*  
*For detailed technical specifications, see: [MetaPilot_HPC_Migration_Plan.md](MetaPilot_HPC_Migration_Plan.md)* 
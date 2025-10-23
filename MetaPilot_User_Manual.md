# MetaPilot User Manual

## Table of Contents
1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Installation and Setup](#installation-and-setup)
4. [Resource Configuration](#resource-configuration)
5. [MAG Database Setup](#mag-database-setup)
6. [Workflow Configuration](#workflow-configuration)
7. [Running Analyses](#running-analyses)
8. [Troubleshooting](#troubleshooting)
9. [Advanced Usage](#advanced-usage)

## Introduction

MetaPilot is a comprehensive metaproteomics data analysis platform that unifies multiple search engines and workflows behind a modern JavaFX GUI and a headless CLI for HPC environments. It supports both Data-Dependent Acquisition (DDA) and Data-Independent Acquisition (DIA) workflows for metaproteomics analysis.

### Key Features
- **Multi-workflow support**: DDA (IGC and MAG) and DIA workflows
- **Multi-engine integration**: MaxQuant, pFind, FragPipe (MSFragger), Sage, DIA-NN
- **Resource management**: Centralized tool/database configuration
- **MAG database management**: Configure MGnify/local MAG catalogs
- **GUI and CLI**: JavaFX desktop app and `--cli` mode for HPC/headless runs

## System Requirements

### Software Requirements
- **Java**: Version 17 or higher
- **Maven**: Version 3.8 or higher (for building from source)
- **JavaFX**: Version 21 runtime (for GUI mode)
- **Operating System**: Windows, Linux, or macOS

### Hardware Requirements
- **RAM**: Minimum 8GB, recommended 16GB or more
- **Storage**: At least 10GB free space for tools and databases
- **CPU**: Multi-core processor recommended for parallel processing

### External Tools
MetaPilot requires several external tools to be installed and configured:
- MaxQuant (for DDA IGC workflows)
- pFind (for DDA workflows)
- FragPipe (for DDA workflows)
- Sage (for DDA workflows)
- DIA-NN (for DIA workflows)
- msconvert (ProteoWizard)
- FlashLFQ (for quantification)
- Python (for scripts)
- FunctionDB (for functional annotation)

## Installation and Setup

### 1. Download and Build

```bash
# Clone the repository (if building from source)
git clone <repository-url>
cd MetaPilot

# Build the application
mvn -q -DskipTests package
```

The shaded JAR will be created in `target/` with the main class `uk.ac.quadram.metapilot.MetaPilotApp`.

### 2. Run MetaPilot

#### GUI Mode (Recommended)
```bash
# Using Maven JavaFX plugin
mvn -q javafx:run

# Or run the JAR directly (requires JavaFX modules)
java --module-path "$PATH_TO_FX" \
     --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base \
     -jar target/*-shaded.jar
```

#### CLI Mode (for HPC)
```bash
java -jar target/*-shaded.jar --cli \
  --params config/last_dda_igc_maxquant_parameters.json \
  --resources config/meta_resources.params
```

## Resource Configuration

### Overview
MetaPilot uses a centralized resource management system to handle all external tool paths and dependencies. Configuration is stored in `config/meta_resources.params` and can be managed through the GUI or manually.

### Setting Up Resources

#### 1. Base Resource Directory
The base resource directory contains all external tools. Default structure:
```
Resources/
├── mq_bin/MaxQuantCmd.exe           (MaxQuant)
├── DBReducer/pFind3/bin/pFind.exe   (pFind)
├── fragpipe/bin/fragpipe.bat        (FragPipe)
├── sage/sage.exe                    (Sage)
├── DIA-NN/2.0/DiaNN.exe            (DIA-NN)
├── ProteoWizard/msconvert.exe       (msconvert)
├── FlashLFQ/CMD.exe                 (FlashLFQ)
├── function/func_def.db             (FunctionDB)
└── Python/python.exe                (Python)
```

#### 2. Configure via GUI
1. Launch MetaPilot
2. Go to the **Resources** tab
3. Set the **Base Resource Path** to your tools directory
4. Configure individual tool paths as needed
5. Click **Save** to persist settings

#### 3. Manual Configuration
Edit `config/meta_resources.params`:
```properties
# MetaPilot Resource Configuration
BASE_PATH=C:\YourPath\Resources

# Tool Paths
PYTHON=C:\YourPath\Resources\python-3.13.7-embed-amd64\python.exe
FRAGPIPE=C:\YourPath\Resources\fragpipe\FragPipe-23.1\bin\fragpipe.bat
FLASHLFQ=C:\YourPath\Resources\FlashLFQ\CMD.exe
FUNCDB=C:\YourPath\Resources\function\func_def.db
MAXQUANT=C:\YourPath\Resources\mq_bin\MaxQuantCmd.exe
SAGE=C:\YourPath\Resources\sage-v0.14.7-x86_64-pc-windows-msvc\sage.exe
DEEPDETECT=C:\YourPath\Resources\DeepDetect\DeepDetect.exe
DIANN=C:\YourPath\Resources\DiaNN\diann.exe
PFIND=C:\YourPath\Resources\DBReducer\pFind3\bin\pFind.exe
DBREDUCER=C:\YourPath\Resources\DBReducer\DBReducer\DBReducer.exe
MSCONVERT=C:\YourPath\Resources\ProteoWizard\msconvert.exe
```

### Tool Download Links
- **MaxQuant**: https://www.maxquant.org/maxquant/
- **FragPipe**: https://fragpipe.nesvilab.org/
- **Sage**: https://github.com/lazear/sage
- **DIA-NN**: https://github.com/vdemichev/DiaNN/releases
- **msconvert**: https://proteowizard.sourceforge.io/download.html

## MAG Database Setup

### Overview
MAG (Metagenome-Assembled Genomes) databases are configured through the `MAGDB_1_0.json` file and managed via the Database Configuration tab in MetaPilot.

### Database Types

#### 1. MGnify Databases
Pre-configured databases from the MGnify resource:
- **Human Gut**: `human-gut` (v2.0.2)
- **Mouse Gut**: `mouse-gut` (v2.0.2)
- **Human Skin**: `human-skin` (v2.0.2)
- **Pig Gut**: `pig-gut` (v2.0.2)

#### 2. Local Databases
Custom MAG databases that you provide:
- Local protein sequence files
- Custom taxonomic annotations
- User-defined functional annotations

### Setting Up MAG Databases

#### 1. Via GUI
1. Launch MetaPilot
2. Go to the **Database Configuration** tab
3. Click **Add Database** for new entries
4. Configure database settings:
   - **Database Name**: Descriptive name
   - **Database Type**: MGnify or Local
   - **Version**: Database version
   - **Local Path**: Path to database files
   - **Download Status**: Available/Downloaded
5. Click **Save** to persist configuration

#### 2. Manual Configuration
Edit `MAGDB_1_0.json` (this file is large, so use the GUI for most changes):

```json
{
  "mgnify_databases": [
    {
      "name": "Human Gut Microbiome",
      "identifier": "MGYG",
      "catalogue_id": "human-gut",
      "catalogue_version": "v2.0.2",
      "versions": [
        {
          "version": "v2.0.2",
          "local_path": "C:/path/to/human-gut/v2.0.2",
          "downloaded": true,
          "size_gb": 15.2
        }
      ]
    }
  ],
  "local_databases": [
    {
      "name": "Custom Environmental Database",
      "identifier": "CUSTOM",
      "versions": [
        {
          "version": "v1.0",
          "local_path": "C:/path/to/custom/database",
          "downloaded": true,
          "size_gb": 8.5
        }
      ]
    }
  ]
}
```

### Database Download
1. Select a database in the Database Configuration tab
2. Click **Download** to fetch the database files
3. Monitor download progress in the console
4. Verify local availability status

## Workflow Configuration

### Supported Workflows

#### DDA (Data-Dependent Acquisition) Workflows

##### IGC Database Workflows
- **MaxQuant DDA IGC**: Label-free quantification using MaxQuant
- **pFind DDA IGC**: Chinese-developed search engine
- **FragPipe DDA IGC**: MSFragger-based pipeline
- **Sage DDA IGC**: Fast Rust-based search engine

##### MAG Database Workflows
- **pFind MAG**: Targeted MAG database searches
- **FragPipe MAG**: MAG-specific FragPipe analysis
- **Sage MAG**: High-speed MAG proteome analysis

#### DIA (Data-Independent Acquisition) Workflows
- **DIA Normal**: Standard DIA analysis with spectral libraries
- **DIA Target**: Targeted DIA analysis for specific proteins
- **DIA Real-time**: Real-time DIA analysis and monitoring

### Workflow Parameters

#### Common Parameters (All Workflows)

##### Input/Output Settings
- **Input Files**: Raw mass spectrometry files (.raw, .mzML, .d)
- **Output Directory**: Results output location
- **Database**: Protein sequence database (IGC or MAG)
- **Host Database**: Optional host protein database
- **Append Host DB**: Whether to include host proteins

##### Sample Metadata
- **Sample Name**: Unique identifier for each sample
- **Condition**: Experimental condition
- **Treatment**: Treatment group
- **Time Point**: Time point in study
- **Age**: Subject age
- **Gender**: Subject gender
- **Custom Fields**: Additional metadata fields

##### Enzyme Settings
- **Enzyme**: Proteolytic enzyme (Trypsin, Lys-C, etc.)
- **Missed Cleavages**: Maximum missed cleavages (0-3)
- **Digestion Mode**: Specific or non-specific

##### Modifications
- **Fixed Modifications**: Always present modifications
  - Carbamidomethyl (C) - standard for cysteine
- **Variable Modifications**: Optional modifications
  - Acetyl (Protein N-term)
  - Oxidation (M)
  - Phosphorylation (S, T, Y)

##### Taxonomic Analysis
- **Bacteria**: Include bacterial proteins
- **Archaea**: Include archaeal proteins
- **Eukaryota**: Include eukaryotic proteins
- **Viruses**: Include viral proteins
- **Unclassified**: Include unclassified sequences
- **Ignore Environmental Samples**: Filter out environmental contaminants
- **Ignore Rank**: Taxonomic rank to ignore (Family, Genus, etc.)

#### DDA-Specific Parameters

##### Search Parameters
- **MS2 Scan Mode**: FTMS or ITMS
- **Threads**: Number of processing threads
- **Generate Sample-Specific DB**: Create custom databases per sample
- **Spectra Clustering**: Enable/disable spectral clustering

##### Quantification
- **Quantitative Mode**: Label-free, Isotopic-labeling, or Isobaric-labeling
- **Label Type**: For labeling experiments
- **Reference Table**: For isobaric labeling

#### DIA-Specific Parameters

##### Basic Parameters
- **MS2 Scan Mode**: DIA
- **Threads**: Number of processing threads
- **Meta Workflow**: Enable metaproteomics-specific processing
- **Core Count**: Number of cores to use

##### DIA Settings
- **Library Search**: Use spectral libraries
- **Self Model Search**: Use self-generated models
- **Models**: Custom model files
- **Libraries**: Spectral library files

### Parameter Configuration Examples

#### DDA IGC MaxQuant Example
```json
{
  "workflow": "DDA_IGC_MAXQUANT",
  "input_output": {
    "host_database": "C:/database/uniprotkb_Human_AND_reviewed_true_AND_m_2025_09_01.fasta",
    "microbiome_database": "C:/database/760MetaHit_139HMP_368PKU_511Bac.uniq.pep.fasta",
    "output_directory": "C:/results/analysis",
    "append_host_db": true,
    "experiments": [
      {
        "file_path": "C:/data/sample1.raw",
        "experiment_name": "Sample1"
      }
    ]
  },
  "enzyme_settings": {
    "enzyme": "Trypsin",
    "missed_cleavages": 2,
    "digestion_mode": "Specific"
  },
  "modifications": {
    "fixed_modifications": ["Carbamidomethyl (C)"],
    "variable_modifications": [
      "Acetyl (Protein N-term)",
      "Oxidation (M)"
    ]
  },
  "taxonomic_analysis": {
    "bacteria": true,
    "archaea": true,
    "eukaryota": true,
    "viruses": false,
    "unclassified": false,
    "ignore_environmental_samples": true,
    "ignore_rank": "Family"
  },
  "additional_parameters": {
    "generate_sample_specific_db": true,
    "ms2_scan_mode": "FTMS",
    "spectra_clustering": false,
    "threads": 16
  }
}
```

#### DIA Normal Example
```json
{
  "workflow": "DiaNNMAG",
  "mag_database": {
    "catalogue_version": "v2.0.2",
    "identifier": "MGYG",
    "catalogue_id": "human-gut"
  },
  "input_output": {
    "micro_db": "C:/database/MAG/human-gut/v2.0.2",
    "output_directory": "C:/results/dia_analysis",
    "experiments": [
      {
        "file_path": "C:/data/dia_sample1.raw",
        "experiment_name": "DIA_Sample1"
      }
    ]
  },
  "basic_parameters": {
    "ms2_scan_mode": "DIA",
    "threads": 20,
    "meta_workflow": true,
    "core_count": 1
  },
  "dia_settings": {
    "library_search": true,
    "self_model_search": false,
    "models": [],
    "use_model": false,
    "libraries": []
  }
}
```

## Running Analyses

### 1. Launch Workflow
1. Start MetaPilot
2. Go to the **Main** tab
3. Select your desired workflow from the workflow panels
4. Click **Launch** to open the workflow configuration window

### 2. Configure Parameters
1. **Input/Output Tab**:
   - Set input files and output directory
   - Configure database paths
   - Add sample metadata

2. **Parameters Tab**:
   - Configure enzyme settings
   - Set modifications
   - Adjust search parameters
   - Configure taxonomic analysis

3. **Advanced Tab** (if available):
   - Set thread counts
   - Configure additional options
   - Set quantitative parameters

### 3. Start Analysis
1. Click **Start Analysis** in the console panel
2. Monitor progress in the console output
3. Analysis will run in the background
4. Results will be saved to the specified output directory

### 4. View Results
1. Click **View Report** when analysis completes
2. Browse generated reports and visualizations
3. Export results as needed

### CLI Usage
For HPC environments, use the CLI mode:

```bash
# Basic CLI usage
java -jar metapilot.jar --cli \
  --params config/my_parameters.json \
  --resources config/meta_resources.params

# CLI help
java -jar metapilot.jar --cli --help
```

## Troubleshooting

### Common Issues

#### 1. JavaFX Runtime Errors
**Problem**: "JavaFX runtime components are missing"
**Solution**: 
- Install JavaFX SDK
- Use `--module-path` and `--add-modules` flags
- Or use Maven JavaFX plugin: `mvn javafx:run`

#### 2. Missing Tools
**Problem**: Workflow launch disabled
**Solution**:
1. Go to **Resources** tab
2. Configure missing tool paths
3. Verify tool availability
4. Save configuration

#### 3. MAG Database Not Loading
**Problem**: Database configuration errors
**Solution**:
1. Check `MAGDB_1_0.json` file exists and is readable
2. Verify database paths are correct
3. Re-download databases if needed
4. Use GUI to reconfigure databases

#### 4. Memory Issues
**Problem**: Out of memory errors
**Solution**:
- Increase JVM heap size: `-Xmx16g`
- Reduce thread count in parameters
- Use smaller databases if possible

#### 5. File Path Issues
**Problem**: File not found errors
**Solution**:
- Use absolute paths
- Check file permissions
- Verify file formats are supported

### Log Files
- Application logs: `logs/metapilot.log`
- Workflow-specific logs in output directories
- Console output in MetaPilot GUI

### Performance Optimization
1. **Thread Configuration**: Set appropriate thread counts based on your system
2. **Database Size**: Use appropriate database sizes for your samples
3. **Memory Allocation**: Allocate sufficient memory for large datasets
4. **Storage**: Use fast storage (SSD) for temporary files

## Advanced Usage

### Parameter File Management
- **Save Parameters**: Use "Save Parameters" to store current settings
- **Load Parameters**: Use "Load Parameters" to restore previous settings
- **Parameter Files**: Stored in `config/` directory
- **Last Used**: Automatically saved as `last_used_*_parameters.json`

### Batch Processing
1. Create parameter files for multiple analyses
2. Use CLI mode for batch processing
3. Script the execution for automation
4. Monitor progress across multiple jobs

### Custom Databases
1. Prepare protein sequence files in FASTA format
2. Add taxonomic annotations if needed
3. Configure in Database Configuration tab
4. Test with small datasets first

### Integration with HPC
1. Use CLI mode for headless execution
2. Submit jobs via job schedulers (SLURM, PBS, etc.)
3. Configure resource requirements
4. Monitor job status and results

### Quality Control
1. Check input file quality
2. Verify database completeness
3. Monitor analysis progress
4. Validate results before interpretation

---

## Support and Documentation

For additional support:
- Check the `docs/` directory for detailed guides
- Review configuration examples in `config/`
- Consult the troubleshooting section above
- Contact the development team for advanced issues

## Version Information

This manual covers MetaPilot v2.0. For the most up-to-date information, refer to the application's built-in help and documentation.

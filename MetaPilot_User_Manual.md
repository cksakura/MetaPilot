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

### DDA Workflows (Standard Proteomics)

#### Quick Start for DDA Workflows
1. **Launch Workflow**: Go to **Main** tab → Select DDA workflow (MaxQuant, pFind, FragPipe, or Sage) → Click **Launch**
2. **Input/Output Tab**:
   - Add raw files (.raw, .mzML, .d)
   - Set output directory
   - Select IGC database (for IGC workflows) or MAG database (for MAG workflows)
   - Optionally add host database
3. **Parameters Tab**:
   - Configure enzyme (default: Trypsin)
   - Set missed cleavages (default: 2)
   - Select modifications
   - Configure taxonomic analysis settings
4. **Console Tab**: Click **Start Analysis**
5. **Monitor Progress**: Watch console output and progress bars
6. **View Results**: Click **View Report** when complete

### DIA Workflows (Advanced Metaproteomics)

DIA workflows are more complex and offer three different analysis strategies. Understanding these strategies is crucial for successful DIA analysis.

#### Prerequisites for DIA Workflows

1. **MAG Database Setup**:
   - Go to **Database Configuration** tab
   - Select a MAG database (e.g., human-gut v2.0.2)
   - Ensure database is downloaded and available locally

2. **Predicted Peptides Generation** (CRITICAL STEP):
   - After selecting MAG database, go to **Input/Output** tab
   - Check the **MAG Library Status** section
   - Click **Check Status** to verify predicted peptide files
   - If not complete, click **Start Prediction** to generate predicted peptides
   - **This step is required before any DIA analysis** - it pre-computes peptides for all genomes in the database
   - Monitor progress in the status console
   - Wait for completion (this can take several hours for large databases)

#### Understanding DIA Analysis Strategies

##### Strategy 1: Library Analysis (Direct Identification)
**Best for**: Samples similar to existing reference libraries

**Process**:
- Uses pre-built spectral libraries for direct identification
- Fastest approach if suitable libraries exist
- No iterative search required
- Recommended for routine analyses with established protocols

**Steps**:
1. Launch **DIA Normal** workflow
2. **Input/Output Tab**:
   - Add DIA raw files
   - Set output directory
   - Select MAG database
   - Verify predicted peptides are available
3. **Parameters Tab**:
   - Set **Analysis Mode** to "Directly identification (no iterative search)"
   - Set thread count (recommended: 50-70% of available cores)
   - **Library Generation Section**:
     - Click **Check Available** to see existing libraries
     - Select libraries to use (check the "Use" checkbox)
     - OR create new library:
       1. Select input type (Existing Library, Peptide List, or Protein Database)
       2. Browse for input file
       3. Enter library name
       4. Click **Generate Library**
       5. Monitor progress in console
       6. Once complete, select the new library for use
4. **Console Tab**: Click **Start Analysis**
5. Results will be in `output_directory/library_analysis/`

##### Strategy 2: Iterative Analysis (Sample-Specific Library Creation)
**Best for**: Novel samples without suitable reference libraries

**Process**:
- Creates sample-specific libraries through 4 rounds of iterative search:
  1. **Round 1**: Initial search using selected broad libraries
  2. **Round 2**: Refine database based on Round 1 identifications
  3. **Round 3**: Further refinement based on Round 2 results
  4. **Round 4**: Final search with optimized sample-specific library
- Each round narrows down the MAG genomes and improves identification
- Takes longer but provides better coverage for complex samples

**Steps**:
1. Launch **DIA Normal** workflow
2. **Input/Output Tab**: Configure as above
3. **Parameters Tab**:
   - Set **Analysis Mode** to "Iterative search"
   - Set thread count (recommended: 50-70% of available cores)
   - **Library Generation Section**:
     - Select initial seed libraries (these guide the first round)
     - Broad, comprehensive libraries work best as seeds
     - Typically 2-5 libraries covering expected taxonomic groups
   - **ML Model Section** (Optional):
     - Can optionally select ML models to enhance iterative search
     - If models are available and selected, they will be used to predict peptide detectability
     - Enhances genome selection in each iterative round
4. **Console Tab**: Click **Start Analysis**
5. Monitor the 4-round iterative process:
   - Round 1: ~15-20% progress
   - Round 2: ~35-50% progress
   - Round 3: ~65-80% progress
   - Round 4: ~80-85% progress
   - Final identification: ~85-100% progress
6. Results will be in `output_directory/iterative_analysis/`
   - Separate folders for each sample with round-by-round results
   - Final combined results in `combined/` folder

##### Strategy 3: Self-Model Analysis (Deep Learning Enhanced)
**Best for**: Complex samples requiring maximum sensitivity

**Process**:
- Combines iterative search with self-generated machine learning models
- **Round 1**: Initial library search to identify peptides
- **Model Building**: Trains ML models on Round 1 results to predict peptide detectability
- **Rounds 2-4**: Iterative search enhanced by ML predictions
- Most computationally intensive but highest sensitivity

**Steps**:
1. **Preparation** - Train initial ML models (if not already available):
   - Go to **Parameters Tab** → **ML Model Generation** section
   - Browse for peptide intensity table (from previous DIA or DDA analysis)
   - Enter model name
   - Click **Generate Model**
   - Wait for training to complete (monitor progress bar)
   - Models are saved to MAG database folder

2. **Analysis Setup**:
   - Launch **DIA Normal** workflow
   - **Input/Output Tab**: Configure as above
   - **Parameters Tab**:
     - Set **Analysis Mode** to "Deep learning mode"
     - Set thread count
     - **Library Generation**: Select seed libraries
     - **ML Model Section**:
       - Click **Check Available** to see trained models
       - Select models to use (check the "Use" checkbox)
       - Multiple models can be selected for ensemble prediction

3. **Console Tab**: Click **Start Analysis**

4. **Process Flow**:
   - Initial library search (Round 1)
   - **Automatic self-model generation** from Round 1 results
   - ML-enhanced iterative search (Rounds 2-4)
   - Final identification and quantification

5. Results will be in `output_directory/self_model_analysis/`
   - Individual sample folders with model subdirectories
   - Models trained specifically for each sample
   - Final combined results with ML-enhanced identifications

#### DIA Workflow Comparison

| Feature | Library Analysis | Iterative Analysis | Self-Model Analysis |
|---------|-----------------|-------------------|---------------------|
| **Speed** | Fastest (~2-4 hours) | Medium (~6-12 hours) | Slowest (~12-24 hours) |
| **Sensitivity** | Good | Better | Best |
| **Requirements** | Suitable libraries | Seed libraries | Seed libraries + intensity data or auto-generated |
| **Best For** | Routine samples | Novel samples | Complex/challenging samples |
| **Disk Space** | Low | Medium | High (stores models) |
| **Technical Skill** | Beginner | Intermediate | Advanced |

#### DIA Target Workflow (Targeted Analysis)
**Purpose**: Focus on specific taxa or functions of interest

**Steps**:
1. Launch **DIA Target** workflow
2. **Input/Output Tab**: Add DIA files, set output directory
3. **Parameters Tab**:
   - **Taxonomic Selection**:
     - Navigate taxonomic tree (Phylum → Class → Order → Family → Genus → Genome)
     - Select specific taxa of interest
     - System will filter to proteins from selected taxa
   - **Functional Selection** (Optional):
     - Choose function type (COG, KEGG, GO, EC)
     - Browse and select specific functions
     - Further narrows protein targets
   - **Library Name**: Enter name for the targeted library
   - Click **Generate Library** to create targeted spectral library
4. **Console Tab**: Click **Start Analysis**
5. Results focused on selected taxa/functions

#### DIA Real-time Workflow (Live Monitoring)
**Purpose**: Real-time monitoring of DIA acquisitions as files are generated

**Steps**:
1. Launch **DIA Real-time** workflow
2. **Input/Output Tab**:
   - Set **Monitor Folder**: Folder where instrument writes .raw files
   - Set **Result Folder**: Where results will be saved
   - Select MAG database
3. **Parameters Tab**:
   - Select spectral libraries for identification
   - Configure thread count
4. **Console Tab**: Click **Start Monitoring**
5. **Live Process**:
   - System watches monitor folder for new .raw files
   - Automatically processes each file as it completes
   - Real-time identification and quantification
   - Results update continuously

### Detailed Guide: Building Spectral Libraries for DIA

Spectral libraries are essential for DIA analysis. MetaPilot provides flexible library generation options.

#### Option 1: From Existing Spectral Library File
**Use when**: You have a pre-built .speclib file from DIA-NN or other tools

1. In **Parameters Tab** → **Spectral Library Generation** section
2. Select **Input Type**: "Existing Library File"
3. Click **Browse** and select your .speclib file
4. Enter a **Library Name** (e.g., "gut_microbiome_lib")
5. Click **Generate Library**
6. Library will be validated and added to available libraries
7. Check the "Use" checkbox to include it in analysis

#### Option 2: From Peptide List
**Use when**: You have a list of peptides from previous analyses

1. Prepare a peptide list file (.txt, .tsv, .csv):
   ```
   PEPTIDEK
   PROTEINR
   SEQUENCEK
   ```
2. In **Spectral Library Generation** section:
   - Select **Input Type**: "Peptide List File"
   - Click **Browse** and select your peptide list file
   - Enter **Library Name**
   - Click **Generate Library**
3. DIA-NN will predict spectra for these peptides
4. Monitor progress in the generation console
5. Once complete, check "Use" to include in analysis

#### Option 3: From Protein Database (FASTA)
**Use when**: You have a FASTA protein database and want to predict peptides

1. Prepare or select a protein database (.fasta file)
2. In **Spectral Library Generation** section:
   - Select **Input Type**: "Protein Database File"
   - Click **Browse** and select your .fasta file
   - Configure **Enzyme Parameters**:
     - **Enzyme**: Trypsin/P, Trypsin, Lys-C, etc.
     - **Max Miss Cleavages**: 0, 1, or 2
     - **Keep protein names**: Check if proteins have meaningful annotations
   - Enter **Library Name**
   - Click **Generate Library**
3. System will:
   - In-silico digest proteins with selected enzyme
   - Predict spectra for resulting peptides
   - Create spectral library
4. Monitor progress (may take 10-60 minutes depending on database size)
5. Check "Use" to include in analysis

#### Library Management Tips
- **Multiple Libraries**: You can select multiple libraries for comprehensive coverage
- **Library Location**: Libraries are stored in `MAG_database_path/libraries/library_name/`
- **Check Available**: Click this button to refresh the list from the database folder
- **Delete**: Remove unwanted libraries (will delete files from disk)
- **Reuse**: Libraries persist across sessions and can be used for multiple analyses

### Detailed Guide: Building ML Models for DIA

ML models predict peptide detectability to enhance identification in iterative searches. **Note: Models are optional** - iterative search works without them.

#### When to Use ML Models
- ✓ You have prior DIA or DDA data from similar sample types
- ✓ You want maximum sensitivity for challenging samples
- ✓ You're using "Iterative search" or "Deep learning mode"
- ✗ Don't use if you have no training data
- ✗ Not needed for simple library-based analysis

#### Building an ML Model

**Prerequisites**:
- Peptide intensity table from previous analysis (DIA-NN output or similar)
- File should contain columns: Protein, Peptide, Intensity values
- Python configured in Resources (for PyTorch model training)

**Steps**:

1. **Prepare Training Data**:
   - Export peptide intensity data from a previous successful analysis
   - File format: TSV/CSV with peptide sequences and intensities
   - Example from DIA-NN: `report.pg_matrix.tsv` or `report.pr_matrix.tsv`

2. **Generate Model** in MetaPilot:
   - Go to **Parameters Tab** → **ML Model Generation** section
   - Click **Browse** next to "Peptide Intensity Table"
   - Select your intensity file
   - Enter a descriptive **Model Name** (e.g., "gut_microbiome_model_2024")
   - Click **Generate Model**

3. **Training Process**:
   - System extracts peptide-intensity relationships
   - Trains neural network model (LSTM + attention mechanism)
   - Training typically takes 30-120 minutes depending on data size
   - Monitor progress in the model generation console
   - Progress bar shows training epochs

4. **Model Validation**:
   - Upon completion, model statistics are displayed:
     - **Train Size**: Number of peptides used for training
     - **Model Size**: File size of trained model
     - **Correlation**: Prediction accuracy (higher is better, >0.7 is good)
     - **Median**: Median prediction value
   - Model files saved to `MAG_database_path/models/model_name/`

5. **Using the Model**:
   - Click **Check Available** to refresh model list
   - Check the "Use" checkbox for models you want to apply
   - Multiple models can be selected for ensemble prediction
   - Selected models will enhance peptide prediction in iterative rounds

#### Understanding Model-Enhanced Iterative Search

When models are selected:
- **Round 1**: Uses traditional predicted peptides based on genome presence
- **Round 2-4**: ML models predict which peptides are most detectable
  - Models score each potential peptide based on learned detectability patterns
  - High-scoring peptides are prioritized in database generation
  - Improves genome selection accuracy
  - Reduces false positive identifications

#### Self-Model Analysis (Automatic Model Generation)

In "Deep learning mode":
1. **Round 1** proceeds with library search (no models yet)
2. **Automatic Model Training**: After Round 1, system automatically:
   - Extracts peptide intensities from Round 1 results
   - Trains a sample-specific ML model
   - No manual intervention required
3. **Rounds 2-4**: Uses the self-generated model for enhanced prediction
4. **Result**: Best possible identification for that specific sample

**Advantages**:
- Sample-specific models capture unique characteristics
- No need for prior training data
- Adapts to your exact experimental conditions

**Disadvantages**:
- Requires successful Round 1 (enough identifications for training)
- Longest analysis time
- Highest computational requirements

#### Model Management
- **Storage**: Models are stored in `MAG_database_path/models/`
- **Reuse**: Models can be used across different analyses
- **Sharing**: Model folders can be copied between MAG databases
- **Deletion**: Remove unwanted models to save space

### DIA Workflow Decision Tree

```
START: Do you have DIA data to analyze?
  │
  ├─ YES → Do you have suitable pre-built spectral libraries?
  │         │
  │         ├─ YES → Use LIBRARY ANALYSIS (Strategy 1)
  │         │        - Fastest approach
  │         │        - Direct identification
  │         │
  │         └─ NO → Is this a novel/complex sample?
  │                  │
  │                  ├─ NO (routine sample) → Use ITERATIVE ANALYSIS (Strategy 2)
  │                  │                         - Sample-specific library creation
  │                  │                         - 4-round refinement
  │                  │
  │                  └─ YES (complex/challenging) → Use SELF-MODEL ANALYSIS (Strategy 3)
  │                                                  - ML-enhanced identification
  │                                                  - Maximum sensitivity
  │
  └─ NO → Start with library/model generation for future analyses
```

### DIA Workflows: Common Questions

#### Q: Which DIA strategy should I use?
**A**: Follow this guideline:
- **New to DIA or routine samples**: Start with Library Analysis
- **No suitable libraries exist**: Use Iterative Analysis
- **Maximum sensitivity needed**: Use Self-Model Analysis
- **Specific taxa/pathways of interest**: Use DIA Target workflow

#### Q: How long does predicted peptide generation take?
**A**: Depends on database size:
- Small database (~1,000 genomes): 2-4 hours
- Medium database (~5,000 genomes): 8-16 hours  
- Large database (~10,000+ genomes): 24-48 hours
- **Important**: This is a one-time process per database version

#### Q: Can I stop and resume predicted peptide generation?
**A**: Yes! The DeepDetect task tracks progress:
- If stopped, resume by clicking "Start Prediction" again
- Already-generated files are skipped
- Only remaining genomes will be processed

#### Q: How do I know if my libraries are good quality?
**A**: Check library generation console output:
- Look for "SUCCESS" messages
- Check file sizes (libraries should be >1MB typically)
- Test with a small dataset first
- Good libraries should give >1000 identifications for complex samples

#### Q: What if iterative search doesn't find enough peptides in Round 1?
**A**: This usually means:
- Libraries don't match your sample type → Try different seed libraries
- MAG database doesn't contain your sample's microbiome → Use broader database
- Technical issues with raw files → Check file quality and format

#### Q: Do I need to generate models for every analysis?
**A**: No! Models are reusable:
- Generate once for a sample type (e.g., human gut samples)
- Use the same model for all similar samples
- Only regenerate if changing to very different sample types

#### Q: Can I use libraries and models from different MAG databases?
**A**: Libraries and models are specific to their MAG database:
- Libraries must be generated from the same MAG database version
- Models can work across databases if trained on similar data
- For best results, keep libraries and models with their source database

#### Q: What's the difference between "library search" and "iterative search"?
**A**: 
- **Library search** (Strategy 1): Direct matching against pre-built libraries - fast but limited to known peptides
- **Iterative search** (Strategy 2): Creates custom libraries from your samples through multiple refinement rounds - slower but more comprehensive

#### Q: How much disk space do I need for DIA analysis?
**A**: Estimate for a typical analysis:
- Predicted peptides: 50-200 GB (one-time, per MAG database)
- Library analysis: 5-20 GB per analysis
- Iterative analysis: 20-100 GB per analysis
- Self-model analysis: 50-200 GB per analysis (includes models)
- **Recommendation**: 500 GB free space for comfortable operation

### General Analysis Tips

#### Before Starting Analysis:
1. ✓ Verify all tools are configured in Resources tab
2. ✓ Check MAG database is downloaded and available
3. ✓ For DIA: Ensure predicted peptides are generated
4. ✓ For DIA: Select appropriate libraries/models
5. ✓ Verify sufficient disk space for results
6. ✓ Set appropriate thread count (don't use all cores)

#### During Analysis:
- Monitor console output for errors
- Check progress bars for workflow stages
- Don't close the application
- Ensure system doesn't go to sleep

#### After Analysis:
- Check result directory for output files
- Review generated reports
- Save parameter file for reproducibility
- Archive results before next analysis

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

### Common Issues - General

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
- Increase JVM heap size: `-Xmx16g -Xms4g`
- Reduce thread count in parameters
- Use smaller databases if possible
- Close other applications

#### 5. File Path Issues
**Problem**: File not found errors
**Solution**:
- Use absolute paths
- Check file permissions
- Verify file formats are supported

### Common Issues - DIA-Specific

#### 1. Predicted Peptides Not Generating
**Problem**: Prediction stalls or fails
**Symptoms**:
- Progress bar stops moving
- Console shows errors
- DeepDetect task fails

**Solutions**:
- **Check Python installation**: Python must be configured in Resources tab
- **Check .NET Runtime**: DIA-NN requires .NET SDK 8.0.407 or later for .raw file reading
  - Download from: https://dotnet.microsoft.com/en-us/download/dotnet/8.0
- **Check DeepDetect path**: Verify DeepDetect.exe is in Resources
- **Memory issues**: Reduce thread count (try 30-40% of cores instead of 60%)
- **Resume partial prediction**: Click "Check Status" then "Start Prediction" again

#### 2. Library Generation Fails
**Problem**: "Generate Library" button doesn't work or fails
**Causes & Solutions**:
- **DIA-NN not configured**: Check Resources tab → DIANN path must be set
- **Invalid input file**: Ensure file format is correct (.speclib, .txt, .fasta)
- **Enzyme parameters**: For protein databases, check enzyme selection
- **Disk space**: Ensure sufficient space in MAG database folder
- **Permissions**: Check write permissions on MAG database directory

**Debug steps**:
1. Check library generation console for specific error messages
2. Verify MAG database is selected
3. Try with a small test file first
4. Check logs in `logs/metapilot.log`

#### 3. Model Training Fails or Hangs
**Problem**: ML model generation doesn't complete
**Symptoms**:
- Training progress stops
- No correlation values shown
- Python errors in console

**Solutions**:
- **Python + PyTorch setup**: Ensure Python configured with PyTorch, pandas, numpy installed
  ```bash
  pip install torch pandas numpy scikit-learn
  ```
- **Training data format**: Check peptide intensity file has correct columns
- **Memory**: Model training needs 4-8 GB RAM, ensure sufficient memory
- **Check training data size**: Need at least 100 peptides, ideally 1000+
- **Timeout issues**: Very large datasets may need longer - check if process is still running

#### 4. Iterative Search Fails in Round 1
**Problem**: "Round 1 library search failed" error
**Causes**:
- **No seed libraries selected**: Must select at least one library for initial search
- **Libraries incompatible**: Libraries don't match MAG database
- **RAW file format**: Cannot read .raw files → install .NET SDK 8.0.407+
- **Thread exhaustion**: Too many threads → reduce to 50% of cores

**Recovery**:
1. Check Console tab for specific error message
2. Verify selected libraries exist (Check Available button)
3. Test with a single file first
4. Check predicted peptides are fully generated

#### 5. "Library panel not available" Message
**Problem**: Library generation panel doesn't appear
**Solution**:
- This happens if MAG database isn't selected yet
- Go to **Input/Output** tab → Select MAG database → Select version
- Return to **Parameters** tab → Panel should now appear
- If still missing, check console for errors

#### 6. Models Not Improving Results
**Problem**: Self-model analysis doesn't perform better than iterative
**Possible Reasons**:
- **Insufficient training data**: Round 1 didn't identify enough peptides
- **Poor model quality**: Check correlation values (<0.5 is problematic)
- **Overfitting**: Model trained on very different samples
- **Sample complexity**: Sometimes traditional approach is sufficient

**Recommendations**:
- Check model correlation values in table (should be >0.7)
- Ensure Round 1 identifies >500 peptides for good training
- For poor models, stick with Strategy 2 (iterative without models)

#### 7. DIA Real-time Not Detecting New Files
**Problem**: New .raw files aren't being processed
**Solutions**:
- **Verify monitor folder path**: Must be exact path where files appear
- **File permissions**: Ensure MetaPilot can read the folder
- **File completion**: Real-time waits for files to finish writing
- **Library selection**: Must have libraries selected before starting
- **Restart monitoring**: Stop and restart if files were added before monitoring started

### DIA-Specific Performance Optimization

1. **Thread Allocation**:
   - **Prediction**: 60% of cores (CPU-intensive)
   - **Library generation**: 50-70% of cores
   - **Analysis**: 50-70% of cores
   - Always leave some cores for system operations

2. **Memory Management**:
   - **Library Analysis**: 8-16 GB recommended
   - **Iterative Analysis**: 16-32 GB recommended
   - **Self-Model Analysis**: 32-64 GB recommended
   - Use `-Xmx` flag to allocate: `java -Xmx32g -jar metapilot.jar`

3. **Storage Optimization**:
   - Use **SSD for temporary files** (in output directory)
   - Keep MAG databases on fast storage
   - Can use slower storage for final results archive
   - Clean up intermediate files after successful completion

4. **Processing Time Estimates** (for 20 samples):
   - Library Analysis: 2-4 hours
   - Iterative Analysis: 8-16 hours
   - Self-Model Analysis: 16-36 hours
   - Scale linearly with sample count

### Log Files
- Application logs: `logs/metapilot.log`
- Workflow-specific logs in output directories
- Console output in MetaPilot GUI
- DIA-NN logs in individual result folders
- Model training logs in model directories

## Advanced Usage

### Complete Example: First-Time DIA Analysis Walkthrough

This example walks through a complete DIA analysis from scratch for a human gut microbiome sample.

#### Scenario
- **Sample Type**: Human gut microbiome
- **Data**: 10 DIA .raw files from a LC-MS/MS run
- **Goal**: Identify and quantify microbial proteins
- **Approach**: Iterative analysis (no suitable libraries available)

#### Step-by-Step Process

**Phase 1: Initial Setup (One-Time, ~4-6 hours)**

1. **Install and Configure Tools** (30 minutes):
   ```bash
   # Install .NET SDK 8.0.407+ for .raw file support
   # Download from: https://dotnet.microsoft.com/download/dotnet/8.0
   
   # Install Python with dependencies
   pip install torch pandas numpy scikit-learn
   
   # Download and install DIA-NN
   # From: https://github.com/vdemichev/DiaNN/releases
   ```

2. **Configure Resources in MetaPilot** (10 minutes):
   - Launch MetaPilot
   - Go to **Resources** tab
   - Set Base Path: `C:\MetaProteomics\Resources`
   - Configure paths:
     - DIANN: `C:\MetaProteomics\Resources\DiaNN\diann.exe`
     - PYTHON: `C:\MetaProteomics\Resources\python\python.exe`
     - DEEPDETECT: `C:\MetaProteomics\Resources\DeepDetect\DeepDetect.exe`
     - FUNCDB: `C:\MetaProteomics\Resources\function\func_def.db`
   - Click **Save**

3. **Setup MAG Database** (20 minutes):
   - Go to **Database Configuration** tab
   - Select **human-gut v2.0.2** from MGnify databases
   - Click **Download** if not available locally
   - Wait for download to complete
   - Verify status shows "Downloaded" and "Local Available"

4. **Generate Predicted Peptides** (4-6 hours):
   - Launch **DIA Normal** workflow
   - **Input/Output** tab:
     - Select MAG database: human-gut v2.0.2
   - Check **MAG Library Status** section
   - Click **Check Status** (should show 0% if first time)
   - Click **Start Prediction**
   - Monitor progress (will process ~13,000 genomes)
   - **Go do something else** - this takes several hours
   - Check back periodically - status should reach 100%

**Phase 2: Create Initial Library** (Optional but recommended, ~30-60 minutes):

5. **Generate Broad Library from IGC Database**:
   - Go to **Parameters** tab
   - **Spectral Library Generation** section:
     - Input Type: "Protein Database File"
     - Click Browse: Select IGC database file (e.g., `IGC.pep.fasta`)
     - Enzyme: Trypsin/P
     - Max Miss Cleavages: 1
     - Library Name: "human_gut_IGC_broad"
     - Click **Generate Library**
   - Wait for generation (~30 minutes for large database)
   - Once complete, library appears in "Available Libraries" table

**Phase 3: Configure and Run Analysis** (~8-12 hours):

6. **Configure Analysis Parameters**:
   - **Input/Output** tab:
     - Click **Add Files** → Select your 10 DIA .raw files
     - For each file, enter experiment name (Sample1, Sample2, etc.)
     - Set **Output Directory**: `C:\Results\HumanGut_DIA_Analysis`
     - MAG Database: human-gut v2.0.2 (should already be selected)
     - Append Host DB: Check (optional)
     - Host Database: Browse to human protein database

7. **Set Analysis Strategy**:
   - **Parameters** tab:
     - **Analysis Mode**: Select "Iterative search"
     - **Thread Count**: 8 (assuming 12-core system)
     - **Spectral Library Generation**:
       - Check "Use" for "human_gut_IGC_broad" library
       - This will seed the iterative search
     - **ML Model Generation**: Leave empty for first analysis
     - **Metadata Configuration**: Fill in sample metadata (optional)

8. **Start Analysis**:
   - Go to **Console** tab
   - Review parameters in the parameter view
   - Click **Start Analysis**
   - **Do NOT close MetaPilot during analysis**

9. **Monitor Progress**:
   - Overall progress bar shows main stages
   - Current task progress bar shows sub-tasks
   - Console shows detailed progress:
     - Converting spectra: 10-15%
     - Round 1 search: 15-35%
     - Round 2 search: 35-50%
     - Round 3 search: 50-65%
     - Round 4 search: 65-80%
     - Final identification: 80-85%
     - Functional annotation: 85-92%
     - Taxonomy analysis: 92-95%
     - Report generation: 95-100%

10. **Analysis Complete**:
    - Console shows "✓ Analysis completed successfully!"
    - Click **View Report** to open interactive dashboard
    - Results in `C:\Results\HumanGut_DIA_Analysis\iterative_analysis\`

**Phase 4: Review Results** (30 minutes):

11. **Explore Reports**:
    - **ID Summary**: Check identification rates per sample
    - **Peptides Summary**: Review peptide identifications and quantification
    - **Proteins Summary**: Explore protein groups
    - **Taxonomy Summary**: View microbial composition
    - **Function Summary**: Analyze functional pathways (COG, KEGG, GO)

12. **Quality Control**:
    - Check `final_summary.tsv`:
      - MS/MS identified should be >60% for good quality
      - Peptide count should be >1000 per sample for complex microbiomes
    - Review `Genome.tsv`:
      - Should see expected microbial taxa
      - Check if p-values and FDR values are acceptable

**Phase 5: Future Analyses** (Much Faster!):

13. **For Next Similar Samples**:
    - **Predicted peptides**: Already done! ✓
    - **Libraries**: Already created! ✓
    - **Just configure**:
      - Load previous parameter file: File → Load Parameter File
      - Update raw files and output directory
      - **Start Analysis** - Will complete much faster (~2-4 hours)

14. **Optional Enhancement - Build ML Model**:
    - Use results from successful analysis to train model
    - Export peptide matrix from results
    - Generate model using peptide intensity data
    - Use model for future analyses to improve sensitivity

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

## Quick Reference Guide

### DIA Workflow Selection Cheat Sheet

**Choose your workflow based on your situation:**

| Your Situation | Recommended Workflow | Key Steps |
|----------------|---------------------|-----------|
| Have suitable spectral libraries | DIA Normal - Library Analysis | Select libraries → Start Analysis |
| Novel sample, no libraries | DIA Normal - Iterative Analysis | Select seed libraries → 4-round search |
| Complex/challenging sample | DIA Normal - Self-Model Analysis | Enable deep learning → Auto-model generation |
| Focus on specific taxa | DIA Target | Select taxa → Generate targeted library |
| Live instrument monitoring | DIA Real-time | Set monitor folder → Watch for files |

### Essential Checklist Before DIA Analysis

```
□ Resources tab: DIA-NN path configured
□ Resources tab: Python path configured (for models/reports)
□ Resources tab: DeepDetect path configured (for prediction)
□ Database tab: MAG database downloaded and available
□ Input/Output tab: MAG database selected and version chosen
□ Input/Output tab: Predicted peptides generated (check status = 100%)
□ Parameters tab: At least one library selected OR iterative mode chosen
□ Parameters tab: Thread count set (50-70% of available cores)
□ Console tab: Output directory has sufficient space (>100 GB)
□ Install .NET SDK 8.0.407+ for .raw file support
```

### DIA File Structure Reference

After a successful DIA analysis, your output directory will contain:

#### Library Analysis Results:
```
output_directory/
└── library_analysis/
    ├── combined/
    │   ├── combined.parquet            # DIA-NN search results
    │   ├── combined.pr_matrix.tsv      # Peptide quantification matrix
    │   ├── combined.pg_matrix.tsv      # Protein quantification matrix
    │   └── combined.stats.tsv          # Search statistics
    ├── final_peptides.tsv              # Final peptide report
    ├── final_proteins.tsv              # Final protein report
    ├── functional_annotation/
    │   ├── functions.xml               # Functional annotations (XML)
    │   ├── functions.html              # Functional report (HTML)
    │   └── [function category reports]
    ├── taxonomy_analysis/
    │   ├── Genome.tsv                  # Genome-level taxonomy
    │   ├── Genome.html                 # Genome report (HTML)
    │   └── Taxa_report.tsv             # Taxonomic abundance
    └── report/
        └── reports/
            ├── index.html              # Main report dashboard
            ├── report_ID_summary.html
            ├── report_peptides_summary.html
            ├── report_proteinGroups_summary.html
            ├── report_taxonomy_summary.html
            └── report_function_summary.html
```

#### Iterative Analysis Results:
```
output_directory/
└── iterative_analysis/
    ├── Sample1/                        # Per-sample folders
    │   ├── firstSearch/
    │   │   ├── firstSearch.parquet
    │   │   ├── firstSearch.fasta
    │   │   └── mags.tsv                # MAGs used in Round 1
    │   ├── secondSearch/
    │   │   ├── secondSearch.parquet
    │   │   ├── secondSearch.fasta
    │   │   └── mags.tsv                # Refined MAGs for Round 2
    │   ├── thirdSearch/
    │   │   └── [similar structure]
    │   └── lastSearch.parquet          # Final Round 4 results
    ├── Sample2/
    │   └── [same structure]
    ├── combined/                       # Combined results across all samples
    │   ├── combined.parquet
    │   ├── combined.pr_matrix.tsv
    │   ├── combined.pg_matrix.tsv
    │   ├── genomes_razor.tsv           # Genome coverage analysis
    │   └── proteins_razor.tsv          # Protein coverage analysis
    └── [final reports same as library analysis]
```

#### Self-Model Analysis Results:
```
output_directory/
└── self_model_analysis/
    ├── Sample1/
    │   ├── firstSearch/
    │   ├── model/                      # Sample-specific ML model
    │   │   ├── model_pep_rank.tsv      # Peptide rankings
    │   │   ├── model_pytorch_predictions.tsv
    │   │   └── best_model_pytorch.pt   # Trained PyTorch model
    │   ├── secondSearch/
    │   ├── thirdSearch/
    │   └── lastSearch.parquet
    └── [similar structure for other samples and combined results]
```

### Key Output Files Explained

| File | Description | Use For |
|------|-------------|---------|
| `final_peptides.tsv` | All identified peptides with quantification | Peptide-level analysis |
| `final_proteins.tsv` | Protein groups with intensities | Protein-level analysis |
| `final_summary.tsv` | Overview of identifications per sample | QC and overview |
| `Genome.tsv` | Genome-level taxonomic assignments | Taxonomic analysis |
| `Taxa_report.tsv` | Abundance at all taxonomic ranks | Microbiome composition |
| `functions.tsv` | Functional annotations (COG, KEGG, GO, EC) | Functional analysis |
| `*.pr_matrix.tsv` | Precursor intensity matrix | Downstream quantification |
| `*.pg_matrix.tsv` | Protein group intensity matrix | Downstream statistics |
| `report/reports/index.html` | Interactive HTML dashboard | Visual exploration |

### Command Reference

#### Launch MetaPilot GUI:
```bash
mvn javafx:run
# OR
java --module-path $PATH_TO_FX --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base -jar target/metapilot-*-shaded.jar
```

#### Launch MetaPilot CLI:
```bash
java -jar target/metapilot-*-shaded.jar --cli --params config/parameters.json --resources config/meta_resources.params
```

#### Allocate More Memory:
```bash
java -Xmx32g -Xms8g -jar target/metapilot-*-shaded.jar
```

### Keyboard Shortcuts (GUI)

| Shortcut | Action |
|----------|--------|
| Ctrl+O | Load parameter file |
| Ctrl+S | Save current parameters |
| Ctrl+R | Check resources configuration |
| Ctrl+D | Open database configuration |
| Ctrl+Q | Quit application |

---

## Version Information

This manual covers MetaPilot v2.0. For the most up-to-date information, refer to the application's built-in help and documentation.

### Version History
- **v2.0**: Added DIA workflows, ML model support, improved MAG database management
- **v1.0**: Initial release with DDA workflows (MaxQuant, pFind, FragPipe, Sage)

### Citation
If you use MetaPilot in your research, please cite:
```
MetaPilot: A Comprehensive Metaproteomics Analysis Platform
Kai Cheng, Quadram Institute Bioscience
Version 2.0 (2024)
```

### Contact and Support
- **GitHub Issues**: [Report bugs or request features]
- **Email**: [Contact information]
- **Documentation**: Check `docs/` directory for additional guides
- **Updates**: Use Help → Check for Updates in the application menu

# MetaPilot

MetaPilot is a next-generation proteomics data analysis platform, building upon the foundation of MetaLab with enhanced features and improved performance.

## Project Structure

```
MetaPilot/
├── src/
│   ├── main/
│   │   ├── java/
│   │   │   └── org/metapilot/
│   │   │       ├── model/        # Deep learning models
│   │   │       ├── data/         # Data processing
│   │   │       ├── utils/        # Utility classes
│   │   │       └── Main.java     # Main entry point
│   │   └── resources/            # Configuration files
│   └── test/
│       ├── java/                 # Test classes
│       └── resources/            # Test resources
├── data/                         # Data files
└── pom.xml                       # Maven configuration
```

## Features

- Advanced peptide intensity prediction using deep learning
- Support for DeepLearning4j framework
- Improved data preprocessing and feature engineering
- Enhanced visualization capabilities
- Comprehensive documentation and examples

## Requirements

- Java 11 or higher
- Maven 3.6 or higher
- DeepLearning4j dependencies (automatically managed by Maven)

## Building the Project

```bash
mvn clean install
```

## Running the Application

```bash
mvn exec:java -Dexec.mainClass="org.metapilot.Main"
```

## Development Roadmap

### 🎯 Current Priority: HPC Migration (2024)
MetaPilot is being enhanced to support high-performance computing (HPC) environments:

- **📋 [HPC Implementation Roadmap](docs/HPC_Implementation_Roadmap.md)** - Quick start guide
- **📖 [Complete Technical Report](docs/MetaPilot_HPC_Migration_Plan.md)** - Comprehensive migration plan

**Key Goals:**
- ✅ **Dual-Mode Architecture**: Desktop GUI + HPC CLI
- ✅ **Linux Server Support**: Headless execution on HPC clusters  
- ✅ **Scalable Performance**: 10x throughput improvement
- ✅ **File Format Compatibility**: Windows .raw → Linux .mzML workflow

### Future Development

1. Core Features
   - [x] Enhanced data preprocessing pipeline
   - [x] Improved model architecture
   - [ ] Advanced visualization tools
   - [ ] Performance optimization for HPC

2. HPC Features (In Progress)
   - [ ] Headless JAR for Linux servers
   - [ ] SLURM/PBS job integration
   - [ ] Parameter file workflow
   - [ ] Web-based result visualization

3. Desktop Features
   - [ ] HPC parameter export wizard
   - [ ] Remote job monitoring
   - [ ] File upload/download tools
   - [ ] Enhanced GUI workflows

4. Documentation
   - [x] HPC migration technical report
   - [x] Linux deployment guides
   - [ ] API documentation
   - [ ] User tutorials

## Contributing

(Contribution guidelines will be added)

## License

(License information will be added) 
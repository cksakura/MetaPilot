# MetaPilot CLI Reference Manual

## Overview

MetaPilot CLI provides command-line interface for running workflows without a GUI. This is useful for:
- Automated batch processing
- HPC job scheduling
- Server-based execution
- Testing and validation

## Available Workflows

### Supported Workflows

Currently supported workflows:

1. **DDA IGC MaxQuant** - Label-free quantification
2. **DDA IGC pFind** - Chinese search engine
3. **DDA IGC FragPipe** - MSFragger pipeline
4. **DDA IGC Sage** - Rust-based fast search
5. **DDA MAG pFind** - MAG database search
6. **DDA MAG FragPipe** - MAG with FragPipe
7. **DDA MAG Sage** - MAG with Sage
8. **DIA Normal** - Standard DIA workflow
9. **DIA Target** - Targeted DIA workflow
10. **DIA Real-Time** - Real-time DIA workflow

**Note**: Only FragPipe, Sage, and DIA-NN workflows are supported on Linux.

## Basic Usage

### Standard Workflow Execution

```bash
java -jar metapilot-windows-gui-1.0.0-beta.2.jar --cli \
    --params /path/to/parameters.json \
    --resources /path/to/resources.params
```

### Parameters

- `--params <file>`: Parameter JSON file
- `--resources <file>`: Resource configuration file
- `--help`: Show help message

### Example

```bash
java -jar metapilot-windows-gui-1.0.0-beta.2.jar --cli \
    --params config/last_dda_igc_maxquant_parameters.json \
    --resources config/meta_resources.params
```

## Parameter File Format

### Structure

Parameter files are JSON format:

```json
{
    "metadata": {
        "rawFiles": [
            "/data/raw/file1.raw",
            "/data/raw/file2.raw"
        ],
        "fastaFile": "/data/database.fasta"
    },
    "workflow": "DDA_IGC_MAXQUANT",
    "result": "/data/results/output"
}
```

### Required Fields

All workflows require:

- `metadata.rawFiles`: Array of raw mass spectrometry files
- `metadata.fastaFile`: Database FASTA file
- `result`: Output directory

### Workflow-Specific Parameters

See individual workflow documentation for specific parameters.

## Resource Configuration

### File Format

Resource file (`meta_resources.params`):

```properties
# Tool Paths
FRAGPIPE_PATH=/opt/tools/fragpipe
SAGE_PATH=/usr/local/bin
DIANN_PATH=/opt/tools/diann

# Database Paths
MAG_DB_DIR=/opt/metapilot/data/mag_databases

# Temporary Directory
TEMP_DIR=/tmp/metapilot

# Python
PYTHON_PATH=/usr/bin/python3
```

### Configuration Keys

| Key | Description | Example |
|-----|-------------|---------|
| `FRAGPIPE_PATH` | FragPipe installation | `/opt/tools/fragpipe` |
| `SAGE_PATH` | Sage binary path | `/usr/local/bin` |
| `DIANN_PATH` | DIA-NN installation | `/opt/tools/diann` |
| `MAG_DB_DIR` | MAG database directory | `/opt/metapilot/data/mag_databases` |
| `TEMP_DIR` | Temporary files | `/tmp/metapilot` |
| `PYTHON_PATH` | Python executable | `/usr/bin/python3` |

## Workflow Examples

### DDA IGC FragPipe

```bash
java -jar metapilot.jar --cli \
    --params config/dda_igc_fragpipe.json \
    --resources config/meta_resources.params
```

Parameter file:
```json
{
    "metadata": {
        "rawFiles": ["sample1.raw", "sample2.raw"],
        "fastaFile": "database.fasta"
    },
    "result": "/data/results/fragpipe",
    "workflow": "DDA_IGC_FRAGPIPE"
}
```

### DDA MAG Sage

```bash
java -jar metapilot.jar --cli \
    --params config/dda_mag_sage.json \
    --resources config/meta_resources.params
```

Parameter file:
```json
{
    "metadata": {
        "rawFiles": ["mag_sample1.raw"],
        "fastaFile": "mgyg-catalogue-v1.fasta"
    },
    "magCatalogue": "mgyg-catalogue-v1",
    "result": "/data/results/mag_sage",
    "workflow": "DDA_MAG_SAGE"
}
```

### DIA Normal

```bash
java -jar metapilot.jar --cli \
    --params config/dia_normal.json \
    --resources config/meta_resources.params
```

Parameter file:
```json
{
    "metadata": {
        "rawFiles": ["dia_sample.raw"],
        "fastaFile": "database.fasta"
    },
    "result": "/data/results/dia",
    "workflow": "DIA_NORMAL"
}
```

## Output Structure

### Results Directory

```
results/
├── metainfo.txt                 # Workflow metadata
├── peptides.txt                # Peptide identifications
├── proteins.txt                # Protein identifications
├── report.html                 # HTML report
└── report_data/                # Report data files
    ├── function_summary.tsv
    ├── taxonomy_summary.tsv
    └── ...
```

### Metadata File

`metainfo.txt` contains:
- Workflow parameters
- Execution timestamp
- Input files
- Tool versions
- Processing statistics

## Batch Processing

### Processing Multiple Files

Create manifest file (`batch_manifest.txt`):

```
# Format: param_file|resource_file|result_dir
config/job1_params.json|config/meta_resources.params|results/job1
config/job2_params.json|config/meta_resources.params|results/job2
config/job3_params.json|config/meta_resources.params|results/job3
```

Run batch:

```bash
java -cp metapilot.jar uk.ac.quadram.metapilot.BatchRunner batch_manifest.txt
```

### HPC Job Submission (SLURM)

Example SLURM script (`submit_job.sh`):

```bash
#!/bin/bash
#SBATCH --job-name=metapilot
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=metapilot_%j.out

module load java/17

java -Xmx12g -jar metapilot.jar --cli \
    --params ${1:-config/params.json} \
    --resources ${2:-config/resources.params}
```

Submit job:

```bash
sbatch submit_job.sh config/my_params.json config/meta_resources.params
```

## Environment Variables

### JVM Options

Set via environment:

```bash
export JAVA_OPTS="-Xms512m -Xmx8g -XX:+UseG1GC"
java -jar metapilot.jar --cli --params ...
```

### Tool Paths

Override via environment:

```bash
export FRAGPIPE_PATH=/usr/local/fragpipe
export SAGE_PATH=/usr/local/bin/sage
java -jar metapilot.jar --cli --params ...
```

## Error Handling

### Common Errors

1. **File Not Found**: Check file paths are absolute or relative to working directory
2. **Tool Not Found**: Verify tool paths in resources file
3. **Out of Memory**: Increase `-Xmx` parameter
4. **Permission Denied**: Check file permissions

### Debugging

Enable verbose logging:

```bash
java -jar metapilot.jar --cli \
    --params params.json \
    --resources resources.params \
    --debug
```

Log file: `metapilot.log` in working directory

## Return Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Invalid parameters |
| 3 | Tool execution failed |
| 4 | Out of memory |

## Performance Tips

### Memory Allocation

```bash
# For small datasets (< 10 files)
java -Xmx4g -jar metapilot.jar --cli

# For medium datasets (10-50 files)
java -Xmx8g -jar metapilot.jar --cli

# For large datasets (> 50 files)
java -Xmx16g -jar metapilot.jar --cli
```

### Parallel Processing

Workflows automatically use available CPUs. Limit via:

```bash
# Set available cores
export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=8"
java -jar metapilot.jar --cli
```

### I/O Optimization

Use fast storage for:
- Input raw files
- Temporary files
- Output results

Prefer SSD over HDD for better performance.

## Integration Examples

### Python Script

```python
import subprocess
import json

def run_metapilot(params_file, resources_file):
    cmd = [
        'java', '-Xmx8g', '-jar', 'metapilot.jar',
        '--cli',
        '--params', params_file,
        '--resources', resources_file
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0

# Usage
success = run_metapilot('config/params.json', 'config/resources.params')
print(f"Workflow {'completed' if success else 'failed'}")
```

### Bash Script

```bash
#!/bin/bash
# run_metapilot.sh

if [ $# -lt 2 ]; then
    echo "Usage: $0 <params_file> <resources_file>"
    exit 1
fi

JAVA_OPTS="-Xmx8g"
PARAMS="$1"
RESOURCES="$2"

java $JAVA_OPTS -jar metapilot.jar --cli \
    --params "$PARAMS" \
    --resources "$RESOURCES"

echo "Exit code: $?"
```

## Troubleshooting

### Tool Paths

Verify paths in resources file:

```bash
# Test FragPipe
ls -l $(grep FRAGPIPE_PATH config/meta_resources.params | cut -d'=' -f2)

# Test Sage
sage --version

# Test DIA-NN
ls -l $(grep DIANN_PATH config/meta_resources.params | cut -d'=' -f2)
```

### File Permissions

Ensure files are readable:

```bash
chmod 644 config/meta_resources.params
chmod 755 /opt/tools/fragpipe
```

### Disk Space

Check available space:

```bash
df -h /tmp
df -h /data/results
```

## Advanced Usage

### Custom Workflows

Create custom parameter files for specific analysis needs:

```json
{
    "metadata": {
        "rawFiles": ["file1.raw", "file2.raw"],
        "fastaFile": "database.fasta"
    },
    "workflow": "DDA_IGC_FRAGPIPE",
    "result": "/results/custom",
    "advancedSettings": {
        "precursorTolerance": 20.0,
        "fragmentTolerance": 0.02
    }
}
```

### Workflow Chaining

Run multiple workflows sequentially:

```bash
#!/bin/bash
# Chain multiple workflows

# Step 1: DDA search
java -jar metapilot.jar --cli --params dda.json --resources resources.params

# Step 2: DIA analysis
java -jar metapilot.jar --cli --params dia.json --resources resources.params

# Step 3: Post-processing
# ... custom scripts
```

## Quick Reference

### Essential Commands

```bash
# Show help
java -jar metapilot.jar --help

# Run workflow
java -jar metapilot.jar --cli --params <params> --resources <resources>

# With custom memory
java -Xmx16g -jar metapilot.jar --cli --params <params> --resources <resources>

# Debug mode
java -jar metapilot.jar --cli --params <params> --resources <resources> --debug
```

### File Locations

| Purpose | Default Location |
|---------|------------------|
| Parameters | `config/` |
| Resources | `config/meta_resources.params` |
| Results | As specified in params |
| Logs | `logs/metapilot.log` |
| Tools | `/opt/tools/` |

## Support

For additional help:

- Review log files: `logs/metapilot.log`
- Check documentation: `docs/`
- Contact: MetaPilot development team

---

**Version**: 1.0.0-beta.3  
**Last Updated**: 2024


# MetaPilot HPC Deployment

This guide describes the current HPC deployment model for MetaPilot `1.0.0-beta.8`.

## Current Deployment Model

Use `metapilot-cli` as the production headless entry point. It validates GUI-exported JSON parameter files, resolves site resources and MAG databases, and dispatches supported core runners.

Supported headless workflows:

- `DIA_MARKER`: Core/Focus, Comprehensive, Deep together, Deep per file.
- `DDA_MAG_FRAGPIPE`.

The Windows GUI remains the primary interface for configuring and running all GUI workflows. The web module is a server wrapper and API/UI surface, not a replacement for the CLI execution engine.

## Recommended Directory Layout

```text
/opt/metapilot/
  metapilot-cli-1.0.0-beta.8-shaded.jar
  metapilot-web-1.0.0-beta.8.jar          # optional

/data/metapilot/
  site/
    meta_resources.params
    MAGDB_1_0.json
  MAG/
    human-gut/v2.0.2/
    mouse-gut/v1.0/
  function/
  studies/
```

Keep raw data, result folders and MAG databases on shared storage visible to compute nodes.

## Tool Paths

Example `meta_resources.params` values:

```properties
DIANN=/opt/tools/diann/diann-linux
FRAGPIPE=/opt/tools/fragpipe-24.0/bin/fragpipe
FRAGPIPE_TOOLS=/opt/tools/fragpipe-24.0/tools
PYTHON=/opt/conda/envs/metapilot/bin/python
FUNCDB=/data/metapilot/function
DEEPDETECT=/opt/tools/deepdetect/DeepDetect
```

Validate before submitting jobs:

```bash
java -jar /opt/metapilot/metapilot-cli-1.0.0-beta.8-shaded.jar resources-check \
  --resources /data/metapilot/site/meta_resources.params
```

## SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=metapilot-dia-marker
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=metapilot-%j.out
#SBATCH --error=metapilot-%j.err

set -euo pipefail

JAR=/opt/metapilot/metapilot-cli-1.0.0-beta.8-shaded.jar
PARAMS=/data/metapilot/studies/study_parameters.json
RESOURCES=/data/metapilot/site/meta_resources.params
MAGDB=/data/metapilot/site/MAGDB_1_0.json

java -Xmx220G -jar "$JAR" params-validate \
  --params "$PARAMS" \
  --resources "$RESOURCES" \
  --magdb "$MAGDB"

java -Xmx220G -jar "$JAR" run \
  --params "$PARAMS" \
  --resources "$RESOURCES" \
  --magdb "$MAGDB"
```

## Deep Per File Controls

DIA Deep per file can run multiple per-raw-file workers. Tune concurrency for node memory:

```bash
export METAPILOT_DEEP_PER_FILE_POOL=2
export METAPILOT_DEEP_PER_FILE_RAM_PER_TASK_GB=80
```

Use a conservative pool first, then increase only after measuring peak memory on your data.

## Web/API Wrapper

The `metapilot-web` module can be deployed when users need a browser/API front end:

```bash
mvn -pl metapilot-web -am -DskipTests package
java -jar metapilot-web/target/metapilot-web-1.0.0-beta.8.jar
```

Configure network access, upload limits and logging through `metapilot-web/src/main/resources/application.properties` or external Spring Boot configuration.

For production execution, prefer submitting CLI jobs through the scheduler rather than running long analyses inside the web JVM.

## Operational Checklist

- Build and record the exact MetaPilot commit.
- Pin DIA-NN, FragPipe and DeepDetect versions.
- Store `meta_resources.params` and `MAGDB_1_0.json` with the study outputs.
- Keep GUI-exported parameter JSON files with results.
- Capture `logging.log` from the result directory.
- Keep Marker4 and predicted peptide artifacts with the MAG database or regenerate them deterministically.
- Ensure each compute node can read all raw files, libraries, databases and tool installations.

## Current Limits

- Generic CLI execution is limited to `DIA_MARKER` and `DDA_MAG_FRAGPIPE`.
- DIA Seed, DIA Target and DIA Real-Time are not production headless scheduler workflows yet.
- MaxQuant is Windows-oriented.
- The web wrapper launches a JAR subprocess and should not be documented as the canonical HPC runner.

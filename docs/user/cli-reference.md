# MetaPilot CLI Reference

This page describes the current headless CLI surface for MetaPilot `1.1.0`.

The CLI module is intended for Linux/HPC operation, but only a subset of workflows has headless runners today. The Windows GUI remains the broadest workflow surface.

## Release Package

Use the CLI JAR included in the MetaPilot release ZIP. Public releases are distributed as ready-to-run artifacts; building from source is not part of the supported user installation path.

Example release artifact:

```text
MetaPilot-1.1.0-linux-cli.zip
```

Example extracted CLI JAR:

```text
metapilot-cli-1.1.0-shaded.jar
```

## Commands

```bash
java -jar metapilot-cli-1.1.0-shaded.jar --help
java -jar metapilot-cli-1.1.0-shaded.jar --version
```

Current commands:

| Command | Purpose | Current status |
| --- | --- | --- |
| `init` | Create the default site configuration files. | Supported |
| `resources-check` | Validate configured tool paths. | Supported |
| `resources-set` | Set or print resource paths. | Supported |
| `mag-register` | Register a local MAG database in the site registry. | Supported |
| `mag-list` | List registered MAG databases. | Supported |
| `mag-download` | Download/register supported MAG catalogues. | Supported where catalogue metadata is available |
| `params-validate` | Validate a GUI-exported parameter JSON before running. | Supported for headless workflows below |
| `run` | Run a GUI-exported parameter JSON through the generic study runner. | Supported for `DIA_MARKER` and `DDA_MAG_FRAGPIPE` |
| `dia-marker` | Direct headless runner for DIA Marker workflows. | Supported |
| `mag-fragpipe` | Direct headless runner for DDA MAG FragPipe. | Supported |
| `deepdetect-predict` | Run DeepDetect prediction for a FASTA. | Supported |
| `deepdetect-magdb` | Prepare DeepDetect predictions for a MAG database. | Supported |

## Headless Workflow Support

The generic `run` and `params-validate` commands currently execute:

- `DIA_MARKER`: Core/Focus, Comprehensive, Deep together, and Deep per file.
- `DDA_MAG_FRAGPIPE`: MAG-aware DDA FragPipe workflow.

Other workflow parameter types may be registered for JSON load/save compatibility but are not yet executable through the headless generic runner.

The GUI supports the wider workflow set:

- DDA IGC: MaxQuant, pFind, FragPipe, Sage.
- DDA MAG: pFind, FragPipe, Sage.
- DIA: Marker, Seed, Target, Real-Time.

Seed Deep per file is exposed in the canonical mode model but is still blocked by GUI validation as unsupported.

## Site Configuration

Initialize site configuration:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar init
```

Useful resource keys include:

| Key | Typical value |
| --- | --- |
| `DIANN` | `/opt/tools/diann/diann-linux` |
| `FRAGPIPE` | `/opt/tools/fragpipe-24.0/bin/fragpipe` |
| `FRAGPIPE_TOOLS` | `/opt/tools/fragpipe-24.0/tools` |
| `PYTHON` | Python executable used for report generation helpers |
| `FUNCDB` | Functional annotation resource directory |
| `DEEPDETECT` | DeepDetect executable |

Set values:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar resources-set \
  --set DIANN=/opt/tools/diann/diann-linux \
  --set FRAGPIPE=/opt/tools/fragpipe-24.0/bin/fragpipe \
  --set FRAGPIPE_TOOLS=/opt/tools/fragpipe-24.0/tools \
  --print
```

Validate values:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar resources-check
```

## Running Parameter Files

Export a parameter JSON from the GUI, then validate it:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar params-validate \
  --params /path/to/study_parameters.json
```

Run it:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar run \
  --params /path/to/study_parameters.json
```

Override resource or MAG database registry paths:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar run \
  --params /path/to/study_parameters.json \
  --resources /path/to/meta_resources.params \
  --magdb /path/to/MAGDB_1_0.json
```

The runner writes `logging.log` in the study output directory when it can resolve one from the parameter file.

## Direct DIA Marker Runner

```bash
java -jar metapilot-cli-1.1.0-shaded.jar dia-marker \
  --params /path/to/dia_marker_parameters.json \
  --resources /path/to/meta_resources.params \
  --magdb /path/to/MAGDB_1_0.json
```

Supported marker modes are persisted as:

- `FOCUS`: Core / Focus.
- `COMPREHENSIVE`: Comprehensive.
- `DEEP`: Deep together.
- `DEEP_PER_FILE`: Deep per file.

`DIRECT` is not valid for DIA Marker; it belongs to DIA Seed.

## Direct DDA MAG FragPipe Runner

```bash
java -jar metapilot-cli-1.1.0-shaded.jar mag-fragpipe \
  --params /path/to/dda_mag_fragpipe_parameters.json \
  --resources /path/to/meta_resources.params \
  --magdb /path/to/MAGDB_1_0.json
```

FragPipe must have the required installed workflows available under the configured installation.

MAG FragPipe per-sample searches default to memory-safe concurrency using a 60 GB per task estimate. Advanced users can override the concurrent task count:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar run \
  --params /path/to/dda_mag_fragpipe_parameters.json \
  --fragpipe-parallel-tasks 4
```

The same override is available on `mag-fragpipe`. If the requested value is above the memory-safe estimate, MetaPilot logs a warning and proceeds with the requested concurrency.

## MAG Database Registry

List registered MAG databases:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar mag-list
```

Download and prepare a supported MAG catalogue before registration:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar mag-download \
  --catalogue human-gut \
  --version v2.0.2 \
  --out /data/metapilot/MAG \
  --resume
```

Register the downloaded or otherwise prepared local database:

```bash
java -jar metapilot-cli-1.1.0-shaded.jar mag-register \
  --catalogue human-gut \
  --version v2.0.2 \
  --identifier MGYG \
  --path /data/metapilot/MAG/human-gut/v2.0.2
```

The database folder should contain the catalogue files expected by MetaPilot, including `original_db`, `genomes-taxa.tsv` or catalogue metadata, and generated Marker4/predicted files where the selected workflow requires them.

## Exit Codes

| Code | Meaning |
| --- | --- |
| `0` | Success |
| `1` | Runtime failure |
| `2` | Invalid arguments, unsupported workflow, or failed validation |

## Current Limits

- Generic headless `run` is intentionally limited to `DIA_MARKER` and `DDA_MAG_FRAGPIPE`.
- MaxQuant remains Windows-oriented.
- DIA Seed, DIA Target and DIA Real-Time are GUI workflows today; do not document them as headless-supported until runners are added.
- The web module invokes a JAR as a subprocess and should be treated as a server wrapper, not the canonical execution engine.

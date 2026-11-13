# MetaPilot Linux Installation

This guide covers the current Linux/headless setup for MetaPilot `1.2.0`.

The canonical Linux entry point is `metapilot-cli`. The Spring Boot `metapilot-web` module is a web/API wrapper that launches a JAR subprocess and should be deployed only when a browser/API front end is needed.

## Prerequisites

- Linux x86_64 system or HPC login/compute environment.
- Java 17+.
- DIA-NN for DIA Marker.
- FragPipe/MSFragger for DDA MAG FragPipe.
- Python where report-generation helpers are enabled.
- Access to a configured MAG database registry and local MAG database folders.

## Install From Release ZIP

Use the MetaPilot Linux CLI release ZIP. Building from source is not part of the supported user installation path.

Example layout after extraction:

```text
/opt/metapilot/
  metapilot-cli-1.2.0-shaded.jar
  metapilot
  config/
    meta_resources.params.template
    MAGDB_1_0.json
  docs/
  examples/
```

## Initialize Site Configuration

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar init
```

Then set tool paths:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar resources-set \
  --set DIANN=/opt/tools/diann/diann-linux \
  --set FRAGPIPE=/opt/tools/fragpipe-24.0/bin/fragpipe \
  --set FRAGPIPE_TOOLS=/opt/tools/fragpipe-24.0/tools \
  --set PYTHON=/opt/conda/envs/metapilot/bin/python \
  --set FUNCDB=/data/metapilot/function \
  --print
```

Validate:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar resources-check
```

## MAG Database Setup

List registered databases:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar mag-list
```

Download and prepare a supported MAG catalogue before registration:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar mag-download \
  --catalogue human-gut \
  --version v2.0.2 \
  --out /data/metapilot/MAG \
  --resume
```

Register the downloaded or otherwise prepared local database:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar mag-register \
  --catalogue human-gut \
  --version v2.0.2 \
  --identifier MGYG \
  --path /data/metapilot/MAG/human-gut/v2.0.2
```

For DIA Marker, the database must be able to provide or generate:

- `Marker4.fasta`
- `genomes-taxa.tsv`
- `Marker4.predicted.tsv` and `Marker4_highpep.fasta` for large Marker4 sets
- `libraries/Marker4/Marker4.speclib`

## Run A Study

Export a parameter JSON from the GUI or generate one through an existing workflow setup.

Validate:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar params-validate \
  --params /data/studies/study_parameters.json
```

Run:

```bash
java -jar /opt/metapilot/metapilot-cli-1.2.0-shaded.jar run \
  --params /data/studies/study_parameters.json
```

Currently supported by the generic headless runner:

- `DIA_MARKER`
- `DDA_MAG_FRAGPIPE`

## Optional Web Wrapper

If a release package includes the optional web wrapper, deploy that packaged artifact separately when a browser/API front end is needed. The CLI remains the preferred execution path for production headless runs.

## Verification Checklist

- `java -version` reports Java 17 or newer.
- `resources-check` passes for the tools needed by the chosen workflow.
- `params-validate` passes for the study JSON.
- MAG database entries resolve from the configured registry.
- The result directory is writable.
- For Deep per file, memory limits are high enough for concurrent DIA-NN workers.

## Troubleshooting

- `workflow is registered but not supported by the headless CLI yet`: use the GUI or add a headless runner before deploying that workflow.
- `DIANN is not set or does not exist`: fix the `DIANN` key in the resource config.
- `FRAGPIPE is not set or does not exist`: fix `FRAGPIPE` and ensure required workflow files exist in the FragPipe installation.
- `MAG database entry could not be resolved`: check `MAGDB_1_0.json`, catalogue id, identifier and version.
- Missing Marker4 library: generate Marker4 prerequisites from the MAG database configuration before running DIA Marker.

# MetaPilot CLI

Headless command-line runner for MetaPilot, intended for Linux/HPC environments.

## Build

From repo root:

```bash
mvn -DskipTests clean package
```

## Run

```bash
java -jar metapilot-cli/target/metapilot-cli-*-shaded.jar --help
java -jar metapilot-cli/target/metapilot-cli-*-shaded.jar --version
```

## Status

This module provides the headless entry point for the upcoming MetaPilot `1.2.0` release.

Current headless execution supports:

- `DIA_MARKER`: Core/Focus, Comprehensive, Deep together, and Deep per file modes.
- `DDA_MAG_FRAGPIPE`: MAG-aware DDA FragPipe workflow.

Use `--help` to list the current command surface, including resource setup, MAG database registration, parameter validation, and workflow execution commands.


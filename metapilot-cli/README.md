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

This module currently provides CLI scaffolding only. Workflow execution will be added as workflow/task code is refactored to be JavaFX-independent (using `metapilot-core`).


# MetaPilot Architecture

MetaPilot uses a Maven multi-module structure. The root project is a parent POM only; production code belongs inside the module directories.

## Repository Layout

```text
MetaPilot/
├── pom.xml                    # Parent POM
├── metapilot-core/            # Shared domain logic and reusable workflow services
├── metapilot-cli/             # Headless command-line and HPC entry point
├── metapilot-windows-gui/     # Windows desktop application
├── metapilot-web/             # Spring Boot web application
└── docs/                      # Documentation
```

Current version target: `1.0.0-beta.8`.

## Modules

### `metapilot-core`

Shared implementation used by the GUI, CLI, and web modules. Core code should not depend on JavaFX, Spring Web UI classes, or platform-specific shell behavior unless it is explicitly isolated behind an adapter.

Typical ownership:

- Workflow orchestration that can run outside the GUI.
- Search engine readers and result parsers.
- MAG database selection, taxonomy, and functional export logic.
- Shared parameter and resource models.

### `metapilot-cli`

Headless entry point for Linux/HPC usage. CLI code should focus on argument parsing, configuration loading, validation, and invoking core runners.

Typical ownership:

- Command definitions and CLI help.
- Batch/HPC study runner commands.
- Resource and parameter file handoff to core services.

Current generic headless execution support is limited to `DIA_MARKER` and `DDA_MAG_FRAGPIPE`. Other workflow parameter readers can exist for persistence compatibility without implying production headless execution.

### `metapilot-windows-gui`

Windows desktop application. GUI code should handle UI state, user interaction, JavaFX panels, and adapting GUI events to core workflow APIs.

Typical ownership:

- JavaFX application entry point and panels.
- GUI-specific validation, dialogs, progress display, and user settings.
- Windows desktop packaging resources.

### `metapilot-web`

Spring Boot web application for server/API deployments.

Typical ownership:

- REST controllers.
- WebSocket progress broadcasting.
- Web service configuration.
- Web UI static assets.

The web module currently launches a JAR subprocess for workflow execution. Treat `metapilot-cli` as the canonical production headless runner.

## Build Commands

Build all modules:

```bash
mvn clean package -DskipTests
```

Build maintained runtime modules with dependencies:

```bash
mvn -pl metapilot-core,metapilot-cli,metapilot-windows-gui -am -DskipTests package
```

Build one module with dependencies:

```bash
mvn -pl metapilot-cli -am -DskipTests package
mvn -pl metapilot-windows-gui -am -DskipTests package
mvn -pl metapilot-web -am -DskipTests package
```

For release-facing CLI/HPC validation, prefer:

```bash
mvn -pl metapilot-core,metapilot-cli -am -DskipTests package
```

## Documentation

Start from [docs/README.md](../README.md). Keep release-facing documentation in `docs/user`, `docs/deployment`, and `docs/developer`. Keep dated progress reports and one-off planning notes in `docs/archive`.

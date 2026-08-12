# MetaPilot Cleanup Plan

## Scope

This cleanup focuses on the maintained product line:

- `metapilot-core`: shared workflow, data, search, quantification, and execution code.
- `metapilot-cli`: headless Linux/HPC entry point.
- `metapilot-windows-gui`: Windows desktop UI and GUI-specific orchestration.

The `uk.ac.ebi.pride` source tree is treated as external vendored code and should not be refactored as part of MetaPilot cleanup.

Current version target: `1.1.1`.

## Design Targets

- Keep workflow/business logic in `metapilot-core`.
- Keep CLI code focused on argument parsing, site configuration, and invoking core runners.
- Keep Windows GUI code focused on UI state, user interaction, and adapting GUI events to core workflows.
- Avoid duplicated source files across modules.
- Avoid wrapper classes that only forward calls without adding a clear boundary or policy.
- Prefer clear names and small, direct methods over defensive indirection.
- Add comments only where they explain domain decisions, non-obvious constraints, or integration contracts.

## Current Findings

- `metapilot-windows-gui` contained byte-identical copies of core quantification/resource classes. The known duplicates listed in [Duplicate Classes Removed](duplicate-classes-removed.md) have been removed so GUI code resolves them from `metapilot-core`.
- Some resource/configuration classes still have GUI-specific and core/headless variants. Consolidate them only after confirming semantic differences and migration impact.
- Generated files and crash logs should remain ignored and out of source control.
- The root Maven build includes `metapilot-web`; for CLI/HPC validation use `-pl metapilot-core,metapilot-cli -am`, and for desktop validation use `-pl metapilot-core,metapilot-cli,metapilot-windows-gui -am`.
- Several large GUI task classes still combine UI, workflow setup, process execution, and result handling. Move reusable behavior into core services/runners only with focused tests around the workflow boundary.
- Headless production execution currently covers `DIA_MARKER` and `DDA_MAG_FRAGPIPE`; do not document wider CLI support until the corresponding core runner exists and is wired through `HpcStudyRunner`.

See [Duplicate Classes Removed](duplicate-classes-removed.md) for the concrete duplicate GUI/core classes already removed.

## Refactor Sequence

1. Remove exact duplicate code and generated artifacts.
2. Consolidate remaining split-package classes after confirming semantic differences.
3. Move reusable workflow execution logic from GUI task classes into core services/runners.
4. Keep CLI and GUI as thin adapters over the same core workflow APIs.
5. Reduce dependency duplication in module POM files.
6. Add focused tests around core workflow configuration and headless runners before behavior-changing refactors.

## Validation Command

Use this command after each cleanup step:

```powershell
mvn -pl metapilot-core,metapilot-cli,metapilot-windows-gui -am -DskipTests package
```

Use this command for the current CLI/headless surface:

```powershell
mvn -pl metapilot-core,metapilot-cli -am -DskipTests package
```

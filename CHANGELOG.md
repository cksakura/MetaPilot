# Changelog

All notable user-facing changes to MetaPilot are documented here.

## 1.1.1 - 2026-08-12

### Fixed

- Fixed the FragPipe DDA open-search workflow so MAG FragPipe open search runs correctly.

## 1.1.0 - 2026-07-23

### Added

- Representative-peptide clustering for high-throughput DIA deep-search workflows, with configurable presets and similarity thresholds to reduce redundant genome searches.
- A workflow-aware result-folder cleaner with preview, Basic, and Deep modes. The cleaner preserves resumability checkpoints and final deliverables according to the selected workflow.
- Expanded CLI and HPC support for DIA Marker and MAG workflows, including additional post-analysis and database-preparation utilities.

### Improved

- Progress reporting now follows completed workflow work more closely across long-running GUI and headless tasks.
- DIA and MAG result layouts, restart behavior, parameter handling, and post-analysis exports are more consistent.
- Report generation and peptide/PSM parsing are more robust when inputs are incomplete or differ between search engines.
- User-facing validation, status messages, and workflow configuration controls have been refined.

### Fixed

- Fixed progress bars that could advance independently of the real workflow state.
- Fixed resume and cleanup edge cases that could remove files still needed by later post-analysis steps.
- Fixed several DIA-NN, FragPipe, pFind, and report-generation compatibility issues.

## 1.0.0

- First formal MetaPilot release.

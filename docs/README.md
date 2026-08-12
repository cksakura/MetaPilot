# MetaPilot Documentation

This folder contains the current release-facing documentation for MetaPilot, plus an archive of historical implementation notes.

Current version target: `1.1.1`.

## User Guides

- [Formal Release Manual Draft](MetaPilot_User_Manual_formal_release_content.md) - slide-ready structure and content for the upgraded user manual, including GUI, CLI, and HPC usage.
- [User Manual](user/user-manual.md) - current applications, workflow families, resource setup, MAG databases, outputs, and troubleshooting.
- [Workflow Design](user/workflows.md) - DDA, DIA Marker, DIA Seed, Deep, Focus, Comprehensive, MAG selection, taxonomy, and functional export behavior.
- [CLI Reference](user/cli-reference.md) - current command-line surface, headless support matrix, parameter files, resource configuration, return codes, and limits.

## Deployment

- [Linux Installation](deployment/linux-installation.md) - current Linux/headless setup using `metapilot-cli`.
- [HPC Deployment](deployment/hpc-deployment.md) - scheduler-oriented deployment for supported headless workflows.

## Developer Notes

- [Architecture](developer/architecture.md) - current Maven module layout and where code belongs.
- [Cleanup Plan](developer/cleanup.md) - cleanup scope, design targets, known findings, and validation command.
- [Duplicate Classes Removed](developer/duplicate-classes-removed.md) - record of duplicate GUI/core classes removed during cleanup.
- [Function Definition Resources](developer/function-definition-resources.md) - how to rebuild `func_def.db` and current bundled database counts.
- [Result Folder Cleaning](developer/result-folder-cleaning.md) - workflow-specific Basic/Deep clean contract, resumability checkpoints, and final deliverables.
- [WebSocket Progress](developer/websocket-progress.md) - current web progress/API status and remaining web hardening work.

## Archive

The [archive](archive/) directory keeps dated implementation summaries and generated planning prompts. These files are useful for project history but should not be treated as current release instructions.

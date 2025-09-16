# MetaPilot

MetaPilot is a comprehensive metaproteomics data analysis platform that unifies multiple search engines and workflows behind a modern JavaFX GUI and a headless CLI for HPC. The main entry point is `uk.ac.quadram.metapilot.MetaPilotApp`, which lets users select and launch workflows tailored to DDA and DIA data.

### Key capabilities
- **Multi-workflow support**: DDA (IGC and MAG) and DIA workflows
- **Multi-engine integration**: MaxQuant, pFind, FragPipe (MSFragger), Sage, DIA-NN
- **Resource management**: Centralized tool/database configuration via `Resources/` and user preferences
- **MAG database management**: Configure MGnify/local MAG catalogs using `MAGDB_1_0.json`
- **GUI and CLI**: JavaFX desktop app and `--cli` mode for HPC/headless runs

## Repository structure

```
MetaPilot/
├── src/
│   ├── main/
│   │   ├── java/uk/ac/quadram/metapilot/         # App sources (GUI, CLI, tasks, UI)
│   │   └── resources/                            # App resources (CSS, icons)
│   └── test/java/                                # JUnit tests
├── config/                                       # Last-used parameters and resource config
│   ├── last_dda_igc_maxquant_parameters.json
│   ├── last_dda_igc_pfind_parameters.json
│   ├── last_used_mag_*_parameters.json
│   ├── last_used_dia*_parameters.json
│   └── meta_resources.params                      # Saved tool paths
├── Resources/                                     # Default base path for external tools/databases
├── docs/                                          # Design docs and guides
├── MAGDB_1_0.json                                 # MAG database configuration
├── pom.xml                                        # Maven build
└── README.md
```

## Architecture overview

- **Entry point**: `MetaPilotApp` (JavaFX `Application`) builds the main window with tabs for Main (workflows), Resources, and Database Configuration
- **Workflow launch**: `MetaPilotApp.WorkflowType` drives workflow selection and opens the respective panels
- **Resources**: `MetaResourceParameters` and `MetaResourceIO` manage tool paths and availability (persisted via Java Preferences or `config/meta_resources.params`)
- **MAG DB**: `MagDbConfigIO` parses/writes `MAGDB_1_0.json`; `MagDbConfigPanelJavaFX` edits entries
- **Tasks**: DDA/DIA panels orchestrate task parameter IO and execution; CPU work runs in background threads
- **Update/registration**: `UpdateChecker` for update notices; registration logic is present but disabled by default

## Supported workflows

- **DDA IGC**
  - MaxQuant (`DDA_IGC_MAXQUANT`)
  - pFind (`DDA_IGC_PFIND`)
  - FragPipe (`DDA_IGC_FRAGPIPE`)
  - Sage (`DDA_IGC_SAGE`)
- **DDA MAG**
  - pFind (`DDA_MAG_PFIND`)
  - FragPipe (`DDA_MAG_FRAGPIPE`)
  - Sage (`DDA_MAG_SAGE`)
- **DIA**
  - Normal (`DIA_NORMAL`)
  - Target (`DIA_TARGET`)
  - Real-time (`DIA_REALTIME`)

## Requirements

- Java 17+
- Maven 3.8+
- JavaFX 21 runtime available at execution time
- External tools as needed by workflows (configure paths under Resources):
  - MaxQuant, pFind, FragPipe, Sage, DIA-NN, msconvert, FlashLFQ, FunctionDB, Python

Notes
- Default tool paths target Windows-style binaries under `Resources/` (e.g., `MaxQuantCmd.exe`, `DiaNN.exe`). Configure actual paths for your OS in the GUI Resources tab or by editing `config/meta_resources.params`.
- On Linux/macOS, ensure the corresponding executables exist, are executable, and update paths accordingly.

## Build

Build a shaded JAR:
```bash
mvn -q -DskipTests package
```
The shaded artifact will be created under `target/` with the main class `uk.ac.quadram.metapilot.MetaPilotApp`.

Run tests:
```bash
mvn -q test
```

## Run (GUI)

Easiest via Maven JavaFX plugin:
```bash
mvn -q javafx:run
```

Run the shaded JAR directly (requires JavaFX modules on the module-path):
```bash
# Set PATH_TO_FX to your JavaFX SDK lib directory
java --module-path "$PATH_TO_FX" \
     --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base \
     -jar target/*-shaded.jar
```

## Run (CLI / HPC)

MetaPilot supports a headless mode via `--cli`, currently wired to the DDA IGC MaxQuant workflow. Defaults read from `config/`.

```bash
java -jar target/*-shaded.jar --cli \
  --params config/last_dda_igc_maxquant_parameters.json \
  --resources config/meta_resources.params
```

Options
- `--params <file>`: Workflow parameter JSON
- `--resources <file>`: Tool/resource parameter file (`meta_resources.params`)
- `--help`: Show CLI help

## Configuration

- **Resources tab** (GUI): configure base resource path and individual tool executables. Backed by Java Preferences and `config/meta_resources.params` when saved via `MetaResourceIO`.
- **Database Configuration tab** (GUI): view/edit MAG catalogs; persists to `MAGDB_1_0.json` on exit.
- **Last-used parameters**: recent workflow parameters are stored under `config/` (e.g., `last_dda_igc_maxquant_parameters.json`).

## Logs

- Application log: `logs/metapilot.log` (if enabled/configured by the app)

## Documentation

- [HPC Implementation Roadmap](docs/HPC_Implementation_Roadmap.md)
- [MetaPilot HPC Migration Plan](docs/MetaPilot_HPC_Migration_Plan.md)
- [DDA IGC UI Unification Guide](docs/DDA_IGC_UI_Unification_Guide.md)
- [MGYG Ftp Task JavaFX Integration](docs/MGYGFtpTaskJavaFX_Integration.md)
- [CDF Resource Issue Resolution](docs/CDF_Resource_Issue_Resolution.md)

## Troubleshooting

- JavaFX errors when launching the JAR: ensure the JavaFX SDK is installed and passed via `--module-path` with the required `--add-modules`.
- Workflow launch disabled: open the Resources tab and configure missing tools; MetaPilot checks availability before launching workflows.
- MAG DB not loading: ensure `MAGDB_1_0.json` exists and is readable; the app will create a blank file if missing.

## License and contributing

Contribution guidelines and license information will be added. For questions or issues, please open an issue or consult the documents in `docs/`.
# MetaPilot

MetaPilot 1.2.0 is a metaproteomics workflow for DDA and DIA. On a MAG catalogue it selects the genomes that explain a sample, builds a sample-specific database, and searches that database. This guide is for installing the release and running it, with the Linux/HPC open-search workflow in full. How the identification works is described in the source repository.

The release contains:

- A Windows desktop application for interactive analysis.
- `metapilot-cli-1.2.0-shaded.jar` for Linux servers and HPC.
- `config/meta_resources.params` for the paths of external tools.

Building from source is not part of this install.

## Requirements

MetaPilot runs on Java 17 or newer. The HPC runs for this version use Java 21.

| Tool | Required for |
| --- | --- |
| FragPipe 24, with MSFragger, Crystal-C, and MSBooster under its `tools/` directory | DDA MAG FragPipe and DDA open search |
| DIA-NN 1.8.x on the cluster; the Windows GUI can use a current DIA-NN | Open-search rescoring, and every DIA workflow |
| FlashLFQ | Peptide quantification. On Linux this is the .NET build, so `dotnet` must be on `PATH` |
| Python 3 | Report pages |
| Functional annotation database | `functional_annotation/` |
| ProteoWizard `msconvert` | Thermo `.raw` input. mzML input does not need it |

Open search reads MSFragger, Crystal-C, and MSBooster from `FRAGPIPE_TOOLS`. It does not launch the FragPipe GUI.

## 1. Install the command line

Copy the jar and the `config` directory to a path that every compute node can read. The two must stay together:

```text
/opt/metapilot/metapilot-cli-1.2.0-shaded.jar
/opt/metapilot/config/meta_resources.params
```

Run from that directory, or pass `--resources` explicitly. Check Java and the tool paths:

```bash
cd /opt/metapilot
java -version
java -jar metapilot-cli-1.2.0-shaded.jar resources-check
```

Set a path that `resources-check` reports as missing:

```bash
java -jar metapilot-cli-1.2.0-shaded.jar resources-set \
  --set FRAGPIPE=/opt/fragpipe-24.0/bin/fragpipe \
  --set FRAGPIPE_TOOLS=/opt/fragpipe-24.0/tools \
  --set DIANN=/opt/fragpipe-24.0/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
  --set FLASHLFQ=/opt/metapilot/FlashLFQ/CMD \
  --set PYTHON=/opt/conda/envs/metapilot/bin/python \
  --set FUNCDB=/data/metapilot/function/func_def.db \
  --print
```

`dotnet --version` must succeed on the compute node before a run that quantifies peptides.

## 2. Prepare the MAG catalogue

Formal studies in this version use the human-gut catalogue v2.0.2. Download it once, then register it:

```bash
java -jar metapilot-cli-1.2.0-shaded.jar mag-download \
  --catalogue human-gut \
  --version v2.0.2 \
  --out /data/metapilot/MAG \
  --resume

java -jar metapilot-cli-1.2.0-shaded.jar mag-register \
  --catalogue human-gut \
  --version v2.0.2 \
  --identifier MGYG \
  --path /data/metapilot/MAG/human-gut/v2.0.2

java -jar metapilot-cli-1.2.0-shaded.jar mag-list
```

The registered folder must contain `original_db/`, `genomes-taxa.tsv`, and `Marker4.fasta`. Open search builds `Marker4_td.fasta` beside it on the first run and reuses that file afterwards. The pangenome check looks for `protein_catalogue/uhgp-100.faa.gz` under the same version. Leave `uhgp100_path` empty to use that default, or set it to the catalogue directory.

## 3. Write the parameter file

Use a parameter file saved from the Windows GUI, or start from the example shipped with the release. These fields decide the run:

```json
{
  "workflow": "DDA_MAG_OPENSEARCH",
  "input_output": {
    "mag_db_version": "v2.0.2",
    "mag_database": "human-gut",
    "append_host_db": false,
    "host_database": "/data/metapilot/host/human_reviewed.fasta",
    "output_directory": "/data/studies/demo",
    "experiments": [
      {
        "file_path": "/data/studies/demo/mzML/sample_a.mzML",
        "experiment_name": "sample_a",
        "fraction": 1
      }
    ]
  },
  "additional_parameters": {
    "identification_version": "v2",
    "genome_selection_strategy": "Comprehensive",
    "threads": 24,
    "pangenome_variant_check": "on"
  }
}
```

Paths in this file are read on the compute node. A drive letter from the Windows machine will not resolve there. `threads` should equal the CPUs you request. `append_host_db` adds the host FASTA to every sample database after genome selection; set it only when host proteins are part of the question. `Focus` selects fewer genomes than `Comprehensive`.

Check the file before submitting:

```bash
java -jar metapilot-cli-1.2.0-shaded.jar params-validate \
  --params /data/studies/demo/mag_opensearch_parameters.json
```

## 4. Submit the open search

Request memory beyond the Java heap. A 100 GB heap needs a job of about 160 GB, because MSFragger, DIA-NN, and the XGBoost library sit outside it. A study of a few dozen mzML files is a one-node job: 24 CPUs and 72 hours is a reasonable start; large cohorts need more of both.

```bash
#!/bin/bash
#SBATCH --job-name=metapilot-open
#SBATCH --cpus-per-task=24
#SBATCH --mem=160G
#SBATCH --time=72:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

export JAVA_HOME=/opt/java/jdk-21
export DOTNET_ROOT=/opt/dotnet
export PATH="$DOTNET_ROOT:$JAVA_HOME/bin:$PATH"

cd /opt/metapilot
java -Xms16g -Xmx100g -XX:+UseParallelGC \
  -jar metapilot-cli-1.2.0-shaded.jar mag-opensearch \
  --params /data/studies/demo/mag_opensearch_parameters.json
```

Submit it with `sbatch run_open.slurm.sh` from the directory where you want the Slurm logs.

The run writes `/data/studies/demo/opensearch_MAG_search/`. The first file to open is `final_summary.tsv`: one row per raw file with the MS/MS count, the number identified, the identification rate, the peptide count, and the PSM count. `combined/validated_psms.tsv` holds the accepted PSMs. `final_peptides.tsv` holds the peptide table, including the shift class. `taxonomy_analysis/Genome.tsv` and `functional_annotation/` are the genome and function reports.

## Running the same study again

Submit the same command. A step whose output is already present and still matches the current parameters is skipped: the marker search, MSBooster, genome selection, and the stage-2 searches.

Genome selection is skipped whenever `search_results/<sample>/marker/<sample>.fasta` exists. After a change to stage-1 rescoring or to the genome-selection settings, delete that FASTA. The stage-2 searches are then older than the new database and run again. A sample that previously produced no stage-1 PSMs has no FASTA and is selected again on the next run without this step.

Do not run two jobs on the same mzML files. MSFragger writes `<sample>.mzBIN_calibrated` next to each mzML, and the second job will overwrite it. Queue the second job with `--dependency=afterany:<jobid>`.

## Entrapment

To measure the error of an open search, run the entrapment command on the same parameter file. It reuses the marker search and genome selection, adds shuffled proteins to the stage-2 database, and writes `opensearch_entrapment_MAG_search/`. It stops before quantification and taxonomy. Read `combined/entrapment_summary.tsv`. This evaluation does not replace the normal open-search result.

```bash
java -Xms16g -Xmx100g -XX:+UseParallelGC \
  -jar metapilot-cli-1.2.0-shaded.jar mag-opensearch-entrapment \
  --params /data/studies/demo/mag_opensearch_parameters.json
```

## Other workflows

| Goal | Command or GUI workflow |
| --- | --- |
| Closed DDA search of the selected MAG genomes | `mag-fragpipe` |
| DIA from the catalogue Marker4 library | `dia-marker` |
| DIA from libraries you already have | `dia-seed` |
| Standard DDA against one database, without genome selection | DDA IGC in the Windows GUI |
| pFind MAG search | Windows GUI or `pfind-mag` |

`dia-marker` and `mag-fragpipe` take the same `--params` and `--resources` form as `mag-opensearch`.

## When a job fails

The Slurm `.err` file is usually empty, because the workflow log is `opensearch_MAG_search/run.log` when the launcher copies it, and otherwise `logging.log` in the result folder. The last lines of that log name the sample and the step.

| Message | What to do |
| --- | --- |
| `resources-check` reports a missing tool | Fix the path in `config/meta_resources.params`. The compute node must see the same path |
| `dotnet` is not found | Export `DOTNET_ROOT` in the job. FlashLFQ needs it |
| Marker search did not start | The mzML path in the parameter file is not visible on the node |
| Stage-1 rescoring accepted 0 targets | The file has no usable marker evidence. Later samples continue. Delete nothing; the next run selects that sample again |
| Genome selection was skipped and the genomes look stale | Delete `marker/<sample>.fasta` and submit again |
| Java was killed near the memory limit | Raise `--mem`. Keep `-Xmx` about 50 GB under the job memory |
| Two jobs searched the same mzML | Rerun them one after another so `.mzBIN_calibrated` is written by one job |

With a report, include the MetaPilot version, the parameter file, `resources-check` output, and the last 80 lines of the workflow log.

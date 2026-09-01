# HLA typing with HLA-LA — a Nextflow pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)

A [Nextflow](https://www.nextflow.io) pipeline that runs **HLA-LA** on a set of BAM/CRAM files and merges the resulting allele calls into two tidy, cohort-level tables that are ready for downstream analysis (e.g. HLA imputation reference panels, association studies).

HLA-LA is a graph-based HLA typing tool for whole-genome and exome sequencing data. See the [paper](https://academic.oup.com/bioinformatics/article/35/21/4394/5426702) (Dilthey *et al.*, 2019) and the [HLA-LA GitHub repository](https://github.com/DiltheyLab/HLA-LA).

---

## Contents

- [What the pipeline does](#what-the-pipeline-does)
- [Requirements](#requirements)
- [Installation](#installation)
- [Preparing your input](#preparing-your-input)
- [Running the pipeline](#running-the-pipeline)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Repository layout](#repository-layout)
- [Support](#support)

---

## What the pipeline does

For every input sample the pipeline:

1. **Types HLA alleles** with `HLA-LA.pl`, using the MHC population reference graph (`PRG_MHC_GRCh38_withIMGT`).
2. **Reads the G-group best-guess calls** from HLA-LA's `R1_bestguess_G.txt` (via `templates/ambigous.R`).
3. **Builds a per-sample `.ped` row** in the format expected by HLA imputation tools (FID, IID, PID, MID, SEX, PHENO, then two columns per gene).
4. **Extracts per-allele coverage** from HLA-LA's best-guess file.

Then, across all samples, it:

5. **Concatenates the `.ped` rows** into a single `.hped` file with a header.
6. **Concatenates the coverage rows** into a single `.coverage` table.

Genes reported (two alleles each): A, B, C, DQA1, DQB1, DRB1, DPA1, DPB1, DRB3, DRB4, E, F, G.

Because it is written in Nextflow, the pipeline is portable across laptops and HPC schedulers, runs samples in parallel, and can be resumed after a failure with `-resume`.

---

## Requirements

| Tool | Notes |
|------|-------|
| [Nextflow](https://www.nextflow.io) ≥ 20.10.0 | Requires Java 11+ |
| [HLA-LA](https://github.com/DiltheyLab/HLA-LA) | Plus its dependencies (bwa, samtools, picard, …) |
| HLA-LA reference graph | `PRG_MHC_GRCh38_withIMGT` (≈ 2.3 GB download; indexing needs up to ~40 GB RAM) |
| R with `dplyr` | Used when parsing HLA-LA's best-guess files |
| samtools | Only needed if you have to extract the MHC region from your BAMs (see below) |

HLA-LA itself is memory-hungry: budget roughly **70 GB of RAM** per sample for the typing step (this is the default in `nextflow.config`).

---

## Installation

### 1. Nextflow

```bash
wget -qO- https://get.nextflow.io | bash
# then move the `nextflow` binary somewhere on your PATH
```

### 2. HLA-LA

The simplest route is Bioconda:

```bash
conda install -c bioconda hla-la
```

Other installation options (Docker, from source) are described in the [HLA-LA README](https://github.com/DiltheyLab/HLA-LA#installation).

### 3. Reference graph

Download the pre-built MHC graph (≈ 2.3 GB):

```bash
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xzf PRG_MHC_GRCh38_withIMGT.tar.gz
```

Then index it. This is a one-off step; it can take **several hours** and up to **40 GB of memory**:

```bash
HLA-LA.pl --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
```

If you installed HLA-LA via conda, the graph should live in `<conda-env>/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT`. Keep note of this path — you will pass it to the pipeline as `--reference_genome`.

### 4. The pipeline

Nothing to install: Nextflow fetches the pipeline straight from GitHub the first time you run it. To update to the latest version later:

```bash
nextflow pull nanjalaruth/hla_typing_using_HLA-LA
```

---

## Preparing your input

The pipeline takes **aligned, sorted, indexed BAM (or CRAM) files**. Each file needs its index (`.bai` / `.crai`) sitting next to it.

To keep run times reasonable, HLA-LA only needs the reads mapping to the MHC region. If your BAMs are whole-genome, extract chromosome 6 from **29 Mb to 34 Mb** (GRCh38) first:

```bash
samtools index sample.sorted.bam
samtools view -b sample.sorted.bam 6:29000000-34000000 > sample.mhc.bam
samtools index sample.mhc.bam
```

> **Note:** use `chr6:29000000-34000000` if your reference uses `chr`-prefixed contig names.

Sample IDs are taken from the file name (everything before the first `.`), so `NA12878.mhc.bam` becomes sample `NA12878`.

---

## Running the pipeline

### Quick test

Runs on a small public BAM (from the nf-core `hlatyping` test dataset):

```bash
nextflow run nanjalaruth/hla_typing_using_HLA-LA \
    -profile test \
    --reference_genome /path/to/PRG_MHC_GRCh38_withIMGT \
    -resume
```

### Your own data

**Option A — pass everything on the command line:**

```bash
nextflow run nanjalaruth/hla_typing_using_HLA-LA \
    -profile slurm \
    --input "/path/to/bams/*.bam" \
    --reference_genome /path/to/PRG_MHC_GRCh38_withIMGT \
    --outdir results \
    -resume
```

> Quote the glob (`"*.bam"`) so that Nextflow, not your shell, expands it.

**Option B — use a config file.** Copy `conf/test.config`, edit the `input` (and any other) parameters to point at your data, then:

```bash
nextflow run nanjalaruth/hla_typing_using_HLA-LA \
    -profile slurm \
    -c /path/to/my.config \
    -resume
```

### Submitting to a cluster

On SLURM, wrap the command in a submission script (see `scripts/trial.sh` for an example) or launch it from a login node inside `screen`/`tmux`. The `slurm` profile submits each sample as its own job.

---

## Parameters

### Pipeline parameters (`--`)

| Parameter | Example | Description |
|-----------|---------|-------------|
| `--input` | `"/project/*.bam"` | Glob pattern matching your BAM/CRAM files. Index files must be alongside. **Required.** |
| `--reference_genome` | `/path/to/PRG_MHC_GRCh38_withIMGT` | Path to the indexed HLA-LA graph directory. **Required.** |
| `--outdir` | `./output` | Where results are written. Default: `./output`. |

### Nextflow options (`-`)

| Option | Values | Description |
|--------|--------|-------------|
| `-profile` | `standard`, `slurm`, `conda`, `debug`, `test` | Execution profile. `standard` runs locally; `slurm` submits jobs to a SLURM scheduler (for other schedulers such as PBS, add a profile in `nextflow.config`). |
| `-c` | `path/to/file.config` | Extra configuration file (overrides defaults). |
| `-resume` | — | Reuse cached results from a previous run. |
| `-r` | branch/tag | Pin a specific revision of the pipeline (default: `master`). |

Resource limits (`max_memory`, `max_cpus`, `max_time`) can be adjusted in `conf/base.config` or your own config file.

---

## Outputs

All results go under `--outdir` (default `./output`):

```
output/
├── hla_types/
│   ├── GGVP.hped        # one row per sample: FID IID PID MID SEX PHENO A.1 A.2 B.1 ... G.2 Pop
│   └── GGVP.coverage    # one row per sample: IID, then (allele, coverage) pairs for every gene
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.png
```

- **`GGVP.hped`** — the merged HLA types in "hped" format (a PLINK-style `.ped` where each pair of columns holds the two alleles of one gene).
- **`GGVP.coverage`** — per-allele read coverage reported by HLA-LA, useful for QC and for filtering low-confidence calls.

The raw HLA-LA output for each sample (including `hla/R1_bestguess_G.txt`) is kept in the HLA-LA working directory.

---

## Repository layout

```
.
├── main.nf              # pipeline definition (DSL2)
├── nextflow.config      # profiles, resources, reports
├── conf/
│   ├── base.config      # default resource limits
│   └── test.config      # settings for the quick test run
├── templates/
│   └── ambigous.R       # parses HLA-LA best-guess G-group calls
└── scripts/             # helper / legacy shell scripts (not used by the pipeline)
```

---

## Support

Questions, bugs and feature requests are tracked in the [GitHub issues](https://github.com/nanjalaruth/hla_typing_using_HLA-LA/issues).

If you use this pipeline, please also cite HLA-LA:

> Dilthey AT, et al. (2019). HLA\*LA — HLA typing from linearly projected graph alignments. *Bioinformatics*, 35(21), 4394–4396. https://doi.org/10.1093/bioinformatics/btz235

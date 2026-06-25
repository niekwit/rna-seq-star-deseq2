# Snakemake workflow: `rna-seq-star-deseq2`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/rna-seq-star-deseq2/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/rna-seq-star-deseq2/actions/workflows/main.yml)

A Snakemake workflow for `rna-seq-star-deseq2`. It will take raw RNA-seq fastq files as input, perform quality control, map the reads to the reference genome using STAR, perform differential gene expression analysis using DESeq2, and generate various plots for data visualization.

Optionally, a viral genome can be included in the analysis.

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and its DOI (see above):

Niek Wit. (2024). niekwit/rna-seq-star-deseq2: v0.6.0 (v0.6.0). Zenodo. https://doi.org/10.5281/zenodo.13693004

## Software dependencies

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [Snakemake > 8.25.5](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Apptainer (recommended)](https://apptainer.org/docs/admin/main/installation.html)

## Installation of software dependencies

Using Conda, install Snakemake in a new environment, and activate the it:

```bash
$ conda create -n star -c bioconda -c defaults snakemake=8.25.5
$ conda activate star
```

If you want to use Apptainer containers for software dependencies (recommended for GNU/Linux systems), install it via Conda into your `star` environment:

```bash
$ conda install -c conda-forge apptainer
```

## Usage

### Fetching the workflow

First clone the repository in your directory of choice using git:

```bash
$ cd /path/to/your/directory
$ git clone https://github.com/niekwit/rna-seq-star-deseq2.git
```

Create a new directory for your analysis and copy the workflow directories `config` and `workflow` to this new directory. It should look like this:

```bash
$ tree
.
├── config
│   ├── config.yaml
│   ├── README.md
│   └── samples.csv
└── workflow
    ├── envs
    │   ├── deseq2.yml
    │   ├── mapping.yml
    │   └── resources.yml
    ├── report
    │   ├── deseq2.rst
    │   ├── gprofiler2.rst
    │   ├── mapping_rates.rst
    │   ├── pca.rst
    │   ├── sample_distance.rst
    │   ├── volcano.rst
    │   └── workflow.rst
    ├── rules
    │   ├── deeptools.smk
    │   ├── deseq2.smk
    │   ├── fastqc.smk
    │   ├── mapping.smk
    │   ├── plotting.smk
    │   ├── resources.smk
    │   └── trim.smk
    ├── schemas
    │   └── config.schema.yaml
    ├── scripts
    │   ├── deseq2.R
    │   ├── general_functions.smk
    │   ├── get_readlength.sh
    │   ├── get_resource.sh
    │   ├── gprofiler2.R
    │   ├── heatmap_sd.R
    │   ├── mapping_rates.R
    │   ├── pca.R
    │   ├── __pycache__
    │   │   └── resources.cpython-312.pyc
    │   ├── resources.py
    │   └── volcano.R
    └── Snakefile

8 directories, 33 files

```

### Configuration

The analysis can be configured by editing the `config/config.yaml` file:

```yaml
# Reference genome parameters
# --------------------------------------
genome: hg38
ensembl_genome_build: 115

# Optional viral genome parameters
# --------------------------------------
# Viral genomes have to be from NCBI Datasets
#https://www.ncbi.nlm.nih.gov/datasets/genome/
viral_genome:
  apply: False # whether to include viral genome in analysis
  name: Human betaherpesvirus 5
  genome_assembly: ViralProj14559
  NCBI_RefSeq_assembly: GCF_000845245.1

# STAR (alignment) parameters
# --------------------------------------
star:
  index:
    extra: ""
  align:
    extra: ""

# DESeq2 parameters
# --------------------------------------
# Check https://htseq.readthedocs.io/en/latest/htseqcount.html for more info on stranded option
stranded: "yes" # QUOTE! unstranded, yes, reverse (as htseq-count -s values)
lfc_shrinkage: TRUE # whether to apply log fold change shrinkage in DESeq2 analysis
fdr_cutoff: 0.05 # cut off for volcano plots
fc_cutoff: 0.5 # log2 fold change cut off for volcano plots

# BigWig parameters for deepTools
# --------------------------------------
deeptools:
  normalisation: RPKM # RPKM, CPM, BPM, RPGC, None
  binsize: 10

# Computing resources
# --------------------------------------
resources:
  trim:
    cpu: 8
    time: 60
  fastqc:
    cpu: 4
    time: 20
  star_index:
    cpu: 32
    time: 60
  mapping:
    cpu: 12
    time: 60
  samtools:
    cpu: 4
    time: 30
  deeptools:
    cpu: 6
    time: 60
  deseq2:
    cpu: 6
    time: 60
  plotting:
    cpu: 2
    time: 10
```

### Preparing the sample sheet

First, prepare a directory `reads/` in your analysis directory and place all your raw fastq files there.

The sample sheet `config/samples.csv` should be edited to include the samples to be analyzed. An example is shown below:

| sample | genotype | treatment | reference |
| :----- | :------- | :-------- | :-------- |
| 5hr_1  | wt       | 5hr       | "yes"     |
| 5hr_2  | wt       | 5hr       | "yes"     |
| 72hr_1 | wt       | 72hr      | no        |
| 72hr_2 | wt       | 72hr      | no        |

Make sure the sample names correspond to the fastq files in the `reads/` directory (e.g. `5hr_1_R1_001.fastq.gz` and `5hr_1_R2_001.fastq.gz` for sample `5hr_1` for paired-end reads (`5hr_1.fastq.gz` for single-end reads)):

```bash
$ tree reads
reads
├── 5hr_1.fastq.gz
├── 5hr_2.fastq.gz
├── 72hr_1.fastq.gz
└── 72hr_2.fastq.gz

0 directories, 4 files
```

### Configuration of Snakemake

Snakemake can be configured to use Apptainer containers for software dependencies. To do so, create a yaml file (`$HOME/.config/snakemake/standard/config.yaml`):

```yaml
cores: 32
latency-wait: 20
use-conda: True
rerun-incomplete: True
printshellcmds: True
cache: False
show-failed-logs: True
use-apptainer: True
```

## Running the workflow

First, run a dry-run to check if everything is set up correctly:

```bash
$ snakemake -np
```

If everything looks good, run the workflow using:

```bash
$ snakemake --profile $HOME/.config/snakemake/standard/
```

## Expected output

After the workflow has finished, the output can be found in the `results/` directory. When run with the test dataset included with this repository, the output should look like this:

```bash
$ tree results/
results/
├── bigwig
│   ├── 5hr_1.bw
│   ├── 5hr_2.bw
│   ├── 72hr_1.bw
│   └── 72hr_2.bw
├── deseq2
│   ├── 72hr_vs_5hr.csv
│   ├── dds.RData
│   └── scale_factors.txt
├── mapped
│   ├── 5hr_1
│   │   ├── 5hr_1.Aligned.sortedByCoord.out.bam
│   │   ├── 5hr_1.Aligned.sortedByCoord.out.bam.bai
│   │   ├── 5hr_1.Log.final.out
│   │   ├── 5hr_1.Log.out
│   │   ├── 5hr_1.Log.progress.out
│   │   ├── 5hr_1.ReadsPerGene.out.tab
│   │   └── 5hr_1.SJ.out.tab
│   ├── 5hr_2
│   │   ├── 5hr_2.Aligned.sortedByCoord.out.bam
│   │   ├── 5hr_2.Aligned.sortedByCoord.out.bam.bai
│   │   ├── 5hr_2.Log.final.out
│   │   ├── 5hr_2.Log.out
│   │   ├── 5hr_2.Log.progress.out
│   │   ├── 5hr_2.ReadsPerGene.out.tab
│   │   └── 5hr_2.SJ.out.tab
│   ├── 72hr_1
│   │   ├── 72hr_1.Aligned.sortedByCoord.out.bam
│   │   ├── 72hr_1.Aligned.sortedByCoord.out.bam.bai
│   │   ├── 72hr_1.Log.final.out
│   │   ├── 72hr_1.Log.out
│   │   ├── 72hr_1.Log.progress.out
│   │   ├── 72hr_1.ReadsPerGene.out.tab
│   │   └── 72hr_1.SJ.out.tab
│   └── 72hr_2
│       ├── 72hr_2.Aligned.sortedByCoord.out.bam
│       ├── 72hr_2.Aligned.sortedByCoord.out.bam.bai
│       ├── 72hr_2.Log.final.out
│       ├── 72hr_2.Log.out
│       ├── 72hr_2.Log.progress.out
│       ├── 72hr_2.ReadsPerGene.out.tab
│       └── 72hr_2.SJ.out.tab
├── plots
│   ├── gprofiler2
│   │   ├── 72hr_vs_5hr
│   │   │   ├── downregulated.pdf
│   │   │   ├── downregulated.txt
│   │   │   ├── upregulated.pdf
│   │   │   └── upregulated.txt
│   │   ├── downregulated.pdf
│   │   ├── downregulated.txt
│   │   ├── upregulated.pdf
│   │   └── upregulated.txt
│   ├── mapping_rates.pdf
│   ├── pca.pdf
│   ├── sample_distance.pdf
│   └── volcano
│       └── 72hr_vs_5hr.pdf
├── qc
│   ├── fastqc
│   │   ├── 5hr_1_fastqc.zip
│   │   ├── 5hr_1.html
│   │   ├── 5hr_2_fastqc.zip
│   │   ├── 5hr_2.html
│   │   ├── 72hr_1_fastqc.zip
│   │   ├── 72hr_1.html
│   │   ├── 72hr_2_fastqc.zip
│   │   └── 72hr_2.html
│   ├── multiqc
│   │   ├── multiqc_data
│   │   │   ├── multiqc_citations.txt
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_fastqc.txt
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc.log
│   │   │   ├── multiqc_software_versions.txt
│   │   │   └── multiqc_sources.txt
│   │   └── multiqc.html
│   └── readlength.txt
└── trimmed
    ├── 5hr_1.fastq.gz_trimming_report.txt
    ├── 5hr_2.fastq.gz_trimming_report.txt
    ├── 72hr_1.fastq.gz_trimming_report.txt
    └── 72hr_2.fastq.gz_trimming_report.txt

16 directories, 68 files
```

## Creating the report

After the workflow has finished, a report can be generated using:

```bash
$ snakemake --report report.html
```

This will create a file `report.html` in your analysis directory containing an overview of the results.

## Integration with NGS Tracker

This workflow can be integrated with [NGS Tracker](https://github.com/niekwit/ngs-tracker).

First install requirements in your `Snakemake` environment:

```bash
$ pip install git+https://github.com/niekwit/ngs-tracker.git
```

To upgrade to the latest version:

```bash
$ pip install --upgrade git+https://github.com/niekwit/ngs-tracker.git
```

If the version number has not changed but you need to pick up a bug fix, use --force-reinstall:

```bash
$ pip install --force-reinstall git+https://github.com/niekwit/ngs-tracker.git
```

An `NGS Tracker` block is already added to the `config/config.yaml` file:

```yaml
# Integration with NGS Tracker for workflow registration and file attachment
# https://github.com/niekwit/ngs-tracker
ngs_tracker:
  enabled: false # set to false to skip registration

  # Connection
  base_url: "http://127.0.0.1:5000/api" # change host/port if needed

  # Run identity  (find project_id on the project detail page)
  project_id: 3
  # run_id: 1 # optional, if not provided a new run will be created
  workflow_name: "rna-seq-star-deseq2"
  workflow_tag: "v0.7.0"
  workflow_system: "snakemake" # snakemake | nextflow | cwl | other
  description: "Testing REST API integration with NGS Tracker"
  tags:
    - Test

  # Files to attach — paths relative to the working directory or absolute
  # type: config | sample_info | qc | results | mapping_rates | other
  files:
    - path: "config/config.yml"
      type: config
      description: "Workflow config"
    - path: "results/qc/multiqc.html"
      type: qc
      description: "MultiQC report"
    - path: "results/qc/mapping_rates.csv"
      type: mapping_rates
      description: "STAR mapping rates"
    - path: "results/qc/*.pdf"
      type: qc
    - path: logs/snakemake/*.log
      type: snakemake_log
```

On a successfull or failed run, the workflow will register itself and attach the specified files to the run in NGS Tracker.

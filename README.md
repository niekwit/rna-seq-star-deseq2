# Snakemake workflow: `rna-seq-star-deseq2`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.13.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/rna-seq-star-deseq2/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/rna-seq-star-deseq2/actions/workflows/main.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13693004.svg)](https://doi.org/10.5281/zenodo.13693004)


A Snakemake workflow for `rna-seq-star-deseq2`. It will 

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Software dependencies

* [Mamba](https://mamba.readthedocs.io/en/stable/installation/mamba-installation.html)
* [Snakemake > 8.13.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
* [Apptainer](https://apptainer.org/docs/admin/main/installation.html)


## Usage

```yaml
genome: hg38
ensembl_genome_build: 110
viral_genome: #https://www.ncbi.nlm.nih.gov/datasets/genome/
  apply: False # whether to include viral genome in analysis
  name: Human betaherpesvirus 5
  genome_assembly: ViralProj14559
  NCBI_RefSeq_assembly: GCF_000845245.1
star:
  index:
    extra: ""
  align:
    extra: ""
stranded: unstranded # unstranded, yes, reverse (as htseq-count -s values)
fdr_cutoff: 0.05 # cut off for volcano plots
fc_cutoff: 0.5 # log2 fold change cut off for volcano plots
deeptools:
  normalisation: RPKM # RPKM, CPM, BPM, RPGC, None
  binsize: 10
resources: # computing resources
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





FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="1eeecd680b0f76e8159465a6ec6217eb8590d47c38c8cb07a54cc72ecf477d8a"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/deeptools/bamcoverage/environment.yaml
#   prefix: /conda-envs/4a098defd6d43db31846a0f3e27479bb
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/4a098defd6d43db31846a0f3e27479bb
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/deeptools/bamcoverage/environment.yaml /conda-envs/4a098defd6d43db31846a0f3e27479bb/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/90ea9a0a5f4de478eb6aa43932774175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/90ea9a0a5f4de478eb6aa43932774175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/fastqc/environment.yaml /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/adec7d399fa9492178d9caa96f93852d
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
RUN mkdir -p /conda-envs/adec7d399fa9492178d9caa96f93852d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/samtools/index/environment.yaml /conda-envs/adec7d399fa9492178d9caa96f93852d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/trim_galore/se/environment.yaml
#   prefix: /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/trim_galore/se/environment.yaml /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5/environment.yaml

# Conda environment:
#   source: workflow/envs/deseq2.yml
#   prefix: /conda-envs/8a458680dda7e0cc2400420358dfd868
#   name: deseq2
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - conda-forge::r-tidyverse=2.0.0
#     - r-rcolorbrewer=1.1_3
#     - r-pheatmap=1.0.12
#     - bioconductor-deseq2=1.42.0
#     - bioconda::bioconductor-tximport=1.30.0
#     - r-ggrepel=0.9.5
#     - r-cowplot=1.1.3
#     - bioconductor-genomicfeatures=1.54.1
#     - r-reshape2=1.4.4
#     - conda-forge::r-gprofiler2=0.2.3
RUN mkdir -p /conda-envs/8a458680dda7e0cc2400420358dfd868
COPY workflow/envs/deseq2.yml /conda-envs/8a458680dda7e0cc2400420358dfd868/environment.yaml

# Conda environment:
#   source: workflow/envs/mapping.yml
#   prefix: /conda-envs/cf282b74abf4216d2af04d350cb5a342
#   name: mapping
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - star=2.7.11a
#     - pigz=2.8
#     - multiqc=1.19
RUN mkdir -p /conda-envs/cf282b74abf4216d2af04d350cb5a342
COPY workflow/envs/mapping.yml /conda-envs/cf282b74abf4216d2af04d350cb5a342/environment.yaml

# Conda environment:
#   source: workflow/envs/resources.yml
#   prefix: /conda-envs/57fcc392392389898086fddaf7d0c929
#   name: deseq2
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.12
#     - pigz=2.8
#     - agat=1.4.0
#     - ncbi-datasets-cli=16.27.1
#     - conda-forge::unzip=6.0
RUN mkdir -p /conda-envs/57fcc392392389898086fddaf7d0c929
COPY workflow/envs/resources.yml /conda-envs/57fcc392392389898086fddaf7d0c929/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/4a098defd6d43db31846a0f3e27479bb --file /conda-envs/4a098defd6d43db31846a0f3e27479bb/environment.yaml && \
    mamba env create --prefix /conda-envs/90ea9a0a5f4de478eb6aa43932774175 --file /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml && \
    mamba env create --prefix /conda-envs/adec7d399fa9492178d9caa96f93852d --file /conda-envs/adec7d399fa9492178d9caa96f93852d/environment.yaml && \
    mamba env create --prefix /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5 --file /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5/environment.yaml && \
    mamba env create --prefix /conda-envs/8a458680dda7e0cc2400420358dfd868 --file /conda-envs/8a458680dda7e0cc2400420358dfd868/environment.yaml && \
    mamba env create --prefix /conda-envs/cf282b74abf4216d2af04d350cb5a342 --file /conda-envs/cf282b74abf4216d2af04d350cb5a342/environment.yaml && \
    mamba env create --prefix /conda-envs/57fcc392392389898086fddaf7d0c929 --file /conda-envs/57fcc392392389898086fddaf7d0c929/environment.yaml && \
    mamba clean --all -y

FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="8dd37d1f99794c69d4501ac75503657aaed0108416d34652b115828c158e5f09"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/deeptools/bamcoverage/environment.yaml
#   prefix: /conda-envs/e84a94c164a7278b8266143f989704c2
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - deeptools =3.5.5
RUN mkdir -p /conda-envs/e84a94c164a7278b8266143f989704c2
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/deeptools/bamcoverage/environment.yaml /conda-envs/e84a94c164a7278b8266143f989704c2/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/fastqc/environment.yaml /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/42a2caff52fc46f2f25454995edb10df
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
RUN mkdir -p /conda-envs/42a2caff52fc46f2f25454995edb10df
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/samtools/index/environment.yaml /conda-envs/42a2caff52fc46f2f25454995edb10df/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/trim_galore/se/environment.yaml
#   prefix: /conda-envs/975a4a78cc02737b050b54b2645b9ecc
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/975a4a78cc02737b050b54b2645b9ecc
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.2.0/bio/trim_galore/se/environment.yaml /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml

# Conda environment:
#   source: workflow/envs/deseq2.yml
#   prefix: /conda-envs/07cbb24dfaf7174bd85a55152fd5c18e
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
#     - conda-forge::r-gprofiler2=0.2.3
#     - bioconda::bioconductor-limma
RUN mkdir -p /conda-envs/07cbb24dfaf7174bd85a55152fd5c18e
COPY workflow/envs/deseq2.yml /conda-envs/07cbb24dfaf7174bd85a55152fd5c18e/environment.yaml

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

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/e84a94c164a7278b8266143f989704c2 --file /conda-envs/e84a94c164a7278b8266143f989704c2/environment.yaml && \
    conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/42a2caff52fc46f2f25454995edb10df --file /conda-envs/42a2caff52fc46f2f25454995edb10df/environment.yaml && \
    conda env create --prefix /conda-envs/975a4a78cc02737b050b54b2645b9ecc --file /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml && \
    conda env create --prefix /conda-envs/07cbb24dfaf7174bd85a55152fd5c18e --file /conda-envs/07cbb24dfaf7174bd85a55152fd5c18e/environment.yaml && \
    conda env create --prefix /conda-envs/cf282b74abf4216d2af04d350cb5a342 --file /conda-envs/cf282b74abf4216d2af04d350cb5a342/environment.yaml && \
    conda env create --prefix /conda-envs/57fcc392392389898086fddaf7d0c929 --file /conda-envs/57fcc392392389898086fddaf7d0c929/environment.yaml && \
    conda clean --all -y
rule get_readlength:  
    input:
        r="results/qc/multiqc/multiqc.html",
        d="results/qc/multiqc/",
        t="results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",
    output:
        "results/qc/readlength.txt",
    conda:
        "../envs/mapping.yml"
    log:
        "logs/mapping/readlength.log"
    script:
        "../scripts/get_readlength.sh"

rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    cache: False
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_genome_fasta.log"
    threads: 1
    resources: 
        runtime=15
    conda:
        "../envs/resources.yml"
    script:
        "../scripts/get_resource.sh"


use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


if config["viral_genome"]["apply"]:
    refseq = config["viral_genome"]["NCBI_RefSeq_assembly"]
    genome_assembly = config["viral_genome"]["genome_assembly"]
    rule get_viral_resources:
        output:
            zip="resources/viral_resources.zip",
        params:
            refseq=refseq,
        log:
            "logs/resources/get_viral_resources.log"
        threads: 1
        resources: 
            runtime=15
        conda:
            "../envs/resources.yml"
        shell:
            "datasets download genome "
            "accession {params.refseq} "
            "--include gff3,rna,cds,protein,genome,seq-report "
            "--filename {output.zip} "
            " > {log} 2>&1"


    rule unpack_resources:
        input:
            zip="resources/viral_resources.zip",
        output:
            fasta=f"resources/{refseq}/ncbi_dataset/data/{refseq}/{refseq}_{genome_assembly}_genomic.fna",
            gff=f"resources/{refseq}/ncbi_dataset/data/{refseq}/genomic.gff",
        log:
            "logs/resources/unpack_resources.log"
        threads: 1
        resources: 
            runtime=5
        conda:
            "../envs/resources.yml"
        shell:
            f"unzip {{input.zip}} -d resources/{refseq} 2> {{log}}"


    rule convert_gff_to_gtf:
        input:
            gff=f"resources/{refseq}/ncbi_dataset/data/{refseq}/genomic.gff",
        output:
            gtf=f"resources/{refseq}/ncbi_dataset/data/{refseq}/genomic.gtf",
        log:
            "logs/resources/convert_gff_to_gtf.log"
        threads: 1
        resources: 
            runtime=5
        conda:
            "../envs/resources.yml"
        shell:
            "agat_convert_sp_gff2gtf.pl --gff {input} -o {output} {log};"
            "mv genomic.agat.log {log}"


    rule merge_fasta_files:
        input:
            host=resources.fasta,
            virus=f"resources/{refseq}/ncbi_dataset/data/{refseq}/{refseq}_{genome_assembly}_genomic.fna",
        output:
            resources.fasta_combined,
        log:
            "logs/resources/merge_fasta_files.log"
        threads: 1
        resources: 
            runtime=5
        conda:
            "../envs/resources.yml"
        shell:
            "cat {input.host} {input.virus} > {output}"


    rule tidy_viral_gtf:
        """
        Remove comments from GTF file (some of them cause trouble later on)
        Rename gene to gene_name for better annotation after DESeq2 analysis
        """
        input:
            gtf=f"resources/{refseq}/ncbi_dataset/data/{refseq}/genomic.gtf",
        output:
            gtf=f"resources/{refseq}/ncbi_dataset/data/{refseq}/genomic_tidy.gtf",
        log:
            "logs/resources/tidy_viral_gtf.log"
        threads: 1
        resources: 
            runtime=5
        conda:
            "../envs/resources.yml"
        shell:
            """grep -v '^#' {input} | sed 's/; gene /; gene_name /' > {output}"""


    rule merge_gtf_files:
        input:
            host=resources.gtf,
            virus=f"resources/{refseq}/ncbi_dataset/data/{refseq}/genomic_tidy.gtf",
        output:
            resources.gtf_combined,
        log:
            "logs/resources/merge_gtf_files.log"
        threads: 1
        resources: 
            runtime=5
        conda:
            "../envs/resources.yml"
        shell:
            "cat {input.host} {input.virus} > {output}"

    
    rule create_star_index:
        input:
            fasta=resources.fasta_combined,
            gtf=resources.gtf_combined,
            rl="results/qc/readlength.txt",
        output:
            directory("resources/index_star/"),
        params:
            sjdbOverhang="$(cat results/qc/readlength.txt)",
            extra=config["star"]["index"]["extra"],
        log:
            "logs/resources/create_star_index.log"
        threads: config["resources"]["star_index"]["cpu"]
        resources:
            runtime=config["resources"]["star_index"]["time"]
        conda:
            "../envs/mapping.yml"
        shell:
            "STAR "
            "--runThreadN {threads} " # Number of threads
            "--runMode genomeGenerate " # Indexation mode
            "--genomeFastaFiles {input.fasta} " # Path to fasta files
            "--sjdbOverhang {params.sjdbOverhang} " # Read-len - 1
            "--sjdbGTFfile {input.gtf} "
            "{params.extra} " # Optional parameters
            "--genomeDir {output} " # Path to output
            "> {log} 2>&1" # Logging

else:
    rule create_star_index:
        input:
            fasta=resources.fasta,
            gtf=resources.gtf,
            rl="results/qc/readlength.txt",
        output:
            directory("resources/index_star/"),
        params:
            sjdbOverhang="$(cat results/qc/readlength.txt)",
            extra=config["star"]["index"]["extra"],
        log:
            "logs/resources/create_star_index.log"
        threads: config["resources"]["star_index"]["cpu"]
        resources:
            runtime=config["resources"]["star_index"]["time"]
        conda:
            "../envs/mapping.yml"
        shell:
            "STAR "
            "--runThreadN {threads} " # Number of threads
            "--runMode genomeGenerate " # Indexation mode
            "--genomeFastaFiles {input.fasta} " # Path to fasta files
            "--sjdbOverhang {params.sjdbOverhang} " # Read-len - 1
            "--sjdbGTFfile {input.gtf} "
            "{params.extra} " # Optional parameters
            "--genomeDir {output} " # Path to output
            "> {log} 2>&1" # Logging
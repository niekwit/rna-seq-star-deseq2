if paired_end:
    rule mapping:
        input:
            fq1="results/trimmed/{sample}_val_1.fq.gz",
            fq2="results/trimmed/{sample}_val_2.fq.gz",
            idx="resources/index_star/",
        output:
            log_final="results/mapped/{sample}/{sample}.Log.final.out",
            aln="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            reads_per_gene="results/mapped/{sample}/{sample}.ReadsPerGene.out.tab",
        params:
            extra=config["star"]["align"]["extra"],
        threads: config["resources"]["mapping"]["cpu"]
        resources:
            runtime=config["resources"]["mapping"]["time"]
        log:
            "logs/mapping/{sample}.log"
        conda:
            "../envs/mapping.yml"
        shell:
            "STAR --runThreadN {threads} "
            "--runMode alignReads "
            "--genomeDir {input.idx} "
            "--readFilesIn {input.fq1} {input.fq2} "
            "--readFilesCommand gunzip -c "
            "--quantMode TranscriptomeSAM GeneCounts "
            "--outSAMtype BAM SortedByCoordinate "
            "--outTmpDir results/mapped/{wildcards.sample}/temp "
            "--outFileNamePrefix results/mapped/{wildcards.sample}/{wildcards.sample}. "
            "{params.extra} "
            "> {log} 2>&1"
else:
    rule mapping:
        input:
            fq="results/trimmed/{sample}.fq.gz",
            idx="resources/index_star/",
        output:
            log_final="results/mapped/{sample}/{sample}.Log.final.out",
            aln="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            reads_per_gene="results/mapped/{sample}/{sample}.ReadsPerGene.out.tab",
        params:
            extra=config["star"]["align"]["extra"],
        threads: config["resources"]["mapping"]["cpu"]
        resources:
            runtime=config["resources"]["mapping"]["time"]
        log:
            "logs/mapping/{sample}.log"
        conda:
            "../envs/mapping.yml"
        shell:
            "STAR --runThreadN {threads} "
            "--runMode alignReads "
            "--genomeDir {input.idx} "
            "--readFilesIn {input.fq} "
            "--readFilesCommand gunzip -c "
            "--quantMode TranscriptomeSAM GeneCounts "
            "--outSAMtype BAM SortedByCoordinate "
            "--outTmpDir results/mapped/{wildcards.sample}/temp "
            "--outFileNamePrefix results/mapped/{wildcards.sample}/{wildcards.sample}. "
            "{params.extra} "
            "> {log} 2>&1"

rule index_bam:
    input:
        bam="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools/index_{sample}.log"
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"]
    wrapper:
        f"{wrapper_version}/bio/samtools/index"
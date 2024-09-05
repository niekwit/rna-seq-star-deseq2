rule bigwig:
    input:
        bam="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
        rl="results/qc/readlength.txt",
    output:
        bw="results/bigwig/{sample}.bw",
    params:
        norm=config["deeptools"]["normalisation"],
        binsize=config["deeptools"]["binsize"],
        genome=resources.name,
        rl="$(cat results/qc/readlength.txt)",
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/bamcoverage_{sample}.log"
    wrapper:
        "v4.2.0/bio/deeptools/bamcoverage"
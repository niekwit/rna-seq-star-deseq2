rule bigwig:
    input:
        bam="results/mapped/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        idx="resources/index_star/",
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
        f"{wrapper_version}/bio/deeptools/bamcoverage"
rule deseq2:
    input:
        counts=expand("results/mapped/{sample}/{sample}.ReadsPerGene.out.tab", sample=SAMPLES),
        gtf=gtf(),
    output:
        csv=report(expand("results/deseq2/{comparison}.csv", comparison=COMPARISONS), caption="../report/deseq2.rst", category="Differential Expression Analysis"),
        rdata="results/deseq2/dds.RData",
    params:
        strand=config["stranded"],
        genome=resources.genome,
        viralgenome=config["viral_genome"]["name"],
    threads: config["resources"]["deseq2"]["cpu"]
    resources:
        runtime=config["resources"]["deseq2"]["time"]
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "../scripts/deseq2.R"
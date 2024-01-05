rule deseq2:
    input:
        counts=expand("results/mapped/{sample}/{sample}ReadsPerGene.out.tab", sample=SAMPLES),
    output:
        deseq2=report("results/deseq2/deseq2.xlsx", caption="report/deseq2.rst", category="Differential Expression Analysis of genes"),
        rdata="results/deseq2/dds.RData",
    retries: 8 # gene annotation may fail due to database failed connection
    params:
        strand=config["strand"],
        genome=resources.genome
    threads: config["resources"]["deseq2"]["cpu"]
    resources:
        runtime=config["resources"]["deseq2"]["time"]
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "../scripts/deseq2.R"
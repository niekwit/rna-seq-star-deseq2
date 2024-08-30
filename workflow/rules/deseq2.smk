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


rule gprofiler2:
    input:
        csv="results/deseq2/{comparison}.csv"
    output:
        pdf_up=report("results/gprofiler2/{comparison}/{comparison}_up.pdf", caption="../report/gprofiler2.rst", category="Gene Ontology Analysis"),
        pdf_down=report("results/gprofiler2/{comparison}/{comparison}_down.pdf", caption="../report/gprofiler2.rst", category="Gene Ontology Analysis"),
        txt_up="results/gprofiler2/{comparison}/{comparison}_up.txt",
        txt_down="results/gprofiler2/{comparison}/{comparison}_down.txt",
    params:
        genome=resources.genome,
        fdr=config["fdr_cutoff"],
        lfc=config["fc_cutoff"],
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/gprofiler2/gprofiler2_{comparison}.log"
    script:
        "../scripts/gprofiler2.R"
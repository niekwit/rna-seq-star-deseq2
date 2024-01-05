rule mapping_rates_plot:
    input:
        expand("results/mapped/{sample}/{sample}Log.final.out", sample=SAMPLES)
    output:
        report("results/plots/mapping_rates.pdf", caption="report/mapping_rates.rst", category="Mapping rates")
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/mapping_rates.log"
    script:
        "../scripts/mapping_rates.R"


rule pca_plot:
    input:
        "results/deseq2/dds.RData",
    output:
        report("results/plots/pca.pdf", caption="report/pca.rst", category="PCA"),
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/pca.log"
    script:
        "../scripts/pca.R"


rule heatmap_sample_distance:
    input:
        "results/deseq2/dds.RData",
    output:
        report("results/plots/sample_distance.pdf", caption="report/sample_distance.rst", category="Sample distances"),
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/sample_distance.log"
    script:
        "../scripts/heatmap_sd.R"


rule volcano_plot:
    input:
        "results/deseq2/deseq2_genes.xlsx",
    output:
        report(directory("results/plots/volcano_genes/"), caption="report/volcano.rst", category="Volcano plots"),
    params:
        fdr=config["fdr_cutoff"],
        fc=config["fc_cutoff"]
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/volcano.log"
    script:
        "../scripts/volcano.R"   



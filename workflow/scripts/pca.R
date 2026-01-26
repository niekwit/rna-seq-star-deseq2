# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(limma)

# Read in data
load(snakemake@input[[1]])

# Select appropriate colour palette
if (length(unique(dds$treatment)) <= 8) {
  palette <- "Dark2"
} else {
  palette <- "Set3"
}

if (length(dds) > 1000) {
  vsd <- vst(dds, blind = FALSE)
} else {
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "mean")
}

if (length(levels(dds$batch)) > 1) {
  print("Removing batch effect from data...")

  # Remove batch variation with limma
  mat <- assay(vsd)
  mm <- model.matrix(~comb, colData(vsd))
  mat <- limma::removeBatchEffect(mat, batch = vsd$batch, design = mm)
  assay(vsd) <- mat

  # Create PCA plot
  pca <- plotPCA(vsd, intgroup = c("genotype", "treatment")) +
    geom_text_repel(aes(label = vsd$sample), size = 6) +
    guides(colour = "none") +
    theme_cowplot(18) +
    scale_color_brewer(palette = palette)
} else {
  pca <- plotPCA(vsd, intgroup = c("genotype", "treatment")) +
    geom_label_repel(aes(label = vsd$sample), size = 5) +
    guides(colour = "none") +
    theme_cowplot(18) +
    scale_color_brewer(palette = palette)
}


# Save plot to file
ggsave(snakemake@output[[1]], pca, width = 10, height = 10)

# Close redirection of output/messages
sink(type = "output")
sink(type = "message")

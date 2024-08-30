# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# Read in data
load(snakemake@input[[1]])

# Log transform data
rld <- rlog(dds)

# Select appropriate colour palette
if (length(unique(rld$treatment)) <= 8) {
  palette <- "Dark2"
} else {
  palette <- "Set3"
}

# Create PCA plot
pca <- plotPCA(rld, intgroup=c("genotype", "treatment")) +
  geom_label_repel(aes(label = rld$sample),
                  size = 5) + 
  guides(colour = "none") +
  theme_cowplot(18) +
  scale_color_brewer(palette = palette)

# Save plot to file
ggsave(snakemake@output[[1]], 
       pca, 
       width=10,
       height=10)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

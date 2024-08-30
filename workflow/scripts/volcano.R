# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(cowplot)

csv <- snakemake@input[["csv"]]
pdf <- snakemake@output[["pdf"]]

# Create output directory
outdir <- dirname(csv)

dir.create(outdir, showWarnings = FALSE)

# Get plotting paramaters
fdr <- -log10(snakemake@params[["fdr"]])
fc <- snakemake@params[["fc"]]

# Load data for contrast
df <- read.csv(csv)

# Get contrast name
contrast <- unique(df$contrast_name)
print(paste0("Generating volcano plot for ", contrast, "..."))

# Log transform padj
df$log.padj <- -log10(df$padj)

# Remove genes with NA in log.padj
df <- df[!(is.na(df$log.padj)), ]

# Add fill colour based on log2FC and pvalue
df <- df %>%
  mutate(colour = case_when(
    log2FoldChange > fc & log.padj > fdr ~ "red",
    log2FoldChange < -fc & log.padj > fdr ~ "navy",
    log2FoldChange < fc & log.padj < fdr ~ "grey40",
    log2FoldChange > -fc & log.padj < fdr ~ "grey40",
  )) 

# Select top 5 down and up regulated for labels
df.up <- df %>%
  filter(log2FoldChange > fc) %>%
  filter(log.padj > fdr)
if (nrow(df.up) > 5){
  df.up <- df.up[1:5,]
}
  
df.down <- df %>%
  filter(log2FoldChange < -fc) %>%
  filter(log.padj > fdr)
  if (nrow(df.down) > 5){
    df.down <- df.down[1:5,]
  }

df.label <- rbind(df.up, df.down)

# Create plot
p <- ggplot(df, aes(x = `log2FoldChange`,
                y = `log.padj`)
        ) + 
  theme_cowplot(18) +
  theme(axis.text=element_text(size = 16),
        axis.title=element_text(size = 16),
        plot.title = element_text(hjust = 0.5,
                                  size = 16),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 16),
        panel.spacing.x = unit(2, "lines")) +
  geom_point(alpha = 0.5,
              shape = 21,
              size = 5,
              colour = "black",
              fill = df$colour) +
  ylab("-log10(adj. p value)") +
  geom_vline(xintercept = c(-fc,fc),
              linetype = "dashed", 
              color = "red", 
              linewidth = 0.5) +
  geom_hline(yintercept = fdr,
              linetype = "dashed", 
              color = "red", 
              linewidth = 0.5) +
  ggtitle(contrast) 

p <- p + geom_label_repel(size = 5,
                    aes(x = `log2FoldChange`,
                        y = `log.padj`,
                        label = `external_gene_name`), 
                    data = df.label,
                    nudge_x = -0.125,
                    nudge_y = 0.05) +
  scale_fill_manual(values = df$`colour`)
  
  # Save plot to file
  ggsave(pdf, p)

sink(log, type = "output")
sink(log, type = "message")
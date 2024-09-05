# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(gprofiler2)
library(tidyverse)

fdr <- snakemake@params[["fdr"]]
lfc <- snakemake@params[["lfc"]]
genome <- snakemake@params[["genome"]]

if (genome == "test") {
  genome <- "hg38"
  test <- TRUE
  fdr <- 1
} else {
  test <- FALSE
}

# Load DESeq2 data
csv <- read.csv(snakemake@input[["csv"]])

genes.up <- csv %>%
  filter(log2FoldChange > lfc,
         padj < fdr,
         str_detect(ensembl_gene_id, "ENS")) %>%
  pull(ensembl_gene_id)
genes.up <- list("Upregulated genes" = genes.up)

genes.down <- csv %>%
  filter(log2FoldChange < -lfc,
         padj < fdr,
         str_detect(ensembl_gene_id, "ENS")) %>%
  pull(ensembl_gene_id) 
genes.down <- list("Downregulated genes" = genes.down)

# Create organism variable
if (grepl("hg", genome, fixed = TRUE)) {
  organism <- "hsapiens"
} else if (grepl("mm", genome, fixed = TRUE)) {
  organism <- "mmusculus"
} else {
  stop("Unknown genome")
}

# Function to run gprofiler
gprofiler <- function(genes, pdf, txt) {
  gostres <- gost(query = genes, 
                organism = organism, 
                ordered_query = FALSE, 
                multi_query = TRUE, 
                significant = ifelse(test == TRUE, FALSE, TRUE), 
                exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, 
                evcodes = FALSE, 
                user_threshold = 0.05, 
                correction_method = "g_SCS", 
                domain_scope = "annotated", 
                custom_bg = NULL, 
                numeric_ns = "", 
                sources = c("GO:BP", "GO:MF", "REAC", "KEGG"), 
                as_short_link = FALSE, 
                highlight = TRUE)
  df <- as.data.frame(gostres$result) %>%
    select(-parents) #misbehaving column
  df <- as.data.frame(lapply(df, function(x) unlist(x)))
  write.table(df, txt, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Select term to highlight
  if (nrow(df) < 10) {
    terms_to_highlight <- df$term_id
  } else {
    terms_to_highlight <- df[1:10,]$term_id
  }
  
  # Create and save plot
  p <- gostplot(gostres, 
              capped = TRUE, 
              interactive = FALSE) 
  #  theme(axis.title = element_text(size = 18),
  #        text = element_text(size = 18),
  #        axis.text.x = element_text(size = 18),
  #        axis.text.y = element_text(size = 18),
  #        plot.title = element_text(size = 18)
  #        ) +
  #  scale_size_continuous(range = c(2.5,7.5))
  
  publish_gostplot(p, 
                   highlight_terms = terms_to_highlight,
                   width = 10,
                   height = 8, 
                   filename = pdf )
  
  #file.remove("Rplots.pdf")
}

gprofiler(genes.up, snakemake@output[["pdf_up"]], snakemake@output[["txt_up"]])
gprofiler(genes.down, snakemake@output[["pdf_down"]], snakemake@output[["txt_down"]])

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
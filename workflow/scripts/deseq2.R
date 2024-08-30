# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load packages
library(DESeq2)
library(GenomicFeatures)
library(tximport)
library(dplyr)
library(stringr)
library(rtracklayer)
library(readr)

# Load Snakemake variables
count.files <- snakemake@input[["counts"]]
genome <- snakemake@params[["genome"]]
gtf <- snakemake@input[["gtf"]]
strand <- snakemake@params[["strand"]]

# Load experiment information
samples <- read.csv("config/samples.csv", header=TRUE)
genotypes <- unique(samples$genotype)
treatments <- unique(samples$treatment)
design <- snakemake@params[["design"]]

if (length(genotypes) > 1 & length(treatments) > 1) {
  samples$comb <- paste0(samples$genotype, "_", samples$treatment)
} else if (length(genotypes) > 1 & length(treatments) == 1) {
  samples$comb <- samples$genotype
} else if (length(genotypes) == 1 & length(treatments) > 1) {
  samples$comb <- samples$treatment
} else {
  stop("Error: Not enough genotypes or treatments for differential expression analysis")
}

# Check if batch column exists
if ("batch" %in% colnames(samples)) {
  batches <- unique(samples$batch)
  samples$batch <- as.factor(samples$batch)
} else {
  # Add batch column with just one batch
  samples$batch <- as.factor("1")
  batches <- 1
}

# Create txdb from GTF
txdb <- makeTxDbFromGFF(gtf)

# Create transcript to gene file
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Gene annotation info
db <- rtracklayer::import(gtf)
gene.info <- data.frame(ensembl_gene_id = db$gene_id, 
                        external_gene_name = db$gene_name) %>%
  filter(!duplicated(ensembl_gene_id))

# Create count matrix
countMatrix <- read.delim(count.files[1], 
                          header = FALSE,
                          skip = 4) %>%
  dplyr::select(V1)
names(countMatrix) <- "index"

for (i in seq_along(count.files)) {
  sample <- sub(".ReadsPerGene.out.tab", "", basename(count.files[i]))
  df <- read.delim(count.files[i],
                   header = FALSE,
                   skip = 4)
  if (strand == "unstranded") {
    df <- df %>% dplyr::select(V1, V2)
  } else if (strand == "yes") {
    df <- df %>% dplyr::select(V1, V3)
  } else if (strand == "reverse") {
    df <- df %>% dplyr::select(V1, V4)
  }
  colnames(df) <- c("index", sample)
  
  countMatrix <- full_join(countMatrix,
                           df,
                           by = "index")
}

# Remove lines with all 0s
countMatrix <- countMatrix[rowSums(countMatrix[,2:ncol(countMatrix)]) > 0,]

# Create named index
rownames(countMatrix) <- countMatrix$index
countMatrix$index <- NULL

# Create DESeqDataSet
if (length(batches) == 1){
  print("Not including batch factor in DESeq2 design...")
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = samples,
                                design = ~comb)
} else {
  print("Including batch factor in DESeq2 design...")
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = samples,
                                design = ~batch + comb)
}

# Save dds to file (input for other scripts)
save(dds, file = snakemake@output[["rdata"]])

# Load reference samples
references <- unique(samples[samples$reference == "yes", ]$comb)

# List to store data from each contrast
resList <- list()

# For each reference sample, perform pairwise comparisons with all the other samples
for (r in seq_along(references)){
  cat(paste0("Setting reference level: ",references[r], " (",r,"/",length(references),")\n"))

  # Copy dds
  dds_relevel <- dds

  # For each reference sample, set it as reference in dds
  dds_relevel$comb <- relevel(dds$comb, ref = references[r])

  # Differential expression analysis
  dds_relevel <- DESeq(dds_relevel)

  # Get comparisons
  comparisons <- resultsNames(dds_relevel)
  comparisons <- strsplit(comparisons, " ")
  comparisons[1] <- NULL

  
  # Create df for each comparison
  for (c in seq_along(comparisons)){
    comparison <- comparisons[[c]]
    comparison <- str_replace(comparison, "comb_", "") # Get name
    print(paste0("Specifiying contrast: ",comparison, " (",c,"/",length(comparisons),")"))

    res <- results(dds_relevel, name = comparisons[[c]])
    df <- as.data.frame(res) %>%
      mutate(ensembl_gene_id = res@rownames, .before = 1)

    # Annotate df
    #df$ensembl_gene_id <- gsub("\\.[0-9]*","",df$ensembl_gene_id) # Tidy up gene IDs
    df <- left_join(df,gene.info,by = "ensembl_gene_id") %>%
      relocate(external_gene_name, .after = ensembl_gene_id)

    # Remove genes with baseMean zero
    df <- df[df$baseMean != 0, ]

    # Add normalised read counts for each sample to df
    temp <- as.data.frame(counts(dds_relevel, normalized = TRUE))
    temp$ensembl_gene_id <- row.names(temp)
    #temp$ensembl_gene_id <- gsub("\\.[0-9]*", "", temp$ensembl_gene_id) # Tidy up gene IDs
    names(temp)[1:length(dds_relevel@colData@listData$sample)] <- dds_relevel@colData@listData$sample

    df <- left_join(df, temp, by = "ensembl_gene_id")

    # Order data for padj
    df <- df[order(df$padj), ]

    # Sheet title can be max 31 characters
    # Change column with contrast name and change it to a number (counter) if too long
    if (nchar(comparison) > 31) {
      # change names to number
      comparison <- c
    }
    df <- df %>%
      mutate(contrast_name = comparison, .before = 1)

    # Create column with whether gene is viral or not
    df <- df %>% 
      mutate(class = case_when(
        grepl("^ENS", ensembl_gene_id, ignore.case = TRUE) ~ "non-viral",
        TRUE ~ "viral"
      )) %>%
      relocate(class, .after = external_gene_name)

    # Save df to resList
    resList[[(length(resList) + 1)]] <- df
  }
}

# Get contrast names from each df in list
names <- lapply(resList, function(x) unique(x$contrast_name))

# Name data frames in list
names(resList) <- names

# Write each df to separate csv file
for (i in seq(resList)){
  write_csv(resList[[i]],
            paste0("results/deseq2/", names(resList)[i], ".csv"))
}

# Close log file
sink(log, type = "output")
sink(log, type = "message")
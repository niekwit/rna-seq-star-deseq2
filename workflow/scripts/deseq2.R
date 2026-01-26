# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load packages
library(tidyverse)
library(DESeq2)
library(GenomicFeatures)
library(tximport)
library(rtracklayer)


# Load Snakemake variables
count_files <- snakemake@input[["counts"]]
genome <- snakemake@params[["genome"]]
viralgenome <- snakemake@params[["viralgenome"]]
gtf <- snakemake@input[["gtf"]]
strand <- snakemake@params[["strand"]]
lfc_shrinkage <- as.logical(snakemake@config[["lfc_shrinkage"]])

# Load experiment information
samples <- read.csv("config/samples.csv", header = TRUE)
genotypes <- unique(samples$genotype)
treatments <- unique(samples$treatment)

# Determine the experimental grouping
if (length(genotypes) > 1 & length(treatments) > 1) {
  samples$comb <- paste0(samples$genotype, "_", samples$treatment)
} else if (length(genotypes) > 1 & length(treatments) == 1) {
  samples$comb <- samples$genotype
} else if (length(genotypes) == 1 & length(treatments) > 1) {
  samples$comb <- samples$treatment
} else {
  stop(
    "Not enough genotypes or treatments for differential expression analysis"
  )
}

# Handle Batch Factor
if ("batch" %in% colnames(samples)) {
  batches <- unique(samples$batch)
  samples$batch <- as.factor(samples$batch)
} else {
  # Add batch column with just one batch
  samples$batch <- as.factor("1")
  batches <- 1
}

# Create txdb from GTF
print("Creating TxDb from GTF...")
txdb <- makeTxDbFromGFF(gtf)

# Create transcript to gene file
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Gene annotation info
db <- rtracklayer::import(gtf)
gene.info <- data.frame(
  ensembl_gene_id = db$gene_id,
  external_gene_name = db$gene_name
) %>%
  filter(!duplicated(ensembl_gene_id))

# Create count matrix
print("Creating count matrix...")
countMatrix <- read.delim(count_files[1], header = FALSE, skip = 4) %>%
  dplyr::select(V1)
names(countMatrix) <- "index"

for (i in seq_along(count_files)) {
  sample <- sub(".ReadsPerGene.out.tab", "", basename(count_files[i]))
  df <- read.delim(count_files[i], header = FALSE, skip = 4)
  if (strand == "unstranded") {
    df <- df %>% dplyr::select(V1, V2)
  } else if (strand == "yes") {
    df <- df %>% dplyr::select(V1, V3)
  } else if (strand == "reverse") {
    df <- df %>% dplyr::select(V1, V4)
  }
  colnames(df) <- c("index", sample)

  countMatrix <- full_join(countMatrix, df, by = "index")
}

# Remove lines with all 0s
print("Removing genes with zero counts across all samples...")
countMatrix <- countMatrix[rowSums(countMatrix[, 2:ncol(countMatrix)]) > 0, ]

# Create named index
rownames(countMatrix) <- countMatrix$index
countMatrix$index <- NULL

# Create DESeqDataSet
print("Creating DESeqDataSet...")
design_formula <- if (length(batches) > 1) ~ batch + comb else ~comb
dds <- DESeqDataSetFromMatrix(
  countData = countMatrix,
  colData = samples,
  design = design_formula
)

# Calculate size factors to normalize for sequencing depth
print("Estimating size factors...")
dds <- estimateSizeFactors(dds)

# Generate Batch-Corrected Data
# VST is better for visualization than raw normalized counts

if (length(dds) > 1000) {
  print("Applying variance stabilizing transformation...")
  vsd <- vst(dds, blind = FALSE)
} else {
  print("Applying variance stabilizing transformation for small datasets...")
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "mean")
}

batch_multipliers <- rep(1, ncol(dds))
names(batch_multipliers) <- colnames(dds)
if (length(batches) > 1) {
  print("Removing batch effect from data...")
  # Create a design matrix for the biological groups of interest.
  # This ensures limma preserves the 'comb' differences while removing 'batch'.
  mod <- model.matrix(~comb, data = colData(vsd))

  # Remove batch effect from the transformed values
  mat <- assay(vsd)
  mat_corrected <- limma::removeBatchEffect(
    mat,
    batch = vsd$batch,
    design = mod
  )

  # Calculate an "effective multiplier" for each sample.
  # We compare the mean intensity after correction vs before correction.
  # This provides a linear factor that approximates the batch shift.
  for (i in seq_along(colnames(dds))) {
    # Ratio of geometric means (effectively)
    batch_multipliers[i] <- mean(2^mat_corrected[, i]) / mean(2^mat[, i])
  }

  assay(vsd) <- mat_corrected
}
batch_corrected_counts <- as.data.frame(assay(vsd))

# Save dds for downstream use
print("Saving DESeqDataSet...")
save(dds, file = snakemake@output[["rdata"]])

# Export comprehensive Scaling Factors for BigWig scaling
print("Calculating scaling factors for BigWig generation...")
sf_df <- data.frame(
  sample = colnames(dds),
  depth_size_factor = sizeFactors(dds),
  batch_multiplier = batch_multipliers,
  # This is the final value you should use for --scaleFactor in deepTools
  final_bigwig_scalefactor = (1 / sizeFactors(dds)) * batch_multipliers
)
print(snakemake@output[["scale_factors"]])
write_csv(sf_df, snakemake@output[["scale_factors"]][1])


# Load reference samples
references <- unique(samples[samples$reference == "yes", ]$comb)

# List to store data from each contrast
resList <- list()

# For each reference sample, perform pairwise comparisons
# with all the other samples
print("Performing differential expression analysis...")
for (r in seq_along(references)) {
  cat(paste0(
    "Setting reference level: ",
    references[r],
    " (",
    r,
    "/",
    length(references),
    ")\n"
  ))

  dds_relevel <- dds
  dds_relevel$comb <- relevel(dds$comb, ref = references[r])
  dds_relevel <- DESeq(dds_relevel)

  # Get comparisons
  comp_names <- resultsNames(dds_relevel)
  comp_names <- comp_names[
    comp_names != "Intercept" & !grepl("batch", comp_names)
  ]

  # Create df for each comparison
  for (comp in comp_names) {
    clean_name <- str_replace(comp, "comb_", "")
    res <- results(dds_relevel, name = comp)

    if (lfc_shrinkage) {
      # Shrink effect size for better ranking of genes
      res <- lfcShrink(
        dds_relevel,
        coef = comp,
        res = res,
        type = "apeglm"
      )
    }

    # Get results
    df <- as.data.frame(res) %>%
      mutate(ensembl_gene_id = rownames(res), .before = 1) %>%
      left_join(gene.info, by = "ensembl_gene_id") %>%
      relocate(external_gene_name, .after = ensembl_gene_id) %>%
      mutate(lfc_shrinkage = lfc_shrinkage) # Add lfc_shrinkage info

    # Add normalized counts (Depth normalized only)
    norm_counts <- as.data.frame(counts(dds_relevel, normalized = TRUE))
    norm_counts$ensembl_gene_id <- rownames(norm_counts)

    # Remove columns with count that are not in the current comparison
    ref_condition <- str_split(clean_name, "_vs_")[[1]][2]
    alt_condition <- str_split(clean_name, "_vs_")[[1]][1]
    remove_cols <- dds$sample[
      !(dds$comb %in% c(ref_condition, alt_condition))
    ]
    norm_counts <- norm_counts %>%
      dplyr::select(
        -all_of(remove_cols)
      )

    # Add gene annotation
    df <- df %>%
      left_join(norm_counts, by = "ensembl_gene_id") %>%
      mutate(
        contrast_name = clean_name,
        class = case_when(
          grepl("^ENS", ensembl_gene_id, ignore.case = TRUE) ~ genome,
          TRUE ~ viralgenome
        )
      ) %>%
      relocate(class, .after = external_gene_name)

    resList[[length(resList) + 1]] <- df
  }
}

# Get contrast names from each df in list
names <- lapply(resList, function(x) unique(x$contrast_name))

# Name data frames in list
names(resList) <- names

# Write each df to separate csv file
print("Saving DESeq2 results...")
for (i in seq(resList)) {
  write_csv(resList[[i]], paste0("results/deseq2/", names(resList)[i], ".csv"))
}

#!/usr/bin/env Rscript

# Usage: Rscript pirna_edger_pairwise.R control_file.uniq.fa treatment_file.uniq.fa

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
  cat("Usage: Rscript pirna_edger_pairwise.R control_file.uniq.fa treatment_file.uniq.fa\n")
  quit(status=1)
}

control_file <- args[1]
treatment_file <- args[2]

# Load required packages
if (!require(edgeR, quietly=TRUE)) {
  install.packages("edgeR")
  library(edgeR)
}

# Install and load Biostrings (for reading FASTA)
if (!require(Biostrings, quietly=TRUE)) {
  if (!require(BiocManager, quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings", update=FALSE, ask=FALSE)
  library(Biostrings)
}

# Read uniq.fa and return a named vector (geneid=sequence, count=reads)
read_fa_counts <- function(fafile) {
  # Read FASTA with Biostrings
  fa <- readDNAStringSet(fafile)

  # Extract sequences and headers
  sequences <- as.character(fa)
  headers <- names(fa)

  # Extract read counts from header
  counts <- integer(length(headers))
  for (i in seq_along(headers)) {
    parts <- unlist(strsplit(headers[i], "-"))
    counts[i] <- as.integer(parts[2])
  }

  # Create named vector with sequences as names
  names(counts) <- sequences
  return(counts)
}

# Read two samples
control_counts <- read_fa_counts(control_file)
treat_counts <- read_fa_counts(treatment_file)

# Merge all gene IDs
all_geneids <- union(names(control_counts), names(treat_counts))
control_vec <- setNames(rep(0, length(all_geneids)), all_geneids)
treat_vec <- control_vec
control_vec[names(control_counts)] <- control_counts
treat_vec[names(treat_counts)] <- treat_counts

count_matrix <- cbind(control=control_vec, treatment=treat_vec)
rownames(count_matrix) <- all_geneids

# Group info
group <- factor(c("control", "treatment"), levels=c("control", "treatment"))

# edgeR analysis
dge <- DGEList(counts = count_matrix, group = group)
dge <- calcNormFactors(dge)

# No biological replicates: set a fixed dispersion (common practice for no-replicate data)
cat("Note: no biological replicates; using fixed dispersion 0.1 for analysis\n")
cat("Results are for reference; experimental validation is recommended\n")

# Set a reasonable dispersion value
dge$common.dispersion <- 0.1

# Exact test
et <- exactTest(dge, pair=c("control", "treatment"))
res <- topTags(et, n=Inf)$table

# Output results
# Create output filenames (based on input names)
control_name <- sub("\\.uniq\\.fa$", "", basename(control_file))
treatment_name <- sub("\\.uniq\\.fa$", "", basename(treatment_file))
output_prefix <- paste0("results/pirna_edger_pairwise_out/", control_name, "_vs_", treatment_name)

# Save count matrix and differential analysis results
result_file <- paste0(output_prefix, "_edger_results.csv")

write.csv(count_matrix, file=count_file, quote=FALSE)
write.csv(res, file=result_file, quote=FALSE)

# Print summary
cat("Analysis completed!\n")
cat("Total sequences:", nrow(count_matrix), "\n")
cat("Significant sequences (FDR < 0.05):", sum(res$FDR < 0.05), "\n")
cat("Upregulated (logFC > 0):", sum(res$FDR < 0.05 & res$logFC > 0), "\n")
cat("Downregulated (logFC < 0):", sum(res$FDR < 0.05 & res$logFC < 0), "\n")
cat("Results saved:\n")
cat("  - Count matrix:", count_file, "\n")
cat("  - DE results:", result_file, "\n")
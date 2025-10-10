### ---------------------------------------------------
### Correlation between RNA expression and methylation
### for whole placenta
### ---------------------------------------------------

# Function to clean sample/column names
# Removes leading zeros from sample IDs (e.g., "W01" -> "W1")
clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}

### --------------------------
### 1. Load gene annotation
### --------------------------

ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = F,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))

# Keep only gene entries
ann <- ann[ann$Class == "gene", ]

# Remove novel genes
ann <- ann[!grepl("novel", ann$Name), ]

# Extract gene IDs from GTF Name field
ann$Name <- sub('.*gene_id ([^;]+);.*', '\\1', ann$Name)

### --------------------------
### 2. Load RNA-seq CPM data
### --------------------------

rna <- read.table("W-CPM.txt", sep = "\t", header = TRUE)

# Clean column names to remove leading zeros
colnames(rna) <- clean_names(colnames(rna))

# For other tissues:
# rna <- read.table("LZ-CPM.txt", sep = "\t", header = TRUE)
# rna <- read.table("JZ-CPM.txt", sep = "\t", header = TRUE)

### --------------------------
### 3. Load methylation summary
### --------------------------

meth <- read.table("methylation_summary.tsv", sep = "\t", header = TRUE, row.names = 1)

# Remove "_CpG" suffix from column names
colnames(meth) <- sub("_CpG", "", colnames(meth))
colnames(meth) <- clean_names(colnames(meth))

# Add gene names as rownames
meth$Gene <- rownames(meth)
rownames(meth) <- meth$Gene

### --------------------------
### 4. Align RNA and methylation matrices
### --------------------------

# Keep only genes present in both datasets
common_genes <- intersect(rownames(rna), rownames(meth))
rna_aligned <- rna[common_genes, ]
meth_aligned <- meth[common_genes, ]

# Subset methylation columns to match RNA (exclude metadata)
meth_aligned <- meth_aligned[, c(4:81)]  # adjust column range for other populations if needed

# Ensure genes match
all(rownames(rna_aligned) == rownames(meth_aligned))  # should be TRUE

# Keep only samples present in both RNA and methylation
common_samples <- intersect(colnames(rna_aligned), colnames(meth_aligned))
rna_aligned <- rna_aligned[, common_samples]
meth_aligned <- meth_aligned[, common_samples]

# Check sample order
all(colnames(rna_aligned) == colnames(meth_aligned))  # should be TRUE

# Convert to matrices for faster computation
rna <- as.matrix(rna_aligned)
meth <- as.matrix(meth_aligned)

### --------------------------
### 5. Compute gene-wise Spearman correlations
### --------------------------

set.seed(123)

gene_names <- intersect(rownames(rna), rownames(meth))
rna_sub <- rna[gene_names, ]
meth_sub <- meth[gene_names, ]

# Data frame to store correlation results
cor_results <- data.frame(Gene = gene_names, Cor = NA, P = NA)

# Loop through genes
for (i in seq_along(gene_names)) {
  expr <- as.numeric(rna_sub[i, ])
  methy <- as.numeric(meth_sub[i, ])
  valid <- complete.cases(expr, methy)  # remove NAs
  if (sum(valid) >= 3) {  # require at least 3 valid points
    res <- suppressWarnings(cor.test(expr[valid], methy[valid], method = "spearman"))
    cor_results$Cor[i] <- res$estimate
    cor_results$P[i] <- res$p.value
  }
}

# Adjust p-values for multiple testing using FDR
cor_results$FDR <- p.adjust(cor_results$P, method = "fdr")

# Count significant correlations
observed_hits <- sum(abs(cor_results$Cor) > 0.3 & cor_results$FDR < 0.05, na.rm = TRUE)
paste("Exon significant correlations in W for ME is", observed_hits, "out of", length(rownames(rna_aligned)))

### --------------------------
### 6. Permutation test for empirical significance
### --------------------------

library(parallel)

n_perm <- 1000
n_cores <- detectCores() - 1
n_genes <- length(gene_names)

# Function to perform one permutation
permute_once <- function(iter) {
  # Shuffle RNA sample columns
  rna_perm <- rna_sub[, sample(ncol(rna_sub), replace = FALSE)]
  cor_vals <- numeric(n_genes)
  p_vals <- numeric(n_genes)

  for (i in seq_len(n_genes)) {
    expr <- as.numeric(rna_perm[i, ])
    methy <- as.numeric(meth_sub[i, ])
    valid <- complete.cases(expr, methy)
    if (sum(valid) >= 3) {
      test <- suppressWarnings(cor.test(expr[valid], methy[valid], method = "spearman"))
      cor_vals[i] <- test$estimate
      p_vals[i] <- test$p.value
    } else {
      cor_vals[i] <- NA
      p_vals[i] <- NA
    }
  }

  # FDR correction
  fdr_vals <- p.adjust(p_vals, method = "fdr")

  # Count significant correlations in this permutation
  sum(abs(cor_vals) > 0.3 & fdr_vals < 0.05, na.rm = TRUE)
}

# Run permutations in parallel
perm_hits <- unlist(mclapply(1:n_perm, permute_once, mc.cores = n_cores))

# Calculate empirical p-value
empirical_p <- (sum(perm_hits >= observed_hits) + 1) / (length(perm_hits) + 1)

# Save permutation hits for later inspection
saveRDS(perm_hits, "W-CPM-perm_hits.rds")

# Output empirical p-value
empirical_p

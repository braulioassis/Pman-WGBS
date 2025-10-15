###############################################################
# Gene-wise Methylation–Expression Association by Population
#   This script tests the relationship between gene expression
#   and DNA methylation across populations using linear models.
#   It aligns WGBS and RNA-seq data by gene and sample, fits
#   population-adjusted regressions per gene, and records
#   effect size, significance, R², and VIF.
###############################################################

library(car)

#--------------------------------------------------------------
# 2. Load and prepare gene annotation data
#--------------------------------------------------------------
# Read GTF annotation and retain only "gene" entries
ann <- read.table(
  "IsoquantAnnotation.gtf",
  sep = "\t", header = FALSE,
  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name")
)
ann <- ann[ann$Class == "gene", ]

# Remove uncharacterized/novel genes
ann <- ann[!grepl("novel", ann$Name), ]

# Extract clean gene IDs from GTF attributes
ann$Name <- sub('.*gene_id ([^;]+);.*', '\\1', ann$Name)

#--------------------------------------------------------------
# 3. Helper function for cleaning sample/gene names
#--------------------------------------------------------------
clean_names <- function(x) {
  # Remove leading zeros after the prefix (e.g., "BW001" → "BW1")
  gsub("^([A-Za-z]+)0*", "\\1", x)
}

#--------------------------------------------------------------
# 4. Load sample metadata and extract population groups
#--------------------------------------------------------------
samples <- read.csv("w.summary.CpG.csv", header = TRUE, sep = "\t")
# For LZ:
# samples <- read.csv("lz.summary.CpG.csv", header = TRUE, sep = "\t")
# For JZ:
# samples <- read.csv("jz.summary.CpG.csv", header = TRUE, sep = "\t")

bw <- samples$FileName[samples$Population == "BW"]
me <- samples$FileName[samples$Population == "ME"]

#--------------------------------------------------------------
# 5. Load RNA-seq and methylation summary tables
#--------------------------------------------------------------
rna <- read.table("W-CPM.txt", sep = "\t", header = TRUE)
# For LZ:
# rna <- read.table("LZ-CPM.txt", sep = "\t", header = TRUE)
# For JZ:
# rna <- read.table("JZ-CPM.txt", sep = "\t", header = TRUE)

colnames(rna) <- clean_names(colnames(rna))

meth <- read.table("w.meth.summary.tsv", sep = "\t", header = TRUE)
# For LZ:
# meth <- read.table("lz.meth.summary.tsv", sep = "\t", header = TRUE)
# For JZ:
# meth <- read.table("jz.meth.summary.tsv", sep = "\t", header = TRUE)
colnames(meth) <- sub("_CpG", "", colnames(meth))
colnames(meth) <- clean_names(colnames(meth))

#--------------------------------------------------------------
# 6. Map genes in methylation data to annotations
#--------------------------------------------------------------
# Match by contig and start position to assign gene names
meth$Gene <- ann$Name[meth$Contig == ann$Contig & meth$Start == ann$Start]
rownames(meth) <- meth$Gene

#--------------------------------------------------------------
# 7. Align RNA and methylation data by common genes and samples
#--------------------------------------------------------------
common_genes <- intersect(rownames(rna), rownames(meth))
rna_aligned <- rna[common_genes, ]
meth_aligned <- meth[common_genes, ]

# Restrict to columns containing methylation data
# (the exact range depends on population)
meth_aligned <- meth_aligned[, 4:81]
# For LZ:
# meth_aligned <- meth_aligned[, 4:71]
# For JZ:
# meth_aligned <- meth_aligned[, 4:66]

# Confirm alignment by gene names
stopifnot(all(rownames(rna_aligned) == rownames(meth_aligned)))

# Align sample order between datasets
common_samples <- intersect(colnames(rna_aligned), colnames(meth_aligned))
rna_aligned <- rna_aligned[, common_samples]
meth_aligned <- meth_aligned[, common_samples]
stopifnot(all(colnames(rna_aligned) == colnames(meth_aligned)))

sample_order <- colnames(meth_aligned)

#--------------------------------------------------------------
# 8. Match population labels to samples
#--------------------------------------------------------------
population <- samples$Population[match(sample_order, samples$FileName)]
population <- factor(population)  # ensure it's a factor
stopifnot(all(!is.na(population)))  # sanity check

# Optional preview of population assignment
head(data.frame(sample = sample_order, population = population), 10)

#--------------------------------------------------------------
# 9. Example: test methylation-expression association for one gene
#--------------------------------------------------------------
g <- rownames(meth_aligned)[1]  # select first gene
meth_values <- as.numeric(meth_aligned[g, ])
expr_values <- as.numeric(rna_aligned[g, ])

# Fit model: expression ~ methylation + population
fit <- lm(expr_values ~ meth_values + population)
summary(fit)

#--------------------------------------------------------------
# 10. Define function to run per-gene regression
#--------------------------------------------------------------
run_model <- function(i) {
  expr_values <- as.numeric(rna_aligned[i, ])
  meth_values <- as.numeric(meth_aligned[i, ])
  
  # Keep complete (non-missing) data
  keep <- !is.na(expr_values) & !is.na(meth_values)
  expr_sub <- expr_values[keep]
  meth_sub <- meth_values[keep]
  pop_sub <- droplevels(population[keep])
  
  # Skip if insufficient data or population diversity
  if (length(expr_sub) < 3 || nlevels(pop_sub) < 2) {
    return(c(beta = NA, pval = NA, r2 = NA, vif.meth_sub = NA))
  }
  
  # Fit model
  fit <- lm(expr_sub ~ meth_sub + pop_sub)
  s <- summary(fit)
  
  # Extract coefficient and p-value for methylation
  coef_table <- s$coefficients
  beta <- if ("meth_sub" %in% rownames(coef_table)) coef_table["meth_sub", "Estimate"] else NA
  pval <- if ("meth_sub" %in% rownames(coef_table)) coef_table["meth_sub", "Pr(>|t|)"] else NA
  
  # Extract model R²
  r2 <- s$r.squared
  
  # Calculate variance inflation factor (collinearity check)
  vif_values <- vif(fit)
  vif_meth <- if ("meth_sub" %in% names(vif_values)) vif_values["meth_sub"] else NA
  
  # Return results as a vector
  c(beta = beta, pval = pval, r2 = r2, vif.meth_sub = vif_meth)
}

#--------------------------------------------------------------
# 11. Apply model across all genes
#--------------------------------------------------------------
results <- t(sapply(1:nrow(meth_aligned), run_model))
rownames(results) <- rownames(meth_aligned)
results <- as.data.frame(results)

#--------------------------------------------------------------
# 12. Summarize multicollinearity diagnostics
#--------------------------------------------------------------
mean(results$vif.meth_sub, na.rm = TRUE)      # average VIF across genes
quantile(results$vif.meth_sub, 0.95, na.rm = TRUE)  # 95th percentile of VIF

#--------------------------------------------------------------
# 13. Multiple testing correction and output
#--------------------------------------------------------------
results$padj <- p.adjust(results$pval, method = "BH")  # FDR adjustment
write.csv(results, "W-popvif-results.csv")
# For LZ:
# write.csv(results, "LZ-popvif-results.csv")
# For JZ:
# write.csv(results, "JZ-popvif-results.csv")

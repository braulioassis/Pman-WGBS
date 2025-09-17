clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}

ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = F,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))
ann <- ann[ann$Class == "gene", ]
ann <- ann[!grepl("novel", ann$Name), ]
ann$Name <- sub('.*gene_id ([^;]+);.*', '\\1', ann$Name)

rna <- read.table("W-CPM.txt", sep = "\t", header = T)
colnames(rna) <- clean_names(colnames(rna))

# For LZ:
# rna <- read.table("LZ-CPM.txt", sep = "\t", header = T)

# For JZ:
# rna <- read.table("JZ-CPM.txt", sep = "\t", header = T)

meth <- read.table("methylation_summary.tsv", sep = "\t", header = T, row.names = 1)
# From appropriate tissue
colnames(meth) <- sub("_CpG", "", colnames(meth))
colnames(meth) <- clean_names(colnames(meth))
meth$Gene <- rownames(meth)
rownames(meth) <- meth$Gene

common_genes <- intersect(rownames(rna), rownames(meth))
rna_aligned <- rna[common_genes, ]
meth_aligned <- meth[common_genes, ]

common_genes <- intersect(rownames(rna), rownames(meth_aligned))
meth_aligned <- meth_aligned[, c(4:81)]
# For LZ, 71:
# meth_aligned <- meth_aligned[, c(4:71)]

# For JZ, 66:
# meth_aligned <- meth_aligned[, c(4:66)]

all(rownames(rna_aligned) == rownames(meth_aligned))

common_samples <- intersect(colnames(rna_aligned), colnames(meth_aligned))
rna_aligned <- rna_aligned[, common_samples]
meth_aligned <- meth_aligned[, common_samples]

all(colnames(rna_aligned) == colnames(meth_aligned))

meth <- as.matrix(meth_aligned)
rna <- as.matrix(rna_aligned)

set.seed(123)

gene_names <- intersect(rownames(rna), rownames(meth))
rna_sub <- rna[gene_names, ]
meth_sub <- meth[gene_names, ]

cor_results <- data.frame(Gene = gene_names, Cor = NA, P = NA)

for (i in seq_along(gene_names)) {
  expr <- as.numeric(rna_sub[i, ])
  methy <- as.numeric(meth_sub[i, ])

  valid <- complete.cases(expr, methy)
  if (sum(valid) >= 3) {
    res <- suppressWarnings(cor.test(expr[valid], methy[valid], method = "spearman"))
    cor_results$Cor[i] <- res$estimate
    cor_results$P[i] <- res$p.value
  }
}

cor_results$FDR <- p.adjust(cor_results$P, method = "fdr")
observed_hits <- sum(abs(cor_results$Cor) > 0.3 & cor_results$FDR < 0.05, na.rm = TRUE)
paste("Exon significant correlations in W for ME is", observed_hits, "out of", length(rownames(rna_aligned)))

library(parallel)

n_perm <- 1000
n_cores <- detectCores() - 1
gene_names <- intersect(rownames(rna), rownames(meth))
rna_sub <- rna[gene_names, ]
meth_sub <- meth[gene_names, ]
n_genes <- length(gene_names)

permute_once <- function(iter) {
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

  fdr_vals <- p.adjust(p_vals, method = "fdr")
  sum(abs(cor_vals) > 0.3 & fdr_vals < 0.05, na.rm = TRUE)
}

perm_hits <- unlist(mclapply(1:n_perm, permute_once, mc.cores = n_cores))

empirical_p <- (sum(perm_hits >= observed_hits) + 1) / (length(perm_hits) + 1)

saveRDS(perm_hits, "W-CPM-perm_hits.rds")
empirical_p

library(edgeR)
library(limma)

# For early gestation, whole placenta
counts <- read.csv("EP_Pman_ExtMMFrac_readcounts_Exon.csv")
samples <- read.csv("w.summary.CpG.csv")
colnames(counts_genes) <- sub(".*_W([^_]+)_.*", "W\\1", colnames(counts_genes))
clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}

# For late gestation, labyrinth
#counts <- read.csv("LP_Pman_ExtMMFrac_readcounts_Exon.csv")
#samples <- read.csv("lz.summary.CpG.csv")
#colnames(counts_genes) <- sub(".*_LZ([^_]+)_.*", "LZ\\1", colnames(counts_genes))
#clean_names <- function(x) {
#  gsub("^([A-Za-z]+)0*", "\\1", x)
#}

# For late gestation, junctional zone
#counts <- read.csv("LP_JZ_Pman_ExtMMFrac_readcounts_Exon.csv")
#samples <- read.csv("jz.summary.CpG.csv")
#colnames(counts_genes) <- sub(".*_JZ([^_]+)_.*", "JZ\\1", colnames(counts_genes))
#clean_names <- function(x) {
#  gsub("^([A-Za-z]+)0*", "\\1", x)
#}

gene_names <- counts[,1]
rownames(counts) <- gene_names
counts_genes <- counts[, -c(1:6)]

colnames(counts_genes) <- clean_names(colnames(counts_genes))

counts_genes <- as.matrix(na.omit(counts_genes))
samples <- samples[samples$FileName %in% colnames(counts_genes), ]
counts_genes <- counts_genes[, colnames(counts_genes) %in% samples$FileName]
y <- DGEList(counts = counts_genes, samples = samples)
y <- calcNormFactors(y)

keep <- rowSums(cpm(y) > 0.5 ) >= (0.75 * ncol(counts_genes))
y <- y[keep, , keep.lib.sizes=FALSE]
design <- model.matrix(~ Population + Treatment, data = y$samples)

# voom transform (uses whole dataset to learn mean-variance)
v <- voom(y, design, plot = TRUE)    # plot helps inspect mean-variance trend

# linear model + empirical Bayes
fit <- lmFit(v, design)
fit <- eBayes(fit)

# extract results
res <- topTable(fit, coef = "PopulationME", number = Inf)
results_subset <- res[rownames(res) %in% c("Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Tet1", "Tet2", "Tet3", "Tdg"), , drop = FALSE]

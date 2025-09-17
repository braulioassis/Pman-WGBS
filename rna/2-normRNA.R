library(edgeR)

# W
counts <- read.csv("EP_Pman_ExtMMFrac_readcounts_Exon.csv")
gene_names <- counts[,1]
rownames(counts) <- gene_names
counts_genes <- counts[, -c(1:6)]

colnames(counts_genes) <- sub(".*_W([^_]+)_.*", "W\\1", colnames(counts_genes))
clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}
colnames(counts_genes) <- clean_names(colnames(counts_genes))

counts_genes <- as.matrix(na.omit(counts_genes))
counts_genes <- DGEList(counts_genes)

weights <- calcNormFactors(counts_genes)
dim(weights)
keep <- rowSums(cpm(weights) > 0.5 ) >= (0.75 * ncol(counts_genes))
rnanorm <- weights[keep,]
dim(rnanorm)
write.table(rnanorm$counts, ("W-CPM.txt"), sep = "\t", row.names = T, quote = F)

# LZ

library(edgeR)

counts <- read.csv("LP_Pman_ExtMMFrac_readcounts_Exon.csv")
gene_names <- counts[,1]
rownames(counts) <- gene_names
counts_genes <- counts[, -c(1:6)]

colnames(counts_genes) <- sub(".*_LZ([^_]+)_.*", "JZ\\1", colnames(counts_genes))
clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}
colnames(counts_genes) <- clean_names(colnames(counts_genes))

counts_genes <- as.matrix(na.omit(counts_genes))
counts_genes <- DGEList(counts_genes)

weights <- calcNormFactors(counts_genes)
dim(weights)
keep <- rowSums(cpm(weights) > 0.5 ) >= (0.75 * ncol(counts_genes))
rnanorm <- weights[keep,]
dim(rnanorm)
write.table(rnanorm$counts, ("LZ-CPM.txt"), sep = "\t", row.names = T, quote = F)

# JZ

library(edgeR)

counts <- read.csv("LP_JZ_Pman_ExtMMFrac_readcounts_Exon.csv")
gene_names <- counts[,1]
rownames(counts) <- gene_names
counts_genes <- counts[, -c(1:6)]

colnames(counts_genes) <- sub(".*_JZ([^_]+)_.*", "JZ\\1", colnames(counts_genes))
clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}
colnames(counts_genes) <- clean_names(colnames(counts_genes))

counts_genes <- as.matrix(na.omit(counts_genes))
counts_genes <- DGEList(counts_genes)

weights <- calcNormFactors(counts_genes)
dim(weights)
keep <- rowSums(cpm(weights) > 0.5 ) >= (0.75 * ncol(counts_genes))
rnanorm <- weights[keep,]
dim(rnanorm)
write.table(rnanorm$counts, ("JZ-CPM.txt"), sep = "\t", row.names = T, quote = F)

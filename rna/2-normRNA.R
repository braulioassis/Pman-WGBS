### ---------------------------------------------------
### edgeR RNA-seq preprocessing for W, LZ, and JZ tissues
### Steps: load counts, clean column names, filter low expression genes,
### normalize counts, and save CPM-filtered data
### ---------------------------------------------------

library(edgeR)


### --------------------------
### 1. Whole placenta, early gestation (W)
### --------------------------

# Load exon read count table
counts <- read.csv("EP_Pman_ExtMMFrac_readcounts_Exon.csv")

# Extract gene names and set as rownames
gene_names <- counts[,1]
rownames(counts) <- gene_names

# Remove non-count columns
counts_genes <- counts[, -c(1:6)]

# Clean column names to standard sample IDs
# Extract sample ID patterns like "_WXXX_" → "WXXX"
colnames(counts_genes) <- sub(".*_W([^_]+)_.*", "W\\1", colnames(counts_genes))

# Function to remove leading zeros from sample IDs
clean_names <- function(x) {
  gsub("^([A-Za-z]+)0*", "\\1", x)
}
colnames(counts_genes) <- clean_names(colnames(counts_genes))

# Convert counts to matrix and remove rows with NA
counts_genes <- as.matrix(na.omit(counts_genes))

# Create DGEList object for edgeR
counts_genes <- DGEList(counts_genes)

# Calculate normalization factors to account for library size differences
weights <- calcNormFactors(counts_genes)

# Filter lowly expressed genes: keep genes expressed >0.5 CPM in ≥75% of samples
keep <- rowSums(cpm(weights) > 0.5 ) >= (0.75 * ncol(counts_genes))
rnanorm <- weights[keep,]

# Save filtered CPM counts
write.table(rnanorm$counts, "W-CPM.txt", sep = "\t", row.names = TRUE, quote = FALSE)


### --------------------------
### 2. Late gestation labyrinth
### --------------------------

# Load exon read count table
counts <- read.csv("LP_Pman_ExtMMFrac_readcounts_Exon.csv")
gene_names <- counts[,1]
rownames(counts) <- gene_names
counts_genes <- counts[, -c(1:6)]

# Clean column names: "_LZXXX_" → "LZXXX"
colnames(counts_genes) <- sub(".*_LZ([^_]+)_.*", "LZ\\1", colnames(counts_genes))
colnames(counts_genes) <- clean_names(colnames(counts_genes))

# Convert to matrix and remove NAs
counts_genes <- as.matrix(na.omit(counts_genes))
counts_genes <- DGEList(counts_genes)

# Normalize library sizes
weights <- calcNormFactors(counts_genes)

# Filter lowly expressed genes
keep <- rowSums(cpm(weights) > 0.5 ) >= (0.75 * ncol(counts_genes))
rnanorm <- weights[keep,]

# Save CPM-filtered counts
write.table(rnanorm$counts, "LZ-CPM.txt", sep = "\t", row.names = TRUE, quote = FALSE)


### --------------------------
### 3. Junctional zone
### --------------------------

# Load exon read count table
counts <- read.csv("LP_JZ_Pman_ExtMMFrac_readcounts_Exon.csv")
gene_names <- counts[,1]
rownames(counts) <- gene_names
counts_genes <- counts[, -c(1:6)]

# Clean column names: "_JZXXX_" → "JZXXX"
colnames(counts_genes) <- sub(".*_JZ([^_]+)_.*", "JZ\\1", colnames(counts_genes))
colnames(counts_genes) <- clean_names(colnames(counts_genes))

# Convert to matrix and remove NAs
counts_genes <- as.matrix(na.omit(counts_genes))
counts_genes <- DGEList(counts_genes)

# Normalize library sizes
weights <- calcNormFactors(counts_genes)

# Filter lowly expressed genes
keep <- rowSums(cpm(weights) > 0.5 ) >= (0.75 * ncol(counts_genes))
rnanorm <- weights[keep,]

# Save CPM-filtered counts
write.table(rnanorm$counts, "JZ-CPM.txt", sep = "\t", row.names = TRUE, quote = FALSE)

### ---------------------------------------------------
### Summarize CpG methylation over gene, promoter, and intron/exon regions
### ---------------------------------------------------

library(data.table)
library(GenomicRanges)
library(future.apply)
library(stringr)

# Use multicore parallel processing
plan(multicore, workers = parallel::detectCores())


### 1. Load gene annotation (GTF file)
ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = FALSE,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))

# Keep only gene-level annotations
ann <- ann[ann$Class == "gene", ]

# Remove novel/uncharacterized genes
ann <- ann[!grepl("novel", ann$Name), ]

# Extract gene_id from GTF 'Name' column
ann$Name <- sub('.*gene_id ([^;]+);.*', '\\1', ann$Name)

# Set rownames for indexing convenience
rownames(ann) <- ann$Name


### 2. Load expression data (filter to genes in annotation)
rna <- read.csv("EP_Pman_ExtMMFrac_readcounts_Exon.csv")
rna <- rna[rna$Geneid %in% ann$Name, ]

# Alternative tissues:
# rna <- read.csv("LP_Pman_ExtMMFrac_readcounts_Exon.csv")  # LZ
# rna <- read.csv("LP_JZ_Pman_ExtMMFrac_readcounts_Exon.csv")  # JZ


### 3. Convert annotation to GRanges for overlap with methylation
setDT(ann, keep.rownames = TRUE)
setnames(ann, old = c(names(ann)[2], names(ann)[5], names(ann)[6]), new = c("Contig", "Start", "End"))

# Expand each gene ±500bp to include promoters/adjacent regions
ann_expanded <- ann[, .(Contig, Start = Start - 500, End = End + 500)]

# Convert to GRanges object
ann_gr <- GRanges(seqnames = ann_expanded$Contig,
                  ranges = IRanges(start = ann_expanded$Start, end = ann_expanded$End))


### 4. Load methylation files
meth_dir <- "/CpG/"
meth_files <- list.files(path = meth_dir, pattern = "\\.meth$", full.names = TRUE)


### 5. Aggregate methylation per gene in parallel
file_results <- future_lapply(meth_files, function(file) {
  
  # Load methylation file
  meth <- fread(file, col.names = c("Contig", "Pos", "End", "Strand", "Score", "Coverage"))
  
  # Keep CpGs with coverage >= 5
  meth <- meth[Coverage >= 5]
  
  # Convert to GRanges
  meth_gr <- GRanges(seqnames = meth$Contig,
                     ranges = IRanges(start = meth$Pos, width = 1),
                     Score = meth$Score)
  
  # Find overlaps between genes and CpG sites
  hits <- findOverlaps(ann_gr, meth_gr)
  
  # Extract overlapping scores
  overlaps_dt <- data.table(queryHits = queryHits(hits),
                            Score = mcols(meth_gr)$Score[subjectHits(hits)])
  
  # Compute average methylation per gene
  avg_by_region <- overlaps_dt[, .(avg = mean(Score, na.rm = TRUE)), by = queryHits]
  
  # Initialize vector for all genes
  avg_vec <- rep(NA_real_, length(ann_gr))
  avg_vec[avg_by_region$queryHits] <- avg_by_region$avg
  
  # Set column name based on file name (remove .meth extension)
  colname <- str_remove(basename(file), "\\.meth$")
  setnames(data.table(avg_vec), colname)
})


### 6. Combine methylation data into one table
methylation_data <- do.call(cbind, file_results)

# Merge methylation with gene coordinates
final <- cbind(ann[, .(Contig, Start, End)], methylation_data)
rownames(final) <- ann$rn

# Save methylation summary
fwrite(final, "/methylation_summary.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


### 7. Promoter methylation summary
# Promoter region defined as 500bp upstream for + strand, downstream for - strand
ann_expanded <- ann[, .(
  Contig,
  Start = pmax(ifelse(Strand == "+", Start - 500, End), 1),
  End   = pmax(ifelse(Strand == "+", Start,       End + 500), 1)
)]


### 8. First intron methylation summary
# Select exon annotations
exons <- ann[Class == "exon"]

# Get first two exons per gene
first_two_exons <- exons[, .SD[1:2], by = Name]

# Split into first and second exon tables
exon1 <- first_two_exons[, .SD[1], by = Name]
exon2 <- first_two_exons[, .SD[2], by = Name]

# Preserve original gene order
gene_order <- unique(ann$Name)

# Merge exon1 and exon2 to define introns (between exons)
setkey(exon1, Name)
setkey(exon2, Name)
introns <- exon2[exon1, nomatch = 0]   # preserves exon1 order

# Define intron coordinates
introns <- introns[, .(
 Name,
 Contig,
 Start = pmin(End + 1, i.Start - 1),
 End   = pmax(End + 1, i.Start - 1)
)]

# Remove incomplete rows (genes missing 2nd exon)
vintrons <- introns[complete.cases(introns)]

# Ensure introns are in same order as ann
introns <- introns[match(gene_order, introns$Name), ]
introns <- introns[!is.na(Name)]

# Convert introns to GRanges
intron_gr <- GRanges(
 seqnames = introns$Contig,
 ranges = IRanges(start = introns$Start, end = introns$End)
)


### 9. First exon methylation summary
# Take first exon per gene
exons <- ann[Class == "exon", .SD[1], by = Name]
genenames <- exons$Name

# Convert to GRanges
exon_gr <- GRanges(
 seqnames = exons$Contig,
 ranges = IRanges(start = exons$Start, end = exons$End)
)

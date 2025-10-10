### ---------------------------------------------------
### Extract gene names overlapping DMRs from genome annotation
### ---------------------------------------------------

# 1. Load chromosome information
chrs <- read.table("sequence_report.tsv", header = TRUE, sep = "\t")

# Adjust GenBank accession numbers
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)


# 2. Load genome annotation
ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = FALSE,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))

# Keep only gene entries
ann <- ann[ann$Class == "gene", ]

# Exclude novel or LOC entries
ann <- ann[!grepl("novel", ann$Name), ]
ann <- ann[!grepl("LOC", ann$Name), ]

# Map Contig to chromosome name using the sequence report
ann$Chrom <- chrs$Sequence.name[match(ann$Contig, chrs$RefSeq.seq.accession)]

# 3. Load DMR (Differentially Methylated Region) data
dmr <- read.table("dmr.bed", sep = "\t", header = FALSE, 
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))

# Map DMR contig to chromosome name
dmr$Chrom <- chrs$Sequence.name[match(dmr$Contig, chrs$RefSeq.seq.accession)]

# Calculate gap per DMR (length divided by count)
dmr$Gap <- (dmr$End - dmr$Start) / dmr$Count

# Filter DMRs:
# - Keep DMRs with more than 29 CpGs (Count > 29)
# - Keep DMRs where average distance between CpGs is less than 50 bp (Gap < 50)
dmr <- dmr[dmr$Count > 29, ]
dmr <- dmr[dmr$Gap < 50, ]
# For hypoxia-specific methylation changes, we use dmr$Count > 7, dmr$Gap > 0


# 4. Identify genes overlapping DMRs
vec <- c()

for (i in 1:nrow(ann)) {
  y <- c()
  
  # Check if any DMR overlaps the gene ±500 bp
  y <- which(
    ((dmr$Chrom == ann$Chrom[i]) & (dmr$Start >= ann$Start[i] - 500) & (dmr$Start <= ann$End[i] + 500)) |
    ((dmr$Chrom == ann$Chrom[i]) & (dmr$End >= ann$Start[i] - 500) & (dmr$End <= ann$End[i] + 500))
  )
  
  # If overlaps exist, extract gene information
  if (length(y) >= 1) {
    x <- ann[i, c(10,9)]  # Columns Name and V8 (gene_id info)
    vec <- rbind(vec, x)
  }
}

# 5. Extract gene IDs from the GTF 'Name' column
# GTF 'gene_id' format: gene_id "XYZ"; transcript_id "ABC";
vec$extracted <- sub('.*gene_id ([^;]+);.*', '\\1', vec$Name)

# Write gene IDs overlapping DMRs to file
write(vec$extracted, "genenames.txt")

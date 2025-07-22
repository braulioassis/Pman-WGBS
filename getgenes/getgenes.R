# Get chromosome names
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)

# Get genome annotation
ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = F,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))
ann <- ann[ann$Class == "gene", ]
ann <- ann[!grepl("novel", ann$Name), ]
ann <- ann[!grepl("LOC", ann$Name), ]
ann$Chrom <- chrs$Sequence.name[match(ann$Contig, chrs$RefSeq.seq.accession)]

# Get DMR's
dmr <- read.table("dmr.bed", sep = "\t", header = F, 
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
dmr$Chrom <- chrs$Sequence.name[match(pop$Contig, chrs$RefSeq.seq.accession)]
dmr$Gap <- (dmr$End - dmr$Start)/dmr$Count
dmr <- dmr[dmr$Count > 29, ]
dmr <- dmr[dmr$Gap < 50, ]
# When filtering for within population methylation changes in hypoxia, dmr$Count > 7, dmr$Gap > 0

vec <- c()
for (i in 1:nrow(ann)) {
  y <- c()
  y <- which(
    ((dmr$Chrom == ann$Chrom[i]) & (dmr$Start >= ann$Start[i] - 500) & (dmr$Start <= ann$End[i] + 500)) |
      ((dmr$Chrom == ann$Chrom[i]) & (dmr$End >= ann$Start[i] - 500) & (dmr$End <= ann$End[i] + 500))
  )
  if (length(y) >= 1) {
    x <- ann[i, c(10,9)]
    vec <- rbind(vec,x)
  }
}

vec$extracted <- sub('.*gene_id ([^;]+);.*', '\\1', vec$Name)
write(vec$extracted, "genenames.txt")

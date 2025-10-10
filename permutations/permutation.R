library(regioneR)

# Get Chr numbers
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)

# Get chromosome sizes
geno <- chrs[, c(14, 11)]
geno <- geno[c(1:24), ]
custom.genome <- toGRanges(geno)

# Get DMRs
dmr <- read.table("dmr.bed", sep = "\t", header = F, 
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
dmr$Chrom <- chrs$Sequence.name[match(dmr$Contig, chrs$RefSeq.seq.accession)]
dmr$Gap <- (dmr$End - dmr$Start)/dmr$Count
dmr <- dmr[dmr$Count > 29, ]
dmr <- dmr[dmr$Gap < 50, ]
# When filtering for within population methylation changes in hypoxia, dmr$Count > 7, dmr$Gap > 0
dmr <- dmr[, c(7, 2, 3)]

# Get inversions
inv <- read.csv("Inversions.csv")
inv$Chr <- paste0("chr",inv$chromosome)
inv <- inv[inv$inversion %in% c("inv6.0", "inv7.0", "inv7.3", "inv15.1", "inv18.0", "inv22.0"), ]
inv <- inv[, c(5, 3, 4)]

A <- toGRanges(dmr)
B <- toGRanges(inv)
numOverlaps(A,B)

pt <- permTest(A = A, ntimes = 1000, alternative = "greater", randomize.function = function(x, ...) randomizeRegions(x, genome = custom.genome),
               evaluate.function = numOverlaps, B = B, verbose = F)
saveRDS(pt, file = "inversions-popdmr-permutation.rds")

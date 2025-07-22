library(regioneR)

# Get Chr numbers
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)

# Get chromosome sizes
genome <- chrs[, c(14, 11)]

# Get DMRs
dmr <- read.table("dmr.bed", sep = "\t", header = F, 
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
dmr$Chrom <- chrs$Sequence.name[match(dmr$Contig, chrs$RefSeq.seq.accession)]
dmr$Gap <- (dmr$End - dmr$Start)/dmr$Count
dmr <- dmr[dmr$Count > 29, ]
dmr <- dmr[dmr$Gap < 50, ]
# For ...
dmr <- dmr[, c(7, 2, 3)]

# Get inversions
inv <- read.csv("Inversions.csv")
inv$Chr <- paste0("chr",inv$chromosome)
inv <- inv[inv$inversion %in% c("inv6.0", "inv7.0", "inv7.3", "inv15.1", "inv18.0", "inv22.0"), ]
inv <- inv[, c(5, 3, 4)]

A <- toGRanges(dmr)
B <- toGRanges(inv)
numOverlaps(A,B)

pt <- permTest(A = A, ntimes = 1000, alternative = "auto", genome = genome, randomize.function = randomizeRegions,
               evaluate.function = numOverlaps, B = B, verbose = F)

# Get coordinates for genes under positive selection
# Two-way selection scan
two <- read.csv("Table_S2_twoWayOutlier_genes.csv")
two <- two[c(1:992), ]

# Get genes from two-way selection scan that are found in the Pman Isoquant annotation
ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = F,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))
ann <- ann[ann$Class == "gene", ]
ann$Name <- sub('.*gene_id ([^;]+);.*', '\\1', ann$Name)
ann_twoway <- ann[ann$Name %in% two$gene_symbol_pman, ]
ann_twoway <- ann_twoway[, c(10, 4, 5)]

B <- toGRanges(ann_twoway)
numOverlaps(A,B)

pt <- permTest(A = A, ntimes = 1000, alternative = "auto", genome = genome, randomize.function = randomizeRegions,
               evaluate.function = numOverlaps, B = B, verbose = F)

# Three-way selection scan
three <- read.csv("Table_S4_threeWayOutliers.csv")
ann_threeway <- ann[ann$Name %in% three$gene_pman, ]
ann_threeway <- ann_threeway[, c(10, 4, 5)]

B <- toGRanges(ann_threeway)
numOverlaps(A,B)

pt <- permTest(A = A, ntimes = 1000, alternative = "auto", genome = genome, randomize.function = randomizeRegions,
               evaluate.function = numOverlaps, B = B, verbose = F)

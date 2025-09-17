library(data.table)
library(GenomicRanges)
library(future.apply)
library(stringr)

plan(multicore, workers = parallel::detectCores())

ann <- read.table("IsoquantAnnotation.gtf", sep = "\t", header = F,
                  col.names = c("Contig", "Source", "Class", "Start", "End", "V6", "Strand", "V8", "Name"))
ann <- ann[ann$Class == "gene", ]
ann <- ann[!grepl("novel", ann$Name), ]
ann$Name <- sub('.*gene_id ([^;]+);.*', '\\1', ann$Name)
rownames(ann) <- ann$Name

# W
rna <- read.csv("EP_Pman_ExtMMFrac_readcounts_Exon.csv")
rna <- rna[rna$Geneid %in% ann$Name, ]

# For LZ:
# rna <- read.csv("LP_Pman_ExtMMFrac_readcounts_Exon.csv")
# For JZ:
# rna <- read.csv("LP_JZ_Pman_ExtMMFrac_readcounts_Exon.csv")

setDT(ann, keep.rownames = T)
setnames(ann, old = c(names(ann)[2], names(ann)[5], names(ann)[6]), new = c("Contig", "Start", "End"))

ann_expanded <- ann[, .(Contig, Start = Start - 500, End = End + 500)]
ann_gr <- GRanges(seqnames = ann_expanded$Contig,
                  ranges = IRanges(start = ann_expanded$Start, end = ann_expanded$End))

meth_dir <- "/CpG/"
meth_files <- list.files(path = meth_dir, pattern = "\\.meth$", full.names = TRUE)

file_results <- future_lapply(meth_files, function(file) {
  meth <- fread(file, col.names = c("Contig", "Pos", "End", "Strand", "Score", "Coverage"))
  # Filter on coverage
  meth <- meth[Coverage >= 5]
  meth_gr <- GRanges(seqnames = meth$Contig,
                     ranges = IRanges(start = meth$Pos, width = 1),
                     Score = meth$Score)
  hits <- findOverlaps(ann_gr, meth_gr)
  overlaps_dt <- data.table(queryHits = queryHits(hits),
                            Score = mcols(meth_gr)$Score[subjectHits(hits)])
  avg_by_region <- overlaps_dt[, .(avg = mean(Score, na.rm = TRUE)), by = queryHits]
  avg_vec <- rep(NA_real_, length(ann_gr))
  avg_vec[avg_by_region$queryHits] <- avg_by_region$avg
  colname <- str_remove(basename(file), "\\.meth$")
  setnames(data.table(avg_vec), str_remove(basename(file), "\\.meth$"))
})

methylation_data <- do.call(cbind, file_results)

final <- cbind(ann[, .(Contig, Start, End)], methylation_data)
final <- data.frame(final_output)
rownames(final) <- ann$rn
fwrite(final, "/methylation_summary.tsv", sep = "\t", quote = F, row.names = T)

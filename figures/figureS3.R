library(ggplot2)

chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)

w <- read.table("radmeth.w.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                 col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
w$Chrom <- chrs$Sequence.name[match(w$Contig, chrs$RefSeq.seq.accession)]
w$Density <- (w$End - w$Start)/w$Count
w <- w[w$Count > 29, ]
w <- w[w$Density < 50, ]


lz <- read.table("radmeth.lz.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
lz$Chrom <- chrs$Sequence.name[match(lz$Contig, chrs$RefSeq.seq.accession)]
lz$Density <- (lz$End - lz$Start)/lz$Count
lz <- lz[lz$Count > 29, ]
lz <- lz[lz$Density < 50, ]

jz <- read.table("radmeth.jz.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                 col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jz$Chrom <- chrs$Sequence.name[match(jz$Contig, chrs$RefSeq.seq.accession)]
jz$Density <- (jz$End - jz$Start)/jz$Count
jz <- jz[jz$Count > 29, ]
jz <- jz[jz$Density < 50, ]

vec <- c()
for (i in 1:nrow(w)) {
  y_lz <- which(
    ((w$Chrom[i] == lz$Chrom) & (w$Start[i] >= lz$Start) & (w$Start[i] <= lz$End)) |
      ((w$Chrom[i] == lz$Chrom) & (w$End[i] >= lz$Start) & (w$End[i] <= lz$End))
  )
  
  y_jz <- which(
    ((w$Chrom[i] == jz$Chrom) & (w$Start[i] >= jz$Start) & (w$Start[i] <= jz$End)) |
      ((w$Chrom[i] == jz$Chrom) & (w$End[i] >= jz$Start) & (w$End[i] <= jz$End))
  )
  
  if (length(y_lz) >= 1 & length(y_jz) >= 1) {
    x <- w[i, c(7, 2, 3)]
    vec <- rbind(vec, x)
  }
}

vec$Chrom <- factor(vec$Chrom, levels = unique(vec$Chrom))
vec$Chrom <- gsub("chr", "", vec$Chrom)
vec$Chrom <- as.numeric(vec$Chrom)
axis <- c(1,3,5,7,9,11,13,15,17,19,21,23)

p1 <- ggplot(vec, aes(x = Chrom)) +
  geom_bar() +
  scale_x_continuous(label = axis, breaks = axis) +
  labs(y = "DMRs conserved across all tissues", x = "Chromosome") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"))

png("FigS3.png", res = 300, width = 9, height = 6, units = "in")
p1
dev.off()

png("FigS3.pdf", width = 9, height = 6)
p1
dev.off()


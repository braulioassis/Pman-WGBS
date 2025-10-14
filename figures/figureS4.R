library(ggplot2)
library(gridExtra)

chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)

w <- read.table("radmeth.w.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
w$Chrom <- chrs$Sequence.name[match(w$Contig, chrs$RefSeq.seq.accession)]
w$Density <- (w$End - w$Start)/w$Count
w <- w[w$Count > 29, ]
w <- w[w$Density < 50, ]


w.bw <- read.table("radmeth.w.bw.trt.adj.dmr.bed", sep = "\t", header = F,
                  col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
w.bw$Chrom <- chrs$Sequence.name[match(w.bw$Contig, chrs$RefSeq.seq.accession)]
w.bw <- w.bw[w.bw$Count > 7, ]

w.me <- read.table("radmeth.w.me.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
w.me$Chrom <- chrs$Sequence.name[match(w.me$Contig, chrs$RefSeq.seq.accession)]
w.me <- w.me[w.me$Count > 7, ]

lz <- read.table("radmeth.lz.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
lz$Chrom <- chrs$Sequence.name[match(lz$Contig, chrs$RefSeq.seq.accession)]
lz$Density <- (lz$End - lz$Start)/lz$Count
lz <- lz[lz$Count > 29, ]
lz <- lz[lz$Density < 50, ]

lz.bw <- read.table("radmeth.lz.bw.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
lz.bw$Chrom <- chrs$Sequence.name[match(lz.bw$Contig, chrs$RefSeq.seq.accession)]
lz.bw <- lz.bw[lz.bw$Count > 7, ]

lz.me <- read.table("radmeth.lz.me.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
lz.me$Chrom <- chrs$Sequence.name[match(lz.me$Contig, chrs$RefSeq.seq.accession)]
lz.me <- lz.me[lz.me$Count > 7, ]

jz <- read.table("radmeth.jz.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                 col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jz$Chrom <- chrs$Sequence.name[match(jz$Contig, chrs$RefSeq.seq.accession)]
jz$Density <- (jz$End - jz$Start)/jz$Count
jz <- jz[jz$Count > 29, ]
jz <- jz[jz$Density < 50, ]

jz.bw <- read.table("radmeth.jz.bw.trt.adj.dmr.bed", sep = "\t", header = F,
                    col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jz.bw$Chrom <- chrs$Sequence.name[match(jz.bw$Contig, chrs$RefSeq.seq.accession)]
jz.bw <- jz.bw[jz.bw$Count > 7, ]

jz.me <- read.table("radmeth.jz.me.trt.adj.dmr.bed", sep = "\t", header = F,
                    col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jz.me$Chrom <- chrs$Sequence.name[match(jz.me$Contig, chrs$RefSeq.seq.accession)]
jz.me <- jz.me[jz.me$Count > 7, ]

p1 <- ggplot(w, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "DMRs", title = "Population-level differential methylation\nEarly gestation, whole placenta") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p4 <- ggplot(w.bw, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "DMRs", title = "Exposure-level differential methylation\nLowland ancestry") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p7 <- ggplot(w.me, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "Differential methylation in hypoxia", y = "DMRs", title = "Exposure-level differential methylation\nHighland ancestry") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p2 <- ggplot(lz, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "", title = "Population-level differential methylation\nLate gestation, labyrinth") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p5 <- ggplot(lz.bw, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "", title = "Exposure-level differential methylation\nLowland ancestry") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p8 <- ggplot(lz.me, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "Differential methylation in hypoxia", y = "", title = "Exposure-level differential methylation\nHighland ancestry") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p3 <- ggplot(jz, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "", title = "Population-level differential methylation\nLate gestation, junctional zone") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p6 <- ggplot(jz.bw, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "", title = "Exposure-level differential methylation\nLowland ancestry") +
  xlim(-1,1) + 
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p9 <- ggplot(jz.me, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "Differential methylation in hypoxia", y = "", title = "Exposure-level differential methylation\nHighland ancestry") +
  theme_bw() +
  xlim(-1,1) + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

png("FigS4.png", res = 300, width = 12, height = 10, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
dev.off()

pdf("FigS4.pdf", width = 12, height = 10)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)
dev.off()

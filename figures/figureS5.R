library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)
library(ggsignif)

# P1
wgbs <- read.csv("jz.summary.CpG.csv")
wgbs <- wgbs[wgbs$FileName != "JZ6" & wgbs$FileName != "JZ52", ]
wgbs$Population[wgbs$Population == "ME"] <- "Highland"
wgbs$Population[wgbs$Population == "BW"] <- "Lowland"
wgbs$Population <- factor(wgbs$Population, levels = c("Lowland", "Highland"))
wgbs$Treatment[wgbs$Treatment == "1N"] <- "Normoxia"
wgbs$Treatment[wgbs$Treatment == "2H"] <- "Hypoxia"
wgbs$Treatment <- factor(wgbs$Treatment, levels = c("Normoxia", "Hypoxia"))

p1 <- ggplot(wgbs[wgbs$MeanDepth > 2 & wgbs$Population == "Lowland", ], aes(x = Treatment, y = MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, fill = "goldenrod1") +
  geom_jitter(color = "black", width = 0.15, alpha = 0.6) +
  geom_signif(comparisons = list(c("Normoxia", "Hypoxia")), annotations = "n.s.", textsize = 3) +
  labs(title = "Global methylation",
       y = "CpG methylation rate",
       x = "") +
  ylim(0.4, 0.65) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        legend.position = c(0.75, 0.8), legend.background = element_rect(color = "black"))

# P2
bw <- read.table("radmeth.jz.bw.trt.adj.dmr.bed")
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)
bw$Chrom <- chrs$Sequence.name[match(bw$V1, chrs$RefSeq.seq.accession)]
bw <- bw[!grepl("scaff|MT", bw$Chrom), ]
bw$Chrom <- gsub("chr", "", bw$Chrom)
bw$Chrom[bw$Chrom == "X"] <- 24
bw$Chrom <- as.numeric(bw$Chrom)
bw$V2 <- as.numeric(bw$V2)
bw$V5 <- as.numeric(bw$V5)
unique(bw$Chrom)
bw <- na.omit(bw[, c("Chrom", "V2", "V5")])

bw <- bw %>% 
  
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(chr_len=max(V2)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(bw, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, V2) %>%
  mutate(BPcum=V2+tot)

axisbw = bw %>%
  group_by(Chrom) %>%
  summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

axisbw <- axisbw[c(1,3,5,7,9,11,13,15,17,19,21,23), ]

p2 <- ggplot(bw, aes(x=BPcum, y=V5)) +
  geom_hline(yintercept = 8, color = "deepskyblue2", size = 1, alpha = 0.7) +
  geom_point(aes(color=as.factor(Chrom)), alpha=0.8, size = 1.3) +
  scale_color_manual(values = rep(c("goldenrod1", "grey"), 23)) +
  scale_x_continuous(label = axisbw$Chrom, breaks = axisbw$center) +
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Lowlanders, late gestation, junctional zone",
       y = "Differential methylation sites",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P3
wgbs <- read.csv("jz.summary.CpG.csv")
wgbs <- wgbs[wgbs$FileName != "JZ6" & wgbs$FileName != "JZ52", ]
wgbs$Population[wgbs$Population == "ME"] <- "Highland"
wgbs$Population[wgbs$Population == "BW"] <- "Lowland"
wgbs$Population <- factor(wgbs$Population, levels = c("Lowland", "Highland"))
wgbs$Treatment[wgbs$Treatment == "1N"] <- "Normoxia"
wgbs$Treatment[wgbs$Treatment == "2H"] <- "Hypoxia"
wgbs$Treatment <- factor(wgbs$Treatment, levels = c("Normoxia", "Hypoxia"))

p3 <- ggplot(wgbs[wgbs$MeanDepth > 2 & wgbs$Population == "Highland", ], aes(x = Treatment, y = MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, fill = "deepskyblue3") +
  geom_jitter(color = "black", width = 0.15, alpha = 0.6) +
  geom_signif(comparisons = list(c("Normoxia", "Hypoxia")), annotations = "n.s.", textsize = 3) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.4, 0.65) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        legend.position = c(0.75, 0.8), legend.background = element_rect(color = "black"))

# P4
me <- read.table("radmeth.jz.me.trt.adj.dmr.bed")
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)
me$Chrom <- chrs$Sequence.name[match(me$V1, chrs$RefSeq.seq.accession)]
me <- me[!grepl("scaff|MT", me$Chrom), ]
me$Chrom <- gsub("chr", "", me$Chrom)
me$Chrom[me$Chrom == "X"] <- 24
me$Chrom <- as.numeric(me$Chrom)
me$V2 <- as.numeric(me$V2)
me$V5 <- as.numeric(me$V5)
unique(me$Chrom)
me <- na.omit(me[, c("Chrom", "V2", "V5")])

me <- me %>% 
  
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(chr_len=max(V2)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(me, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, V2) %>%
  mutate(BPcum=V2+tot)

axisme = me %>%
  group_by(Chrom) %>%
  summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

axisme <- axisme[c(1,3,5,7,9,11,13,15,17,19,21,23), ]

p4 <- ggplot(me, aes(x=BPcum, y=V5)) +
  geom_hline(yintercept = 8, color = "goldenrod1", size = 1, alpha = 0.7) +
  geom_point(aes(color=as.factor(Chrom)), alpha=0.8, size = 1.3) +
  scale_color_manual(values = rep(c("deepskyblue3", "grey"), 23)) +
  scale_x_continuous(label = axisme$Chrom, breaks = axisme$center) +
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.05))) + # remove space between plot area and x axis
  labs(title = "Highlanders, late gestation, junctional zone",
       y = "Differential methylation sites",
       x = "Chromosome") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P5
jz.bw <- read.table("radmeth.jz.bw.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jz.bw$Chrom <- chrs$Sequence.name[match(jz.bw$Contig, chrs$RefSeq.seq.accession)]
jz.bw <- jz.bw[jz.bw$Count > 7, ]

p5 <- ggplot(jz.bw, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "Number of significant DMRs", title = "DMRs in hypoxia") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  ylim(0, 10) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P6
jz.me <- read.table("radmeth.jz.me.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jz.me$Chrom <- chrs$Sequence.name[match(jz.me$Contig, chrs$RefSeq.seq.accession)]
jz.me <- jz.me[jz.me$Count > 7, ]

p6 <- ggplot(jz.me, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "Methylation diff. in hypoxia", y = "Number of significant DMRs") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  ylim(0, 10) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# Layout
layout <- rbind(
  c(1, 2, 2, 5),
  c(3, 4, 4, 6)
)

# Render
png("Figure S5.png", res = 300, width = 8.5, height = 6, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)
dev.off()

pdf("Figure S5.pdf", width = 8.5, height = 6)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)
dev.off()

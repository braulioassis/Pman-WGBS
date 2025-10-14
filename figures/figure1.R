library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)
library(ggsignif)

# P1
wgbs <- read.csv("w.summary.CpG.csv")
wgbs$Population[wgbs$Population == "ME"] <- "Highland"
wgbs$Population[wgbs$Population == "BW"] <- "Lowland"
wgbs$Population <- factor(wgbs$Population, levels = c("Lowland", "Highland"))
wgbs$Treatment[wgbs$Treatment == "1N"] <- "Normoxia"
wgbs$Treatment[wgbs$Treatment == "2H"] <- "Hypoxia"
wgbs$Treatment <- factor(wgbs$Treatment, levels = c("Normoxia", "Hypoxia"))

p1 <- ggplot(wgbs[wgbs$MeanDepth > 2 & wgbs$Population == "Lowland", ], aes(x = Treatment, y = MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 1.3, fill = "goldenrod1") +
  geom_jitter(color = "black", width = 0.15, alpha = 0.6) +
  geom_signif(comparisons = list(c("Normoxia", "Hypoxia")), annotations = "n.s.", textsize = 5) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.35, 0.65) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 20),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        legend.position = c(0.75, 0.8), legend.background = element_rect(color = "black"))

# P2
bw <- read.table("radmeth.w.bw.trt.adj.dmr.bed")
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
  scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.05))) + # remove space between plot area and x axis
  labs(title = "Lowlanders, early gestation, whole placenta",
       y = "Differential methylation sites",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P3
wgbs <- read.csv("w.summary.CpG.csv")
wgbs$Population[wgbs$Population == "ME"] <- "Highland"
wgbs$Population[wgbs$Population == "BW"] <- "Lowland"
wgbs$Population <- factor(wgbs$Population, levels = c("Lowland", "Highland"))
wgbs$Treatment[wgbs$Treatment == "1N"] <- "Normoxia"
wgbs$Treatment[wgbs$Treatment == "2H"] <- "Hypoxia"
wgbs$Treatment <- factor(wgbs$Treatment, levels = c("Normoxia", "Hypoxia"))

p3 <- ggplot(wgbs[wgbs$MeanDepth > 2 & wgbs$Population == "Highland", ], aes(x = Treatment, y = MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 1.3, fill = "deepskyblue3") +
  geom_jitter(color = "black", width = 0.15, alpha = 0.6) +
  geom_signif(comparisons = list(c("Normoxia", "Hypoxia")), annotations = "n.s.", textsize = 5) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.35, 0.65) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 20),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        legend.position = c(0.75, 0.8), legend.background = element_rect(color = "black"))

# P4
me <- read.table("radmeth.w.me.trt.adj.dmr.bed")
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
  labs(title = "Highlanders, early gestation, whole placenta",
       y = "Differential methylation sites",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

layout <- rbind(
  c(1, 2, 2, 2),
  c(3, 4, 4, 4)
  )

png("Fig1.png", res = 300, width = 12, height = 8, units = "in")
grid.arrange(p1, p2, p3, p4, layout_matrix = layout)
dev.off()

pdf("Fig1.pdf", width = 12, height = 8)
grid.arrange(p1, p2, p3, p4, layout_matrix = layout)
dev.off()

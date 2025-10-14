library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggsignif)

# P1
wgbs <- read.csv("w.summary.CpG.csv")
wgbs$Population[wgbs$Population == "ME"] <- "Highland"
wgbs$Population[wgbs$Population == "BW"] <- "Lowland"
wgbs$Population <- factor(wgbs$Population, levels = c("Lowland", "Highland"))
wgbs$Treatment[wgbs$Treatment == "1N"] <- "Normoxia"
wgbs$Treatment[wgbs$Treatment == "2H"] <- "Hypoxia"
wgbs$Treatment <- factor(wgbs$Treatment, levels = c("Normoxia", "Hypoxia"))

p1 <- ggplot(wgbs[wgbs$MeanDepth > 2, ], aes(Population, MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 1.3, fill = "#C0AED9FF") +
  geom_signif(comparisons = list(c("Lowland", "Highland")), annotations = "p < 2.2e-16", textsize = 5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.35, 0.7) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"))

# P2
w <- read.table("radmeth.w.pop.trt.filteredCTGA.adj.dmr.bed")
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)
w$Chrom <- chrs$Sequence.name[match(w$V1, chrs$RefSeq.seq.accession)]
w <- w[!grepl("scaff|MT", w$Chrom), ]
w$Chrom <- gsub("chr", "", w$Chrom)
w$Chrom[w$Chrom == "X"] <- 24
w$Chrom <- as.numeric(w$Chrom)
w$V2 <- as.numeric(w$V2)
w$V5 <- as.numeric(w$V5)
w$Gap <- (w$V3 - w$V2)/w$V5
w <- w[w$Gap < 50, ]
unique(w$Chrom)
w <- na.omit(w[, c("Chrom", "V2", "V5")])

w <- w %>% 
  
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(chr_len=max(V2)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(w, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, V2) %>%
  mutate(BPcum=V2+tot)

axisw = w %>%
  group_by(Chrom) %>%
  summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

axisw <- axisw[c(1,3,5,7,9,11,13,15,17,19,21,23), ]

p2 <- ggplot(w, aes(x=BPcum, y=V5)) +
  geom_point( aes(color=as.factor(Chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#C0AED9FF", "grey"), 23)) +
  scale_x_continuous(label = axisw$Chrom, breaks = axisw$center) +
  scale_y_continuous(limits = c(0, 201), expand = expansion(mult = c(0, 0.05))) + # remove space between plot area and x axis
  geom_hline(yintercept = 30, color = "forestgreen", size = 1) +
  labs(title = "Early gestation, whole placenta",
       y = "Differential methylation sites",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P3
LZ <- read.csv("lz.summary.CpG.csv")
LZ$Population[LZ$Population == "ME"] <- "Highland"
LZ$Population[LZ$Population == "BW"] <- "Lowland"
LZ$Population <- factor(LZ$Population, levels = c("Lowland", "Highland"))

p3 <- ggplot(LZ[LZ$MeanDepth > 2, ], aes(Population, MeanMeth)) +
  geom_boxplot(linewidth = 1.3, fill = "#C0AED9FF") +
  geom_signif(comparisons = list(c("Lowland", "Highland")), annotations = "p < 2.2e-16", textsize = 5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.35, 0.7) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"))

# P4
lzman <- read.table("radmeth.lz.pop.trt.filteredCTGA.adj.dmr.bed")
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$RefSeq.seq.accession)
lzman$Chrom <- chrs$Sequence.name[match(lzman$V1, chrs$RefSeq.seq.accession)]
lzman <- lzman[!grepl("scaff|MT", lzman$Chrom), ]
lzman$Chrom <- gsub("chr", "", lzman$Chrom)
lzman$Chrom[lzman$Chrom == "X"] <- 24
lzman$Chrom <- as.numeric(lzman$Chrom)
lzman$V2 <- as.numeric(lzman$V2)
lzman$V5 <- as.numeric(lzman$V5)
lzman$Gap <- (lzman$V3 - lzman$V2)/lzman$V5
lzman <- lzman[lzman$Gap < 50, ]

unique(lzman$Chrom)
lzman <- na.omit(lzman[, c("Chrom", "V2", "V5")])

lzman <- lzman %>% 
  
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(chr_len=max(V2)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(lzman, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, V2) %>%
  mutate(BPcum=V2+tot)

axislzman = lzman %>%
  group_by(Chrom) %>%
  summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

axislzman <- axislzman[c(1,3,5,7,9,11,13,15,17,19,21,23), ]

p4 <- ggplot(lzman, aes(x=BPcum, y=V5)) +
  geom_point( aes(color=as.factor(Chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#C0AED9FF", "grey"), 23)) +
  scale_x_continuous(label = axislzman$Chrom, breaks = axislzman$center) +
  scale_y_continuous(limits = c(0, 201), expand = expansion(mult = c(0, 0.05))) + # remove space between plot area and x axis
  geom_hline(yintercept = 30, color = "forestgreen", size = 1) +
  labs(title = "Late gestation, labyrinth",
       y = "Differential methylation sites",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P5
JZ <- read.table("jz.summary.CpG.csv")
JZ <- JZ[JZ$FileName != "JZ134", ]
JZ$Population[JZ$Population == "ME"] <- "Highland"
JZ$Population[JZ$Population == "BW"] <- "Lowland"
JZ$Population <- factor(JZ$Population, levels = c("Lowland", "Highland"))
JZ$Treatment[JZ$Treatment == "1N"] <- "Normoxia"
JZ$Treatment[JZ$Treatment == "2H"] <- "Hypoxia"
JZ$Treatment <- factor(JZ$Treatment, levels = c("Normoxia", "Hypoxia"))

p5 <- ggplot(wgbs[wgbs$MeanDepth > 2, ], aes(Population, MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 1.3, fill = "#C0AED9FF") +
  geom_signif(comparisons = list(c("Lowland", "Highland")), annotations = "p = 1.4e-13", textsize = 5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(y = "CpG methylation rate",
       x = "Population") +
  ylim(0.35, 0.7) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"))

# P6
jzman <- read.table("radmeth.jz.pop.trt.filteredCTGA.adj.dmr.bed")
chrs <- read.table("sequence_report.tsv", header = T, sep = "\t")
chrs$GenBank.seq.accession <- gsub("\\.2$", ".1", chrs$GenBank.seq.accession)
jzman$Chrom <- chrs$Sequence.name[match(jzman$V1, chrs$RefSeq.seq.accession)]
jzman <- jzman[!grepl("scaff|MT", jzman$Chrom), ]
jzman$Chrom <- gsub("chr", "", jzman$Chrom)
jzman$Chrom[jzman$Chrom == "X"] <- 24
jzman$Chrom <- as.numeric(jzman$Chrom)
jzman$V2 <- as.numeric(jzman$V2)
jzman$V5 <- as.numeric(jzman$V5)
jzman$Gap <- (jzman$V3 - jzman$V2)/jzman$V5
jzman <- jzman[jzman$Gap < 50, ]
unique(jzman$Chrom)
jzman <- na.omit(jzman[, c("Chrom", "V2", "V5")])

jzman <- jzman %>% 
  
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(chr_len=max(V2)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(jzman, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, V2) %>%
  mutate(BPcum=V2+tot)

axisjz = jzman %>%
  group_by(Chrom) %>%
  summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

axisjz <- axisjz[c(1,3,5,7,9,11,13,15,17,19,21,23), ]

p6 <- ggplot(jzman, aes(x=BPcum, y=V5)) +
  geom_point( aes(color=as.factor(Chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#C0AED9FF", "grey"), 23)) +
  scale_x_continuous(label = axisjz$Chrom, breaks = axisjz$center) +
  scale_y_continuous(limits = c(0, 201), expand = expansion(mult = c(0, 0.05))) + # remove space between plot area and x axis
  geom_hline(yintercept = 30, color = "forestgreen", size = 1) +
  labs(title = "Late gestation, junctional zone",
       y = "Differential methylation sites",
       x = "Chromosome") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

layout <- rbind(
  c(1, 2, 2, 2), # Row 1: p1 (col 1), p2 (cols 2–4)
  c(3, 4, 4, 4), # Row 2: p3 (col 1), p4 (cols 2–4)
  c(5, 6, 6, 6)  # Row 3: p5 (col 1), p6 (cols 2–4)
)

png("Fig1.png", res = 300, width = 14, height = 12, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)
dev.off()

pdf("Fig1.pdf", width = 14, height = 12)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)
dev.off()

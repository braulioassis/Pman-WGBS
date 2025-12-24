library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggsignif)
library(viridis)

# P1
wgbs <- read.csv("w.summary.CpG.csv")
wgbs$Population[wgbs$Population == "ME"] <- "Highland"
wgbs$Population[wgbs$Population == "BW"] <- "Lowland"
wgbs$Population <- factor(wgbs$Population, levels = c("Lowland", "Highland"))
wgbs$Treatment[wgbs$Treatment == "1N"] <- "Normoxia"
wgbs$Treatment[wgbs$Treatment == "2H"] <- "Hypoxia"
wgbs$Treatment <- factor(wgbs$Treatment, levels = c("Normoxia", "Hypoxia"))

p1 <- ggplot(wgbs[wgbs$MeanDepth > 2, ], aes(Population, MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, fill = "#C0AED9FF") +
  geom_signif(comparisons = list(c("Lowland", "Highland")), annotations = "p < 2.2e-16", textsize = 3) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(y = "CpG methylation rate",
       x = "",
       title = "Global methylation") +
  ylim(0.35, 0.7) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
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
  geom_boxplot(linewidth = 0.8, fill = "#C0AED9FF") +
  geom_signif(comparisons = list(c("Lowland", "Highland")), annotations = "p < 2.2e-16", textsize = 3) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.35, 0.7) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P5
JZ <- read.csv("jz.summary.CpG.csv")
JZ <- JZ[JZ$FileName != "JZ6" & JZ$FileName != "JZ52" & JZ$FileName != "JZ134", ]
JZ$Population[JZ$Population == "ME"] <- "Highland"
JZ$Population[JZ$Population == "BW"] <- "Lowland"
JZ$Population <- factor(JZ$Population, levels = c("Lowland", "Highland"))
JZ$Treatment[JZ$Treatment == "1N"] <- "Normoxia"
JZ$Treatment[JZ$Treatment == "2H"] <- "Hypoxia"
JZ$Treatment <- factor(JZ$Treatment, levels = c("Normoxia", "Hypoxia"))

p5 <- ggplot(JZ[JZ$MeanDepth > 2, ], aes(Population, MeanMeth)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, fill = "#C0AED9FF") +
  geom_signif(comparisons = list(c("Lowland", "Highland")), annotations = "p = 2.2e-15", textsize = 3) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(y = "CpG methylation rate",
       x = "Population") +
  ylim(0.35, 0.7) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P7
wdiff <- read.table("radmeth.w.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
wdiff$Chrom <- chrs$Sequence.name[match(wdiff$Contig, chrs$RefSeq.seq.accession)]
wdiff$Density <- (wdiff$End - wdiff$Start)/wdiff$Count
wdiff <- wdiff[wdiff$Count > 29, ]
wdiff <- wdiff[wdiff$Density < 50, ]

p7 <- ggplot(wdiff, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "Number of significant DMRs", title = "DMRs in highlanders") +
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  ylim(0, 250) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P8
lzdiff <- read.table("radmeth.lz.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                 col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
lzdiff$Chrom <- chrs$Sequence.name[match(lzdiff$Contig, chrs$RefSeq.seq.accession)]
lzdiff$Density <- (lzdiff$End - lzdiff$Start)/lzdiff$Count
lzdiff <- lzdiff[lzdiff$Count > 29, ]
lzdiff <- lzdiff[lzdiff$Density < 50, ]

p8 <- ggplot(lzdiff, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "Number of significant DMRs") +
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  ylim(0, 250) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P9
jzdiff <- read.table("radmeth.jz.pop.trt.filteredCTGA.adj.dmr.bed", sep = "\t", header = F, 
                 col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
jzdiff$Chrom <- chrs$Sequence.name[match(jzdiff$Contig, chrs$RefSeq.seq.accession)]
jzdiff$Density <- (jzdiff$End - jzdiff$Start)/jzdiff$Count
jzdiff <- jzdiff[jzdiff$Count > 29, ]
jzdiff <- jzdiff[jzdiff$Density < 50, ]

p9 <- ggplot(jzdiff, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "Methylation diff. in highlanders", y = "Number of significant DMRs") +
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  ylim(0, 250) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P10
df <- read.csv("GO table.csv")
df$Tissue <- factor(df$Tissue, levels = c("W", "LZ", "JZ"))
df$term_name[df$term_name == "Growth hormone synthesis, secretion and action"] <- "Growth hormone synthesis[...]"
df$term_name[df$term_name == "Cortisol synthesis and secretion"] <- "Cortisol synthesis[...]"
df$logF <- log(df$Fold.enrichment)
#df$term_name <- paste0(df$source, ": ", df$term_name)
df$term_name <- factor(df$term_name, levels = rev(c("Growth hormone synthesis[...]",
                                                    "Longevity regulating pathway",
                                                    "GnRH secretion",
                                                    "Cortisol synthesis[...]",
                                                    "Insulin secretion",
                                                    "Estrogen signaling",
                                                    "FoxO signaling",
                                                    "ZBTB14",
                                                    "KLF6",
                                                    "CTCF",
                                                    "ZNF354C",
                                                    "E2F1")))

kegg <- df[df$source == "KEGG",]
tf <- df[df$source == "TF",]

p10 <- ggplot(kegg, aes(x = Tissue, y = term_name)) +
  geom_point(aes(fill = negative_log10_of_adjusted_p_value, size = logF), color = "black", shape = 21, stroke = 1) +
  scale_fill_viridis(option = "plasma", limits=c(0, 2.5), breaks = c(0, 1.25, 2.5)) +
  scale_size_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  labs(title = "KEGG enrichments",
       y = "",
       x = "",
       fill = "-log10(FDR)",
       size = "log Fold-enrichment") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.key.size = unit(0.15, "cm"),
        legend.title = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(
    fill = guide_colorbar(order = 2),
    size = guide_legend(order = 1)) +
  coord_cartesian()


# P11
p11 <- ggplot(tf, aes(x = Tissue, y = term_name)) +
  geom_point(aes(fill = negative_log10_of_adjusted_p_value, size = logF), color = "black", shape = 21, stroke = 1) +
  scale_fill_viridis(option = "plasma", limits=c(0, 30), breaks = c(0, 15, 30)) +
  labs(title = "Transcription factor motif enrichments",
       y = "",
       x = "",
       fill = "-log10(FDR)",
       size = "log Fold-enrichment") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.key.size = unit(0.15, "cm"),
        legend.title = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(
    fill = guide_colorbar(order = 2),
    size = guide_legend(order = 1)) +
  coord_cartesian()

# Layout
layout <- rbind(
  c(1, 2, 2, 7),
  c(3, 4, 4, 8),
  c(5, 6, 6, 9),
  c(10, 10, 11, 11)
)

png("Figure 1.png", res = 300, width = 8.5, height = 11, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, layout_matrix = layout)
dev.off()

pdf("Figure 1.pdf", width = 8.5, height = 11)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, layout_matrix = layout)
dev.off()

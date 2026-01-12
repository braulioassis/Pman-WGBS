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
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, fill = "goldenrod1") +
  geom_jitter(color = "black", width = 0.15, alpha = 0.6) +
  geom_signif(comparisons = list(c("Normoxia", "Hypoxia")), annotations = "n.s.", textsize = 3) +
  labs(y = "CpG methylation rate",
       x = "",
       title = "Global methylation") +
  ylim(0.35, 0.65) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
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
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, fill = "deepskyblue3") +
  geom_jitter(color = "black", width = 0.15, alpha = 0.6) +
  geom_signif(comparisons = list(c("Normoxia", "Hypoxia")), annotations = "n.s.", textsize = 3) +
  labs(y = "CpG methylation rate",
       x = "") +
  ylim(0.35, 0.65) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
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

# P5
w.bw <- read.table("radmeth.w.bw.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
w.bw$Chrom <- chrs$Sequence.name[match(w.bw$Contig, chrs$RefSeq.seq.accession)]
w.bw <- w.bw[w.bw$Count > 7, ]

p5 <- ggplot(w.bw, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "", y = "Number of significant DMRs", title = "DMRs in hypoxia") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  ylim(0, 45) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P6
w.me <- read.table("radmeth.w.me.trt.adj.dmr.bed", sep = "\t", header = F,
                   col.names = c("Contig", "Start", "End", "DMR", "Count", "Stat"))
w.me$Chrom <- chrs$Sequence.name[match(w.me$Contig, chrs$RefSeq.seq.accession)]
w.me <- w.me[w.me$Count > 7, ]

p6 <- ggplot(w.me, aes(x = Stat)) +
  geom_histogram() +
  labs(x = "Methylation diff. in hypoxia", y = "Number of significant DMRs") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  ylim(0, 45) +
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
df <- read.csv("w.summary.hypoxia.regions.csv")

colnames(df) <- gsub("_CpG", "", colnames(df))

sample_info <- read.csv("w.summary.CpG.csv")
stopifnot(all(colnames(meth_matrix) %in% sample_info$FileName))

sample_info <- sample_info[
  match(colnames(meth_matrix), sample_info$FileName),
]
sample_info <- sample_info[sample_info$Population == "ME" | sample_info$Population == "BW", ]
sample_info$Treatment[sample_info$Treatment == "1N"] <- "Normoxia"
sample_info$Treatment[sample_info$Treatment == "2H"] <- "Hypoxia"
sample_info$Population[sample_info$Population == "ME"] <- "Highland"
sample_info$Population[sample_info$Population == "BW"] <- "Lowland"

sample_info$Group <- with(sample_info,
                          interaction(Population, Treatment, sep = "_"))

sample_info$Group <- factor(sample_info$Group,
                            levels = c("Lowland_Normoxia", "Lowland_Hypoxia",
                                       "Highland_Normoxia", "Highland_Hypoxia")
)

levels(sample_info$Group) <- c(
  "Lowland, normoxia",
  "Lowland, hypoxia",
  "Highland, normoxia",
  "Highland, hypoxia"
)

meth_matrix <- as.matrix(df[ , -(1:3)])
rownames(meth_matrix) <- paste(df$Contig, df$Start, df$End, sep = "_")

new_order <- order(sample_info$Group)

meth_matrix2 <- meth_matrix[, new_order]
sample_info2 <- sample_info[new_order, ]


meth_matrix2 <- t(apply(meth_matrix2, 1, function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

meth_matrix2 <- meth_matrix2[rowSums(is.na(meth_matrix2)) == 0, ]

row_var <- apply(meth_matrix2, 1, var)
meth_matrix2 <- meth_matrix2[row_var > 0, ]

meth_z <- t(scale(t(meth_matrix2), center = TRUE, scale = TRUE))
hm <- pheatmap(
  meth_z,
  silent = T,
  )
row_order <- hm$tree_row$order
col_order <- hm$tree_col$order

df <- as.data.frame(meth_z)
df$row <- rownames(df)

long <- pivot_longer(
  df,
  -row,
  names_to = "column",
  values_to = "value"
)
long$Group <- sample_info2$Group[match(long$column, colnames(meth_z))]
long$Group <- factor(
  long$Group,
  levels = c("Lowland, normoxia", "Lowland, hypoxia",
             "Highland, normoxia", "Highland, hypoxia")
)
long$row <- factor(long$row, levels = rownames(meth_z)[row_order])
long$column <- factor(long$column, levels = colnames(meth_z)[col_order])

p7 <- ggplot(long, aes(column, row, fill = value)) +
  geom_tile() +
  facet_grid(
    . ~ Group,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_viridis_c(
    name   = "Z-scored\nmethylation rate",
    limits = c(-5, 5),
    breaks = c(-5, 0, 5)
  ) +
  labs(title = "Methylation rate, hypoxia-sensitive DMRs") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    panel.spacing.x = unit(0, "pt"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

# Layout
layout <- rbind(
  c(1, 2, 2, 5),
  c(3, 4, 4, 6),
  c(7),
  c(7)
)

# Render
png("Figure 2.png", res = 300, width = 8.5, height = 11, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, layout_matrix = layout)
dev.off()

pdf("Figure 2.pdf", width = 8.5, height = 11)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, layout_matrix = layout)
dev.off()

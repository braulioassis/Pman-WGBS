library(ggplot2)
library(gridExtra)
library(dplyr)
library(Cairo)

# P1
cor_results <- read.csv("W.genebody_cor_results.csv")
cor_results <- na.omit(cor_results)
cor_results$bin <- cut(cor_results$Cor, breaks = seq(floor(min(cor_results$Cor)), ceiling(max(cor_results$Cor)), by = 0.1), include.lowest = T)
cor_counts <- as.data.frame(table(cor_results$bin))
colnames(cor_counts) <- c("bin", "count")
cor_counts$midpoint <- sapply(strsplit(as.character(cor_counts$bin), ","), function(x) {
  as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.05
})
cor_counts$highlight <- with(cor_counts, ifelse(midpoint >= 0.3 | midpoint <= -0.3, "Highlight", "Normal"))

sig_results <- filter(cor_results, abs(Cor) >= 0.3, FDR < 0.05)

p1 <- ggplot(cor_counts, aes(x = midpoint, y = count, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Highlight" = "#C0AED9FF", "Normal" = "grey70")) +
  labs(x = "Spearman \u03C1", y = "Genes", title = "Early gestation, whole placenta") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P2
cor_results <- read.csv("LZ.genebody_cor_results.csv")
cor_results <- na.omit(cor_results)
cor_results$bin <- cut(cor_results$Cor, breaks = seq(floor(min(cor_results$Cor)), ceiling(max(cor_results$Cor)), by = 0.1), include.lowest = T)
cor_counts <- as.data.frame(table(cor_results$bin))
colnames(cor_counts) <- c("bin", "count")
cor_counts$midpoint <- sapply(strsplit(as.character(cor_counts$bin), ","), function(x) {
  as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.05
})
cor_counts$highlight <- with(cor_counts, ifelse(midpoint >= 0.3 | midpoint <= -0.3, "Highlight", "Normal"))

p2 <- ggplot(cor_counts, aes(x = midpoint, y = count, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Highlight" = "#C0AED9FF", "Normal" = "grey70")) +
  labs(x = "Spearman \u03C1", y = "", title = "Late gestation, labyrinth") +
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
cor_results <- read.csv("JZ.genebody_cor_results.csv")
cor_results <- na.omit(cor_results)
cor_results$bin <- cut(cor_results$Cor, breaks = seq(floor(min(cor_results$Cor)), ceiling(max(cor_results$Cor)), by = 0.1), include.lowest = T)
cor_counts <- as.data.frame(table(cor_results$bin))
colnames(cor_counts) <- c("bin", "count")
cor_counts$midpoint <- sapply(strsplit(as.character(cor_counts$bin), ","), function(x) {
  as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.05
})
cor_counts$highlight <- with(cor_counts, ifelse(midpoint >= 0.3 | midpoint <= -0.3, "Highlight", "Normal"))

p3 <- ggplot(cor_counts, aes(x = midpoint, y = count, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Highlight" = "#C0AED9FF", "Normal" = "grey70")) +
  labs(x = "Spearman \u03C1", y = "", title = "Late gestation, junctional zone") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# P4
df <- read.csv("RNA-meth-summary-with-direction.csv")
df$Domain <- factor(df$Domain, levels = c("Whole gene", "Promoter", "Intron", "Exon"))
df$Tissue <- factor(df$Tissue, levels = c("W", "LZ", "JZ"))
df$Correlation <- factor(df$Correlation, levels = c("Positive", "Negative"))

p4 <- ggplot(df[df$Tissue %in% c("W"),], aes(Domain, PercentSignif, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("seagreen3", "brown1")) +
  ylim(0, 55) +
  labs(x = "", y = "Significant correlations (%)") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))


p5 <- ggplot(df[df$Tissue %in% c("LZ"),], aes(Domain, PercentSignif, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("seagreen3", "brown1")) + 
  ylim(0, 55) +
  labs(x = "Gene regions", y = "") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

p6 <- ggplot(df[df$Tissue %in% c("JZ"),], aes(Domain, PercentSignif, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("seagreen3", "brown1")) + 
  ylim(0, 55) +
  labs(x = "", y = "") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.8, 0.8),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

# Layout
layout <- rbind(
  c(1, 2, 3),
  c(4, 5, 6)
  )

# Render
png("Figure 3.png", res = 300, width = 8.5, height = 5, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)
dev.off()

#pdf("Figure 3.pdf", width = 15, height = 8)
#grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)
#dev.off()

ggsave("Figure 3.pdf", width = 8.5, height = 5,
       plot = (grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = layout)), device = cairo_pdf)

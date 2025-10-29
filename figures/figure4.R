library(ggplot2)
library(gridExtra)
library(dplyr)

cor_results <- read.csv("W.genebody_cor_results.csv")
cor_results <- na.omit(cor_results)
cor_results$bin <- cut(cor_results$Cor, breaks = seq(floor(min(cor_results$Cor)), ceiling(max(cor_results$Cor)), by = 0.1), include.lowest = T)
cor_counts <- as.data.frame(table(cor_results$bin))
colnames(cor_counts) <- c("bin", "count")
cor_counts$midpoint <- sapply(strsplit(as.character(cor_counts$bin), ","), function(x) {
  as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.05
})
cor_counts$highlight <- with(cor_counts, ifelse(midpoint > 0.3 | midpoint < -0.3, "Highlight", "Normal"))

sig_results <- filter(cor_results, abs(Cor) >= 0.3, FDR < 0.05)

p1 <- ggplot(cor_counts, aes(x = midpoint, y = count, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Highlight" = "#C0AED9FF", "Normal" = "grey70")) +
  labs(x = "Spearman \u03C1", y = "Genes", title = "Early gestation, whole placenta") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

cor_results <- read.csv("LZ.genebody_cor_results.csv")
cor_results <- na.omit(cor_results)
cor_results$bin <- cut(cor_results$Cor, breaks = seq(floor(min(cor_results$Cor)), ceiling(max(cor_results$Cor)), by = 0.1), include.lowest = T)
cor_counts <- as.data.frame(table(cor_results$bin))
colnames(cor_counts) <- c("bin", "count")
cor_counts$midpoint <- sapply(strsplit(as.character(cor_counts$bin), ","), function(x) {
  as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.05
})
cor_counts$highlight <- with(cor_counts, ifelse(midpoint > 0.3 | midpoint < -0.3, "Highlight", "Normal"))

p2 <- ggplot(cor_counts, aes(x = midpoint, y = count, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Highlight" = "#C0AED9FF", "Normal" = "grey70")) +
  labs(x = "Spearman \u03C1", y = "", title = "Late gestation, labyrinth") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

cor_results <- read.csv("JZ.genebody_cor_results.csv")
cor_results <- na.omit(cor_results)
cor_results$bin <- cut(cor_results$Cor, breaks = seq(floor(min(cor_results$Cor)), ceiling(max(cor_results$Cor)), by = 0.1), include.lowest = T)
cor_counts <- as.data.frame(table(cor_results$bin))
colnames(cor_counts) <- c("bin", "count")
cor_counts$midpoint <- sapply(strsplit(as.character(cor_counts$bin), ","), function(x) {
  as.numeric(gsub("[^0-9.-]", "", x[1])) + 0.05
})
cor_counts$highlight <- with(cor_counts, ifelse(midpoint > 0.3 | midpoint < -0.3, "Highlight", "Normal"))

p3 <- ggplot(cor_counts, aes(x = midpoint, y = count, fill = highlight)) +
  geom_col() +
  scale_fill_manual(values = c("Highlight" = "#C0AED9FF", "Normal" = "grey70")) +
  labs(x = "Spearman \u03C1", y = "", title = "Late gestation, junctional zone") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

png("Fig4.png", res = 300, width = 18, height = 5, units = "in")
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

ggsave("Fig4.pdf", width = 18, height = 5, plot = (grid.arrange(p1, p2, p3, nrow = 1)), device = cairo_pdf)


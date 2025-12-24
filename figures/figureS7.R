library(ggplot2)
library(gridExtra)

df <- read.csv("RNA-meth-summary-with-direction.csv")
df$Domain <- factor(df$Domain, levels = c("Whole gene", "Promoter", "Intron", "Exon"))
df$Tissue <- factor(df$Tissue, levels = c("W", "LZ", "JZ", "W-BW", "LZ-BW", "JZ-BW", "W-ME", "LZ-ME", "JZ-ME"))
df$Correlation <- factor(df$Correlation, levels = c("Positive", "Negative"))

p1 <- ggplot(df[df$Tissue %in% c("W-BW", "LZ-BW", "JZ-BW"),], aes(Domain, PercentSignif, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("seagreen3", "brown1")) + 
  facet_grid(. ~ Tissue) +
  labs(x = "", y = "Significant correlations (%)", title = "Methylation-expression correlations\nLowlanders only") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

p2 <- ggplot(df[df$Tissue %in% c("W-ME", "LZ-ME", "JZ-ME"),], aes(Domain, PercentSignif, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("seagreen3", "brown1")) + 
  facet_grid(. ~ Tissue) +
  ylim(0, 30) +
  labs(x = "", y = "Significant correlations (%)", title = "Methylation-expression correlations\nHighlanders only") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.8, 0.7),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

# Render
png("Figure S7.png",
    res = 300, height = 4, width = 8.5, units = "in")
grid.arrange(p1, p2, ncol = 2)
dev.off()

pdf("Figure S7.pdf",
    height = 4, width = 8.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

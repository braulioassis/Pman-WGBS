library(ggplot2)
library(viridis)
library(cowplot)

df <- read.csv("GO table.csv")
df$Tissue <- factor(df$Tissue, levels = c("W", "LZ", "JZ"))
df$logF <- log(df$Fold.enrichment)
df$term_name <- paste0(df$source, ": ", df$term_name)
df$term_name <- factor(df$term_name, levels = rev(c("KEGG: Growth hormone synthesis, secretion and action",
                                                    "KEGG: Longevity regulating pathway",
                                                    "KEGG: GnRH secretion",
                                                    "KEGG: Cortisol synthesis and secretion",
                                                    "KEGG: Insulin secretion",
                                                    "KEGG: Estrogen signaling",
                                                    "KEGG: FoxO signaling",
                                                    "TF: ZBTB14",
                                                    "TF: KLF6",
                                                    "TF: CTCF",
                                                    "TF: ZNF354C",
                                                    "TF: E2F1")))

kegg <- df[df$source == "KEGG",]
tf <- df[df$source == "TF",]

p1 <- ggplot(kegg, aes(x = Tissue, y = term_name)) +
  geom_point(aes(fill = negative_log10_of_adjusted_p_value, size = logF), color = "black", shape = 21, stroke = 1.5) +
  scale_fill_viridis(option = "plasma", limits=c(0, 2.5), breaks = c(0, 1.25, 2.5)) +
  scale_size_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  labs(title = "KEGG enrichments",
       y = "",
       x = "",
       fill = "-log10(FDR)",
       size = "log Fold-enrichment") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 1),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(
    fill = guide_colorbar(order = 2),
    size = guide_legend(order = 1)) +
  coord_cartesian()

p2 <- ggplot(tf, aes(x = Tissue, y = term_name)) +
  geom_point(aes(fill = negative_log10_of_adjusted_p_value, size = logF), color = "black", shape = 21, stroke = 1.5) +
  scale_fill_viridis(option = "plasma", limits=c(0, 30), breaks = c(0, 15, 30)) +
  labs(title = "Transcription factor motif enrichments",
       y = "",
       x = "",
       fill = "-log10(FDR)",
       size = "log Fold-enrichment") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 1),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(
    fill = guide_colorbar(order = 2),
    size = guide_legend(order = 1)) +
  coord_cartesian()

png("Fig3.png", res = 300, width = 16, height = 8, units = "in")
plot_grid(p1, p2, nrow = 1, align = "hv", axis = "tb")
dev.off()

pdf("Fig3.pdf", width = 16, height = 8)
plot_grid(p1, p2, nrow = 1, align = "hv", axis = "tb")
dev.off()

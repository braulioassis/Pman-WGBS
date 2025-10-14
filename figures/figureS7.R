library(ggplot2)
library(gridExtra)

df <- read.csv("RNA-meth-summary.csv")
df$Domain <- factor(df$Domain, levels = c("Genebody", "Promoter", "Intron", "Exon"))

all <- df[df$Tissue %in% c("W", "LZ", "JZ"), ]
all$Tissue <- factor(all$Tissue, levels = c("W", "LZ", "JZ"))

p1 <- ggplot(all, aes(Domain, Percent, fill = Tissue)) +
  geom_col(position = "dodge") +
  ylim(0, 50) +
  scale_fill_manual(values = c("#003366", "steelblue", "skyblue")) +
  labs(y = "Percent genes significantly correlated", x = "", title = "Across lowlanders and highlanders") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

bw <- df[df$Tissue %in% c("W-BW", "LZ-BW", "JZ-BW"), ]
bw$Tissue <- factor(bw$Tissue, levels = c("W-BW", "LZ-BW", "JZ-BW"))

p2 <- ggplot(bw, aes(Domain, Percent, fill = Tissue)) +
  geom_col(position = "dodge") +
  ylim(0, 50) +
  scale_fill_manual(values = c("#003366", "steelblue", "skyblue")) +
  labs(y = "", x = "", title = "Lowlanders only") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

me <- df[df$Tissue %in% c("W-ME", "LZ-ME", "JZ-ME"), ]
me$Tissue <- factor(me$Tissue, levels = c("W-ME", "LZ-ME", "JZ-ME"))

p3 <- ggplot(me, aes(Domain, Percent, fill = Tissue)) +
  geom_col(position = "dodge") +
  ylim(0, 50) +
  scale_fill_manual(values = c("#003366", "steelblue", "skyblue")) +
  labs(y = "", x = "", title = "Highlanders only") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

png("FigS7.png", res = 300, width = 18, height = 5, units = "in")
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

pdf("FigS7.pdf", width = 18, height = 5)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

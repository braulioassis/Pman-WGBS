library(ggplot2)
library(gridExtra)

df <- read.csv("w.eigen.csv")
pc <- read.csv("w.PCs.csv")

p1 <- ggplot(df, aes(PC1, PC2, fill = Population)) +
  geom_point(size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = c("goldenrod1", "deepskyblue3")) +
  labs(x = "PC1 (82.4% variance explained)", 
       y = "PC2 (0.4% variance explained)", 
       title = "Early gestation") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

df <- read.csv("lz.eigen.csv")
pc <- read.csv("lz.PCs.csv")

p2 <- ggplot(df, aes(PC1, PC2, fill = Population)) +
  geom_point(size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = c("goldenrod1", "deepskyblue3")) +
  labs(x = "PC1 (82.9% variance explained)", y = "PC2 (0.6% variance explained)", title = "Labyrinth zone") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

df <- read.csv("jz.eigen.csv")
rownames(df) <- df$FileName
df <- df[df$FileName != "JZ6" & df$FileName != "JZ52", ]
pc <- read.csv("jz.PCs.csv")

p3 <- ggplot(df, aes(PC1, PC2, fill = Population)) +
  geom_point(size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = c("goldenrod1", "deepskyblue3")) +
  labs(x = "PC1 (87.1% variance explained)", y = "PC2 (0.5% variance explained)", title = "Junctional zone") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "gray85"), panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(color = "#FFFFFF", fill = "#FFFFFF"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"))

png("Figure S1.png",
    res = 300, height = 2.8, width = 8.5, units = "in")
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

pdf("Figure S1.pdf",
    height = 2.8, width = 8.5)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

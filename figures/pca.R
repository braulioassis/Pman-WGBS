library(ggplot2)
library(gridExtra)

df <- read.csv("w.eigen.csv")
pc <- read.csv("w.PCs.csv")

p1 <- ggplot(df, aes(PC1, PC2, fill = Population)) +
  geom_point(size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = c("goldenrod1", "deepskyblue3")) +
  labs(x = "PC1 (82.4% variance explained)", y = "PC2 (0.4% variance explained)", title = "Early gestation") +
  theme_bw()

df <- read.csv("lz.eigen.csv")
pc <- read.csv("lz.PCs.csv")

p2 <- ggplot(df, aes(PC1, PC2, fill = Population)) +
  geom_point(size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = c("goldenrod1", "deepskyblue3")) +
  labs(x = "PC1 (82.9% variance explained)", y = "PC2 (0.6% variance explained)", title = "Labyrinth zone") +
  theme_bw()

df <- read.csv("jz.eigen.csv")
rownames(df) <- df$FileName
df <- df[df$FileName != "JZ6" & df$FileName != "JZ52", ]
pc <- read.csv("jz.PCs.csv")

p3 <- ggplot(df, aes(PC1, PC2, fill = Population)) +
  geom_point(size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = c("goldenrod1", "deepskyblue3")) +
  labs(x = "PC1 (87.1% variance explained)", y = "PC2 (0.5% variance explained)", title = "Junctional zone") +
  theme_bw()

png("PCA.png", res = 300, height = 4, width = 15, units = "in")
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

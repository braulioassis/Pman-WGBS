library(ggplot2)
library(gridExtra)

perms <- readRDS("W-genebody-CPM-perm_hits.rds")
hist(perms)

perms <- data.frame(perms)
perms$Hits <- perms$perms

p1 <- ggplot(perms, aes(x = Hits)) +
  geom_histogram(bins = 30, 
                 fill = "slategray3", 
                 color = "black", 
                 alpha = 0.8) +
  geom_vline(xintercept = quantile(perms$Hits, 0.95), color = "firebrick3", size = 1) +
  geom_vline(xintercept = 6264, color = "seagreen3", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(title = "Early gestation, whole placenta",
       x = "", y = "Count")

perms <- readRDS("LZ-genebody-CPM-perm_hits.rds")
hist(perms)

perms <- data.frame(perms)
perms$Hits <- perms$perms


p2 <- ggplot(perms, aes(x = Hits)) +
  geom_histogram(bins = 30, 
                 fill = "slategray3", 
                 color = "black", 
                 alpha = 0.8) +
  geom_vline(xintercept = quantile(perms$Hits, 0.95), color = "firebrick3", size = 1) +
  geom_vline(xintercept = 4261, color = "seagreen3", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(title = "Late gestation, labyrinth",
       x = "Significant methylation-expression correlations", y = "")

perms <- readRDS("JZ-genebody-CPM-perm_hits.rds")
hist(perms)

perms <- data.frame(perms)
perms$Hits <- perms$perms

p3 <- ggplot(perms, aes(x = Hits)) +
  geom_histogram(bins = 30, 
                 fill = "slategray3", 
                 color = "black", 
                 alpha = 0.8) +
  geom_vline(xintercept = quantile(perms$Hits, 0.95), color = "firebrick3", size = 1) +
  geom_vline(xintercept = 2367, color = "seagreen3", size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(title = "Late gestation, junctional zone",
       x = "", y = "")

png("FigS6.png", res = 300, width = 15, height = 5, units = "in")
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

pdf("FigS6.pdf", width = 15, height = 5)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

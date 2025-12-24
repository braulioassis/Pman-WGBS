library(ggplot2)
library(gridExtra)

# P1
perms <- readRDS("W-genebody-CPM-perm_hits.rds")
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8)) +
  labs(x = "", y = "Count", title = "Early gestation, whole placenta")

# P2
perms <- readRDS("LZ-genebody-CPM-perm_hits.rds")
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8)) +
  labs(x = "Permuted methylation-expression correlations", y = "", title = "Late gestation, labyrinth")

# P3
perms <- readRDS("JZ-genebody-CPM-perm_hits.rds")
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
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8)) +
  labs(x = "", y = "", title = "Late gestation, junctional zone")

# Render
png("Figure S6.png", res = 300, width = 8.5, height = 2.8, units = "in")
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

pdf("Figure S6.pdf", width = 8.5, height = 2.8)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

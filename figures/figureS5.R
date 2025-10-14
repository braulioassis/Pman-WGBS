library(ggplot2)

# P1
winv <- readRDS("inversions-popdmrCTGA-permutation.rds")
w.inv <- data.frame(winv$numOverlaps$permuted)
mu <- mean(w.inv$winv.numOverlaps.permuted)
sigma <- sd(w.inv$winv.numOverlaps.permuted)

observed <- 68
signif <- quantile(w.inv$winv.numOverlaps.permuted, 0.95)

p1 <- ggplot(w.inv, aes(x = winv.numOverlaps.permuted)) +
  geom_histogram(aes(y = ..count..), 
                 bins = 30, 
                 fill = "slategray3", 
                 alpha = 0.3) +
  geom_vline(xintercept = signif, color = "firebrick3", size = 1) +
  geom_vline(xintercept = observed, color = "seagreen3", size = 1) +
  annotate("text", x = observed - 1, y = 75, 
           label = "", angle = 90, vjust = -0.5, 
           color = "seagreen3", size = 6) +
  annotate("text", x = signif - 1, y = 75, 
           label = "\u03B1 = 0.05", angle = 90, vjust = -0.5, 
           color = "firebrick3", size = 6) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(title = "Early gestation, whole placenta",
       x = "Overlaps", y = "")

# P2
lzinv <- readRDS("inversions-popdmrCTGA-permutation.rds")
lz.inv <- data.frame(lzinv$numOverlaps$permuted)
mu <- mean(lz.inv$lzinv.numOverlaps.permuted)
sigma <- sd(lz.inv$lzinv.numOverlaps.permuted)

observed <- 42
signif <- quantile(lz.inv$lzinv.numOverlaps.permuted, 0.95)
p2 <- ggplot(lz.inv, aes(x = lzinv.numOverlaps.permuted)) +
  geom_histogram(aes(y = ..count..), 
                 bins = 30, 
                 fill = "slategray3", 
                 alpha = 0.3) +
  geom_vline(xintercept = signif, color = "firebrick3", size = 1) +
  geom_vline(xintercept = observed, color = "seagreen3", size = 1) +
  annotate("text", x = observed - 0.2, y = 75, 
           label = "", angle = 90, vjust = -0.5, 
           color = "seagreen3", size = 6) +
  annotate("text", x = signif - 0.2, y = 75, 
           label = "\u03B1 = 0.05", angle = 90, vjust = -0.5, 
           color = "firebrick3", size = 6) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(title = "Late gestation, labyrinth",
       x = "Overlaps", y = "")

# P3
jzinv <- readRDS("inversions-popdmrCTGA-permutation.rds")
jz.inv <- data.frame(jzinv$numOverlaps$permuted)
mu <- mean(jz.inv$jzinv.numOverlaps.permuted)
sigma <- sd(jz.inv$jzinv.numOverlaps.permuted)

observed <- 36
signif <- quantile(jz.inv$jzinv.numOverlaps.permuted, 0.95)
p3 <- ggplot(jz.inv, aes(x = jzinv.numOverlaps.permuted)) +
  geom_histogram(aes(y = ..count..), 
                 bins = 30, 
                 fill = "slategray3", 
                 alpha = 0.3) +
  geom_vline(xintercept = signif, color = "firebrick3", size = 1) +
  geom_vline(xintercept = observed, color = "seagreen3", size = 1) +
  annotate("text", x = observed - 0.2, y = 75, 
           label = "", angle = 90, vjust = -0.5, 
           color = "seagreen3", size = 6) +
  annotate("text", x = signif - 0.2, y = 75, 
           label = "\u03B1 = 0.05", angle = 90, vjust = -0.5, 
           color = "firebrick3", size = 6) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(title = "Late gestation, junctional zone",
       x = "Overlaps", y = "")


png("FigS5.png", res = 300, width = 15, height = 5, units = "in")
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

pdf("FigS5.pdf", width = 15, height = 5)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

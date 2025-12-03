library(pheatmap)

# W
df <- read.csv("w.summary.hypoxia.regions.csv")

colnames(df) <- gsub("_CpG", "", colnames(df))

sample_info <- read.csv("w.summary.CpG.csv")
sample_info <- sample_info[sample_info$Population == "ME" | sample_info$Population == "BW", ]
sample_info$Treatment[sample_info$Treatment == "1N"] <- "Normoxia"
sample_info$Treatment[sample_info$Treatment == "2H"] <- "Hypoxia"

sample_info$Group <- with(sample_info,
                          interaction(Population, Treatment, sep = "_"))

sample_info$Group <- factor(sample_info$Group,
                            levels = c("BW_Normoxia", "BW_Hypoxia",
                                       "ME_Normoxia", "ME_Hypoxia")
)

meth_matrix <- as.matrix(df[ , -(1:3)])
rownames(meth_matrix) <- paste(df$Contig, df$Start, df$End, sep = "_")

new_order <- order(sample_info$Group)

meth_matrix2 <- meth_matrix[, new_order]
sample_info2 <- sample_info[new_order, ]

meth_matrix2 <- t(apply(meth_matrix2, 1, function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

meth_matrix2 <- meth_matrix2[rowSums(is.na(meth_matrix2)) == 0, ]

row_var <- apply(meth_matrix2, 1, var)
meth_matrix2 <- meth_matrix2[row_var > 0, ]

annotation <- data.frame(Group = sample_info2$Group)
rownames(annotation) <- colnames(meth_matrix2)

png("w.heatmap.png",
    res = 300, width = 6, height = 6, units = "in")
pheatmap(
  meth_matrix2,
  cluster_cols = F,
  cluster_rows = T,
  annotation_col = annotation,
  annotation_names_col = T,
  show_rownames = F,
  show_colnames = F,
  main = "Methylation rate, hypoxia-sensitive DMRs\nEarly gestation, whole placenta"
)
dev.off()

# LZ
df <- read.csv("lz.summary.hypoxia.regions.csv")

colnames(df) <- gsub("_CpG", "", colnames(df))

sample_info <- read.csv("lz.summary.CpG.csv")
sample_info <- sample_info[sample_info$Population == "ME" | sample_info$Population == "BW", ]
sample_info$Treatment[sample_info$Treatment == "1N"] <- "Normoxia"
sample_info$Treatment[sample_info$Treatment == "2H"] <- "Hypoxia"

sample_info$Group <- with(sample_info,
                          interaction(Population, Treatment, sep = "_"))

sample_info$Group <- factor(sample_info$Group,
                            levels = c("BW_Normoxia", "BW_Hypoxia",
                                       "ME_Normoxia", "ME_Hypoxia")
)

meth_matrix <- as.matrix(df[ , -(1:3)])
rownames(meth_matrix) <- paste(df$Contig, df$Start, df$End, sep = "_")

new_order <- order(sample_info$Group)

meth_matrix2 <- meth_matrix[, new_order]
sample_info2 <- sample_info[new_order, ]

meth_matrix2 <- t(apply(meth_matrix2, 1, function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

meth_matrix2 <- meth_matrix2[rowSums(is.na(meth_matrix2)) == 0, ]

row_var <- apply(meth_matrix2, 1, var)
meth_matrix2 <- meth_matrix2[row_var > 0, ]

annotation <- data.frame(Group = sample_info2$Group)
rownames(annotation) <- colnames(meth_matrix2)

png("lz.heatmap.png",
    res = 300, width = 6, height = 6, units = "in")
pheatmap(
  meth_matrix2,
  cluster_cols = F,
  cluster_rows = T,
  annotation_col = annotation,
  annotation_names_col = T,
  show_rownames = F,
  show_colnames = F,
  main = "Methylation rate, hypoxia-sensitive DMRs\nLate gestation, labyrinth"
)
dev.off()

# JZ
df <- read.csv("jz.summary.hypoxia.regions.csv")

colnames(df) <- gsub("_CpG", "", colnames(df))

sample_info <- read.csv("jz.summary.CpG.csv")
sample_info <- sample_info[sample_info$Population == "ME" | sample_info$Population == "BW", ]
sample_info$Treatment[sample_info$Treatment == "1N"] <- "Normoxia"
sample_info$Treatment[sample_info$Treatment == "2H"] <- "Hypoxia"

sample_info$Group <- with(sample_info,
                          interaction(Population, Treatment, sep = "_"))

sample_info$Group <- factor(sample_info$Group,
                            levels = c("BW_Normoxia", "BW_Hypoxia",
                                       "ME_Normoxia", "ME_Hypoxia")
)

meth_matrix <- as.matrix(df[ , -(1:3)])
rownames(meth_matrix) <- paste(df$Contig, df$Start, df$End, sep = "_")

new_order <- order(sample_info$Group)

meth_matrix2 <- meth_matrix[, new_order]
sample_info2 <- sample_info[new_order, ]

meth_matrix2 <- t(apply(meth_matrix2, 1, function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

meth_matrix2 <- meth_matrix2[rowSums(is.na(meth_matrix2)) == 0, ]

row_var <- apply(meth_matrix2, 1, var)
meth_matrix2 <- meth_matrix2[row_var > 0, ]

annotation <- data.frame(Group = sample_info2$Group)
rownames(annotation) <- colnames(meth_matrix2)

png("jz.heatmap.png",
    res = 300, width = 6, height = 6, units = "in")
pheatmap(
  meth_matrix2,
  cluster_cols = F,
  cluster_rows = T,
  annotation_col = annotation,
  annotation_names_col = T,
  show_rownames = F,
  show_colnames = F,
  main = "Methylation rate, hypoxia-sensitive DMRs\nLate gestation, junctional zone"
)
dev.off()

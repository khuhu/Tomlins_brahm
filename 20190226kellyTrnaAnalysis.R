library(pheatmap)

tRnaDat <- read_xlsx("/home/kevhu/data/tRNA_Mod_Analysis_Dec2018.xlsx", sheet = 5,skip = 8, col_names = FALSE)

raw_trna <- tRnaDat[1:42,1:26]
normalized_trna <- tRnaDat[1:42,c(1,2,27:50)]

colnames(raw_trna) <- raw_trna[1,]
raw_trna <- raw_trna[-1,]

colnames(normalized_trna) <- normalized_trna[1,]
normalized_trna <- normalized_trna[-1,]

raw_trna <- raw_trna[-which(raw_trna$HEK1 == "ND"),]
normalized_trna <- normalized_trna[-which(normalized_trna$HEK1 == "ND"),]

### 20190227: below not using because doing averages of triplicaes and log2 ratios
group1 <- c("H358 Puro", "H358 TPRKB shJ17", "H358 TPRKB sh46")

g1Idx <- NULL
for (i in group1) {
  tmpIdx <- grep(i,colnames(normalized_trna) , fixed = FALSE)
  g1Idx <- c(g1Idx, tmpIdx)
}

g1Idx <- c(1, g1Idx)

tmpSeq <- seq_along(colnames(normalized_trna))
g2Idx <- tmpSeq[-g1Idx]
g2Idx[1] <- 1 
  
g3Idx <- seq_along(colnames(normalized_trna))
g3Idx <- g3Idx[-2]


normalized_trna_g1 <- normalized_trna[,g1Idx]
normalized_trna_g2 <- normalized_trna[,g2Idx]
normalized_trna_g3 <- normalized_trna[,g3Idx]

normalized_trna_g1_mat <- as.matrix(data.frame(lapply(normalized_trna_g1[,2:ncol(normalized_trna_g1)], as.numeric)))
rownames(normalized_trna_g1_mat) <- normalized_trna_g1$Nucleoside
normalized_trna_g1_mat_scaled <- log10(normalized_trna_g1_mat)
normalized_trna_g1_mat_scaled <- apply(normalized_trna_g1_mat_scaled, 1, scale)
rownames(normalized_trna_g1_mat_scaled) <- colnames(normalized_trna_g1)[2:ncol(normalized_trna_g1)]

heatMapCol <- colorRampPalette(c("blue","black","yellow"))(1000)
quantile(normalized_trna_g1_mat_scaled, probs = seq(0, 1, 0.01))
normalized_trna_g1_mat_scaled[normalized_trna_g1_mat_scaled > 2.1] <- 2.1
colors.breaks <- seq(-2.1,2.1, 4.2/1000)

pdf("/home/kevhu/data/20190226group1Heatmap.pdf", useDingbats = TRUE)
pheatmap(t(normalized_trna_g1_mat_scaled), color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 7,
         breaks = colors.breaks,border_color = "black",
         cellwidth = 14, cellheight = 9)
dev.off()


### group 2

normalized_trna_g2_mat <- as.matrix(data.frame(lapply(normalized_trna_g2[,2:ncol(normalized_trna_g2)], as.numeric)))
rownames(normalized_trna_g2_mat) <- normalized_trna_g2$Nucleoside
normalized_trna_g2_mat_scaled <- log10(normalized_trna_g2_mat)
normalized_trna_g2_mat_scaled <- apply(normalized_trna_g2_mat_scaled, 1, scale)
rownames(normalized_trna_g2_mat_scaled) <- colnames(normalized_trna_g2)[2:ncol(normalized_trna_g2)]

quantile(normalized_trna_g2_mat_scaled, probs = seq(0, 1, 0.01))
normalized_trna_g2_mat_scaled[normalized_trna_g2_mat_scaled > 3] <- 3
colors.breaks2 <- seq(-3,3, 6/1000)

pdf("/home/kevhu/data/20190226group2Heatmap.pdf", useDingbats = TRUE)
pheatmap(t(normalized_trna_g2_mat_scaled), color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 7,
         breaks = colors.breaks2,border_color = "black",
         cellwidth = 14, cellheight = 9)
dev.off()

### group 3 i.e all
normalized_trna_g3_mat <- as.matrix(data.frame(lapply(normalized_trna_g3[,2:ncol(normalized_trna_g3)], as.numeric)))
rownames(normalized_trna_g3_mat) <- normalized_trna_g3$Nucleoside
normalized_trna_g3_mat_scaled <- log10(normalized_trna_g3_mat)
normalized_trna_g3_mat_scaled <- apply(normalized_trna_g3_mat_scaled, 1, scale)
rownames(normalized_trna_g3_mat_scaled) <- colnames(normalized_trna_g3)[2:ncol(normalized_trna_g3)]

quantile(normalized_trna_g3_mat_scaled, probs = seq(0, 1, 0.01))
colors.breaks3 <- seq(-3,3, 6/1000)

pdf("/home/kevhu/data/20190226allSamplesHeatmap.pdf", useDingbats = TRUE, width = 10)
pheatmap(t(normalized_trna_g3_mat_scaled), color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 7,
         breaks = colors.breaks3,border_color = "black",
         cellwidth = 14, cellheight = 9)
dev.off()




### below shoould be the final version
###
###

normalized_trna2 <- data.frame(lapply(normalized_trna[,3:ncol(normalized_trna)], as.numeric))

normalized_trna3 <- NULL
for (i in seq(3, ncol(normalized_trna2), 3)) {
  tmpRow <- rowMeans(normalized_trna2[,c(i - 2, i - 1 , i)])
  normalized_trna3 <- rbind(normalized_trna3, tmpRow)
}

normalized_trna3 <- t(normalized_trna3)

tmpNames <- colnames(normalized_trna2)[seq(3, ncol(normalized_trna2), 3)]
colnames(normalized_trna3) <- substr(tmpNames, 1, nchar(tmpNames) - 1)
rownames(normalized_trna3) <- normalized_trna$Nucleoside

group1 <- c("H358.Puro", "H358.TPRKB.shJ17", "H358.TPRKB.sh46")

g1Idx <- NULL
for (i in group1) {
  tmpIdx <- grep(i,colnames(normalized_trna3) , fixed = FALSE)
  g1Idx <- c(g1Idx, tmpIdx)
}

tmpSeq <- seq_along(colnames(normalized_trna3))
g2Idx <- tmpSeq[-g1Idx]

normalized_trna3.g1 <- data.frame(normalized_trna3[,g1Idx])
normalized_trna3.g2 <- data.frame(normalized_trna3[,g2Idx])

normalized_trna3.combinedRatio <- data.frame(normalized_trna3.g1[,2:3]/normalized_trna3.g1[,1],
           normalized_trna3.g2[,2:3]/normalized_trna3.g2[,1],
           normalized_trna3.g2[,5]/normalized_trna3.g2[,4])

colnames(normalized_trna3.combinedRatio)[5] <- colnames(normalized_trna3.g2)[5]
normalized_trna3.combinedRatio.log2 <- log2(normalized_trna3.combinedRatio)



quantile(as.matrix(normalized_trna3.combinedRatio.log2), probs = seq(0, 1, 0.01))
colors.breaks4 <- seq(-1,1, 2/1000)
heatMapCol <- colorRampPalette(c("blue","black","yellow"))(1000)

pdf("/home/kevhu/data/20190227allSamplesHeatmap.pdf", useDingbats = TRUE)
pheatmap(normalized_trna3.combinedRatio.log2, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 7,
         breaks = colors.breaks4,border_color = "black",
         cellwidth = 30, cellheight = 10)
dev.off()


normalized_trna3.combinedRatio.log2_g1 <- normalized_trna3.combinedRatio.log2[,1:2]
normalized_trna3.combinedRatio.log2_g2 <- normalized_trna3.combinedRatio.log2[,3:5]

pdf("/home/kevhu/data/20190227shSamplesHeatmap.pdf", useDingbats = TRUE)
pheatmap(normalized_trna3.combinedRatio.log2_g1, color = heatMapCol,
         cluster_cols = FALSE, cluster_rows = TRUE, fontsize = 7,
         breaks = colors.breaks4,border_color = "black",
         cellwidth = 30, cellheight = 10)
dev.off()

pdf("/home/kevhu/data/20190227koSamplesHeatmap.pdf", useDingbats = TRUE)
pheatmap(normalized_trna3.combinedRatio.log2_g2, color = heatMapCol,
         cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 7,
         breaks = colors.breaks4,border_color = "black",
         cellwidth = 30, cellheight = 10)
dev.off()




library(argyle)
library(data.table)
load("/br_z1/kevin_storage/misc/snps.gigamuga.Rdata")
load("/br_z1/kevin_storage/misc/clusters.gigamuga.Rdata")

# data(ex)
# qcplot(ex)

snpRes <- read.beadstudio(prefix = "",
                          in.path = "/br_z1/kevin_storage/advancedGenomicsCore/tmp/",
                          snps = snps)

normedSnpRes <- tQN(snpRes, adjust.lrr = FALSE, clusters = clusters)
baf_normed <- get.baf(normedSnpRes)
baf_mat <- dcast(baf_normed, marker + chr + pos ~ iid, value.var = c("BAF"))
logR_mat <- dcast(baf_normed, marker + chr + pos ~ iid, value.var = c("LRR"))
baf_mat$chr <- stringr::str_remove(baf_mat$chr, "chr")
logR_mat$chr <- stringr::str_remove(logR_mat$chr, "chr")

tmp <- tQN(snpRes, clusters = clusters)
baf_normed2 <- get.baf(tmp)
baf_mat2 <- dcast(baf_normed2, marker + chr + pos ~ iid, value.var = c("BAF"))
logR_mat2 <- dcast(baf_normed2, marker + chr + pos ~ iid, value.var = c("LRR"))
baf_mat2$chr <- stringr::str_remove(baf_mat2$chr, "chr")
logR_mat2$chr <- stringr::str_remove(logR_mat2$chr, "chr")

noNormSnpRes  <- tQN(snpRes, adjust.lrr = TRUE, clusters = clusters, prenorm = FALSE, xynorm = FALSE)
baf_normed3 <- get.baf(noNormSnpRes)
baf_mat3 <- dcast(baf_normed3, marker + chr + pos ~ iid, value.var = c("BAF"))
logR_mat3 <- dcast(baf_normed3, marker + chr + pos ~ iid, value.var = c("LRR"))
baf_mat3$chr <- stringr::str_remove(baf_mat3$chr, "chr")
logR_mat3$chr <- stringr::str_remove(logR_mat3$chr, "chr")


### different scaling factors - doesn't change anything between same within and between params
###
###

# scaledSnpRes <- tQN(snpRes, adjust.lrr = TRUE, clusters = clusters, thresholds = c(2.5, 2.5))
# baf_scaled2.5 <- get.baf(scaledSnpRes)
# baf_mat_s_2.5 <- dcast(baf_scaled2.5, marker + chr + pos ~ iid, value.var = c("BAF"))
# logR_mat_s_2.5 <- dcast(baf_scaled2.5, marker + chr + pos ~ iid, value.var = c("LRR"))
# baf_mat_s_2.5$chr <- stringr::str_remove(baf_mat_s_2.5$chr, "chr")
# logR_mat_s_2.5$chr <- stringr::str_remove(logR_mat_s_2.5$chr, "chr")
# 
# 
# scaledSnpRes3 <- tQN(snpRes, adjust.lrr = TRUE, clusters = clusters, thresholds = c(3, 3))
# baf_scaled3 <- get.baf(scaledSnpRes3)
# baf_mat_s_3 <- dcast(baf_scaled3, marker + chr + pos ~ iid, value.var = c("BAF"))
# logR_mat_s_3 <- dcast(baf_scaled3, marker + chr + pos ~ iid, value.var = c("LRR"))
# baf_mat_s_3$chr <- stringr::str_remove(baf_mat_s_3$chr, "chr")
# logR_mat_s_3$chr <- stringr::str_remove(logR_mat_s_3$chr, "chr")
# 
# 
# scaledSnpRes3.5 <- tQN(snpRes, adjust.lrr = TRUE, clusters = clusters, thresholds = c(3.5, 3.5))
# baf_scaled3.5 <- get.baf(scaledSnpRes3.5)
# baf_mat_s_3.5 <- dcast(baf_scaled3.5, marker + chr + pos ~ iid, value.var = c("BAF"))
# logR_mat_s_3.5 <- dcast(baf_scaled3.5, marker + chr + pos ~ iid, value.var = c("LRR"))
# baf_mat_s_3.5$chr <- stringr::str_remove(baf_mat_s_3.5$chr, "chr")
# logR_mat_s_3.5$chr <- stringr::str_remove(logR_mat_s_3.5$chr, "chr")




write.table(test_ascat_tumor.logR, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(test_ascat_tumor.baf, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

test_ascat_tumor.logR2 <- logR_mat2
test_ascat_tumor.baf2 <- baf_mat2

write.table(test_ascat_tumor.logR2, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR2.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(test_ascat_tumor.baf2, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf2.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


test_ascat_tumor.logR3 <- logR_mat3
test_ascat_tumor.baf3 <- baf_mat3

write.table(test_ascat_tumor.logR3, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(test_ascat_tumor.baf3, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


library(ASCAT)
source("/home/kevhu/scripts/20220331ascat.predictGermlineGenotypes.R")

ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/unadjLrr/")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/unadjLrr/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/unadjLrr/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/unadjLrr/")

save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220405ascsatbcLrrUn.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220405ascsatoutLrrUn.rds")



ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR2.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf2.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/adjLrr/")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/adjLrr/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/adjLrr/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/adjLrr/")


save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220405ascsatbcLrrAdj.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220405ascsatoutLrrAdj.rds")



ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/noNorm/")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/noNorm/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/noNorm/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/noNorm/")


save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220405ascsatbcLrrAdjNoNorm.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220405ascsatoutLrrAdjNoNorm.rds")







nameStripper <- function(df){
  require(stringr)
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}

tcDf <- read.table("/br_z1/kevin_storage/misc/20220321hgscTc.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
sampleMap <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/Sample_Map.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)
sampleMap$stripped <- nameStripper(sampleMap$ID)
sampleMap$stripped[21] <- "13085lt"
sampleMap$stripped[32] <- "14154lt"
# sampleMap$stripped[42] <- "15774rt"
# sampleMap$stripped[43] <- "15774lt"

library(foreach)
library(doParallel)


runGridAscat <- function(i, sampleMap, tcDf, outdir, gamma  = 0.15, illu = NULL){
  setwd(outdir)
  source("/home/kevhu/scripts/20220404ascatCnBafCustom.R", local = TRUE)
  
  if (is.null(illu)) {
    matchedName <- sampleMap$stripped[which(sampleMap$Name == colnames(ascat.bc[[1]])[i])]
  } else{
    matchedName <- sampleMap$stripped[which(sampleMap$stripped == colnames(ascat.bc[[1]])[i])]
    if (length(matchedName) == 0) {
      matchedName <- sampleMap$stripped[which(sampleMap$ID == colnames(ascat.bc[[1]])[i])]
    }
  }
  
  ### hard coding sample mismatch
  if (matchedName == "15774lt" | matchedName == "15774rt") {
    if (matchedName == "15774lt") {
      matchedName <- "15774rt"
    } else if (matchedName == "15774rt") {
      matchedName <- "15774lt"
    }
  }
  
  if (length(tcDf$tc[which(tcDf$sample == matchedName)]) == 0) {
    tmpTc <- 0.96
  } else{
    tmpTc <- signif(tcDf$tc[which(tcDf$sample == matchedName)][1], 2)
  }

  
  
  if (tmpTc > 0.94) {
    purityVec <- seq(0.88, 1, 0.04)
  } else{
    purityVec <- seq(tmpTc - 0.02 * 3, tmpTc + 0.02 * 3, 0.04)
  }
  
  ploidyVec <- seq(1.6, 5.4, 0.4)
  ploidyVec2 <- rep(ploidyVec, 4)
  purityVec2 <- NULL
  for (j in purityVec) {
    purityVec2 <- c(purityVec2, rep(j, 10))
  }
  
  
  tmp.bc <- ascat.bc
  tmp.bc[[1]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[1]][,i]), 40))
  tmp.bc[[2]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[2]][,i]), 40))
  tmp.bc[[7]] <- paste0(names(ascat.bc[[1]][i]), "_rho", purityVec2, "_psi", ploidyVec2)
  tmp.bc[[8]] <- rep("XX", 40)
  tmp.bc[[11]] <- rep(unlist(ascat.bc[[11]][i]), 40)
  tmp.bc[[12]] <- rep(unlist(ascat.bc[[12]][i]), 40)
  tmp.bc[[13]] <- rep(unlist(ascat.bc[[13]][i]), 40)
  tmp.bc[[14]] <- rep(unlist(ascat.bc[[14]][i]), 40)
  tmp.bc[[15]] <- Rfast::rep_col(unlist(ascat.bc[[15]][,i]), 40)
  for (j in 1:40) {
    tmp.bc[[16]][[j]] <- ascat.bc[[16]][[i]]
  }
  
  colnames(tmp.bc[[1]]) <- tmp.bc[[7]]
  colnames(tmp.bc[[2]]) <- tmp.bc[[7]]
  names(tmp.bc[[11]]) <- tmp.bc[[7]]
  names(tmp.bc[[12]]) <- tmp.bc[[7]]
  names(tmp.bc[[13]]) <- tmp.bc[[7]]
  names(tmp.bc[[14]]) <- tmp.bc[[7]]
  colnames(tmp.bc[[15]]) <- tmp.bc[[7]]
  
  rownames(tmp.bc[[1]]) <- rownames(ascat.bc[[1]])
  rownames(tmp.bc[[2]]) <- rownames(ascat.bc[[2]])
  rownames(tmp.bc[[15]]) <- rownames(ascat.bc[[15]])
  if (!file.exists(paste0("./", names(ascat.bc[[1]][i])))) {
    dir.create(paste0("./", names(ascat.bc[[1]][i])))
  }
  setwd(paste0("./", names(ascat.bc[[1]][i])))
  ascat.output = ascat.runAscat(tmp.bc, rho_manual = purityVec2, psi_manual = ploidyVec2, gamma = gamma)
  ascatDf <- data.frame("sample" = names(ascat.output$ploidy), "ploidy" = ascat.output$ploidy, 
                        "purity" = ascat.output$purity, "goodnessOfFit" = ascat.output$goodnessOfFit)
  ascatDf <- ascatDf[order(ascatDf$goodnessOfFit, decreasing = TRUE), ]
  ascatDf <- ascatDf[1:15,]
  rownames(ascatDf) <- NULL
  
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpMatch <- match(rownames(tmpBafDf), rownames(ascat.bc$Tumor_BAF_segmented[[i]]))
  tmpMatch <- which(!is.na(tmpMatch))
  tmpBafDf$AF[-tmpMatch] <- NA
  tmpBafGraph <- freqPlot_baf(tmpBafDf)
  tmpBafGrob <- ggplotGrob(tmpBafGraph)
  # ascat_cn(df = ascat.output$segments)
  ascat_cn2(df = ascat.output$segments)
  
  write.table(ascat.output$segments, file = paste0(outdir, "allSegResults.txt"), sep = "\t", quote = FALSE,
              col.names = TRUE, row.names = FALSE, append = TRUE)
  
  ascatDf
}




load("/br_z1/kevin_storage/ASCAT/20220405ascsatbcLrrUn.rds")
load("/br_z1/kevin_storage/ASCAT/20220405ascsatoutLrrUn.rds")
no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/unadjLrr/")
stopCluster(cl)
unadjRes <- res



load("/br_z1/kevin_storage/ASCAT/20220405ascsatbcLrrAdj.rds")
load("/br_z1/kevin_storage/ASCAT/20220405ascsatoutLrrAdj.rds")
no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/adjLrr/")
stopCluster(cl)

adjRes <- res


### look at the homo snps again ... and see where those homozygous snps on the edges lie  within all samples


tmpDev <- NULL
for (i in 1:48) {
  tmpDev <- c(tmpDev, rownames(ascat.bc$Tumor_BAF_segmented[[i]])[which(ascat.bc$Tumor_BAF_segmented[[i]] < 0.35)])
}

possibleArtifacts05 <- names(table(tmpDev))[which(table(tmpDev) > 3)]
possibleArtifacts10 <- names(table(tmpDev))[which(table(tmpDev) > 5)]
possibleArtifacts20 <- names(table(tmpDev))[which(table(tmpDev) > 10)]

bafData <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorBaf.txt", sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
bafHom <- bafData[,4:ncol(bafData)]
meanBafHom2 <- apply(bafHom, 1, mean)
sdBafHom2 <- apply(bafHom, 1, sd)
plot(meanBafHom2, sdBafHom2)

sd30mean30 <- rownames(ascat.bc$SNPpos)[intersect(which(meanBafHom2 > 0.7 | meanBafHom2 < 0.3),
                                                  which(sdBafHom2 < 0.3))]
sd20mean30 <- rownames(ascat.bc$SNPpos)[intersect(which(meanBafHom2 > 0.7 | meanBafHom2 < 0.3),
                                                  which(sdBafHom2 < 0.2))]

for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% possibleArtifacts20)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% possibleArtifacts20)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpArt20/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% possibleArtifacts10)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% possibleArtifacts10)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpArt10/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}


for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% possibleArtifacts05)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% possibleArtifacts05)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpArt05/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% sd30mean30)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% sd30mean30)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/sd30mean30/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}


for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% sd20mean30)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% sd20mean30)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/sd20mean30/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}


### trying new baf calling from non-normalized snp data

load("/br_z1/kevin_storage/ASCAT/20220405ascsatbcLrrAdjNoNorm.rds")
load("/br_z1/kevin_storage/ASCAT/20220405ascsatoutLrrAdjNoNorm.rds")


oldBafData <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf.txt", sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)
bafData <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt", sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
bafHom <- bafData[,4:ncol(bafData)]
meanBafHom2 <- apply(bafHom, 1, mean)
sdBafHom2 <- apply(bafHom, 1, sd)
plot(meanBafHom2, sdBafHom2)

sd30mean30 <- rownames(ascat.bc$SNPpos)[intersect(which(meanBafHom2 > 0.7 | meanBafHom2 < 0.3),
                                                  which(sdBafHom2 < 0.3))]
sd20mean30 <- rownames(ascat.bc$SNPpos)[intersect(which(meanBafHom2 > 0.7 | meanBafHom2 < 0.3),
                                                  which(sdBafHom2 < 0.2))]


tmpDev <- NULL
for (i in 1:48) {
  tmpDev <- c(tmpDev, rownames(ascat.bc$Tumor_BAF_segmented[[i]])[which(ascat.bc$Tumor_BAF_segmented[[i]] < 0.35)])
}

possibleArtifacts05 <- names(table(tmpDev))[which(table(tmpDev) > 3)]
possibleArtifacts10 <- names(table(tmpDev))[which(table(tmpDev) > 5)]
possibleArtifacts20 <- names(table(tmpDev))[which(table(tmpDev) > 10)]



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% sd20mean30)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% sd20mean30)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/noNormSnpsSd20Mean30/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% sd30mean30)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% sd20mean30)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/noNormSnpsSd30Mean30/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% possibleArtifacts20)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% possibleArtifacts20)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/noNormSnpArt20/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% possibleArtifacts10)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% possibleArtifacts10)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/noNormSnpArt10/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}


for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf$AF[which(tmpBafDf$AF > .90 | tmpBafDf$AF < .10 )] <- NA
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% possibleArtifacts05)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% possibleArtifacts05)] <- NA
  # tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/noNormSnpArt05/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}

save(sd30mean30, file = "/br_z1/kevin_storage/ASCAT/20220410noNorm3030MarkersToFilt.rds")


###  filtering prior with no norm

load("/br_z1/kevin_storage/ASCAT/20220410noNorm3030MarkersToFilt.rds")
tmpTable <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
tmpTable[which(tmpTable$marker %in% sd30mean30), 4:ncol(tmpTable)] <- NA
write.table(tmpTable, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3Filted3030.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3Filted3030.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat10/")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat10/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat10/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat10/")


save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220410ascsatbcLrrAdjNoNorm303010.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220410ascsatoutLrrAdjNoNorm303010.rds")




ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3Filted3030.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat15/")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat15/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat15/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat15/")


save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220410ascsatbcLrrAdjNoNorm303015.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220410ascsatoutLrrAdjNoNorm303015.rds")

### doing custom ascat on new filters


load("/br_z1/kevin_storage/ASCAT/20220410ascsatbcLrrAdjNoNorm303010.rds")
load("/br_z1/kevin_storage/ASCAT/20220410ascsatoutLrrAdjNoNorm303010.rds")
no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf,"/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat10/")
stopCluster(cl)
res10 <- res

write.table(res10, "/br_z1/kevin_storage/ASCAT/20220410ascatTableNoNormRes10.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)


load("/br_z1/kevin_storage/ASCAT/20220410ascsatbcLrrAdjNoNorm303015.rds")
load("/br_z1/kevin_storage/ASCAT/20220410ascsatoutLrrAdjNoNorm303015.rds")
no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/noNorm3030FiltAscat15/")
stopCluster(cl)
res15 <- res
write.table(res15, "/br_z1/kevin_storage/ASCAT/20220410ascatTableNoNormRes15.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)




tmp <- cbind(ascat.bc$SNPpos, ascat.bc$Tumor_LogR)
mean(tmp$`5516-KH-47`[which(tmp$chr == 19)], na.rm = TRUE)


### no norm, the snps looks better as the het snps revolve more around 50% 
### testing different gamma

ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/tmp/")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/tmp/", ascat.gg = ascat.gg, penalty = 50)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/tmp/")

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/tmp/")
stopCluster(cl)
res


### instead of multiplying the logRs can I just alter gamma and penalty ...l
### make sure to only use smaller sample size ... i.e look at segmentation  group plots
library(foreach)
library(doParallel)
library(ASCAT)
source("/home/kevhu/scripts/20220331ascat.predictGermlineGenotypes.R")

tcDf <- read.table("/br_z1/kevin_storage/misc/20220321hgscTc.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
load("/br_z1/kevin_storage/ASCAT/20220402ascatbcRes24NoMatched.rds")
sampleMap <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/Sample_Map.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)
sampleMap$stripped <- nameStripper(sampleMap$ID)
sampleMap$stripped[21] <- "13085lt"
sampleMap$stripped[32] <- "14154lt"


tmpSampleMap <- sampleMap
tmpSampleMap$tc  <- tcDf$tc[match(tmpSampleMap$ID, tcDf$sample)]

tmp.logr <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt", sep = "\t",
                     header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
tmp.baf <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

illu.logr <- read.table("/br_z1/kevin_storage/ASCAT/20220410illuminaLogRr.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
illu.baf <- read.table("/br_z1/kevin_storage/ASCAT/20220410illuminaBaf.txt", sep = "\t",
                      header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

setOfSamples <- c("5516-KH-1", "5516-KH-13", "5516-KH-45", "5516-KH-47", "5516-KH-48")
setOfSamples2 <- c("1628lt", "6786lt", "kc01", "kc07", "kc10")

tmp.logr <- cbind(tmp.logr[,1:3], tmp.logr[, setOfSamples])
tmp.baf <- cbind(tmp.baf[,1:3], tmp.baf[, setOfSamples])


illu.logr <- cbind(illu.logr[,1:3], illu.logr[, setOfSamples2])
illu.baf <- cbind(illu.baf[,1:3], illu.baf[, setOfSamples2])


write.table(tmp.logr, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411TumorLogrSmall.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(tmp.baf, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411tumorBafSmall.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(illu.logr, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220418IlluTumorLogrSmall.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(illu.baf, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220418IlluTumorBafSmall.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)



runGridAscatGam <- function(i, sampleMap, tcDf, gam){
  source("/home/kevhu/scripts/20220404ascatCnBafCustom.R", local = TRUE)
  
  matchedName <- sampleMap$stripped[which(sampleMap$Name == colnames(ascat.bc[[1]])[i])]
  if (length(tcDf$tc[which(tcDf$sample == matchedName)]) == 0) {
    tmpTc <- 0.96
  } else{
    tmpTc <- signif(tcDf$tc[which(tcDf$sample == matchedName)][1], 2)
  }
  
  if (tmpTc > 0.94) {
    purityVec <- seq(0.88, 1, 0.04)
  } else{
    purityVec <- seq(tmpTc - 0.02 * 3, tmpTc + 0.02 * 3, 0.04)
  }
  
  ploidyVec <- seq(1.6, 5.4, 0.4)
  ploidyVec2 <- rep(ploidyVec, 4)
  purityVec2 <- NULL
  for (j in purityVec) {
    purityVec2 <- c(purityVec2, rep(j, 10))
  }
  
  
  tmp.bc <- ascat.bc
  tmp.bc[[1]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[1]][,i]), 40))
  tmp.bc[[2]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[2]][,i]), 40))
  tmp.bc[[7]] <- paste0(names(ascat.bc[[1]][i]), "_rho", purityVec2, "_psi", ploidyVec2)
  tmp.bc[[8]] <- rep("XX", 40)
  tmp.bc[[11]] <- rep(unlist(ascat.bc[[11]][i]), 40)
  tmp.bc[[12]] <- rep(unlist(ascat.bc[[12]][i]), 40)
  tmp.bc[[13]] <- rep(unlist(ascat.bc[[13]][i]), 40)
  tmp.bc[[14]] <- rep(unlist(ascat.bc[[14]][i]), 40)
  tmp.bc[[15]] <- Rfast::rep_col(unlist(ascat.bc[[15]][,i]), 40)
  for (j in 1:40) {
    tmp.bc[[16]][[j]] <- ascat.bc[[16]][[i]]
  }
  
  colnames(tmp.bc[[1]]) <- tmp.bc[[7]]
  colnames(tmp.bc[[2]]) <- tmp.bc[[7]]
  names(tmp.bc[[11]]) <- tmp.bc[[7]]
  names(tmp.bc[[12]]) <- tmp.bc[[7]]
  names(tmp.bc[[13]]) <- tmp.bc[[7]]
  names(tmp.bc[[14]]) <- tmp.bc[[7]]
  colnames(tmp.bc[[15]]) <- tmp.bc[[7]]
  
  rownames(tmp.bc[[1]]) <- rownames(ascat.bc[[1]])
  rownames(tmp.bc[[2]]) <- rownames(ascat.bc[[2]])
  rownames(tmp.bc[[15]]) <- rownames(ascat.bc[[15]])
  if (!file.exists(paste0("./", names(ascat.bc[[1]][i])))) {
    dir.create(paste0("./", names(ascat.bc[[1]][i])))
  }
  setwd(paste0("./", names(ascat.bc[[1]][i])))
  ascat.output = ascat.runAscat(tmp.bc, rho_manual = purityVec2, psi_manual = ploidyVec2, gamma = gam)
  ascatDf <- data.frame("sample" = names(ascat.output$ploidy), "ploidy" = ascat.output$ploidy, 
                        "purity" = ascat.output$purity, "goodnessOfFit" = ascat.output$goodnessOfFit)
  ascatDf <- ascatDf[order(ascatDf$goodnessOfFit, decreasing = TRUE), ]
  ascatDf <- ascatDf[1:15,]
  rownames(ascatDf) <- NULL
  
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpMatch <- match(rownames(tmpBafDf), rownames(ascat.bc$Tumor_BAF_segmented[[i]]))
  tmpMatch <- which(!is.na(tmpMatch))
  tmpBafDf$AF[-tmpMatch] <- NA
  tmpBafGraph <- freqPlot_baf(tmpBafDf)
  tmpBafGrob <- ggplotGrob(tmpBafGraph)
  # ascat_cn(df = ascat.output$segments, prefix = gam)
  ascat_cn2(df = ascat.output$segments, prefix = gam)
  
  ascatDf
}



listOfPenalty <- seq(10, 100, 10)
listOfGamma <- 0.15
for (q in listOfPenalty) {
  dir <- "/br_z1/kevin_storage/ASCAT/testPenaltyGamma/"
  newDir <- paste0(dir, q)
  dir.create(newDir)
  setwd(newDir)
  
  ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411TumorLogrSmall.txt",
                             "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411tumorBafSmall.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
  ascat.plotRawData(ascat.bc)
  ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = ascat.gg, penalty = q)
  ascat.plotSegmentedData(ascat.bc)
  optimalRes <- NULL
  for (j in listOfGamma) {
    no_cores <- 5
    cl <- makeCluster(no_cores, type="FORK")  
    registerDoParallel(cl) 
    res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
                   .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscatGam(i, sampleMap, tcDf, j)
    stopCluster(cl)
    res$sample <- paste0(res$sample, "_pen", i, "_gam", j)
    optimalRes <- rbind(optimalRes, res)
  }
  write.table(optimalRes, "./optimalRes.txt",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

### from results above, the penalty doesn't seem to matter much from low to high
### gamma which defines what a one copy loss is does. i.e difference from 0.15 to 0.25 is where one copy losses become cn-loh
### test between 0.10 and 0.20 with default penalty for 70 - one last thing to test is the segment size .. this will change accuracy
### of certain segments lacking het snps i think



listOfPenalty <- 70
listOfGamma <- seq(0.10,  0.20, 0.01)
for (i in listOfPenalty) {
  dir <- "/br_z1/kevin_storage/ASCAT/testPenaltyGamma2/"
  newDir <- paste0(dir, i)
  dir.create(newDir)
  setwd(newDir)
  
  ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411TumorLogrSmall.txt",
                             "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411tumorBafSmall.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
  ascat.plotRawData(ascat.bc)
  ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = ascat.gg)
  ascat.plotSegmentedData(ascat.bc)
  optimalRes <- NULL
  for (j in listOfGamma) {
    no_cores <- 5
    cl <- makeCluster(no_cores, type="FORK")  
    registerDoParallel(cl) 
    res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
                   .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscatGam(i, sampleMap, tcDf, j)
    stopCluster(cl)
    res$sample <- paste0(res$sample, "_pen", i, "_gam", j)
    optimalRes <- rbind(optimalRes, res)
  }
  write.table(optimalRes, "./optimalRes.txt",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


### testing different window segments


listOfPenalty <- 70
listOfGamma <- seq(0.10,  0.20, 0.01)
for (i in listOfPenalty) {
  dir <- "/br_z1/kevin_storage/ASCAT/testPenaltyGamma2/"
  newDir <- paste0(dir, i)
  dir.create(newDir)
  setwd(newDir)
  
  ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411TumorLogrSmall.txt",
                             "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411tumorBafSmall.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
  ascat.plotRawData(ascat.bc)
  ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15Seg100")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = ascat.gg)
  ascat.plotSegmentedData(ascat.bc)
  optimalRes <- NULL
  for (j in listOfGamma) {
    no_cores <- 5
    cl <- makeCluster(no_cores, type="FORK")  
    registerDoParallel(cl) 
    res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
                   .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscatGam(i, sampleMap, tcDf, j)
    stopCluster(cl)
    res$sample <- paste0(res$sample, "_pen", i, "_gam", j)
    optimalRes <- rbind(optimalRes, res)
  }
  write.table(optimalRes, "./optimalRes.txt",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}



listOfPenalty <- 70
listOfGamma <- seq(0.10,  0.20, 0.01)
for (i in listOfPenalty) {
  dir <- "/br_z1/kevin_storage/ASCAT/testPenaltyGamma3/"
  newDir <- paste0(dir, i)
  dir.create(newDir)
  setwd(newDir)
  
  ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411TumorLogrSmall.txt",
                             "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220411tumorBafSmall.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
  ascat.plotRawData(ascat.bc)
  ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15Seg30")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = ascat.gg)
  ascat.plotSegmentedData(ascat.bc)
  optimalRes <- NULL
  for (j in listOfGamma) {
    no_cores <- 5
    cl <- makeCluster(no_cores, type="FORK")  
    registerDoParallel(cl) 
    res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
                   .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscatGam(i, sampleMap, tcDf, j)
    stopCluster(cl)
    res$sample <- paste0(res$sample, "_pen", i, "_gam", j)
    optimalRes <- rbind(optimalRes, res)
  }
  write.table(optimalRes, "./optimalRes.txt",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


### trying various gamma pen 7 for raw illumina


listOfPenalty <- 70
listOfGamma <- seq(0.17,  0.31, 0.02)
for (i in listOfPenalty) {
  dir <- "/br_z1/kevin_storage/ASCAT/testPenaltyGammaIllu15/"
  newDir <- paste0(dir, i)
  dir.create(newDir)
  setwd(newDir)
  
  ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220418IlluTumorLogrSmall.txt",
                             "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220418IlluTumorBafSmall.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
  ascat.plotRawData(ascat.bc)
  ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15Seg100")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = ascat.gg)
  ascat.plotSegmentedData(ascat.bc)
  optimalRes <- NULL
  for (j in listOfGamma) {
    no_cores <- 5
    cl <- makeCluster(no_cores, type="FORK")  
    registerDoParallel(cl) 
    res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
                   .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscatGam(i, sampleMap, tcDf, j)
    stopCluster(cl)
    res$sample <- paste0(res$sample, "_pen", i, "_gam", j)
    optimalRes <- rbind(optimalRes, res)
  }
  write.table(optimalRes, "./optimalRes.txt",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}




listOfPenalty <- 70
listOfGamma <- seq(0.11,  0.20, 0.02)
for (i in listOfPenalty) {
  dir <- "/br_z1/kevin_storage/ASCAT/testPenaltyGammaIllu05/"
  newDir <- paste0(dir, i)
  dir.create(newDir)
  setwd(newDir)
  
  ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220418IlluTumorLogrSmall.txt",
                             "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220418IlluTumorBafSmall.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
  ascat.plotRawData(ascat.bc)
  ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes05")
  ascat.bc = ascat.aspcf(ascat.bc, ascat.gg = ascat.gg)
  ascat.plotSegmentedData(ascat.bc)
  optimalRes <- NULL
  for (j in listOfGamma) {
    no_cores <- 5
    cl <- makeCluster(no_cores, type="FORK")  
    registerDoParallel(cl) 
    res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
                   .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscatGam(i, sampleMap, tcDf, j)
    stopCluster(cl)
    res$sample <- paste0(res$sample, "_pen", i, "_gam", j)
    optimalRes <- rbind(optimalRes, res)
  }
  write.table(optimalRes, "./optimalRes.txt",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}








### final sets run. one with 15 cutoff, penalty 70, gamma 0.15 regular segment size - no other filets
### second set is 3030 cutoff - 
### second to last is penalty 70, gamma 0.15, regular segment size. then filter SNPs by proportion of het i.e AF > 0.8 | AF < 0.2

bafData <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt", sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
bafHom <- bafData[,4:ncol(bafData)]
propHom <- (bafHom > 0.9 | bafHom < 0.1)
propHom2 <- apply(propHom, 1, sum)/ncol(propHom)
quantile(propHom2, seq(0, 1, 0.05), na.rm = TRUE)
meanBafHom2 <- apply(bafHom, 1, mean)
sdBafHom2 <- apply(bafHom, 1, sd)
sd30mean30 <- rownames(ascat.bc$SNPpos)[intersect(which(meanBafHom2 > 0.7 | meanBafHom2 < 0.3),
                                                  which(sdBafHom2 < 0.3))]
filtBaf3030 <- bafData
filtBafProp90 <- bafData

propHomSnps100 <- bafData$marker[which(propHom2 == 1)]
propHomSnps90 <- bafData$marker[which(propHom2 > 0.90)]
propHomSnps80 <- bafData$marker[which(propHom2 > 0.80)]

filtBaf3030[which(filtBaf3030$marker %in% sd30mean30), 4:ncol(filtBaf3030)] <- NA
filtBafProp90[which(filtBafProp90$marker %in% propHomSnps90), 4:ncol(filtBafProp90)] <- NA



write.table(filtBaf3030, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3Filted3030.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(filtBafProp90, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3FiltedProp90.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



### (1) res15 no filt

tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res15NoFilt/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### (1.5) res10 no filt

tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res10NoFilt/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### (1.7) res05 no filt

tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res05NoFilt/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405
                           testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes05")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


### (2) res05 sd 3030


tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res05Sd3030/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3Filted3030.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes05")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


### (2.5) res10 sd 3030


tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res10Sd3030/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3Filted3030.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



### (3) res15 sd 3030


tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res15Sd3030/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3Filted3030.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### 4 res 05 prop90


tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res05Prop90/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3FiltedProp90.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes05")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

if (names(ascat.bc)[11] == "failedarrays") {
  ascat.bc[["failedarrays"]] <- NULL
}

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### 5 res 10 prop90


tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res10Prop90/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3FiltedProp90.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

if (names(ascat.bc)[11] == "failedarrays") {
  ascat.bc[["failedarrays"]] <- NULL
}

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


### 6 res15 prop90

tmpDir <- "/br_z1/kevin_storage/ASCAT/20220414final/res15Prop90/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorLogR3.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220414testTumorBaf3FiltedProp90.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

if (names(ascat.bc)[11] == "failedarrays") {
  ascat.bc[["failedarrays"]] <- NULL
}

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




### res 05 seg 100 gamm 0.29

tmpDir <- "/br_z1/kevin_storage/ASCAT/illuRes05Seg100/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/ASCAT/20220410illuminaLogRr.txt",
                           "/br_z1/kevin_storage/ASCAT/20220410illuminaBaf.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes05Seg100")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

if (names(ascat.bc)[11] == "failedarrays") {
  ascat.bc[["failedarrays"]] <- NULL
}

no_cores <- 20
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir, gamma = 0.29, illu = TRUE)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### res 10 seg 100 gamm 0.29

tmpDir <- "/br_z1/kevin_storage/ASCAT/illuRes10Seg100/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/ASCAT/20220410illuminaLogRr.txt",
                           "/br_z1/kevin_storage/ASCAT/20220410illuminaBaf.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10Seg100")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

if (names(ascat.bc)[11] == "failedarrays") {
  ascat.bc[["failedarrays"]] <- NULL
}

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir, gamma = 0.29, illu = TRUE)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



### raw Illu, res 15 seg 10 0 gamma 0.29


tmpDir <- "/br_z1/kevin_storage/ASCAT/finalIlluRes15Seg100/"
setwd(tmpDir)
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/ASCAT/20220410illuminaLogRr.txt",
                           "/br_z1/kevin_storage/ASCAT/20220410illuminaBaf.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, tmpDir)
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes15Seg100")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = tmpDir, ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = tmpDir)

if (names(ascat.bc)[11] == "failedarrays") {
  ascat.bc[["failedarrays"]] <- NULL
}

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat(i, sampleMap, tcDf,tmpDir, gamma = 0.29, illu = TRUE)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "optimalRes.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


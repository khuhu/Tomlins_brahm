### alternatively, I need to see between all samples, which are always homozygous - these I can possibly rule out for segmentation

library(grid)
library(ggplot2)
library(stringr)

freqPlot_baf2 <- function(df, main = "no title", chromTextSpec = NULL){
  require(ggplot2)
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/br_z1/kevin_storage/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  ### divide positions by megabase locations
  
  df$CHROM <- as.numeric(str_remove(df$CHROM, "chr"))
  df$point <- df$POS/1e6
  
  for (i in unique(df$CHROM)) {
    df$point[which(df$CHROM == i)] <- df$point[which(df$CHROM == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  df <- df[,c("CHROM", "point", "AF", "color")]
  colnames(df) <- c("chrom", "pos", "af", "color")
  
  
  if (nrow(df) == 0) {
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.025, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  } else{
    # df$color <- "#000000"
    # df$color[which(df$af < 0.3)] <- "#FF0000"
    # df$color[which(df$af > 0.7)] <- "#FF0000"
    
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.025, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) + ggtitle(main)
  }
}


# bafData <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorBaf.txt", sep = "\t", header = TRUE,
#                       stringsAsFactors = FALSE, check.names = FALSE)

bafData <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405testTumorBaf3.txt", sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)


load("/br_z1/kevin_storage/ASCAT/20220402ascatbcRes05NoMatched.rds")

bafHom <- bafData[,4:ncol(bafData)]
# meanAf <- apply(bafHom, 1, mean)
# medianAf <- apply(bafHom, 1, median)
# 
# homMarkers <- bafData$marker[intersect(which(meanAf > 0.2 | meanAf < 0.8), which(medianAf > 0.2 | medianAf < 0.8))]

propHom <- (bafHom > 0.9 | bafHom < 0.1)
propHom2 <- apply(propHom, 1, sum)/ncol(propHom)

coefVar <- apply(bafHom, 1, sd)/apply(bafHom, 1, mean)
meanBafHom2 <- apply(bafHom, 1, mean)
sdBafHom2 <- apply(bafHom, 1, sd)
plot(meanBafHom2, sdBafHom2)


# i <- "5516-KH-48"
# normalBafs <- ascat.bc$Tumor_BAF[,i]
# normalHomSnps <- rownames(ascat.bc$SNPpos)[which(normalBafs < 0.3 | normalBafs > 0.7)]
# normalHetSnps <- rownames(ascat.bc$SNPpos)[which(normalBafs > 0.3 & normalBafs < 0.7)]

### Snps left after filtering by proportion: 100  is 60% of snps, 0.90 is 33%, 0.80  ~25%, 0.50 is 4%, 0.25 is 0.8%, 0.10 is 0.04%, 0.05 is .025%

quantile(propHom2, seq(0, 1, 0.05), na.rm = TRUE)

propHomSnps100 <- bafData$marker[which(propHom2 == 1)]
propHomSnps90 <- bafData$marker[which(propHom2 > 0.90)]
propHomSnps80 <- bafData$marker[which(propHom2 > 0.80)]

propHomSnps50 <- bafData$marker[which(propHom2 > 0.50)]
propHomSnps25 <- bafData$marker[which(propHom2 > 0.25)]
propHomSnps10  <- bafData$marker[which(propHom2 > 0.10)]
propHomSnps05  <- bafData$marker[which(propHom2 > 0.05)]

meanSdHet <- rownames(ascat.bc$SNPpos)[intersect(which(meanBafHom2 > 0.4 & meanBafHom2 < 0.6),
                                                 which(sdBafHom2 < 0.4))]


quantile(propHom2, seq(0, 1, 0.05), na.rm = TRUE)



for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% homMarkers)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% homMarkers)] <- NA
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% homMarkers)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/allHomSnpFilt/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  dev.off()
}




for (i in colnames(bafData)[4:ncol(bafData)]) {
  tmpBafDf <- bafData[,2:3]
  rownames(tmpBafDf) <- bafData$marker
  tmpBafDf$AF <- bafData[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  # tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% propHomSnps100)] <- NA
  tmpBafDf2 <- tmpBafDf2[which(rownames(tmpBafDf2) %in% propHomSnps100),]
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  # tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% propHomSnps100)] <- NA
  tmpBafDf3 <- tmpBafDf3[-which(rownames(tmpBafDf3) %in% propHomSnps100),]
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpsProp100/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(bafData)[4:ncol(bafData)]) {
  tmpBafDf <- bafData[,2:3]
  rownames(tmpBafDf) <- bafData$marker
  tmpBafDf$AF <- bafData[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  # tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% propHomSnps100)] <- NA
  tmpBafDf2 <- tmpBafDf2[which(rownames(tmpBafDf2) %in% propHomSnps90),]
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  # tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% propHomSnps100)] <- NA
  tmpBafDf3 <- tmpBafDf3[-which(rownames(tmpBafDf3) %in% propHomSnps90),]
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpsProp90/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  # grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  grid::grid.draw(rbind(tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



for (i in colnames(bafData)[4:ncol(bafData)]) {
  tmpBafDf <- bafData[,2:3]
  rownames(tmpBafDf) <- bafData$marker
  tmpBafDf$AF <- bafData[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  # tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% propHomSnps100)] <- NA
  tmpBafDf2 <- tmpBafDf2[which(rownames(tmpBafDf2) %in% propHomSnps80),]
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  # tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% propHomSnps100)] <- NA
  tmpBafDf3 <- tmpBafDf3[-which(rownames(tmpBafDf3) %in% propHomSnps80),]
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpsProp80/", i, ".png"),
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
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% propHomSnps50)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% propHomSnps50)] <- NA
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% homMarkers)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpsProp50/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  dev.off()
}


for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% propHomSnps25)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% propHomSnps25)] <- NA
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% homMarkers)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpsProp25/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  dev.off()
}

for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[-which(rownames(tmpBafDf2) %in% propHomSnps10)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[which(rownames(tmpBafDf3) %in% propHomSnps10)] <- NA
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% homMarkers)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpsProp10/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  dev.off()
}


for (i in colnames(ascat.bc$Tumor_BAF)) {
  tmpBafDf <- ascat.bc$SNPpos
  tmpBafDf$AF <- ascat.bc$Tumor_BAF[,i]
  colnames(tmpBafDf) <- c("CHROM", "POS", "AF")
  tmpBafDf$color <- "#000000"
  tmpBafDf2 <- tmpBafDf
  tmpBafDf2$color <- "#0000FF"
  tmpBafDf2$AF[which(rownames(tmpBafDf2) %in% meanSdHet)] <- NA
  tmpBafDf3 <- tmpBafDf
  tmpBafDf3$color <- "#FF0000"
  tmpBafDf3$AF[-which(rownames(tmpBafDf3) %in% meanSdHet)] <- NA
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% homMarkers)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  # tmpBafDf$color[which(rownames(tmpBafDf) %in% normalHomSnps)] <- "#0000FF"
  tmpBafGraph <- ggplotGrob(freqPlot_baf2(tmpBafDf, main = paste0("all_", i)))
  tmpBafGraph2 <- ggplotGrob(freqPlot_baf2(tmpBafDf2, main = paste0("hom_", i)))
  tmpBafGraph3 <- ggplotGrob(freqPlot_baf2(tmpBafDf3, main = paste0("het_", i)))
  png(filename = paste0("/br_z1/kevin_storage/ASCAT/snpHetMean4060/", i, ".png"),
      width = 2000, height = 1000, res = 200)
  grid::grid.newpage()
  grid::grid.draw(rbind(tmpBafGraph, tmpBafGraph2, tmpBafGraph3))
  dev.off()
}



### using the lists above to filter out homozygous SNPs, then run as bulk with no matched normal


library(argyle)
load("/br_z1/kevin_storage/misc/snps.gigamuga.Rdata")
load("/br_z1/kevin_storage/misc/clusters.gigamuga.Rdata")
snpRes <- read.beadstudio(prefix = "",
                          in.path = "/br_z1/kevin_storage/advancedGenomicsCore/tmp/",
                          snps = snps)
normedSnpRes <- tQN(snpRes, clusters = clusters)
baf_normed <- get.baf(normedSnpRes)
baf_mat <- dcast(baf_normed, marker + chr + pos ~ iid, value.var = c("BAF"))

test_intense <- get.intensity(normedSnpRes)
logR_mat <- dcast(test_intense, marker + chr + pos ~ iid, value.var = c("si"))
baf_mat$chr <- stringr::str_remove(baf_mat$chr, "chr")
logR_mat$chr <- stringr::str_remove(logR_mat$chr, "chr")

test_ascat_tumor.logR <- logR_mat
test_ascat_tumor.logR[ , 4:51] <- log2(test_ascat_tumor.logR[ , 4:51])
test_ascat_tumor.baf <- baf_mat
test_ascat_tumor.baf_allFilt <- test_ascat_tumor.baf 
test_ascat_tumor.baf_normFilt <- test_ascat_tumor.baf 

test_ascat_tumor.baf_allFilt[which(test_ascat_tumor.baf_allFilt$marker %in% homMarkers), 4:ncol(test_ascat_tumor.baf_allFilt)] <- NA_integer_
test_ascat_tumor.baf_normFilt[which(test_ascat_tumor.baf_normFilt$marker %in% normalHomSnps), 4:ncol(test_ascat_tumor.baf_allFilt)] <- NA_integer_

write.table(test_ascat_tumor.logR, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405allLogR.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(test_ascat_tumor.baf_allFilt, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405allFiltBaf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(test_ascat_tumor.baf_normFilt, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405normFiltBaf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


source("/home/kevhu/scripts/20220331ascat.predictGermlineGenotypes.R")
# norm filt
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405allLogR.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405normFiltBaf.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/normFiltRho/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/normFiltRho/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/normFiltRho/", gamma = 0.25)

save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220405ascatbcNormFilt.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220405ascatoutputNormFilt.rds")


# all filt
ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405allLogR.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220405allFiltBaf.txt")
ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes10")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/allFiltRho/", ascat.gg = ascat.gg)
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/allFiltRho/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/allFiltRho/")

save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220405ascatbcAllFilt.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220405ascatoutputAllFilt.rds")


nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}


ccpAmps <- read.table("/br_z1/kevin_storage/misc/20211115allCcpAmps.txt", sep = "\t",
                      header = TRUE,  stringsAsFactors = FALSE)
colnames(ccpAmps)[7:ncol(ccpAmps)] <- nameStripper(colnames(ccpAmps)[7:ncol(ccpAmps)])
ccpAmps[,7:ncol(ccpAmps)] <- log2(ccpAmps[,7:ncol(ccpAmps)])
tmpMat <- ccpAmps[,7:ncol(ccpAmps)] 
tmpMat[tmpMat > 3] <- 3
tmpMat[tmpMat < -3] <- -3
ccpAmps[,7:ncol(ccpAmps)] <- tmpMat

ccpTable <- read.table("/br_z1/kevin_storage/misc/20211115ccpSegTableFilt_1SD.txt", sep = "\t",
                       header = TRUE,  stringsAsFactors = FALSE)

ccpTable$ID <- nameStripper(ccpTable$ID)

sumfiles <- read.table("/br_z1/kevin_storage/misc/2021115ccpSummaryFile.txt", header = TRUE,
                       sep = "\t", stringsAsFactors = FALSE)

sumfiles$Sample.Name <- nameStripper(sumfiles$Sample.Name)
sumfiles$Uniformity <- as.numeric(str_remove(sumfiles$Uniformity,"%"))
sumfiles$On.Target <- as.numeric(str_remove(sumfiles$On.Target,"%"))

sumfiles <- sumfiles[which(sumfiles$Uniformity > 80 & sumfiles$On.Target > 80 & sumfiles$Mean.Depth > 300),]
sumfiles <- sumfiles[-which(sumfiles$Sample.Name == "nne"),]

ccpTable2 <- ccpTable[which(ccpTable$ID %in% unique(sumfiles$Sample.Name)),]


# dirs of interest Auto_user_AUS5-45-PR959toPR970_CCP_218_116
dirsOfInterest <- c("Auto_user_AUS5-45-PR959toPR970_CCP_218_116",
                    "Auto_user_AUS5-18-MC1104X49_CCP_191_060",
                    "Auto_user_AUS5-50-BDC55_CCP_223_126",
                    "Auto_user_AUS5-75-BDC71_CCP_253_183")

combinedAnnoTable <- NULL
annoFiles <- system('find /br_z1/kevin_storage/mouseData/humanCcpAnno/ -type f -name "*_anno.txt"', intern = TRUE)
for (i in annoFiles) {
  tmpAnno <- read.table(i, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  combinedAnnoTable <- rbind(combinedAnnoTable, tmpAnno)
}

combinedAnnoTable <- combinedAnnoTable[which(combinedAnnoTable$mm10_mpgpv6_Indels == "."),]


snpFilters <- c("A", "T", "G", "C")
combinedAnnoTable <- combinedAnnoTable[which(combinedAnnoTable$Ref %in% snpFilters),]
combinedAnnoTable <- combinedAnnoTable[which(combinedAnnoTable$Alt %in% snpFilters),]
combinedAnnoTable <- combinedAnnoTable[-which(combinedAnnoTable$AF < 0.1 | combinedAnnoTable$AF > 0.9),]


combinedAnnoTable$POS <- combinedAnnoTable$Start
combinedAnnoTable$AF <- as.numeric(combinedAnnoTable$AF)

fdpFilt <- which(combinedAnnoTable$FDP > 150)
faoFilt <- which(combinedAnnoTable$FAO > 5)
freqFilt <- which(combinedAnnoTable$AF > 0.05)
hrunFilt <- which(combinedAnnoTable$HRUN < 4)
qualFilt <- which(combinedAnnoTable$QUAL > 30)
gqFilt <- which(combinedAnnoTable$GQ > 20)
strandRatio <- intersect(which(combinedAnnoTable$FSAF/combinedAnnoTable$FSAR > 0.2),
                         which(combinedAnnoTable$FSAF/combinedAnnoTable$FSAR < 5))
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, hrunFilt, gqFilt, qualFilt))
combinedAnnoTable <- combinedAnnoTable[goodSamps,]
combinedAnnoTable$Sample <- nameStripper(combinedAnnoTable$Sample)

colnames(combinedAnnoTable)[2] <- "CHROM"
combinedAnnoTable$CHROM <- str_replace(combinedAnnoTable$CHROM, "chrX", "chr23")
### BAFs prior to additional filtering
sampleNames <- unique(ccpTable2$ID)

# samples don't have all chromosome segs


# for (i in sampleNames) {
#   testDf_cn <- ccpTable2[which(ccpTable2$ID == i),]
#   testDf_amp <- cbind(ccpAmps[,c(3:6)], ccpAmps[[i]])
#   b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = i, chromTextSpec = "hg19")
#   
#   testDf <- combinedAnnoTable[which(combinedAnnoTable$Sample == i),]
#   a <- freqPlot_baf(testDf, chromTextSpec = "hg19")
#   
#   gA <- ggplotGrob(a)
#   gB_2 <- ggplotGrob(b_2)
#   
#   dev.off()
#   pdf(file = paste0("/br_z1/kevin_storage/misc/20211115ccpBafs/20211115_ccp_noAddFilt", i, ".pdf"), width = 15, height = 7)
#   grid::grid.newpage()
#   grid::grid.draw(rbind(gB_2, gA))
#   dev.off()
#   
# }


### BAFs after filtering

combinedAnnoTable$string <- paste0(combinedAnnoTable$CHROM, ":", combinedAnnoTable$POS)

tableAf <- NULL
for (i in unique(combinedAnnoTable$string)) {
  tmpTable <- combinedAnnoTable[which(combinedAnnoTable$string ==  i),]
  count <- nrow(tmpTable)
  medianAf <- signif(median(tmpTable$AF) * 100, digits = 3)
  meanAf <- signif(mean(tmpTable$AF) * 100, digits = 3)
  varAf <- signif(var(tmpTable$AF) * 100, digits = 3)
  madAf <- signif(mad(tmpTable$AF) * 100, digits = 3)
  tableAf <- rbind(tableAf, c(i, count, medianAf, meanAf, varAf, madAf))
  
}

tableAf <- as.data.frame(tableAf, stringsAsFactors = FALSE)
tableAf[,2:6] <- lapply(tableAf[,2:6], as.numeric)
colnames(tableAf) <- c("position", "count", "medianAf", "meanAf", "varAf", "madAf")
tableAf$varAf[which(is.na(tableAf$varAf))] <- 0
filt <- which(tableAf$medianAf > 45 & tableAf$medianAf < 55)
table_filt <- tableAf[filt,]

goodPosString_hg19  <- tableAf$position[filt]

goodPosString_hg19  <- tableAf$position[which(tableAf$count >= 1)]
combinedAnnoTable_medFilt <- combinedAnnoTable[which(combinedAnnoTable$string %in% goodPosString_hg19),]



variants_grange <- unique(GRanges(seqnames = combinedAnnoTable_medFilt$CHROM,
                                  IRanges(start = combinedAnnoTable_medFilt$POS, end = combinedAnnoTable_medFilt$POS)))
segs_grange <- GRanges(seqnames = paste0("chr", ccpTable2$chrom),
                       IRanges(start = ccpTable2$loc.start,  end = ccpTable2$loc.end))

combinedAnnoTable_medFilt$absAf <- abs(combinedAnnoTable_medFilt$AF -  0.5)
ccpTable2$dipCount <- 2^ccpTable2$seg.mean * 2
ccpTable2$absDipCount <- abs(ccpTable2$dipCount - 2)
### how to deal with losses and gains in terms of correlation? convert the segments
### into absolute counts i.e reference of diploid

corTable <- NULL
for (i in 1:length(variants_grange)) {
  tmpGrange <- variants_grange[i]
  tmpSegs <- ccpTable2[queryHits(findOverlaps(segs_grange, tmpGrange)),]
  
  
  tmpPos <- paste0(variants_grange@seqnames[i], ":", variants_grange@ranges@start[i])
  tmpVcf <- combinedAnnoTable_medFilt[which(combinedAnnoTable_medFilt$string == tmpPos),]
  
  allNames <- intersect(tmpVcf$Sample, tmpSegs$ID)
  
  tmpSegs <- tmpSegs[which(tmpSegs$ID %in% allNames),]
  tmpVcf <- tmpVcf[which(tmpVcf$Sample %in% allNames),]
  
  tmpSegs <- tmpSegs[match(tmpVcf$Sample, tmpSegs$ID),]
  
  tmpCor <- cor(tmpSegs$absDipCount, tmpVcf$absAf)
  meanSeg <- median(tmpSegs$absDipCount)
  sdSeg <- mad(tmpSegs$absDipCount)
  meanAbsVf <- median(tmpVcf$absAf) * 100
  sdAbsVf <- mad(tmpVcf$absAf) * 100
  tmpCount <- nrow(tmpVcf)
  medVf <- median(tmpVcf$AF) * 100
  madVf <- mad(tmpVcf$AF) * 100
  
  corTable <- rbind(corTable, c(tmpPos, tmpCount, tmpCor, meanSeg, sdSeg, meanAbsVf, sdAbsVf, medVf, madVf))
}

corTable <- data.frame(corTable, stringsAsFactors = FALSE)
colnames(corTable) <- c("Position", "Count.hgsc", "Cor.p", "Med.Seg",
                        "Mad.Seg", "Med.AbsVf", "Mad.AbsVf", "Med.Vf", "Mad.Vf")
corTable[,2:ncol(corTable)] <- lapply(corTable[,2:ncol(corTable)], function(x) signif(as.numeric(x), digits = 3))

corTable.na <- corTable[which(is.na(corTable$Cor.p)),]
corTable <- corTable[-which(is.na(corTable$Cor.p)),]

badPos_CorNa <- corTable.na$Position[which(corTable.na$Med.Vf > 55 | corTable.na$Med.Vf < 45)]


# badPos_Cor <- corTable$Position[which(corTable$Cor.p < 0.30)]
badPos_Cor <- corTable$Position[which(corTable$Cor.p < 0.20)]
addBadPos <- c(badPos_CorNa, badPos_Cor)
combinedAnnoTable_allFilt <- combinedAnnoTable_medFilt
combinedAnnoTable_allFilt <- combinedAnnoTable_allFilt[-which(combinedAnnoTable_allFilt$string %in% addBadPos),]

# mc1110, pr94020ccp, bdc52
i <- "mc1110"
i <- "pr94020ccp"
i <- "bdc52"

for (i in sampleNames) {
  
  testDf_cn <- ccpTable2[which(ccpTable2$ID == i),]
  testDf_amp <- cbind(ccpAmps[,c(3:6)], ccpAmps[[i]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = paste0(i, "CCP, Baf: no filt, Baf w/ filts Cor>= 0.30"),
                        chromTextSpec = "hg19")
  
  testDf <- combinedAnnoTable[which(combinedAnnoTable$Sample == i),]
  a_nofilt <- freqPlot_baf(testDf, chromTextSpec = "hg19")
  
  testDf <- combinedAnnoTable_allFilt[which(combinedAnnoTable_allFilt$Sample == i),]
  a_0 <- freqPlot_baf(testDf, chromTextSpec = "hg19")
  
  gnoFilt <- ggplotGrob(a_nofilt)
  g0 <- ggplotGrob(a_0)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211115ccpBafs/20211115_ccp_filtComp_cor20", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, gnoFilt, g0))
  dev.off()
  
}


# mc1110, pr94020ccp, bdc52

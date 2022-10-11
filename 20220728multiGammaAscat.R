### res 05 seg 100 gamm 0.29 - calculating tables for rmse
library(foreach)
library(doParallel)
library(ASCAT)
source("/home/kevhu/scripts/20220331ascat.predictGermlineGenotypes.R")


runGridAscat_gamma <- function(i, sampleMap, tcDf, outdir, gamma  = 0.15, illu = NULL){
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
  
  ploidyVec <- seq(1.6, 5.4, 0.2)
  ploidyVec2 <- rep(ploidyVec, 4)
  purityVec2 <- NULL
  for (j in purityVec) {
    purityVec2 <- c(purityVec2, rep(j, 20))
  }
  
  
  tmp.bc <- ascat.bc
  tmp.bc[[1]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[1]][,i]), 80))
  tmp.bc[[2]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[2]][,i]), 80))
  tmp.bc[[7]] <- paste0(names(ascat.bc[[1]][i]), "_rho", purityVec2, "_psi", ploidyVec2)
  tmp.bc[[8]] <- rep("XX", 80)
  tmp.bc[[11]] <- rep(unlist(ascat.bc[[11]][i]), 80)
  tmp.bc[[12]] <- rep(unlist(ascat.bc[[12]][i]), 80)
  tmp.bc[[13]] <- rep(unlist(ascat.bc[[13]][i]), 80)
  tmp.bc[[14]] <- rep(unlist(ascat.bc[[14]][i]), 80)
  tmp.bc[[15]] <- Rfast::rep_col(unlist(ascat.bc[[15]][,i]), 80)
  for (j in 1:80) {
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
  
  ascatRes <- NULL
  varGamma <- gamma[1]
  for (varGamma in gamma) {
    if (!file.exists(paste0(outdir, "/", names(ascat.bc[[1]][i]), "_", varGamma))) {
      dir.create(paste0(outdir, "/", names(ascat.bc[[1]][i]), "_", varGamma))
    }
    setwd(paste0(outdir, "/", names(ascat.bc[[1]][i]), "_", varGamma))
    ascat.output = ascat.runAscat(tmp.bc, rho_manual = purityVec2, psi_manual = ploidyVec2, gamma = varGamma)
    
    ### I think goodness of fit is still good information to have based, but need to grab segments
    ascatDf <- data.frame("sample" = names(ascat.output$ploidy), "ploidy" = ascat.output$ploidy,
                          "purity" = ascat.output$purity, "goodnessOfFit" = ascat.output$goodnessOfFit)
    ascatDf <- ascatDf[order(ascatDf$goodnessOfFit, decreasing = TRUE), ]
    # ascatDf <- ascatDf[1:15,]
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
    
    write.table(ascat.output$segments, file = paste0(outdir, "/", "allSegResults_", varGamma, ".txt"), sep = "\t", quote = FALSE,
                col.names = TRUE, row.names = FALSE, append = TRUE)
    
    ascatDf$gamma <- varGamma
    ascatRes <- rbind(ascatRes, ascatDf)
  }
  
  
  ascatRes
}

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

tmpDir <- "/br_z1/kevin_storage/ASCAT/testingGammaRmseRes05Illu"
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

save(ascat.bc, file = "/br_z1/kevin_storage/misc/20221004ascastSeg.Robj")

### instead of above make sure I choose just 5-6 of high tumor content diploid samples
### so the foreach function just runs on those

tmpDir <- "/br_z1/kevin_storage/ASCAT/testingGammaRmseRes05IlluV2"

no_cores <- 20
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind",
               .packages = c("ASCAT", "ggplot2", "grid", "stringr")) %dopar% runGridAscat_gamma(i, sampleMap, tcDf,tmpDir, gamma = seq(0.05,1, 0.05), illu = TRUE)
stopCluster(cl)
write.table(res, file = paste0(tmpDir, "/gammaOptimalResAll.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

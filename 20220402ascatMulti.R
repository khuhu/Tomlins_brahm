library(ASCAT)
library(stringr)
source("/home/kevhu/scripts/20220331ascat.predictGermlineGenotypes.R")

# ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/testTumorLogR.txt",
#                            "/br_z1/kevin_storage/advancedGenomicsCore/tmp/testTumorBaf.txt")
# ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
# ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/res24/")
# ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrent")
# ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/res24/", ascat.gg = ascat.gg)
# ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/res24/")
# ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/res24/")
# save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220402ascatouputRes24NoMatched.rds")
# save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220402ascatbcRes24NoMatched.rds")


# ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/testTumorLogR.txt",
#                            "/br_z1/kevin_storage/advancedGenomicsCore/tmp/testTumorBaf.txt")
# ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
# ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/20220402allNoMatched05/")
# ascat.gg <- ascat.predictGermlineGenotypes(ascat.bc, "IonTorrentRes05")
# ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/20220402allNoMatched05/", ascat.gg = ascat.gg)
# save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220402ascatbcRes05NoMatched.rds")

# (1) save ascat.output. (2) make function to run the outputs at 60 ploidy anbd purity combinations (3) make custom graph so I can scale to genome size and not markers
# (4) similar to 3, but a bit easier, I need to graph BAfs to genome scale too

### from the results, I can reformat ascat.output into something useable for specific ploidies

freqPlot_baf <- function(df, main = "no title", chromTextSpec = NULL){
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
  
  df <- df[,c("CHROM", "point", "AF")]
  colnames(df) <- c("chrom", "pos", "af")
  
  
  if (nrow(df) == 0) {
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  } else{
    df$color <- "#000000"
    df$color[which(df$af < 0.3)] <- "#FF0000"
    df$color[which(df$af > 0.7)] <- "#FF0000"
    
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  }
}

ascat_cn <- function(df, chromTextSpec = NULL){
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
  for (i in unique(df$sample)) {
    main <- i
    tmpDf <- df[which(df$sample == i),]
    tmpDf_major <- tmpDf[,1:5]
    tmpDf_major$col <- "red"
    tmpDf_minor <- tmpDf[,c(1:4, 6)]
    tmpDf_minor$col <- "blue"
    colnames(tmpDf_major)[5] <- "cn"
    colnames(tmpDf_minor)[5] <- "cn"
    tmpDf_major$cn <- ifelse(tmpDf_major$cn == 0, tmpDf_major$cn, tmpDf_major$cn + 0.1)
    tmpDf_minor$cn <- ifelse(tmpDf_minor$cn == 0, tmpDf_minor$cn + 0.1, tmpDf_minor$cn - 0.4)
    tmpDf2 <- rbind(tmpDf_major, tmpDf_minor)
    tmpDf2 <- tmpDf2[which(tmpDf2$chr %in% c(1:19)),]
    tmpDf2 <- tmpDf2[order(as.numeric(tmpDf2$chr)),]
    colnames(tmpDf2) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean", "col")
    
    tmpDf2$loc.start <- tmpDf2$loc.start/1e6
    tmpDf2$loc.end <- tmpDf2$loc.end/1e6
    
    for (j in unique(tmpDf2$chrom)) {
      tmpDf2$loc.start[which(tmpDf2$chrom == j)] <- tmpDf2$loc.start[which(tmpDf2$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
      tmpDf2$loc.end[which(tmpDf2$chrom == j)] <- tmpDf2$loc.end[which(tmpDf2$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
    }
    
    tmpDf2 <- tmpDf2[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
    colnames(tmpDf2) <- c("chrom", "start", "end","cn", "col")
    
    cnplot <- ggplot(tmpDf2) + geom_hline(yintercept = c(0, 1, 2, 3,4, 5), color = "#000000") + 
      geom_rect(aes(xmin = start, xmax = end, ymin = cn, ymax = cn + 0.3), fill = tmpDf2$col) + 
      geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 5, by = 1),
                                                             limits = c(0, 5.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 4.3, label = chrom), size = 2.5) + ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5))
    
    cnplot_grob <- ggplotGrob(cnplot)
    
    dev.off()
    # pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211122_filtComps_20_", i, ".pdf"), width = 15, height = 7)
    # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
    grid::grid.newpage()
    grid::grid.draw(rbind(cnplot_grob, tmpBafGrob))
    dev.off()
  }
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
load("/br_z1/kevin_storage/ASCAT/20220402ascatbcRes24NoMatched.rds")
sampleMap <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/Sample_Map.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)
sampleMap$stripped <- nameStripper(sampleMap$ID)
sampleMap$stripped[21] <- "13085lt"



### to save the outputs to sort the n best results and vectors, I can use doParallel and save the ascat.output
### could also sort before so less filtering downstream and graph

for (i in 1:(dim(ascat.bc[[1]])[2])) {
  setwd("/br_z1/kevin_storage/ASCAT/multiRhoPsi/")
  # pain in the butt because I need to replace the object of list of 16
  # lists I can keep 3 (SNPpos), 4 (ch), 5 (chr), 6 (chrs), 9 (sexchromosomes), 10 (X_nonPAR)
  # other columns need to be pulled and replicated to the number of test ploidies
  
  matchedName <- sampleMap$stripped[which(sampleMap$Name == colnames(ascat.bc[[1]])[i])]
  if (length(tcDf$tc[which(tcDf$sample == matchedName)]) == 0) {
    tmpTc <- 0.96
  } else{
    tmpTc <- signif(tcDf$tc[which(tcDf$sample == matchedName)][1], 2)
  }
  
  if (tmpTc > 0.94) {
    purityVec <- seq(0.88, 1, 0.02)
  } else{
    purityVec <- seq(tmpTc - 0.02 * 3, tmpTc + 0.02 * 3, 0.02)
  }
  
  ploidyVec <- seq(1.6, 5.4, 0.4)
  ploidyVec2 <- rep(ploidyVec, 7)
  purityVec2 <- NULL
  for (j in purityVec) {
    purityVec2 <- c(purityVec2, rep(j, 10))
  }
  
  
  tmp.bc <- ascat.bc
  tmp.bc[[1]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[1]][,i]), 70))
  tmp.bc[[2]] <- data.frame(Rfast::rep_col(unlist(ascat.bc[[2]][,i]), 70))
  tmp.bc[[7]] <- paste0(names(ascat.bc[[1]][i]), "_rho", purityVec2, "_psi", ploidyVec2)
  tmp.bc[[8]] <- rep("XX", 70)
  tmp.bc[[11]] <- rep(unlist(ascat.bc[[11]][i]), 70)
  tmp.bc[[12]] <- rep(unlist(ascat.bc[[12]][i]), 70)
  tmp.bc[[13]] <- rep(unlist(ascat.bc[[13]][i]), 70)
  tmp.bc[[14]] <- rep(unlist(ascat.bc[[14]][i]), 70)
  tmp.bc[[15]] <- Rfast::rep_col(unlist(ascat.bc[[15]][,i]), 70)
  for (j in 1:70) {
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
  ascat.output = ascat.runAscat(tmp.bc, rho_manual = purityVec2, psi_manual = ploidyVec2)
}


### take the ascat output for the baf graph - make that and save it into an object. then 
### comment out tmp match in order to all snps



library(foreach)
library(doParallel)


runGridAscat <- function(i, sampleMap, tcDf, outdir){
  setwd(outdir)
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
  ascat.output = ascat.runAscat(tmp.bc, rho_manual = purityVec2, psi_manual = ploidyVec2)
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
  ascat_cn(df = ascat.output$segments)
  
  ascatDf
}

load("/br_z1/kevin_storage/ASCAT/20220402ascatbcRes24NoMatched.rds")

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/multiRhoPsiV2/")
stopCluster(cl)

# write.table(res, "/br_z1/kevin_storage/ASCAT/20220404res24ASCATgridTable.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


load("/br_z1/kevin_storage/ASCAT/20220402ascatbcRes05NoMatched.rds")
no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/multiRhoPsiV3/")
stopCluster(cl)

res24Table <- read.table("/br_z1/kevin_storage/ASCAT/20220404res24ASCATgridTable.txt", sep = "\t", stringsAsFactors = FALSE , header = TRUE)


# normFilt
load("/br_z1/kevin_storage/ASCAT/20220405ascatbcNormFilt.rds")
load("/br_z1/kevin_storage/ASCAT/20220405ascatoutputNormFilt.rds")
ascat.bc  <- ascat.bc[-which(names(ascat.bc)  == "failedarrays")]

no_cores <- 12
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/normFiltRho/")
stopCluster(cl)

normFiltRes <- res

# allFilt
load("/br_z1/kevin_storage/ASCAT/20220405ascatbcAllFilt.rds")
load("/br_z1/kevin_storage/ASCAT/20220405ascatoutputAllFilt.rds")
no_cores <- 18
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/allFiltRho/")
stopCluster(cl)

allFiltRes <- res


library(vcfR)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)


allVcfs <- system('find /br_z1/kevin_storage/mouseData/tvcOut/ -type f -name "*.vcf" | grep -v "TSVC" | grep -v "small_variants" | grep -v "indel" | grep -v "geno"', intern = TRUE)

genoVcfs <- system('find /br_z1/kevin_storage/mouseData/tvcOut/ -type f -name "*_geno.vcf"', intern = TRUE)

#small subset
#allVcfs <- system('find /br_z1/kevin_storage/mouseData/tvcOut/Auto_user_AUS5-142-MG_cho_20210701_357_353/ -type f -name "*.vcf" | grep -v "TSVC" | grep -v "small_variants" | grep -v "indel" | grep -v "geno"', intern = TRUE)
#genoVcfs <- system('find /br_z1/kevin_storage/mouseData/tvcOut/Auto_user_AUS5-142-MG_cho_20210701_357_353/ -type f -name "*_geno.vcf"', intern = TRUE)

allAmps <- read.table("/br_z1/kevin_storage/misc/20211129all_mouse_amps.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)

allSegs <- read.table("/br_z1/kevin_storage/misc/20211129allMMouseSegResNoTc.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)
pcfSegs <-  read.table("/br_z1/kevin_storage/misc/20211025allPcfRnaseq_zscore.txt", sep = "\t", stringsAsFactors = FALSE,
                       header = TRUE)

allRnaCbs <- read.table("/br_z1/kevin_storage/misc/20211109allCbsRnaseq_zscore.txt", sep = "\t", stringsAsFactors = FALSE,
                        header = TRUE)

amps_normal <- read.table("/br_z1/kevin_storage/misc/20211101normals_amps.txt", sep = "\t", stringsAsFactors = FALSE,
                          header = TRUE, check.names = FALSE)

segs_normal <- read.table("/br_z1/kevin_storage/misc/20211101allsegRes_normal.txt", sep = "\t", stringsAsFactors = FALSE,
                          header = TRUE)


genoBed <- read.table("/br_z1/kevin_storage/misc/20211120snpsOnIAD202670_167_positions.bed", sep = "\t",
                      stringsAsFactors = FALSE, header = FALSE)

positionString <- paste0(genoBed$V1,":" ,genoBed$V3)

allSegs$seg.mean[which(abs(allSegs$seg.mean) < 0.2)] <- 0
#allSegs$seg.mean[which(allSegs$seg.mean < (log2(5/4)) & allSegs$seg.mean > (log2(3/4)))] <- 0
allSegs$seg.mean[which(allSegs$q.val2 > 0.05)] <- 0
allSegs$seg.mean[which(allSegs$seg.mean < -3)] <- -3
allSegs$seg.mean[which(allSegs$seg.mean > 3)] <- 3
#allSegs$seg.mean[which((allSegs$end.pos - allSegs$start.pos) < (10 * 1e6))] <- 0


pcfSegs$seg.mean[which(abs(pcfSegs$seg.mean) < 0.2)] <- 0
#pcfSegs$seg.mean[which(pcfSegs$seg.mean < (log2(5/4)) & pcfSegs$seg.mean > (log2(3/4)))] <- 0
pcfSegs$seg.mean[which(pcfSegs$seg.mean < -3)] <- -3
pcfSegs$seg.mean[which(pcfSegs$seg.mean > 3)] <- 3
pcfSegs$seg.mean[which((pcfSegs$end.pos - pcfSegs$start.pos) < (10 * 1e6))] <- 0
pcfSegs$seg.mean[which(pcfSegs$q.val2 > 0.05)] <- 0


allRnaCbs$seg.mean[which(abs(allRnaCbs$seg.mean) < 0.2)] <- 0
#allRnaCbs$seg.mean[which(allRnaCbs$seg.mean < (log2(5/4)) & allRnaCbs$seg.mean > (log2(3/4)))] <- 0
allRnaCbs$seg.mean[which(allRnaCbs$q.val2 > 0.05)] <- 0
allRnaCbs$seg.mean[which(allRnaCbs$seg.mean < -3)] <- -3
allRnaCbs$seg.mean[which(allRnaCbs$seg.mean > 3)] <- 3
allRnaCbs$seg.mean[which((allRnaCbs$loc.end - allRnaCbs$loc.start) < (10 * 1e6))] <- 0
allRnaCbs$seg.mean[which(allRnaCbs$q.val2 > 0.05)] <- 0

allNames <- colnames(allAmps)[7:ncol(allAmps)]
for (i in allNames) {
  tmpVector <- allAmps[[i]]
  tmpVector[tmpVector > 3] <- 3
  tmpVector[tmpVector < -3] <- -3
  allAmps[[i]] <- tmpVector
}


# all_vcf_processed <- NULL
# for (i in allVcfs[10]) {
#   sampleName <- str_remove(unlist(str_split(i, "/"))[7], ".vcf")
#   tmpVcf <- vcfR::read.vcfR(i)
#   tmpVcf2 <- data.frame(tmpVcf@fix, stringsAsFactors = FALSE)
#   pass_filt <- which(tmpVcf2$FILTER == "PASS")
# 
#   GQ <- unlist(lapply(strsplit(tmpVcf@gt[,2], ":"), '[[', 2))[pass_filt]
#   tmpVcf2 <- tmpVcf2[pass_filt,]
# 
#   tmpVcf3 <- NULL
#   strSplitRes <- strsplit(tmpVcf2$INFO, ";")
#   for (j in 1:length(strSplitRes)) {
#     tmpVector <- unlist(strSplitRes[[j]])[grep("^AF|AO|DP|FR|FAO|FDP|FSAF|FSAR|HRUN",unlist(strSplitRes[[j]]))]
#     tmpVector <- tmpVector[-grep("FRO", tmpVector)]
#     tmpVector <- str_remove(tmpVector, ".*=")
#     tmpVector2 <- unlist(c(tmpVcf2[j,1:7], tmpVector))
#     tmpVcf3 <- rbind(tmpVcf3, tmpVector2)
#   }
#   tmpVcf3 <- data.frame(tmpVcf3, stringsAsFactors = FALSE)
#   tmpVcf3$GQ <- unlist(GQ)
#   rownames(tmpVcf3) <- NULL
#   colnames(tmpVcf3) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
#                          "FILT", "AF", "AO", "DP", "FAO", "FDP", "FR",
#                          "FSAF", "FSAR", "HRUN", "GQ")
#   tmpVcf3[,c(6, 8:12, 14:17)] <- lapply(tmpVcf3[,c(6, 8:12, 14:17)], as.numeric)
#   tmpVcf3 <- cbind("sample" = sampleName, tmpVcf3)
#   all_vcf_processed <- rbind(all_vcf_processed, tmpVcf3)
# }


## run for mgp
# cl <- makeCluster(30)
# registerDoParallel(cl)
# 
# all_vcf_processed <- NULL
# all_vcf_processed <- foreach(i = allVcfs, .combine = 'rbind', .packages = c('vcfR', 'stringr')) %dopar% {
#   sampleName <- str_remove(unlist(str_split(i, "/"))[7], ".vcf")
#   tmpVcf <- vcfR::read.vcfR(i)
#   tmpVcf2 <- data.frame(tmpVcf@fix, stringsAsFactors = FALSE)
#   pass_filt <- which(tmpVcf2$FILTER == "PASS")
# 
#   GQ <- unlist(lapply(strsplit(tmpVcf@gt[,2], ":"), '[[', 2))[pass_filt]
#   tmpVcf2 <- tmpVcf2[pass_filt,]
# 
#   tmpVcf3 <- NULL
#   strSplitRes <- strsplit(tmpVcf2$INFO, ";")
#   for (j in 1:length(strSplitRes)) {
#     tmpVector <- unlist(strSplitRes[[j]])[grep("^AF|AO|DP|FR|FAO|FDP|FSAF|FSAR|HRUN",unlist(strSplitRes[[j]]))]
#     tmpVector <- tmpVector[-grep("FRO", tmpVector)]
#     tmpVector <- str_remove(tmpVector, ".*=")
#     tmpVector2 <- unlist(c(tmpVcf2[j,1:7], tmpVector))
#     tmpVcf3 <- rbind(tmpVcf3, tmpVector2)
#   }
#   tmpVcf3 <- data.frame(tmpVcf3, stringsAsFactors = FALSE)
#   tmpVcf3$GQ <- unlist(GQ)
#   rownames(tmpVcf3) <- NULL
#   colnames(tmpVcf3) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
#                          "FILT", "AF", "AO", "DP", "FAO", "FDP", "FR",
#                          "FSAF", "FSAR", "HRUN", "GQ")
#   tmpVcf3[,c(6, 8:12, 14:17)] <- lapply(tmpVcf3[,c(6, 8:12, 14:17)], as.numeric)
#   tmpVcf3 <- cbind("sample" = sampleName, tmpVcf3)
#   return(tmpVcf3)
# }
# 
# stopCluster(cl)
# 
# # ### run for seprate geno

# cl <- makeCluster(30)
# registerDoParallel(cl)
# 
# geno_vcf_processed <- NULL
# geno_vcf_processed <- foreach(i = genoVcfs, .combine = 'rbind', .packages = c('vcfR', 'stringr')) %dopar% {
#   sampleName <- str_remove(unlist(str_split(i, "/"))[7], ".vcf")
#   tmpVcf <- vcfR::read.vcfR(i)
#   tmpVcf2 <- data.frame(tmpVcf@fix, stringsAsFactors = FALSE)
#   # pass_filt <- which(tmpVcf2$FILTER == "PASS")
#   
#   GQ <- unlist(lapply(strsplit(tmpVcf@gt[,2], ":"), '[[', 2))
#   # GQ <- unlist(lapply(strsplit(tmpVcf@gt[,2], ":"), '[[', 2))[pass_filt]
#   # tmpVcf2 <- tmpVcf2[pass_filt,]
# 
#   tmpVcf2_posString <- paste0(tmpVcf2$CHROM, ":", tmpVcf2$POS)
#   GQ <- GQ[which(tmpVcf2_posString %in% positionString)]
#   tmpVcf2 <- tmpVcf2[which(tmpVcf2_posString %in% positionString),]
#   
#   tmpVcf3 <- NULL
#   strSplitRes <- strsplit(tmpVcf2$INFO, ";")
#   for (j in 1:length(strSplitRes)) {
#     tmpVector <- unlist(strSplitRes[[j]])[grep("^AF|AO|DP|FR|FAO|FDP|FSAF|FSAR|HRUN",unlist(strSplitRes[[j]]))]
#     tmpVector <- tmpVector[-grep("FRO", tmpVector)]
#     tmpVector <- str_remove(tmpVector, ".*=")
#     tmpVector2 <- unlist(c(tmpVcf2[j,1:7], tmpVector))
#     tmpVcf3 <- rbind(tmpVcf3, tmpVector2)
#   }
#   tmpVcf3 <- data.frame(tmpVcf3, stringsAsFactors = FALSE)
#   tmpVcf3$GQ <- unlist(GQ)
#   rownames(tmpVcf3) <- NULL
#   colnames(tmpVcf3) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
#                          "FILT", "AF", "AO", "DP", "FAO", "FDP", "FR",
#                          "FSAF", "FSAR", "HRUN", "GQ")
#   tmpVcf3[,c(6, 8:12, 14:17)] <- lapply(tmpVcf3[,c(6, 8:12, 14:17)], as.numeric)
#   tmpVcf3 <- cbind("sample" = sampleName, tmpVcf3)
#   return(tmpVcf3)
# }
# 
# stopCluster(cl)
# 
# 
# geno_vcf_processed$sample <- str_remove(geno_vcf_processed$sample, "_geno")
# 
# 
# write.table(geno_vcf_processed, "/br_z1/kevin_storage/misc/20211124Genopvars.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# write.table(all_vcf_processed, "/br_z1/kevin_storage/misc/20211020mgpvars.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# write.table(geno_vcf_processed, "/br_z1/kevin_storage/misc/2021102Genopvars.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

all_vcf_processed <- read.table("/br_z1/kevin_storage/misc/20211020mgpvars.txt", sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE)
geno_vcf_processed <- read.table("/br_z1/kevin_storage/misc/20211124Genopvars.txt", sep = "\t", header = TRUE,
                                 stringsAsFactors = FALSE)
geno_vcf_processed <- geno_vcf_processed[which(geno_vcf_processed$FILT == "PASS"), ]

all_vcf_processed2 <- rbind(all_vcf_processed, geno_vcf_processed)

snpFilters <- c("A", "T", "G", "C")
all_vcf_processed2 <- all_vcf_processed2[which(all_vcf_processed2$REF %in% snpFilters),]
all_vcf_processed2 <- all_vcf_processed2[which(all_vcf_processed2$ALT %in% snpFilters),]


non_zero <- all_vcf_processed2[-which(all_vcf_processed2$AF < 0.1 | all_vcf_processed2$AF > 0.9),]

# ggplot(data = non_zero) + geom_histogram(aes(x = AF), binwidth = 0.05, color = "white") + xlab("allele frequency") +
#   ggtitle("Histogram of het. AF") + theme(plot.title = element_text(hjust = 0.5))
# 
# hist(as.numeric(table(non_zero$sample)), breaks = seq(from=0, to=1500, by=10),
#      main = "Histogram of per sample het. SNPs (n = 278)", xlab = "# of het SNPs", ylab = "counts")
# 
# summary(as.numeric(table(non_zero$sample)))

non_zero$POS <- as.numeric(non_zero$POS)
non_zero$AF <- as.numeric(non_zero$AF)
non_zero$CHROM <- str_replace(non_zero$CHROM, "chrX", "chr20")

fdpFilt <- which(non_zero$FDP > 150)
faoFilt <- which(non_zero$FAO > 5)
freqFilt <- which(non_zero$AF > 0.05)
hrunFilt <- which(non_zero$HRUN < 4)
qualFilt <- which(non_zero$QUAL > 30)
gqFilt <- which(non_zero$GQ > 20)
strandRatio <- intersect(which(non_zero$FSAF/non_zero$FSAR > 0.2),
                         which(non_zero$FSAF/non_zero$FSAR < 5))
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, hrunFilt, gqFilt, qualFilt))
non_zero <- non_zero[goodSamps,]


### function to graph BAF + cn overlap 
###

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

freqPlot_cn <- function(df, main = NULL, chromTextSpec = NULL){
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
  
  df$loc.start <- df$loc.start/1e6
  df$loc.end <- df$loc.end/1e6
  df$col <- "#000000"
  df$col[which(df$seg.mean > 0.2)] <- "#FF0000"
  df$col[which(df$seg.mean < -0.2)] <- "#0000FF"
  
  for (i in unique(df$chrom)) {
    df$loc.start[which(df$chrom == i)] <- df$loc.start[which(df$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    df$loc.end[which(df$chrom == i)] <- df$loc.end[which(df$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df <- df[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(df) <- c("chrom", "start", "end","cn", "col")
  
  #df2 <- rbind(df2, dummyPoints)
  
  ggplot(df) + geom_hline(yintercept = c(-3, -2 , -1 , 0, 1, 2, 3), color = "#D4D4D4") + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = df$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1),
                                                           limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}


freqPlot_cn_amp <- function(df_cn, df_amp, main = NULL, chromTextSpec = NULL){
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
  
  df_cn$loc.start <- df_cn$loc.start/1e6
  df_cn$loc.end <- df_cn$loc.end/1e6
  df_cn$col <- "#000000"
  df_cn$col[which(df_cn$seg.mean > 0.2)] <- "#FF0000"
  df_cn$col[which(df_cn$seg.mean < -0.2)] <- "#0000FF"
  
  for (i in unique(df_cn$chrom)) {
    df_cn$loc.start[which(df_cn$chrom == i)] <- df_cn$loc.start[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    df_cn$loc.end[which(df_cn$chrom == i)] <- df_cn$loc.end[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df_cn <- df_cn[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(df_cn) <- c("chrom", "start", "end","cn", "col")
  
  
  df_amp$Pos <- ((df_amp$StartPos + df_amp$EndPos)/2)/1e6
  df_amp$col <- "#64e5d1"
  df_amp$col[grep("SNP", df_amp$Gene)] <- "#006400"
  
  for (i in unique(df_amp$ChromNum)) {
    df_amp$Pos[which(df_amp$ChromNum == i)] <- df_amp$Pos[which(df_amp$ChromNum == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df_amp <- df_amp[,c(1,6,5,7)]
  colnames(df_amp) <- c("chrom", "pos", "cn", "col")
  
  
  ggplot(df_cn) + geom_hline(yintercept = c(-3, -2 , -1 , 0, 1, 2, 3), color = "#D4D4D4") +
    geom_point(data = df_amp, aes(x = pos, y = cn),
               colour = df_amp$col, size = 0.05, inherit.aes = FALSE) +
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = df_cn$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1),
                                                           limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}


nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}


non_zero$sample <- nameStripper(non_zero$sample)
non_zero$sample <- str_remove(str_remove(non_zero$sample, "x.*"), "-")
non_zero$dupeString <- paste0(non_zero$sample, non_zero$CHROM, ":", non_zero$POS)
non_zero <- non_zero[-which(duplicated(non_zero$dupeString)),]
non_zero$string <- paste0(non_zero$CHROM, ":", non_zero$POS)
non_zero_noPosFilt <- non_zero

tableAf <- NULL
for (i in unique(non_zero$string)) {
  tmpTable <- non_zero[which(non_zero$string ==  i),]
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

filt <- which(tableAf$medianAf > 45 & tableAf$medianAf < 55 & tableAf$count > 1)
goodPosString  <- tableAf$position[filt]


non_zero <- non_zero[which(non_zero$string %in% goodPosString),]
non_zero$CHROM <- str_replace(non_zero$CHROM, "chrX", "chr20")
non_zero_noPosFilt$CHROM <- str_replace(non_zero_noPosFilt$CHROM, "chrX", "chr20")

samplesVcf <- c("2611lt", "2942rt", "2285lt", "2519lt", "10879lt")
#sampleList <- c("2611lt", "2611lt", "2942rt", "2285lt", "2519lt","10879lt") 
sampleList <- c("2611T", "2942T","2285T","2519T","10879T")
samplesSeg <- c("2611lt", "2942rt", "2285lt", "2519lt", "10879lt")

colnames(pcfSegs) <- c("ID", "chrom", "arm", "loc.start", "loc.end", "num.marks", "seg.mean")
# allSegs$chrom <- str_replace(allSegs$chrom, "23", "20")
pcfSegs$chrom <- str_replace(pcfSegs$chrom, "X", "20")
allRnaCbs$chrom <- str_replace(allRnaCbs$chrom, "X", "20")


for (i in seq_along(samplesVcf)) {
  testDf <- non_zero_16[which(non_zero_16$sample == samplesVcf[i]),]
  a <- freqPlot_baf(testDf)
  
  testDf_cn <- allSegs[which(allSegs$ID == samplesSeg[i]),]
  # b <- freqPlot_cn(testDf_cn, main = samplesSeg[i])
  
  testDf_pcf <- pcfSegs[which(pcfSegs$ID == sampleList[i]),]
  c <- freqPlot_cn(testDf_pcf)
  
  testDf_cbs <- allRnaCbs[which(allRnaCbs$ID == samplesSeg[i]),]
  d <- freqPlot_cn(testDf_cbs)
  
  testDf_amp <- cbind(allAmps[,c(3:6)], allAmps[[samplesVcf[i]]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = paste(samplesSeg[i], "MousePanel, RnaPcf, RnaSeg, Baf"))
  
  
  gA <- ggplotGrob(a)
  gC <- ggplotGrob(c)
  gD <- ggplotGrob(d)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211109_rnaComps", samplesSeg[i], ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, gC, gD, gA))
  dev.off()
  
}



### graph normals - check normal variability 

for (i in unique(segs_normal$ID)) {
  
  testDf_cn <- segs_normal[which(segs_normal$ID == i),]
  testDf_amp <- cbind(amps_normal[,c(3:6)], amps_normal[[i]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = i)
  
  testDf <- non_zero[which(non_zero$sample == i),]
  a <- freqPlot_baf(testDf)
  
  
  gA <- ggplotGrob(a)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211101_normals_amp", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, gA))
  dev.off()
  
}

for (i in unique(allSegs$ID)) {
  
  testDf_cn <- allSegs[which(allSegs$ID == i),]
  testDf_amp <- cbind(allAmps[,c(3:6)], allAmps[[i]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = i)
  
  testDf <- non_zero[which(non_zero$sample == i),]
  a <- freqPlot_baf(testDf)
  
  
  gA <- ggplotGrob(a)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211101_bafFilt_", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, gA))
  dev.off()
  
}


### getting correlations to changes

variants_grange <- unique(GRanges(seqnames = non_zero$CHROM,
                                  IRanges(start = non_zero$POS, end = non_zero$POS)))
segs_grange <- GRanges(seqnames = paste0("chr", allSegs$chrom),
                       IRanges(start = allSegs$loc.start,  end = allSegs$loc.end))

non_zero$absAf <- abs(non_zero$AF -  0.5)
allSegs$dipCount <- 2^allSegs$seg.mean * 2
allSegs$absDipCount <- abs(allSegs$dipCount - 2)

### how to deal with losses and gains in terms of correlation? convert the segments
### into absolute counts i.e reference of diploid


corTable <- NULL
for (i in 1:length(variants_grange)) {
  tmpGrange <- variants_grange[i]
  tmpSegs <- allSegs[queryHits(findOverlaps(segs_grange, tmpGrange)),]
  
  
  tmpPos <- paste0(variants_grange@seqnames[i], ":", variants_grange@ranges@start[i])
  tmpVcf <- non_zero[which(non_zero$string == tmpPos),]
  
  allNames <- intersect(tmpVcf$sample, tmpSegs$ID)
  
  tmpSegs <- tmpSegs[which(tmpSegs$ID %in% allNames),]
  tmpVcf <- tmpVcf[which(tmpVcf$sample %in% allNames),]
  
  tmpSegs <- tmpSegs[match(tmpVcf$sample, tmpSegs$ID),]
  
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
corTable <- corTable[-which(is.na(corTable$Med.Vf)),]

corTable.na <- corTable[which(is.na(corTable$Cor.p)),]
corTable <- corTable[-which(is.na(corTable$Cor.p)),]


### expect things to not shift if the seg mean is 0 for all samps i.e 
### cor is na
ggplot(data = corTable.na) + geom_point(aes(x = Med.Vf, y = Count.hgsc, color = Mad.Vf)) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, 5)) + 
  ggtitle("Stats on BAFs (n = 159) all samples with filter") +
  theme(plot.title = element_text(hjust = 0.5))



ggplot(data = corTable) + geom_point(aes(x = Med.AbsVf, y = Mad.AbsVf, color = Cor.p)) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 26, 1)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 26, 1)) + 
  ggtitle("Stats on BAFs (n = 159) all samples with filter") +
  theme(plot.title = element_text(hjust = 0.5))



ggplot(data = corTable) + geom_point(aes(x = Med.Vf, y = Mad.Vf, color = Cor.p)) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, 5)) + 
  ggtitle("Stats on BAFs (n = 159) all samples with filter") +
  theme(plot.title = element_text(hjust = 0.5))


badPos_CorNa <- corTable.na$Position[which(corTable.na$Med.Vf > 55 | corTable.na$Med.Vf < 45)]

badPos_Cor <- corTable$Position[which(corTable$Cor.p < 0.20)]
addBadPos <- c(badPos_CorNa, badPos_Cor)
non_zero_20 <- non_zero
non_zero_20 <- non_zero_20[-which(non_zero_20$string %in% addBadPos),]


badPos_Cor <- corTable$Position[which(corTable$Cor.p < 0.05)]
addBadPos <- c(badPos_CorNa, badPos_Cor)
non_zero_0 <- non_zero
non_zero_0 <- non_zero_0[-which(non_zero_0$string %in% addBadPos),]


badPos_Cor <- corTable$Position[which(corTable$Cor.p < 0.30)]
addBadPos <- c(badPos_CorNa, badPos_Cor)
non_zero_30 <- non_zero
non_zero_30 <- non_zero_30[-which(non_zero_30$string %in% addBadPos),]



for (i in unique(allSegs$ID)) {
  
  testDf_cn <- allSegs[which(allSegs$ID == i),]
  testDf_amp <- cbind(allAmps[,c(3:6)], allAmps[[i]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = paste0(i, ":cor cutoffs 0, 20 , 30"))
  
  testDf <- non_zero_30[which(non_zero_30$sample == i),]
  a_30 <- freqPlot_baf(testDf)
  
  testDf <- non_zero_20[which(non_zero_20$sample == i),]
  a_20 <- freqPlot_baf(testDf)
  
  testDf <- non_zero_0[which(non_zero_0$sample == i),]
  a_0 <- freqPlot_baf(testDf)
  
  
  g30 <- ggplotGrob(a_30)
  g20 <- ggplotGrob(a_20)
  g0 <- ggplotGrob(a_0)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211115_bafFilt_3_", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, g0, g20, g30))
  dev.off()
  
}


### bafs with the two types of filtering - (1) median, mad, min count (2) correlation 

i  <- "6941lt"
i <- "3867lt"
i <- "6786lt"
i <- "6260lt"
i <- "5438lt"
i <- "efd57"

allSegs$ID <- str_remove(allSegs$ID, "x.*")
colnames(allAmps)[7:ncol(allAmps)] <- str_remove(colnames(allAmps)[7:ncol(allAmps)] , "x.*")

for (i in unique(allSegs$ID)) {
  print(i)
  
  testDf_cn <- allSegs[which(allSegs$ID == i),]
  testDf_amp <- cbind(allAmps[,c(3:6)], allAmps[[i]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = paste0(i, "MousePanel, Baf: no filt, Baf w/ filts"))
  
  testDf <- non_zero_noPosFilt[which(non_zero_noPosFilt$sample == i),]
  a_nofilt <- freqPlot_baf(testDf)
  
  testDf <- non_zero_20[which(non_zero_20$sample == i),]
  a_20 <- freqPlot_baf(testDf)
  
  gnoFilt <- ggplotGrob(a_nofilt)
  g20 <- ggplotGrob(a_20)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211122_filtComps_20_", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, gnoFilt, g20))
  dev.off()
  
}



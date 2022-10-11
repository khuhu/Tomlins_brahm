library(gridExtra)


nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}

all_vcf_processed <- read.table("/br_z1/kevin_storage/misc/20211020mgpvars.txt", sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE)
geno_vcf_processed <- read.table("/br_z1/kevin_storage/misc/2021102Genopvars.txt", sep = "\t", header = TRUE,
                                 stringsAsFactors = FALSE)

bedFile <- read.table("/br_z1/kevin_storage/misc/IAD202670_167_Designed.gc.bed", sep = "\t", stringsAsFactors = FALSE, header = FALSE)

all_vcf_processed2 <- rbind(all_vcf_processed, geno_vcf_processed)

snpFilters <- c("A", "T", "G", "C")
all_vcf_processed2 <- all_vcf_processed2[which(all_vcf_processed2$REF %in% snpFilters),]
all_vcf_processed2 <- all_vcf_processed2[which(all_vcf_processed2$ALT %in% snpFilters),]


non_zero <- all_vcf_processed2[-which(all_vcf_processed2$AF < 0.1 | all_vcf_processed2$AF > 0.9),]

non_zero$POS <- as.numeric(non_zero$POS)
non_zero$AF <- as.numeric(non_zero$AF)

non_zero$sample <- nameStripper(non_zero$sample)
non_zero$sample <- str_remove(str_remove(non_zero$sample, "x.*"), "-")


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
non_zero$string <- paste0(non_zero$CHROM, ":", non_zero$POS)
non_zero$dupeString <- paste0(non_zero$sample, non_zero$CHROM, ":", non_zero$POS)
non_zero <- non_zero[-which(duplicated(non_zero$dupeString)),]


mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

mouseNormal <- nameStripper(mouseNormal)


# normal_non_zero <- non_zero[which(non_zero$sample %in% mouseNormal),]
# badPos <- names(table(normal_non_zero$string[which(normal_non_zero$AF > 0.6 | normal_non_zero$AF < 0.4)]))
# non_zero2 <- non_zero[-which(non_zero$string %in% badPos),]


### removing the more variable SNPs from the normal positions doesn't seem to shift the AF
### there is a bias towards shifts in higher AF - should double check see two graphs
# a <- ggplot(data = non_zero) + geom_histogram(aes(x = AF), binwidth = 0.05, color = "white") + xlab("allele frequency") +
#   ggtitle("All Het. AF") + theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_x_continuous(limits = c(0,1), breaks = seq(0,1, 0.05))
# 
# b <- ggplot(data = non_zero2) + geom_histogram(aes(x = AF), binwidth = 0.05, color = "white") + xlab("allele frequency") +
#   ggtitle("Het. AF (- bad variable position in normals)") + theme(plot.title = element_text(hjust = 0.5)) +
#   scale_x_continuous(limits = c(0,1), breaks = seq(0,1, 0.05))
# 
# grid.arrange(a,b, ncol = 1)

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

### cross-referencing with large changes: making df of it so i can look back
samp2519lt <- non_zero$string[which(non_zero$sample == "2519lt" & non_zero$CHROM == "chr15")]
samp2942rt <- non_zero$string[which(non_zero$sample == "2942rt" & non_zero$CHROM == "chr16")]
samp1531lt <- non_zero$string[which(non_zero$sample == "1531lt" & non_zero$CHROM == "chr7")]
samp1549rt <- non_zero$string[which(non_zero$sample == "1549rt" & non_zero$CHROM == "chr1")]
samp1628lt <- non_zero$string[which(non_zero$sample == "1628lt" & non_zero$CHROM == "chr4")]

samp2007rt <- non_zero$string[which((non_zero$sample == "2007rt" & non_zero$CHROM == "chr8") | 
                                      (non_zero$sample == "2007rt" & non_zero$CHROM == "chr11") |
                                      (non_zero$sample == "2007rt" &non_zero$CHROM == "chr13"))]
samp2027lts <- non_zero$string[which((non_zero$sample == "2027lts" & non_zero$CHROM == "chr4") | 
                                      (non_zero$sample == "2027lts" & non_zero$CHROM == "chr7") |
                                      (non_zero$sample == "2027lts" &non_zero$CHROM == "chr13"))]

samp2796rt <- non_zero$string[which((non_zero$sample == "2796rt" & non_zero$CHROM == "chr18") | 
                                       (non_zero$sample == "2796rt" & non_zero$CHROM == "chr19"))]

samp5293lt <- non_zero$string[which((non_zero$sample == "5293lt" & non_zero$CHROM == "chr3") | 
                                      (non_zero$sample == "5293lt" & non_zero$CHROM == "chr18"))]

samp6786lt <- non_zero$string[which((non_zero$sample == "6786lt" & non_zero$CHROM == "chr2") | 
                                      (non_zero$sample == "6786lt" & non_zero$CHROM == "chr12") |
                                      (non_zero$sample == "6786lt" & non_zero$CHROM == "chr13"))]

varPosFilt <- c(samp2519lt, samp2942rt, samp1531lt, samp1549rt, samp1628lt, samp2007rt, samp2027lts,
                samp2796rt, samp5293lt, samp6786lt)
tableAf_filt <- tableAf[which(tableAf$position %in% varPosFilt),]

filt <- which(tableAf$medianAf > 45 & tableAf$medianAf < 55 & tableAf$count > 30 & tableAf$madAf > 5)
table_filt2 <- tableAf[filt,]



a <- ggplot(data = tableAf) + geom_point(aes(x = medianAf, y = count, color = madAf)) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,250, 10)) +
  ggtitle("Stats on BAFs (n = 1574) all samples") +
  theme(plot.title = element_text(hjust = 0.5))

b <- ggplot(data = tableAf_filt) + geom_point(aes(x = medianAf, y = count, color = madAf)) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,250, 10)) +
  ggtitle("Stats on BAFs (n = 182) from 10 samples") +
  theme(plot.title = element_text(hjust = 0.5))

c <- ggplot(data = table_filt2) + geom_point(aes(x = medianAf, y = count, color = madAf)) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100, 10)) +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,250, 10)) + 
  ggtitle("Stats on BAFs (n = 802) all samples with filter") +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(a, b, c, ncol = 1)


dev.off()
pdf("/br_z1/kevin_storage/misc/20211102bafMedianProfiles.pdf", height = 10)
grid.arrange(a, b, c, ncol = 1)
dev.off()

### mad does seemm like a good measurement 
### what I want to get rid of: 
### (1) SNPs that are called less than 5 times 
### (2) more than 5 instances, higher or low AF but low MAD - systematic error 
### (3) high counts, within het range with low mad - meaning no matter what. it's always het mad of 5 - based on quantiles


goodPosString  <- tableAf$position[filt]

non_zero_filt <- non_zero[which(non_zero$string %in% goodPosString), ]


ampGrange <- GRanges(seqnames = bedFile$V1,
                     IRanges(start = bedFile$V2, end = bedFile$V3))
varGrange <- GRanges(seqnames = non_zero_filt$CHROM,
                     IRanges(start = non_zero_filt$POS, end = non_zero_filt$POS))

table(queryHits(findOverlaps(ampGrange, varGrange)))
subjectHits(findOverlaps(ampGrange, varGrange))

View(non_zero_filt[subjectHits(findOverlaps(ampGrange, varGrange)), ])

non_zero_filt$snpClusterAmp <- "0"
non_zero_filt$snpClusterAmp <- queryHits(findOverlaps(ampGrange, varGrange))


which(duplicated(subjectHits(findOverlaps(ampGrange, varGrange))))

View(non_zero_filt[which(duplicated(subjectHits(findOverlaps(ampGrange, varGrange)))),])


overlapRes <- findOverlaps(ampGrange, varGrange)
overlapRes <- overlapRes[-which(duplicated(subjectHits(overlapRes)))]

non_zero_filt$snpClusterAmp <- "0"
non_zero_filt$snpClusterAmp <- queryHits(overlapRes)[order(subjectHits(overlapRes))]

non_zero_filt$ampSample <- paste0(non_zero_filt$sample, ":",non_zero_filt$snpClusterAmp)


### no it doesn't seem to fluctuate much for SNPs on the same amplicon, the clusters though are one
### chromosomes 3, 4 , 8 - obvious ones - i.e clusters with variabilty 



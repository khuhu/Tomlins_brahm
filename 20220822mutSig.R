# install.packages("mutSignatures")
# install.packages("kableExtra")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(dplyr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSignatures)
library(vcfR)
library(stringr)
library(BSgenome.Mmusculus.UCSC.mm10)

hg19 <- BSgenome.Hsapiens.UCSC.hg19
mm10 <- BSgenome.Mmusculus.UCSC.mm10

### get AF and other stats to filter out below
### only .156 Mb covering CDS regions of genes (excluding 400 markers for SNPs)

singleBase <- c("A", "T", "G", "C")

combinedVars  <- read.table("/avatar_data6/mouseData/reportAnno/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562_anno.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

combinedVars_reseq  <- read.table("/avatar_data6/mouseData/reportAnno/Auto_user_AUS5-260-BBN_mouse_bladder_MG_2_514_605_anno.txt", sep = "\t",
                                  header = TRUE, stringsAsFactors = FALSE)

combinedVars$Sample2 <- str_remove(combinedVars$Sample, "\\_.*")
combinedVars_reseq$Sample2 <- str_remove(combinedVars_reseq$Sample, "\\_.*")

combinedVars$string <- paste(combinedVars$Sample2, combinedVars$Chr, combinedVars$Start, combinedVars$End)
combinedVars_reseq$string <- paste(combinedVars_reseq$Sample2, combinedVars_reseq$Chr, combinedVars_reseq$Start, combinedVars_reseq$End)

combinedVars_concord <- combinedVars[which(combinedVars$string %in%  combinedVars_reseq$string),]
combinedVars_concord <- combinedVars_concord[which(combinedVars_concord$mm10_mpgpv6_Indels == ""),]
combinedVars_concord$dupeString <- paste(combinedVars_concord$Chr, combinedVars_concord$Start, combinedVars_concord$End)
combinedVars_concord <- combinedVars_concord[-which(duplicated(combinedVars_concord$dupeString)),]

fdpFilt <- which(combinedVars_concord$FDP > 100)
faoFilt <- which(combinedVars_concord$FAO > 10)
freqFilt <- which(combinedVars_concord$AF > 0.05)
hrunFilt <- which(combinedVars_concord$HRUN < 4)
strandRatio <- intersect(which(combinedVars_concord$FSAF/combinedVars_concord$FSAR > 0.2),
                         which(combinedVars_concord$FSAF/combinedVars_concord$FSAR < 5))
qualFilt <- which(combinedVars_concord$QUAL >= 50)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt, hrunFilt))
combinedVars_concord_goodsamps <- combinedVars_concord[goodSamps,]
combinedVars_concord_goodsamps <- combinedVars_concord_goodsamps[which(combinedVars_concord_goodsamps$Ref %in% singleBase),]
combinedVars_concord_goodsamps <- combinedVars_concord_goodsamps[which(combinedVars_concord_goodsamps$Alt %in% singleBase),]

combinedVars_concord_goodsamps_exon <- combinedVars_concord_goodsamps[which(combinedVars_concord_goodsamps$Func.refGene == "exonic"), ]
combinedVars_concord_goodsamps_exon$type <- paste0(combinedVars_concord_goodsamps_exon$Ref, ">", combinedVars_concord_goodsamps_exon$Alt)
# combinedVars_concord_goodsamps_exon_put <- combinedVars_concord_goodsamps_exon[-which(combinedVars_concord_goodsamps_exon$ExonicFunc.refGene == "synonymous SNV"),]

combinedVars_concord_goodsamps_exon$string <- paste0(combinedVars_concord_goodsamps_exon$Sample2,
                                                     combinedVars_concord_goodsamps_exon$Chr,
                                                     combinedVars_concord_goodsamps_exon$Start)

combinedVars_concord_goodsamps_exon$type <- paste(combinedVars_concord_goodsamps_exon$Ref, ">", combinedVars_concord_goodsamps_exon$Alt)


combinedVars_concord_goodsamps$string  <- paste0(combinedVars_concord_goodsamps$Sample2,
                                                 combinedVars_concord_goodsamps$Chr,
                                                 combinedVars_concord_goodsamps$Start)

# combinedVars_concord_goodsamps_nonsyn <- combinedVars_concord_goodsamps[which(combinedVars_concord_goodsamps$ExonicFunc.refGene == "nonsynonymous SNV"), ]
# combinedVars_concord_goodsamps_nonsyn$string <- paste0(combinedVars_concord_goodsamps_nonsyn$Sample, combinedVars_concord_goodsamps_nonsyn$Chr,
#                                                combinedVars_concord_goodsamps_nonsyn$Start)
# combinedVars_concord_goodsamps_nonsyn$type <- paste0(combinedVars_concord_goodsamps_nonsyn$Ref, ">", combinedVars_concord_goodsamps_nonsyn$Alt)


### discordant variants

combinedVars_discord <- rbind(combinedVars[-which(combinedVars$string %in%  combinedVars_reseq$string),],
                              combinedVars_reseq[-which(combinedVars_reseq$string %in% combinedVars$string),])
combinedVars_discord <- combinedVars_discord[which(combinedVars_discord$mm10_mpgpv6_Indels == ""),]
combinedVars_discord$dupeString <- paste(combinedVars_discord$Chr, combinedVars_discord$Start, combinedVars_discord$End)
combinedVars_discord <- combinedVars_discord[-which(duplicated(combinedVars_discord$dupeString)),]

fdpFilt <- which(combinedVars_discord$FDP > 100)
faoFilt <- which(combinedVars_discord$FAO > 10)
freqFilt <- which(combinedVars_discord$AF > 0.05)
hrunFilt <- which(combinedVars_discord$HRUN < 4)
strandRatio <- intersect(which(combinedVars_discord$FSAF/combinedVars_discord$FSAR > 0.2),
                         which(combinedVars_discord$FSAF/combinedVars_discord$FSAR < 5))
qualFilt <- which(combinedVars_discord$QUAL >= 50)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt, hrunFilt))
combinedVars_discord_goodsamps <- combinedVars_discord[goodSamps,]
combinedVars_discord_goodsamps <- combinedVars_discord_goodsamps[which(combinedVars_discord_goodsamps$Ref %in% singleBase),]
combinedVars_discord_goodsamps <- combinedVars_discord_goodsamps[which(combinedVars_discord_goodsamps$Alt %in% singleBase),]

combinedVars_discord_goodsamps_exon <- combinedVars_discord_goodsamps[which(combinedVars_discord_goodsamps$Func.refGene == "exonic"), ]
combinedVars_discord_goodsamps_exon$type <- paste0(combinedVars_discord_goodsamps_exon$Ref, ">", combinedVars_discord_goodsamps_exon$Alt)
combinedVars_discord_goodsamps_exon$string <- paste0(combinedVars_discord_goodsamps_exon$Sample, combinedVars_discord_goodsamps_exon$Chr,
                                                     combinedVars_discord_goodsamps_exon$Start)

# combinedVars_discord_goodsamps_exon_put <- combinedVars_discord_goodsamps_exon[-which(combinedVars_discord_goodsamps_exon$ExonicFunc.refGene == "synonymous SNV"),]
# combinedVars_discord_goodsamps_nonsyn <- combinedVars_discord_goodsamps[which(combinedVars_discord_goodsamps$ExonicFunc.refGene == "nonsynonymous SNV"), ]
# combinedVars_discord_goodsamps_nonsyn$string <- paste0(combinedVars_discord_goodsamps_nonsyn$Sample, combinedVars_discord_goodsamps_nonsyn$Chr,
#                                                        combinedVars_discord_goodsamps_nonsyn$Start)
# combinedVars_discord_goodsamps_nonsyn$type <- paste0(combinedVars_discord_goodsamps_nonsyn$Ref, ">", combinedVars_discord_goodsamps_nonsyn$Alt)

### redid how I went about initial filtering

combinedVarsV2 <- combinedVars
combinedVarsV2 <- combinedVarsV2[which(combinedVarsV2$mm10_mpgpv6_Indels == ""),]
combinedVarsV2$dupeString <- paste(combinedVarsV2$Chr, combinedVarsV2$Start, combinedVarsV2$End)
combinedVarsV2 <- combinedVarsV2[-which(duplicated(combinedVarsV2$dupeString)),]

fdpFilt <- which(combinedVarsV2$FDP > 100)
faoFilt <- which(combinedVarsV2$FAO > 10)
freqFilt <- which(combinedVarsV2$AF > 0.05)
hrunFilt <- which(combinedVarsV2$HRUN < 4)
strandRatio <- intersect(which(combinedVarsV2$FSAF/combinedVarsV2$FSAR > 0.2),
                         which(combinedVarsV2$FSAF/combinedVarsV2$FSAR < 5))
qualFilt <- which(combinedVarsV2$QUAL >= 50)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt, hrunFilt))
combinedVarsV2_goodsamps <- combinedVarsV2[goodSamps,]
combinedVarsV2_goodsamps <- combinedVarsV2_goodsamps[which(combinedVarsV2_goodsamps$Ref %in% singleBase),]
combinedVarsV2_goodsamps <- combinedVarsV2_goodsamps[which(combinedVarsV2_goodsamps$Alt %in% singleBase),]


combinedVarsV2_goodsamps_exon <- combinedVarsV2_goodsamps[which(combinedVarsV2_goodsamps$Func.refGene == "exonic"), ]
combinedVarsV2_goodsamps_exon$type <- paste0(combinedVarsV2_goodsamps_exon$Ref, ">", combinedVarsV2_goodsamps_exon$Alt)
combinedVarsV2_goodsamps_exon$string <- paste0(combinedVarsV2_goodsamps_exon$Sample2, combinedVarsV2_goodsamps_exon$Chr,
                                               combinedVarsV2_goodsamps_exon$Start)

# combinedVarsV2_goodsamps_exon_put <- combinedVarsV2_goodsamps_exon[-which(combinedVarsV2_goodsamps_exon$ExonicFunc.refGene == "synonymous SNV"),]

# combinedVarsV2_goodsamps_nonsyn <- combinedVarsV2_goodsamps[which(combinedVarsV2_goodsamps$ExonicFunc.refGene == "nonsynonymous SNV"), ]
# combinedVarsV2_goodsamps_nonsyn$string <- paste0(combinedVarsV2_goodsamps_nonsyn$Sample, combinedVarsV2_goodsamps_nonsyn$Chr,
#                                                  combinedVarsV2_goodsamps_nonsyn$Start)
# combinedVarsV2_goodsamps_nonsyn$type <- paste0(combinedVarsV2_goodsamps_nonsyn$Ref, ">", combinedVarsV2_goodsamps_nonsyn$Alt)


combinedVars_reseqV2 <- combinedVars_reseq
combinedVars_reseqV2 <- combinedVars_reseqV2[which(combinedVars_reseqV2$mm10_mpgpv6_Indels == ""),]
combinedVars_reseqV2$dupeString <- paste(combinedVars_reseqV2$Chr, combinedVars_reseqV2$Start, combinedVars_reseqV2$End)
combinedVars_reseqV2 <- combinedVars_reseqV2[-which(duplicated(combinedVars_reseqV2$dupeString)),]

fdpFilt <- which(combinedVars_reseqV2$FDP > 100)
faoFilt <- which(combinedVars_reseqV2$FAO > 10)
freqFilt <- which(combinedVars_reseqV2$AF > 0.05)
hrunFilt <- which(combinedVars_reseqV2$HRUN < 4)
strandRatio <- intersect(which(combinedVars_reseqV2$FSAF/combinedVars_reseqV2$FSAR > 0.2),
                         which(combinedVars_reseqV2$FSAF/combinedVars_reseqV2$FSAR < 5))
qualFilt <- which(combinedVars_reseqV2$QUAL >= 50)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt, hrunFilt))
combinedVars_reseqV2_goodsamps <- combinedVars_reseqV2[goodSamps,]
combinedVars_reseqV2_goodsamps <- combinedVars_reseqV2_goodsamps[which(combinedVars_reseqV2_goodsamps$Ref %in% singleBase),]
combinedVars_reseqV2_goodsamps <- combinedVars_reseqV2_goodsamps[which(combinedVars_reseqV2_goodsamps$Alt %in% singleBase),]


combinedVars_reseqV2_goodsamps_exon <- combinedVars_reseqV2_goodsamps[which(combinedVars_reseqV2_goodsamps$Func.refGene == "exonic"), ]
combinedVars_reseqV2_goodsamps_exon$type <- paste0(combinedVars_reseqV2_goodsamps_exon$Ref, ">", combinedVars_reseqV2_goodsamps_exon$Alt)
combinedVars_reseqV2_goodsamps_exon_put <- combinedVars_reseqV2_goodsamps_exon[-which(combinedVars_reseqV2_goodsamps_exon$ExonicFunc.refGene == "synonymous SNV"),]

combinedVars_reseqV2_goodsamps_exon$string <- paste0(combinedVars_reseqV2_goodsamps_exon$Sample2,
                                                     combinedVars_reseqV2_goodsamps_exon$Chr, 
                                                     combinedVars_reseqV2_goodsamps_exon$Start)

# combinedVars_reseqV2_goodsamps_nonsyn <- combinedVars_reseqV2_goodsamps[which(combinedVars_reseqV2_goodsamps$ExonicFunc.refGene == "nonsynonymous SNV"), ]
# # combinedVars_reseqV2_goodsamps_nonsyn$string <- paste0(combinedVars_reseqV2_goodsamps_nonsyn$Sample, combinedVars_reseqV2_goodsamps_nonsyn$Chr,
#                                                  combinedVars_reseqV2_goodsamps_nonsyn$Start)
# 
# combinedVars_reseqV2_goodsamps_nonsyn$type <- paste0(combinedVars_reseqV2_goodsamps_nonsyn$Ref, ">", combinedVars_reseqV2_goodsamps_nonsyn$Alt)
# 


### first seq

tmpVcf <- read.table("/avatar_data3/eros_tmp/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562/plugin_out/variantCaller_out.1152/R_2022_07_01_11_06_33_user_AUS5-239-BBN_mouse_bladder_MG.xls",
                     sep = "\t", header = TRUE)
tmpVcf$Sample.Name2 <- str_remove(tmpVcf$Sample.Name, "\\_.*")
tmpVcf$string <- paste0(tmpVcf$Sample.Name2, tmpVcf$Chrom, tmpVcf$Position)
tmpVcf <- tmpVcf[which(tmpVcf$string %in% combinedVarsV2_goodsamps_exon$string), ]
vcfRate <- tmpVcf[, c("Chrom", "Position", "Ref", "Variant", "Allele.Call", "Filter", "Sample.Name", "Quality", "Frequency")]
colnames(vcfRate) <- c("CHROM", "POS", "REF", "ALT", "XTR1", "FILTER", "SAMPLEID", "QUAL", "FREQ")
vcfRate$XTR1 <- str_replace_all(vcfRate$XTR1, "Heterozygous", "0/1")
vcfRate$XTR1 <- str_replace_all(vcfRate$XTR1, "Homozygous", "1/1")
vcfRate$ID <- "."
vcfRate$INFO <- "."
vcfRate$FILTER <- "."
vcfRate$FORMAT <- "GT:PL"
vcfRate2 <- vcfRate[, c("CHROM", "POS", "ID", "REF", "ALT",
                        "QUAL", "FILTER", "INFO", "FORMAT",
                        "XTR1", "SAMPLEID", "FREQ")]
vcfRate3 <- filterSNV(dataSet = vcfRate2,  seq_colNames = c("REF", "ALT"))
vcfRate4 <- attachContext(mutData = vcfRate3,
                          chr_colName = "CHROM",
                          start_colName = "POS",
                          end_colName = "POS",
                          nucl_contextN = 3,
                          BSGenomeDb = mm10)
vcfRate4 <- removeMismatchMut(mutData = vcfRate4,
                              refMut_colName = "REF",
                              context_colName = "context",
                              refMut_format = "N")
vcfRate4 <- attachMutType(mutData = vcfRate4,
                          ref_colName = "REF",
                          var_colName = "ALT",
                          context_colName = "context")
vcfRate4$mutTypeSingle <- str_remove(str_remove(vcfRate4$mutType, ".*\\["), "\\].*")
vcfRate4$SAMPLEID <- str_remove(vcfRate$SAMPLEID, "\\_.*")
vcfRate4$string <- paste(vcfRate4$SAMPLEID, vcfRate4$CHROM, vcfRate4$POS, vcfRate4$POS)
vcfRate4$string2 <- paste0(vcfRate4$SAMPLEID, vcfRate4$CHROM, vcfRate4$POS)

### reseq

tmpVcf2 <- read.table("/avatar_data3/eros_tmp/Auto_user_AUS5-260-BBN_mouse_bladder_MG_2_514_605/plugin_out/variantCaller_out.1254/R_2022_09_13_14_13_30_user_AUS5-260-BBN_mouse_bladder_MG_2.xls",
                     sep = "\t", header = TRUE)
tmpVcf2$Sample.Name2 <- str_remove(tmpVcf2$Sample.Name, "\\_.*")
tmpVcf2$string <- paste0(tmpVcf2$Sample.Name2, tmpVcf2$Chrom, tmpVcf2$Position)
tmpVcf2 <- tmpVcf2[which(tmpVcf2$string %in% combinedVars_reseqV2_goodsamps_exon$string), ]
vcfRateReseq <- tmpVcf2[, c("Chrom", "Position", "Ref", "Variant", "Allele.Call", "Filter", "Sample.Name", "Quality", "Frequency")]
colnames(vcfRateReseq) <- c("CHROM", "POS", "REF", "ALT", "XTR1", "FILTER", "SAMPLEID", "QUAL", "FREQ")
vcfRateReseq$XTR1 <- str_replace_all(vcfRateReseq$XTR1, "Heterozygous", "0/1")
vcfRateReseq$XTR1 <- str_replace_all(vcfRateReseq$XTR1, "Homozygous", "1/1")
vcfRateReseq$ID <- "."
vcfRateReseq$INFO <- "."
vcfRateReseq$FILTER <- "."
vcfRateReseq$FORMAT <- "GT:PL"
vcfRateReseq2 <- vcfRateReseq[, c("CHROM", "POS", "ID", "REF", "ALT",
                        "QUAL", "FILTER", "INFO", "FORMAT",
                        "XTR1", "SAMPLEID", "FREQ")]
vcfRateReseq3 <- filterSNV(dataSet = vcfRateReseq2,  seq_colNames = c("REF", "ALT"))
vcfRateReseq4 <- attachContext(mutData = vcfRateReseq3,
                          chr_colName = "CHROM",
                          start_colName = "POS",
                          end_colName = "POS",
                          nucl_contextN = 3,
                          BSGenomeDb = mm10)
vcfRateReseq4 <- removeMismatchMut(mutData = vcfRateReseq4,
                              refMut_colName = "REF",
                              context_colName = "context",
                              refMut_format = "N")
vcfRateReseq4 <- attachMutType(mutData = vcfRateReseq4,
                          ref_colName = "REF",
                          var_colName = "ALT",
                          context_colName = "context")
vcfRateReseq4$mutTypeSingle <- str_remove(str_remove(vcfRateReseq4$mutType, ".*\\["), "\\].*")
vcfRateReseq4$SAMPLEID <- str_remove(vcfRateReseq4$SAMPLEID, "\\_.*")
vcfRateReseq4$string <- paste(vcfRateReseq4$SAMPLEID, vcfRateReseq4$CHROM, vcfRateReseq4$POS, vcfRateReseq4$POS)
vcfRateReseq4$string2 <- paste(vcfRateReseq4$SAMPLEID, vcfRateReseq4$CHROM, vcfRateReseq4$POS, vcfRateReseq4$POS)

### concord

vcfRateConcord <- vcfRate4[which(vcfRate4$string2 %in% combinedVars_concord_goodsamps_exon$string),]

### discord

vcfRateDiscord <- rbind(vcfRate4[which(vcfRate4$string2 %in% combinedVars_discord$string),],
                        vcfRateReseq4[which(vcfRateReseq4$string2 %in% combinedVars_discord$string),])


mutTypeCount <- NULL
for (i in unique(vcfRate4$SAMPLEID)) {
  tmpDf <-  vcfRate4[which(vcfRate4$SAMPLEID == i),]
  tmpDf2 <- data.frame("sample" = i, table(tmpDf$mutTypeSingle))
  mutTypeCount <- rbind(mutTypeCount, tmpDf2)
}

mutTypeCount2 <- NULL
for (i in unique(vcfRateReseq4$SAMPLEID)) {
  tmpDf <-  vcfRateReseq4[which(vcfRateReseq4$SAMPLEID == i),]
  tmpDf2 <- data.frame("sample" = i, table(tmpDf$mutTypeSingle))
  mutTypeCount2 <- rbind(mutTypeCount2, tmpDf2)
}

mutTypeCount3 <- NULL
for (i in unique(vcfRateConcord$SAMPLEID)) {
  tmpDf <-  vcfRateConcord[which(vcfRateConcord$SAMPLEID == i),]
  tmpDf2 <- data.frame("sample" = i, table(tmpDf$mutTypeSingle))
  mutTypeCount3 <- rbind(mutTypeCount3, tmpDf2)
}


dummyFactor <- data.frame("sample" = rep("dummy", 6), "Var1" = as.factor(levels(mutTypeCount$Var1)), "Freq" = rep(1,6))
mutTypeCount4 <- dummyFactor
for (i in unique(vcfRateDiscord$SAMPLEID)) {
  tmpDf <-  vcfRateDiscord[which(vcfRateDiscord$SAMPLEID == i),]
  tmpDf2 <- data.frame("sample" = i, table(tmpDf$mutTypeSingle))
  mutTypeCount4 <- rbind(mutTypeCount4, tmpDf2)
}

levels(mutTypeCount$Var1)
levels(mutTypeCount2$Var1)
levels(mutTypeCount3$Var1)
levels(mutTypeCount4$Var1)

mutTypeCount4$Var1 <- factor(mutTypeCount4$Var1, levels = c("C>A", "C>T", "T>A", "T>C", "C>G", "T>G"))

levels(mutTypeCount4$Var1)


mutTypeCount$sample <- str_remove(mutTypeCount$sample, "\\_.*")
mutTypeCount2$sample <- str_remove(mutTypeCount2$sample, "\\_.*")
mutTypeCount3$sample <- str_remove(mutTypeCount3$sample, "\\_.*")
mutTypeCount4$sample <- str_remove(mutTypeCount4$sample, "\\_.*")


a <- ggplot(mutTypeCount, aes(fill=Var1, y=Freq, x=sample)) +
  geom_bar(position="stack", stat="identity") + ylab("# of mutations") +
  ggtitle("nonsyn SBS seq1 (0.156Mb)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) + ylim(0, 120)



b <- ggplot(mutTypeCount2, aes(fill=Var1, y=Freq, x=sample)) +
  geom_bar(position="stack", stat="identity") + ylab("# of mutations") +
  ggtitle("nonsyn SBS seq2 (0.156Mb))") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) + ylim(0, 120)



c <- ggplot(mutTypeCount3, aes(fill=Var1, y=Freq, x=sample)) +
  geom_bar(position="stack", stat="identity") + ylab("# of mutations") +
  ggtitle("nonsyn SBS Concord (0.156Mb))") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) + ylim(0, 120)



d <- ggplot(mutTypeCount4, aes(fill=Var1, y=Freq, x=sample)) +
  geom_bar(position="stack", stat="identity") + ylab("# of mutations") +
  ggtitle("nonsyn SBS Discord (0.156Mb))") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) + ylim(0, 120)



grid.arrange(a,b,c,d, ncol = 2)

### actual context
tmpVcf <- read.table("/avatar_data3/eros_tmp/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562/plugin_out/variantCaller_out.1152/R_2022_07_01_11_06_33_user_AUS5-239-BBN_mouse_bladder_MG.xls",
                     sep = "\t", header = TRUE)

tmpVcf$Sample.Name <- str_remove(tmpVcf$Sample.Name, "\\_.*")
tmpVcf$string <- paste0(tmpVcf$Sample.Name, tmpVcf$Chrom, tmpVcf$Position)
tmpVcf <- tmpVcf[which(tmpVcf$string %in% combinedVarsV2_goodsamps_exon$string), ]


tmpVcf2 <- tmpVcf[, c("Chrom", "Position", "Ref", "Variant", "Allele.Call", "Filter", "Sample.Name", "Quality")]
colnames(tmpVcf2) <- c("CHROM", "POS", "REF", "ALT", "XTR1", "FILTER", "SAMPLEID", "QUAL")

tmpVcf2$XTR1 <- str_replace_all(tmpVcf2$XTR1, "Heterozygous", "0/1")
tmpVcf2$XTR1 <- str_replace_all(tmpVcf2$XTR1, "Homozygous", "1/1")
tmpVcf2$ID <- "."
tmpVcf2$INFO <- "."
tmpVcf2$FILTER <- "."
tmpVcf2$FORMAT <- "GT:PL"

tmpVcf3 <- tmpVcf2[, c("CHROM", "POS", "ID", "REF", "ALT",
                       "QUAL", "FILTER", "INFO", "FORMAT", 
                       "XTR1", "SAMPLEID")]

tmpVcf3 <- filterSNV(dataSet = tmpVcf3,  seq_colNames = c("REF", "ALT"))

tmpVcf4 <- attachContext(mutData = tmpVcf3,
                         chr_colName = "CHROM",
                         start_colName = "POS",
                         end_colName = "POS",
                         nucl_contextN = 3,
                         BSGenomeDb = mm10)


tmpVcf4 <- removeMismatchMut(mutData = tmpVcf4,                  
                       refMut_colName = "REF",
                       context_colName = "context",
                       refMut_format = "N")

tmpVcf4 <- attachMutType(mutData = tmpVcf4,            
                   ref_colName = "REF",              
                   var_colName = "ALT",             
                   context_colName = "context") 

head(tmpVcf4) %>% kable() %>% kable_styling(bootstrap_options = "striped")

allVcfCounts <- countMutTypes(mutTable = tmpVcf4,
                              mutType_colName = "mutType",
                              sample_colName = "SAMPLEID")

# print(allVcfCounts)
# 
# 
# vcf.params <- 
#   mutSignatures::setMutClusterParams( 
#     num_processesToExtract = 7,
#     num_totIterations = 1000,
#     num_parallelCores = 28) 
# 


# vcf.analysis <- 
#   decipherMutationalProcesses(input = allVcfCounts,
#                               params = vcf.params)


# Retrieve signatures (results)
vcf.sig <- vcf.analysis$Results$signatures

# Retrieve exposures (results)
vcf.exp <- vcf.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(vcf.sig, signature = 1, ylim = c(0, 0.10))

msigPlot(vcf.exp, main = "BLCA samples") + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))



cosmixSigs <- mutSigData$blcaSIGS %>% 
  dplyr::select(starts_with("COSMIC")) %>% 
  as.mutation.signatures()

blcaKnwnSigs <- mutSigData$blcaSIGS %>% 
  dplyr::select(starts_with("BLCA")) %>% 
  as.mutation.signatures()


msig1 <- matchSignatures(mutSign = vcf.sig, reference = cosmixSigs, 
                         threshold = 0.45, plot = TRUE) 
msig2 <- matchSignatures(mutSign = vcf.sig, reference = blcaKnwnSigs, 
                         threshold = 0.45, plot = TRUE)

msigPlot(vcf.sig, signature = 2, ylim = c(0, 0.20))

hm1 <- msig1$plot + ggtitle("Match to COSMIC signs.")
hm2 <- msig2$plot + ggtitle("Match to known BLCA signs.")

# Show
grid.arrange(hm1, hm2, ncol = 2)

### for Scott have the profiles + just raw variant counts  - be sure to also filter these variants
### by filters we normally use on the IonTorrent data


# vcfRate <- tmpVcf[, c("Chrom", "Position", "Ref", "Variant", "Allele.Call", "Filter", "Sample.Name", "Quality", "Frequency")]
# colnames(vcfRate) <- c("CHROM", "POS", "REF", "ALT", "XTR1", "FILTER", "SAMPLEID", "QUAL", "FREQ")
# vcfRate$XTR1 <- str_replace_all(vcfRate$XTR1, "Heterozygous", "0/1")
# vcfRate$XTR1 <- str_replace_all(vcfRate$XTR1, "Homozygous", "1/1")
# vcfRate$ID <- "."
# vcfRate$INFO <- "."
# vcfRate$FILTER <- "."
# vcfRate$FORMAT <- "GT:PL"
# 
# vcfRate2 <- vcfRate[, c("CHROM", "POS", "ID", "REF", "ALT",
#                        "QUAL", "FILTER", "INFO", "FORMAT",
#                        "XTR1", "SAMPLEID", "FREQ")]
# 
# vcfRate3 <- filterSNV(dataSet = vcfRate2,  seq_colNames = c("REF", "ALT"))
# 
# vcfRate4 <- attachContext(mutData = vcfRate3,
#                          chr_colName = "CHROM",
#                          start_colName = "POS",
#                          end_colName = "POS",
#                          nucl_contextN = 3,
#                          BSGenomeDb = mm10)
# 
# 
# vcfRate4 <- removeMismatchMut(mutData = vcfRate4,
#                              refMut_colName = "REF",
#                              context_colName = "context",
#                              refMut_format = "N")
# 
# vcfRate4 <- attachMutType(mutData = vcfRate4,
#                          ref_colName = "REF",
#                          var_colName = "ALT",
#                          context_colName = "context")
# 
# vcfRate4$mutTypeSingle <- str_remove(str_remove(vcfRate4$mutType, ".*\\["), "\\].*")
# 
# 
# mutTypeCount2 <- NULL
# for (i in unique(vcfRate4$SAMPLEID)) {
#   tmpDf <-  vcfRate4[which(vcfRate4$SAMPLEID == i),]
#   tmpDf2 <- data.frame("sample" = i, table(tmpDf$mutTypeSingle))
#   mutTypeCount2 <- rbind(mutTypeCount2, tmpDf2)
# }
# 
# mutTypeCount2$Var1
# 
# c <- ggplot(mutTypeCount2, aes(fill=Var1, y=Freq, x=sample)) +
#   geom_bar(position="stack", stat="identity") + ylab("# of mutations") +
#   ggtitle("exonic single base mutations types (0.156Mb)") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         plot.title = element_text(hjust = 0.5))
# 
# c

# ggplot(vcfRate4, aes(x = FREQ/100, color = mutTypeSingle)) + geom_density() + facet_wrap(vars(SAMPLEID)) +
#   ylim(c(0,100)) + ggtitle("Per sample density plots of mutation type allele frequency") + xlab("Frequency") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(vcfRate4, aes(x = FREQ/100, color = mutType)) + geom_density() + facet_wrap(vars(SAMPLEID)) +
#   ylim(c(0,100)) + ggtitle("Per sample density plots of mutation type allele frequency") + xlab("Frequency") +
#   theme(plot.title = element_text(hjust = 0.5))


# grid.arrange(c, d, ncol = 2)

msigPlot(vcf.sig, signature = 2, ylim = c(0, 0.20))
allSigs <- mutSignatures::getCosmicSignatures()

allSigsJune2022 <- read.table("/br_z1/kevin_storage/cosmic/mutSig/COSMIC_v3.3_SBS_GRCh38.txt",
                              header = TRUE, stringsAsFactors = FALSE)
allSigsJune2022Df <- allSigsJune2022[, 2:ncol(allSigsJune2022)]
rownames(allSigsJune2022Df) <- allSigsJune2022$Type
allSigsJune2022s3 <- as.mutation.signatures(allSigsJune2022Df)


allSigs2015 <- read.table("/br_z1/kevin_storage/cosmic/mutSig/COSMIC_v2_SBS_GRCh38.txt",
                              header = TRUE, stringsAsFactors = FALSE)
allSigs2015Df <- allSigs2015[, 2:ncol(allSigs2015)]
rownames(allSigs2015Df) <- allSigs2015$Type
allSigs2015s3 <- as.mutation.signatures(allSigs2015Df)


cosmixSigs <- mutSigData$blcaSIGS %>% 
  dplyr::select(starts_with("COSMIC")) %>% 
  as.mutation.signatures()

blcaKnwnSigs <- mutSigData$blcaSIGS %>% 
  dplyr::select(starts_with("BLCA")) %>% 
  as.mutation.signatures()


msigAll <- matchSignatures(mutSign = vcf.sig, reference = allSigsJune2022s3, 
                         threshold = 0.45, plot = TRUE)


msigAllOld <- matchSignatures(mutSign = vcf.sig, reference = allSigs2015s3, 
                           threshold = 0.45, plot = TRUE)

testOutDf <- msigAll$distanceDataFrame
testOutDf$Diff <- 1- testOutDf$dist
msigPlot(vcf.sig, signature = 2, ylim = c(0, 0.20))


testOutDf2 <- msigAllOld$distanceDataFrame
testOutDf2$Diff <- 1- testOutDf2$dist



### reduced for easy viewing
allSigsJune2022Df_red <- allSigsJune2022Df
abovePointSix <- c("30", "11", "23", "7b", "32")
abovePointSix <- paste0("SBS", abovePointSix)

allSigsJune2022Df_red <- allSigsJune2022Df_red[,which(colnames(allSigsJune2022Df_red) %in% abovePointSix)]
allSigsJune2022_reds3 <- as.mutation.signatures(allSigsJune2022Df_red)


matchSignatures(mutSign = vcf.sig, reference = allSigsJune2022_reds3, 
                threshold = 0.45, plot = TRUE)



### redoing the mutational signature analysis but for only concordant features

concordVcf <- read.table("/avatar_data3/eros_tmp/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562/plugin_out/variantCaller_out.1152/R_2022_07_01_11_06_33_user_AUS5-239-BBN_mouse_bladder_MG.xls",
                     sep = "\t", header = TRUE)
concordVcf$Sample2 <- str_remove(concordVcf$Sample.Name, "\\_.*")
concordVcf$string <- paste0(concordVcf$Sample2, concordVcf$Chrom, concordVcf$Position)
concordVcf <- concordVcf[which(concordVcf$string %in% combinedVars_concord_goodsamps$string), ]
# concordVcf <- concordVcf[which(concordVcf$string %in% combinedVars_concord_goodsamps_exon$string), ]

concordVcf2 <- concordVcf[, c("Chrom", "Position", "Ref", "Variant", "Allele.Call", "Filter", "Sample.Name", "Quality")]
colnames(concordVcf2) <- c("CHROM", "POS", "REF", "ALT", "XTR1", "FILTER", "SAMPLEID", "QUAL")

concordVcf2$XTR1 <- str_replace_all(concordVcf2$XTR1, "Heterozygous", "0/1")
concordVcf2$XTR1 <- str_replace_all(concordVcf2$XTR1, "Homozygous", "1/1")
concordVcf2$ID <- "."
concordVcf2$INFO <- "."
concordVcf2$FILTER <- "."
concordVcf2$FORMAT <- "GT:PL"

concordVcf3 <- concordVcf2[, c("CHROM", "POS", "ID", "REF", "ALT",
                       "QUAL", "FILTER", "INFO", "FORMAT", 
                       "XTR1", "SAMPLEID")]

concordVcf3 <- filterSNV(dataSet = concordVcf3,  seq_colNames = c("REF", "ALT"))

concordVcf4 <- attachContext(mutData = concordVcf3,
                         chr_colName = "CHROM",
                         start_colName = "POS",
                         end_colName = "POS",
                         nucl_contextN = 3,
                         BSGenomeDb = mm10)


concordVcf4 <- removeMismatchMut(mutData = concordVcf4,                  
                             refMut_colName = "REF",
                             context_colName = "context",
                             refMut_format = "N")

concordVcf4 <- attachMutType(mutData = concordVcf4,            
                         ref_colName = "REF",              
                         var_colName = "ALT",             
                         context_colName = "context") 

head(concordVcf4) %>% kable() %>% kable_styling(bootstrap_options = "striped")

concordVcfCounts <- countMutTypes(mutTable = concordVcf4,
                              mutType_colName = "mutType",
                              sample_colName = "SAMPLEID")


### used for alexandrov paper
sigprofile_df <- data.frame(concordVcfCounts@counts)
colnames(sigprofile_df) <- unlist(concordVcfCounts@sampleId)
sigprofile_df$`Mutation Types` <- unlist(concordVcfCounts@mutTypes)
sigprofile_df2 <- sigprofile_df[, c(colnames(sigprofile_df)[14], colnames(sigprofile_df)[1:13])]
# write.table(sigprofile_df2, "/br_z1/kevin_storage/misc/20221004concordMutsMatrix.txt", sep = "\t",
#             col.names = TRUE, row.names = FALSE, quote = FALSE)

which(apply(sigprofile_df2[, 2:14], 2, sum) < 10)
which(apply(sigprofile_df2[, 2:14], 2, function(x) sum(x)/0.156) < 13.8)
# sigprofile_df2_tmbMinMut <- sigprofile_df2[, -which(apply(sigprofile_df2[, 2:14], 2, sum) < 10)]
sigprofile_df2_tmbMinMut <- sigprofile_df2
sigprofile_df2_tmbMinMut_red <- sigprofile_df2_tmbMinMut[, -which(colnames(sigprofile_df2_tmbMinMut) %in% c("PP22-001_Dx81", "PP22-009_Dx89", "PP22-010_Dx90"))]

sigprofile_df2_tmbMinMut_final <- data.frame("Mutation Types" = sigprofile_df2$`Mutation Types`, sigprofile_df2_tmbMinMut)
sigprofile_df2_tmbMinMut_red_final <- data.frame("Mutation Types" = sigprofile_df2$`Mutation Types`, sigprofile_df2_tmbMinMut_red)

write.table(sigprofile_df2_tmbMinMut_final, "/br_z1/kevin_storage/misc/20221004concordMutMat_AllMut.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(sigprofile_df2_tmbMinMut_red_final, "/br_z1/kevin_storage/misc/20221004concordMutMat_AllMut_red.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)


# write.table(sigprofile_df2_tmbMinMut_final, "/br_z1/kevin_storage/misc/20221004concordMutMat_TmbMinMut.txt", sep = "\t",
#             col.names = TRUE, row.names = FALSE, quote = FALSE)
# write.table(sigprofile_df2_tmbMinMut_red_final, "/br_z1/kevin_storage/misc/20221004concordMutMat_TmbMinMut_red.txt", sep = "\t",
#             col.names = TRUE, row.names = FALSE, quote = FALSE)

# vcf.params2 <- 
#   mutSignatures::setMutClusterParams( 
#     num_processesToExtract = 5,
#     num_totIterations = 1000,
#     num_parallelCores = 28) 


# 
# vcf.analysis.concord <- 
#   decipherMutationalProcesses(input = concordVcfCounts,
#                               params = vcf.params2)


# Retrieve signatures (results)
vcf.sig.concord <- vcf.analysis.concord$Results$signatures

# Retrieve exposures (results)
vcf.exp.concord <- vcf.analysis.concord$Results$exposures


msigAllConcord <- matchSignatures(mutSign = vcf.sig.concord, reference = allSigsJune2022s3, 
                           threshold = 0.45, plot = TRUE)


testOutDfConcord <- msigAllConcord$distanceDataFrame
testOutDfConcord$Diff <- 1 - testOutDfConcord$dist

# 
# msigPlot(vcf.sig, signature = 1, ylim = c(0, 0.50))
# msigPlot(vcf.sig.concord, signature = 4, ylim = c(0, 0.5))


reducedSignConcord <- paste0("SBS", c(16, 32, "10b", 11, 30))
allSigsJune2022Df_concord <- allSigsJune2022[, which(colnames(allSigsJune2022) %in% reducedSignConcord)]
rownames(allSigsJune2022Df_concord) <- allSigsJune2022$Type
allSigsJune2022s3_concord <- as.mutation.signatures(allSigsJune2022Df_concord)

matchSignatures(mutSign = vcf.sig.concord, reference = allSigsJune2022s3_concord, 
                threshold = 0.45, plot = TRUE)

msigPlot(vcf.sig.concord, signature = 1, ylim = c(0, 0.5))

# save.image("/br_z1/kevin_storage/misc/20220920mutSig.RData")
load("/br_z1/kevin_storage/misc/20220920mutSig.RData")




### from previous count data, combine to see of noise vs less noise for patterns - also get var data to see certain VAF dist for groups
### then can also combine to see low vs high tmb

noise_samps <- c("PP22-007", "PP22-009", "PP22-012", "PP22-013")


vcfRateConcord$noiseGroup <- "no"
vcfRateDiscord$noiseGroup <- "yes"

allVarsLabeled <- rbind(vcfRateConcord, vcfRateDiscord)
allVarsLabeled$SAMPLEID2 <- str_remove(allVarsLabeled$SAMPLEID, "\\_.*")

allVarsLabeled_onlyC <- allVarsLabeled[which(allVarsLabeled$mutTypeSingle == "C>T"),]

ggplot(allVarsLabeled, aes(y=FREQ, x= mutTypeSingle)) +
  geom_boxplot() + facet_wrap(vars(SAMPLEID2))

ggplot(allVarsLabeled_onlyC, aes(y=FREQ, x= mutType)) +
  geom_boxplot() + geom_jitter() + 
  facet_wrap(vars(noiseGroup), nrow = 2)

### below only looks at sampels with FFPE noise since, there are real vars within each set. still no enough muts to tell ...
###

allVarsLabeled_onlyC_noise <- allVarsLabeled_onlyC[which(allVarsLabeled_onlyC$SAMPLEID2 %in% noise_samps), ]

ggplot(allVarsLabeled_onlyC_noise, aes(y=FREQ, x= mutType)) +
  geom_boxplot() + geom_jitter() + 
  facet_wrap(vars(noiseGroup), nrow = 2)

### looking at concord exon-put mutations
combinedVars_reseq_goodsamps_exon_put$Sample2 <- str_remove(combinedVars_reseq_goodsamps_exon_put$Sample, "\\_.*")
combinedVars_reseq_goodsamps_exon_put$string <- paste0(combinedVars_reseq_goodsamps_exon_put$Sample2,
                                                       combinedVars_reseq_goodsamps_exon_put$Chr,
                                                       combinedVars_reseq_goodsamps_exon_put$Start,
                                                       combinedVars_reseq_goodsamps_exon_put$End)

combinedVars_goodsamps_exon_put$Sample2 <- str_remove(combinedVars_goodsamps_exon_put$Sample, "\\_.*")
combinedVars_goodsamps_exon_put$string <- paste0(combinedVars_goodsamps_exon_put$Sample2,
                                                       combinedVars_goodsamps_exon_put$Chr,
                                                       combinedVars_goodsamps_exon_put$Start,
                                                       combinedVars_goodsamps_exon_put$End)

combinedVars_concord_exon <- combinedVars_goodsamps_exon_put[which(combinedVars_goodsamps_exon_put$string %in%
                                                                     combinedVars_reseq_goodsamps_exon_put$string), ]


### getting optimal k-proceeses

# prelimProcessAssess (this may take a while)
# https://github.com/dami82/mutSignatures/issues/15
test <- prelimProcessAssess(input = concordVcfCounts, maxProcess = 14, plot = FALSE)

# Build a plot
ggplot(test, aes(x=numProcess, y=percentErr)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(1, 14, by = 1)) +
  theme_bw() + xlab('number of signatures') + ylab('% Error')


# increase the number of signatures, k=5, k=6, ... at each cycle

my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 2,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = concordVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
a <- match_my_sigs_ki$plot

b <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

# msigPlot(my_sigs_ki, signature = 1, ylim = c(0, 0.15))
# msigPlot(my_sigs_ki, signature = 2, ylim = c(0, 0.15))
# msigPlot(my_sigs_ki, signature = 3, ylim = c(0, 0.15))


gridExtra::grid.arrange(a, b)


mean8 <- mean(my_silhouette_ki$silhouette_value)

### k = 7 0.647;  k = 9 0.637; k=10 ;


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 7,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = concordVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
c <- match_my_sigs_ki$plot

d <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean7 <- mean(my_silhouette_ki$silhouette_value)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 6,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = concordVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
e <- match_my_sigs_ki$plot

f <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean6 <- mean(my_silhouette_ki$silhouette_value)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 5,
    num_totIterations = 50,
    num_parallelCores = 20) 
tmp_analysis_ki <- decipherMutationalProcesses(input = concordVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
g <- match_my_sigs_ki$plot

h <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean5 <- mean(my_silhouette_ki$silhouette_value)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 4,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = concordVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
i <- match_my_sigs_ki$plot

j <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean4 <- mean(my_silhouette_ki$silhouette_value)


grid.arrange(b,a, d, c, f, e, h, g, j, i, ncol = 2)

# means 8 to 4; 0.63, 0.645, 0.664, 0.67, 0.61 respectively
# winner is k = 5



### running for original sequencing run 1

test <- prelimProcessAssess(input = allVcfCounts, maxProcess = 14, plot = FALSE)

# Build a plot
ggplot(test, aes(x=numProcess, y=percentErr)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(1, 14, by = 1)) +
  theme_bw() + xlab('number of signatures') + ylab('% Error')



my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 8,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
a <- match_my_sigs_ki$plot

b <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 


gridExtra::grid.arrange(a, b)


mean8 <- mean(my_silhouette_ki$silhouette_value)

### k = 7 0.647;  k = 9 0.637; k=10 ;


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 7,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
c <- match_my_sigs_ki$plot

d <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean7 <- mean(my_silhouette_ki$silhouette_value)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 6,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
e <- match_my_sigs_ki$plot

f <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean6 <- mean(my_silhouette_ki$silhouette_value)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 5,
    num_totIterations = 50,
    num_parallelCores = 20) 
tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
g <- match_my_sigs_ki$plot

h <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean5 <- mean(my_silhouette_ki$silhouette_value)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 4,
    num_totIterations = 50,
    num_parallelCores = 20) 

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
i <- match_my_sigs_ki$plot

j <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean4 <- mean(my_silhouette_ki$silhouette_value)


grid.arrange(b,a, d, c, f, e, h, g, j, i, ncol = 2)


my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 9,
    num_totIterations = 50,
    num_parallelCores = 20)

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
k <- match_my_sigs_ki$plot

l <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean9 <- mean(my_silhouette_ki$silhouette_value)



my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 3,
    num_totIterations = 50,
    num_parallelCores = 20)

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
m <- match_my_sigs_ki$plot

n <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean3 <- mean(my_silhouette_ki$silhouette_value)

my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 2,
    num_totIterations = 50,
    num_parallelCores = 20)

tmp_analysis_ki <- decipherMutationalProcesses(input = allVcfCounts,
                                               params = my_params_ki)

my_sigs_ki <- tmp_analysis_ki$Results$signatures
my_silhouette_ki <- tmp_analysis_ki$Supplementary$silhouette
match_my_sigs_ki <- matchSignatures(my_sigs_ki, my_sigs_ki)
o <- match_my_sigs_ki$plot

p <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 

mean2 <- mean(my_silhouette_ki$silhouette_value)

### means for seq1 are 9, 8, 7, 6, 5, 4, 3; 0.64, 0.73, 0.71, 0.66 - winner 7; although lower, everything but last basis is strong

grid.arrange(l, k, b,a, d, c, f, e, h, g, ncol = 2)


### after doing the testing for different groups k
# save.image("/br_z1/kevin_storage/misc/20220925mutSig.RData")

load("/br_z1/kevin_storage/misc/20220925mutSig.RData")
dev.off()
par(mfrow=c(2,2))
msigPlot(vcf.sig.concord, signature = 2, ylim = c(0, 1))
msigPlot(allSigsJune2022s3, signature = 22, ylim = c(0, 1))
msigPlot(allSigsJune2022s3, signature = 76, ylim = c(0, 1))


dev.off()
par(mfrow=c(2,2))
msigPlot(vcf.sig.concord, signature = 5, ylim = c(0, 1))
msigPlot(allSigsJune2022s3, signature = 14, ylim = c(0, 1))
msigPlot(allSigsJune2022s3, signature = 37, ylim = c(0, 1))
msigPlot(allSigsJune2022s3, signature = 39, ylim = c(0, 1))





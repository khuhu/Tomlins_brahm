library(mutSignatures)
library(vcfR)
library(stringr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(gridExtra)

mm10 <- BSgenome.Mmusculus.UCSC.mm10

variantTable <- readxl::read_xlsx("/br_z1/kevin_storage/misc/20220927FinalDNAseqvariantdata.xlsx")
variantTable2 <- variantTable[-which(is.na(variantTable$Sample)), ]
variantTable2$string <- paste0(variantTable2$Sample, variantTable2$Chr, variantTable2$Start)

tmpVcf <- read.table("/avatar_data3/eros_tmp/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562/plugin_out/variantCaller_out.1152/R_2022_07_01_11_06_33_user_AUS5-239-BBN_mouse_bladder_MG.xls",
                     sep = "\t", header = TRUE)
tmpVcf$Sample.Name <- str_remove(tmpVcf$Sample.Name, "\\_.*")
tmpVcf$string <- paste0(tmpVcf$Sample.Name, tmpVcf$Chrom, tmpVcf$Position)
tmpVcf <- tmpVcf[which(tmpVcf$string %in% variantTable2$string), ]


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


allVcfCounts <- countMutTypes(mutTable = tmpVcf4,
                              mutType_colName = "mutType",
                              sample_colName = "SAMPLEID")



sigprofile_df <- data.frame(allVcfCounts@counts)
colnames(sigprofile_df) <- unlist(allVcfCounts@sampleId)
sigprofile_df$`Mutation Types` <- unlist(allVcfCounts@mutTypes)
sigprofile_df2 <- sigprofile_df[, c(colnames(sigprofile_df)[11], colnames(sigprofile_df)[1:10])]

which(apply(sigprofile_df2[, 2:11], 2, sum) < 10)
which(apply(sigprofile_df2[, 2:11], 2, function(x) sum(x)/0.156) < 13.8)

write.table(sigprofile_df2, "/br_z1/kevin_storage/misc/20221004concordMutMat_aaronList.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

### finding optimal k
###
###

test <- prelimProcessAssess(input = allVcfCounts, maxProcess = 14, plot = FALSE)

# Build a plot
ggplot(test, aes(x=numProcess, y=percentErr)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(1, 14, by = 1)) +
  theme_bw() + xlab('number of signatures') + ylab('% Error')



my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 11,
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

mean11 <- mean(my_silhouette_ki$silhouette_value)

my_params_ki <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 10,
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



mean10 <- mean(my_silhouette_ki$silhouette_value)


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
e <- match_my_sigs_ki$plot

f <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 



mean9 <- mean(my_silhouette_ki$silhouette_value)


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
g <- match_my_sigs_ki$plot

h <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 



mean8 <- mean(my_silhouette_ki$silhouette_value)



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
i <- match_my_sigs_ki$plot

j <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
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
k <- match_my_sigs_ki$plot

l <- ggplot(my_silhouette_ki, aes(y = silhouette_value, x = id, fill = signature)) +
  geom_bar(stat = 'identity', width = 1) + ylim(c(-1, 1)) + 
  coord_flip() + theme_bw() 



mean6 <- mean(my_silhouette_ki$silhouette_value)

grid.arrange(d, c, f, e, h, g, ncol = 2)



### running nmf after getting best k = 9
###
###



vcf.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 9,
    num_totIterations = 1000,
    num_parallelCores = 28) 



vcf.analysis <- 
  decipherMutationalProcesses(input = allVcfCounts,
                              params = vcf.params)


# Retrieve signatures (results)
vcf.sig <- vcf.analysis$Results$signatures

# Retrieve exposures (results)
vcf.exp <- vcf.analysis$Results$exposures

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

msigAll <- matchSignatures(mutSign = vcf.sig, reference = allSigsJune2022s3, 
                           threshold = 0.45, plot = TRUE)

testOutDf <- msigAll$distanceDataFrame
testOutDf$Diff <- 1- testOutDf$dist


msigPlot(vcf.sig, signature = 1, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 2, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 3, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 4, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 5, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 6, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 7, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 8, ylim = c(0, 0.20))
msigPlot(vcf.sig, signature = 9, ylim = c(0, 0.20))



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

# save.image("/br_z1/kevin_storage/misc/20220928aaronSetBbn.Rdata")

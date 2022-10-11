#install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", version = "3.8")

source("/home/kevhu/scripts/fastReadFile.R")

library(stringr)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(AnnotationHub)
ah <- AnnotationHub()
edb <- ah[["AH64944"]]


###
###
### first may try to aggregate all the different genes and see the concordance between the lists

ccpOncominev3 <- as.character(unlist(read.table("/home/kevhu/data/oncomineV3")))
mskGeneList <- read_xlsx("/home/kevhu/data/mskImpactGeneListSuppTab4.xlsx", skip = 2)
mskGeneList2 <- unique(mskGeneList$Gene)

mskHotspots <- read.table("/home/kevhu/data/hotspots_msk.txt",sep = "\t", header = TRUE, stringsAsFactors = FALSE)

tumorSupressors <- c("TP53","RB1", "WT1", "NF1", "NF2", "APC", "TSC1", "TSC2", "SMAD4",
                     "BRCA1","BRCA2", "PTEN", "STK11", "MSH2", "MLH1", "CDH1", 
                     "VHL", "CDKN2A", "PTCH1", "MEN1")

tumorSupressors[which(tumorSupressors %in% mskGeneList2)]

combinedGenesList <- unique(c(mskGeneList2, ccpOncominev3))

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  #humanx <- unique(genesV2[, 2])
  #return(humanx)
  return(genesV2)
}

mouseCancerGenes <- convertHumanGeneList(combinedGenesList)
### FAM123B is Amer1, PAK7 is p21, MLL2 is KMT2D, MLL is Kmt2a, MLL3 is Kmt2c, PARK2 is prkn,
### FAM175A is Abraxas1, RFWD2 is Cop1, FAM46C is Tent5c, MYCL1 is Mcyl

googledConversions <- c("H3f3a", "Cdkn2a", "Xiap", "Amer1", "P21", "Bbc3", "Kmt2d", 
                        "Kmt2a", "Mre11a", "Gnaq", "Kmt2c", "Prkn", "Notch4", "Hist1h3b",
                        "Npm1", "Abraxas1", "Cop1", "Tent5c", "Mycl")

tmp1 <- combinedGenesList[-which(combinedGenesList %in% mouseCancerGenes$HGNC.symbol)]
tmp1 <- tmp1[-which(tmp1 == "-")]

### my googled list is in reverse to what tmp1 is
tmp2 <- cbind(tmp1, rev(googledConversions))
colnames(tmp2) <- colnames(mouseCancerGenes)
mouseCancerGenesV2 <- rbind(mouseCancerGenes, tmp2)

mouseCancerGenesV3 <- mouseCancerGenesV2[!duplicated(mouseCancerGenesV2$HGNC.symbol),]
googleSet2 <- mouseCancerGenesV2[duplicated(mouseCancerGenesV2$HGNC.symbol),]
mouseCancerGenesV3$MGI.symbol[which(mouseCancerGenesV3$HGNC.symbol == "MYC")] <- "Myc"
mouseCancerGenesV3$MGI.symbol[which(mouseCancerGenesV3$HGNC.symbol == "SPOP")] <- "Spop"
mouseCancerGenesV3$MGI.symbol[which(mouseCancerGenesV3$HGNC.symbol == "EIF1AX")] <- "Eif1ax"
mouseCancerGenesV3$MGI.symbol[which(mouseCancerGenesV3$HGNC.symbol == "MYCN")] <- "Mycn"
mouseCancerGenesV3$MGI.symbol[which(mouseCancerGenesV3$HGNC.symbol == "PMS2")] <-"Pms2"



### only keeping hotspots found in COSMIC.need to get coordinates for genomic hotspots too:
### after limiting the list. all hotspots are within already targeted genes

cosmicTable <- faster.readfile(x = "/home/kevhu/data/CosmicMutantExport.tsv", y= 0,sep = '\t')
cosmicTable2 <- cosmicTable[-which(cosmicTable$Mutation_genome_position == ""),]
cosmicTableAAsubset <- cosmicTable2$Mutation_AA
cosmicTableAAsubset <- str_remove_all(cosmicTableAAsubset, "p\\.*")
cosmicTableAAsubset <- gsub(".{1}$", "", cosmicTableAAsubset)
tableToMatch <- paste0(cosmicTable2$Gene_name, cosmicTableAAsubset)
keysToMatch <- paste0(mskHotspots$Gene, mskHotspots$Residue)


cosmicTable3 <- cosmicTable2[match(keysToMatch, tableToMatch),]

mskHotspots$GenomicLocation <- cosmicTable3$Mutation_genome_position
mskHotspots$MutStrand <- cosmicTable3$Mutation_strand_SNP
mskHotspots$cosmicID <- cosmicTable3$Mutation_ID
mskHotspots.Cosmic <- mskHotspots[-which(is.na(cosmicTable3$Gene_name)),]

#write.table(mskHotspots.Cosmic,"/home/kevhu/data/20190201mouseHotspots.txt",
#            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


### this shows all hotspot mutations found on gene panel
hotspotOnly <- unique(mskHotspots.Cosmic$Gene)[-which(unique(mskHotspots.Cosmic$Gene) %in% unique(c(mouseCancerGenesV3$HGNC.symbol, toupper(mouseCancerGenesV3$MGI.symbol))))]


### adding annotations for the targets

targetAnno <- rep("Oncogene", nrow(mouseCancerGenesV3))
targetAnno[which(mouseCancerGenesV3$HGNC.symbol %in% tumorSupressors)] <- "TumorSupressorGene"
mouseCancerGenesV3$Annotation <- targetAnno

#which(mouseCancerGenesV3$HGNC.symbol %in% unique(mskHotspots.Cosmic$Gene))
annotatedAmplicon <- rep("CopyNumber", nrow(mouseCancerGenesV3))
annotatedAmplicon[which(mouseCancerGenesV3$HGNC.symbol %in% unique(mskHotspots.Cosmic$Gene))] <- "CopyNumber/Hotspot"
annotatedAmplicon[which(toupper(mouseCancerGenesV3$MGI.symbol) %in% unique(mskHotspots.Cosmic$Gene))] <- "CopyNumber/Hotspot"
annotatedAmplicon[which(mouseCancerGenesV3$HGNC.symbol %in% tumorSupressors)] <- "AllExons"
mouseCancerGenesV3$ampliconPurpose <- annotatedAmplicon



#write.table(mouseCancerGenesV3,"/home/kevhu/data/20190129preliminaryMousePanelListCosmic.txt", sep = "\t", col.names = TRUE, row.names = FALSE,
#            quote = FALSE)


####: 20190204 reading in BLAST results

uniProtCanonical <- read.table("/home/kevhu/data/uniprotCanonicalList.txt", skip = 1, header = FALSE,
                               stringsAsFactors = FALSE, sep = "\t")

blast.res <- read.table("/home/kevhu/data/20190204blastMouse.txt", header = FALSE,
                        stringsAsFactors = FALSE, sep = "\t")
colnames(blast.res) <- c("queryID", "subjectID", "length", "queryLength", "subjectLength",
                         "queryStart", "queryEnd", "subjectStart", "subjectEnd", "percentIdentical",
                         "percentPosScoreMatches", "numberOfGaps","e-value")

### only really interested in genes where we are looking at hotspots ..
tmpTable <- mouseCancerGenesV3
tmpTable$uniProt <- uniProtCanonical$V2[match(mouseCancerGenesV3$HGNC.symbol, uniProtCanonical$V1)]

tmpTable$HGNC.symbol[which(is.na(tmpTable$uniProt))]
tmpTable$uniProt[22] <- "P12524"
tmpTable$uniProt[305] <- "P46100"
tmpTable$uniProt[351] <- "Q6UWZ7"
tmpTable <- tmpTable[-348,]

tmpTable$uniProt[which(tmpTable$ampliconPurpose == "CopyNumber/Hotspot")]

matchIdxTmp <- NULL
for(i in tmpTable$uniProt[which(tmpTable$ampliconPurpose == "CopyNumber/Hotspot")]){
  print(i)
  a <- which(grepl(i, blast.res$queryID))
  matchIdxTmp <- c(matchIdxTmp, a)
}

blast.res2 <- blast.res[matchIdxTmp,]

boxplot(blast.res2$percentIdentical)


### 20190208 : after blast results
### finding only regions where we aren't tiling a TSG because at that point we don't need a blast match
### also editing total list after Anj and Eric's comments

mouseCancerGenesV4 <- mouseCancerGenesV3
mouseCancerGenesV4$ampliconPurpose[which(mouseCancerGenesV4$ampliconPurpose == 'AllExons')] <- "CDS"

oncoToTsg <- c("Smad2", "Tgfbr2", "Arid1a", "Sox9", "Chek2", "Pms1", "Pms2", "Notch1", "Pole", "Axin1")

mouseCancerGenesV4$ampliconPurpose[which(mouseCancerGenesV4$MGI.symbol %in% oncoToTsg)] <- "CDS"
mouseCancerGenesV4$Annotation[which(mouseCancerGenesV4$MGI.symbol %in% oncoToTsg)] <- "TumorSupressorGene"

### adding other tsgs not on mouse list
tsg.hgnc <- c("DCC", "CYLD", "PTCH2")
tsg.mgi <- c("Dcc", "Cyld", "Ptch2")
tsg.purpose <- rep("CDS",3)
tsg.anno <- rep("TumorSurpressorGene", 3)
tsg.tmp <- cbind(tsg.hgnc, tsg.mgi, tsg.anno, tsg.purpose)
colnames(tsg.tmp) <- colnames(mouseCancerGenesV4)
mouseCancerGenesV4 <- rbind(mouseCancerGenesV4, tsg.tmp)

### same for oncogenes
mouseCancerGenesV4 <- rbind(mouseCancerGenesV4, c("GLI1", "Gli1", "TumorSupressorGene", "CDS"))


### now I filter out the cosmic list by only oncogenes -> helps reduce amount of things not found from blast results 
### copying code from above so I can just blast using the smaller list

mouseHotspotGenes <- mouseCancerGenesV4$HGNC.symbol[which(mouseCancerGenesV4$ampliconPurpose == "CopyNumber/Hotspot")]
mskHotspots.Cosmic.onlyOncogenes <- mskHotspots.Cosmic[which(mskHotspots.Cosmic$Gene %in% mouseHotspotGenes),]


### above is mostly not working ... instead maybe just try to pull from primary amino acid structure
### so below needs to be tweaked just a little bit with tryCatch in order to NULL things out of bounds
library(Biostrings)

uniprot.Mouse <- readAAStringSet("/home/kevhu/data/20190211mouseUniprotOncogeneCanonical.fa")


blastList.human <- NULL
for(i in 1:nrow(sortedVars)){
  tmpAAseq <- which(grepl(sortedVars$Gene[i],toupper(names(uniprot.Mouse))))
  for (k in seq_along(tmpAAseq)) {
    for (j in 1:3) {
      if (j == 1) {
        tmpAA <- tryCatch(subseq(uniprot.Mouse[[tmpAAseq[k]]], start = (sortedVars$position[i] - 6), end = sortedVars$position[i]), error = function(x) return(NULL))
        if (is.null(tmpAA)) {
          next()
        }
        tmpAA <- as.character(tmpAA)
        blastList.human <- c(blastList.human, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList.human <- c(blastList.human, tmpAA)
      }
      if (j == 2) {
        tmpAA <- tryCatch(subseq(uniprot.Mouse[[tmpAAseq[k]]], start = sortedVars$position[i], end = (sortedVars$position[i] + 6)), error = function(x) return(NULL))
        if (is.null(tmpAA)) {
          next()
        }
        tmpAA <- as.character(tmpAA)
        blastList.human <- c(blastList.human, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList.human <- c(blastList.human, tmpAA)
      }
      if (j == 3) {
        tmpAA <- tryCatch(subseq(uniprot.Mouse[[tmpAAseq[k]]], start = (sortedVars$position[i] - 3), end = (sortedVars$position[i] + 3)), error = function(x) return(NULL))
        if (is.null(tmpAA)) {
          next()
        }
        tmpAA <- as.character(tmpAA)
        blastList.human <- c(blastList.human, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList.human <- c(blastList.human, tmpAA)
      }
    }
  }
}

blastList.human
con.human <-file("/home/kevhu/data/20190213customAaTblastnHuman.fa")
writeLines(blastList.human, con.human)
close(con.human)

### alright going back to the original idea .... because it looks way easier .....  I'm real dumn if this works ....
### below checked .. and it seems to grab the genomic positions
### most interesting part is I may need to grab unreviewed proteins because the ensdb identifies both
### some of the unreviewed ones but not reviewed ones i.e NRAS


library(GenomicRanges)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)

uniProtCanonical <- read.table("/home/kevhu/data/20190214wholeCanonicalMouseProt.txt", skip = 1, header = FALSE,
                               stringsAsFactors = FALSE, sep = "\t")
load("/home/kevhu/data/20190214mskHotspots.onlyOnco.df.Robj")
mouseCancerGenesV4 <- read.table("/home/kevhu/data/20190214mouseCancerGenesV4.txt", header = TRUE,
                                 stringsAsFactors = FALSE, sep = "\t")

#canProtName <- uniProtCanonical$V2[match(mskHotspots.Cosmic.onlyOncogenes$Gene,uniProtCanonical$V1)]


matchedVars <- mskHotspots.Cosmic.onlyOncogenes$Residue[match(uniProtCanonical$V1, mskHotspots.Cosmic.onlyOncogenes$Gene)]
uniProtCanonical$residueLocation <- matchedVars


#### below is example of it working 
#### 20190630 KH: going to make a working example of a better version. before the script only looked at the Nth amino acid
#### now I will pull the 3/5 AA sequence in human as noted before and match it with the mouse - if there is a match
#### I can grab the coordinates

testNRAS <- IRanges(start = 59, end = 61, names = "Q9D091")
testNRAS.out <- proteinToGenome(testNRAS, EnsDb.Mmusculus.v79, idType = "uniprot_id")
as.character(testNRAS.out[[1]]@ranges)
as.character(testNRAS.out[[1]]@seqnames@values)

testNRAS2 <- GRanges(seqnames = testNRAS.out[[1]]@seqnames@values, ranges = testNRAS.out[[1]]@ranges)
testNRAS.out2 <- genomeT(testNRAS2, EnsDb.Mmusculus.v79)
NRAS <- proteins(EnsDb.Mmusculus.v79, filter = ~ genename == "Nras")
NRAS.human <- proteins(EnsDb.Hsapiens.v86, filter = ~ genename == "NRAS", order.type = "descending")
NRAS.human$protein_sequence[1]
AAtrio <- substr(NRAS.human$protein_sequence[1], 59, 63)
testRes <- matchPattern(AAtrio, NRAS$protein_sequence[1])
### to get exact location of the match
(testRes@ranges@width - 1 + 2 * testRes@ranges@start)/2

### after this create the range and search against 

justResPosition <- as.numeric(substr(uniProtCanonical$residueLocation,
                                     2, nchar(mskHotspots.Cosmic.onlyOncogenes$Residue)))

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


###
### so I couldn't figure out which ones weren't getting entered but they'd be "NA" anyways.
### so create empty vector of NA's first to get right final length


chromConvert2 <- rep("NA", nrow(uniProtCanonical))
seqConvert2 <- rep("NA", nrow(uniProtCanonical))


for(i in seq_along(uniProtCanonical$V2)) {
  tmpChrom <- NULL
  tmpSeq <- NULL
  print(paste("start:", i))
  tmpHumanProt <- proteins(EnsDb.Hsapiens.v86, filter = ~ genename == uniProtCanonical$V1[i], return.type="AAStringSet")
  longestHumProt <- which(tmpHumanProt@ranges@width == max(tmpHumanProt@ranges@width))[1]
  tmpHumanProt2 <- tmpHumanProt[[longestHumProt]]
  tmpAAtrio <- tryCatch(substr(tmpHumanProt2, justResPosition[i] - 2, justResPosition[i] + 2),
                        error = function(x) return(NULL))
  
  ###same thing with mouse
  tmpMouseProt <- proteins(EnsDb.Mmusculus.v79, filter = ~ genename == simpleCap(tolower(uniProtCanonical$V1[i])), return.type="AAStringSet")
  longestMouseProt <- which(tmpMouseProt@ranges@width == max(tmpMouseProt@ranges@width))[1]
  tmpMouseProt2 <- tmpMouseProt[[longestMouseProt]]
  
  if(is.null(tmpAAtrio) || length(tmpMouseProt) == 0){
    tmpConv <- tryCatch(proteinToGenome(IRanges(start = justResPosition[i], end = justResPosition[i], names = uniProtCanonical$V2[i]),
                                        EnsDb.Mmusculus.v79, idType = "uniprot_id"), warning = function(x) return(NULL))
    if(is.null(tmpConv)) {
      tmpChrom <- "NA"
      tmpSeq <- "NA"
      chromConvert2[i] <- tmpChrom
      seqConvert2[i] <- tmpSeq
      print(paste("end:", i))
      next()
    }
    if(length(tmpConv[[1]]) > 1) {
      tmpChrom <- tryCatch(as.character(tmpConv[[1]]@unlistData@seqnames@values[1]),
                           error = function(x) return(NULL))
      tmpSeq <- tryCatch(as.character(tmpConv[[1]]@unlistData@ranges[1]),
                         error = function(x) return(NULL))
      if(is.null(tmpChrom)) {
        tmpChrom <- as.character(tmpConv[[1]]@seqnames@values[1])
        tmpSeq <- as.character(tmpConv[[1]]@ranges[1])
      }
      else{
        tmpChrom <- as.character(tmpConv[[1]]@unlistData@seqnames@values[1])
        tmpSeq <- as.character(tmpConv[[1]]@unlistData@ranges[1])
      }
    }
    else{
      tmpChrom <- as.character(tmpConv[[1]]@seqnames@values)
      tmpSeq <- as.character(tmpConv[[1]]@ranges)
    }
    chromConvert2[i] <- tmpChrom
    seqConvert2[i] <- tmpSeq
    print(paste("end:", i))
    next()
  }

  #tmpConv <- tryCatch(proteinToGenome(tmpRange, EnsDb.Mmusculus.v79, idType = "uniprot_id"),
  #                    error = function(x) return(NULL))
  tmpQuery <- tryCatch(matchPattern(tmpAAtrio, tmpMouseProt2),
                      error = function(x) return(NULL))
  tmpLoc <- tryCatch(tmpQuery@ranges@width[1] - 1 + 2 * tmpQuery@ranges@start[1]/2,
                     error = function(x) return(NULL))
  #if(length(tmpLoc) == 0 | is.na(tmpLoc) | length(tmpLoc) == 0 | is.null(tmpLoc)) {
  if(is.null(tmpQuery) || is.null(tmpLoc) || is.na(tmpLoc)) {
    tmpChrom <- "NA"
    tmpSeq <- "NA"
    chromConvert2[i] <- tmpChrom
    seqConvert2[i] <- tmpSeq
    print(paste("end:", i))
    next()
  }
  print(tmpQuery)
  ### to get exact location of the match
  tmpConv <- tryCatch(proteinToGenome(IRanges(start = tmpLoc, end = tmpLoc, names = uniProtCanonical$V2[i]),
                                      EnsDb.Mmusculus.v79, idType = "uniprot_id"), warning = function(x) return(NULL))
  if(is.null(tmpConv)) {
    tmpChrom <- "NA"
    tmpSeq <- "NA"
    chromConvert2[i] <- tmpChrom
    seqConvert2[i] <- tmpSeq
    print(paste("end:", i))
    next()
  }
  if(length(tmpConv[[1]]) > 1) {
    tmpChrom <- tryCatch(as.character(tmpConv[[1]]@unlistData@seqnames@values[1]),
                       error = function(x) return(NULL))
    tmpSeq <- tryCatch(as.character(tmpConv[[1]]@unlistData@ranges[1]),
                       error = function(x) return(NULL))
    if(is.null(tmpChrom)) {
      tmpChrom <- as.character(tmpConv[[1]]@seqnames@values[1])
      tmpSeq <- as.character(tmpConv[[1]]@ranges[1])
    }
    else{
      tmpChrom <- as.character(tmpConv[[1]]@unlistData@seqnames@values[1])
      tmpSeq <- as.character(tmpConv[[1]]@unlistData@ranges[1])
    }
  }
  else{
    tmpChrom <- as.character(tmpConv[[1]]@seqnames@values)
    tmpSeq <- as.character(tmpConv[[1]]@ranges)
  }
 chromConvert2[i] <- tmpChrom
 seqConvert2[i] <- tmpSeq
 print(paste("end:", i))
}


uniProtCanonical$locationChr <- chromConvert2
uniProtCanonical$locationRange <- seqConvert2

uniProtCanonical.reduced <- uniProtCanonical[-which(uniProtCanonical$locationChr == "NA"),]
uniProtCanonical.reduced <- uniProtCanonical.reduced[!duplicated(uniProtCanonical.reduced$locationRange),]

tmpLocation <- strsplit(uniProtCanonical.reduced$locationRange, "-")
oncoHotspotBedF1 <- paste0("chr",uniProtCanonical.reduced$locationChr)
oncoHotspotBedF2 <- unlist(lapply(tmpLocation, '[[',1))
#oncoHotspotBedF3 <- unlist(lapply(tmpLocation, '[[',2))
oncoHotspotBedF4 <- paste0(uniProtCanonical.reduced$V1, "|",uniProtCanonical.reduced$residueLocation)

mouseHotspotBed <- data.frame("V1" = oncoHotspotBedF1, "V2" = oncoHotspotBedF2, "V3" = oncoHotspotBedF2,
                              "V4" = oncoHotspotBedF4, stringsAsFactors = FALSE)

#mouseHotspotBed3Col <- data.frame("Chr" = oncoHotspotBedF1,"Start" = oncoHotspotBedF2, stringsAsFactors = FALSE)


#write.table(mouseHotspotBed3Col, "/home/kevhu/data/20190630mouseOncoGeneHotspots3Col.bed",
#            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)


### need to combine these scopes with other ones
###
###


snpPanelTmp <- read.table("/home/kevhu/data/20190628snpCompPanel.bed", sep = "\t", header = FALSE)

allPanelSnps  <- rbind(snpPanelTmp, mouseHotspotBed)

#write.table(allPanelSnps, "/home/kevhu/data/20190701allPanelSnps.bed", quote = FALSE, row.names = FALSE, col.names = FALSE,
#            sep = "\t")

### next is figuring out how many amplicons of the hotspots cover ... also seems like there is some overlap of the hotspots ... mainly just fgfr2 ()
###
###

table(uniProtCanonical.reduced$V1)

### pseudo code for getting the 10 amplicons that tile the gene
### (1) grab genic/exonix sequence. (2) linearize if exonic (3) grab 10 individually spaced locations 


### which genes are non tsgs/full-gene coverage

Cn_Hotspot_genes <- mouseCancerGenesV4$MGI.symbol[which(mouseCancerGenesV4$ampliconPurpose == "CopyNumber/Hotspot")]
Cn_Hotspot_genes <- c(Cn_Hotspot_genes, "Cdkn1a")

non_uniProtCanonical.reduced <- Cn_Hotspot_genes[-which(toupper(Cn_Hotspot_genes) %in% unique(uniProtCanonical.reduced$V1))]
non_uniProt_tmp <- uniProtCanonical$V2[which(uniProtCanonical$V1 %in% toupper(non_uniProtCanonical.reduced))]

### so not all of the proteins from uniprot are in our reduced list, so I added ones that aren't to it through google
### also P21 is cdkn1a
Cn_Hotspot_genes_proteins <- c(uniProtCanonical.reduced$V2,non_uniProt_tmp, "H3BKG7", "O88898", "Q8BRH4",
                               "P55200", "P39689")

Cn_Hotspot_genes_proteins2 <- c(uniProtCanonical$V2, "H3BKG7", "O88898", "Q8BRH4",
                                "P55200", "P39689")
###tried to reduce number of transcripts filtered 

allCds <- exonsBy(edb, filter = list(GeneNameFilter(Cn_Hotspot_genes)))

allCds.sortedExons <- allCds[order(unlist(lapply(allCds, length)),decreasing = TRUE)]

dummyList2 <- NULL
for (i in seq_along(allCds.sortedExons)) {
  dummyList2 <- c(dummyList2, allCds.sortedExons[[i]]$gene_name[1])
}

### seems cdk1n1a is listed twice along with P21 ... so we can reudce it down to that
allCds.longestExons <- allCds.sortedExons[!duplicated(dummyList2)]


## iteratively making 10 bed locations per oncogene
finalBed <- NULL
for (i in seq_along(names(allCds.longestExons))) {
  print(i)
  listOfIntervalIdx <- NULL
  tmpGeneOfInt <- allCds.longestExons[[i]]
  tmpTotalLength <- sum(allCds.longestExons[[i]]@ranges@width)
  tmpIncrement <- floor(tmpTotalLength/13)
  tmpInterval <- cumsum(allCds.longestExons[[i]]@ranges@width)
  for (k in 1:12) {
    tmpIdx <- k * tmpIncrement
    #tmpPost <- findInterval(tmpIdx, tmpInterval)
    tmpPost <- findInterval(tmpIdx, tmpInterval) + 1
    listOfIntervalIdx <- c(listOfIntervalIdx, tmpPost)
  }
  tmpTableExons <- names(table(listOfIntervalIdx))
  tmpTableCounts <- table(listOfIntervalIdx)
  tmpBed <- NULL
  for (j in seq_along(tmpTableExons)) {
    print(j)
    tmpRangeStart <- tmpGeneOfInt@ranges[j]@start
    tmpDist <-  floor(tmpGeneOfInt@ranges[j]@width/3)
    tmpChr <- as.numeric(as.character(tmpGeneOfInt@seqnames)[1])
    tmpGeneName <- tmpGeneOfInt@elementMetadata$gene_name[1]
    for (l in 1:tmpTableCounts[j]) {
      tmpPos <- tmpRangeStart + tmpDist * l
      tmpBedStart <- tmpPos - 50
      tmpBedEnd <- tmpPos + 50
      tmpBedLine <- c(tmpChr, tmpBedStart, tmpBedEnd, tmpGeneName )
      finalBed <- rbind(finalBed, tmpBedLine)
    }
  }
}

rownames(finalBed) <- NULL
finalBed <- data.frame(finalBed, stringsAsFactors = FALSE)
colnames(finalBed) <- c("Chr", "Start", "Stop", "Gene")
finalBed$Chr[which(is.na(finalBed$Chr))] <- "X"
finalBed$Chr <- paste0("chr",finalBed$Chr)

#write.table(finalBed,"/home/kevhu/data/20190219oncogeneMousePositions_12AmpPer.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


### maybe it would be better to (1) remove the amplicons from a list that cover a hotspot and then use closest on the nonOverlap list
### in order to find a replacement amplicon so I don't remove ones that are effectively in use already
### afterwards add ATM to TSGs. grab all exons for tsgs to add to bed -> then I can make a bed track layer for everyone to see


bedIntersectV2Overlap <- read.table("/home/kevhu/data/20190701intersectMouseSnps.bed", header = FALSE,
                                       stringsAsFactors = FALSE, sep = "\t")

bedIntersectV2Overlap_sorted <- bedIntersectV2Overlap[order(bedIntersectV2Overlap$V9),]
bedIntersectV2Overlap_sorted_red <- bedIntersectV2Overlap_sorted[!duplicated(paste0(bedIntersectV2Overlap_sorted$V4,
                                                                                          bedIntersectV2Overlap_sorted$V5,
                                                                                          bedIntersectV2Overlap_sorted$V6)),]

bedIntersect_dupliacted <- bedIntersectV2Overlap_sorted[duplicated(paste0(bedIntersectV2Overlap_sorted$V6,
                                                                             bedIntersectV2Overlap_sorted$V7,
                                                                             bedIntersectV2Overlap_sorted$V8)),]

#bedIntersect_duplicated2 <- bedIntersect_dupliacted[,c(1:4)]
#colnames(bedIntersect_duplicated2) <- c("Chr","Start","End","Gene")
finalBed.key <- paste0(finalBed$Chr,finalBed$Start,finalBed$Stop)

interchangeableAmps <- finalBed
interchangeableAmps <- interchangeableAmps[-which(finalBed.key %in% paste0(bedIntersectV2Overlap_sorted_red$V1,
                                                         bedIntersectV2Overlap_sorted_red$V2,
                                                         bedIntersectV2Overlap_sorted_red$V3)),] 


colnames(bedIntersectV2nonOverlap_sorted_red)[c(5:8)] <- colnames(finalBedV2)
finalBedV2 <- rbind(finalBedV2, bedIntersectV2nonOverlap_sorted_red[,c(5:8)])


### point of the fgfr2 is a single amplicon can cover almost all called variants so I decided to just create amplicons that spanned them
### rather than have it mess up the algorithm below
interchangeableAmps_key <- paste0(interchangeableAmps$Chr, interchangeableAmps$Start, interchangeableAmps$Stop)

nonOverlap <- read.table("/home/kevhu/data/20190701intersectNoOverlapMouseSnps.bed", header = FALSE,
                         stringsAsFactors = FALSE, sep = "\t")

fgfr2Sites <- data.frame("V1" = rep("chr7",4),
                         "V2" = c(130167668, 130172770, 130172850, 130261703),
                         "V3" = c(130167770, 130172872, 130172950, 130261805),
                         "V4" = rep("Fgfr2", 4), stringsAsFactors = FALSE)
ptprdSites <- data.frame("V1" = rep("chr4",2),
                         "V2" = c(76128690, 76084458),
                         "V3" = c(76128790, 76084458),
                         "V4" = rep("Ptprd", 2), stringsAsFactors = FALSE)

nonOverlap <- nonOverlap[-grep("FGFR2", nonOverlap$V4),]
nonOverlap <- nonOverlap[-grep("PTPRD", nonOverlap$V4),]

nonOverlap <- rbind(nonOverlap, fgfr2Sites, ptprdSites)


nonOverlap_genes <- nonOverlap[-which(nonOverlap$V4 == "snp:LOH"),]
nonOverlap_genes <- nonOverlap_genes[-which(nonOverlap_genes$V4 == "snp:StrainEst"),]
nonOverlap_genes$V5 <- toupper(str_remove(nonOverlap_genes$V4, "\\|.*"))
nonOverlap_genes$V2 <- as.numeric(nonOverlap_genes$V2)

interchangeableAmps$Gene <- toupper(interchangeableAmps$Gene)
interchangeableAmps$Start <- as.numeric(interchangeableAmps$Start)

library(DescTools)



interchangeableAmps_notUsed <- interchangeableAmps[-which(interchangeableAmps$Gene %in% unique(nonOverlap_genes$V5)),]
interchangeableAmps_pre <- interchangeableAmps[which(interchangeableAmps$Gene %in% unique(nonOverlap_genes$V5)),]

interchangeableAmps_post <- NULL
for (i in seq_along(unique(nonOverlap_genes$V5))) {
  tmpDf <- nonOverlap_genes[which(nonOverlap_genes$V5 == unique(nonOverlap_genes$V5)[i]),]
  tmpDf2 <- interchangeableAmps_pre[which(interchangeableAmps_pre$Gene == unique(nonOverlap_genes$V5)[i]),]
  tmpDf3 <- interchangeableAmps_pre[which(interchangeableAmps_pre$Gene == unique(nonOverlap_genes$V5)[i]),]
  listOfReps <- NULL
  for (j in 1:nrow(tmpDf)) {
    tmpIdx <- which(tmpDf2$Start == Closest(tmpDf2$Start, tmpDf$V2[j]))
    listOfReps <- c(listOfReps, tmpIdx)
    tmpDf2[-tmpIdx,]
  }
  tmpDf3[listOfReps,] <- tmpDf
  interchangeableAmps_post <- rbind(interchangeableAmps_post, tmpDf3)
}



### combines amps that were already intersected, amps that were replaceable and changed, and amps that were replaceable but not changed

finalBedV3 <- rbind(interchangeableAmps_notUsed, interchangeableAmps_post,
                    finalBed[which(finalBed.key %in% paste0(bedIntersectV2Overlap_sorted_red$V1,
                                                             bedIntersectV2Overlap_sorted_red$V2,
                                                             bedIntersectV2Overlap_sorted_red$V3)),] )

finalBedV3 <- finalBedV3[-which(finalBedV3$Gene == "ATM"),]

tsgGenes <- mouseCancerGenesV4$MGI.symbol[which(mouseCancerGenesV4$ampliconPurpose == "CDS")]
tsgGenes <- c(tsgGenes, "Atm")
tsg.cds <- cdsBy(edb, filter = list(GeneNameFilter(tsgGenes)))
tsg.cds_sorted <- tsg.cds[order(unlist(lapply(tsg.cds, length)), decreasing = TRUE)]

dummyList <- NULL
for (i in seq_along(tsg.cds_sorted)) {
  dummyList <- c(dummyList, tsg.cds_sorted[[i]]$gene_name[1])
}


tsg.cds_sorted_red <- tsg.cds_sorted[!duplicated(dummyList)]

tsg.bed <- NULL
for (i in seq_along(tsg.cds_sorted_red)) {
  tmpObj <- tsg.cds_sorted_red[[i]]
  tmpRange <- as.character(tmpObj@ranges)
  tmpChr <- rep(as.character(tmpObj@seqnames)[1], length(tmpRange))
  tmpGeneName <- rep(as.character(tmpObj@elementMetadata$gene_name)[1], length(tmpRange))
  tmpTable <- cbind(tmpChr, tmpRange,tmpGeneName)
  tsg.bed <- rbind(tsg.bed, tmpTable)
}


tsg.chr <- paste0("chr", tsg.bed[,1])
tsg.start <- unlist(lapply(strsplit(tsg.bed[,2], "-"), '[[',1))
tsg.end <- unlist(lapply(strsplit(tsg.bed[,2], "-"), '[[',2))
tsg.tmpTable <- data.frame("Chr" = tsg.chr,
                           "Start" = tsg.start, 
                           "Stop" = tsg.end,
                           "Gene" = tsg.bed[,3], stringsAsFactors = FALSE)

finalBedV4 <- rbind(finalBedV3, tsg.tmpTable)
finalBedV4$Start <- as.numeric(finalBedV4$Start)
finalBedV4$Stop <- as.numeric(finalBedV4$Stop)
finalBedV4$Gene <- toupper(finalBedV4$Gene)

### total targets from lifetech  below
###
### added ATM, Cycld + Ptch2 to TSG gene list. need to edit the mouse table for this. change labels of some of the oncogenes.
### need to add gli1 too ...

gli1.cds <- cdsBy(edb, filter = list(GeneNameFilter("Gli1")))

longest.gli <- gli1.cds[[1]]@ranges

gli_longest_cds <- NULL
for (i in seq_along(longest.gli)) {
  tmpChr <- "chr10"
  tmpGene <- "Gli1"
  tmpCds <- round(longest.gli@start[i]) + round(longest.gli@width[i]/2)
  tmpCds_start <- tmpCds - 50 
  tmpCds_end <- tmpCds + 50
  if (i == 5) {
    next()
  }
  tmpTable <- c(tmpChr, as.numeric(tmpCds_start), as.numeric(tmpCds_end), tmpGene)
  gli_longest_cds <- rbind(gli_longest_cds, tmpTable)
}

colnames(gli_longest_cds) <- colnames(finalBedV4)
finalBedV4 <- rbind(finalBedV4 , gli_longest_cds)


### V4 is everything and V3 is just tsgs
finalBedV3 <- rbind(finalBedV3, gli_longest_cds)
finalBedV3$Start <- as.numeric(finalBedV3$Start)
finalBedV3$Stop <- as.numeric(finalBedV3$Stop)
finalBedV3$Stop[which(finalBedV3$Start == finalBedV3$Stop)] <- finalBedV3$Stop[which(finalBedV3$Start == finalBedV3$Stop)] + 1


finalBedV3$Gene  <- toupper(finalBedV3$Gene)

finalBedV3 <- finalBedV3[order(finalBedV3$Chr, finalBedV3$Start, finalBedV3$Stop),]

#write.table(finalBedV3, "/home/kevhu/data/20190701oncogenesContig.bed", col.names = FALSE, row.names = FALSE, quote = FALSE,
#            sep = "\t")


### separate bed for variants
nonOverlap_snps <- nonOverlap[c(which(nonOverlap$V4 == "snp:LOH"), which(nonOverlap$V4 == "snp:StrainEst")),]
nonOverlap_snps$V3 <- nonOverlap_snps$V3 + 1

#write.table(nonOverlap_snps, "/home/kevhu/data/20190701snps.bed", col.names = FALSE, row.names = FALSE, quote = FALSE,
#            sep = "\t")


### 20190225: remaking the figure of amplicons spanning all the chromosomes
###
###


mouseCancerGenesV4$Annotation[which(mouseCancerGenesV4$MGI.symbol == "Atm")] <- "TumorSupressorGene"
mouseCancerGenesV4$ampliconPurpose[which(mouseCancerGenesV4$MGI.symbol == "Atm")] <- "CDS"



#write.table(mouseCancerGenesV4, "/home/kevhu/data/20190225mousePanelAnnotationsV2.txt", sep = "\t", quote = FALSE,
#            row.names = FALSE, col.names = TRUE)


#### edit two lines below so I get variants found theoretically on panel by tiling tsg. combine with list
###
###

mouseConvertedTargets <- read.table("/home/kevhu/data/20190215mouseOncoGeneHotspots.bed", sep = "\t",
                                    header = FALSE, stringsAsFactors = FALSE)
### 148 unique areas
listOfDetectedVars <- c(unique(mouseConvertedTargets$V4), unique(paste0(mskHotspots.Cosmic.onlyTsgs$Gene,
                                                                        "|", mskHotspots.Cosmic.onlyTsgs$Residue)))
mskHotspots.Cosmic_in <- mskHotspots.Cosmic[which(mskHotspots.Cosmic$Gene %in% mouseCancerGenesV4$HGNC.symbol),]

### number of hotspots we can detect vs numbers listed from analysis
431/892


mouseHotspotGenes_Tsgs <- mouseCancerGenesV4$HGNC.symbol[-which(mouseCancerGenesV4$ampliconPurpose == "CopyNumber/Hotspot")]
mskHotspots.Cosmic.onlyTsgs <- mskHotspots.Cosmic[which(mskHotspots.Cosmic$Gene %in% mouseHotspotGenes_Tsgs),]
mskHotspots.Cosmic.onlyOncogenes2 <- mskHotspots.Cosmic.onlyOncogenes
mskHotspots.Cosmic.onlyOncogenes2 <- mskHotspots.Cosmic.onlyOncogenes2[-which(mskHotspots.Cosmic.onlyOncogenes2$Gene == "ATM"),]

mouseConvertedTargetsV2 <- data.frame("Hotspot Residue" = mouseConvertedTargets$V4, "Chromosome" = mouseConvertedTargets$V1,
                                      "Start" = mouseConvertedTargets$V2, "End" = mouseConvertedTargets$V3,stringsAsFactors = FALSE)

### re-final price
4047 * 8
4047 * 9
4047 * 10


#write.table(mouseConvertedTargetsV2, "/home/kevhu/data/20190225mouseConvertedOncogeneHotspots.txt", col.names = TRUE,
#            row.names = FALSE, quote = FALSE, sep = "\t")


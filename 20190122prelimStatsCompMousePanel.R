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


mskGeneList <- read_xlsx("/home/kevhu/data/mskImpactGeneListSuppTab4.xlsx", skip = 2)
mskGeneList2 <- unique(mskGeneList$Gene)
genome <- TxDb.Hsapiens.UCSC.hg19.knownGene
conversionTab <- as.list(org.Hs.egALIAS2EG)

mskGeneList.entrezID <- unlist(conversionTab[which(names(conversionTab) %in% mskGeneList2)])

genome.subset <-  genes(genome)[which(genes(genome)$gene_id %in% mskGeneList.entrezID),]

boxplot(mskGeneList$stop - mskGeneList$start)

### things to think about - we don't know where the hotspots would be analogous for mice as compared to humans right?
### so making a hotspot for another species seems increasingly hard to do 
### look at a few of the genes homologues and see how likely the hotspots are to be in the same space
### my general thoughts on this is firs tto explore a subset of panels i.e kathy's shows one of two things
### targeted panels can detect genomic instability in mice similar to what is found in humans
### 


bccGeneList <- c("FBXW7", "PPP6C", "STK19", "CASP8","RB1","KNSTRN","ERRB2",
                 "NOTCH2", "NOTCH1", "ARID1A", "SMO", "TP53", "PTCH1", "SUFU")
bccGeneList.entrez <- unlist(conversionTab[which(names(conversionTab) %in% bccGeneList)])
genome.bcc <- genes(genome)[which(genes(genome)$gene_id %in% bccGeneList.entrez),]
bccGenes.exons <- subsetByOverlaps(exons(genome), genome.bcc)
sum(bccGenes.exons@ranges@width)
sum(bccGenes.exons@ranges@width)/120

(707 + 305 + 366 + 479 + 928 + 316 + 1255 + 2471 + 2555 + 2285 + 787 + 393 + 1447 + 484) * 3 / 120 

### converting to mouse and doing similar thing
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
biomartFilters <- listFilters(ensembl)

getBM(attributes = c("ensembl_gene_id",
                     "external_gene_name"),
      filters = "with_entrezgene",
      values = names(bccGeneList.entrez),
      mart = ensembl)

###
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


### ensembl method below was trying to identify the genomic locatuions of the mutations. however, since there is
### multiple isofroms of proteins ... it made it difficult. Pulled locations via COSMIC annotations seen below

#unique(mskHotspots$Gene)[-which(unique(mskHotspots$Gene) %in% toupper(mouseCancerGenesV3$MGI.symbol))]
#ensembl.hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#queryHotspots <- mouseCancerGenesV3$HGNC.symbol[which(grepl("Hotspot", mouseCancerGenesV3$ampliconPurpose))]

#ensemblHotspots <- getBM(attributes = c("ensembl_peptide_id","cdna_coding_start", "cdna_coding_end"),
#                         filters = "hgnc_symbol", values = queryHotspots, mart = ensembl.hs)

#ensemblHotspots2 <- getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id","hgnc_symbol"),
#                          filters = "hgnc_symbol", values = queryHotspots, mart = ensembl.hs)

#ensemblHotspotsCombined <- merge(x = ensemblHotspots, ensemblHotspots2, by = "ensembl_peptide_id")


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

write.table(mskHotspots.Cosmic,"/home/kevhu/data/20190201mouseHotspots.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


### this shows all hotspot mutations found on gene panel
hotspotOnly <- unique(mskHotspots.Cosmic$Gene)[-which(unique(mskHotspots.Cosmic$Gene) %in% unique(c(mouseCancerGenesV3$HGNC.symbol, toupper(mouseCancerGenesV3$MGI.symbol))))]

### below was only used because some of the non-Cosmic genes were not on the panel
#tmp3 <- tolower(hotspotOnly)
#substr(tmp3,1,1) <- toupper(substr(tmp3,1,1))
#tmp3[4] <- "Eloc"
#tmp3[26] <- "Bat1b"
#tmp3[27] <- "Gm8909"
#tmp3[40] <- "Ccnq"
#tmp3[42] <- "Nsd2"
#tmp4 <- cbind(hotspotOnly, tmp3, rep("Oncogene", length(hotspotOnly)),
#                               rep("HotspotOnly", length(hotspotOnly)))
#colnames(tmp4) <- colnames(mouseCancerGenesV3)
#mouseCancerGenesV3 <- rbind(mouseCancerGenesV3, tmp4)

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



write.table(mouseCancerGenesV3,"/home/kevhu/data/20190129preliminaryMousePanelListCosmic.txt", sep = "\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE)


### getting the accession numbers so I can run the protein blast - i have a bookmark on a page that does this
### need to try and grab canonical seqeunces

ensembl.hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomartFilters <- listFilters(ensembl.hs)
ensemblHotspots <- getBM(attributes = c("hgnc_symbol","entrezgene","ensembl_gene_id"),
                         filters = "hgnc_symbol", values = mouseCancerGenesV3$HGNC.symbol, mart = ensembl.hs)

res <- select(org.Hs.eg.db, keys = as.character(ensemblHotspots$entrezgene), columns = c("ACCNUM", "UNIPROT",
                                                          "SYMBOL"), keytype = "ENTREZID")

### best to to find canonical isoform + protein produced by it: end result only gave us 56 unique genes.
### pretty bad .... somehow doing it directly from ucscKG give it worst of 31 ...
#genCode <- read.table("/home/kevhu/data/knownToGencodeV20.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#refSeq.canonical <- read.table("/home/kevhu/data/knownCanonical.txt", sep = "\t", header = FALSE,
#                               stringsAsFactors = FALSE)

#res2 <- res[which(res$UNIPROT %in% uniProtCanonical$V2),]
###
###will take uniprot list and then get appropriate fasta files for it to blast against

### get list of gene names then put it into the retrieve unitprot .. for canonical accessioin or structure



### quick calculation @ $8 dollars per amplicon
numberAmpsCN <- length(which(mouseCancerGenesV3$ampliconPurpose == "CopyNumber")) * 10 

#mouseCancerGenesV3[which(mouseCancerGenesV3$Annotation == "TumorSupressorGene"),]
numberAmpsTSGs <- 701
numberAmpsTSGsUTRs <- 1059
newFearonAmps <- 1275

numberOfHotspotsInGenes <- table(mskHotspots.Cosmic$Gene[which(mskHotspots.Cosmic$Gene %in% mouseCancerGenesV3$HGNC.symbol[which(mouseCancerGenesV3$ampliconPurpose == "CopyNumber/Hotspot")])])
numberOfHotspotsInGenes2 <- numberOfHotspotsInGenes 
numberOfHotspotsInGenes2[which(numberOfHotspotsInGenes2 < 10)] <- 10
numberCnHotspot <- sum(numberOfHotspotsInGenes2)

totalAmps <- numberAmpsCN + numberAmpsTSGs + numberCnHotspot
totalAmps * 8
totalAmps2 <- numberAmpsCN + newFearonAmps + numberCnHotspot
totalAmps2 * 8



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

### need to expand the single base regions to ~100bp regions in order to use the liftover tools - allow can be used to
### collapse multiple hotspots onto one contig; need to convert string to genomic coordinates - in the end
### didn't need to convert because the input needed was in format of what cosmic gave

chrom <- NULL
startCoord <- NULL
endCoord <- NULL
for(i in mskHotspots.Cosmic$GenomicLocation){
  tmp <- strsplit(i, ":")
  chrom <- c(chrom, tmp[[1]][1])
  tmp2 <- strsplit(tmp[[1]][2], "-")
  startCoord <- c(startCoord, tmp2[[1]][1])
  endCoord <- c(endCoord, tmp2[[1]][2])
}


startCoord <- as.numeric(startCoord)
startCoord <- startCoord - 50
endCoord <- as.numeric(endCoord)
endCoord <- endCoord + 50

chrom2 <- chrom
chrom2[which(chrom2 == "23")] <- "X"
chrom3 <- paste0("chr", chrom2)

newBed <- data.frame("Chr" = chrom, "Start" = startCoord, "End" = endCoord,  "COSMIC" = mskHotspots.Cosmic$cosmicID,stringsAsFactors = FALSE )
newBed2<- data.frame("chr" = chrom3, "start" = startCoord, "stop" = endCoord)

write.table(newBed2,"/home/kevhu/data/20190205mouseHotspotHumanV2.bed", quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")

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

chrom <- NULL
startCoord <- NULL
endCoord <- NULL
for(i in mskHotspots.Cosmic.onlyOncogenes$GenomicLocation){
  tmp <- strsplit(i, ":")
  chrom <- c(chrom, tmp[[1]][1])
  tmp2 <- strsplit(tmp[[1]][2], "-")
  startCoord <- c(startCoord, tmp2[[1]][1])
  endCoord <- c(endCoord, tmp2[[1]][2])
}


startCoord <- as.numeric(startCoord)
startCoord <- startCoord - 25
endCoord <- as.numeric(endCoord)
endCoord <- endCoord + 25

chrom2 <- chrom
chrom2[which(chrom2 == "23")] <- "X"
chrom3 <- paste0("chr", chrom2)

bedReducedHotspot <- data.frame(chrom3, startCoord, endCoord, mskHotspots.Cosmic.onlyOncogenes$cosmicID)

write.table(bedReducedHotspot,"/home/kevhu/data/20190211mouseCosmicOncogenes.bed", quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")


###
newBlast.res <- read.table("/home/kevhu/data/20190211blastSeqAlignMouseV2.txt",sep = "\t", stringsAsFactors = FALSE)

newBlast.name <- NULL
for (i in newBlast.res$V1) {
  newBlast.name <- c(newBlast.name, strsplit(i, "\\|")[[1]][3])
}
newBlast.res$V5 <- newBlast.name

kgXref <- read.table("/home/kevhu/data/kgXref.txt", sep = "\t", stringsAsFactors = FALSE)

tmpNewBlastV6 <- uniProtCanonical$V1[match(newBlast.res$V5,uniProtCanonical$V3)]

newBlast.res$V6 <- tmpNewBlastV6
### idea might be to just iterate through mskHotspot list and create 3 sets of mouse amino acids for each ....
### i.e aa1 aa2 aa3 target aa4 aa5 aa6, aa1-aa6 target, target aa1-aa6, then blastx these guys back to mouse nucleotides from 
### ther amino acid sequence

### if this doesn't work, get same sequence as above, but for humans and tblastn it agaisnt mouse
### since the output for mouse is whatever mapped, there are gaps ... so ... need to do reversed or
### just look at the fasta file 


sortedVars <- mskHotspots.Cosmic.onlyOncogenes[,c("Gene","Residue")]
sortedVars.position <- substr(sortedVars$Residue, 2, nchar(sortedVars$Residue))
sortedVars.AA <- substr(sortedVars$Residue, 1, 1)

sortedVars$aminoAcid <- sortedVars.AA
sortedVars$position <- as.numeric(sortedVars.position)

rownames(sortedVars) <- NULL

blastList <- NULL
for(i in 1:nrow(sortedVars)){
  tmpAAseq <- newBlast.res$V4[which(newBlast.res$V6  == sortedVars$Gene[i])]
  for (k in seq_along(tmpAAseq)) {
    tmpAAseq.1 <- tmpAAseq[k]
    for (j in 1:3) {
      if (j == 1) {
        tmpAA <- substr(tmpAAseq.1, (sortedVars$position - 6), sortedVars$position)
        blastList <- c(blastList, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList <- c(blastList, tmpAA)
      }
      if (j == 2) {
        tmpAA <- substr(tmpAAseq.1, sortedVars$position,  (sortedVars$position + 6))
        blastList <- c(blastList, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList <- c(blastList, tmpAA)
      }
      if (j == 3) {
        tmpAA <- substr(tmpAAseq.1, (sortedVars$position - 3), (sortedVars$position + 3))
        blastList <- c(blastList, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList <- c(blastList, tmpAA)
      }
    }
  }
}

con<-file("/home/kevhu/data/20190212customAABlastXMouse.fa")
writeLines(blastList, con)
close(con)

### doing this for human

blastList.human <- NULL
for(i in 1:nrow(sortedVars)){
  tmpAAseq <- newBlast.res$V3[which(newBlast.res$V6  == sortedVars$Gene[i])]
  for (k in seq_along(tmpAAseq)) {
    tmpAAseq.1 <- tmpAAseq[k]
    for (j in 1:3) {
      if (j == 1) {
        tmpAA <- tryCatch(substr(tmpAAseq.1, (sortedVars$position - 6), sortedVars$position), error = function(x) return(NULL))
        if (is.null(tmpAA)) {
          next()
        }
        blastList.human <- c(blastList.human, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList.human <- c(blastList.human, tmpAA)
      }
      if (j == 2) {
        tmpAA <- tryCatch(substr(tmpAAseq.1, sortedVars$position,  (sortedVars$position + 6)), error = function(x) return(NULL))
        if (is.null(tmpAA)) {
          next()
        }
        blastList.human <- c(blastList.human, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList.human <- c(blastList.human, tmpAA)
      }
      if (j == 3) {
        tmpAA <- tryCatch(substr(tmpAAseq.1, (sortedVars$position - 3), (sortedVars$position + 3)), error = function(x) return(NULL))
        if (is.null(tmpAA)) {
          next()
        }
        blastList.human <- c(blastList.human, paste0(">", sortedVars$Gene[i], ".",sortedVars$Residue[i]))
        blastList.human <- c(blastList.human, tmpAA)
      }
    }
  }
}

#con.human <-file("/home/kevhu/data/20190212customAABlastXHuman.fa")
#writeLines(blastList.human, con.human)
#close(con.human)


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


uniProtCanonical <- read.table("/home/kevhu/data/20190214wholeCanonicalMouseProt.txt", skip = 1, header = FALSE,
                               stringsAsFactors = FALSE, sep = "\t")
load("/home/kevhu/data/20190214mskHotspots.onlyOnco.df.Robj")
mouseCancerGenesV4 <- read.table("/home/kevhu/data/20190214mouseCancerGenesV4.txt", header = TRUE,
                                 stringsAsFactors = FALSE, sep = "\t")

#canProtName <- uniProtCanonical$V2[match(mskHotspots.Cosmic.onlyOncogenes$Gene,uniProtCanonical$V1)]


matchedVars <- mskHotspots.Cosmic.onlyOncogenes$Residue[match(uniProtCanonical$V1, mskHotspots.Cosmic.onlyOncogenes$Gene)]

uniProtCanonical$residueLocation <- matchedVars
#mskHotspots.Cosmic.onlyOncogenes$protName[316] <-"P46100"
#mskHotspots.Cosmic.onlyOncogenes$protName[560] <-"P46100"
#mskHotspots.Cosmic.onlyOncogenes$protName[350] <- "Q9QXA3"

# below is example of it working 
testNRAS <- IRanges(start = 59, end = 61, names = "Q9D091")
testNRAS.out <- proteinToGenome(testNRAS, EnsDb.Mmusculus.v79, idType = "uniprot_id")
as.character(testNRAS.out[[1]]@ranges)
as.character(testNRAS.out[[1]]@seqnames@values)

justResPosition <- as.numeric(substr(uniProtCanonical$residueLocation,
                          2, nchar(mskHotspots.Cosmic.onlyOncogenes$Residue)))

chromConvert <- NULL
seqConvert <- NULL
for (i in seq_along(uniProtCanonical$V2)) {
  print(i)
  tmpRange <- tryCatch(IRanges(start = (justResPosition[i]), end = (justResPosition[i]),names = uniProtCanonical$V2[i]),
                       warning = function(x) return(NULL))
  if (is.null(tmpRange)) {
    tmpChrom <- "NA"
    tmpSeq <- "NA"
    chromConvert <- c(chromConvert, tmpChrom)
    seqConvert <- c(seqConvert, tmpSeq)
    next()
  }
  print(tmpRange)
  tmpConv <- tryCatch(proteinToGenome(tmpRange, EnsDb.Mmusculus.v79, idType = "uniprot_id"),
                      warning = function(x) return(NULL))
  if (is.null(tmpConv)) {
    tmpChrom <- "NA"
    tmpSeq <- "NA"
    chromConvert <- c(chromConvert, tmpChrom)
    seqConvert <- c(seqConvert, tmpSeq)
    next()
  }
  print(tmpConv)
  
  if (length(tmpConv[[1]]) > 1) {
    tmpChrom <- as.character(tmpConv[[1]]@unlistData@seqnames@values[1])
    tmpSeq <- as.character(tmpConv[[1]]@unlistData@ranges[1])
  }
#  if (length(tmpConv[[1]]@ranges) > 1) {
#    tmpMax <- max(tmpConv[[1]]@ranges@width)
#    tmpMaxIdx <- which(tmpConv[[1]]@ranges@width == tmpMax)
#    tmpChrom <- as.character(tmpConv[[1]]@seqnames@values)
#    tmpSeq <- as.character(tmpConv[[1]]@ranges)[tmpMaxIdx]
#    chromConvert <- c(chromConvert, tmpChrom)
#    seqConvert <- c(seqConvert, tmpSeq)
#    next()
#  }
  else{
    tmpChrom <- as.character(tmpConv[[1]]@seqnames@values)
    tmpSeq <- as.character(tmpConv[[1]]@ranges)
  }
  #tmpSeq <- as.character(tmpConv[[1]]@ranges)
  chromConvert <- c(chromConvert, tmpChrom)
  seqConvert <- c(seqConvert, tmpSeq)
}


#tmpRange <- IRanges(start = (justResPosition[24] - 2), end = (justResPosition[24] + 2),names = "Q9JJU4")
#tmpConv <- try(proteinToGenome(tmpRange, EnsDb.Mmusculus.v79, idType = "uniprot_id"))
#as.character(tmpConv[[1]]@unlistData@ranges[1])

uniProtCanonical$locationChr <- chromConvert
uniProtCanonical$locationRange <- seqConvert

uniProtCanonical.reduced <- uniProtCanonical[-which(uniProtCanonical$locationChr == "NA"),]
uniProtCanonical.reduced <- uniProtCanonical.reduced[!duplicated(uniProtCanonical.reduced$locationRange),]

tmpLocation <- strsplit(uniProtCanonical.reduced$locationRange, "-")
oncoHotspotBedF1 <- paste0("chr",uniProtCanonical.reduced$locationChr)
oncoHotspotBedF2 <- unlist(lapply(tmpLocation, '[[',1))
oncoHotspotBedF3 <- unlist(lapply(tmpLocation, '[[',2))
oncoHotspotBedF4 <- paste0(uniProtCanonical.reduced$V1, "|",uniProtCanonical.reduced$residueLocation)

mouseHotspotBed <- data.frame("Chr" = oncoHotspotBedF1, "Start" = oncoHotspotBedF2,
                              "End" = oncoHotspotBedF3, "Residue" = oncoHotspotBedF4,
                              stringsAsFactors = FALSE)

mouseHotspotBed3Col <- data.frame("Chr" = oncoHotspotBedF1, "Start" = oncoHotspotBedF2,
                              "End" = oncoHotspotBedF3,
                              stringsAsFactors = FALSE)

write.table(mouseHotspotBed3Col, "/home/kevhu/data/20190215mouseOncoGeneHotspots3Col.bed",
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

### may try process again but with a newer ens database
### tried and nothign changes really ... exactly same variatns are re-converted despite 3-4 years update ..

library(AnnotationHub)
ah <- AnnotationHub()
#query(ah, c("Mus musculus", "release-93"))
edb <- ah[["AH64944"]]
#ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb", 94))
#edb <- EnsDb(ahDb)

testNRAS <- IRanges(start = 61, end = 61, names = "Q9D091")
testNRAS.out <- proteinToGenome(testNRAS, edb, idType = "uniprot_id")
as.character(testNRAS.out[[1]]@unlistData@seqnames)
as.character(testNRAS.out[[1]]@unlistData@ranges)



chromConvert2 <- NULL
seqConvert2 <- NULL
for (i in seq_along(uniProtCanonical$V2)) {
  tmpRange <- tryCatch(IRanges(start = (justResPosition[i]), end = (justResPosition[i]),names = uniProtCanonical$V2[i]),
                       warning = function(x) return(NULL))
  if (is.null(tmpRange)) {
    tmpChrom <- "NA"
    tmpSeq <- "NA"
    chromConvert2 <- c(chromConvert2, tmpChrom)
    seqConvert2 <- c(seqConvert2, tmpSeq)
    next()
  }
  tmpConv <- tryCatch(proteinToGenome(tmpRange, edb, idType = "uniprot_id"),
                      warning = function(x) return(NULL))
  if (is.null(tmpConv)) {
    tmpChrom <- "NA"
    tmpSeq <- "NA"
    chromConvert2 <- c(chromConvert2, tmpChrom)
    seqConvert2 <- c(seqConvert2, tmpSeq)
    next()
  }
  if (length(tmpConv[[1]]) > 1) {
    tmpChrom <- as.character(tmpConv[[1]]@unlistData@seqnames@values[1])
    tmpSeq <- as.character(tmpConv[[1]]@unlistData@ranges[1])
  }
  else{
    tmpChrom <- as.character(tmpConv[[1]]@seqnames@values)
    tmpSeq <- as.character(tmpConv[[1]]@ranges)
  }
  chromConvert2 <- c(chromConvert2, tmpChrom)
  seqConvert2 <- c(seqConvert2, tmpSeq)
}

uniProtCanonicalV89 <- uniProtCanonical
uniProtCanonicalV89$locationChr <- chromConvert2
uniProtCanonicalV89$locationRange <- seqConvert2

uniProtCanonicalV89.reduced <- uniProtCanonicalV89[-which(uniProtCanonicalV89$locationChr == "NA"),]
uniProtCanonicalV89.reduced <- uniProtCanonicalV89.reduced[!duplicated(uniProtCanonicalV89.reduced$locationRange),]



### next is figuring out how many amplicons of the hotspots cover ... also seems like there is some overlap of the hotspots ... mainly just fgfr2 ()
###
###

table(uniProtCanonical.reduced$V1)

### pseudo code for getting the 10 amplicons that tile the gene
### (1) grab genic/exonix sequence. (2) linearize if exonic (3) grab 10 individually spaced locations 
### 2019/02/18 - so I did some more thinking and research into the whole LOH thing ... originally I though I could just add SNPs into amplicons
### where we were solely detecting copy number ... i.e did not lie on a hotspot. this design was seemingly deployed in on the panGU
### only probelm with this is SNPs are generally polymorphic in humans and thus allow allele-specific loss to be done or theoretically done.
### For mice, if strains are isogenic or there is low positions of polymorphism especially in genes we target i.e look at trp53 or nras ..
### then doing so is not feasible with SNPS or other types of genetic variations
### exon rank is just the index system they employ ...


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

#allCds_transcriptsWithRankOneExon <- names(allCds)
#allCds2 <- cdsBy(edb, filter = list(GeneNameFilter(Cn_Hotspot_genes),
#                                   ExonRankFilter(10, "<")))

#allCds2 <- allCds2[which(names(allCds2) %in% allCds_transcriptsWithRankOneExon)]

dummyList <- NULL
for (i in seq_along(allCds)) {
  dummyList <- c(dummyList, allCds[[i]]$gene_name[1])
}

### sanity check to see if everything is there ... P21 is Cdkn1a + not listed in order of most exons. quick fix below that
dummyList1.1 <- unique(dummyList)
Cn_Hotspot_genes[-which(Cn_Hotspot_genes %in% dummyList1.1)]

allCds[[996]]
allCds[[997]]
allCds[[998]]

allCds.sortedExons <- allCds[order(unlist(lapply(allCds, length)),decreasing = TRUE)]

dummyList2 <- NULL
for (i in seq_along(allCds.sortedExons)) {
  dummyList2 <- c(dummyList2, allCds.sortedExons[[i]]$gene_name[1])
}

### seems cdk1n1a is listed twice along with P21 ... so we can reudce it down to that
allCds.longestExons <- allCds.sortedExons[!duplicated(dummyList2)]

floor(sum(allCds.longestExons[[1]]@ranges@width)/10)
cumsum(allCds.longestExons[[1]]@ranges@width)

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

write.table(finalBed,"/home/kevhu/data/20190219oncogeneMousePositions_12AmpPer.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


### read in results from intersect. find out which ones didn't cut it and then add amplicons for them
### need to add ATM to genes listed as TSGs

bedIntersect <- read.table("/home/kevhu/data/20190219mouseBedsOverlap_12AmpPer.txt", header = FALSE,
                           stringsAsFactors = FALSE, sep = "\t")

bedIntMatched <- paste0(bedIntersect$V5,bedIntersect$V6, bedIntersect$V7, bedIntersect$V8)
allMouseHotspotsKeys <- paste0(mouseHotspotBed$Chr, mouseHotspotBed$Start,
                               mouseHotspotBed$End,mouseHotspotBed$Residue)

mouseHotspotBed.notOverlap <- mouseHotspotBed[-which(allMouseHotspotsKeys %in% bedIntMatched),]

###padding the boundaries os 
newStart100 <- as.numeric(mouseHotspotBed.notOverlap$Start) - 50
newEnd100 <- as.numeric(mouseHotspotBed.notOverlap$End) + 50

hotspotBedNonOverlap <- data.frame("chr" = mouseHotspotBed.notOverlap$Chr, "start" = newStart100, 
                                   "end" = newEnd100, stringsAsFactors = FALSE)

write.table(hotspotBedNonOverlap,"/home/kevhu/data/20190219oncogeneHotspotsNonOverlap.bed",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


### maybe it would be better to (1) remove the amplicons from a list that cover a hotspot and then use closest on the nonOverlap list
### in order to find a replacement amplicon so I don't remove ones that are effectively in use already
### afterwards add ATM to TSGs. grab all exons for tsgs to add to bed -> then I can make a bed track layer for everyone to see

finalBed.key <- paste0(finalBed$Chr,finalBed$Start,finalBed$Stop)
bedIntersect.oncogeneBed <- paste0(bedIntersect$V1, bedIntersect$V2, bedIntersect$V3)
finalBed.nonOnverlap <- finalBed[-which(finalBed.key %in% bedIntersect.oncogeneBed),]
write.table(finalBed.nonOnverlap,"/home/kevhu/data/20190220oncogeneMousePos12AmpPerNonoverlap.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


bedIntersectV2nonOverlap <- read.table("/home/kevhu/data/20190220oncogeneHotspotNonoverlapClosestV2.txt", header = FALSE,
                                       stringsAsFactors = FALSE, sep = "\t")

bedIntersectV2nonOverlap <- bedIntersectV2nonOverlap[-which(bedIntersectV2nonOverlap$V7 == "Fgfr2"),]

#bedIntersectV2_keys <- paste0(bedIntersectV2nonOverlap$V1, bedIntersectV2nonOverlap$V2, bedIntersectV2nonOverlap$V3)

which(duplicated(paste0(bedIntersectV2nonOverlap$V4, bedIntersectV2nonOverlap$V5, bedIntersectV2nonOverlap$V6)))


#bedIntersectV3nonOverlap <- read.table("/home/kevhu/data/20190220oncogeneHotspotNonoverlapClosestV3.txt", header = FALSE,
#                                       stringsAsFactors = FALSE, sep = "\t")

#
#going to remove fgfr2 ones as they cause an error and will make custom amplicons for them
#bedIntersectV3nonOverlap <- bedIntersectV3nonOverlap[-which(bedIntersectV3nonOverlap$V7 == "Fgfr2"),]
#hotspotBedNonOverlap_keys <- paste0(hotspotBedNonOverlap$chr, hotspotBedNonOverlap$start, hotspotBedNonOverlap$end)
#tmpInt_keys <- paste0(bedIntersectV3nonOverlap)
#tmpBedTable <- bedIntersectV3nonOverlap
#bedInUniqueTable <- NULL

#for (i in seq_along(hotspotBedNonOverlap_keys)){
#  print(i)
#  tmpIntHotspot_keys <- paste0(tmpBedTable$V1, tmpBedTable$V2, tmpBedTable$V3)
#  tmpIntBed_keys <- paste0(tmpBedTable$V4, tmpBedTable$V5, tmpBedTable$V6)
#  tmpIdx <- which(tmpIntHotspot_keys %in% hotspotBedNonOverlap_keys[i])[1]
#  bedInUniqueTable <- rbind(bedInUniqueTable, tmpBedTable[tmpIdx,])
#  print(bedInUniqueTable)
#  tmpIdx2 <- which(tmpIntBed_keys %in% paste0(bedInUniqueTable$V4[i], bedInUniqueTable$V5[i], bedInUniqueTable$V6[i]))
#  tmpBedTable <- tmpBedTable[-tmpIdx2,]
#}


### might just be going to do some of these manually to make the process easier
### first keep ones that are closest and then manually change the further duplications


bedIntersectV2nonOverlap_sorted <- bedIntersectV2nonOverlap[order(bedIntersectV2nonOverlap$V8),]
bedIntersectV2nonOverlap_sorted_red <- bedIntersectV2nonOverlap_sorted[!duplicated(paste0(bedIntersectV2nonOverlap_sorted$V4,
                                                                                          bedIntersectV2nonOverlap_sorted$V5,
                                                                                          bedIntersectV2nonOverlap_sorted$V6)),]

bedIntersect_dupliacted <- bedIntersectV2nonOverlap_sorted[duplicated(paste0(bedIntersectV2nonOverlap_sorted$V4,
                                                                              bedIntersectV2nonOverlap_sorted$V5,
                                                                              bedIntersectV2nonOverlap_sorted$V6)),]

bedIntersect_duplicated2 <- bedIntersect_dupliacted[,c(1:3,7)]
colnames(bedIntersect_duplicated2) <- c("Chr","Start","End","Gene")
finalBedV2 <- finalBed
  
finalBedV2 <- finalBedV2[-which(finalBed.key %in% paste0(bedIntersectV2nonOverlap_sorted$V4,
                                         bedIntersectV2nonOverlap_sorted$V5,
                                         bedIntersectV2nonOverlap_sorted$V6)),] 

colnames(bedIntersectV2nonOverlap_sorted_red)[c(1:3,7)] <- colnames(finalBedV2)
finalBedV2 <- rbind(finalBedV2, bedIntersectV2nonOverlap_sorted_red[,c(1:3,7)])

finalBedV2_key <- paste0(finalBedV2$Chr, finalBedV2$Start, finalBedV2$Stop)
interchangeableAmps <- finalBedV2[-which(finalBedV2_key %in% c(bedIntersect.oncogeneBed,
                                                               paste0(bedIntersectV2nonOverlap_sorted$V1,
                                                                      bedIntersectV2nonOverlap_sorted$V2,
                                                                      bedIntersectV2nonOverlap_sorted$V3))),]

fgfr2Sites <- data.frame("Chr" = rep("chr7",4),
                         "Start" = c(130167668, 130172770, 130172850, 130261703),
                         "End" = c(130167770, 130172872, 130172950, 130261805),
                         "Gene" = rep("Fgfr2", 4), stringsAsFactors = FALSE)

interchangeableAmps[which(interchangeableAmps$Gene == "Fgfr2"),]

interchangeableAmps[c(613:616),] <- fgfr2Sites

for (i in seq_along(unique(bedIntersect_duplicated2$Gene))) {
  print(i)
  tmpGene <- unique(bedIntersect_duplicated2$Gene)[i]
  tmpReplacement <- bedIntersect_duplicated2[which(bedIntersect_duplicated2$Gene == tmpGene),]
  print(tmpReplacement)
  tmpIdx <- which(interchangeableAmps$Gene == tmpGene)
  tmpReplacementIdx <- tmpIdx[(length(tmpIdx) - nrow(tmpReplacement) + 1):(length(tmpIdx))]
  print(tmpReplacementIdx)
  interchangeableAmps[tmpReplacementIdx,] <- tmpReplacement
}



finalBedV3 <- rbind(interchangeableAmps, 
                    finalBedV2[which(finalBedV2_key %in% c(bedIntersect.oncogeneBed,
                                                            paste0(bedIntersectV2nonOverlap_sorted$V1,
                                                                   bedIntersectV2nonOverlap_sorted$V2,
                                                                   bedIntersectV2nonOverlap_sorted$V3))),])

finalBedV3 <- finalBedV3[-which(finalBedV3$Gene == "Atm"),]
  
write.table(finalBedV3, "/home/kevhu/data/20190221finalSetOfOncogeneContigs.bed", quote = FALSE,
            col.names = FALSE, row.names = FALSE, sep = "\t")

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
write.table(finalBedV4, "/home/kevhu/data/20190221mouseWholeBedV4.bed", col.names = FALSE,
            row.names = FALSE, quote = FALSE, sep = "\t")


### need to check coverage per chromosome and see whether enough of it is covered per chromosome in order to calculate 
### arm loss - trying with ideoViz package. can be done using granges so I need to convert granges 
finalBedV4 <- read.table("/home/kevhu/data/20190221mouseWholeBedV4.bed", header = FALSE,
                         sep = "\t", stringsAsFactors = FALSE)
finalBedV4$V2 <- as.numeric(finalBedV4$V2)
finalBedV4$V3 <- as.numeric(finalBedV4$V3)

library(RColorBrewer)
library(IdeoViz)

finalBedV4.granges <- GRanges(finalBedV4$V1, IRanges(finalBedV4$V2, finalBedV4$V3))

pallette <- rep("#ff0000", 5)
ideo_mm10 <- getIdeo("mm10")
chromosomes <- unique(finalBedV4$V1)
chrom_bins <- getBins(chromosomes, ideo_mm10, stepSize = 5*100*1000)

valueGrDat <- list(finalBedV4.granges,finalBedV4.granges,finalBedV4.granges,finalBedV4.granges,finalBedV4.granges)

pdf(file = "/home/kevhu/data/20190222finalBedV4Ideogram.pdf", useDingbats = TRUE, width = 25, height = 60)
plotOnIdeo(chrom = seqlevels(chrom_bins),
           ideoTable = ideo_mm10, 
           values_GR = valueGrDat,
           plotType = "seg_tracks", value_cols = "value",
           vertical = FALSE, col = pallette)
dev.off()


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

rownames(gli_longest_cds) <- NULL

tsgOncoGenesBed <- read.table("/home/kevhu/data/20190221finalSetOfOncogeneContigs.bed", header = FALSE, stringsAsFactors = FALSE,
           sep = "\t")

tsgOncoGenesBed <- rbind(tsgOncoGenesBed, gli_longest_cds)
write.table(tsgOncoGenesBed, "/home/kevhu/data/20190225oncogeneContigs.bed", col.names = FALSE, row.names = FALSE, quote = FALSE,
            sep = "\t")

### 20190225: remaking the figure of amplicons spanning all the chromosomes
###
###

finalBedV4 <- rbind(finalBedV4, gli_longest_cds)
finalBedV4$V2 <- as.numeric(finalBedV4$V2)
finalBedV4$V3 <- as.numeric(finalBedV4$V3)
finalBedV4.granges <- GRanges(finalBedV4$V1, IRanges(finalBedV4$V2, finalBedV4$V3))

pallette <- rep("#ff0000", 5)
ideo_mm10 <- getIdeo("mm10")
chromosomes <- unique(finalBedV4$V1)
chrom_bins <- getBins(chromosomes, ideo_mm10, stepSize = 5*100*1000)

valueGrDat <- list(finalBedV4.granges,finalBedV4.granges,finalBedV4.granges,finalBedV4.granges,finalBedV4.granges)

pdf(file = "/home/kevhu/data/20190225finalBedV4Ideogram.pdf", useDingbats = TRUE, width = 25, height = 60)
plotOnIdeo(chrom = seqlevels(chrom_bins),
           ideoTable = ideo_mm10, 
           values_GR = valueGrDat,
           plotType = "seg_tracks", value_cols = "value",
           vertical = FALSE, col = pallette)
dev.off()

mouseCancerGenesV4$Annotation[which(mouseCancerGenesV4$MGI.symbol == "Atm")] <- "TumorSupressorGene"
mouseCancerGenesV4$ampliconPurpose[which(mouseCancerGenesV4$MGI.symbol == "Atm")] <- "CDS"



write.table(mouseCancerGenesV4, "/home/kevhu/data/20190225mousePanelAnnotationsV2.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)


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

### re-final price with Gli10  
(3884 + 11) * 8

write.table(mouseConvertedTargetsV2, "/home/kevhu/data/20190225mouseConvertedOncogeneHotspots.txt", col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = "\t")

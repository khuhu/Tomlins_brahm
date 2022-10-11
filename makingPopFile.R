# tmpTable <- read.table("/br_z1/kevin_storage/ADMIXTURE/merge48ForAdmixture.fam", sep = " ", 
#                        stringsAsFactors = FALSE, header = FALSE)

tmpTable <- read.table("/br_z1/kevin_storage/ADMIXTURE/20220912merge48Autosomes.fam", sep = " ", 
                       stringsAsFactors = FALSE, header = FALSE)


popFile <- tmpTable[,1:2]
popFile$V2[grep("5516", popFile$V2)] <- "-"
popFile$V2[grep("129", popFile$V2)] <- "129"
popFile$V2[grep("BALB", popFile$V2)] <- "BALB"
popFile$V2[grep("C57B", popFile$V2)] <- "C57BL"
popFile$V2[grep("C3H", popFile$V2)] <- "C3H"

unique(popFile$V2)


# write.table(popFile[,2], "/br_z1/kevin_storage/ADMIXTURE/merge48ForAdmixture.pop", col.names = FALSE,
#             row.names = FALSE, quote = FALSE)


write.table(popFile[,2], " admixtureFilt.fam", col.names = FALSE,
            row.names = FALSE, quote = FALSE)

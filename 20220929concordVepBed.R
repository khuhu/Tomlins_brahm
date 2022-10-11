concordBed <- combinedVars_concord_goodsamps[grep("exonic|splicing", combinedVars_concord_goodsamps$Func.refGene),]
concordBed2 <- concordBed[, c("Chr", "Start", "End")]
colnames(concordBed2) <- c("CHROM", "POS", "END")
concordBed2$POS <- concordBed2$POS - 1

write.table(concordBed2, "/br_z1/kevin_storage/misc/20220930bbnConcordMutBed.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


variantTable3 <- variantTable2[grep("exonic|splicing", variantTable2$Func.refGene),]
variantTable4 <- variantTable3[, c("Chr", "Start", "End")]
colnames(variantTable4) <- c("CHROM", "POS", "END")
variantTable4$POS <- variantTable4$POS - 1

write.table(variantTable4, "/br_z1/kevin_storage/misc/20220930aaronBbnConcordMutBed.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


write.table(variantTable3, "/br_z1/kevin_storage/misc/20221004aaronVarList.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

library(argyle)
load("/br_z1/kevin_storage/misc/snps.gigamuga.Rdata")
load("/br_z1/kevin_storage/misc/clusters.gigamuga.Rdata")

snpRes <- read.beadstudio(prefix = "",
                          in.path = "/br_z1/kevin_storage/advancedGenomicsCore/tmp/",
                          snps = snps)


sampleQc <- argyle::run.sample.qc(snpRes, max.N = 0.0838047 * 143259)
sampleQc <- argyle::run.marker.qc(sampleQc, max.N = 0.0838047 * 143259)
sampleQc <- argyle::run.sample.qc(snpRes, max.N = 0.0838047 * 143259, ref.intensity = clusters$Rmean)


normedSnpRes <- tQN(sampleQc, clusters = clusters)
baf_normed <- get.baf(normedSnpRes)
baf_mat <- dcast(baf_normed, marker + chr + pos ~ iid, value.var = c("BAF"))
test_intense <- get.intensity(normedSnpRes)
logR_mat <- dcast(test_intense, marker + chr + pos ~ iid, value.var = c("si"))



baf_mat$chr <- stringr::str_remove(baf_mat$chr, "chr")
logR_mat$chr <- stringr::str_remove(logR_mat$chr, "chr")


test_ascat_tumor.logR <- logR_mat
test_ascat_tumor.logR[ , 4:51] <- log2(test_ascat_tumor.logR[ , 4:51])
test_ascat_tumor.baf <- baf_mat

norm_ascat_tumor.logR <- logR_mat[,c("marker", "chr", "pos", "5516-KH-48")]
norm_ascat_tumor.logR$`5516-KH-48` <- log2(norm_ascat_tumor.logR$`5516-KH-48`)
norm_ascat_tumor.baf <- baf_mat[,c("marker", "chr", "pos", "5516-KH-48")]

dupeLogR <- Rfast::rep_col(norm_ascat_tumor.logR$`5516-KH-48`, 47)
norm_ascat_tumor.logR <- cbind(norm_ascat_tumor.logR, dupeLogR)
dupeBaf <- Rfast::rep_col(norm_ascat_tumor.baf$`5516-KH-48`, 47)
norm_ascat_tumor.baf <- cbind(norm_ascat_tumor.baf, dupeBaf)


# test_ascat_tumor.logR2 <- test_ascat_tumor.logR[, c("marker", "chr", "pos", "5516-KH-46", "5516-KH-48")]
# test_ascat_tumor.baf <- baf_mat[, c("marker", "chr", "pos", "5516-KH-46", "5516-KH-48")]
### kc6 and kc10 are the normals

write.table(test_ascat_tumor.logR, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorLogR.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(test_ascat_tumor.baf, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorBaf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(norm_ascat_tumor.logR, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh48NormLogR.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(norm_ascat_tumor.baf, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh48NormBaf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


norm_ascat_tumor.logR <- logR_mat[,c("marker", "chr", "pos", "5516-KH-46")]
norm_ascat_tumor.logR$`5516-KH-46` <- log2(norm_ascat_tumor.logR$`5516-KH-46`)
norm_ascat_tumor.baf <- baf_mat[,c("marker", "chr", "pos", "5516-KH-46")]

dupeLogR <- Rfast::rep_col(norm_ascat_tumor.logR$`5516-KH-46`, 1)
norm_ascat_tumor.logR <- cbind(norm_ascat_tumor.logR, dupeLogR)
dupeBaf <- Rfast::rep_col(norm_ascat_tumor.baf$`5516-KH-46`, 1)
norm_ascat_tumor.baf <- cbind(norm_ascat_tumor.baf, dupeBaf)

write.table(norm_ascat_tumor.logR, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh46NormLogR.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(norm_ascat_tumor.baf, "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh46NormBaf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


# 48

ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorLogR.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorBaf.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh48NormLogR.txt",
                           "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh48NormBaf.txt")


ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/norm2_48/")
ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/norm2_48/")
ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/norm2_48/")
ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/norm2_48/")

save(ascat.bc, file = "/br_z1/kevin_storage/ASCAT/20220404scatbcNorm48.rds")
save(ascat.output, file = "/br_z1/kevin_storage/ASCAT/20220404acsatoutputNorm48.rds")

no_cores <- 8
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl) 
res <- foreach(i = 1:(dim(ascat.bc[[1]])[2]), .combine = "rbind", .packages = c("ASCAT", "ggplot2", "grid")) %dopar% runGridAscat(i, sampleMap, tcDf, "/br_z1/kevin_storage/ASCAT/norm2_48/")
stopCluster(cl)

# 46


# ascat.bc <- ascat.loadData("/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorLogR.txt",
#                            "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404testTumorBaf.txt",
#                            "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh46NormLogR.txt",
#                            "/br_z1/kevin_storage/advancedGenomicsCore/tmp/20220404kh46NormBaf.txt")
# 
# 
# ascat.bc = ascat.correctLogR(ascat.bc, "/br_z1/kevin_storage/ASCAT/GCcontent_SNPloci.txt")
# ascat.plotRawData(ascat.bc, "/br_z1/kevin_storage/ASCAT/norm2_46/")
# ascat.bc = ascat.aspcf(ascat.bc,out.dir = "/br_z1/kevin_storage/ASCAT/norm2_46/")
# ascat.plotSegmentedData(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/norm2_46/")
# ascat.output = ascat.runAscat(ascat.bc, img.dir = "/br_z1/kevin_storage/ASCAT/norm2_46/")
# 





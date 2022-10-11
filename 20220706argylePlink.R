### prepping data for ADMIXTUERE
### for giga muga following paper: https://pubmed.ncbi.nlm.nih.gov/29410804/

library(reshape2)
library(argyle)
load("/br_z1/kevin_storage/misc/snps.gigamuga.Rdata")
load("/br_z1/kevin_storage/misc/clusters.gigamuga.Rdata")

load("/br_z1/kevin_storage/gigamugaFiles/GM_geno.Rdata")
load("/br_z1/kevin_storage/gigamugaFiles/GM_code.Rdata")
load("/br_z1/kevin_storage/gigamugaFiles/GM_sex.Rdata")
load("/br_z1/kevin_storage/gigamugaFiles/GM_x.Rdata")
load("/br_z1/kevin_storage/gigamugaFiles/GM_y.Rdata")


geno <- read.beadstudio(prefix = "", snps, in.path = "/br_z1/kevin_storage/advancedGenomicsCore/tmp/")


ref_map <- attr(geno, "map")
GM_geno_red <- GM_geno[which(rownames(GM_geno) %in% ref_map$marker),]
GM_x_red <- GM_x[which(rownames(GM_x) %in% ref_map$marker),]
GM_y_red <- GM_y[which(rownames(GM_y) %in% ref_map$marker),]
ref_fam <- data.frame(fid = GM_code, iid = GM_code, sex = GM_sex)
ref_geno <- genotypes(GM_geno_red, map = ref_map, alleles = "native",
                      intensity = list(GM_x_red, GM_y_red))


geno2 <- argyle::run.marker.qc(geno)
geno2 <- apply.filters(geno2)


### printing thresholds and names of bad samples

geno3 <- run.sample.qc(geno2, max.N = 0.05 * 143259)
print(paste("0.05", names(which(attr(geno3, "filter.samples") == "N"))))

geno3 <- run.sample.qc(geno2, max.N = 0.04 * 143259)
print(paste("0.04", names(which(attr(geno3, "filter.samples") == "N"))))

geno3 <- run.sample.qc(geno2, max.N = 0.03 * 143259)
print(paste("0.03", names(which(attr(geno3, "filter.samples") == "N"))))

geno3 <- run.sample.qc(geno2, max.N = 0.02 * 143259)
print(paste("0.02", names(which(attr(geno3, "filter.samples") == "N"))))

geno3 <- run.sample.qc(geno2, max.N = 0.01 * 143259)
print(paste("0.01", names(which(attr(geno3, "filter.samples") == "N"))))



write.plink(ref_geno, "/br_z1/kevin_storage/gigamugaFiles/20220705refGeno")
write.plink(geno, "/br_z1/kevin_storage/gigamugaFiles/20220705gigaMuga48Geno")
write.plink(geno2, "/br_z1/kevin_storage/gigamugaFiles/20220705gigaMuga48MarkerFiltGeno")


refBim <- read.table("/br_z1/kevin_storage/gigamugaFiles/test/refsamples.bim", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
refBim_filt <- refBim[-which(refBim$V5 == refBim$V6), ]
autosomes <- as.character(1:19)
refBim_filt <- refBim_filt[which(refBim_filt$V1 %in% autosomes), ]
write.table(refBim_filt, "/br_z1/kevin_storage/gigamugaFiles/test/refsamplesFilt.bim", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


ref_geno <- read.plink("/br_z1/kevin_storage/gigamugaFiles/test/refsamples", markers = refBim_filt$V2)
ref_geno_native <- recode(ref_geno, "native")
write.plink(ref_geno_native, "/br_z1/kevin_storage/gigamugaFiles/test/refsamplesAutosomes")

geno48 <- read.plink("/br_z1/kevin_storage/gigamugaFiles/test/20220705gigaMuga48Geno", markers = refBim_filt$V2)
geno48_native <- recode(geno48, "native")
write.plink(geno48_native, "/br_z1/kevin_storage/gigamugaFiles/test/geno48Autosomes")


### can further filter which reference samples I use by getting geno , ped and map matrix and filter by samples accordingly
ref_geno_autosomes <- read.plink("/br_z1/kevin_storage/gigamugaFiles/test/refsamplesAutosomes")
ref_geno2 <- get.call(gty = ref_geno_autosomes)
ref_geno2_ped <- attr(ref_geno_autosomes, "ped")
ref_geno2_map <- attr(ref_geno_autosomes, "map")

### filter by fam name grep to BALB, C57BL ... etc

interestedStrains <- c("129", "C57", "BALB", "C3H", "FVB")
ref_geno2_ped_filt <- ref_geno2_ped[grep(paste0(interestedStrains, collapse = "|"), ref_geno2_ped$iid),]
ref_geno2_ped_filt <- ref_geno2_ped_filt[which(ref_geno2_ped_filt$fid == "NILG"), ]

ref_geno2_filt <- ref_geno2[which(ref_geno2$iid %in% ref_geno2_ped_filt$iid),]
ref_geno2_filt$iid <- as.character(ref_geno2_filt$iid)
ref_geno2_filt$marker <- as.character(ref_geno2_filt$marker)

ref_geno2_filt_tall <- as.matrix(dcast(ref_geno2_filt, iid ~ marker, value.var = "call"), nrow = 37, ncol = 137476)
ref_geno2_filt_tall <- ref_geno2_filt_tall[, 2:ncol(ref_geno2_filt_tall)]
ref_geno2_filt_tall <- apply(ref_geno2_filt_tall, 2, as.numeric)
dimnames(ref_geno2_filt_tall) <- list(unique(ref_geno2_ped_filt$iid), unique(ref_geno2_filt$marker))
ref_geno2_filt_tall <- ref_geno2_filt_tall[, match(ref_geno2_map$marker, colnames(ref_geno2_filt_tall))]
ref_geno2_filt_tall <- t(ref_geno2_filt_tall)

ref_geno3 <- genotypes(G = ref_geno2_filt_tall, map = ref_geno2_map, ped = ref_geno2_ped_filt, alleles = "01")
ref_geno3_native <- recode(ref_geno3, "native")
write.plink(ref_geno3_native, "/br_z1/kevin_storage/gigamugaFiles/test/refSamplesAuto5strains")

# colnames(ref_geno2_filt_tall)[which(colnames(ref_geno2_filt_tall) %in% ref_geno2_map$marker)]
# ref_geno2_map$marker[which(ref_geno2_map$marker %in% colnames(ref_geno2_filt_tall))]

save.image(file = "/br_z1/kevin_storage/gigamugaFiles/processingPlankFiles.RData")


### trying other ways of calling. beause using multiple strains may mess it up i.e multiple C57Bl as a founder C57Bl

load("/br_z1/kevin_storage/gigamugaFiles/processingPlankFiles.RData")


interestedStrains <- c("129X1/SvJ", "C57BL/6NJm39418", "BALB/cJ", "C3H/HeJm39336", "FVB001")
ref_geno2_ped_filt <- ref_geno2_ped[grep(paste0(interestedStrains, collapse = "$|"), ref_geno2_ped$iid),]
ref_geno2_ped_filt <- ref_geno2_ped_filt[which(ref_geno2_ped_filt$fid == "NILG"), ]

ref_geno2_filt <- ref_geno2[which(ref_geno2$iid %in% ref_geno2_ped_filt$iid),]
ref_geno2_filt$iid <- as.character(ref_geno2_filt$iid)
ref_geno2_filt$marker <- as.character(ref_geno2_filt$marker)

ref_geno2_filt_tall <- as.matrix(dcast(ref_geno2_filt, iid ~ marker, value.var = "call"), nrow = 37, ncol = 137476)
ref_geno2_filt_tall <- ref_geno2_filt_tall[, 2:ncol(ref_geno2_filt_tall)]
ref_geno2_filt_tall <- apply(ref_geno2_filt_tall, 2, as.numeric)
dimnames(ref_geno2_filt_tall) <- list(unique(ref_geno2_ped_filt$iid), unique(ref_geno2_filt$marker))
ref_geno2_filt_tall <- ref_geno2_filt_tall[, match(ref_geno2_map$marker, colnames(ref_geno2_filt_tall))]
ref_geno2_filt_tall <- t(ref_geno2_filt_tall)

ref_geno3 <- genotypes(G = ref_geno2_filt_tall, map = ref_geno2_map, ped = ref_geno2_ped_filt, alleles = "01")
ref_geno3_native <- recode(ref_geno3, "native")
write.plink(ref_geno3_native, "/br_z1/kevin_storage/gigamugaFiles/test/20220912refSamplesAuto5strains")



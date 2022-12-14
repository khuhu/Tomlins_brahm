ascat.predictGermlineGenotypes <- function (ASCATobj, platform = "AffySNP6", img.dir = ".", img.prefix = "") 
{
  Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
  rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)
  if (platform == "Custom10k") {
    maxHomozygous = 0.05
    proportionHetero = 0.59
    proportionHomo = 0.38
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform == "IlluminaASA") {
    maxHomozygous = 0.05
    proportionHetero = 0.15
    proportionHomo = 0.82
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform == "IlluminaGSAv3") {
    maxHomozygous = 0.05
    proportionHetero = 0.16
    proportionHomo = 0.8
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform == "Illumina109k") {
    maxHomozygous = 0.05
    proportionHetero = 0.35
    proportionHomo = 0.6
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform == "IlluminaCytoSNP") {
    maxHomozygous = 0.05
    proportionHetero = 0.28
    proportionHomo = 0.62
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform == "IlluminaCytoSNP850k") {
    maxHomozygous = 0.05
    proportionHetero = 0.23
    proportionHomo = 0.72
    proportionOpen = 0.01
    segmentLength = 60
  }
  else if (platform == "Illumina610k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform == "Illumina660k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform == "Illumina700k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform == "Illumina1M") {
    maxHomozygous = 0.05
    proportionHetero = 0.22
    proportionHomo = 0.74
    proportionOpen = 0.02
    segmentLength = 100
  }
  else if (platform == "Illumina2.5M") {
    maxHomozygous = 0.05
    proportionHetero = 0.21
    proportionHomo = 0.745
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform == "IlluminaOmni5") {
    maxHomozygous = 0.05
    proportionHetero = 0.13
    proportionHomo = 0.855
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform == "Affy10k") {
    maxHomozygous = 0.04
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 20
  }
  else if (platform == "Affy100k") {
    maxHomozygous = 0.05
    proportionHetero = 0.27
    proportionHomo = 0.62
    proportionOpen = 0.09
    segmentLength = 30
  }
  else if (platform == "Affy250k_sty") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform == "Affy250k_nsp") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform == "AffySNP6") {
    maxHomozygous = 0.05
    proportionHetero = 0.25
    proportionHomo = 0.67
    proportionOpen = 0.04
    segmentLength = 100
  }
  else if (platform == "AffyOncoScan") {
    maxHomozygous = 0.04
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 30
  }
  else if (platform == "AffyCytoScanHD") {
    maxHomozygous = 0.04
    proportionHetero = 0.32
    proportionHomo = 0.6
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform == "HumanCNV370quad") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 20
  }
  else if (platform == "HumanCore12") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 20
  }
  else if (platform == "HumanCoreExome24") {
    maxHomozygous = 0.05
    proportionHetero = 0.175
    proportionHomo = 0.79
    proportionOpen = 0.02
    segmentLength = 100
  }
  else if (platform == "HumanOmniExpress12") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 100
  }
  else if (platform == "IlluminaOmniExpressExome") {
    maxHomozygous = 0.05
    proportionHetero = 0.35
    proportionHomo = 0.6
    proportionOpen = 0.03
    segmentLength = 100
  } else if(platform == "IonTorrent"){
    maxHomozygous = 0.24
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 50
  } else if(platform == "IonTorrentRes05"){
    maxHomozygous = 0.05
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 50
  } else if(platform == "IonTorrentRes05_2"){
    maxHomozygous = 0.05
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 50
  }else if(platform == "IonTorrentRes10"){
    maxHomozygous = 0.10
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 50
  } else if(platform == "IonTorrentRes15"){
    maxHomozygous = 0.15
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 50
  } else if(platform == "IonTorrentRes05Seg100"){
    maxHomozygous = 0.05
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 100
  } else if(platform == "IonTorrentRes10Seg100"){
    maxHomozygous = 0.10
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 100
  } else if(platform == "IonTorrentRes15Seg100"){
    maxHomozygous = 0.15
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 100
  } else if(platform == "IonTorrentRes15Seg30"){
    maxHomozygous = 0.15
    proportionHetero = 0.11
    proportionHomo = 0.86
    proportionOpen = 0.04
    segmentLength = 30
  }else {
    print("Error: platform unknown")
  }

  failedarrays = NULL
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[, 
                                                                  i]), i]
    names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[, 
                                                                                   i])]
    Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA), 
                                          i]
    names(Tumor_LogR_noNA) = names(Tumor_BAF_noNA)
    chr_noNA = list()
    prev = 0
    for (j in 1:length(ASCATobj$chr)) {
      chrke = ASCATobj$chr[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke, 
                                                   i]))
      chr_noNA[[j]] = (prev + 1):next2
      prev = next2
    }
    ch_noNA = list()
    prev = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrke = ASCATobj$ch[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke, 
                                                   i]))
      ch_noNA[[j]] = (prev + 1):next2
      prev = next2
    }
    tbsam = Tumor_BAF_noNA
    bsm = ifelse(tbsam < 0.5, tbsam, 1 - tbsam)
    homoLimit = max(sort(bsm)[round(length(bsm) * proportionHomo)], 
                    maxHomozygous)
    if (homoLimit > 0.25) {
      failedarrays = c(failedarrays, ASCATobj$samples[i])
    }
    Hom = ifelse(bsm < homoLimit, T, NA)
    Homo = sum(Hom == T, na.rm = T)
    Undecided = sum(is.na(Hom))
    extraHetero = round(min(proportionHetero * length(Tumor_BAF_noNA), 
                            Undecided - proportionOpen * length(Tumor_BAF_noNA)))
    Hetero = 0
    if (extraHetero > 0) {
      allProbes = 1:length(Tumor_BAF_noNA)
      nonHomoProbes = allProbes[is.na(Hom) | Hom == F]
      lowestDist = NULL
      bsmHNA = bsm
      bsmHNA[!is.na(Hom) & Hom] = NA
      for (chrke in chr_noNA) {
        chrNonHomoProbes = intersect(nonHomoProbes, chrke)
        if (length(chrNonHomoProbes) > 5) {
          segmentLength2 = min(length(chrNonHomoProbes) - 
                                 1, segmentLength)
          chrNonHomoProbesStartWindowLeft = c(rep(NA, 
                                                  segmentLength2), chrNonHomoProbes[1:(length(chrNonHomoProbes) - 
                                                                                         segmentLength2)])
          chrNonHomoProbesEndWindowLeft = c(NA, chrNonHomoProbes[1:(length(chrNonHomoProbes) - 
                                                                      1)])
          chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)], 
                                               NA)
          chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2 + 
                                                                 1):length(chrNonHomoProbes)], rep(NA, segmentLength2))
          chrNonHomoProbesStartWindowMiddle = c(rep(NA, 
                                                    segmentLength2/2), chrNonHomoProbes[1:(length(chrNonHomoProbes) - 
                                                                                             segmentLength2/2)])
          chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2 + 
                                                                  1):length(chrNonHomoProbes)], rep(NA, segmentLength2/2))
          chrLowestDist = NULL
          for (probeNr in 1:length(chrNonHomoProbes)) {
            probe = chrNonHomoProbes[probeNr]
            if (!is.na(chrNonHomoProbesStartWindowLeft[probeNr]) & 
                !is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
              medianLeft = median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], 
                                  na.rm = T)
            }
            else {
              medianLeft = NA
            }
            if (!is.na(chrNonHomoProbesStartWindowRight[probeNr]) & 
                !is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
              medianRight = median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]], 
                                   na.rm = T)
            }
            else {
              medianRight = NA
            }
            if (!is.na(chrNonHomoProbesStartWindowMiddle[probeNr]) & 
                !is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
              medianMiddle = median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], 
                                      bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]), 
                                    na.rm = T)
            }
            else {
              medianMiddle = NA
            }
            chrLowestDist[probeNr] = min(abs(medianLeft - 
                                               bsm[probe]), abs(medianRight - bsm[probe]), 
                                         abs(medianMiddle - bsm[probe]), Inf, na.rm = T)
          }
        }
        else {
          chrLowestDist = NULL
          if (length(chrNonHomoProbes) > 0) {
            chrLowestDist[1:length(chrNonHomoProbes)] = 1
          }
        }
        lowestDist = c(lowestDist, chrLowestDist)
      }
      lowestDistUndecided = lowestDist[is.na(Hom[nonHomoProbes])]
      names(lowestDistUndecided) = names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]
      sorted = sort(lowestDistUndecided)
      Hom[names(sorted[1:min(length(sorted), extraHetero)])] = F
      Hetero = sum(Hom == F, na.rm = T)
      Homo = sum(Hom == T, na.rm = T)
      Undecided = sum(is.na(Hom))
    }
    png(filename = file.path(img.dir, paste(img.prefix, "tumorSep", 
                                            colnames(ASCATobj$Tumor_LogR)[i], ".png", sep = "")), 
        width = 2000, height = 500, res = 200)
    title = paste(paste(colnames(ASCATobj$Tumor_BAF)[i], 
                        Hetero), Homo)
    ascat.plotGenotypes(ASCATobj, title, Tumor_BAF_noNA, 
                        Hom, ch_noNA)
    dev.off()
    Hom[is.na(Hom)] = T
    Homozygous[names(Hom), i] = Hom
  }
  return(list(germlinegenotypes = Homozygous, failedarrays = failedarrays))
}

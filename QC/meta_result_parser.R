library(data.table)

# test destination
# '../../trans_eQTL_meta_analysis/Script_for_QC_05012015/eQTLProbesFDR0.05-ProbeLevel_21_datasets.txt'

parseMetaResult <- function(x){

  require(data.table)
  require(stringr)
  require(dplyr)
  
  andRaw <- fread(x)
  
  zScores <- andRaw$DatasetsZScores
  zScores <- t(as.data.frame(str_split(zScores, ';')))
  rownames(zScores) <- paste(andRaw$SNPName, andRaw$ProbeName)
  Names <- t(as.data.frame(str_split(andRaw$DatasetsWhereSNPProbePairIsAvailableAndPassesQC, ';')))
  
  uniqueColNames <- function(y){unique(y[!y == '-'])}
  Names <- apply(Names, 2, uniqueColNames)
  colnames(zScores) <- Names
  zScores[zScores == "-"] <- NA
  zScores <- apply(zScores, 2, as.numeric)
  rownames(zScores) <- paste(andRaw$SNPName, andRaw$ProbeName)
  # parse samples sizes
  Samples <- t(as.data.frame(str_split(andRaw$DatasetsNrSamples, ';')))
  
  Samples[Samples == '-'] <- NA
  
  Samples <- apply(Samples, 2, as.numeric)
  colnames(Samples) <- Names
  rownames(Samples) <- paste(andRaw$SNPName, andRaw$ProbeName)
  
  metaZ <- andRaw$OverallZScore

  annotation <- data.frame(SNPName_meta = andRaw$SNPName, SNPChr_meta = andRaw$SNPChr, SNPChrPos_meta = andRaw$SNPChrPos, ProbeChr_meta = andRaw$ProbeChr, ProbeCenterPos_meta = andRaw$ProbeCenterChrPos, SNPType_meta = andRaw$SNPType, AlleleAssessed_meta = andRaw$AlleleAssessed)
  
  return(list(annotation = annotation, zScoreMatrix = zScores, sampleSizes = Samples, metaZ = metaZ))
  
}

#meta <- parseMetaResult('../../trans_eQTL_meta_analysis/Script_for_QC_05012015/eQTLProbesFDR0.05-ProbeLevel_21_datasets.txt')

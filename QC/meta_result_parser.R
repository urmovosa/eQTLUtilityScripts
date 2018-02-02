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
  
  Sum_samples <- rowSums(Samples, na.rm = T)
  Sum_samples <- data.frame(eQTL = names(Sum_samples), N = Sum_samples)
  SNP <- str_replace(Sum_samples$eQTL, ' .*', '')
  Probe <- str_replace(Sum_samples$eQTL, '.* ', '')
  
  Sum_samples <- data.frame(SNP = SNP, Probe = Probe, N = Sum_samples$N)
  
  # Find in how many datasets the eQTL was tested:
  
  cohort_count <- apply(Samples, 1, function(x) sum(!is.na(x)))
  cohort_count <- data.frame(eQTL = names(cohort_count), N = cohort_count)
  SNP <- str_replace(cohort_count$eQTL, ' .*', '')
  Probe <- str_replace(cohort_count$eQTL, '.* ', '')
  
  cohort_count <- data.frame(SNP = SNP, Probe = Probe, N = cohort_count$N)
  
  metaZ <- andRaw$OverallZScore
  P_value <- andRaw$PValue
  
  
  annotation <- data.frame(SNPName_meta = andRaw$SNPName, SNPChr_meta = andRaw$SNPChr, SNPChrPos_meta = andRaw$SNPChrPos, ProbeName_meta = andRaw$ProbeName, ProbeChr_meta = andRaw$ProbeChr, ProbeCenterPos_meta = andRaw$ProbeCenterChrPos, SNPType_meta = andRaw$SNPType, AlleleAssessed_meta = andRaw$AlleleAssessed)
  
  return(list(annotation = annotation, 
              zScoreMatrix = zScores, 
              sampleSizes = Samples, 
              sumSampleSize = Sum_samples,
              nrOfCohorts = cohort_count,
              metaZ = metaZ,
              uncorrected_P = P_value))
  
}

#meta <- parseMetaResult('../../trans_eQTL_meta_analysis/Script_for_QC_05012015/eQTLProbesFDR0.05-ProbeLevel_21_datasets.txt')

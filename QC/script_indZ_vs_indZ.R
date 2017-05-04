
if(!exists("parseMetaResult", mode = "function")) source("meta_result_parser.R")

library(limma)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(rtracklayer)
library(stringr)
library(data.table)

##########################
# 1. Visualize Z-scores ##
##########################

## read in the meta-analysis summary

and <- parseMetaResult('../eQTLsFDR0.05-ProbeLevel.txt')


allSetZ <- as.data.frame(and$zScoreMatrix)

# NB! change according the final list of studies!
#names[c(12, 13)] <- c('SIGN', 'Rotterdam')

allSetZ$metaZ <- as.numeric(and$metaZ)
# NB! change according the final list of studies!
allSetZ2 <- gather(allSetZ, "study", "Z_score", 1:21)
allSetZ2$concordance <- 'concordant'


allSetZ2[allSetZ2$Z_score < 0 & allSetZ2$metaZ > 0 | 
           allSetZ2$Z_score > 0 & allSetZ2$metaZ < 0 | 
           is.na(allSetZ2$metaZ) | 
           is.na(allSetZ2$Z_score), ]$concordance <- 'unconcordant'

# pairwise Z-score comparison

library(GGally)

png('z_scores_pairwise_interim.png', width = 40, height = 40, units = 'in', res = 300)
p <- ggpairs(allSetZ,
             lower = list(continuous = wrap("points", alpha = 0.1))) + theme_bw()

#p <- ggpairs(allSetZPair2, alpha = 0.1, params = c(fill = "white", color = "black"), lower = list(continuous = "points", combo = "dot")) + theme_bw()
p

dev.off()

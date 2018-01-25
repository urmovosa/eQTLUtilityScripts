
if(!exists("parseMetaResult", mode = "function")) source("meta_result_parser.R")

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

and <- parseMetaResult('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/trans_PRS_meta_analysis_20180125/eQTLsFDR-Significant-0.05.txt')


allSetZ <- as.data.frame(and$zScoreMatrix)

## Replace the names with newest versions:

# Use mapping file: 1. col old dataset name (name_old), 2. column new standardized name (name_new)
# Order the names so that they should emerge on the plot
name_mapping <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/trans_PRS_meta_analysis_20180125/name_mapping.txt')

# works only where there is some eQTL which is present in all datasets

names_orig <- colnames(allSetZ)
names_orig <- unlist(str_split(names_orig, ';'))

names <- data.frame(name_old = names_orig)
names <- merge(names, name_mapping, by = 'name_old', all.x = T)
names$name_new <- as.character(names$name_new)
names$name_old <- as.character(names$name_old)

names <- names[match(names_orig, names$name_old), ]

if (length(names$name_new[is.na(names$name_new)] > 0)){print(paste0('Issue with old name: ', as.character(names[is.na(names$name_new), ]$name_old), '. Check your name mapping file!'))}

####
colnames(allSetZ) <- names$name_new

allSetZ$metaZ <- as.numeric(and$metaZ)



allSetZ2 <- gather(allSetZ, "study", "Z_score", 1 : ncol(and$sampleSizes))
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

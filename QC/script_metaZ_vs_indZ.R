######################################################################################################################################
# This script takes eQTLMappingPipeline output and does basic QC: compares each individual Z-score against the meta-analysed Z-score #
######################################################################################################################################

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

## read in the meta-analysis summary

and <- fread('../Interpretation/eQTLsFDR0.05-ProbeLevel.txt')


allSetZ <- and$DatasetsZScores
allSetZ <- read.table(text = allSetZ, sep = ';', na.strings = '-')

# works only where there are some eQTL which is present in all datasets
names <- and$DatasetsWhereSNPProbePairIsAvailableAndPassesQC[!str_detect(and$DatasetsWhereSNPProbePairIsAvailableAndPassesQC, '-')][1]
names <- unlist(str_split(names, ';'))

N <- and$DatasetsNrSamples[!str_detect(and$DatasetsNrSamples, '-')][1]
N <- unlist(str_split(N, ';'))

colnames(allSetZ) <- names
allSetZ$eQTL <- paste(and$SNPName, and$ProbeName, sep = '_')
allSetZ$metaZ <- and$OverallZScore

abi1 <- data.frame(study = names, N = N)

allSetZ2 <- gather(allSetZ, "study", "Z_score", 1:(ncol(allSetZ)-2))
allSetZ2$concordance <- 'concordant'


allSetZ2[allSetZ2$Z_score < 0 & allSetZ2$metaZ > 0 | 
           allSetZ2$Z_score > 0 & allSetZ2$metaZ < 0 | 
           is.na(allSetZ2$metaZ) | 
           is.na(allSetZ2$Z_score), ]$concordance <- 'unconcordant'

allSetZ2 <- allSetZ2[!is.na(allSetZ2$Z_score), ]

# calculate % of concordance per dataset

allSetZ2 <- allSetZ2 %>%
  group_by(study) %>%
  mutate(nr_of_assoc = n(), nr_of_concordant = length(study[concordance == 'concordant']), nr_of_unconcordant = length(study[concordance == 'unconcordant']), 
         perc_of_concordant = round((length(study[concordance == 'concordant'])/n()) * 100, digits = 1), y_pos = (max(Z_score)))

abi2 <- unique(allSetZ2[, c(3, 6, 7, 8, 9, 10)])

ann_text <- data.frame(metaZ = (min(allSetZ2$metaZ) + 0.1 * min(allSetZ2$metaZ)), Z_score = abi2$y_pos,
                       lab = paste("nr. of assoc: ", abi2$nr_of_assoc, "\nconcordant: ", abi2$nr_of_concordant, "\nunconcordant: ", abi2$nr_of_unconcordant, "\nconcordant: ", abi2$perc_of_concordant, '%', sep = ''),
                       study = factor(unique(allSetZ2$study)), concordance = factor('concordant', levels = c('concordant', 'unconcordant')))

ann_text <- merge(ann_text, abi1, by = "study")
ann_text$lab <- paste('N: ', ann_text$N, '\n', ann_text$lab, sep = '')

# visualize

p <- ggplot(allSetZ2, aes(x = metaZ, y = Z_score, colour = concordance)) + 
  geom_point(alpha = 0.2) + 
  theme_bw() + 
  facet_wrap(~study, scales = 'free_y') + 
  geom_hline(aes(yintercept = 0), linetype = "longdash", colour = 'forestgreen') + 
  geom_vline(aes(xintercept = 0), linetype = "longdash", colour = 'forestgreen') + 
  scale_color_manual(values = c('black', 'red')) + 
  ylab('Study Z-score') + 
  xlab('Meta-analysis Z-score') + geom_text(data = ann_text, aes(label = lab), size = 2, hjust=0, vjust = 1) + theme(legend.position = "none")
#p

ggsave(p, filename = 'z_scores_interim.png', width = 8, height = 8)

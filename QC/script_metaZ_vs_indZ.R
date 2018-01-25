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

and <- fread('eQTLsFDR-Significant-0.05.txt')

## Replace the names with newest versions:

# Use mapping file: 1. col old dataset name (name_old), 2. column new standardized name (name_new)
# Order the names so that they should emerge on the plot
name_mapping <- fread('name_mapping.txt')

# works only where there is some eQTL which is present in all datasets

names_orig <- and$DatasetsWhereSNPProbePairIsAvailableAndPassesQC[!str_detect(and$DatasetsWhereSNPProbePairIsAvailableAndPassesQC, ';-;')][1]
names_orig <- unlist(str_split(names_orig, ';'))

names <- data.frame(name_old = names_orig)
names <- merge(names, name_mapping, by = 'name_old', all.x = T)
names$name_new <- as.character(names$name_new)
names$name_old <- as.character(names$name_old)

names <- names[match(names_orig, names$name_old), ]

if (length(names$name_new[is.na(names$name_new)] > 0)){print(paste0('Issue with old name: ', as.character(names[is.na(names$name_new), ]$name_old), '. Check your name mapping file!'))}

and2 <- and[, c(2, 5, 13), with = F] %>%
  separate(DatasetsZScores, as.character(names$name_old), ";")

allSetZ <- and2[, -c(1, 2), with = F]

N <- and$DatasetsNrSamples[!str_detect(and$DatasetsNrSamples, '-')][1]
N <- unlist(str_split(N, ';'))

colnames(allSetZ) <- names$name_new
allSetZ$eQTL <- paste(and$SNPName, and$ProbeName, sep = '_')
allSetZ$metaZ <- and$OverallZScore

abi1 <- data.frame(study = names$name_new, N = N)

allSetZ2 <- gather(allSetZ, "study", "Z_score", 1:(ncol(allSetZ) - 2))
allSetZ2$concordance <- 'concordant'

#allSetZ2[allSetZ2$metaZ == '-', ]$metaZ <- NA
allSetZ2[allSetZ2$Z_score == '-', ]$Z_score <- NA

allSetZ2$metaZ <- as.numeric(allSetZ2$metaZ)
allSetZ2$Z_score <- as.numeric(allSetZ2$Z_score)

allSetZ2[allSetZ2$Z_score < 0 & allSetZ2$metaZ > 0 | 
           allSetZ2$Z_score > 0 & allSetZ2$metaZ < 0 | 
           is.na(allSetZ2$metaZ) | 
           is.na(allSetZ2$Z_score), ]$concordance <- 'unconcordant'

allSetZ2 <- allSetZ2[!is.na(allSetZ2$Z_score), ]

# calculate % of concordance per dataset

allSetZ2 <- allSetZ2 %>%
  group_by(study) %>%
  mutate(nr_of_assoc = n(), perc_overlapping_with_meta = round((n()/nrow(and)) * 100, digits = 1), nr_of_concordant = length(study[concordance == 'concordant']), nr_of_unconcordant = length(study[concordance == 'unconcordant']), 
         perc_of_concordant = round((length(study[concordance == 'concordant'])/n()) * 100, digits = 1), y_pos = (max(Z_score)))

abi2 <- unique(allSetZ2[, c(3, 6, 7, 8, 9, 10, 11)])

ann_text <- data.frame(metaZ = (min(allSetZ2$metaZ) + 0.1 * min(allSetZ2$metaZ)), Z_score = abi2$y_pos,
                       lab = paste("nr. of assoc: ", abi2$nr_of_assoc, "\ntested in dataset: ", abi2$perc_overlapping_with_meta, '%', "\nconcordant: ", abi2$nr_of_concordant, "\nunconcordant: ", abi2$nr_of_unconcordant, "\nconcordant: ", abi2$perc_of_concordant, '%', sep = ''),
                       study = factor(unique(allSetZ2$study)), concordance = factor('concordant', levels = c('concordant', 'unconcordant')))

ann_text <- merge(ann_text, abi1, by = "study")
ann_text$lab <- paste('N: ', ann_text$N, '\n', ann_text$lab, sep = '')

min_Z <- min(as.numeric(allSetZ2$metaZ), as.numeric(allSetZ2$Z_score))
max_Z <- max(as.numeric(allSetZ2$metaZ), as.numeric(allSetZ2$Z_score))

ann_text$study <- as.character(ann_text$study)

# Order the studies by the platform and name:
allSetZ2 <- allSetZ2 %>% ungroup()

allSetZ2$study <- factor(allSetZ2$study, 
                         levels = as.character(name_mapping$name_new))

ann_text$study <- factor(ann_text$study, 
                         levels = as.character(name_mapping$name_new))

# visualize

p <- ggplot(allSetZ2, aes(x = metaZ, y = Z_score, colour = concordance)) + 
  geom_point(alpha = 0.2) + 
  theme_bw() + 
  facet_wrap(~ study, scales = 'free_y') + 
  geom_hline(aes(yintercept = 0), linetype = "longdash", colour = 'forestgreen') + 
  geom_vline(aes(xintercept = 0), linetype = "longdash", colour = 'forestgreen') + 
  scale_color_manual(values = c('black', 'red')) + 
  ylab('Study Z-score') + 
  xlab('Meta-analysis Z-score') + 
  geom_text(data = ann_text, aes(label = lab, x = metaZ, y = Z_score), size = 2, hjust = 0, vjust = 1, inherit.aes = FALSE) + 
  theme(legend.position = "none")

ggsave(p, filename = 'Compare_meta_Z_vs_cohort_Z_', Sys.Date(),'.png', width = 10 * 1.7 * 1.4, height = 10 * 1.5 * 1.2)

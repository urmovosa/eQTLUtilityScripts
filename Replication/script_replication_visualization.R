#############################################################################################################
# This script uses discovery meta-analysis output, replication output files and constructs replication plot #
#############################################################################################################

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

discovery <- args[1]
replication <- args[2]
replication_sig <- args[3]
plotname <- args[4]

## read in the discovery meta-analysis results:

discovery <- fread(args[1])

##  replication:

replication <- fread(args[2])
replication_sig <- fread(args[3])

# merge together 

discovery$pair <- paste(discovery$SNPName, discovery$ProbeName, sep = '_')
discovery <- discovery[, c(22, 10, 11), with = F]
colnames(discovery) <- paste0('Discovery_', colnames(discovery))

colnames(replication)[c(10, 11, 22)] <- c('AlleleAssessed', 'OverallZScore', 'pair')
replication <- replication[, c(22, 10, 11), with = F]
colnames(replication) <- paste0('Replication_', colnames(replication))

colnames(replication_sig)[c(10, 11, 22)] <-c('AlleleAssessed', 'OverallZScore', 'pair')

replication$sig_replication <- 'no'
replication[replication$Replication_pair %in% replication_sig$pair, ]$sig_replication <- 'yes'

all_merged <- merge(discovery, replication, by.x = 'Discovery_pair', by.y = 'Replication_pair', all.x = T)
all_merged[!all_merged$Discovery_AlleleAssessed == all_merged$Replication_AlleleAssessed, ]$Replication_OverallZScore <- -all_merged[!all_merged$Discovery_AlleleAssessed == all_merged$Replication_AlleleAssessed, ]$Replication_OverallZScore
all_merged$same_dir <- 'no'

all_merged[(all_merged$Discovery_OverallZScore > 0 & all_merged$Replication_OverallZScore > 0) |
          (all_merged$Discovery_OverallZScore < 0 & all_merged$Replication_OverallZScore < 0) ]$same_dir <- 'yes'

all_merged$replicated <- 'no'
all_merged[all_merged$sig_replication == 'yes' & all_merged$same_dir == 'yes', ]$replicated <- 'yes'


# calculate statistics:

statistics_overview <- all_merged %>%
  summarize(nr_of_disc_assoc = length(Replication_OverallZScore), 
         nr_of_tested_in_replication = sum(!is.na(Replication_OverallZScore)),
         same_direction = sum(same_dir == 'yes'),
         nr_of_sig_in_replication = sum(!is.na(Replication_OverallZScore) & sig_replication == 'yes'),
         nr_of_replicated_same_dir = sum(!is.na(Replication_OverallZScore) & sig_replication == 'yes' & same_dir == 'yes')) %>%
  mutate(tested_eqtl_perc = round(nr_of_tested_in_replication/nr_of_disc_assoc, digits = 4) * 100,
         same_dir_perc = round(same_direction/nr_of_tested_in_replication, digits = 4) * 100,
         replication_perc = round(nr_of_replicated_same_dir/nr_of_tested_in_replication, digits = 4) * 100)
         

ann_text <- data.frame(Discovery_OverallZScore = (min(all_merged$Discovery_OverallZScore, na.rm = T) + 0.1 * min(all_merged$Discovery_OverallZScore, na.rm = T)), 
                       Replication_OverallZScore = max(all_merged$Replication_OverallZScore, na.rm = T),
                       lab = paste0("Nr. of discovery assoc: ", statistics_overview$nr_of_disc_assoc, 
                                   "\nTested in replication cohort: ", statistics_overview$nr_of_tested_in_replication, ' (', statistics_overview$tested_eqtl_perc, '%)',
                                   "\nSame dir. in replication cohort: ", statistics_overview$same_direction, ' (', statistics_overview$same_dir_perc, '%)', 
                                    "\nReplicated in replication cohort: ", statistics_overview$nr_of_replicated_same_dir, ' (', statistics_overview$replication_perc, '%)'),
                       `Replication\nFDR<0.05` = factor('yes', levels = c('yes', 'no')))
colnames(ann_text)[ncol(ann_text)] <- "Replication\nFDR<0.05"

###

colnames(all_merged)[8] <- 'Replication\nFDR<0.05'
all_merged <- all_merged[!is.na(all_merged$`Replication\nFDR<0.05`), ]

all_merged$`Replication\nFDR<0.05` <- factor(all_merged$`Replication\nFDR<0.05`, 
                                            levels = c('yes', 'no'))

p <- ggplot(all_merged, aes(x = Discovery_OverallZScore, y = Replication_OverallZScore, colour = `Replication\nFDR<0.05`)) + 
  geom_point(alpha = 0.35) + 
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) + 
  scale_color_manual(values = c("yes" = "red", "no" = "darkgrey")) + 
  xlab('Discovery analysis Z-score') + 
  ylab('Replication analysis Z-score') +
  geom_text(data = ann_text, aes(label = lab), size = 3, hjust=0, vjust = 1, col = 'black')

#scale_y_continuous(limits = c(min(all_merged$Discovery_OverallZScore, all_merged$Replication_OverallZScore, na.rm = T), max(all_merged$Discovery_OverallZScore, all_merged$Replication_OverallZScore, na.rm = T))) + 

ggsave(paste0(arg[4], '.png'), height = 7 * 0.9, width = 8 * 0.9, dpi = 600, units = 'in')

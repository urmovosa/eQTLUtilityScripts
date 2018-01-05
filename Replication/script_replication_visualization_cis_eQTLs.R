#############################################################################################################
# This script uses discovery meta-analysis output, replication output files and constructs replication plot #
#############################################################################################################

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

## read in the discovery meta-analysis results:

discovery <- fread(args[1])

##  replication:

replication <- fread(args[3])
replication_sig <- fread(args[2])

# calculate the sample sizes for discovery and replication


Ndisc <- t(as.data.frame(str_split(unique(discovery$DatasetsNrSamples), ';')))
rownames(Ndisc) <- paste(discovery$SNPName, discovery$ProbeName)
Ndisc[Ndisc == "-"] <- NA

if (class(Ndisc) == 'matrix' & ncol(Ndisc) > 1){
Ndisc <- apply(Ndisc, 2, as.numeric)
Ndisc <- Ndisc[complete.cases(Ndisc), ]
Ndisc <- max(rowSums(Ndisc))
} else {Ndisc <- unique(as.numeric(Ndisc))}


#if(class(Ndisc) == 'vector' & length(Ndisc) < 1){stop("None of the discovery eQTLs is present in all of the cohorts!")}

# if(class(Ndisc) == 'matrix' & ncol(Ndisc) > 1){
# Ndisc <- max(rowSums(Ndisc))
# } else {Ndisc <- unique(Ndisc)}

colnames(replication)[c(2, 5, 14)] <- c('SNPName', 'ProbeName', 'DatasetsNrSamples')

Nrep <- t(as.data.frame(str_split(unique(replication$DatasetsNrSamples), ';')))
dim(Nrep)
dim(replication)
rownames(Nrep) <- paste(replication$SNPName, replication$ProbeName)
Nrep[Nrep == "-"] <- NA


if (class(Nrep) == 'matrix' & ncol(Nrep) > 1){
Nrep <- apply(Nrep, 2, as.numeric)
Nrep <- Nrep[complete.cases(Nrep), ]
Nrep <- max(rowSums(Nrep))
} else {Nrep <- unique(as.numeric(Nrep))}

print(Nrep)

#if(class(Nrep) == 'vector' & length(Nrep) < 1){stop("None of the replication eQTLs is present in all of the cohorts!")}

# if(class(Nrep) == 'matrix'){
#   Nrep <- max(rowSums(Nrep))
# } else {Nrep <- unique(Nrep)}

# merge together 
dim(discovery)
discovery$pair <- paste(discovery$SNPName, discovery$ProbeName, sep = '_')
dim(discovery)
discovery <- discovery[, c(23, 10, 11), with = F]
dim(discovery)
colnames(discovery) <- paste0('Discovery_', colnames(discovery))
replication$pair <- paste(replication$SNPName, replication$ProbeName, sep = '_')
colnames(replication)[c(10, 11, 23)] <- c('AlleleAssessed', 'OverallZScore', 'pair')
replication <- replication[, c(23, 10, 11), with = F]
colnames(replication) <- paste0('Replication_', colnames(replication))

replication_sig$pair <- paste(replication_sig$SNPName, replication_sig$ProbeName, sep = '_')
colnames(replication_sig)[c(10, 11, 23)] <-c('AlleleAssessed', 'OverallZScore', 'pair')

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
         replication_perc = round(nr_of_replicated_same_dir/nr_of_tested_in_replication, digits = 4) * 100,
         replicated_same_dir_perc = round(nr_of_replicated_same_dir/nr_of_sig_in_replication, digits = 4) * 100)

# define limits for x and y axes
x_axis <- max(abs(min(all_merged$Discovery_OverallZScore)), abs(max(all_merged$Discovery_OverallZScore, na.rm = T)), na.rm = T)
x_axis <- x_axis + 0.1 * x_axis
y_axis <- max(abs(min(all_merged$Replication_OverallZScore)), abs(max(all_merged$Replication_OverallZScore, na.rm = T)), na.rm = T)
y_axis <- y_axis + 0.1 * y_axis
         
ann_text <- data.frame(Discovery_OverallZScore = -x_axis, 
                       Replication_OverallZScore = y_axis,
                       lab = paste0("Nr. of discovery associations: \n", statistics_overview$nr_of_disc_assoc, 
                                   "\nTested in replication cohort: \n", statistics_overview$nr_of_tested_in_replication, ' (', statistics_overview$tested_eqtl_perc, '%)',
                                   "\nSig. replicated in replication cohort: \n", statistics_overview$nr_of_sig_in_replication, ' (', round(statistics_overview$nr_of_sig_in_replication/statistics_overview$nr_of_tested_in_replication * 100, 2), '%)',
                                   "\nSig. replicated with same direction: \n", statistics_overview$nr_of_replicated_same_dir, ' (', round(statistics_overview$nr_of_replicated_same_dir/statistics_overview$nr_of_sig_in_replication * 100, 2), '%)'),
                       `Replication\nFDR<0.05` = factor('yes', levels = c('yes', 'no')))

colnames(ann_text)[ncol(ann_text)] <- "Sig. replication"

###
colnames(all_merged)[8] <- 'Sig. replication'
colnames(all_merged)[6] <- 'FDR<0.05'
all_merged <- all_merged[!is.na(all_merged$`FDR<0.05`), ]

all_merged$`FDR<0.05` <- factor(all_merged$`FDR<0.05`, 
                                            levels = c('yes', 'no'))
all_merged$`Effect in replication\ncohort:` <- 'non sig.'
all_merged[all_merged$`FDR<0.05` == 'yes' & all_merged$same_dir == 'no', ]$`Effect in replication\ncohort:` <- 'FDR<0.05 and opp. dir.'
all_merged[all_merged$`FDR<0.05` == 'yes' & all_merged$same_dir == 'yes', ]$`Effect in replication\ncohort:` <- 'FDR<0.05 and same dir.'
all_merged$`Effect in replication\ncohort:` <- factor(all_merged$`Effect in replication\ncohort:`, levels = c('non sig.', 'FDR<0.05 and opp. dir.', 'FDR<0.05 and same dir.'))

p <- ggplot(all_merged, aes(x = Discovery_OverallZScore, y = Replication_OverallZScore, colour = `Effect in replication\ncohort:`)) + 
  geom_point(alpha = 0.35) + 
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.5) + 
  scale_color_manual(values = c("FDR<0.05 and same dir." = "red", "FDR<0.05 and opp. dir." = "orange", "non sig." = "darkgrey")) + 
  xlab('Discovery analysis Z-score') + 
  ylab('Replication analysis Z-score') +
  geom_text(data = ann_text, aes(label = lab), size = 4, hjust=0, vjust = 1, col = 'black') +
  ggtitle(paste0('Discovery cohort: ', args[5], ' (N=', Ndisc, ')\nReplication cohort: ', args[6], ' (N=', Nrep, ')')) +
  scale_y_continuous(limits = c(-y_axis, y_axis)) + 
  scale_x_continuous(limits = c(-x_axis, x_axis))

ggsave(paste0(args[4], '_', Sys.Date(), '.png'), height = 7 * 1.2, width = 8 * 1.2, dpi = 600, units = 'in')

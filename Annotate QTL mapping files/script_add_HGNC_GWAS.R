######################################################################################################################################################
# This scripts does following:                                                                                                                       #
# - Adds HGNC (ENSEMBL v75) identifiers for the ouput file of eQTLMapping pipeline, given that ENSEMBL gene names are given in the ProbeMapping file #
# - Adds annotation from separate GWAS SNP table.                                                                                                    #
######################################################################################################################################################

library(data.table)
library(stringr)
library(biomaRt)
library(dplyr)
library(tidyr)

# annotation:
annot <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/ProbeAnnotation2_CorrectlyMapped_UniquelyMapping_AdditionalInfo.txt')

# QTL file:
qtl <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/trans_PRS_meta_analysis_20170115/trans_meta_analysis_20170115/eQTLsFDR0.05-ProbeLevel_PC_corrected.txt')

# Get probe to one SNP associations:

probe_gene <- unique(annot[, c(24, 1), with = F])

probe_gene <- probe_gene %>% 
  mutate(HGNCName_single = strsplit(as.character(GeneName), ";")) %>% 
  unnest(HGNCName_single)

probe_gene <- unique(probe_gene[, c(2, 3)])

# get HGNC gene symbols (ENSEMBL v75)

ensembl_names <- probe_gene$HGNCName_single
ensembl_names <- ensembl_names[!ensembl_names == '-']

ensembl75 <- useMart('ENSEMBL_MART_ENSEMBL', host = "feb2014.archive.ensembl.org", dataset = 'hsapiens_gene_ensembl')
genes <- getBM(ensembl75, attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = ensembl_names)

probe_gene <- merge(probe_gene, genes, by.x = 'HGNCName_single', by.y = 'ensembl_gene_id', all.x = T)

probe_gene[is.na(probe_gene$HGNCName_single) | probe_gene$HGNCName_single %in% '', ]$HGNCName_single  <- '-'
probe_gene[is.na(probe_gene$hgnc_symbol) | probe_gene$hgnc_symbol %in% '', ]$hgnc_symbol  <- '-'

probe_gene <- probe_gene %>% 
  group_by(Probe) %>% summarise(ENSEMBL_75_name = paste(unique(HGNCName_single), collapse="; "), ENSEMBL_75_symbol = paste(unique(as.character(hgnc_symbol)), collapse="; "))

probe_gene <- as.data.frame(probe_gene)

probe_gene[probe_gene$ENSEMBL_75_name == 'NA' | probe_gene$ENSEMBL_75_name == '', ]$ENSEMBL_75_name  <- '-'
probe_gene[probe_gene$ENSEMBL_75_symbol == 'NA' | probe_gene$ENSEMBL_75_symbol == ' ', ]$ENSEMBL_75_symbol  <- '-'


qtl2 <- merge(qtl, probe_gene, by.x = 'ProbeName', by.y = 'Probe', all.x = T)
qtl2 <- qtl2[, c(2:5, 1, 6:ncol(qtl2)), with = F]
qtl2[is.na(qtl2$ENSEMBL_75_symbol), ]$ENSEMBL_75_symbol <- '-'
qtl2 <- qtl2[order(qtl2$PValue), ]

# add GWAS annotation

snp_annot <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/genetic_risk_factors_19042016/merged_files_for_analysis/annotated_genetic_risk_factors_21112016.txt')
snp_annot <- unique(snp_annot[, c(9, 4, 7, 5), with = F])

snp_annot <- snp_annot %>% 
  group_by(ID) %>% summarise(traits = paste(trait, collapse="; "), P = paste(P, collapse="; "), PUBIDs = paste(pubid, collapse="; "))

snp_annot <- as.data.table(snp_annot)
qtl3 <- merge(qtl2, snp_annot, by.x = 'SNPName', by.y = 'ID', all.x = T)
qtl3 <- qtl3[, c(2, 1, 3:ncol(qtl3)), with = F]
colnames(qtl3)[25:27] <- c('GWAS_traits', 'GWAS_Pvalues', 'GWAS_PUBIDs')
qtl3 <- qtl3[order(abs(qtl3$OverallZScore), decreasing = T), ]

write.table(qtl3, '/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/trans_PRS_meta_analysis_20170115/trans_meta_analysis_20170115/Interpretation/trans_eQTLsFDR0.05-ProbeLevel_PC_corrected_ENSEMBL75_HGNC_GWAS_information_added.txt', sep = '\t', quote = F, row.names = F)


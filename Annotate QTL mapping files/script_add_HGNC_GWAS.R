######################################################################################################################################################
# This scripts does following:                                                                                                                       #
# - Adds HGNC (ENSEMBL v75) identifiers for the ouput file of eQTLMapping pipeline, given that ENSEMBL gene names are given in the ProbeMapping file #
# - Adds annotation from separate GWAS SNP table.                                                                                                    #
######################################################################################################################################################

library(data.table)
library(stringr)
library(biomaRt)
library(dplyr)

# annotation:
annot <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/ProbeAnnotation2_CorrectlyMapped_UniquelyMapping_AdditionalInfo.txt')

# QTL file
qtl <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/trans_PRS_meta_analysis_20170115/trans_meta_analysis_20170115/eQTLsFDR0.05-ProbeLevel_PC_corrected.txt')

# get HGNC gene symbols (ENSEMBL v75)

ensembl_names <- annot$GeneName
ensembl_names <- ensembl_names[!ensembl_names == '-']

ensembl75 <- useMart('ENSEMBL_MART_ENSEMBL', host = "feb2014.archive.ensembl.org", dataset = 'hsapiens_gene_ensembl')
genes <- getBM(ensembl75, attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = ensembl_names)

annot <- merge(annot, genes, by.x = "GeneName", by.y = "ensembl_gene_id", all.x = T)
annot_a <- annot[, c(2, 1, ncol(annot)), with = F]
qtl2 <- merge(qtl, annot_a, by.x = 'ProbeName', by.y = 'Probe', all.x = T)
qtl2 <- qtl2[, c(2:5, 1, 6:ncol(qtl2)), with = F]
qtl2[is.na(qtl2$hgnc_symbol), ]$hgnc_symbol <- '-'
qtl2 <- qtl2[order(qtl2$PValue), ]
colnames(qtl2)[ncol(qtl2)] <- 'HGNC_symbol_ENSEMBL75'

# add GWAS annotation

snp_annot <- fread('/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/genetic_risk_factors_19042016/merged_files_for_analysis/annotated_genetic_risk_factors_21112016.txt')
snp_annot <- unique(snp_annot[, c(9, 4, 7, 5), with = F])

snp_annot <- snp_annot %>% 
  group_by(ID) %>% summarise(traits = paste(trait, collapse="; "), P = paste(P, collapse="; "), PUBIDs = paste(pubid, collapse="; "))

snp_annot <- as.data.table(snp_annot)
prs_unc3 <- merge(prs_unc2, snp_annot, by.x = 'SNPName', by.y = 'ID')
prs_unc3 <- prs_unc3[, c(2, 1, 3:ncol(prs_unc3)), with = F]
colnames(prs_unc3)[25:27] <- c('GWAS_traits', 'GWAS_Pvalues', 'GWAS_PUBIDs')
prs_unc3 <- prs_unc3[order(prs_unc3$PValue), ]

write.table(prs_unc3, '/Users/urmovosa/Documents/move_to_mac/trans_eQTL_meta_analysis/trans_PRS_meta_analysis_20170115/trans_meta_analysis_20170115/Interpretation/trans_eQTLsFDR0.05-ProbeLevel_PC_corrected_ENSEMBL75_HGNC_GWAS_information_added.txt', sep = '\t', quote = F, row.names = F)


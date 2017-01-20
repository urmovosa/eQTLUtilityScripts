#################################################################################################
# This script takes the standard output of the eQTLMappingPipeline and constructs the dot-plot  #
# SNP and probe positions have to be in hg19 coordinates                                        #
# This script also needs the file chr_sizes, which was extracted from UCSC (GRCh37)             #
# ###############################################################################################


library(data.table)
library(ggplot2)

# Read in chromosome sizes
chromosome <- fread('chr_sizes.txt')
xy <- chromosome[chromosome$chr %in% c('X', 'Y'), ]
chromosome <- chromosome[chromosome$chr %in% c(1:22), ]
chromosome$chr <- as.numeric(chromosome$chr)
chromosome$abs_coord <- 0
chromosome <- chromosome[order(chromosome$chr), ]

chromosome[1, ]$abs_coord <- chromosome[1, ]$size 

for (i in 2:nrow(chromosome)){
  
  chromosome[i, ]$abs_coord <- chromosome[i - 1, ]$abs_coord + chromosome[i, ]$size
  
}



# chrom tick marks

chromosome$tick_mark <- chromosome$abs_coord - chromosome$size/2
chromosome$tick_mark[1] <- chromosome$abs_coord[1]/2
max_chr <- chromosome[chromosome$chr == 'Y', ]$abs_coord + chromosome[chromosome$chr == 'Y', ]$size

# visualize

p <- ggplot(chromosome, aes(x = abs_coord, y = abs_coord)) + theme_bw() + 
  geom_vline(xintercept = c(chromosome$abs_coord, 0, max_chr), colour = 'lightgrey') + geom_hline(yintercept = c(chromosome$abs_coord, 0, max_chr), colour = 'lightgrey') + 
  scale_x_continuous(breaks = chromosome$tick_mark, labels = paste("chr", chromosome$chr, sep = ''), limits = c(0, max(chromosome$abs_coord)), expand = c(0.01, 0.01)) + 
  scale_y_continuous(breaks = chromosome$tick_mark, labels = paste("chr", chromosome$chr, sep = ''), limits = c(0, max(chromosome$abs_coord)), expand = c(0.01, 0.01)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + xlab('SNP position (hg19)') + ylab('Probe/gene center position (hg19)')


# sort and tweak the data (trans or cis + trans eQTLsFDR0.05-Probelevel.txt)
trans <- fread('eQTLsFDR0.05-ProbeLevel.txt')

trans_f <- trans[, c(1, 2, 3, 4, 5, 6, 7, 11, 25:27), with = F]

# Calculate absolute chromosome coordinates for SNPs and probes

#
trans_f$abs_SNPChrPos <- 0
trans_f[trans_f$SNPChr == 1, ]$abs_SNPChrPos <- trans_f[trans_f$SNPChr == 1, ]$SNPChrPos


for (i in 2:22){
  
  trans_f[trans_f$SNPChr == i, ]$abs_SNPChrPos <- trans_f[trans_f$SNPChr == i, ]$SNPChrPos + chromosome[chromosome$chr == i - 1, ]$abs_coord
  
}

#
trans_f$abs_ProbeChrPos <- 0
trans_f[trans_f$ProbeChr == 1, ]$abs_ProbeChrPos <- trans_f[trans_f$ProbeChr == 1, ]$ProbeCenterChrPos

for (i in 2:22){
  
  trans_f[trans_f$ProbeChr == i, ]$abs_ProbeChrPos <- trans_f[trans_f$ProbeChr == i, ]$ProbeCenterChrPos + chromosome[chromosome$chr == i - 1, ]$abs_coord
  
}

# visualize locations on dot-plot

p + geom_point(data = as.data.frame(trans_f), aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = -log10(PValue)), alpha = 0.15, colour = 'black') + 
  scale_size_continuous(breaks = c(10, 50, 100, 200, 300), range = c(1, 5), guide = guide_legend(title = expression(paste(-log[10]("P-value")))))

ggsave('dot_plot_BIOS.png', height = 9, width = 9 * 1.2, dpi = 400)

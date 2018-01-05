library(data.table)
library(ggplot2)

## Setup the graph

# Read in chromosome sizes
chromosome <- data.frame(chr = c("1", "2", "3", "4", "5", "6", "7", "X", "8", "9", "10",
                                 "11", "12", "13", "14", "15", "16", "17", "18", "20", "Y", "19", "22", "21"),
size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 155270560, 146364022,
         141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210,
         78077248,  63025520,  59373566,  59128983, 51304566, 48129895))
chromosome$chr <- as.character(chromosome$chr)

xy <- chromosome[chromosome$chr %in% c('X', 'Y'), ]
chromosome <- chromosome[chromosome$chr %in% c(1:22), ]
chromosome$chr <- as.numeric(chromosome$chr)
chromosome$abs_coord <- 0
chromosome <- chromosome[order(chromosome$chr), ]

chromosome[1, ]$abs_coord <- chromosome[1, ]$size 

for (i in 2:nrow(chromosome)){
  
  chromosome[i, ]$abs_coord <- chromosome[i - 1, ]$abs_coord + chromosome[i, ]$size
  
}

xy$abs_coord <- NA
chromosome <- rbind(chromosome, xy)

chromosome$abs_coord[23] <- chromosome$abs_coord[22] + chromosome$size[23]
chromosome$abs_coord[24] <- chromosome$abs_coord[23] + chromosome$size[24]

# chrom tick marks

chromosome$tick_mark <- chromosome$abs_coord - chromosome$size/2
chromosome$tick_mark[1] <- chromosome$abs_coord[1]/2
max_chr <- chromosome[chromosome$chr == 'Y', ]$abs_coord + chromosome[chromosome$chr == 'Y', ]$size

# visualize

p <- ggplot(chromosome, aes(x = abs_coord, y = abs_coord)) + theme_bw() + 
  geom_vline(xintercept = c(chromosome$abs_coord, 0, max_chr), colour = 'lightgrey') + geom_hline(yintercept = c(chromosome$abs_coord, 0, max_chr), colour = 'lightgrey') + 
  scale_x_continuous(breaks = chromosome$tick_mark, labels = paste("chr", chromosome$chr, sep = ''), limits = c(0, max(chromosome$abs_coord[1:22])), expand = c(0.01, 0.01)) + 
  scale_y_continuous(breaks = chromosome$tick_mark, labels = paste("chr", chromosome$chr, sep = ''), limits = c(0, max(chromosome$abs_coord)), expand = c(0.01, 0.01)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + xlab('SNP position (hg19)') + ylab('Probe/gene center position (hg19)')


# sort and tweak the data (trans or cis + trans eQTLsFDR0.05-Probelevel.txt)
trans <- fread('eQTLsFDR0.05-ProbeLevel.txt')

# Calculate absolute chromosome coordinates for SNPs and probes

trans$abs_SNPChrPos <- 0
trans[trans$SNPChr == 1, ]$abs_SNPChrPos <- trans[trans$SNPChr == 1, ]$SNPChrPos


for (i in 2:22){
  
  trans[trans$SNPChr == i, ]$abs_SNPChrPos <- trans[trans$SNPChr == i, ]$SNPChrPos + chromosome[chromosome$chr == i - 1, ]$abs_coord
  
}

#
trans$abs_ProbeChrPos <- 0
trans[trans$ProbeChr == 1, ]$abs_ProbeChrPos <- trans[trans$ProbeChr == 1, ]$ProbeCenterChrPos
trans$abs_ProbeChrPos <- as.numeric(trans$abs_ProbeChrPos)

for (i in c(2:24)){
  
  trans[trans$ProbeChr == chromosome$chr[i], ]$abs_ProbeChrPos <- as.numeric(trans[trans$ProbeChr == chromosome$chr[i], ]$ProbeCenterChrPos) + as.numeric(chromosome[chromosome$chr == chromosome$chr[i - 1], ]$abs_coord)
  
}

# visualize locations on dot-plot

# Size is -log10(P-value)
p + geom_point(data = as.data.frame(trans), aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = -log10(PValue)), alpha = 0.15, colour = 'black') + 
  scale_size_continuous(breaks = c(10, 50, 100, 200, 300), range = c(1, 5), guide = guide_legend(title = expression(paste(-log[10]("P-value"))))) + 
  theme(axis.title=element_text(size = 18, face = "bold")) + ylab('Gene position (hg19)')

ggsave('dot_plot.png', height = 9, width = 9 * 1.1, dpi = 400)

# Size is Z-score instead of P-value
p + geom_point(data = as.data.frame(trans), aes(x = abs_SNPChrPos, y = abs_ProbeChrPos, size = abs(OverallZScore)), alpha = 0.15, colour = 'black') + 
  scale_size_continuous( range = c(1, 5), guide = guide_legend(title = expression(paste("meta-analysis Z-score"))))

ggsave('dot_plot_Z.png', height = 9, width = 9, dpi = 400)

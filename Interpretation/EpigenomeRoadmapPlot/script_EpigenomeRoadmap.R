###################################################################################
# This script makes a chrohmm plot for given region, using Epigenome Roadmap data #
###################################################################################

library(AnnotationHub)
library(stringr)
library(tidyr)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

# Specify the positions
chr <- args[1]
pos_start <- as.numeric(args[2])
pos_end <- as.numeric(args[3])
graphics_device <- args[4]

print(paste0('Chromosome: ', chr))
print(paste0('Start position: ', pos_start))
print(paste0('End position: ', pos_end))
print(paste0('Will write file into format: ', graphics_device))

# Color scheme for annotations:

ah <- AnnotationHub()

epiFiles <- query(ah, "EpigenomeRoadMap")
epiFiles <- epiFiles[epiFiles$description == "15 state chromatin segmentations from EpigenomeRoadMap Project", ]

sourcetype <- as.data.frame(epiFiles$tags)
colnames(sourcetype) <- "tag"
sourcetype$tag <- as.character(sourcetype$tag)
sourcetype <- str_split_fixed(sourcetype$tag, ", ", 8)

full_marks <- as.data.frame(matrix(NA, nrow = 2, ncol = 10))
colnames(full_marks) <- c(
  "seqnames", "start", "end", "width", "strand", "abbr",
  "name", "color_name", "color_code", "tissue"
)

full_marks <- full_marks[-c(1, 2), ]

for (i in c(1:127)) {
  marks <- epiFiles[[c(epiFiles$ah_id[i])]]
  marks$sample <- sourcetype[, 5][i]
  marks$tissue <- sourcetype[, 6][i]
  # marks <- marks[seqnames(marks) == 'chr11' & ranges(marks)]
  marks <- as.data.frame(marks)
  # marks <- marks[marks$start > pos_start - 1000 & marks$end < pos_end + 1000, ]

  full_marks <- rbind(full_marks, marks)
  print(i)
}

full_marks <- full_marks[full_marks$seqnames == chr & full_marks$start > pos_start - 1000 & full_marks$end < pos_end + 1000, ]

full_marks$end <- as.numeric(full_marks$end)
full_marks <- full_marks[order(full_marks$tissue, full_marks$sample, full_marks$seqnames, full_marks$start), ]
full_marks$start <- as.numeric(full_marks$start)

# Add colors for tissues:

tissue_colors <- unique(full_marks[, c(10, 11)])
tissue_colors$tissue <- as.factor(tissue_colors$tissue)

tissue_colors$x_pos <- 5
tissue_colors$y_pos <- c(1:nrow(tissue_colors))
tissue_colors$color <- "unknown"

tissue_colors[tissue_colors$tissue == "ES-deriv", ]$color <- "goldenrod"
tissue_colors[tissue_colors$tissue == "ESC", ]$color <- "darkgrey"
tissue_colors[tissue_colors$tissue == "iPSC", ]$color <- "lightgreen"
tissue_colors[tissue_colors$tissue == "Mesench", ]$color <- "salmon"
tissue_colors[tissue_colors$tissue == "IMR90", ]$color <- "#D1CE77"
tissue_colors[tissue_colors$tissue == "Blood & T-cell", ]$color <- "darkgreen"
tissue_colors[tissue_colors$tissue == "Epithelial", ]$color <- "navy"
tissue_colors[tissue_colors$tissue == "HSC & B-cell", ]$color <- "orange"
tissue_colors[tissue_colors$tissue == "Adipose", ]$color <- "darkgoldenrod1"
tissue_colors[tissue_colors$tissue == "Brain", ]$color <- "red"
tissue_colors[tissue_colors$tissue == "Digestive", ]$color <- "lightblue"
tissue_colors[tissue_colors$tissue == "ENCODE2012", ]$color <- "brown"
tissue_colors[tissue_colors$tissue == "Heart", ]$color <- "blue"
tissue_colors[tissue_colors$tissue == "HSC & B-cell", ]$color <- "firebrick"
tissue_colors[tissue_colors$tissue == "Muscle", ]$color <- "skyblue"
tissue_colors[tissue_colors$tissue == "Myosat", ]$color <- "pink"
tissue_colors[tissue_colors$tissue == "Neurosph", ]$color <- "sienna1"
tissue_colors[tissue_colors$tissue == "Other", ]$color <- "aquamarine"
tissue_colors[tissue_colors$tissue == "Sm. Muscle", ]$color <- "grey"
tissue_colors[tissue_colors$tissue == "Thymus", ]$color <- "black"

############################
# Here starts the plotting #
############################

if (graphics_device == "png") {
  png(paste(paste(chr, pos_start, pos_end, sep = "_"), "_", Sys.Date(), "_.png", sep = ""), width = 11, height = 4, units = "in", res = 600)
}
if (graphics_device == "pdf") {
  pdf(paste(paste(chr, pos_start, pos_end, sep = "_"), "_", Sys.Date(), "_.pdf", sep = ""), width = 11, height = 4)
}

layout(matrix(c(1, 2, 3), 1, 3),
  widths = c(15, 40, 15), heights = c(1, 1, 1)
)

par(mar = c(0, 0, 0, 0), oma = c(4, 2, 2, 3))

plot(NULL, xlim = c(1, 6), ylim = c(0, length(unique(full_marks$sample))), axes = F)

rect(tissue_colors$x_pos,
  tissue_colors$y_pos - 1,
  tissue_colors$x_pos + 0.5,
  tissue_colors$y_pos,
  border = NA,
  col = tissue_colors$color
)

# This is for lines separating the tissues

lines_help <- data.frame(
  x_pos_start = rep(min(full_marks$start), times = length(table(tissue_colors$tissue))),
  x_pos_end = rep(max(full_marks$end), times = length(table(tissue_colors$tissue))),
  y_pos = NA
)

start_pos <- 0
for (i in 1:length(table(tissue_colors$tissue))) {
  lines_help$y_pos[i] <- start_pos + as.numeric(table(tissue_colors$tissue)[i])
  start_pos <- lines_help$y_pos[i]
}

lines_help <- rbind(lines_help, data.frame(x_pos_start = min(full_marks$start), x_pos_end = max(full_marks$end), y_pos = 0))


# add text annotation for the tissues

y_position <- seq(
  from = 0,
  to = length(unique(full_marks$sample)),
  by = length(unique(full_marks$sample)) / (length(unique(full_marks$tissue)) - 1)
)

y_position <- y_position
x_position <- c(min(lines_help$x_pos_start) - ((min(lines_help$x_pos_end) - max(lines_help$x_pos_start))) * 0.7)

text(labels = unique(tissue_colors$tissue), x = 4, y = y_position, pos = 2, col = unique(tissue_colors$color), cex = 1.5)

# add lines for annotations
abi <- 0
y_positions <- NA

for (i in 1:(nrow(lines_help) - 1)) {
  y_positions[i] <- ((lines_help$y_pos[i] - abi) - (lines_help$y_pos[i] - abi) / 2) + abi
  abi <- lines_help$y_pos[i]
}

segments(x0 = 4, y0 = y_position, x1 = 5, y1 = y_positions, col = unique(tissue_colors$color), cex = 2)

# plot itself

plot(NULL, xlim = c(pos_start / 1000000, pos_end / 1000000), ylim = c(0, length(unique(full_marks$sample))), axes = F, xlab = "", ylab = "")
mtext(paste(chr, " position (Mb)(hg19)", sep = ""), side = 1, line = 2)

j <- 0
for (i in rev(unique(full_marks$sample))) {
  rect(full_marks[full_marks$sample == i, ]$start / 1000000,
    rep(j, nrow(full_marks[full_marks$sample == i, ])),
    full_marks[full_marks$sample == i, ]$end / 1000000,
    rep(j + 1, nrow(full_marks[full_marks$sample == i, ])),
    col = full_marks[full_marks$sample == i, ]$color_code,
    border = NA
  )

  j <- j + 1
}

axis(1)
box()

# add lines sparating the tissues

for (i in 1:nrow(lines_help)) {
  lines(c(lines_help$x_pos_start[i] / 1000000 - 1000, lines_help$x_pos_end[i] / 1000000 + 1000), c(lines_help$y_pos[i], lines_help$y_pos[i]), lwd = 1)
}

# Add annotation for chromatin regions

col_abi <- unique(full_marks[, c(7, 9)])
col_abi <- col_abi[match(rev(c(
  "Active TSS", "Flanking Active TSS", "Transcr. at gene 5' and 3'",
  "Strong transcription", "Weak transcription", "Genic enhancers", "Enhancers",
  "ZNF genes & repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed PolyComb", "Weak Repressed PolyComb",
  "Quiescent/Low"
)), col_abi$name), ]


plot(NULL, xlim = c(1, 6), ylim = c(0, length(unique(full_marks$sample))), axes = F)

rect(1,
  seq(
    from = 0,
    to = length(unique(full_marks$sample)),
    by = length(unique(full_marks$sample)) / 15
  )[-16],
  1.5,
  seq(
    from = 0,
    to = length(unique(full_marks$sample)),
    by = length(unique(full_marks$sample)) / 15
  )[-1],
  border = "black",
  col = col_abi$color_code
)

text(
  x = 1.5, y = seq(from = 0, to = length(unique(full_marks$sample)), by = length(unique(full_marks$sample)) / 15)[-1] -
    (seq(from = 0, to = length(unique(full_marks$sample)), by = length(unique(full_marks$sample)) / 15)[-1] - seq(from = 0, to = length(unique(full_marks$sample)), by = length(unique(full_marks$sample)) / 15)[-16]) / 2,
  label = col_abi$name,
  pos = 4,
  cex = 1.25
)


dev.off()
# Qiwen Hu 2018
# Generate peak region and count matrix for each cell
# Input: aligned bam file with barcode tag
# Output: peak regions and count matrix
# Usage: Rscript sc.count.R bamfile

source("../src/core.R")
library(dplyr)
library(GenomicRanges)

# Function to convert from ranges to Granges
# From Jean Fan
range2GRanges <- function(df){
# Input: df Dataframe with columns as sequence name, start, and end
# Output: GRanges version 
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3])
  )
  return(gr)
}

# Command line arguments
options <- commandArgs(trailingOnly = T)
output.prefix <- options[[1]]

# read bam files
filename <- paste(output.prefix, ".bam", sep = "")
rl <- read.bam.tags(filename)

# get valid chromosomes
chrl <- paste0("chr", c(1:19, "X", "Y"))
pos <- rl$rl$tags[which(names(rl$rl$tags) %in% chrl)]
bindex <- rl$rl$bindex[which(names(rl$rl$tags) %in% chrl)]

bcode.index <- rl$index
bcode.index <- bcode.index[order(bcode.index$ID), ]

# get peaks for each chromosome
bandwidth <- 300
step = 100
thr = 5
span = 10
peaks <- lapply(pos, t.get.density.peaks, step=step, thr=thr, span=span)

# get peak regions 
peaks.region <- do.call(rbind, lapply(chrl, function(chr){
	df <- peaks[[chr]]
	data.frame(chr, df[,1] - bandwidth/2, df[,1] + bandwidth/2)
}))

peaks.region <- range2GRanges(peaks.region)
names(peaks.region) <- paste(peaks.region)

# write out bed file
gr <- peaks.region
bed.df <- data.frame(seqnames = seqnames(gr),
                     starts = start(gr) - 1,
                     ends = end(gr),
                     names = c(rep(".", length(gr))),
                     scores = c(rep(".", length(gr))),
                     strands = strand(gr))
write.table(bed.df, file = paste0(output.prefix, "_peaks.1", ".bed"), 
            quote = F, sep = "\t", row.names = F, col.names = F)

# generate count matrix

nbcode <- max(unlist(lapply(bindex, max)))
t <- mapply(ccwindow_n_tags_barcodes, pos, bindex, nbcode, peaks, bandwidth)
count.matrix <- do.call(rbind, t)
colnames(count.matrix) <- bcode.index$barcode
#rownames(count.matrix) <- names(peaks.region)
save(count.matrix, file = paste0(output.prefix, ".counts.RData"))

#save column
write.table(bcode.index$barcode, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t",
              paste0(output.prefix, ".cmatrix.col.txt"))

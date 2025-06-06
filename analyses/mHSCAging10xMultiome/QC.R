# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(Seurat)
library(SummarizedExperiment)
library(ggplot2)
library(ArchR)
source("../../code/utils.R")

########################
# Make depth-FRIP plot #
########################

dataset <- "mHSCAging10xMultiome"

# Load ATAC peaks
peakBed <- read.table(paste0("../../data/", dataset, "/peaks.bed"), sep = "\t", header = F)
peaks <- GRanges(seqnames = peakBed$V1,
                 ranges = IRanges(start = peakBed$V2, end = peakBed$V3))

# Get peak-by-cell count matrix
scATAC <- getCountsFromFrags(paste0("../../data/", dataset, "/all.frags.tsv.gz"), peaks)

# Plot Depth-FRIP plot
plotData <- data.frame(
  FRIP = scATAC$FRIP,
  depth = scATAC$depth,
  sample = scATAC$sample,
  density = get_density(log10(scATAC$depth), scATAC$FRIP, n = 100)
)
pdf(paste0("../../data/", dataset, "/plots/QC/FRIPDepth.pdf"),
    height = 5, width = 5)
ggplot(plotData) +
  geom_point(aes(x = depth, y = FRIP, color = density)) +
  scale_x_continuous(trans = "log10")
dev.off()

############################
# Visualize fragment sizes #
############################

# Load fragment coordinates
frags <- data.table::fread(paste0("../../data/", dataset, "/all.frags.filt.tsv.gz"), 
                           sep = "\t", showProgress = TRUE, nThread = 16) 

# Calculate fragment lengths
fragLen <- frags$V3 - frags$V2

# Calculate histogram
histogram <- table(factor(fragLen, levels = 1:max(fragLen)))
histogram <- as.numeric(histogram)

# Visualize
plotData <- data.frame(
  length = 1:length(histogram),
  freq = histogram / sum(histogram)
)
plotData <- plotData[plotData$length < 1000, ]

pdf(paste0("../../data/", dataset, "/plots/QC/fragSizes.pdf"),
    width = 7, height = 5.5)
ggplot(plotData) +
  geom_line(aes(x = length, y = freq), size = 1, color = "#2171B5") +
  xlab("Fragment size") + ylab("Frequency") +
  theme_classic()
dev.off()

######################
# Get TSS enrichment #
######################

TSSEnrichment <- getTSSEnrichment(frags, "mm10")

pdf(paste0("../../data/", dataset, "/plots/QC/TSSPlot.pdf"),
    width = 7, height = 5.5)
ggplot(TSSEnrichment) +
  geom_line(aes(x = position, y = normValue), size = 1) +
  theme_classic()
dev.off()


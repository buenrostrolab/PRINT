# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(Seurat)
library(SummarizedExperiment)
library(ggplot2)
library(ArchR)
library(BuenColors)
source("../../code/utils.R")

########################
# Make depth-FRIP plot #
########################

dataset <- "K562"

# Create folder to store result plots
system(paste0("mkdir ../../data/", dataset, "/plots/QC"))

# Load ATAC peaks
peakBed <- read.table(paste0("../../data/", dataset, "/hg19Peaks.bed"), sep = "\t", header = F)
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
  geom_point(aes(x = depth, y = FRIP, color = density), size = 0.1) +
  scale_x_continuous(trans = "log10") +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  xlab("Depth") + ylab("FRIP") + ggtitle(dataset) +
  theme_classic()
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

TSSEnrichment <- getTSSEnrichment(frags, "hg38")

pdf(paste0("../../data/", dataset, "/plots/QC/TSSPlot.pdf"),
    width = 7, height = 5.5)
ggplot(TSSEnrichment) +
  geom_line(aes(x = position, y = normValue), size = 1) +
  theme_classic()
dev.off()


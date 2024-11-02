# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(BuenRTools)
library(data.table)
source("../../code/utils.R")

###############################################
# Get genomic ranges for running footprinting #
###############################################

# Load chromosome sizes
chrSizes <- read.table("../../data/shared/sacCer3.chrom.sizes", sep = "\t")

# Create GRanges for each chromosome
chrRanges <- GRanges(paste0(chrSizes$V1, ":1-", chrSizes$V2))

# Divide the whole genome in to bins
binSize <- 5000
chrBins <- Reduce(c, slidingWindows(chrRanges, width = binSize, step = binSize))

# Only keep bins with the correct size
regions <- chrBins[width(chrBins) == binSize]

# Convert to bed format
regionBed <- data.frame(
  chr = seqnames(regions),
  start = start(regions),
  end = end(regions)
)

# Save to file
saveRDS(regions, "../../data/yeast/regionRanges.rds")
write.table(regionBed, "../../data/yeast/regions.bed", sep = "\t",
            quote = F, row.names = F, col.names = F)

###################################################################
# Add barcode to fragments file so it matches the standard format #
###################################################################

frags <- data.table::fread("../../data/yeast/merged.frags.gz", sep = "\t")
if(dim(frags)[2] == 4){
  frags$V5 <- "yeast"
  frags <- frags[, c(1,2,3,5,4)]
  data.table::fwrite(frags, "../../data/yeast/merged.frags.gz", sep = "\t", quote = F, 
                     row.names = F, col.names = F, nThread = 16, showProgress = T)
}

##################
# Clean up peaks #
##################

# Path to summit file(s)
peakCallingFiles <- list.files("../../data/yeast/ATAC/peakCalling")
summitFile <- peakCallingFiles[stringr::str_detect(peakCallingFiles, "summits")]
chrSizes <- "../../data/shared/refGenomes/sacCer3.chrom.sizes"

# Call summitsToCleanPeaks function for peak clean-up
summitsToCleanPeaks(summitFiles = paste0("../../data/yeast/ATAC/peakCalling/", summitFile),
                    peakWidth = 800, # Window to use for peak summit padding
                    blackList = "../../data/shared/BlacklistFiles/placeHolder.bed",
                    chromSizes = chrSizes,
                    topNPeaks = NULL, # Filter top N peaks after?
                    useQuantileRanks = FALSE, # Use normalized ranks instead of raw MACS score?
                    FDRcutoff = 0.01, # MACS peak FDR cut-off used for all peaks
                    outDir = "../../data/yeast/ATAC/peakCalling/", # Make sure this directory exists first
                    expName= "filtered", # Name to give as prefix to resulting filtered peak bed file
                    resizeTo = 1000 # Re-size peaks after filtering?
)

system("mv ../../data/yeast/ATAC/peakCalling/filtered.fixedwidthpeaks_800bp_1000reSized.bed ../../data/yeast/peaks.bed")

# Get counts in peaks
peakBed <- read.table("../../data/yeast/peaks.bed", sep = "\t")
peaks <- GRanges(seqnames = peakBed$V1,
                 ranges = IRanges(start = peakBed$V2, end = peakBed$V3))
counts <- BuenRTools::getCountsFromFrags(
  fragFile = "../../data/yeast/merged.frags.gz",
  peaks = peaks
)
peaks$score <- assay(counts)[, 1]

# Merger overlapping peaks
peaks <- mergeRegions(peaks)

# Save to file
saveRDS(peaks, "../../data/yeast/ATACPeaks.rds")

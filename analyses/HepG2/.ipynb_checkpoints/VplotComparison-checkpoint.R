# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(dplyr)
library(GenomicRanges)
library(ggplot2)
source("../../code/utils.R")
source("../../code/getBias.R")

###########################################################
# Compare V-plot and multi-scale footprints in cell lines #
###########################################################

# Load genomic regions of CREs
regions <- readRDS("../../data/HepG2/regionRanges.rds")

# Load fragments and convert to genomic ranges
fragFile <- "../../data/HepG2/all.frags.filt.tsv.gz"
frags <- data.table::fread(fragFile, sep = "\t", showProgress = TRUE, nrows=Inf) %>% 
  data.frame() %>% GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
                                                           start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, 
                                                           starts.in.df.are.0based = T)

# Load Unibind TF binding sites
UnibindChIPRanges <- readRDS("../../data/shared/unibind/HepG2ChIPRanges.rds")
UnibindChIPRangesAll <- Reduce(c, sapply(
  names(UnibindChIPRanges), 
  function(TF){
    ranges <- UnibindChIPRanges[[TF]]
    if(length(ranges) > 0){
      ranges$TF <- TF
    }
    ranges
  })
)

# For a TF-of-interest, pick a CRE region bound by this TF
TF <- "BHLHE40"
Reduce(intersect, list(
  findOverlaps(regions, motifPositions[[TF]])@from,
  findOverlaps(regions, UnibindChIPRanges[[TF]])@from
  )
)
regionInd <- 140856

# Visualize V-plot in a CRE region
regionFrags <- subsetByOverlaps(frags, regions[regionInd])
midPoints <- start(resize(regionFrags, 1, fix = "center"))
pltRegion <- resize(regions[regionInd], 800, fix = "center")
fragLengths <- width(regionFrags)
pltData <- data.frame(
  fragL = fragLengths,
  midPoint = midPoints
)

# Find the relative position of this TF
TFSite <- subsetByOverlaps(UnibindChIPRanges[[TF]], regions[regionInd])
start(TFSite) - start(pltRegion)

# Visualize results
ggplot(pltData) + 
  geom_point(aes(x = midPoint, y = fragL), size = 0.5) +
  #ylim(0, 300) +
  ylab("Fragment size") + xlab("Mid-point location") +
  scale_x_continuous(limits = c(start(pltRegion), end(pltRegion)), 
                     breaks = seq(start(pltRegion) + 100, end(pltRegion), b=200), 
                     expand = c(0,0)) +
  theme_classic()
  
# List all unibind TF sites in thie region
regionOfInterest <- resize(shift(resize(pltRegion, 1, fix = "start"), 200), 100, fix = "center")
ovTFs <- subsetByOverlaps(UnibindChIPRangesAll, pltRegion)
ovTFs
start(ovTFs) - start(pltRegion)

##########################################################
# Compare V-plot and multi-scale footprints on naked DNA #
##########################################################

# Initialize a footprintingProject object
projectName <- "BAC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
dataDir(project) <- "../../data/BAC"
mainDir(project) <- "../../"

# Load region ranges
regions <- readRDS(paste0(projectDataDir, "tileRanges.rds"))

# Set the regionRanges slot
regionRanges(project) <- regions

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getPrecomputedBias(project, nCores = 24)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load fragments and convert to genomic ranges
fragFile <- "../../data/BAC/rawData/all.fragments.tsv.gz"
frags <- data.table::fread(fragFile, sep = "\t", showProgress = TRUE, nrows=Inf) %>% 
  data.frame() %>% GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
                                                           start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, 
                                                           starts.in.df.are.0based = T)

# Visualize V-plot in a CRE region
regionInd <- 3500
regionFrags <- subsetByOverlaps(frags, regions[regionInd])
regionFrags <- regionFrags[sample(1:length(regionFrags), 5000)]
midPoints <- start(resize(regionFrags, 1, fix = "center"))
pltRegion <- resize(regions[regionInd], 800, fix = "center")
fragLengths <- width(regionFrags)
pltData <- data.frame(
  fragL = fragLengths,
  midPoint = midPoints
)

# Visualize results
ggplot(pltData) + 
  geom_point(aes(x = midPoint, y = fragL), size = 0.5, color = "#FF6B39") +
  ylab("Fragment size") + xlab("Mid-point location") +
  scale_x_continuous(limits = c(start(pltRegion), end(pltRegion)), 
                     breaks = seq(start(pltRegion) + 100, end(pltRegion), b=200), 
                     expand = c(0,0)) +
  theme_classic()

# Plot Tn5 bias
positions <- start(pltRegion):end(pltRegion)
barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                      y1 = 0, y2 = regionBias(project)[regionInd, 101:900])
ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Tn5 bias") +
  theme_classic() 

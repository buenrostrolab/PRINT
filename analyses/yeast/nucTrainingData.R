# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getFootprints.R")
source("../../code/getAggregateFootprint.R")
library(hdf5r)
library(rtracklayer)

###################################
# Load multi-scale footprint data #
###################################

# Load region ranges
regions <- readRDS("../../data/yeast/regionRanges.rds")

projectName <- "yeast"

# Load footprints
footprintRadii <- c(10, 20, 30, 50, 80)
multiScaleFootprints <- list()
for(footprintRadius in footprintRadii){
  print(paste0("Loading data for footprint radius = ", footprintRadius))  
  chunkDir <- paste0("../../data/", projectName, "/chunkedFootprintResults/", footprintRadius, "/")
  chunkFiles <- gtools::mixedsort(list.files(chunkDir))
  scaleFootprints <- pbmcapply::pbmclapply(
    chunkFiles,
    function(chunkFile){
      chunkData <- readRDS(paste0(chunkDir, chunkFile))
      as.data.frame(t(sapply(
        chunkData,
        function(regionData){
          regionData$aggregateScores
        }
      )))
    },
    mc.cores = 16
  )
  scaleFootprints <- data.table::rbindlist(scaleFootprints)
  multiScaleFootprints[[as.character(footprintRadius)]] <- as.matrix(scaleFootprints)
}

# Determine whether or not to keep only data in the CRE region
CREOnly <- T

####################################################
# Load chemically mapped nucleosome occupancy data #
####################################################

# Load bigwig file. Data is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2561059
nucOccupancy <- lapply(
  1:3,
  function(replicate){
    import.bw(paste0("../../data/yeast/nucleosomeH3Q85C/GSM256105", replicate + 6, 
                     "_Occupancy_H3_CC_rep_", replicate, ".bw"))
  }
)

# Load CRE ranges
CREs <- readRDS("../../data/yeast/ATACPeaks.rds")

# Retrieve nucleosome occupancy data for each base pair 
regions$regionInd <- 1:length(regions)
region1bpDataFrame <- data.table::rbindlist(
  pbmcapply::pbmclapply(
    1:length(regions),
    function(regionInd){
      region <- regions[regionInd]
      data.frame(
        chr = as.character(seqnames(region)),
        start = start(region):end(region),
        pos = 1:width(region),
        regionInd = regionInd
      )
    },
    mc.cores = 16
  )
)

# Convert 1bp bins to GRanges object
region1bpPos <- GRanges(
  seqnames = region1bpDataFrame$chr,
  ranges = IRanges(start = region1bpDataFrame$start,
                   end = region1bpDataFrame$start)
)
region1bpPos$pos <- region1bpDataFrame$pos
region1bpPos$regionInd <- region1bpDataFrame$regionInd

# Generate region-by-position matrix of nucleosome occupancy data
nucTracks <- pbmcapply::pbmclapply(
  1:3,
  function(replicate){
    nucOv <- findOverlaps(nucOccupancy[[replicate]], region1bpPos)
    nucTracks <- Matrix::sparseMatrix(
      i = region1bpPos[nucOv@to]$regionInd,
      j = region1bpPos[nucOv@to]$pos,
      x = nucOccupancy[[replicate]][nucOv@from]$score,
      dims = c(length(regions), unique(width(regions)))
    )
  },
  mc.cores = 3
)

# Average over the replicates
nucTracks <- Reduce("+", nucTracks) / 3

# Save resulst
saveRDS(nucTracks, "../../data/yeast/nucTracks.rds")

##########################
# Generate training data #
##########################

regionWid <- unique(width(regions))
contextRadius <- 100
stepSize <- 10

footprintContext <- pbmcapply::pbmclapply(
  1:length(regions),
  function(regionInd){
    
    # Iterate over the position relative to the start of the region
    relativePos <- seq(contextRadius + 1, regionWid - contextRadius, by = stepSize)
    
    regionFootprints <- sapply(
      # Iterate over the position relative to the start of the region
      relativePos,
      function(pos){
        # Retrieve multi-scale footprints 
        contextInds <- (pos - contextRadius) : (pos + contextRadius)
        footprints <- lapply(
          names(multiScaleFootprints ),
          function(kernelSize){
            kernelFootprints <- multiScaleFootprints[[as.character(kernelSize)]][regionInd, contextInds]
          }
        )
        footprints <- Reduce(c, footprints)
      }
    )
    
    regionFootprints <- as.data.frame(t(regionFootprints))
    
    if(CREOnly){
      singlebpPos <- GenomicRanges::slidingWindows(regions[regionInd], width = 1, step = 1)[[1]][relativePos]
      filter <- findOverlaps(singlebpPos, CREs)@from
    }else{
      filter <- 1:length(coverage)
    }
    
    regionFootprints <- regionFootprints[filter, ]
  },
  mc.cores = 16
)
footprintContext <- data.table::rbindlist(footprintContext)
gc()

#################################################
# Generate training targets along with metadata #
#################################################

# Load footprintingProject object
project <- readRDS("../../data/yeast/project.rds")

# Calculate region-by-position matrix of ATAC counts
project <- getATACTracks(project)

# Get nucleosome occupancy label and metadata
nucLabel <- data.table::rbindlist(pbmcapply::pbmclapply(
  1:length(regions),
  function(regionInd){
    
    # Iterate over the position relative to the start of the region
    relativePos <- seq(contextRadius + 1, regionWid - contextRadius, by = stepSize)
    
    # Get nucleosome chemical mapping signal at each position
    regionNucLabel <- nucTracks[regionInd, relativePos]
    
    # Get local ATAC coverage at the same positions
    coverage <- sapply(
      relativePos,
      function(relativePos){
        sum(ATACTracks(project)[regionInd, (relativePos - contextRadius) : (relativePos + contextRadius)])
      }
    )
    
    if(CREOnly){
      singlebpPos <- GenomicRanges::slidingWindows(regions[regionInd], width = 1, step = 1)[[1]][relativePos]
      filter <- findOverlaps(singlebpPos, CREs)@from
    }else{
      filter <- 1:length(coverage)
    }
    
    # Combine results
    data.frame(
      nucOccupancy = regionNucLabel, 
      chr = as.character(seqnames(regions[regionInd])),
      coverage = coverage)[filter, ]
  },
  mc.cores = 16
))

#####################
# Save data to file #
#####################

# Write TFBS training data to a file
if(CREOnly){
  h5_path <- "../../data/yeast/nucDataCRE.h5"
}else{
  h5_path <- "../../data/yeast/nucData.h5"
}

if(file.exists(h5_path)){
  system(paste0("rm ", h5_path))
}
h5file <- H5File$new(h5_path, mode="w")
h5file[["footprints"]] <- as.matrix(footprintContext)
h5file[["nucLabel"]] <- as.data.frame(nucLabel)
h5file[["scales"]] <- footprintRadii
h5file$close_all()

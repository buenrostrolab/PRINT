# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getTFBS.R")
source("../../code/getAggregateFootprint.R")
library(GenomicRanges)
library(patchwork)

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "inVitroFootprints"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regionBed <- read.table(paste0(projectDataDir, "BAC.bed"))
breaks <- seq(regionBed$V2 + 1, regionBed$V3, 1000)
regions <- GRanges(paste0(regionBed$V1, ":", breaks[1:(length(breaks) - 1)], "-", breaks[2:length(breaks)] - 1))
saveRDS(regions, paste0(projectDataDir, "regionRanges.rds"))

# Set the regionRanges slot
regionRanges(project) <- regions

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getPrecomputedBias(project)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load barcodes for each pseudo-bulk
barcodeGroups <- data.frame(
  barcode = c(
    "aligned/FT6_11/FT6_11_S190_001",  "aligned/FT7_11/FT7_11_S200_001",  "aligned/FT7_5/FT7_5_S195_001",
    "aligned/FT5_5/FT5_5_S185_001",  "aligned/FT6_9/FT6_9_S188_001",  "aligned/FT7_2/FT7_2_S192_001",
    "aligned/FT7_8/FT7_8_S198_001"
    ),
  group = c(
    "Control", "100 nM CEBPA", "100 nM FOS", "100 nM MYC/MAX", "100 nM IRF3", "100 nM JUN", "100 nM JUN/FOS"
    )
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))
groupCellType(project) <- c("BAC")

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- "/data/PRINT/multiScaleFootprinting/data/inVitroFootprints/fragments/all.frags.tsv.gz"
project <- getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T)

# Load background dispersion model
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <-
    readRDS(paste0(projectMainDir,"data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

###################################
# Find candidate TF binding sites #
###################################

# Load unibind ChIP data
unibindCellTypes <- c("HepG2", "K562", "A549", "GM12878")
unibindList <- lapply(
  unibindCellTypes,
  function(cellType){
    readRDS(paste0("/data/PRINT/multiScaleFootprinting/data/shared/unibind/", cellType, "ChIPRanges.rds"))
  }
)
names(unibindList) <- unibindCellTypes

# Merge Unibin data of the same TF across cell types and get candiaate TF binding sites
TFSites <- list()
for(TF in c("CEBPA", "MAX", "MYC", "IRF3", "JUN", "FOS")){
  TFSites[[TF]] <- NULL
  for(cellType in unibindCellTypes){
    if(TF %in% names(unibindList[[cellType]])){
      cellTypeSites <- unibindList[[cellType]][[TF]]
      cellTypeSites$cellType <- cellType
      if(is.null(TFSites[[TF]])){
        TFSites[[TF]] <- cellTypeSites
      }else{
        TFSites[[TF]] <- GenomicRanges::union(TFSites[[TF]], cellTypeSites)
      }
    }
  }
}

###############################
# Visualize a specific region #
###############################

TF <- "MAX"

# Select a region that overlaps with candidate sites
boundRegions <- findOverlaps(TFSites[[TF]], resize(regions, 800, fix = "center"))@to
regionInd <- sample(boundRegions, 1)
regionInd

# Get differential footprints between treated and untreated
groupList <- list("CEBPA" = "CEBPA", "MYC" =  "MYC/MAX", "MAX" = "MYC/MAX", "IRF3" = "IRF3")
pControl <- plotMultiScale(project = project, regionInd = regionInd, lineageGroups = c("Control"), vmax = 3)
pTreated <- plotMultiScale(project = project, regionInd = regionInd, lineageGroups = c(paste0("100 nM ", groupList[[TF]])), vmax = 3)
pDiff <- pTreated
pDiff@matrix <- pTreated@matrix - pControl@matrix

# Label the positions of TFs
labelSites <- F
if(labelSites){
  region <- regions[regionInd]
  regionMotifSites <- subsetByOverlaps(TFSites[[TF]], region)
  relativePos <- start(resize(regionMotifSites, 1, fix = "center")) - start(region) + 1
  relativePos <- relativePos - 100
  relativePos <- relativePos[(relativePos > 0) & (relativePos < 800)]
  for (pos in relativePos){
    pDiff@matrix[, (pos-1):(pos + 1)] <- 2
  }
}
pDiff

plots <- list()
positions <- 100:900
treated <- c(paste0("100 nM ", groupList[[TF]]))
for(condition in c("Control", treated)){
  
  Tn5Insertion <- getFeatureMatrix(
    project = project,
    regionInd = regionInd,
    feature = "insertion",
    groupIDs = condition
  )[, 1]
  
  # Plot Tn5 insertion counts
  barData <- data.frame(x1 = positions - 1 + start(region) - 1, 
                        x2 = positions + 1 + start(region) - 1, 
                        y1 = 0, y2 = Tn5Insertion[positions])
  plots[[condition]] <- ggplot(barData) + 
    geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
              fill = "#7F3F98") +
    scale_y_continuous(expand = c(0, 0)) + ylab(paste0(condition, "\nTn5 insertion")) +
    theme_classic() 
}

# Plot predicted Tn5 bias
Tn5Bias <- regionBias(project)[regionInd, ]
barData <- data.frame(x1 = positions - 1 + start(region) - 1, 
                      x2 = positions + 1 + start(region) - 1, 
                      y1 = 0, y2 = Tn5Bias[positions])
plots[["bias"]] <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#E21F26") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Predicted\nTn5 bias") +
  theme_classic() 

pdf(paste0("/data/PRINT/multiScaleFootprinting/data/inVitroFootprints/plots/", 
    stringr::str_replace_all(as.character(region), "[:-]", "_"), "_tracks.pdf"), width = 12, height = 6)
Reduce("/", plots) + xlab(paste0("Coordinate on ", as.character(seqnames(region))))
dev.off()

pdf(paste0("/data/PRINT/multiScaleFootprinting/data/inVitroFootprints/plots/", 
           stringr::str_replace_all(as.character(region), "[:-]", "_"), "_diff_multiscale.pdf"), width = 12, height = 6)
pDiff@matrix[,-(positions - 100)] = 2
pDiff
dev.off()

pdf(paste0("/data/PRINT/multiScaleFootprinting/data/inVitroFootprints/plots/", 
           stringr::str_replace_all(as.character(region), "[:-]", "_"), "_control_multiscale.pdf"), width = 12, height = 6)
pControl@matrix[,-(positions - 100)] = 2
pControl
dev.off()

pdf(paste0("/data/PRINT/multiScaleFootprinting/data/inVitroFootprints/plots/", 
           stringr::str_replace_all(as.character(region), "[:-]", "_"), "_treated_multiscale.pdf"), width = 12, height = 6)
pTreated@matrix[,-(positions - 100)] = 2
pTreated
dev.off()

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

###############################################
# Get genomic ranges for running footprinting #
###############################################
  
# Initialize a footprintingProject object
projectName <- "yeast"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "sacCer3")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Set the regionRanges slot
regionRanges(project) <- readRDS("../../data/yeast/regionRanges.rds")

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getPrecomputedBias(project, nCores = 24)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load barcodes for each pseudo-bulk
barcodeGroups <- data.frame(
  barcode = "yeast",
  group = 1L
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))
groupCellType(project) <- c(projectName)

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- paste0(projectDataDir, "merged.frags.gz")
#project <- getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T)

# Load background dispersion model
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <-
    readRDS(paste0(projectMainDir,"data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

if(!file.exists("../../data/yeast/project.rds")){
  saveRDS(project, "../../data/yeast/project.rds")
}

####################
# Run footprinting #
####################

# Specify which programs to run
args <- commandArgs(trailingOnly=TRUE)
options <- rep(F, 2)
options[as.integer(args[1])] <- T
runTFBS <- options[1]
runFootprinting <- options[2]

# Run footprinting
if(runFootprinting){
  
  # Set footprint radius
  if(is.na(args[2])){
    footprintRadius <- 20
  }else{
    footprintRadius <- as.integer(args[2])
  }
  
  project <- getFootprints(
    project, 
    mode = as.character(footprintRadius),
    nCores = 16, 
    footprintRadius = footprintRadius,
    flankRadius = footprintRadius)
  
  system(paste0("mkdir ../../data/", projectName, "/footprints/"))
  saveRDS(footprints(project, as.character(footprintRadius)), 
          paste0(projectDataDir, "footprints/", footprintRadius, "bp_footprints.rds"))
}

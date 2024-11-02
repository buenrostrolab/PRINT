# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../../code/utils.R")
source("../../../code/getCounts.R")
source("../../../code/getBias.R")
source("../../../code/getFootprints.R")
source("../../../code/visualization.R")
source("../../../code/getGroupData.R")
source("../../../code/footprintTracking.R")
source("../../../code/getAggregateFootprint.R")
source("../../../code/getTFBS.R")
library(patchwork)

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "saturatedMut/BMMC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regionOrig <- GRanges(c("chr22:19723266-19723650", "chr11:5227022-5227208", "chr11:5249805-5250078",
                        "chr1:155301395-155301864", "chr2:60494940-60495539"))
regions <- resize(regionOrig, 1000, fix = "center")

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

# Load barcodes for each pseudo-bulk
pathToFragGrouping <- paste0(projectDataDir, "barcodeGrouping.txt")
barcodeGroups <- read.table(pathToFragGrouping, header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- paste0(projectMainDir, "data/", projectName, "/merged.fragments.tsv.gz")
#getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = F)

for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <- 
    readRDS(paste0("../../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Load pseudo-bulk metadata
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")
groupUMAP(project) <- as.matrix(cbind(groupInfo$UMAP1, groupInfo$UMAP2))
groupCellType(project) <- groupInfo$cellType
groupPseudoTime(project) <- groupInfo$Pseudotime

# Get color mapping for cell types
cellTypes <- unique(groupInfo$cellType)
cellTypeColors <- groupInfo$color[match(cellTypes, groupInfo$cellType)]
names(cellTypeColors) <- cellTypes

# Save the footprintingProject object with loaded data
projectPath <- paste0(projectDataDir, "footprintingProject.rds")
if(!file.exists(projectPath)){
  saveRDS(project, projectPath)
}

# Get region-by-position matrix of Tn5 insertion
ATACTrks <- getATACTracks(project)@ATACTracks

##################################
# Compare results with MPRA data #
##################################

system("mkdir ../../../data/saturatedMut/plots/")

# Select a region for visualization
regionInd <- 5
region <- resize(regionRanges(project)[regionInd], 800, fix = "center")
startPos <- start(regionOrig[regionInd]) - start(region) + 1
endPos <- end(regionOrig[regionInd]) - start(region) + 1

# Make multi-scale plot
selectedGroups <- which(groupCellType(project) %in% c("HSC/MPP", "CMP", "MEP", "CLP", "GMP", "early-Ery", "late-Ery"))
p <- plotMultiScale(project, regionInd, vmax = 1.5)

# Subset to the region with MPRA data
pdf(paste0("../../../data/saturatedMut/plots/BMMC_region_", regionInd, "_subset_multiscale.pdf"),
    width = 5, height = 5)
Heatmap(p@matrix[, startPos:endPos], 
        col = p@matrix_color_mapping, 
        cluster_rows = F, 
        cluster_columns = F, border = T)
dev.off()

# Plot whole region
pdf(paste0("../../../data/saturatedMut/plots/BMMC_region_", regionInd, "_full_multiscale.pdf"),
    width = 5, height = 5)
p@matrix[, outer(c(startPos,endPos), 0:1, "+")] <- 5
p
dev.off()

# Plot Tn5 insertion
positions <- start(regionOrig[regionInd]):end(regionOrig[regionInd])
plotData <- data.frame(
  x1 = positions - 1, x2 = positions + 1, y1 = 0,
  y2 = ATACTrks[regionInd,startPos:endPos]
)
pIns <- ggplot(plotData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("BAC Tn5\ninsertion") +
  theme_classic() 
pdf(paste0("../../../data/saturatedMut/plots/BMMC_region_", regionInd, "_Tn5_Insertion.pdf"),
    width = 5, height = 5)
pIns
dev.off()

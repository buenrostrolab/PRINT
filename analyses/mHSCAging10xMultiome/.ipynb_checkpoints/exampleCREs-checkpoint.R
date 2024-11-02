# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getGroupData.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getModelPrediction.R")
source("../../code/getTFBS.R")
source("../../code/getSubstructures.R")

###################
# Load input data #
###################

# Save the project object with pre-loaded slots
project <- readRDS("../../data/mHSCAging10xV3/project.rds")

# Get CRE regions
regions <- regionRanges(project)

# Get TSS positions
TSSs <- BuenRTools::mm10TSSRanges
names(TSSs) <- as.character(TSSs$gene_name)

# Load pseudobulk clustering
pbulkClusters <- read.table("../../data/mHSCAging10xV3/pbulkClusters.txt")
cellStates <- rep("", length(pbulkClusters$V1))
cellStates[pbulkClusters$V1 %in% c("Old_1", "Old_3")] <- "Old Mk-biased"
cellStates[pbulkClusters$V1 %in% c("Old_2")] <- "Young-like Old"
cellStates[pbulkClusters$V1 %in% c("Old_4")] <- "Old Multi-lineage"
cellStates[pbulkClusters$V1 %in% c("Young_3")] <- "Young Mk-biased"
cellStates[pbulkClusters$V1 %in% c("Young_1", "Young_2")] <- "Young Multi-lineage"
groupCellType(project) <- cellStates

# Define the order of cell states for plotting
cellStateOrder <- c("Old Mk-biased", "Old Multi-lineage",
                    "Young-like Old", "Young Mk-biased", "Young Multi-lineage")
groupOrder <- order(match(cellStates, cellStateOrder))

system("mkdir ../../data/mHSCAging10xV3/plots/exampleCREs")

# Find the CRE that overlaps with the promoter of interest
gene <- "Psat1"
regionInd <- findOverlaps(regions, TSSs[gene])@from

# Load pseudobulked RNA data
groupRNA(project) <- readRDS("../../data/mHSCAging10xV3/pseudobulkedRNANormed.rds")

####################################
# Visualize multi-scale footprints #
####################################

plots <- list()
for(cellState in unique(groupCellType(project))){
  stateGroups <- groups(project)[groupCellType(project) == cellState]
  p <- plotMultiScale(project, regionInd, lineageGroups = stateGroups)
  pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
             "_", gene, "_", cellState, "_multiscale.pdf"), width = 7, height = 5)
  print(p)
  dev.off()
  plots[[cellState]] <- p
}

# Make color palette for differential multi-scale
footprintColors = brewer.pal(9, "RdBu")
footprintColors[5] <- "white"
footprintColors <- colorRamp2(
  seq(-2, 2,length.out=9),
  colors = rev(footprintColors))

# Plot differential multi-scale
state1 <- "Young-like Old"
state2 <- "Young Multi-lineage"
diffMtx <- plots[[state1]]@matrix - plots[[state2]]@matrix
pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
           "_", gene, "_diff_", state1, "_vs_", state2, "_multiscale.pdf"), width = 7, height = 5)
Heatmap(diffMtx,
        use_raster = TRUE,
        col = footprintColors,
        cluster_rows = F,
        show_row_dend = F,
        cluster_columns = F,
        name = "Footprint\nScore",
        border = TRUE,
        column_title = "Position (bp)",
        column_title_side = "bottom",
        column_names_rot = 0)
dev.off()

########################################
# Visualize TF and nucleosome heatmaps #
########################################

# Plot 50 bp scale footprint scores as a proxy of nucleosomes
pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
           "_", gene, "_50bp_footprint.pdf"), width = 7, height = 5)
plotFeatureHeatmap(
  project = project,
  regionInd = regionInd,
  feature = "footprint",
  rowOrder = groupOrder,
  cellTypeOrder = cellStateOrder,
  gene = gene, 
  footprintRadius = 50)
dev.off()

# Plot TF binding scores (model II)
pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
           "_", gene, "_TF_model_II.pdf"), width = 7, height = 5)
TFBindingModel(project) <- loadTFBSModel("../../data/TFBSPrediction/TFBS_model.h5")
plotFeatureHeatmap(
  project = project,
  regionInd = regionInd,
  feature = "TFBS",
  rowOrder = groupOrder,
  cellTypeOrder = cellStateOrder,
  gene = gene, 
  vmax = 0.6,
)
dev.off()

# Plot TF binding scores (model I)
pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
           "_", gene, "_TF_model_I.pdf"), width = 7, height = 5)
TFBindingModel(project) <- loadTFBSModel("../../data/TFBSPrediction/TFBS_model_cluster_I.h5")
plotFeatureHeatmap(
  project = project,
  regionInd = regionInd,
  rowOrder = groupOrder,
  cellTypeOrder = cellStateOrder,
  feature = "TFBS",
  gene = gene, 
  vmax = 0.5)
dev.off()

# Plot nucleosome scores
pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
           "_", gene, "_nucleosome.pdf"), width = 7, height = 5)
predictionModel(project) <- loadModel("../../data/yeast/nucleosome_model.h5")
plotFeatureHeatmap(
  project = project,
  regionInd = regionInd,
  feature = "nucleosome",
  rowOrder = groupOrder,
  cellTypeOrder = cellStateOrder,
  gene = gene, 
  vmax = 1
)
dev.off()

######################################
# Visualize CRE segmentation results #
######################################

pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
           "_", gene, "_segmentation.pdf"), width = 7, height = 5)
plotSegmentHeatmap(
  project = project,
  regionInd = regionInd
)
dev.off()

#################################
# Visualize ATAC and RNA levels #
#################################

colnames(groupRNA(project)) <- colnames(groupATAC(project))
featureSignal <- list(
  "ATAC" = groupATAC(project)[regionInd, ],
  "RNA" = groupRNA(project)[gene, ]
)

for(feature in names(featureSignal)){
  
  plotData <- data.table::rbindlist(
    lapply(
      unique(groupCellType(project)),
      function(cellState){
        data.frame(
          signal = featureSignal[[feature]][groupCellType(project) == cellState],
          state = cellState
        )
      }
    )
  )
  
  pdf(paste0("../../data/mHSCAging10xV3/plots/exampleCREs/region_", regionInd,
             "_", gene, "_", feature, ".pdf"), width = 10, height = 6)
  print(ggpubr::ggviolin(
    plotData, x = "state", y = "signal", fill = "state", width = 0.5,
    ylab = feature, xlab = "") + 
      ggpubr::stat_compare_means(
        comparisons = list(c("Old Multi-lineage", "Young Multi-lineage"), 
                           c("Old Mk-biased", "Young Mk-biased")),
        label.y = max(plotData$signal) * 1.1))
  dev.off()
  
}


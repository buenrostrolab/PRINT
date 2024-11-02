# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getGroupData.R")
source("../../code/getTFBS.R")

###################
# Load input data #
###################

# Read footprint project
project <- readRDS("../../data/BMMCDonorCellType/project.rds")

# Get color mapping for cell types
cellTypeColors <- c("#0E4421","#00AF99","#9AD9E9","#A896C8","#FAA31B","#F28238",
                    "#F26625","#D13C27","#46A147","#EF3741","#CC1F39","#901838",
                    "#A896C8","#B280B9","#C757A1","#854199","#005F9F","#1481C4","#492264","#67B545")
names(cellTypeColors) <- c("HSC/MPP", "LMPP","CLP", "pro/pre-B","GMP", "CD14mono", "CD16mono", "DC",
                           "CMP", "MEP", "early-Ery", "late-Ery", "NaiveB", "MemoryB", "plasmaB", 
                           "pDC", "CD4", "CD8", "NK", "Baso")

# Find the region overlapping with the MPRA region
regions <- regionRanges(project)
regionInd <- findOverlaps(GRanges("chr11:5227050-5227150"), regions)@to

# Get donor labels
donors <- stringr::str_split_fixed(groups(project), "_", 2)[, 1] 

# Find the region with MPRA data
MPRARegion <- GRanges("chr11:5227022-5227208")
region <- resize(regions[regionInd], 800, fix = "center")
MPRAStart <- start(MPRARegion) - start(region) + 1
MPRAEnd <- end(MPRARegion) - start(region) + 1

####################################
# Visualize multi-scale footprints #
####################################

system("mkdir -p ../../data/BMMCDonorCellType/plots/HBB/")
multiScalePlots <- list()
for(cellType in c("early-Ery", "late-Ery")){
  
  # Plot multi-scale footprints in a specific cell type
  pdf(paste0("../../data/BMMCDonorCellType/plots/HBB/region_", regionInd, "_multiscale_", cellType, ".pdf"),
      width = 12, height = 5)
  p <- plotMultiScale(
    project = project,
    regionInd = regionInd,
    lineageGroups = groups(project)[groupCellType(project) %in% cellType],
    vmax = 1.5
  )
  print(p)
  multiScalePlots[[cellType]] <- p
  dev.off()
  
  # Make the same plot but this time label the boundary of the MPRA region
  p@matrix[, MPRAStart] <- 1
  p@matrix[, MPRAEnd] <- 1
  pdf(paste0("../../data/BMMCDonorCellType/plots/HBB/region_", regionInd, "_multiscale_", cellType, "_MPRA_range.pdf"),
      width = 12, height = 5)
  print(p)
  dev.off()
}

# Calculate differeces in footprints
diffMtx <- multiScalePlots$`late-Ery`@matrix- multiScalePlots$`early-Ery`@matrix
diffMtx <- diffMtx[, MPRAStart:MPRAEnd]

# Make color palette for differential multi-scale
footprintColors = brewer.pal(9, "RdBu")
footprintColors[5] <- "white"
footprintColors <- colorRamp2(
  seq(-2, 2,length.out=9),
  colors = rev(footprintColors))

# Plot differential multi-scale
pdf(paste0("../../data/BMMCDonorCellType/plots/HBB/region_", regionInd, "_diff_multiscale.pdf"), 
    width = 15, height = 5)
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

####################################
# Compare footrpints across donors #
####################################

# Calculate 20bp-diameter footprints across donors
plotData <- NULL
for(cellType in c("early-Ery", "late-Ery")){
  
  fpTracks <- plotFeatureHeatmap(
    project = project,
    regionInd = regionInd,
    lineageGroups = groups(project)[groupCellType(project) %in% cellType],
    footprintRadius = 10,
  )@matrix[, MPRAStart:MPRAEnd]
  
  plotData <- rbind(
    plotData, 
    data.frame(
      position = start(MPRARegion):end(MPRARegion),
      footprint = colMeans(fpTracks),
      sd = colSds(fpTracks),
      celltype = cellType
    ))
  
}

# Visualize 20bp-diameter footprints across donors
pdf(paste0("../../data/BMMCDonorCellType/plots/HBB/region_", regionInd, "_10bp_footprint.pdf"),
    width = 15, height = 5)
ggplot(data = plotData) +
  geom_line(aes(x = position, y = footprint, color = celltype)) +
  geom_ribbon(
    aes(x = position, ymin = footprint - sd, ymax = footprint + sd, fill = celltype),
    alpha = 0.2)+
  ylab("Footprint score (20 bp scale)") +
  #scale_color_manual(values = c("early-Ery" = "#CC1F39", "late-Ery" = "#901838")) +
  #scale_fill_manual(values = c("early-Ery" = "#CC1F39", "late-Ery" = "#901838")) +  
  xlab("Position") +
  theme_classic()
dev.off()
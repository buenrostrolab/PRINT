# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(ggplot2)
library(Seurat)
library(GenomicRanges)
library(BuenColors)
library(ArchR)
source("../../code/utils.R")

###################
# Load input data #
###################

system("mkdir ../../data/BMMC/plots/QC")

# Load CRE ranges
regions <- readRDS("../../data/BMMC/regionRanges.rds")

# Load single-cell ATAC data
scATACSeurat <- readRDS("../../data/BMMC/atac.se.rds")

# Load single-cell RNA Seurat object
scRNASeurat <- readRDS("../../data/BMMC/scRNA.rds")
scRNASeurat <- scRNASeurat[, colnames(scATACSeurat)]
scRNASeurat$cellType <- scATACSeurat$cistopic.assign.l2.rank

# Load barcodes for each pseudobulk
barcodeGroups <- read.table("../../data/BMMC/barcodeGrouping.txt",
                            header = T)

# Load pseudo-bulk metadata
groupInfo <- read.table("../../data/BMMC/groupInfo.txt", header = T,
                        comment.char = "")

############################
# Visualize fragment sizes #
############################

# Load fragment coordinates
frags <- data.table::fread("../../data/BMMC/merged.fragments.tsv.gz", 
                           sep = "\t", showProgress = TRUE) 

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

pdf("../../data/BMMC/plots/QC/fragSizes.pdf",
    width = 7, height = 5.5)
ggplot(plotData) +
  geom_line(aes(x = length, y = freq), size = 1, color = "#2171B5") +
  theme_classic()
dev.off()

###########################################
# Visualize marker RNA for each cell type #
###########################################

dotPlot <- function(dataMtx, markers, groupLabels,
                    scaleFeature = T,
                    palette = NULL,
                    vmax = NULL,
                    vmin = NULL){
  
  dataMtx <- dataMtx[markers, ]
  
  groups <- unique(groupLabels)
  pctMtx <- pbmcapply::pbmcmapply(
    function(group){
      groupData <- dataMtx[markers, groupLabels == group]
      percentage <- rowMeans(groupData > 0)
    },
    groups,
    mc.cores = 8
  ) * 100
  colnames(pctMtx) <- groups
  
  expMtx <- pbmcapply::pbmcmapply(
    function(group){
      groupData <- dataMtx[markers, groupLabels == group]
      exp <- rowMeans(groupData)
    },
    groups,
    mc.cores = 8
  )
  colnames(expMtx) <- groups
  
  if(scaleFeature){
    expMtx <- expMtx / rowMeans(expMtx)
  }
  
  plotData <- as.data.frame(expand.grid(markers, groups))
  colnames(plotData) <- c("markers", "groups")
  
  plotData$exp <- sapply(1:dim(plotData)[1], function(i){expMtx[as.character(plotData$markers)[i], 
                                                                as.character(plotData$groups)[i]]})
  
  plotData$percent <- sapply(1:dim(plotData)[1], function(i){pctMtx[as.character(plotData$markers)[i], 
                                                                    as.character(plotData$groups)[i]]})
  
  if(!is.null(vmax)){ plotData$exp <- pmin(plotData$exp, vmax) }
  if(!is.null(vmin)){ plotData$exp <- pax(plotData$exp, vmin) }
  
  if(is.null(palette)){palette <- jdb_palette("brewer_blue")}
  
  p <- ggplot(plotData) +
    geom_point(aes(x = markers, y = groups, color = exp, size = percent)) +
    scale_color_gradientn(colors = palette) +
    theme_classic() + 
    RotatedAxis()
  
  p
}

##################
# TSS enrichment #
##################

TSSEnrichment <- getTSSEnrichment(frags, "hg38")

pdf("../../data/BMMC/plots/QC/TSSPlot.pdf",
    width = 7, height = 5.5)
ggplot(plotData) +
  geom_line(aes(x = position, y = normValue), size = 1) +
  theme_classic()
dev.off()

#########################################################
# Visualize marker gene scores (RNA) for each cell type #
#########################################################

scRNACounts <- scRNASeurat@assays$RNA@counts
scRNACounts <- centerCounts(scRNACounts, chunkSize = 1e5)

markers <- c("CD34", "TCL1A", "CD79A", "CD4", "CD3D", "GIMAP5",
             "CD8A", "CD8B", "S100A4", "CCR7", "CCL5", "NKG7",
             "GNLY", "TSPAN13", "JCHAIN", "CD14", "LYZ", "S100A9", "IL1B",
             "HBA2", "GNG11", "CD14", "CD16", "IL3RA", "FCER1A")
markers <- intersect(markers, rownames(scRNASeurat))

pdf("../../data/BMMC/plots/QC/RNADotPlot.pdf",
    width = 7, height = 5.5)
p <- dotPlot(scRNACounts, sort(markers), scRNASeurat$cellType,
             scaleFeature = T, vmax = 5,
             palette = jdb_palette("brewer_red"))
p
dev.off()

##########################################################
# Visualize marker gene scores (ATAC) for each cell type #
##########################################################

# Create scATAC SummarizedExperiment object
scATACSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(count = scATACSeurat@assays@data$counts),
  rowRanges = regions
)

# Calculate gene scores
geneScores <- getGeneScoresFromPeaks(scATACSE, 
                                     genome = "hg38",
                                     TSSwindow = 10000, 
                                     getWeightsOnly = FALSE)

# Normalize gene scores
geneScores <- centerCounts(geneScores, chunkSize = 1e5)

pdf("../../data/BMMC/plots/QC/ATACDotPlot.pdf",
    width = 7, height = 5.5)
markers <- sort(intersect(rownames(geneScores), markers))
p <- dotPlot(as.matrix(geneScores[markers,]), markers, scRNASeurat$cellType,
             scaleFeature = T, vmax = 5)
p
dev.off()

############################
# Visualize donors on UMAP #
############################

plotData <- data.frame(
  UMAP1 = scATACSeurat$cisTopic.umap1,
  UMAP2 = scATACSeurat$cisTopic.umap2,
  Donor = scATACSeurat$donor
)
donorColors <- c("#0E4421","#00AF99","#9AD9E9","#A896C8","#FAA31B","#F28238", "#F26625")
names(donorColors) <- sort(unique(scATACSeurat$donor))

system("mkdir ../../data/BMMC/plots/QC")
png("../../data/BMMC/plots/QC/UMAPDonor.png",
    width = 1500, height = 1300)
ggplot(plotData[sample(1:dim(plotData)[1]),]) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = Donor), size = 0.001) +
  scale_color_manual(values = donorColors) +
  theme_classic()
dev.off()

#################################
# Visualize pseudobulks on UMAP #
#################################

# Assign a color to each cell type
cellTypeColors <- c("#0E4421","#00AF99","#9AD9E9","#A896C8","#FAA31B","#F28238",
                    "#F26625","#D13C27","#46A147","#EF3741","#CC1F39","#901838",
                    "#A896C8","#B280B9","#C757A1","#854199","#005F9F","#1481C4","#492264","#67B545")
names(cellTypeColors) <- c("HSC/MPP", "LMPP","CLP", "pro/pre-B","GMP", "CD14mono", "CD16mono", "DC",
                           "CMP", "MEP", "early-Ery", "late-Ery", "NaiveB", "MemoryB", "plasmaB", 
                           "pDC", "CD4", "CD8", "NK", "Baso")

plotData <- data.frame(
  UMAP1 = scATACSeurat$cisTopic.umap1,
  UMAP2 = scATACSeurat$cisTopic.umap2,
  cellType = scATACSeurat$cistopic.assign.l2.rank
)

groupInd <- 825
groupBarcodes <- barcodeGroups$barcode[barcodeGroups$group == groupInd]
png(paste0("../../data/BMMC/plots/QC/Pseudobulk_", groupInd, "_", groupInfo$cellType[groupInd], ".png"),
    width = 1500, height = 1300)
ggplot() +
  geom_point(data = plotData, aes(x = UMAP1, y = UMAP2, color = cellType),
             size = 0.001) +
  scale_colour_manual(values = cellTypeColors) +
  geom_point(data = plotData[colnames(scATACSeurat) %in% groupBarcodes, ], aes(x = UMAP1, y = UMAP2),
             size = 0.1, color = "Black") +
  ggtitle(paste0("Pseudobulk ", groupInd, " ", groupInfo$cellType[groupInd])) +
  theme_classic()
dev.off()

# Visualize pseudo-bulk centers
centerBarcodes <- pbmcapply::pbmcmapply(
  function(groupInd){
    barcodeGroups[barcodeGroups$group == groupInd, ]$barcode[1]
  },
  unique(barcodeGroups$group),
  mc.cores = 8
)
png(paste0("../../data/BMMC/plots/QC/PseudobulkCenters.png"),
    width = 1500, height = 1300)
ggplot() +
  geom_point(data = plotData, aes(x = UMAP1, y = UMAP2, color = cellType),
             size = 0.001) +
  scale_colour_manual(values = cellTypeColors) +
  geom_point(data = plotData[colnames(scATACSeurat) %in% centerBarcodes, ], aes(x = UMAP1, y = UMAP2),
             size = 3, color = "Black") +
  theme_classic()
dev.off()

# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getGroupData.R")
library(BuenColors)
library(ggrepel)
library(Seurat)

###################
# Load input data #
###################

# Load single cell ATAC data
scATAC <- readRDS("../../data/mHSCAging10xV3/scATACSeurat.rds")

# Load single cell RNA data
scRNA <- readRDS("../../data/mHSCAging10xV3/scRNA.rds")

# Load LSI embedding of single cells
cellEmbedding <- readRDS("../../data/mHSCAging10xV3/LSIEmbedding.rds")

# Load barcodes for each pseudobulk
barcodeGroups <- read.table("../../data/mHSCAging10xV3/barcodeGrouping.txt",
                            header = T)

# Load pseudobulk centers
pseudobulkCenters <- read.table("../../data/mHSCAging10xV3/pseudobulkCenters.txt",
                                header = T)
rownames(pseudobulkCenters) <- pseudobulkCenters$group

# Load differential RNA testing results
diffRNA <- read.table("../../data/mHSCAging10xV3/diffRNA.tsv")

##################################
# Normalize single cell RNA data #
##################################

# Only keep cells in the HSC cluster 
HSCBarcodes <- colnames(scATAC)[scATAC$seurat_clusters %in% c(0, 6, 10)]
scRNAFilt <- scRNA[, intersect(HSCBarcodes, colnames(scRNA)),]

# Remove top highly expressed genes
geneExpLevel <- rowMeans(scRNAFilt@assays$RNA@counts)
scRNAFilt <- scRNAFilt[geneExpLevel < quantile(geneExpLevel, 0.95),]

# Blacklist ribosomal genes and mito genes
blacklist <- c(rownames(scRNAFilt)[stringr::str_detect(rownames(scRNAFilt), "Rpl")],
               rownames(scRNAFilt)[stringr::str_detect(rownames(scRNAFilt), "Rps")],
               rownames(scRNAFilt)[stringr::str_detect(rownames(scRNAFilt), "mt-")])
scRNAFilt <- scRNAFilt[which(!(rownames(scRNAFilt) %in% blacklist)),]

# Filter cells with low RNA reads
cellFilter <- colSums(scRNAFilt@assays$RNA@counts) > 100
scRNAFilt <- scRNAFilt[, cellFilter]

# Normalize RNA data
scRNAFilt <- SCTransform(
  scRNAFilt, verbose = T, 
  variable.features.n = 5000, 
  return.only.var.genes = F
)

# Normalized scRNA data
scRNAMat <- scRNAFilt[["SCT"]]@scale.data

#########################################
# Visualize gene scores of marker genes #
#########################################

scATACSE <- SummarizedExperiment(
  assays=list(counts=scATAC@assays$ATAC@counts), 
  rowRanges = GRanges(rownames(scATAC)))

# Calculate gene scores
geneScores <- getGeneScoresFromPeaks(scATACSE, 
                                     genome = "mm10", 
                                     TSSwindow = 10000, 
                                     getWeightsOnly = FALSE)

# Normalize gene scores
geneScores <- centerCounts(geneScores, chunkSize = 1e4)

# Marker genes to plot
markers <- c("Cd3d", "Gata1", "Cd79a", "Selp",
             "Nupr1", "Elane", "Hlf", "Clu", "Fos",
             "Stat3", "Gata2")
markers <- intersect(markers, rownames(geneScores))

# Smooth the gene scores of marker genes
KNN <- FNN::get.knn(scATAC@reductions$lsi@cell.embeddings[, 2:20], k = 500)$nn.index
geneScoresFilt <- geneScores[markers, ]
smoothedMarker <- pbmcapply::pbmcmapply(
  function(cellInd){
    rowMeans(geneScoresFilt[, KNN[cellInd, ]])
  },
  1:dim(scATAC)[2],
  mc.cores = 12
)

# Visualize data on UMAP
marker <- "Cd3d"
pltData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2],
  marker = smoothedMarker[marker, ]
)

p <- ggplot(pltData[sample(1:dim(pltData)[1]),]) +
  geom_point(aes(x =  UMAP1, y = UMAP2, color = marker), size = 0.5, alpha = 1) +
  geom_rect(aes(xmin = -5, xmax = 5, ymin = -4.5, ymax = 4), 
            fill = "blue", alpha = 0, color = "black") +
  scale_color_gradientn(colors = c("#FFEDB0", "#FFDF5F", "#FEC510", "#FA8E24", 
                                   "#F14C2B", "#DA2828", "#BE2222", "#A31D1D")) +
  theme_classic() +
  ggtitle(marker)

pdf(paste0("../../data/mHSCAging10xV3/plots/UMAP_marker_", marker, ".pdf"), width = 7, height = 7)
p
dev.off()

png(paste0("../../data/mHSCAging10xV3/plots/UMAP_marker_", marker, ".png"), 
    width = 1000, height = 1000, units = "px")
p
dev.off()

# Visualize data on HSC UMAP
marker <- "Clu"
pltData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2],
  marker = smoothedMarker[marker, ]
)
pltData <- pltData[(pltData$UMAP1 > -5) & (pltData$UMAP1 < 0) &
                     (pltData$UMAP2 > -1) & (pltData$UMAP2 < 4), ]

p <- ggplot(pltData[sample(1:dim(pltData)[1]),]) +
  geom_point(aes(x =  UMAP1, y = UMAP2, color = marker), size = 0.5, alpha = 1) +
  scale_color_gradientn(colors = c("#FFEDB0", "#FFDF5F", "#FEC510", "#FA8E24", 
                                   "#F14C2B", "#DA2828", "#BE2222", "#A31D1D")) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic() +
  ggtitle(marker)

pdf(paste0("../../data/mHSCAging10xV3/plots/UMAP_marker_", marker, ".pdf"), width = 7, height = 7)
p
dev.off()

png(paste0("../../data/mHSCAging10xV3/plots/UMAP_marker_", marker, ".png"), 
    width = 700, height = 700, units = "px")
p
dev.off()

#################################
# Visualize RNA of marker genes #
#################################

markers <- intersect(markers, rownames(scRNAMat))

# Smooth the gene scores of marker genes
embeddingFilt <- scATAC@reductions$lsi@cell.embeddings[colnames(scRNAMat), ]
KNN <- FNN::get.knn(embeddingFilt[, 2:20], k = 500)$nn.index
RNAFilt <- scRNAMat[markers, ]
smoothedMarker <- pbmcapply::pbmcmapply(
  function(cellInd){
    rowMeans(RNAFilt[, KNN[cellInd, ]])
  },
  1:dim(scRNAMat)[2],
  mc.cores = 12
)

# Visualize data on HSC UMAP
marker <- "Selp"
pltData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[colnames(scRNAMat), 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[colnames(scRNAMat), 2],
  marker = smoothedMarker[marker, ]
)
pltData <- pltData[(pltData$UMAP1 > -5) & (pltData$UMAP1 < 0) &
                     (pltData$UMAP2 > -1) & (pltData$UMAP2 < 4), ]

p <- ggplot(pltData[sample(1:dim(pltData)[1]),]) +
  geom_point(aes(x =  UMAP1, y = UMAP2, color = marker), size = 0.5, alpha = 1) +
  scale_color_gradientn(colors = c("#FFEDB0", "#FFDF5F", "#FEC510", "#FA8E24", 
                                   "#F14C2B", "#DA2828", "#BE2222", "#A31D1D")) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic() +
  ggtitle(marker)
p

pdf(paste0("../../data/mHSCAging10xV3/plots/UMAP_marker_", marker, "_RNA.pdf"), width = 7, height = 7)
p
dev.off()

png(paste0("../../data/mHSCAging10xV3/plots/UMAP_marker_", marker, "_RNA.png"), 
    width = 700, height = 700, units = "px")
p
dev.off()

#############################################
# Visualize differential expression results #
#############################################

diffRNA$logP <- -log10(diffRNA$pvalue)
diffRNA$logP[diffRNA$logP > 100] <- 100
selected <- rownames(diffRNA) %in% c("Clu", "Selp", "Nupr1", "Gstm2", "Fos", "Gata2", "Jun",
                                     "Igf1")
labels <- rep("", dim(diffRNA)[1])
labels[selected] <- rownames(diffRNA)[selected]
pdf(paste0("../../data/mHSCAging10xV3/plots/diffRNAVolcano.pdf"), width = 7, height = 7)
ggplot() +
  geom_point(data = diffRNA[diffRNA$log2FoldChange > 1.5,], 
             aes(x = log2FoldChange, y = logP), color = "red") +
  geom_point(data = diffRNA[abs(diffRNA$log2FoldChange) < 1.5,], 
             aes(x = log2FoldChange, y = logP), color = "grey") +
  geom_point(data = diffRNA[diffRNA$log2FoldChange < -1.5,], 
             aes(x = log2FoldChange, y = logP), color = "blue") +
  ylim(0, 100) +
  geom_text_repel(data = diffRNA,
                  aes(x= log2FoldChange, y = logP), 
                  label = labels, color = "black", 
                  max.overlaps = 1e5) +
  theme_classic()
dev.off()

##################################
# Visualize cell type annotation #
##################################

plotData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2],
  cellType = scATAC$cellType
)

cellTypeColors <- c("#0E4421", "#FAA31B", "#AD8E5C", "#EF3741", "#F28238", "#1481C4", "#8F7853", 
                    "#67B545", "#C757A1", "#A896C8")
names(cellTypeColors) <- c("HSC","GMP","Neu-GMP","MEP","Mono-GMP","T cells", 
                           "Granulocytes","Baso/Mast","Plasma cells", "B cells")

png("../../data/mHSCAging10xV3/plots/UMAPCellTypes.png",
    width = 1500, height = 1500)
ggplot(plotData[sample(1:dim(plotData)[1]),]) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = cellType), size = 1, alpha = 1) +
  scale_color_manual(values = cellTypeColors) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic()
dev.off()

###################################
# Visualize donor age on the UMAP #
###################################

age <- scATAC$age == "Old"

plotData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2],
  age = age
)
plotData <- plotData[(plotData$UMAP1 > -5) & (plotData$UMAP1 < 0) &
                       (plotData$UMAP2 > -1) & (plotData$UMAP2 < 4), ]

png("../../data/mHSCAging10xV3/plots/UMAP_donor_age.png",
    width = 1500, height = 1500)
ggplot(plotData[sample(1:dim(plotData)[1]),]) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = age), size = 1) +
  scale_color_manual(values = c("#3361A5", "#A31D1D")) +
  theme_classic()
dev.off()

########################
# Visualize HSC states #
########################

groupLabels <- as.character(scATAC$seurat_clusters)
groupLabels[!(groupLabels %in% c("0", "6", "10"))] <- "others"
colors <- c("#00441B", "#4154A5", "#FFA300", "grey")
names(colors) <- c("0", "6", "10", "others")

plotData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2],
  groupLabels = groupLabels
)

png("../../data/mHSCAging10xV3/plots/threeHSCStates.png",
    width = 1500, height = 1500)
ggplot(plotData[sample(1:dim(plotData)[1]),]) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = groupLabels), size = 1, alpha = 0.5) +
  scale_color_manual(values = colors) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic()
dev.off()

#########################
# Re-number pseudobulks #
#########################

pseudobulkNames <- gtools::mixedsort(unique(barcodeGroups$group))
pseudobulkCenters <- pseudobulkCenters[pseudobulkNames,]
pseudobulkAges <- sapply(pseudobulkCenters$barcode, function(x){strsplit(x, "-")[[1]][2]})
pseudobulkLSI <- cellEmbedding[pseudobulkCenters$barcode, 2]

# Within each age-group, re-order pseudobulks by their LSI-2 coordinate
# We don't use the first LSI because it highly correlates with depth
reOrder <- c(order(pseudobulkLSI[pseudobulkAges == "Old"]), 
             order(pseudobulkLSI[pseudobulkAges == "Young"]) + sum(pseudobulkAges == "Old"))

# Re-label pseudobulk membership of single cells
barcodeGroups <- data.table::rbindlist(lapply(
  1:length(reOrder),
  function(i){
    x <- barcodeGroups[barcodeGroups$group == pseudobulkNames[reOrder[i]],]
    x$group <- pseudobulkNames[i]
    x
  }
))

# Re-number pseudobulk centers
pseudobulkCenters <- pseudobulkCenters[pseudobulkNames,]
pseudobulkCenters <- pseudobulkCenters[reOrder, ]
pseudobulkCenters$group <- pseudobulkNames

################################
# Visualize pseudobulk centers #
################################

plotData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2]
)
rownames(plotData) <- colnames(scATAC)

plotDataCenter <- plotData[pseudobulkCenters$barcode, ]
oldLabels <- sapply(pseudobulkCenters$group[pseudobulkAges == "Old"],
                    function(s){strsplit(s, "_")[[1]][2]})
youngLabels <- sapply(pseudobulkCenters$group[pseudobulkAges == "Young"],
                      function(s){strsplit(s, "_")[[1]][2]})
png("../../data/mHSCAging10xV3/plots/pseudobulkCenters.png",
    width = 1500, height = 1500)
ggplot(plotData) +
  geom_point(aes(x = UMAP1, y = UMAP2), 
             color = "grey", size = 0.5, alpha = 0.3) +
  geom_point(data = plotDataCenter[pseudobulkAges == "Young", ],
             aes(x = UMAP1, y = UMAP2), color = "#3361A5", size = 5) +
  geom_text_repel(data = plotDataCenter[pseudobulkAges == "Young", ],
                  aes(x = UMAP1, y = UMAP2), 
                  label = youngLabels, color = "#3361A5",
                  max.overlaps = 100, size = 12) +
  geom_point(data = plotDataCenter[pseudobulkAges == "Old", ],
             aes(x = UMAP1, y = UMAP2), color = "#A31D1D", size = 5) +
  geom_text_repel(data = plotDataCenter[pseudobulkAges == "Old", ],
                  aes(x = UMAP1, y = UMAP2), 
                  label = oldLabels, color = "#A31D1D",
                  max.overlaps = 100, size = 12) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic()  
dev.off()

################################
# Visualize pseudobulk members #
################################

plotData <- data.frame(
  UMAP1 = scATAC@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATAC@reductions$umap@cell.embeddings[, 2]
)

system("mkdir ../../data/mHSCAging10xV3/plots")
for(ID in unique(barcodeGroups$group)){
  p <- ggplot(plotData[sample(1:dim(plotData)[1]), ]) +
    geom_point(aes(x = UMAP1, y = UMAP2), color = "grey", size = 0.5) +
    geom_point(data = plotData[barcodeGroups$barcode[barcodeGroups$group == ID],],
               aes(x = UMAP1, y = UMAP2), color = "black", size = 1) +
    geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
              fill = "blue", alpha = 0, color = "black") +
    ggtitle(ID) +
    theme_classic()
  
  # Save plot
  png(paste0("../../data/mHSCAging10xV3/plots/pseudobulk_", ID, ".png"),
      width = 1500, height = 1500)
  print(p)
  dev.off()
}

# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(ggplot2)
library(ggrepel)
library(BuenColors)
source("../../code/getPseudoBulks.R")

###################
# Load input data #
###################

# Load single cell ATAC
scATACSeurat <- readRDS("../../data/mHSCAging10xMultiome/scATACSeurat.rds")
seuratClusters <- scATACSeurat$seurat_clusters
cellType <- scATACSeurat$group

# Get LSI embedding
cellEmbedding <- scATACSeurat@reductions$lsi@cell.embeddings[, 2:20]

# Only keep HSCs
HSCfilter <- (cellType == "HSC") & (seuratClusters %in% c(0,6,10))
cellEmbedding <- cellEmbedding[HSCfilter, ]
scATACSeurat <- scATACSeurat[, HSCfilter]

# Get read count per cell
scATACCounts <- scATACSeurat@assays$ATAC@counts
cellCounts <- colSums(scATACCounts)

###########################################
# Calculate geodesic distance among cells #
###########################################

# Find KNN of cells
nCells <- dim(cellEmbedding)[1]
LSIKNN <- FNN::get.knn(cellEmbedding, k = 3)$nn.index

# Construct adjacency matrix of KNN graph
KNNAdjMatrix <- as.matrix(data.table::rbindlist(pbmcapply::pbmclapply(
  1:nCells, 
  function(i){
    as.data.frame(cbind(i, LSIKNN[i, ]))
  }
)))
KNNAdjMatrix <- sparseMatrix(
  i = KNNAdjMatrix[, 1], j = KNNAdjMatrix[, 2], x = 1)
KNNAdjMatrix <- as.matrix(KNNAdjMatrix)
KNNAdjMatrix <- pmin(KNNAdjMatrix + t(KNNAdjMatrix), 1)

# Construct KNN graph object
KNNGraph <- igraph::graph_from_adjacency_matrix(KNNAdjMatrix, mode = "undirected")

# Compute shortest paths between all cell pairs. We checked that this returns a symmetric matrix
# The length of this path is the geodesic distance between two cells
cellDist <- igraph::shortest.paths(KNNGraph)
rownames(cellDist) <- rownames(cellEmbedding)

########################################
# Load and re-order pseudobulk centers #
########################################

# Load SEACells results
SEACells <- read.table("../../data/mHSCAging10xMultiome/SEACells.tsv",
                       sep = "\t")
colnames(SEACells) <- c("barcode", "center")
centerCells <- unique(SEACells$center)
centerCells <- data.frame(
  barcode = centerCells,
  age = sapply(centerCells, function(x){strsplit(x, "-")[[1]][2]})
)

# Order pseudobulk centers
for(ageGroup in c("Young", "Old")){
  ageInds <- centerCells$age == ageGroup
  centerLSI <- cellEmbedding[centerCells$barcode[ageInds], 1]
  centerCells[ageInds, ] <- centerCells[ageInds, ][order(centerLSI), ]
}

########################
# Pseudobulk the cells #
########################

barcodeGroups <- NULL
pseudobulkCenters <- NULL
countThreshold <- 5e6
for(ageGroup in c("Old", "Young")){
  
  # Select HSC cells in one age group
  mask <- scATACSeurat$age == ageGroup
  centerBarcodes <- centerCells[centerCells$age == ageGroup,]$barcode
  
  # For each center cell, rank all cells in the age group by distance in LSI space
  groupMembers <- t(pbmcapply::pbmcmapply(
    function(barcode){
      order(cellDist[barcode, mask])
    },
    centerBarcodes
  ))
  
  # Generate pseudobulks by keep adding the next nearest cell from the center cell
  # Stop until we reach a certain total number of reads
  scBarcodes <- rownames(cellEmbedding)[mask]
  barcodeGp <- pbmcapply::pbmclapply(
    1:dim(groupMembers)[1],
    function(groupInd){
      scCounts <- cellCounts[mask][groupMembers[groupInd,]]
      nKeptCell <- min(which(cumsum(scCounts) > countThreshold))
      data.frame(barcode = scBarcodes[groupMembers[groupInd, 1:nKeptCell]],
                 group = paste0(ageGroup, "_", groupInd))
    }
  )
  barcodeGp <- data.table::rbindlist(barcodeGp)
  barcodeGroups <- rbind(barcodeGroups, barcodeGp)
  
  centerBarcodes <- data.frame(
    barcode = centerBarcodes,
    group = paste0(ageGroup, "_", 1:dim(groupMembers)[1])
  )
  pseudobulkCenters <- rbind(pseudobulkCenters, centerBarcodes)
}
groupIDs <- unique(barcodeGroups$group)

write.table(barcodeGroups, "../../data/mHSCAging10xMultiome/barcodeGrouping.txt",
            row.names = F, sep = "\t", quote = F)
write.table(pseudobulkCenters, "../../data/mHSCAging10xMultiome/pseudobulkCenters.txt",
            row.names = F, sep = "\t", quote = F)

################################
# Visualize pseudobulk centers #
################################

plotData <- data.frame(
  UMAP1 = scATACSeurat@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATACSeurat@reductions$umap@cell.embeddings[, 2]
)
rownames(plotData) <- colnames(scATACSeurat)
pseudobulkAges <- stringr::str_split_fixed(pseudobulkCenters$group, "_", 2)[, 1]

plotDataCenter <- plotData[pseudobulkCenters$barcode, ]
oldLabels <- sapply(pseudobulkCenters$group[pseudobulkAges == "Old"],
                    function(s){strsplit(s, "_")[[1]][2]})
youngLabels <- sapply(pseudobulkCenters$group[pseudobulkAges == "Young"],
                      function(s){strsplit(s, "_")[[1]][2]})
system("mkdir ../../data/mHSCAging10xMultiome/plots/")
png("../../data/mHSCAging10xMultiome/plots/pseudobulkCenters.png",
    width = 1000, height = 1000)
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
  UMAP1 = scATACSeurat@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATACSeurat@reductions$umap@cell.embeddings[, 2]
)

system("mkdir ../../data/mHSCAging10xMultiome/plots/pseudobulkMembers/")
for(ID in unique(barcodeGroups$group)){
  p <- ggplot(plotData[sample(1:dim(plotData)[1]), ]) +
    geom_point(aes(x = UMAP1, y = UMAP2), color = "grey", size = 1) +
    geom_point(data = plotData[barcodeGroups$barcode[barcodeGroups$group == ID],],
               aes(x = UMAP1, y = UMAP2), color = "black", size = 1) +
    geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
              fill = "blue", alpha = 0, color = "black") +
    ggtitle(ID) +
    theme_classic()
  
  # Save plot
  png(paste0("../../data/mHSCAging10xMultiome/plots/pseudobulkMembers/pseudobulk_", ID, ".png"),
      width = 1000, height = 1000)
  print(p)
  dev.off()
}

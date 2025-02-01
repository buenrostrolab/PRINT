# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getAggregateFootprint.R")
source("../../code/getTFBS.R")
source("../../code/getGroupData.R")
source("../../code/getSubstructures.R")
source("../../code/getGeneCorr.R")

set.seed(42)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(ggpubr)

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "BMMC"
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")

# Load the footprintingProject object
project <- readRDS(paste0(projectDataDir, "footprintingProject.rds"))

# Get color mapping for cell types
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")
cellTypes <- unique(groupInfo$cellType)
cellTypeColors <- groupInfo$color[match(cellTypes, groupInfo$cellType)]
names(cellTypeColors) <- cellTypes

# Get subCRE-by-pseudobulk SummarizedExperiment object
subCREPath <- "../../data/BMMC/substructureSE.rds"
if(file.exists(subCREPath)){
  subCRESE <- readRDS("../../data/BMMC/substructureSE.rds")
}else{
  subCRESE <- getSubstructureSE(project)
  saveRDS(subCRESE, subCREPath)
}

# Filter low signal subCREs
subCRESignal <- rowMaxs(subCRESE@assays@data$substructures)
subCREFilter <- which(subCRESignal > 0.5)
subCRESE <- subCRESE[subCREFilter, ]
subCRERanges <- rowRanges(subCRESE)

# Load motif match positions
motifPositions <- readRDS("../../data/BMMC/motifPositions.rds")

###########################################################
# Only keep CREs with more than one non-isolated sub-CREs #
###########################################################

# First remove isoated sub-CREs
# Extend all sub-CREs by 1 bp on each side and see whether they overlap with the next sub-CRE
paddedSubCRERanges <- resize(subCRERanges, width(subCRERanges) + 2, fix = "center")
adjacencyFilter <- which(table(findOverlaps(paddedSubCRERanges, paddedSubCRERanges)@to) > 1)

# Next filter out CREs with only 1 sub-CRE
CREs <- regionRanges(project)
numSubCRE <- table(findOverlaps(subCRERanges, CREs)@to)
keptCREs <- as.integer(names(numSubCRE[numSubCRE > 1]))
numSubCREFilter <- findOverlaps(subCRERanges, CREs[keptCREs])@from

# Combine the above two filters
subCREFilter <- intersect(adjacencyFilter, numSubCREFilter)
subCRESE <- subCRESE[subCREFilter, ]
subCRERanges <- rowRanges(subCRESE)

########################################
# Find sub-CREs around erythroid genes #
########################################

# Get pseudo-bulks along the erythroid lineage
lineageInd <- 2
lineageProbs <- groupInfo[, c("MyeloidProbs", "ErythroidProbs", "BLymphoidProbs", "TLymphoidProbs")]
lineageIdentities <- sapply(1:dim(lineageProbs)[1], 
                            function(x){which(lineageProbs[x,] == max(lineageProbs[x,]))})
lineageGroups <- which(lineageIdentities %in% lineageInd)

# Re-order pseudo-bulks by pseudo-time
lineagePseudoTime <- groupPseudoTime(project)[lineageGroups]
lineageGroups <- lineageGroups[order(lineagePseudoTime)]

# Get pseudo-bulked RNA ordered by pseudotime
RNAMtx <- groupRNA(project)[, lineageGroups]

# Get pseudo-bulked ATAC ordered by pseudotime
ATACMtx <- groupATAC(project)[, lineageGroups]

# Find genes that are strongly activated during erythroid differentiation
ptimeCor <- cor(t(RNAMtx), 1:length(lineagePseudoTime))
ptimeCor[is.na(ptimeCor)] <- 0
eryGenes <- rownames(RNAMtx)[ptimeCor > 0.5]

# Find sub-CREs within 50kb from these genes and calculate correlation with pseudo-time
TSS <- BuenRTools::hg38TSSRanges
TSSFilt <- TSS[TSS$gene_name %in% eryGenes]
erySubCREOv <- findOverlaps(resize(TSSFilt, 100000, fix = "center"), subCRERanges)
erySubCREMtx <- assay(subCRESE)[erySubCREOv@to, lineageGroups]
erySubCRECor <- cor(t(erySubCREMtx), 1:length(lineagePseudoTime))

############################################
# Cluster sub-CREs by timing of activation #
############################################

# Select sub-CREs that are up-regulated during erythropoiesis
# Cluster them by their dynamics
erySubCREFilter <- which(erySubCRECor > 0.5)
erySubCREMtxFilt <- erySubCREMtx[erySubCREFilter, ]
erySubCRERanges <- subCRERanges[erySubCREOv@to[erySubCREFilter]]
erySubCREMtxScaled <- t(apply(erySubCREMtxFilt, 1, function(x){(x - min(x))/(max(x) - min(x))}))

# Rank sub-CREs by their timing of activation. Timing is determined by area under rescaled activity curve (AUC)
# (Higher AUC means earlier activation and vice versa)
smoothRadius <- 10
auc <- sapply(
  1:length(erySubCREFilter),
  function(i){
    
    # Smoothe curve
    curve <- conv(erySubCREMtxFilt[i,], smoothRadius)/(smoothRadius * 2)
    
    # Remove edges that are affected by out-of-bound smoothing
    curve <- curve[smoothRadius:(length(curve) - smoothRadius)]
    
    # Re-scale the curve to 0-1, otherwise AUC won't be comparable
    curve <- (curve - min(curve)) / (max(curve) - min(curve))
    
    # Calculate AUC
    auc <- mean(curve)
    auc
  }
)

# Divide into 5 clusters
clusterLabels <- floor(rank(-auc)/(length(auc) + 1)/0.2) + 1
clusters <- sort(unique(clusterLabels))

# Prepare sub-CRR-by-pseudotime matrix for plotting
smoothRadius <- 5
plotMtx <- erySubCREMtxScaled[order(clusterLabels), ] # Re-order clusters
plotMtx <- t(apply(plotMtx, 1, function(x){conv(x, smoothRadius)/(2 * smoothRadius)})) # Smooth tracks
plotMtx <- plotMtx[,smoothRadius:(dim(plotMtx)[2] - smoothRadius)] 
colnames(plotMtx) <- NULL

# Generate color map for heatmap colors
dataColors <- colorRamp2(seq(quantile(plotMtx, 0.01), quantile(plotMtx, 0.99),length.out=9),
                         colors = jdb_palette("solar_extra"))

# Generate side color bar indicating cluster membership of each sub-CRE
library(RColorBrewer)
colorPalette = brewer.pal.info[brewer.pal.info$category == 'qual',]
colorVector = unique(unlist(mapply(brewer.pal, colorPalette$maxcolors, rownames(colorPalette))))
clusterColors <- colorVector[clusters + 1]
names(clusterColors) <- clusters
rightAnno <- rowAnnotation("cluster" = factor(clusterLabels[order(clusterLabels)]),
                           col = list("cluster" = clusterColors))

# Generate top color bar indicating pseudotime
ptimeColors <- colorRamp2(seq(0, dim(plotMtx)[1],length.out=9),
                          colors = jdb_palette("solar_rojos"))
topAnno <- HeatmapAnnotation("Pseudotime" = (1:(dim(plotMtx)[2])) + smoothRadius,
                             col = list("Pseudotime" = ptimeColors)) # Clip edge region that are affected by smoothing

# Plot Heatmap
pdf("../../data/BMMC/plots/subCREDynamics/subCREClustering.pdf", width = 10, height = 8)
ComplexHeatmap::Heatmap(plotMtx, 
                        col = dataColors,
                        cluster_rows = F,
                        cluster_columns = F,
                        right_annotation = rightAnno,
                        top_annotation = topAnno,
                        name = "sub-cCRE\nactivity")
dev.off()

##############################################
# Calculate motif enrichment in each cluster #
##############################################

# Add GC bias
subCRESEFilt <- subCRESE[, lineageGroups]
subCRESEFilt <- chromVAR::addGCBias(subCRESEFilt, genome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

# Impute NA values
bias <- rowRanges(subCRESEFilt)$bias
rowRanges(subCRESEFilt)$bias[is.na(rowRanges(subCRESEFilt)$bias)] <- mean(bias[!is.na(bias)])

# Get background subCREs or CREs with matched GC bias and average signal
assayNames(subCRESEFilt) <- "counts"
bgSites <- chromVAR::getBackgroundPeaks(subCRESEFilt, niterations = 100)

# Get subCRE-by-TF motif matching matrix
motifMtx <- pbmcapply::pbmcmapply(
  function(TF){
    matchVec <- rep(0, length(rowRanges(subCRESE)))
    motifCenter <- resize(motifPositions[[TF]], 1, fix = "center")
    matchVec[findOverlaps(motifCenter, rowRanges(subCRESE))@to] <- 1
    matchVec
  },
  names(motifPositions),
  mc.cores = 16
)

# Calculate motif enrichment in each cluster
erySubCREInds <- erySubCREOv@to[erySubCREFilter]
motifEnrichment <- data.table::rbindlist(
  lapply(
    clusters,
    function(cluster){
      enrichment <- motifRegionZTest(regionSet = erySubCREInds[clusterLabels == cluster], 
                                     bgRegions = bgSites, tfMat = motifMtx)
      enrichment$fdr <- p.adjust(enrichment$pval.z, method = "fdr")
      enrichment$cluster <- cluster
      enrichment
    }
  )
)
motifEnrichment <- as.data.frame(motifEnrichment)
motifEnrichment <- motifEnrichment[!is.na(motifEnrichment$z_test),]

# Visualize enrichment of top TFs using heatmap
plotTFs <- Reduce(union, lapply(
  clusters,
  function(cluster){
    motifEnrichment[motifEnrichment$cluster == cluster,][1:10,]$motifID
  })
)
plotData <- motifEnrichment[motifEnrichment$motifID %in% plotTFs,]
plotMtx <- as.matrix(Matrix::sparseMatrix(
  i = match(plotData$motifID, plotTFs),
  j = plotData$cluster,
  x = -log10(plotData$fdr)
))
rownames(plotMtx) <- plotTFs
colnames(plotMtx) <- paste("Cluster", clusters)
plotMtx <- pmin(plotMtx, 5)
heat.cols <- colorRamp2(seq(0,3,length.out=9),colors = jdb_palette("solar_extra"))
pdf("../../data/BMMC/plots/subCREDynamics/motifEnrichment.pdf", width = 8, height = 8)
Heatmap(matrix = as.matrix(plotMtx),
        col = heat.cols,
        cluster_columns = F, cluster_rows = F)
dev.off()

########################################################
# Further characterize sub-CRE arrangement within CREs #
########################################################

# Calculate distance of each sub-CRE to CRE centers
CREs <- regionRanges(project)
ov <- findOverlaps(erySubCRERanges, CREs)
CRECenters <- resize(CREs, 1, fix = "center")
erySubCRECenters <- resize(erySubCRERanges, 1, fix = "center")
CRECenterDist <- abs(start(erySubCRECenters[ov@from]) - start(CRECenters[ov@to]))

# Compare the above distances across clusters
plotData <- data.frame(
  dist = CRECenterDist,
  cluster = clusterLabels[ov@from]
)
comparisons <- lapply(2:5, function(x){c("1", as.character(x))}) 
pdf("../../data/BMMC/plots/subCREDynamics/CREDist.pdf", width = 5, height = 5)
ggpubr::ggviolin(data = plotData, x = "cluster", y = "dist", fill = "cluster",
                 xlab = "Cluster", ylab = "Distance to CRE center",add = "median_q1q3") + 
  stat_compare_means(comparisons = comparisons)
dev.off()

# Construct a CRE-by-cluster matrix. This records the number of sub-CREs from each cluster in every CRE
# For example if row 10 in this matrix is 0 1 0 0 2, it means the 10th CRE has 1 sub-CRE in cluster 2
# and 2 sub-CREs in cluster 5
eryCREInds <- ov@to
CREClusterMtx <- as.matrix(Matrix::sparseMatrix(
  i = match(eryCREInds, sort(unique(eryCREInds))), 
  j = clusterLabels[ov@from],
  x = 1,
  dims = c(length(unique(eryCREInds)), length(clusters))
))

# For each CRE, find sub-CREs in it and count how many clusters are they from
nClusters <- rowSums(CREClusterMtx > 0)

# Calculate percentage cases whether sub-CREs of different clusters share the same CRE
sum(CREClusterMtx[nClusters > 1,]) / sum(CREClusterMtx)

# Find specific examples for the above situation
sort(sort(unique(eryCREInds))[which(rowSums(CREClusterMtx > 0) > 1)])
CREInd <- 119174
subRanges <- subsetByOverlaps(erySubCRERanges, CREs[CREInd])
start(resize(subRanges, 1, fix = "center")) - start(CREs[CREInd])

# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

myPath <- .libPaths()
myPath <- c(myPath,'/packages')
.libPaths(myPath)

library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(hdf5r)
library(sctransform)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(patchwork)
source("../../code/utils.R")
source("../../code/getGroupData.R")

###################
# Read input data #
###################

# Load single cell RNA data
scRNA <- readRDS("../../data/mHSCAging10xMultiome/scRNA.rds")

# Load single cell ATAC data
scATACSeurat <- readRDS("../../data/mHSCAging10xMultiome/scATACSeurat.rds")

# Only keep cells in the HSC cluster 
HSCBarcodes <- colnames(scATACSeurat)[scATACSeurat$seurat_clusters %in% c(0, 6, 10)]
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
scATACFilt <- scATACSeurat[, colnames(scRNAFilt)]

# Normalize RNA data
scRNAFilt <- SCTransform(
  scRNAFilt, verbose = T, 
  variable.features.n = 5000, 
  return.only.var.genes = F
)

# Load barcodes for each pseudobulk
barcodeGroups <- read.table("../../data/mHSCAging10xMultiome/barcodeGrouping.txt", header = T, sep = "\t")

# Normalized scRNA data
scRNAMat <- scRNAFilt[["SCT"]]@scale.data

########################
# Save data to h5 file #
########################

# Write input scRNA data to a file
h5File <- H5File$new("../../data/mHSCAging10xMultiome/scRNA.h5", mode="w")
h5File[["RNA_matrix"]] <- as.matrix(scRNAFilt[["SCT"]]@data)
barcodes <- colnames(scRNAFilt[["SCT"]]@data)
h5File[["barcodes"]] <- barcodes
h5File[["genes"]] <- rownames(scRNAFilt[["SCT"]]@data)
h5File[["UMAP"]] <- scATACSeurat@reductions$umap@cell.embeddings[barcodes,]
h5File$close_all()

###################################
# Normalize pseudobulked RNA data #
###################################

# Get gene-by-pseudobulk raw count matrix
scRNACounts <- scRNA@assays$RNA@counts
pseudobulkRNA <- getGroupRNA(scRNACounts, barcodeGroups)

# Use DESeq to estimate size factor
metadata <- data.frame(age = stringr::str_split_fixed(colnames(pseudobulkRNA), "_", 2)[,1])
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = pseudobulkRNA,
  colData = metadata,
  design = ~age)
dds <- DESeq2::DESeq(dds)

# Normalize RNA counts
# See: https://rdrr.io/bioc/DESeq2/man/sizeFactors.html
# "The sizeFactors vector assigns to each column of the count matrix a value, the size factor, 
#  such that count values in the columns can be brought to a common scale by dividing by the corresponding size factor"
RNAMatNormed <- t(t(pseudobulkRNA) / dds$sizeFactor)

# Save normalized data
saveRDS(RNAMatNormed, "../../data/mHSCAging10xMultiome/pseudobulkedRNANormed.rds")

######################################
# Score and cluster spectra programs #
######################################

# Run runSpectra.py to get spectra programs and then load them in
spectraPrograms <- read.table("../../data/mHSCAging10xMultiome/spectra/spectra_markers_lam_10.tsv", 
                              sep = "\t", header = T)
programNames <- colnames(spectraPrograms)

# For each gene, find background genes with matched overall expression
# Even though we will scale the expression of each gene, this is still helpful because
# expression level is associated with biological function
geneExpLevel <- rowMeans(RNAMatNormed)
expKNN <- FNN::get.knn(geneExpLevel, k = 100)$nn.index

# Generate cell-by-program matrix of program scores
programMat <- pbmcapply::pbmcmapply(
  function(programInd){
    
    # Calculate spectra program scores as foreground
    programGenes <- spectraPrograms[, programInd]
    programGenes <- intersect(programGenes, rownames(RNAMatNormed))
    programGenesInds <- match(programGenes, rownames(RNAMatNormed))
    programGeneMat <- RNAMatNormed[programGenes, ]
    programGeneMat <- t(scale(t(programGeneMat)))
    programScores <- colMeans(programGeneMat)
    
    # Generate background programs consisting of expression-matched genes
    # Use them to calculate background scores
    bgScores <- sapply(
      1:dim(expKNN)[2],
      function(i){
        bgGeneInds <- expKNN[programGenesInds, i]
        bgGeneMat <- RNAMatNormed[bgGeneInds, ]
        bgGeneMat <- t(scale(t(bgGeneMat)))
        colMeans(bgGeneMat)
      }
    )
    
    # Get z-scores
    zscores <- (programScores - rowMeans(bgScores)) / rowSds(bgScores)
  },
  1:ncol(spectraPrograms),
  mc.cores = 8
)
colnames(programMat) <- stringr::str_split_fixed(colnames(spectraPrograms), "\\.", 5)[, 5]

# Save matrix to table
write.table(
    programMat,
    "../../data/mHSCAging10xMultiome/pbulk_by_spectra_program_mat.tsv", sep = "\t",
    row.names = T, col.names = T, quote = F)

# Cluster pseudobulks
youngPbulkClusters <- kmeans(cor(t(programMat[metadata$age == "Young", ])), 
                             centers = 3, iter.max = 100, nstart = 10)$cluster
oldPbulkClusters <- kmeans(cor(t(programMat[metadata$age == "Old", ])), 
                           centers = 4, iter.max = 100, nstart = 10)$cluster
pbulkClusters <- rep("", dim(programMat)[1])
pbulkClusters[metadata$age == "Young"] <- youngPbulkClusters
pbulkClusters[metadata$age == "Old"] <- oldPbulkClusters
pbulkClusters <- paste(metadata$age, pbulkClusters, sep = "_")
write.table(pbulkClusters, "../../data/mHSCAging10xMultiome/pbulkClusters.txt", 
            quote = F, row.names = F, col.names = F)

# Cluster programs
set.seed(123)
programClusterLabels <- kmeans(cor(programMat), centers = 4, iter.max = 100, nstart = 10)$cluster
programClusters <- sort(unique(programClusterLabels))
programClusterOrder <- order(-sapply(
    programClusters,
    function(cluster){
        clusterScore <- rowMeans(programMat[, programClusterLabels == cluster])
        ratio <- mean(clusterScore[metadata$age == "Young"]) - mean(clusterScore[metadata$age == "Old"])
    }
))

# Group pseudo-bulk clusters into sub-populations
subpopAnno <- list(
  "Old_1" = "Old Mk-biased", 
  "Old_2" = "Old intermediate",
  "Old_3" = "Old Mk-biased",
  "Old_4" = "Old multi-lineage",
  "Young_1" = "Young multi-lineage",
  "Young_2" = "Young multi-lineage",
  "Young_3" = "Young Mk-biased")
subpopLabels <- sapply(pbulkClusters, function(x){subpopAnno[[x]]})

# Order subpopulations
subpopOrder <- c("Young multi-lineage", "Young Mk-biased", 
                 "Old intermediate", "Old multi-lineage", "Old Mk-biased")

# Plot program-by-pseudobulk clustermap
colors <- circlize::colorRamp2(seq(quantile(programMat, 0.05), quantile(programMat, 0.95),length.out=9),
                               colors = BuenColors::jdb_palette("solar_rojos"))
pdf("../../data/mHSCAging10xMultiome/plots/pseudobulk_cluster_heatmap.pdf",
    width = 20, height = 30)
Heatmap(t(programMat), 
        col = colors,
        column_split = subpopLabels, 
        row_split = programClusterLabels)
dev.off()

# Plot a filtered version by keeping top variable programs
pdf("../../data/mHSCAging10xMultiome/plots/pseudobulk_cluster_heatmap_filt.pdf",
    width = 20, height = 30)
filter <- rank(-colSds(programMat)) < 30
options(repr.plot.width = 15, repr.plot.height = 10)
p <- Heatmap(programMat[, filter], 
        col = colors,
        row_split = factor(subpopLabels, levels = subpopOrder), 
        cluster_row_slices = FALSE,
        column_split = factor(programClusterLabels[filter], levels = programClusters[programClusterOrder]),
        cluster_column_slices = FALSE)
dev.off()

pdf("../../data/mHSCAging10xMultiome/plots/pseudobulk_cluster_correlation_map.pdf")
colors <- circlize::colorRamp2(seq(-1, 1,length.out=9),
                               colors = BuenColors::jdb_palette("solar_extra"))
corMat <- cor(t(programMat))
dimnames(corMat) <- NULL
Heatmap(corMat, 
        col = colors,
        row_split = subpopLabels,
        column_split = subpopLabels)
dev.off()

###############################################
# Plot pathway z-scores in each subpopulation #
###############################################

plotPathways <- c(
  "Rodriguez_Fraticelli_et_al_mkBiased", "all_unfolded.protein.response", 
  "all_MHC.I.presentation", "all_G2M.transition", "all_G2M.transition",
  "all_VAL.LEU.ILE_metabolism", "all_oxidative.phosphorylation",
  "all_hypoxia.response")
for(pathway in plotPathways){
  plotData <- data.frame(
    subpop = subpopLabels,
    expression = programMat[, pathway]
  )
  pdf(paste0("../../data/mHSCAging10xMultiome/plots/subpop_", pathway, ".pdf"), width = 8, height = 4)
  print(ggpubr::ggboxplot(
    plotData, x = "subpop", y = "expression", fill = "subpop", width = 0.3,
    ylab = "Program z-score", xlab = "", palette = "npg", 
    add.params = list(fill = "white"), add = "jitter") +
      ylim(c(min(plotData$expression), max(plotData$expression) * 2)) +
      ggpubr::stat_compare_means(
        comparisons = list(c("Old multi-lineage", "Young multi-lineage"), 
                           c("Old Mk-biased", "Young Mk-biased"),
                           c("Old multi-lineage", "Old Mk-biased"),
                           c("Young multi-lineage", "Young Mk-biased")),
        label = "p.signif"))
  dev.off()
}

##########################################
# Visualize pseudo-bulk clusters on UMAP #
##########################################

# Load pseudobulk center cel barcodes
pseudobulkCenters <- read.table("../../data/mHSCAging10xMultiome/pseudobulkCenters.txt", sep = "\t", header = T)
pseudobulkAges <- sapply(pseudobulkCenters$barcode, function(x){strsplit(x, "-")[[1]][2]})

# Get UMAP coordinates of cells
plotData <- data.frame(
  UMAP1 = scATACSeurat@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATACSeurat@reductions$umap@cell.embeddings[, 2]
)
rownames(plotData) <- colnames(scATACSeurat)
plotData <- plotData[(plotData$UMAP1 > -5) & (plotData$UMAP1 < 0) &
                       (plotData$UMAP2 > -1) & (plotData$UMAP2 < 4), ]

# Get the UMAP coordinate and subpopulation labels of young and old pseudo-bulks
plotDataPbulk <- plotData[pseudobulkCenters$barcode, ]
plotDataPbulk$subpopulation <- subpopLabels

# Label young and old pseudo-bulks separately
oldLabels <- sapply(pseudobulkCenters$group[pseudobulkAges == "Old"],
                    function(s){strsplit(s, "_")[[1]][2]})
youngLabels <- sapply(pseudobulkCenters$group[pseudobulkAges == "Young"],
                      function(s){strsplit(s, "_")[[1]][2]})
png("../../data/mHSCAging10xMultiome/plots/pseudobulkClusters.png",
    width = 1000, height = 1000)
ggplot(plotData) +
  ggrastr::rasterise(geom_point(aes(x = UMAP1, y = UMAP2), color = "grey", size = 0.5, alpha = 0.3)) +
  geom_point(data = plotDataPbulk, aes(x = UMAP1, y = UMAP2, color = subpopulation), size = 5) +
  geom_text_repel(data = plotDataPbulk[pseudobulkAges == "Young", ], 
                  aes(x = UMAP1, y = UMAP2), label = youngLabels, color = "#3361A5",
                  max.overlaps = 100, size = 12) +
  geom_text_repel(data = plotDataPbulk[pseudobulkAges == "Old", ],
                  aes(x = UMAP1, y = UMAP2), label = oldLabels, color = "#A31D1D",
                  max.overlaps = 100, size = 12) +
  theme_classic()  
dev.off()

###########################################
# Score each gene program in single cells #
###########################################

# Get cell embedding in LSI and UMAP space
cellEmbedding <- scATACSeurat@reductions$lsi@cell.embeddings[colnames(scRNAFilt), 2:20]
cellUMAP <- scATACSeurat@reductions$umap@cell.embeddings[colnames(scRNAFilt), ]

# Get KNN for smoothing.
KNN <- FNN::get.knn(cellEmbedding, k = 50)$nn.index

# Select a program to visualize
programInd <- which(stringr::str_detect(programNames, "Fraticelli_et_al_mkBiased"))
programGenes <- spectraPrograms[, programInd]
programName <- colnames(programMat)[programInd]
programScores <- colMeans(t(scale(t(scRNAMat[programGenes, ]))))
#write.table(programGenes, "test.txt", sep = "\t", quote = F, row.names = F, col.names = F)

# Smooth the RNA of marker genes
programScores <- pbmcapply::pbmcmapply(
  function(cellInd){
    mean(programScores[KNN[cellInd, ]])
  },
  1:dim(scRNAFilt)[2],
  mc.cores = 12
)

# Visualize data on HSC UMAP
pltData <- data.frame(
  UMAP1 = cellUMAP[, 1],
  UMAP2 = cellUMAP[, 2],
  marker = programScores
)
pltData <- pltData[(pltData$UMAP1 > -5) & (pltData$UMAP1 < 0) &
                     (pltData$UMAP2 > -1) & (pltData$UMAP2 < 4), ]
p <- ggplot(pltData[sample(1:dim(pltData)[1]),]) +
  geom_point(aes(x =  UMAP1, y = UMAP2, color = marker), size = 1.5, alpha = 1) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_rojos")) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic() +
  ggtitle(programName)
png(paste0("../../data/mHSCAging10xMultiome/plots/UMAP_", programName, ".png"),
    width = 1000, height = 1000)
p
dev.off()

##############################################################
# Visualize previously published aging signature on the UMAP #
##############################################################

# Get cell embedding in LSI and UMAP space
cellEmbedding <- scATACSeurat@reductions$lsi@cell.embeddings[colnames(scRNAFilt), 2:20]
cellUMAP <- scATACSeurat@reductions$umap@cell.embeddings[colnames(scRNAFilt), ]

# Get KNN for smoothing.
KNN <- FNN::get.knn(cellEmbedding, k = 50)$nn.index

# Get aging signature score using the aging genes published in this paper doi: 10.1038/ncomms11075 (2016).
agingGenes <- c(
  "Sdpr", "Nupr1", "Btg2", "Ddx3y", "Mt1", "Eif2s3y", "Sult1a1", 
  "Egr1", "Selp", "Aldh1a1", "Clca1", "Tmem181b-ps", "H2-Ab1", 
  "Gstm1", "Gda", "H2-Aa", "Clec1a", "Plek", "Plscr2", "Uty", 
  "Tmem181c-ps", "Nr4a1", "Id2", "Iigp1", "Ier2", "Dhrs3", 
  "Klhl4", "Cd74", "Oxr1", "Ube2h", "Nfkbiz", "Itgb3", "Casp12", 
  "5430417L22Rik", "Ier3", "Enpp5", "Sat1", "Egr3", "Itga6", 
  "Clu", "Prtn3", "Exoc6b", "Fhl1", "H2-Eb1", "Prune2", 
  "Gata2", "M6pr", "Ccnt1", "Ndrg1", "Tgm2", "Vmp1", 
  "Rnf167", "Tsc22d1", "Gstm2", "Gpr183", "Cers2", 
  "Muc13", "Il2rg", "Art4", "Gstm3"
)
agingGenes <- intersect(agingGenes, rownames(scRNAMat))
agingScores <- colMeans(t(scale(t(scRNAMat[agingGenes, ]))))

# Smooth the RNA of marker genes
agingScores <- pbmcapply::pbmcmapply(
  function(cellInd){
    mean(agingScores[KNN[cellInd, ]])
  },
  1:dim(scRNAFilt)[2],
  mc.cores = 12
)

# Visualize data on HSC UMAP
pltData <- data.frame(
  UMAP1 = cellUMAP[, 1],
  UMAP2 = cellUMAP[, 2],
  score = agingScores
)
pltData <- pltData[(pltData$UMAP1 > -5) & (pltData$UMAP1 < 0) &
                     (pltData$UMAP2 > -1) & (pltData$UMAP2 < 4), ]
p <- ggplot(pltData[sample(1:dim(pltData)[1]),]) +
  geom_point(aes(x =  UMAP1, y = UMAP2, color = score), size = 1.5, alpha = 1) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_rojos")) +
  geom_rect(aes(xmin = -4.5, xmax = -1, ymin = -0.5, ymax = 2), 
            fill = "blue", alpha = 0, color = "black") +
  theme_classic() +
  ggtitle("Grover et. al. Aging signature")
png(paste0("../../data/mHSCAging10xMultiome/plots/UMAP_Grover_et_al_aging_signature.png"),
    width = 1000, height = 1000)
p
dev.off()

###########################################
# Color single cells by HSC state on UMAP #
###########################################

# Get cell embedding in LSI and UMAP space
cellEmbedding <- scATACSeurat@reductions$lsi@cell.embeddings[, 2:20]

# Only keep HSCs
seuratClusters <- scATACSeurat$seurat_clusters
cellType <- scATACSeurat$group
HSCfilter <- (cellType == "HSC") & (seuratClusters %in% c(0,6,10))
cellEmbedding <- cellEmbedding[HSCfilter, ]

# For each cell, find its nearest pseudobulk center cell with annotated subpopulation
KNN <- FNN::get.knnx(cellEmbedding[pseudobulkCenters$barcode, ], cellEmbedding, k = 1)$nn.index
scSubpopLabels <- unname(subpopLabels[KNN])

# Visualize results
colors <- c("#efc86e", "#5c66a8", "#6f9969", "#808fe1", "#97c684")
names(colors) <- unique(subpopLabels)
plotData <- data.frame(
  UMAP1 = scATACSeurat@reductions$umap@cell.embeddings[HSCfilter, 1],
  UMAP2 = scATACSeurat@reductions$umap@cell.embeddings[HSCfilter, 2],
  subpop = scSubpopLabels
)
plotData <- plotData[(plotData$UMAP1 > -5) & (plotData$UMAP1 < 0) &
                       (plotData$UMAP2 > -1) & (plotData$UMAP2 < 4), ]
png("../../data/mHSCAging10xMultiome/plots/HSC_subpop_single_cell_UMAP.png",
    width = 800, height = 800)
ggplot(plotData) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = subpop), 
             size = 1, alpha = 0.5) +
  scale_color_manual(values = colors) +
  theme_classic()  
dev.off()

#############################################################
# Plot genomic tracks of marker genes in each subpopulation #
#############################################################

# Load ATAC fragments
frags <- data.table::fread(
  "../../data/mHSCAging10xMultiome/all.frags.filt.tsv.gz",
  nThread = 16, showProgress = T, sep = "\t"
)
colnames(frags) <- c("chr", "start", "end", "barcode")

# Only keep HSC frags
HSCBarcodes <- colnames(scATACSeurat)[HSCfilter]
frags <- frags[frags$barcode %in% HSCBarcodes, ]

# Map each cell to a subpopulation
names(scSubpopLabels) <- HCSBarcodes
subpops <- unique(scSubpopLabels)

# Calculate coverage in a region-of-interest in each subpopulation
markerGene <- "Socs2"
TSSs <- FigR::mm10TSSRanges
names(TSSs) <- as.character(TSSs$gene_name)
locus <- resize(TSSs[markerGene], 1e4, fix = "center")
locusBins <- tile(locus, width = 100)[[1]]
tracks <- lapply(
  subpops,
  function(subpop){
    subpopFrags <- frags[scSubpopLabels[frags$barcode] %in% subpop,]
    subpopFragRanges <- GRanges(seqnames = subpopFrags$chr, ranges = IRanges(start = subpopFrags$start, end = subpopFrags$end))
    nMillionFrags <- dim(subpopFrags)[1] / 1e6
    coverage <- as.integer(table(factor(findOverlaps(locusBins, subpopFragRanges)@from, levels = 1:length(locusBins))))
    coverageRPM <- coverage / nMillionFrags # Normalize to get read-per-million
  }
)
names(tracks) <- subpops

# Visualize results
pdf(paste0("../../data/mHSCAging10xMultiome/plots/genomic_tracks/", markerGene, ".pdf"), width = 3, height = 7)
ymax <- max(sapply(tracks, max))
plots <- lapply(
  subpops,
  function(subpop){
    plotData <- data.frame(
      position = start(resize(locusBins, 1, fix = "center")),
      coverage = tracks[[subpop]]
    )
    ggplot(plotData) +
      geom_ribbon(aes(x = position, ymin=0, ymax=coverage)) +
      xlab("") + ylim(0, ymax) + ylab(subpop) +
      theme_classic()
  }
)
Reduce("/", plots)
dev.off()

Reduce("/", plots)

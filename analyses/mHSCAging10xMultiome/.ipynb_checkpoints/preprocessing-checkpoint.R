# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(ggplot2)
source("../../code/utils.R")

#########################
# Merge fragments files #
#########################

metadata <- data.frame(
  fragsFiles = c("HHVCLBGXJ_Multiome_GEX_ATAC_Output_OldHSC_atac_fragments.tsv",
                 "HHVCLBGXJ_Multiome_GEX_ATAC_Output_OldLIN_atac_fragments.tsv",
                 "HHVCLBGXJ_Multiome_GEX_ATAC_Output_YoungHSC_atac_fragments.tsv",
                 "HHVCLBGXJ_Multiome_GEX_ATAC_Output_YoungLIN_atac_fragments.tsv"),
  fragsFolder = rep("../../data/mHSCAging10xMultiome/frags/", 4),
  age = c("Old", "Old", "Young", "Young"),
  group = c("HSC", "LinNeg", "HSC", "LinNeg"),
  replicate = c("rep1", "rep1", "rep1", "rep1"),
  dataset = rep("Multiome", 4)
)

if(!file.exists("../../data/mHSCAging10xMultiome/all.frags.tsv.gz")){
  
  frags <- pbmcapply::pbmclapply(
    1:dim(metadata)[1],
    function(i){
      fgs <- data.table::fread(paste0(metadata$fragsFolder[i], metadata$fragsFiles[i]))
      fgs$V4 <- paste0(metadata$dataset[i], "-", metadata$age[i], "-", metadata$group[i], "-", 
                       metadata$replicate[i], "-", fgs$V4)
      fgs[, 1:4]
    }
  ) 
  
  # Combine fragments from different samples and save to a file
  fragsCombined <- data.table::rbindlist(frags)
  data.table::fwrite(fragsCombined, "../../data/mHSCAging10xMultiome/all.frags.tsv.gz", 
                     compress = "gzip", sep = "\t")
  
}

##################
# Get ATAC peaks #
##################

# Get peak ranges and resize to 1kb
peakBed <- read.table("../../data/mHSCAging10xMultiome/peaks.bed")
peaks <- GRanges(seqnames = peakBed$V1,
                 ranges = IRanges(start = peakBed$V2,
                                  end = peakBed$V3))
peaks <- IRanges::resize(peaks, width = 1000, fix = "center")
saveRDS(peaks, "../../data/mHSCAging10xMultiome/regionRanges.rds")
peakBed <- data.frame(chr = as.character(seqnames(peaks)),
                      start = start(peaks),
                      end = end(peaks))
write.table(peakBed, "../../data/mHSCAging10xMultiome/peaks.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")

############################################
# Generating count matrix and filter cells #
############################################

# Here go to run data/mHSCAging10xMultiome/pipeline.sh. It will generate a single cell SummarizedExperiment object and filter cells based on QC

# Then we run ArchR/doubletScoring.R to compute doublet scores

######################################################
# Filter cells based on depth and doublet enrichment #
######################################################

# Load scATAC
if(!file.exists("../../data/mHSCAging10xMultiome/scATACSE.rds")){
  
  scATACSE <- readRDS("../../data/mHSCAging10xMultiome/scATAC.rds")
  
  # Filter by depth
  scATACSE <- scATACSE[, scATACSE@colData$depth > 300]
  
  # Filter by doublet scores
  doubletEnrichment <- readRDS("../../data/mHSCAging10xMultiome/doubletEnrichment.rds")
  doubScores <- rep(0, dim(scATACSE)[2])
  names(doubScores) <- colnames(scATACSE)
  ovCells <- intersect(names(doubScores), names(doubletEnrichment))
  doubScores[ovCells] <- doubletEnrichment[ovCells]
  scATACSE <- scATACSE[, doubScores < quantile(doubScores, 0.95)]
  
  # Save to file
  saveRDS(scATACSE, "../../data/mHSCAging10xMultiome/scATACSE.rds")
  
}else{
  scATACSE <- readRDS("../../data/mHSCAging10xMultiome/scATACSE.rds")
}

#####################
# preprocess scATAC #
#####################

peaks <- rowRanges(scATACSE)

# Get metadata of single cells
metadata <- as.data.frame(t(sapply(colnames(scATACSE), function(x){strsplit(x, "-")[[1]][1:4]})))
colnames(metadata) <- c("dataset", "age", "group", "replicate")

# Create Seurat object
ATACCounts <- assay(scATACSE)
rownames(ATACCounts) <- as.character(peaks)
scATACSeurat <- CreateSeuratObject(counts = ATACCounts, assay = "ATAC")
scATACSeurat$age <- metadata$age
scATACSeurat$group <- metadata$group
scATACSeurat$dataset <- metadata$dataset
scATACSeurat$replicate <- metadata$replicate
scATACSeurat$depth <- scATACSE$depth
scATACSeurat$FRIP <- scATACSE$FRIP

# Remove peaks with zero counts
scATACSeurat <- scATACSeurat[rowSums(ATACCounts) > 0,]

# LSI embedding
DefaultAssay(scATACSeurat) <- "ATAC"
scATACSeurat <- Signac::RunTFIDF(scATACSeurat)
scATACSeurat <- Signac::FindTopFeatures(scATACSeurat, min.cutoff = 'q80')
scATACSeurat <- Signac::RunSVD(scATACSeurat)
LSIEmbedding <- scATACSeurat@reductions$lsi@cell.embeddings
saveRDS(LSIEmbedding, "../../data/mHSCAging10xMultiome/LSIEmbedding.rds")

# Examine LSI-depth correlation
cor(LSIEmbedding, scATACSE@colData$depth)

# Clustering
scATACSeurat <- FindNeighbors(scATACSeurat, dims = 2:20, reduction = "lsi") # Don't use the first PC because it's highly correlated with sequencing depth
scATACSeurat <- FindClusters(scATACSeurat, resolution = 0.5)

# UMAP embedding
set.seed(42)
scATACSeurat <- RunUMAP(scATACSeurat, 
                        reduction = "lsi", 
                        dims = 2:20,
                        n.neighbors = 20,
                        min.dist = 0.01,
                        spread = 0.5,
                        repulsion.strength = 1)

# Visualize embedding
scATACSeurat@active.ident <- factor(scATACSeurat$seurat_clusters)
DimPlot(scATACSeurat, reduction = "umap")

plot <- DimPlot(scATACSeurat, reduction = "umap") 
system("mkdir ../../data/mHSCAging10xMultiome/plots/")
pdf("../../data/mHSCAging10xMultiome/plots/UMAPClusters.pdf", width = 8, height = 7)
LabelClusters(plot = plot, id = "ident")
dev.off()

#################
# Visualization #
#################

# Calculate gene scores
geneScores <- getGeneScoresFromPeaks(scATACSE[, colnames(scATACSeurat)], 
                                     genome = "mm10", 
                                     TSSwindow = 10000, 
                                     getWeightsOnly = FALSE)

# Normalize gene scores
geneScores <- centerCounts(geneScores)

# Calculate gene score of cluster markers
seuratClusters <- scATACSeurat$seurat_clusters
markerList <- read.csv("../../data/mHSCAging10xMultiome/hemeMarkers.txt")
markers <- Reduce(c, lapply(
  colnames(markerList), 
  function(cellType){
    cellTypeMarkers <- markerList[, cellType]
    cellTypeMarkers <- cellTypeMarkers[cellTypeMarkers != ""]
    names(cellTypeMarkers) <- rep(cellType, length(cellTypeMarkers))
    cellTypeMarkers
  }))
markers <- markers[markers %in% rownames(geneScores)]
clusterMarkerSignal <- sapply(
  as.character(sort(unique(seuratClusters))),
  function(cluster){
    rowMeans(geneScores[markers, seuratClusters == cluster])
  }
)

set.seed(42)
clusterMarkerSignal <- (clusterMarkerSignal - rowMeans(clusterMarkerSignal)) / rowSds(clusterMarkerSignal)
pdf("../../data/mHSCAging10xMultiome/plots/clusterMarkers.pdf", width = 15, height = 20)
library(ComplexHeatmap)
leftAnno <- rowAnnotation(cellType=factor(names(markers), levels = unique(names(markers))),
                          border=TRUE,show_annotation_name=FALSE)
ComplexHeatmap::Heatmap(
  clusterMarkerSignal,
  left_annotation = leftAnno,
  cluster_rows = F
)
dev.off()

# Annotate cell types
cellTypeMapping <-  c("HSC", "Neu-GMP", "GMP", "GMP", "Mono-GMP",
                      "Baso/Mast", "HSC", "MEP", "Plasma cells",
                      "T cells", "HSC", "MEP", "B cells", "Granulocytes")
names(cellTypeMapping) <- as.character(0:13)
scATACSeurat$cellType <- cellTypeMapping[as.character(scATACSeurat$seurat_clusters)]

# Visualize data on UMAP
pltData <- data.frame(
  UMAP1 = scATACSeurat@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scATACSeurat@reductions$umap@cell.embeddings[, 2],
  cluster = seuratClusters %in% c(0,6,10),
  cellType = scATACSeurat$cellType,
  group = scATACSeurat$group,
  depth = log10(scATACSE$depth),
  age = scATACSeurat$age
)

clusterColors <- as.character(MetBrewer::met.brewer("Derain", 14))
pdf("../../data/mHSCAging10xMultiome/plots/UMAPSamples.pdf", width = 8, height = 7)
ggplot(pltData[sample(1:dim(pltData)[1]),]) +
  geom_point(aes(x =  UMAP1, y = UMAP2, color = group), size = 0.1, alpha = 0.5) +
  scale_color_manual(values = c("#E64B35FF", "#91D1C2FF")) +
  #scale_color_manual(values = colors) +
  theme_classic()
dev.off()

# Save results
saveRDS(scATACSeurat, "../../data/mHSCAging10xMultiome/scATACSeurat.rds")

####################
# preprocess scRNA #
####################

# Load single cell RNA count matrix
youngHSCRNA <- Seurat::Read10X("../../data/mHSCAging10xMultiome/rawFeatureMatrix/YoungHSC_raw_feature_bc_matrix")[["Gene Expression"]]
oldHSCRNA <- Seurat::Read10X("../../data/mHSCAging10xMultiome/rawFeatureMatrix/OldHSC_raw_feature_bc_matrix")[["Gene Expression"]]

# Add age tag to barcodes
colnames(youngHSCRNA) <- paste0("Multiome-Young-HSC-rep1-", colnames(youngHSCRNA))
colnames(oldHSCRNA) <- paste0("Multiome-Old-HSC-rep1-", colnames(oldHSCRNA))

# Generate seurat object
scRNACounts <- cbind(youngHSCRNA, oldHSCRNA)
scRNA <- CreateSeuratObject(counts = scRNACounts)
scRNA$age <- c(rep("Young", dim(youngHSCRNA)[2]), rep("Old", dim(oldHSCRNA)[2]))

# Only keep cells overlapping with our filtered scATAC dataset
scRNA <- scRNA[, colnames(scATACSeurat)]

# Normalize the data by total counts per cell (followed by log transform)
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# Find top variable features
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 1000)

# Scaling the data
all_genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all_genes)

# PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA, ndims = 50)

# Clustering
scRNA <- FindNeighbors(scRNA, dims = 2:10)
scRNA <- FindClusters(scRNA, resolution = 0.1)

# UMAP dimensionality reduction
set.seed(42)
scRNA <- RunUMAP(scRNA, dims = 2:10, n.neighbors = 50)
plotData <- data.frame(
  UMAP1 = scRNA@reductions$umap@cell.embeddings[, 1],
  UMAP2 = scRNA@reductions$umap@cell.embeddings[, 2],
  Cluster = scRNA$seurat_clusters,
  Age = scRNA$age 
)
ggplot(plotData) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = Age)) +
  theme_classic()

FeaturePlot(scRNA, c("Hlf", "Gata1", "Fos", "Clu", "Selp", "Nupr1"))

# Save results
saveRDS(scRNA, "../../data/mHSCAging10xMultiome/scRNA.rds")

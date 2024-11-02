# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(ggrepel)

###################
# Load input data #
###################

# Load ATAC peak ranges
regions <- readRDS("../../data/mHSCAging10xMultiome/regionRanges.rds")

# Load peak-by-pseudobulk count matrix
pbulkCounts <- readRDS("../../data/mHSCAging10xMultiome/pseudobulkATAC.rds")

# Load pseudobulk clustering results
pbulkClusters <- read.table("../../data/mHSCAging10xMultiome/pbulkClusters.txt", sep = "\t")$V1

# Group pseudo-bulk clusters into sub-populations
subpopAnno <- list(
  "Old_1" = "Old Mk-biased", 
  "Old_2" = "Old intermediate",
  "Old_3" = "Old Mk-biased",
  "Old_4" = "Old multi-lineage",
  "Young_1" = "Young multi-lineage",
  "Young_2" = "Young multi-lineage",
  "Young_3" = "Young Mk-biased")
subpopLabels <- unname(sapply(pbulkClusters, function(x){subpopAnno[[x]]}))

# Load TF motifs
cisBPMotifs <- readRDS("../../data/shared/cisBP_mouse_pwms_2021.rds")

# Load differential RNA results
diffRNA <- read.table("../../data/mHSCAging10xMultiome/diffRNA.tsv", sep = "\t")

##########################
# Calculate motif scores #
##########################

# Convert to RangedSummarizedExperiment object
pbulkSE <- SummarizedExperiment(list(counts = pbulkCounts), rowRanges = regions)

# Remove peaks with 0 counts (This could happen because we are using a consensus peak set that is not called on this dataset)
pbulkSE <- pbulkSE[rowSums(assay(pbulkSE)) > 0,]

# Add GC bias
pbulkSE <- chromVAR::addGCBias(pbulkSE, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)

# Remove NA values
bias <- rowRanges(pbulkSE)$bias
rowRanges(pbulkSE)$bias[is.na(rowRanges(pbulkSE)$bias)] <- mean(bias[!is.na(bias)])

# Get background substructures or CREs with matched GC bias and average signal
bgSites <- chromVAR::getBackgroundPeaks(pbulkSE, niterations = 100) 

# Get pseudobulk-by-TF motif score matrix
motifScores <- pbmcapply::pbmcmapply(
  function(TF){
    
    # Find motif matches
    motifMatches <- motifmatchr::matchMotifs(cisBPMotifs[TF], 
                                             pbulkSE, 
                                             genome = "mm10")
    
    # Compute motif scores (deviation score from background)
    deviation <- chromVAR::computeDeviations(object = pbulkSE, 
                                             annotations = motifMatches, 
                                             background_peaks = bgSites)
    deviation <- assay(deviation)
  },
  names(cisBPMotifs),
  mc.cores = 8
)
rownames(motifScores) <- colnames(pbulkSE)
motifScores[is.na(motifScores)] <- 0
saveRDS(motifScores, "../../data/mHSCAging10xMultiome/motifScores.rds")

#####################
# Visualize results #
#####################

selectedTFs <- c("Stat3", "Jun", "Egr1", "Pbx3", "Smad1", "Mycn", "Spi1", "Nfe2", "Fli1", 
                 "Klf6", "Xbp1", "Bcl11a", "Nfe2l2", "Atf4", "Atf6b", "Tfec", "Ets1", "Clock")
RNAColors <- colorRamp2(seq(-2,2,length.out=9), colors = jdb_palette("solar_extra"))
topAnno <- HeatmapAnnotation(
  "RNAlog2FC" = diffRNA[selectedTFs, ]$log2FoldChange,
  col = list(RNAlog2FC=RNAColors))
pdf("../../data/mHSCAging10xMultiome/plots/motifScores.pdf", width = 10, height = 15)
Heatmap(
  motifScores[, selectedTFs], 
  row_split = subpopLabels,
  top_annotation = topAnno)
dev.off()

####################################
# Differential motif score testing #
####################################

age <- stringr::str_split_fixed(rownames(motifScores), "_", 2)[, 1]
diffMean <- colMeans(motifScores[age == "Old", ]) - colMeans(motifScores[age == "Young", ]) 
pvals <- twoSampTTest(t(motifScores[age == "Old", ]), t(motifScores[age == "Young", ]))
diffResults <- data.frame(
  TF = colnames(motifScores),
  diffMean = diffMean, 
  pvals = pvals
)

# Also save as tsv file
write.table(diffResults, "../../data/mHSCAging10xMultiome/diffMotifScores.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

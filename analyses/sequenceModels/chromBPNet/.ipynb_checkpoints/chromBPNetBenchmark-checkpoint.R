# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../../code/utils.R")
library(rhdf5)
library(matrixStats)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(patchwork)

###########################
# Load chromBPNet results #
###########################

# Load CRE region ranges
regions <- readRDS("../../../data/HepG2/regionRanges.rds")

# Load chromBPNet contribution scores
modes <- c("profile", "count")
h5Paths <- list(
  "profile" = "../../../data/sequenceModels/chromBPNet/K562/fold_0_nobias_contribs_liftoverpeak.profile_scores.h5",
  "count" = "../../../data/sequenceModels/chromBPNet/K562/fold_0_nobias_contribs_liftoverpeak.counts_scores.h5"
)
scores <- list()
for(mode in modes){
  h5f <- H5Fopen(h5Paths[[mode]])
  scores[[mode]] <- h5f$projected_shap$seq
  h5closeAll()
}

# Keep the center region where we have footprinting results
regionRadius <- 500
paddedRadius <- dim(scores$profile)[1] / 2
scores <- lapply(scores, function(x){x[(paddedRadius - regionRadius + 1):((paddedRadius + regionRadius)), , ]})

# Get region-by-position score matrix
scoreMtx <- lapply(
  scores,
  function(x){
    t(sapply(
      1:dim(x)[3],
      function(i){rowMaxs(abs(x[, ,i]))}
    ))
  }
)

###################
# Load input data #
###################

projectName <- "K562"

# Load TF ChIP ranges
TFChIPRanges <- readRDS(paste0("../../../data/shared/unibind/", projectName, "ChIPRanges.rds"))

# Load chain file required for lift-over
pathToChain <- "../../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Load TF cluster membership
TFClustering <- read.csv("../../../data/TFBSPrediction/clusterLabelsAllTFs.txt", sep = "\t")

# Load CRE regions
regions <- readRDS("../../../data/K562/regionRanges.rds")

# Load multi-scale TFBS predictions
multiScalePred <- read.table(paste0("../../../data/TFBSPrediction/", projectName, "_pred_data.tsv"))
multiScaleRanges <- GRanges(multiScalePred$range)
strand(multiScaleRanges) <- "*"
multiScaleRanges$score <- multiScalePred$predScore
multiScaleRanges$TF <- multiScalePred$TF
TFs <- unique(multiScaleRanges$TF)

# Map chromBPNet prediction scores to the same sites
mode <- "count"
chromBPNetRanges <- multiScaleRanges
chromBPNetRanges$score <- 0
motifCREOv <- findOverlaps(chromBPNetRanges, regions, type = "within")
relativePos <- start(chromBPNetRanges[motifCREOv@from]) - start(regions[motifCREOv@to]) + 1
chromBPNetRanges[motifCREOv@from]$score <- scoreMtx[[mode]][cbind(motifCREOv@to, relativePos)]    

###############################
# Evalutate model performance #
###############################

benchmarkResults <- t(pbmcapply::pbmcmapply(
  function(TF){
    
    footprintRanges <- list()
    
    # Retrieve multi-scale prediction for the current TF
    footprintRanges[["multiScale"]] <- multiScaleRanges[multiScaleRanges$TF == TF]
    motifRanges <- footprintRanges[["multiScale"]]
    
    # Retrive chromBPNet predictions
    footprintRanges[["chromBPNet"]] <- subsetByOverlaps(chromBPNetRanges, motifRanges)
    
    # Map prediction scores by different methods to the same set of sites
    footprintRangesRemap <- lapply(
      names(footprintRanges),
      function(method){
        remap <- resize(motifRanges, 1, fix = "center")
        remap$score <- 0
        ov <- findOverlaps(remap, footprintRanges[[method]])
        remap[ov@from]$score <- footprintRanges[[method]][ov@to]$score
        remap$score
      }
    )
    names(footprintRangesRemap) <- names(footprintRanges)
    
    # Get ground truth binding label 
    bindingLabels <- rep(0, length(motifRanges))
    bindingLabels[findOverlaps(motifRanges, TFChIPRanges[[TF]])@from] <- 1
    
    nPred <- as.integer(length(bindingLabels) * 0.1)
    multiScalePrecision <- mean(bindingLabels[rank(-footprintRangesRemap$multiScale) <= nPred])
    chromBPnetPrecision <- mean(bindingLabels[rank(-footprintRangesRemap$chromBPNet) <= nPred])
    motifPrecision <- mean(bindingLabels)
    
    c(multiScalePrecision, chromBPnetPrecision, motifPrecision, nPred) 
    
  },
  TFs,
  mc.cores = 16
))
colnames(benchmarkResults) <- c("MultiScale", "chromBPNet", "Motif", "nPred")
cluster1TFs <- intersect(TFClustering$TF[TFClustering$cluster == 1], TFs)
colMeans(benchmarkResults)
colMeans(benchmarkResults[cluster1TFs,])
colMedians(benchmarkResults)
colMedians(benchmarkResults[cluster1TFs,])

#####################
# Visualize results #
#####################

system("mkdir ../../../data/sequenceModels/plots")
methods = c("MultiScale", "chromBPNet", "Motif")
plotData <- data.frame(
  precision = colMedians(as.matrix(benchmarkResults[, methods])),
  method = factor(methods, levels = methods))
pdf(paste0("../../../data/sequenceModels/plots/K562_chormBPNet_benchmark_", mode, ".pdf"), width = 3.5, height = 5)
ggplot(plotData) +
  geom_col(aes(x = method, y = precision), width = 0.75,
           position = position_stack(reverse = TRUE)) +
  xlab("") + ylab("Precision") + ggtitle("All TFs") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
dev.off()

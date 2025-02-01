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
regions <- readRDS("../../../data/K562/regionRanges.rds")

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

#################################################################
# Compare chromBPNet scores at bound and unbound TF motif sites #
#################################################################

# Load TF motif matches
motifPositions <- readRDS("../../../data/K562/TFMotifRanges.rds")

# Load unibind ChIP validated binding sites
TFChIPRanges <- readRDS("../../../data/shared/unibind/K562ChIPRanges.rds")

# Determine the radius around TF motifs we want to plot
contextRadius <- 50

# Select a TF for visualization
TF <- "YY1"
mode <- "profile"

# Overlap motif sites with CRE regions
TFMotifSites <- motifPositions[[TF]]
motifRegionOv <- findOverlaps(
  resize(TFMotifSites, 2 * contextRadius + 2, fix = "center"), 
  regions, type = "within")

# Calculate the relative position of each motif site in the corresponding CRE region
motifCenters <- resize(TFMotifSites[motifRegionOv@from], 1, fix = "center")
pos <- start(motifCenters) - start(regions[motifRegionOv@to]) + 1

# Retrieve chormBPNet scores within the +/- contextRadius window centered around each motif site
motifWindow <- c(outer((-contextRadius + 1):contextRadius, pos, "+"))
centeredScores <- as.matrix(Matrix::sparseMatrix(
  i = rep(1:length(motifRegionOv), each = 2 * contextRadius),
  j = rep(1:(2*contextRadius), length(motifRegionOv)),
  x = scoreMtx[[mode]][cbind(rep(motifRegionOv@to, each = 2 * contextRadius), motifWindow)]
))

# Randomly sample 500 bound sites and 500 unbound sites
boundSiteInds <- findOverlaps(TFMotifSites[motifRegionOv@from], TFChIPRanges[[TF]])@from
boundLabels <- 1:length(motifRegionOv) %in% boundSiteInds
sampleInds <- c(sample(which(boundLabels == 1), 500, replace = T),
               sample(which(boundLabels == 0), 500, replace = T))
boundLabels <- c("Unbound", "Bound")[(boundLabels + 1)[sampleInds]]

# Visualize results on a heatmap
plotMtx <- centeredScores[sampleInds, ]
colors <- colorRamp2(seq(0, quantile(plotMtx, 0.95), length.out=9),
                     colors = colorRampPalette(c(rep("white", 2), "#9ECAE1", "#08519C", "#08306B"))(9))
rowOrder <- order(rowMaxs(plotMtx[, (contextRadius - 5):(contextRadius + 5)]), decreasing = T)
pdf(paste("../../../data/sequenceModels/plots/K562_chromBPNet", TF, mode, "scores.pdf", sep = "_"),
    width = 3, height = 5)
Heatmap(
  plotMtx,
  col = colors,
  row_split = factor(boundLabels[rowOrder], levels = c("Unbound", "Bound")),
  border = TRUE,
  cluster_columns = F,
  column_title = TF,
  column_title_side = "top",
  row_order = rowOrder)
dev.off()



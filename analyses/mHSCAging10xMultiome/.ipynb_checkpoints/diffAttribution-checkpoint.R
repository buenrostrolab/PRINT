# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getTFBS.R")
source("../../code/getTargetPathway.R")
library("ComplexHeatmap")

###################
# Load input data #
###################

# Load region ranges
regions <- readRDS("../../data/mHSCAging10xV3/regionRanges.rds")

# Merge overlapping regions
regionsMerged <- reduce(resize(regions, 1000, fix = "center"))

# Load sequence attribution scores and reformat into a position-by-sample matrix
cellStates <- c("Old_Mk_biased", "Young_Mk_biased")
attScores <- NULL
bwDir <- "/data/PRINT/multiScaleFootprinting/data/sequenceModels/scPrinter/mHSCAging10x/"
regionWidth <- 1000
for(cellState in cellStates){
  for(fold in 0:0){
    bw <- rtracklayer::import.bw(paste0(bwDir, cellState, "_fold_", fold, "_TFscore.bigwig"))
    ov <- findOverlaps(bw, regionsMerged)
    scoreVec <- rep(0, sum(width(regionsMerged))) # Create a n_region * region_width empty vector
    relativePos <- start(bw[ov@from]) - start(regionsMerged[ov@to]) + 1 # Relative position of each base in each region
    cumulativeWidth <- c(0, cumsum(width(regionsMerged))) # Cumulative width up until each region in regionsMerged
    vecInd <- cumulativeWidth[ov@to] + relativePos # Find index of each base in scoreVec
    scoreVec[vecInd] <- bw$score[ov@from]
    attScores <- cbind(attScores, scoreVec)
    rm(bw)
    rm(scoreVec)
    rm(relativePos)
    rm(vecInd)
    gc()
  }
} 
colnames(attScores) <- cellStates

# Normalize within each sample. This is because the overall signal level could be 
# affected by sequencing depth etc.
attScores <- scale(attScores)

# Get a GRanges object of all single base-pair positions in CREs
positions <- data.table::rbindlist(pbmcapply::pbmclapply(
  1:length(regionsMerged),
  function(regionInd){
    region <- regionsMerged[regionInd]
    data.frame(
      chr = as.character(seqnames(region)), 
      pos = start(region):end(region))
  },
  mc.cores = 8
))
posRanges <- GRanges(seqnames = positions$chr, ranges = IRanges(start = positions$pos, end = positions$pos))
gc()

# Generate a position-by-sample matrix of sequence attribution scores
attSE <- SummarizedExperiment(
  assays = list(attribution = attScores),
  rowRanges = posRanges
)

diffAtt <- attScores[, 1] - attScores[, 2]
meanSignal <- rowMeans(attScores)

sampleInd <- sample(1:dim(attScores)[1], 5000)
qplot(meanSignal[sampleInd], diffAtt[sampleInd])

##########################
# Get TF motif positions #
##########################

# Load PWM data
cisBPMotifs <- readRDS("/data/PRINT/multiScaleFootprinting/data/shared/cisBP_mouse_pwms_2021.rds")

# Find motif matches for all TFs
motifPath <- "/data/PRINT/multiScaleFootprinting/data/mHSCAging10xV3/TFMotifRanges.rds"
if(!file.exists(motifPath)){
  TFMotifRanges <- pbmcapply::pbmclapply(
    names(cisBPMotifs),
    function(TF){
      motifmatchr::matchMotifs(cisBPMotifs[TF], 
                               regions, 
                               genome = "mm10",
                               out = "positions",
                               p.cutoff = 1e-4)[[1]]
    },
    mc.cores = 16
  )
  names(TFMotifRanges) <- names(cisBPMotifs)
  saveRDS(TFMotifRanges, motifPath)
}else{
  TFMotifRanges <- readRDS(motifPath)
}

##########################
# Get TF motif positions #
##########################

bwOldPath <- "/data/PRINT/multiScaleFootprinting/data/sequenceModels/scPrinter/mHSCAging10x/Old_Mk_biased_fold_0_TFscore.bigwig"
bwYoungPath <- "/data/PRINT/multiScaleFootprinting/data/sequenceModels/scPrinter/mHSCAging10x/Young_Mk_biased_fold_0_TFscore.bigwig"
testResults <- data.table::rbindlist(pbmcapply::pbmclapply(
  c("Atf4", "Atf6b", "Bcl11a", "Clock", "Ctcf", "E2f1", 
    "Egr1", "Ets1", "Fli1", "Jun", "Klf6", "Mycn", 
    "Nfe2", "Nfe2l2", "Nrf1", "Pbx3", "Smad1", "Spi1", 
    "Tfec", "Xbp1"), 
  function(TF){
    
    # Load young and old TF scores at motif sites
    selection <- rtracklayer::BigWigSelection(TFMotifRanges[[TF]])
    oldScores <- rtracklayer::import.bw(bwOldPath, selection = selection)$score
    youngScores <- rtracklayer::import.bw(bwYoungPath, selection = selection)$score
    
    # Calculate difference between young and old
    diff <- oldScores - youngScores
    
    # Only keep top differential sites
    diff <- diff[abs(diff) > quantile(abs(diff), 0.99)]
    
    # One-sample t-test to see if the scores are going up or down
    test <- t.test(diff)
    data.frame(
      TF = TF,
      meanDiff = test$estimate,
      pval = test$p.value
    )
  },
  mc.cores = 8
))

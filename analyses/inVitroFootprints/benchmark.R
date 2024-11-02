# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(GenomicRanges)

##################################################
# First find motif matches within our BAC region #
##################################################

# Load TF motifs
cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")

# Define BAC region
BACRange <- GRanges("chr2:238036072-238243952")

# Match motifs
motifSites <- motifmatchr::matchMotifs(
  pwms = cisBPMotifs[c("CEBPA", "MAX", "MYC")],
  subject = BACRange,
  out = "position",
  genome = "hg38",
  p.cutoff = 5e-5
)

# Find sites that are both MYC and MAX motif matches
motifSites$MYCMAX <- subsetByOverlaps(motifSites$MAX, motifSites$MYC)

# Resize to 1 bp
motifSites$CEBPACenters <- resize(motifSites$CEBPA, 1, fix = "center")
motifSites$MYCMAXCenters <- resize(motifSites$MYCMAX, 1, fix = "center")
centerPos <- list(
  "CEBPA" = start(motifSites$CEBPACenters) - start(BACRange) + 1,
  "MYCMAX" = start(motifSites$MYCMAXCenters) - start(BACRange) + 1
)

# Find background sites
bgSites <- list()
bgSites[["CEBPA"]] <- setdiff(
  1:width(BACRange), 
  findOverlaps(tile(BACRange, width = 1)[[1]], motifSites$CEBPA)@from
)
bgSites[["MYCMAX"]] <- setdiff(
  1:width(BACRange), 
  findOverlaps(tile(BACRange, width = 1)[[1]], c(motifSites$MAX, motifSites$MYC))@from
)

#######################
# Benchmark with HINT #
#######################

method <- "HINT"
samples <- c("Control", "CEBPA", "MYCMAX")
modelPred <- lapply(
  samples, 
  function(sample){
    if(method == "HINT"){
      pred <- import.bed(paste0("../../data/inVitroFootprints/HINT/", sample, ".bed"))
      pred <- resize(bed, 1, fix = "center")
    }else if(method == "TOBIAS"){
      pred <- import.bw(paste0("../../data/inVitroFootprints/TOBIAS/", sample, "/", sample, "_footprints.bw"))
      pred <- subsetByOverlaps(pred, BACRange)
    }
    scoreVec <- rep(0, width(BACRange))
    scoreVec[start(pred) - start(BACRange) + 1] <- pred$score
    scoreVec
  }
)
names(modelPred) <- samples

# Visualize results
TF <- "CEBPA"
plotData <- data.table::rbindlist(list(
  data.frame(group = paste0(TF, "\nmotif TF(+)"), scores = modelPred[[TF]][centerPos[[TF]]]),
  data.frame(group = paste0(TF, "\nbg TF(+)"), scores = modelPred[[TF]][bgSites[[TF]]]),
  data.frame(group = paste0(TF, "\nmotif TF(-)"), scores = modelPred[["Control"]][centerPos[[TF]]]),
  data.frame(group = paste0(TF, "\nbg TF(-)"), scores = modelPred[["Control"]][bgSites[[TF]]])
))
ggpubr::ggboxplot(plotData, x = "group", y = "scores", title = method, xlab = "")

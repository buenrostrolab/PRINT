# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(dplyr)
library(GenomicRanges)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(patchwork)
library(RColorBrewer)
library(patchwork)
source("../../code/utils.R")
source("../../code/getFootprints.R")

###################
# Load input data #
###################

# Load background dispersion model
dispersionModel <- list()
for(kernelSize in 2:100){
  dispersionModel[[as.character(kernelSize)]] <-
    readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Load PWM data
cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")

#################################################################################################
# For TFs where the bias will cause a footprint, examine whether we can detect real footprints  #
#################################################################################################

# Determine parameters
windowRadius <- 1000
localRadius <- 500
footprintRadius <- 20

TFList <- stringr::str_match(list.files("../../data/HepG2/scFootprints/"),"Bias(.+?).rds")[, 2]
TFList <- TFList[!is.na(TFList)]

results <- pbmcapply::pbmclapply(
  TFList,
  function(TF){
    
    # Load single cell aggregate Tn5 insertion and bias
    aggCuts <- as.matrix(t(readRDS(paste0("../../data/HepG2/scFootprints/scAggCuts", TF, ".rds"))))
    aggBias <- as.matrix(t(readRDS(paste0("../../data/HepG2/scFootprints/scAggBias", TF, ".rds"))))
    
    # Only analyze the region around the motif
    subsetInds <- (windowRadius - localRadius):(windowRadius + localRadius)
    aggCuts <- colSums(aggCuts[, subsetInds])
    aggBias <- colMeans(aggBias[, subsetInds])
    
    obsFootprint <- -log10(footprintScoring(
      Tn5Insertion = aggCuts,
      Tn5Bias = aggBias,
      dispersionModel = dispersionModel[[as.character(footprintRadius)]],
      footprintRadius = footprintRadius,
      flankRadius = footprintRadius
    ))
    
    biasFootprint <- -log10(footprintScoring(
      Tn5Insertion = rmultinom(1, size = sum(aggCuts), prob = aggBias / sum(aggBias)),
      Tn5Bias = rep(1, length(aggCuts)),
      dispersionModel = dispersionModel[[as.character(footprintRadius)]],
      footprintRadius = footprintRadius,
      flankRadius = footprintRadius
    ))
    
    # Retrieve score at the middle point
    midPoint <- as.integer(length(obsFootprint)/2)
    obsScore <- max(obsFootprint[(midPoint - footprintRadius):(midPoint + footprintRadius)])
    biasScore <- max(biasFootprint[(midPoint - footprintRadius):(midPoint + footprintRadius)])
   
    list(
      insertion = aggCuts,
      bias = aggBias,
      obsFootprint = obsFootprint,
      biasFootprint = biasFootprint,
      obsScore = obsScore,
      biasScore = biasScore,
      TF = TF
    ) 
  },
  mc.cores = 16
)
names(results) <- TFList

#####################
# Visualize results #
#####################

plotData <- as.data.frame(data.table::rbindlist(
  lapply(
    results,
    function(x){
      data.frame(
        obsScore = x$obsScore,
        biasScore = x$biasScore,
        TF = x$TF
      )
    }
  )
))

# Visualize results
labels <- rep("", length(TFList))
labelInds <- (plotData$obsScore > 1) | (plotData$biasScore > 1)
labels[labelInds] <- TFList[labelInds]
system("mkdir ../../data/HepG2/plots/motifBias")
pdf("../../data/HepG2/plots/motifBias/comparison.pdf", width = 7, height = 7)
ggplot(plotData) +
  geom_point(aes(x = obsScore, y = biasScore), size = 0.5, color = "brown") +
  xlab("Observed footprint score") + ylab("Bias footprint score") +
  geom_vline(xintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  ggrepel::geom_text_repel(
    aes(x = obsScore, y = biasScore), 
    label = labels, max.overlaps = 30) +
  theme_classic()
dev.off()

############################
# Visualize individual TFs #
############################

TF <- "ZNF613"
positions = -localRadius:localRadius
plotData <- data.frame(
  x1 = positions - 2,
  x2 = positions + 2,
  positions = positions,
  insertion = results[[TF]]$insertion,
  bias = results[[TF]]$bias,
  obsFootprint = results[[TF]]$obsFootprint,
  biasFootprint = results[[TF]]$biasFootprint
)
plotData <- plotData[abs(plotData$positions) <= 250,]

# Plot Tn5 insertion counts
p1 <- ggplot(plotData) + 
  geom_line(aes(x = positions, y = insertion), 
            color = "#006838", size = 1) +
  ggtitle(TF) +
  ylim(0, max(plotData$insertion) * 1.1) + xlab("") + ylab("Tn5 insertion") +
  theme_classic() 

p2 <- ggplot(plotData) + 
  geom_line(aes(x = positions, y = bias), 
            color = "#006838", size = 1) +
  ylim(0, max(plotData$bias) * 1.1) + xlab("") + ylab("predicted bias") +
  theme_classic() 

# Plot observed - expected deviation of Tn5 insertion
expectedRatio <- conv(plotData$bias, 20) / conv(plotData$bias, 40)
observedRatio <- conv(plotData$insertion, 20) / conv(plotData$insertion, 40)
deviation <- observedRatio - expectedRatio
p3 <- ggplot() +
  geom_ribbon(aes_string(x = plotData$positions, ymin = 0, 
                         ymax = deviation), 
              alpha = 0.5) + 
  ylab("HepG2\nObs-Exp\ndeviation") +
  xlab("") +
  theme_classic()

p4 <- ggplot(plotData) +
  geom_ribbon(aes(x = positions, ymin = 0, ymax = obsFootprint), 
              alpha = 0.5) + 
  scale_y_continuous(limits = c(0, max(plotData$obsFootprint)), expand = c(0, 0)) + 
  ylab("HepG2\nground truth\ncorrected footprints") +
  xlab("") +
  theme_classic()

pdf(paste0("../../data/HepG2/plots/motifBias/", TF, "_tracks.pdf"), width = 7, height = 7)
p1 / p2 / p3 / p4
dev.off()

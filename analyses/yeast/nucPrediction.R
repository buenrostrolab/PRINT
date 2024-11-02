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
source("../../code/getModelPrediction.R")

######################################
# Predict nuclesome occupancy signal #
######################################

# Initialize a footprintingProject object
projectName <- "yeast"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "sacCer3")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Set the regionRanges slot
regionRanges(project) <- readRDS("../../data/yeast/regionRanges.rds")
regions <- regionRanges(project)

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getPrecomputedBias(project, nCores = 24)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load barcodes for each pseudo-bulk
barcodeGroups <- data.frame(
  barcode = "yeast",
  group = 1L
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))
groupCellType(project) <- c(projectName)

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- paste0(projectDataDir, "merged.frags.gz")
#project <- getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T)

# Load background dispersion model
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <-
    readRDS(paste0(projectMainDir,"data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Load nucleosome prediction model
h5Path <- "../../data/yeast/nucleosome_model.h5"
predictionModel(project) <- loadModel(h5Path)

# Load chemically mapping data
nucTracks <- readRDS("../../data/yeast/nucTracks.rds")

# Calculate region-by-position matrix of ATAC counts
project <- getATACTracks(project)

#########################################
# Compute genome-wide nucleosome scores #
#########################################

# Predict nucleosome binding
getPrediction(project,
              nCores = 16)

######################################
# Load NucleoATAC prediction results #
######################################

nucleoATACSignal <- rtracklayer::import("../../data/yeast/nucleoATAC/yeastCRE.nucleoatac_signal.bedgraph.gz",
                                        format = "bedgraph")
nucleoATACTracks <- t(pbmcapply::pbmcmapply(
  function(regionInd){
    region <- regions[regionInd]
    regionTrack <- rep(0, width(region))
    regionPos <- slidingWindows(region, width = 1, step = 1)[[1]]
    regionSignal <- subsetByOverlaps(nucleoATACSignal, region)
    ov <- findOverlaps(regionPos, regionSignal)
    regionTrack[ov@from] <- regionSignal[ov@to]$score
    regionTrack
  },
  1:length(regions),
  mc.cores = 16
))

###############################
# Visualize a specific region #
###############################

library(patchwork)

regionInd <- 1370

# Compute predicted nucleosome occupancy score
scores <- getRegionPredMatrix(project, regionInd)

positions <- start(regions[regionInd]):end(regions[regionInd])
plotData <- data.frame(
  Tn5Insertion = ATACTracks(project)[regionInd, ],
  x1 = positions - 1,
  x2 = positions + 2,
  baseline = 0,
  scores = conv(scores, 10)/20,
  nucleoATAC = nucleoATACTracks[regionInd, ],
  occupancy = nucTracks[regionInd, ]
)

p1 <- ggplot(plotData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = baseline, ymax = Tn5Insertion), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Tn5\ninsertion") +
  theme_classic() + 
  ggtitle(regions[regionInd])

p2 <- ggplot(plotData) + 
  geom_line(aes_string(x = "x1", y = "nucleoATAC"), 
            color = "#69ACD5") + xlab("") + ylab("nucleoATAC\nscores") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

p3 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "x1", ymin = min(plotData$scores), ymax = "scores"), 
              fill = "#69ACD5") + xlab("") + ylab("ATAC\nnucleosme\nfootprints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

p4 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "x1", ymin = 0, ymax = "occupancy"), 
              fill = "#69ACD5") + xlab("") + ylab("Chemically mapped\nnucleosome\nfootprints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

pdf(paste0("../../data/yeast/plots/region_", regionInd, ".pdf"), width = 10, height = 5)
p1 / p2 / p3 / p4
dev.off()

###########################
# Examine all CRE regions #
###########################

CRERanges <- readRDS("../../data/yeast/ATACPeaks.rds")

# Load whole-genome nucleosome predictions
footprintResults <- Reduce(c, lapply(
  gtools::mixedsort(list.files("../../data/yeast/chunkedPredResults/")),
  function(file){readRDS(paste0("../../data/yeast/chunkedPredResults/", file))}
))

# Get region-by-single bp matrix of predicted nucleosome scores
footprintScores <- t(sapply(
  footprintResults,
  function(x){rep(x$predScores, each = 10)}
))

# Resize regions to remove 100 bp edges where we can't make predictions
regionsResized <- resize(regions, 4800, fix = "center")
nucTracksResized <- nucTracks[, 101:4900]
nucleoATACResized <- nucleoATACTracks[, 101:4900]
ov <- findOverlaps(CRERanges, regionsResized)

# Only evaluate performance on the test set
testChr <- c("chrVI", "chrXVI", "chrXIII", "chrIII", "chrXIV", "chrXV")
ov <- ov[as.character(seqnames(regions)[ov@to]) %in% testChr]

# Make a list of prediction results by the two methods
scores <- list(
  "footprint" = footprintScores,
  "nucleoATAC" = nucleoATACResized
)

# Calculate precision-recall curve for each method
nBreaks <- 30
distThreshold <- 75
performance <- lapply(
  names(scores),
  function(method){
    sapply(
      quantile(scores[[method]][scores[[method]] != 0], seq(0.03, 0.97, length.out = nBreaks)),
      function(sigThreshold){
        
        # For each threshold, examine CRE regions
        # Identify predicted and ground truth score summits and calculate their distance
        nearestNucDist <- pbmcapply::pbmclapply(
          1:length(ov),
          function(ovInd){
            
            # Find the region of interest
            CREInd <- ov@from[ovInd]
            regionInd <- ov@to[ovInd]
            CRE <- CRERanges[CREInd]
            region <- regionsResized[regionInd]
            CREStart <- max(start(region) - start(CRE) + 1, 1)
            regionStart <- max(start(CRE) - start(region) + 1, 1)
            CREEnd <- min(width(CRE) - end(CRE) + end(region), width(CRE))
            regionEnd <- min(width(region) - end(region) + end(CRE), width(region))
            
            # Extract predicted and ground truth scores in the region
            chemMapTrack <-  nucTracksResized[regionInd, regionStart:regionEnd]
            predTrack <- scores[[method]][regionInd, regionStart:regionEnd]
            realSummits <- findSummits(chemMapTrack, r = 50, threshold = 1)
            predSummits <- findSummits(predTrack, r = 50, threshold = sigThreshold)
            
            # Calculate distance
            if(length(predTrack) > 200){
              nearestNucDist <- list(
                # For each predicted nucleosome, calculate distance to nearest ground truth nucleosome
                "predDist" = colMins(abs(outer(realSummits, predSummits, FUN ='-'))),
                # For each ground truth nucleosome, calculate distance to nearest predicted nucleosome
                "targetDist" = rowMins(abs(outer(realSummits, predSummits, FUN ='-'))))
            }else{
              NULL
            }
          },
          mc.cores = 16
        )
        
        # Integrate results across CREs
        predDist <- Reduce(c, lapply(nearestNucDist, function(x){x[["predDist"]]}))
        targetDist <- Reduce(c, lapply(nearestNucDist, function(x){x[["targetDist"]]}))
        predDist <- predDist[is.finite(predDist) & (!is.na(predDist))]
        targetDist <- targetDist[is.finite(targetDist) & (!is.na(targetDist))]
        
        # Calculate precision and recall
        c(mean(predDist < distThreshold),
          mean(targetDist < distThreshold))
        
      }
    )
  }
)
names(performance) <- names(scores)

# Visualize results
plotData <- as.data.frame(rbind(t(performance$footprint[, 1:nBreaks]),
                                t(performance$nucleoATAC[, 1:nBreaks])))
colnames(plotData) <- c("Precision", "Recall")
plotData$Method <- c(rep("PRINT", nBreaks), rep("NucleoATAC", nBreaks))
pdf(paste0("../../data/yeast/plots/precisionRecall", distThreshold,"bpThreshold.pdf"),
    width = 6, height = 5)
ggplot(plotData) +
  geom_line(aes(x = Recall, y = Precision, group = Method, color = Method)) +
  theme_classic()
dev.off()

#############################################################
# Evaluate how depth and distance to CRE affect performance #
#############################################################

sigThreshold <- 0.6
distThreshold <- 75
method <- "footprint"

# Find CRE centers on each chromosome
chrCRECenters <- lapply(
  unique(seqnames(regions)),
  function(chr){
    chrCRE <- CRERanges[as.character(seqnames(CRERanges)) == chr]
    start(resize(chrCRE, 1,, fix = "center"))
  }
)
names(chrCRECenters) <- unique(seqnames(regions))

# For each threshold, examine CRE regions
# Identify predicted and ground truth score summits and calculate their distance
nearestNucDist <- data.table::rbindlist(pbmcapply::pbmclapply(
  1:length(regionsResized),
  function(regionInd){
    
    # Extract predicted and ground truth scores in the region
    chemMapTrack <-  nucTracksResized[regionInd, ]
    predTrack <- scores[[method]][regionInd, ]
    targetSummits <- findSummits(chemMapTrack, r = 50, threshold = 1)
    predSummits <- findSummits(predTrack, r = 50, threshold = sigThreshold)
    
    # Calculate ATAC coverage
    ATACTrk <- conv(ATACTracks(project)[regionInd, 101:4900], 200)/400
    
    #plt(ATACTracks(project)[regionInd, 101:4900]) / plt(predTrack) / plt(chemMapTrack)
    
    # Calculate distance to closest CRE
    region <- regions[regionInd]
    predCoords <- start(region) + predSummits - 1
    predCREDist <- colMins(abs(outer(chrCRECenters[[as.character(seqnames(region))]], predCoords, FUN ='-')))
    targetCoords <- start(region) + targetSummits - 1
    targetCREDist <- colMins(abs(outer(chrCRECenters[[as.character(seqnames(region))]], targetCoords, FUN ='-')))
    
    # For each predicted nucelosome, record the distance to nearest ground truth nucleosoem, 
    # distance to nearest CRE, and local coverage
    if(length(predSummits) > 0){
      predDf <- data.frame(
        dist = colMins(abs(outer(targetSummits, predSummits, FUN ='-'))),
        coverage = ATACTrk[predSummits],
        regionInd = regionInd,
        chr = as.character(seqnames(regions[regionInd])),
        CREDist = predCREDist,
        position = predSummits,
        type = "pred"
      )
    }else{
      predDf <- NULL
    }
    
    # For each ground truth nucelosome, record the distance to nearest predicted nucleosoem, 
    # distance to nearest CRE, and local coverage
    if(length(targetSummits) > 0){
      targetDf <- data.frame(
        dist = rowMins(abs(outer(targetSummits, predSummits, FUN ='-'))),
        coverage = ATACTrk[targetSummits],
        regionInd = regionInd,
        chr = as.character(seqnames(regions[regionInd])),
        CREDist = targetCREDist,
        position = targetSummits,
        type = "target"
      )
    }else{
      targetDf <- NULL
    }
    
    rbind(predDf, targetDf)
  },
  mc.cores = 16
))

nearestNucDist <- nearestNucDist[nearestNucDist$chr %in% testChr,]

################################################
# Examine how local coverage affects precision #
################################################

# Integrate results across CREs
comparison <- nearestNucDist[nearestNucDist$type == "pred",]

depthBins <- quantile(comparison$coverage, seq(0.0,1, 0.1))
plotData <- data.table::rbindlist(lapply(
  1:(length(depthBins) - 1),
  function(i){
    filter <- (comparison$coverage >= depthBins[i]) & (comparison$coverage < depthBins[i + 1])
    data.frame(
      depthBins = mean(comparison$coverage[filter]),
      precision = mean(comparison$dist[filter] < distThreshold)
    )
  }
))

pdf(paste0("../../data/yeast/plots/precision_coverage_relationship.pdf"),
    width = 5, height = 5)
ggplot(plotData) +
  geom_line(aes(x = depthBins, y = precision), color = "grey", size = 1.5) +
  scale_x_continuous(trans = "log10") +
  xlab("Average local coverage") +
  ylab("Precision") +
  theme_classic()+
  ylim(0.4, 1)
dev.off()

###################################################################
# Examine how the distance to the nearest CRE affects performance #
###################################################################

comparison <- nearestNucDist[nearestNucDist$type == "pred",]

distBins <- c(100,200,300,400,500,600,700,1000,1500,2000,3000)
plotData <- data.table::rbindlist(lapply(
  1:(length(distBins) - 1),
  function(i){
    filter <- (comparison$CREDist >= distBins[i]) & (comparison$CREDist < distBins[i + 1])
    data.frame(
      distBins = mean(comparison$CREDist[filter]),
      precision = mean(comparison$dist[filter] < distThreshold)
    )
  }
))

pdf(paste0("../../data/yeast/plots/precision_CREDist_relationship.pdf"),
    width = 5, height = 5)
ggplot(plotData) +
  geom_line(aes(x = distBins, y = precision), color = "grey", size = 1.5) +
  xlab("Distance to nearest CRE (bp)") +
  ylab("Precision") +
  theme_classic() +
  ylim(0.4, 0.95)
dev.off()

#############################################
# Examine how local coverage affects recall #
#############################################

# Integrate results across CREs
comparison <- nearestNucDist[nearestNucDist$type == "target",]

depthBins <- quantile(comparison$coverage, seq(0.0,1, 0.1))
plotData <- data.table::rbindlist(lapply(
  1:(length(depthBins) - 1),
  function(i){
    filter <- (comparison$coverage >= depthBins[i]) & (comparison$coverage < depthBins[i + 1])
    data.frame(
      depthBins = mean(comparison$coverage[filter]),
      recall = mean(comparison$dist[filter] < distThreshold)
    )
  }
))

pdf(paste0("../../data/yeast/plots/recall_coverage_relationship.pdf"),
    width = 5, height = 5)
ggplot(plotData) +
  geom_line(aes(x = depthBins, y = recall), color = "grey", size = 1.5) +
  scale_x_continuous(trans = "log10") +
  xlab("Local ATAC coverage") +
  ylab("Recall") +
  ylim(0.4, 0.95) +
  theme_classic() 
dev.off()

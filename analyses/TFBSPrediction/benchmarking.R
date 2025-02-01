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
library(hdf5r)
require(PRROC)
library(GenomicRanges)

###################
# Load input data #
###################

projectName <- "K562"

# Either "Union" or "Intersect". This is how multiple datasets of the same TF
# are integrated into a consensus set
integration <- "Union" 

# Whether to test model I or model II
model <- "I"

# Whether to evaluate precision at top 10% sites (metric = "Precision") or AUPRC (metric = "AUPRC")
metric <- "AUPRC"

# Load TF ChIP ranges
if(integration == "Intersect"){
  TFChIPRanges <- readRDS(paste0("../../data/shared/unibind/", projectName, "ChIPRanges.rds"))
}else if(integration == "Union"){
  TFChIPRanges <- readRDS(paste0("../../data/shared/unibind/", projectName, "UnionChIPRanges.rds"))
}

# Load chain file required for lift-over
pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Load TF cluster membership
TFClustering <- read.csv("../../data/TFBSPrediction/clusterLabelsAllTFs.txt", sep = "\t")

# Load multi-scale TFBS predictions
if(model == "I"){
  multiScalePred <- read.table(paste0("../../data/TFBSPrediction/", projectName, "_pred_data_cluster_I.tsv"))
}else if(model == "II"){
  multiScalePred <- read.table(paste0("../../data/TFBSPrediction/", projectName, "_pred_data.tsv"))
}
multiScaleRanges <- GRanges(multiScalePred$range)
strand(multiScaleRanges) <- "*"
multiScaleRanges$score <- multiScalePred$predScore
multiScaleRanges$TF <- multiScalePred$TF

# Load DNaseI footprints
DNaseIFootprints <- list()
DNaseIFootprints[["ENCLB253REF"]] <- rtracklayer::import.bw("../../data/DNaseIFootprint/vierstra/interval.all.winlnpval.K562-DS15363.bw")
DNaseIFootprints[["ENCLB253REF"]] <- subsetByOverlaps(DNaseIFootprints[["ENCLB253REF"]], multiScaleRanges)
DNaseIFootprints[["ENCLB843GMH"]] <- rtracklayer::import.bw("../../data/DNaseIFootprint/vierstra/interval.all.winlnpval.K562-DS16924.bw")
DNaseIFootprints[["ENCLB843GMH"]] <- subsetByOverlaps(DNaseIFootprints[["ENCLB843GMH"]], multiScaleRanges)
DNaseIFootprints[["ENCLB096YUZ"]] <- rtracklayer::import.bw("../../data/DNaseIFootprint/vierstra/interval.all.winlnpval.h.K562-DS52908.bw")
DNaseIFootprints[["ENCLB096YUZ"]] <- subsetByOverlaps(DNaseIFootprints[["ENCLB096YUZ"]], multiScaleRanges)

# Load HINT prediction results
HINTBed <- read.table(paste0("../../data/", projectName, "/HINT/footprints.bed"))
HINTRanges <- GRanges(seqnames = HINTBed$V1,
                      ranges = IRanges(start = HINTBed$V2, end = HINTBed$V3))
HINTRanges$score <- HINTBed$V5
seqlevelsStyle(HINTRanges) <- "UCSC"  # necessary
HINTRanges <- unlist(rtracklayer::liftOver(HINTRanges, ch))

# Get the list of TFs
TOBIASFolder <- paste0(paste0("../../data/", projectName, "/Tobias/prediction"))
TOBIASTFs <- unname(sapply(list.files(TOBIASFolder), function(x){strsplit(x, "_")[[1]][1]}))
TFs <- Reduce(intersect, list(TOBIASTFs, unique(multiScaleRanges$TF), names(TFChIPRanges)))
cluster1TFs <- intersect(TFClustering$TF[TFClustering$cluster == 1], TFs)

##############################
# Evaluate model performance #
##############################

benchmarkResults <- t(pbmcapply::pbmcmapply(
  function(TF){
    
    footprintRanges <- list()
    
    # Retrieve multi-scale prediction for the current TF
    footprintRanges[["multiScale"]] <- multiScaleRanges[multiScaleRanges$TF == TF]
    motifRanges <- footprintRanges[["multiScale"]]
    
    # Load TOBIAS prediction scores
    TOBIASBed <- read.table(paste0("../../../multiScaleFootprinting/data/", projectName, "/Tobias/prediction/", TF, 
                                   "_motif/beds/", TF, "_motif_all.bed"))
    footprintRanges[["TOBIAS"]] <- GRanges(seqnames = TOBIASBed$V1,
                                           ranges = IRanges(start = TOBIASBed$V2, end = TOBIASBed$V3))
    footprintRanges[["TOBIAS"]]$score <- TOBIASBed$V10
    
    # Lift over from hg19 to hg38
    seqlevelsStyle(footprintRanges[["TOBIAS"]]) <- "UCSC"  # necessary
    footprintRanges[["TOBIAS"]] <- unlist(rtracklayer::liftOver(footprintRanges[["TOBIAS"]], ch))
    
    # Retrieve DNaseI-based prediction for the current TF
    for(dataset in names(DNaseIFootprints)){
      footprintRanges[[paste0("DNaseI-", dataset)]] <- DNaseIFootprints[[dataset]]
    }
    
    # Retrive HINT footprint ranges
    footprintRanges[["HINT"]] <- subsetByOverlaps(HINTRanges, motifRanges)
    
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
    
    bindingLabels <- rep(0, length(motifRanges))
    bindingLabels[findOverlaps(motifRanges, TFChIPRanges[[TF]])@from] <- 1
    
    if(metric == "Precision"){
      
      nPred <- as.integer(length(bindingLabels) * 0.1)
      multiScalePrecision <- mean(bindingLabels[rank(-footprintRangesRemap$multiScale) <= nPred])
      TOBIASPrecision <- mean(bindingLabels[rank(-footprintRangesRemap$TOBIAS) <= nPred])
      HINTPrecision <- mean(bindingLabels[rank(-footprintRangesRemap$HINT) <= nPred])
      DNaseIPrecision1 <- mean(bindingLabels[rank(-footprintRangesRemap$`DNaseI-ENCLB253REF`) <= nPred])
      DNaseIPrecision2 <- mean(bindingLabels[rank(-footprintRangesRemap$`DNaseI-ENCLB843GMH`) <= nPred])
      DNaseIPrecision3 <- mean(bindingLabels[rank(-footprintRangesRemap$`DNaseI-ENCLB096YUZ`) <= nPred])
      motifPrecision <- mean(bindingLabels)
      performance <- c(multiScalePrecision, DNaseIPrecision1, DNaseIPrecision2, DNaseIPrecision3,
                       HINTPrecision, TOBIASPrecision, motifPrecision, nPred) 
      
    }else if(metric == "AUPRC"){
      
      require(PRROC)
      AUPRC <- sapply(
        names(footprintRanges),
        function(method){
          fg <- footprintRangesRemap[[method]][bindingLabels == 1]
          bg <- footprintRangesRemap[[method]][bindingLabels == 0]
          pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
          pr$auc.integral
        }
      )
      performance <- AUPRC[c("multiScale", "DNaseI-ENCLB253REF", "DNaseI-ENCLB843GMH", 
                             "DNaseI-ENCLB096YUZ", "HINT", "TOBIAS")]
      
    }
    
    performance
  },
  TFs,
  mc.cores = 16
))
if(metric == "Precision"){
  colnames(benchmarkResults) <- c("PRINT", "DNaseI-ENCLB253REF", "DNaseI-ENCLB843GMH", 
                                  "DNaseI-ENCLB096YUZ", "HINT", "TOBIAS", "Motif", "nPred")
}else if(metric == "AUPRC"){
  colnames(benchmarkResults) <- c("PRINT", "DNaseI-ENCLB253REF", "DNaseI-ENCLB843GMH", 
                                  "DNaseI-ENCLB096YUZ", "HINT", "TOBIAS")
}

performance <- data.frame(
  AllTFMean = colMeans(benchmarkResults),
  AllTFMedian = colMedians(benchmarkResults),
  Cluster1Mean = colMeans(benchmarkResults[cluster1TFs,]),
  Cluster1Median = colMedians(benchmarkResults[cluster1TFs,])
)
performance 

#####################
# Visualize results #
#####################

# Get average performance for each method
if(metric == "Precision"){
  methods <- c("PRINT", "DNaseI-ENCLB253REF", "DNaseI-ENCLB843GMH", 
               "DNaseI-ENCLB096YUZ", "HINT", "TOBIAS", "Motif")
}else if(metric == "AUPRC"){
  methods <- c("PRINT", "DNaseI-ENCLB253REF", "DNaseI-ENCLB843GMH", 
               "DNaseI-ENCLB096YUZ", "HINT", "TOBIAS")
}
if(model == "I"){
  plotData <- data.frame(
    performance = colMeans(as.matrix(benchmarkResults[cluster1TFs, methods])),
    method = factor(methods, levels = methods))
}else if(model == "II"){
  plotData <- data.frame(
    performance = colMeans(as.matrix(benchmarkResults[, methods])),
    method = factor(methods, levels = methods))
}
plotData <- plotData[order(plotData$performance, decreasing = T),]
plotData$method <- factor(plotData$method, levels = plotData$method)

pdf(paste0("../../data/TFBSPrediction/plots/benchmark_", metric, "_Model_", model, "_", integration, "_TFBS.pdf"), 
    width = 3.5, height = 5)
ggplot(plotData) +
  geom_col(aes(x = method, y = performance), width = 0.75,
           position = position_stack(reverse = TRUE)) +
  xlab("") + ylab(metric) + ggtitle(list("II" = "All TFs", "I" = "Cluster 1 TFs")[[model]]) +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
dev.off()

# Save results to a file
write.table(benchmarkResults,
            paste0("../../data/TFBSPrediction/benchmark_", metric, "_Model_", model, "_", integration, "_TFBS.txt"),
            quote = F, sep = "\t")

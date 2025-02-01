# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../../../code/utils.R")
source("../../../../code/getCounts.R")
source("../../../../code/getBias.R")
source("../../../../code/getFootprints.R")
source("../../../../code/visualization.R")
source("../../../../code/getTFBS.R")
library(hdf5r)
require(PRROC)
library(GenomicRanges)

###################
# Load input data #
###################

projectName <- "K562"

# Either "Union" or "Intersect". This is how multiple datasets of the same TF
# are integrated into a consensus set
integration <- "Intersect" 

# Whether to test model I or model II
model <- "II"

# Whether to evaluate precision at top 10% sites (metric = "Precision") or AUPRC (metric = "AUPRC")
metric <- "Precision"

# Load TF ChIP ranges
if(integration == "Intersect"){
  TFChIPRanges <- readRDS(paste0("../../../../data/shared/unibind/", projectName, "ChIPRanges.rds"))
}else if(integration == "Union"){
  TFChIPRanges <- readRDS(paste0("../../../../data/shared/unibind/", projectName, "UnionChIPRanges.rds"))
}

# Load multi-scale TFBS predictions
if(model == "I"){
  multiScalePred <- read.table(paste0("../../../../data/TFBSPrediction/", projectName, "_pred_data_cluster_I.tsv"))
}else if(model == "II"){
  multiScalePred <- read.table(paste0("../../../../data/TFBSPrediction/", projectName, "_pred_data.tsv"))
}
multiScaleRanges <- GRanges(multiScalePred$range)
strand(multiScaleRanges) <- "*"
multiScaleRanges$score <- multiScalePred$predScore
multiScaleRanges$TF <- multiScalePred$TF

# Load scPrinter sequence model predictions
seqModelPred <- rtracklayer::import.bw("../../../../data/sequenceModels/scPrinter/K562/ground_truth_1bp_231112.bigwig")

# Filter to keep only test chromosomes
testChrs <- unique(as.character(seqnames(seqModelPred)))
multiScaleRanges <- multiScaleRanges[as.character(seqnames(multiScaleRanges)) %in% testChrs]

##############################
# Evaluate model performance #
##############################

TFs <- unique(multiScaleRanges$TF)
benchmark <- t(pbmcapply::pbmcmapply(
  function(TF){
    
    footprintRanges <- list()
    
    # Retrieve multi-scale prediction for the current TF
    footprintRanges[["TFModel"]] <- multiScaleRanges[multiScaleRanges$TF == TF]
    motifRanges <- footprintRanges[["TFModel"]]
    
    # Load TOBIAS prediction scores
    footprintRanges[["seqModel"]] <- seqModelPred

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
      TFModelPrecision <- mean(bindingLabels[rank(-footprintRangesRemap$TFModel) <= nPred])
      seqModelPrecision <- mean(bindingLabels[rank(-footprintRangesRemap$seqModel) <= nPred])
      motifPrecision <- mean(bindingLabels)
      performance <- c(TFModelPrecision, seqModelPrecision, motifPrecision, nPred) 
      
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
      performance <- AUPRC[c("TFModel", "seqModel")]
      
    }
    
    performance
  },
  TFs,
  mc.cores = 16
))
benchmark <- as.data.frame(benchmark)

if(metric == "Precision"){
  colnames(benchmark) <- c("TF Model", "sequenceModel", "Motif", "nPred")
  benchmark <- benchmark[, colnames(benchmark) != "nPred"]
}else if(metric == "AUPRC"){
  colnames(benchmark) <- c("TF Model", "sequenceModel")
}

performance <- data.frame(
  AllTFMean = colMeans(benchmark),
  AllTFMedian = colMedians(as.matrix(benchmark))
)
performance 

if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(rtracklayer)

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getFootprints.R")
source("../../code/getBias.R")
source("../../code/visualization.R")

#################################
# Load motif and footprint data #
#################################

# Load TF motifs
cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")

# Find motif matches
regions <- readRDS("../../data/K562/regionRanges.rds")
selectedTFs <- c("USF1", "NRF1", "YY1")
motifPositions <- motifmatchr::matchMotifs(
  pwms = cisBPMotifs[selectedTFs],
  subject = regions,
  genome = "hg38", 
  out = "positions")

# Load footprintingProject object
project <- readRDS("../../data/K562/project.rds")

# Load file for liftOver
pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Load pre-computed multi-scale footprints
footprintRadii <- c(10, 20, 30, 50, 80, 100)
multiScaleFootprints <- list()
for(footprintRadius in footprintRadii){
  
  print(paste0("Loading data for footprint radius = ", footprintRadius))  
  chunkDir <- paste0("../../data/K562/chunkedFootprintResults/", footprintRadius, "/")
  chunkFiles <- gtools::mixedsort(list.files(chunkDir))
  scaleFootprints <- pbmcapply::pbmclapply(
    chunkFiles,
    function(chunkFile){
      chunkData <- readRDS(paste0(chunkDir, chunkFile))
      as.data.frame(t(sapply(
        chunkData,
        function(regionData){
          regionData$aggregateScores
        }
      )))
    },
    mc.cores = 16
  )
  scaleFootprints <- data.table::rbindlist(scaleFootprints)
  multiScaleFootprints[[as.character(footprintRadius)]] <- as.matrix(scaleFootprints)
}

######################
# Load ChIP-exo data #
######################

TF <- "NRF1"

# Load TF ChIP-exo bigwig data
TFChIPExoDir <- "../../data/TFBSPrediction/ChIPExo/"
TFChIPExoFiles <- list.files(TFChIPExoDir)
TFChIPExoFiles <- TFChIPExoFiles[stringr::str_detect(TFChIPExoFiles, paste0(TF, "[A-Z\\d_]+?K562"))]

# If there are multiple files, merge the peak list
TFChIPExo <- mergeRegions(Reduce(c, lapply(
  TFChIPExoFiles,
  function(TFChIPExoFile){
    TFChIPExo <- import.bw(paste0(TFChIPExoDir, TFChIPExoFile))
    TFChIPExo <- TFChIPExo[TFChIPExo$score == 1]
  }
)))

# LiftOver from hg19 to hg38
seqlevelsStyle(TFChIPExo) <- "UCSC"  # necessary
TFChIPExo <- unlist(rtracklayer::liftOver(TFChIPExo, ch))

# Only focus on binding sites in CREs
TFChIPExo <- subsetByOverlaps(TFChIPExo, regions)

#############################################################
# Visualize ChIPExo-footprint comparison in a single region #
#############################################################

# Create folder for storing results
plotDir <- "../../data/TFBSPrediction/plots/ChIPExoValidation/"
if(!dir.exists(plotDir)){
  system(paste("mkdir", plotDir))
}

# Focus on sites with a motif match
TFBS <- subsetByOverlaps(motifPositions[[TF]], TFChIPExo)

# Select an example site
regionChIPOv <- findOverlaps(resize(regions, 800, fix = "center"), TFBS)
ovInd <- 70
regionInd <- regionChIPOv@from[ovInd]
ChIPInd <- regionChIPOv@to[ovInd]
region <- regions[regionInd]
ChIPSite <- TFBS[ChIPInd]

# Compute multi-scale footprints
p <- plotMultiScale(project, regionInd = regionInd)

# Label the position of the motif
pos <- start(resize(ChIPSite, 1, fix = "center")) - start(region) - 100
p@matrix[round(1:dim(p@matrix)[1]/5) %% 2 == 0, pos:(pos+2)] <- 5
pdf(paste0(plotDir, TF, "_region_", regionInd, ".pdf"), width = 7, height = 5)
p
dev.off()

############################################################
# Visualize ChIPExo-footprint comparison across the genome #
############################################################

plotRadius <- 100

# Extract multi-scale footprints centered around each TF binding site
TFBS <- motifPositions[[TF]][findOverlaps(motifPositions[[TF]], TFChIPExo)@from]
ChIPRegionOv <- findOverlaps(TFBS, resize(regions, 800, fix = "center"), type = "within")
singleSiteMultiScale <- t(pbmcapply::pbmcmapply(
  function(ovInd){
    regionInd <- ChIPRegionOv@to[ovInd]
    ChIPInd <- ChIPRegionOv@from[ovInd]
    region <- regions[regionInd]
    ChIPSite <- TFBS[ChIPInd]
    pos <- start(resize(ChIPSite, 1, fix = "center")) - start(region)
    Reduce(c, sapply(
      multiScaleFootprints,
      function(x){
        x[regionInd, (pos - plotRadius):(pos + plotRadius)]
      }
    ))
  },
  1:length(ChIPRegionOv),
  mc.cores = 16
))

# Split heatmap columns by kernel size
footprintRadii <- c(10, 20, 30, 50, 80, 100)
colGroups <- Reduce(c, lapply(
  footprintRadii,
  function(r){
    paste(rep(r, plotRadius * 2 + 1), "bp")
  }
))

# Visualize using heatmap
colors <- colorRamp2(seq(0, quantile(singleSiteMultiScale, 0.95), length.out=9),
                     colors = colorRampPalette(c(rep("white", 2),  
                                                 "#9ECAE1", "#08519C", "#08306B"))(9))
sampleInd <- sample(1:dim(singleSiteMultiScale)[1], min(500, length(ChIPRegionOv)))
rowOrder <- order(singleSiteMultiScale[sampleInd, 5 * plotRadius + 2], decreasing = T)

pdf(paste0(plotDir, TF, "_sites.pdf"), width = 8, height = 5)
Heatmap(singleSiteMultiScale[sampleInd[rowOrder],],
        col = colors,
        cluster_rows = F,
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        cluster_column_slices = F,
        column_split = factor(colGroups, levels = paste(footprintRadii, "bp")),
        name = paste(TF, "\nfootprint\nscore"),
        border = TRUE,
        column_title = "Kernel sizes",
        column_title_side = "bottom"
)
dev.off()



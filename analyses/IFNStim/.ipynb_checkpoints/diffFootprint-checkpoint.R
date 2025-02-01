# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/getGroupData.R")
source("../../code/getAggregateFootprint.R")
source("../../code/getTFBS.R")
source("../../code/visualization.R")

library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

###################
# Load input data #
###################

# Load object
project <- readRDS("../../data/IFNStim/project.rds")

# Load PWM data
motifs <- readRDS("../../data/shared/cisBP_mouse_pwms_2021.rds")

# Get genomic ranges of CREs
regions <- regionRanges(project)

# Get motif matched sites
motifSites <- motifmatchr::matchMotifs(
  motifs$Stat2, subject = regions, genome = "mm10", out = "positions")[[1]]

# Load diff ChIP data
ChIPInfo <- read.csv("../../data/IFNStim/ChIP/Supplement_5_S2E.txt", sep = "\t")
ChIPRanges <- GRanges(seqnames = ChIPInfo$chr, ranges = IRanges(start = ChIPInfo$start, end = ChIPInfo$end))
ChIPRanges$control <- ChIPInfo$Stat2.Signal.control.cells
ChIPRanges$IFN <- ChIPInfo$Stat2.signal.IFN.treated.cells
ChIPRanges$FC <- as.numeric(ChIPInfo$Stat2.FoldChange)
ChIPRanges$FC[is.na(ChIPRanges$FC)] <- 1

##########################
# Aggregate footprinting #
##########################

# Get region-by-position Tn5 insertion count matrix
samples <- groups(project)
ATACTrackList <- lapply(
  samples,
  function(group){
    project <- getATACTracks(project, groupIDs = group)
    ATACTracks(project)
  }
)
names(ATACTrackList) <- samples

for(group in samples){
  
  ATACTracks(project) <- ATACTrackList[[group]]
  
  selectedSites <- subsetByOverlaps(
    motifSites, resize(ChIPRanges, 100, fix = "center")[(ChIPRanges$FC > 10) & (ChIPRanges$IFN > 5)])
  
  # Get aggregate Tn5 insertion profile centered around motif sites for a specific TF
  aggregateFP <- getAggregateFootprint(
    project, 
    selectedSites
  )
  
  # For each scale size, calculate footprint scores for the aggregate profile
  footprintRadii <- 2:100
  width <- unique(regionWidth(project))
  plotRadius <- 500
  plotRange <- (width - plotRadius + 1):(width + plotRadius)
  aggregateMultiScale <- pbmcapply::pbmcmapply(
    function(footprintRadius){
      footprintPvals <- footprintScoring(
        Tn5Insertion = aggregateFP$uncorrectedATAC[plotRange],
        Tn5Bias = aggregateFP$Tn5Bias[plotRange],
        dispersionModel = dispModel(project, as.character(footprintRadius)),
        footprintRadius = footprintRadius,
        flankRadius = footprintRadius
      )
      # Remove regions affected by edge effect
      # For edges, the sum of bias and counts on the left / right flank might
      # be much lower than what the model has seen, since we're adding many zeros.
      # As a result, the model prediction of ratioSD will not be accurate
      footprintScores <- -log10(footprintPvals)
      footprintScores <- caTools::runmax(footprintScores, 5)
      footprintScores <- conv(footprintScores, 5) / 10
    },
    footprintRadii,
    mc.cores = 16
  )
  
  smoothedMultiScale <- aggregateMultiScale
  for(i in 1:(dim(aggregateMultiScale)[2] - 2)){
    smoothedMultiScale[, i] <- rowMaxs(aggregateMultiScale[, i:(i+2)])
  }
  
  footprintColors <- colorRamp2(seq(quantile(smoothedMultiScale, 0.01),1,length.out=9),
                                colors = jdb_palette("brewer_blue"))
  
  colnames(smoothedMultiScale) <- rep("", length(footprintRadii))
  colnames(smoothedMultiScale)[footprintRadii %% 10 == 0] <- 2 * footprintRadii[footprintRadii %% 10 == 0]
  
  positions <- 1:dim(smoothedMultiScale)[1]
  rownames(smoothedMultiScale) <- rep("", dim(smoothedMultiScale)[1])
  rownames(smoothedMultiScale)[positions %% 100 == 0] <- positions[positions %% 100 == 0] - plotRadius
  
  pdf(paste0("../../data/IFNStim/plots/aggregate_mulitscale_", group, ".pdf"), width = 8, height = 6)
  print(Heatmap(t(smoothedMultiScale[(plotRadius-300):(plotRadius+300), rev(1:length(footprintRadii))]),
          use_raster = TRUE,
          col = footprintColors,
          cluster_rows = F,
          show_row_dend = F,
          cluster_columns = F,
          name = paste(group, "Footprint\nscore"),
          border = TRUE,
          column_title = "Position (bp)",
          column_title_side = "bottom",
          column_names_rot = 0))
  dev.off()

}


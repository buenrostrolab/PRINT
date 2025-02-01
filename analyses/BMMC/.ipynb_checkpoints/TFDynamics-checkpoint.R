# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

myPath <- .libPaths()
myPath <- c(myPath,'/packages')
.libPaths(myPath)

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getAggregateFootprint.R")
source("../../code/getTFBS.R")
source("../../code/getGroupData.R")
source("../../code/getSubstructures.R")
source("../../code/getGeneCorr.R")

set.seed(42)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(ggpubr)

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "BMMC"
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")

# Load the footprintingProject object
project <- readRDS(paste0(projectDataDir, "footprintingProject.rds"))

# Load motif match positions
motifPositions <- readRDS("../../data/BMMC/motifPositions.rds")

# Get color mapping for cell types
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")

cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")

###########################################
# Find CREs around lineage-specific genes #
###########################################

# Get pseudo-bulks along the erythroid lineage
lineageInd <- 2
lineage <- c("Myeloid", "Erythroid", "BLymphoid", "TLymphoid")[lineageInd]
lineageProbs <- groupInfo[, c("MyeloidProbs", "ErythroidProbs", "BLymphoidProbs", "TLymphoidProbs")]
lineageIdentities <- sapply(1:dim(lineageProbs)[1], 
                            function(x){which(lineageProbs[x,] == max(lineageProbs[x,]))})
lineageGroups <- which(lineageIdentities %in% lineageInd)

# Re-order pseudo-bulks by pseudo-time
lineagePseudoTime <- groupPseudoTime(project)[lineageGroups]
lineageGroups <- lineageGroups[order(lineagePseudoTime)]

# Get pseudo-bulked RNA ordered by pseudotime
RNAMtx <- groupRNA(project)[, lineageGroups]

# Get pseudo-bulked ATAC ordered by pseudotime
ATACMtx <- groupATAC(project)[, lineageGroups]

# Find genes that are strongly activated during erythroid differentiation
ptimeCor <- cor(t(RNAMtx), 1:length(lineagePseudoTime))
ptimeCor[is.na(ptimeCor)] <- 0
lineageGenes <- rownames(RNAMtx)[ptimeCor > 0.5]

# Find CREs within 50kb from these genes
TSS <- FigR::hg38TSSRanges
TSSFilt <- TSS[as.character(TSS$gene_name) %in% lineageGenes]
CREs <- regionRanges(project)
lineageCREInds <- findOverlaps(CREs, resize(TSSFilt, 100000, fix = "center"))@from

# Save CREs passing filter
keptCREs  <- CREs[lineageCREInds]
eryCREBed <- data.frame(
  chr = as.character(seqnames(keptCREs)),
  start = start(keptCREs) - 1,
  end = end(keptCREs)
)
write.table(eryCREBed, "../../data/BMMC/eryCREs.bed", sep = "\t", row.names = F, col.names = F, quote = F)

# Load TF binding scores
TFBindingSE <- getTFBindingSE(project, nCores = 16, regionList = lineageCREInds)
TFBindingSE <- TFBindingSE[, lineageGroups]

# Filter out low signal sites
signal <- rowMaxs(assay(TFBindingSE))
sigFilter <- signal > 0.25

# Only keep sites up-regulated during differentiation
corFilter <- cor(t(assay(TFBindingSE)), 1:length(lineageGroups))[,1] > 0.5

# Remove sites overlapping with more than one CREs
CREOv <- findOverlaps(rowRanges(TFBindingSE), CREs)
nOvCRE <- as.integer(table(factor(CREOv@from, levels = 1:length(signal))))
CREOvFilter <- nOvCRE == 1

# Combine and apply the above filters
filterInds <- which(sigFilter & corFilter & CREOvFilter)
TFBindingSEFilt <- TFBindingSE[filterInds, ]
TFBSScores <- assay(TFBindingSEFilt)
TFBSRanges <- rowRanges(TFBindingSEFilt)

###########################################
# Track binding dynamics of different TFs #
###########################################

# Function to calculate standard error of mean
stderror <- function(x) sd(x)/sqrt(length(x))

smoothRadius <- 5
TFBSDynamics <- pbmcapply::pbmclapply(
  names(cisBPMotifs),
  function(TF){
    
    # Calculate TF binding signal tracks
    TFBSMotifOv <- findOverlaps(TFBSRanges, motifPositions[[TF]])
    if(length(TFBSMotifOv) < 100){return(NULL)}
    TFBSKept <- TFBSRanges[TFBSMotifOv@from]
    TFTracks <- TFBSScores[TFBSMotifOv@from, ]
    TFTracks <- t(sapply(1:nrow(TFTracks), function(i){processTrack(TFTracks[i,])}))
    TFStde <- t(sapply(1:ncol(TFTracks), function(i){stderror(TFTracks[,i])}))[1,]
    
    # Calculate CRE ATAC signal tracks in the same regions
    TFBSCREOv <- findOverlaps(TFBSKept, CREs)
    CRETracks <- ATACMtx[TFBSCREOv@to,]
    CRETracks <- t(sapply(1:nrow(CRETracks), function(i){processTrack(CRETracks[i,])}))
    
    # Calculate RNA tracks of nearby genes
    TFBSPromoterOv <- findOverlaps(TFBSKept, resize(TSSFilt, 100000, fix = "center"))
    RNATracks <- RNAMtx[as.character(TSSFilt$gene_name)[TFBSPromoterOv@to],]
    RNATracks <- t(sapply(1:nrow(RNATracks), function(i){processTrack(RNATracks[i,])}))
    
    # After smoothing and rescaling the curve, we can use the area under the curve (AUC) to
    # quantify the timing of activation. Since the curve goes from 0 to 1, a higher AUC means 
    # earlier activation. If the curve rises very late then the AUC will be low.
    
    # We take the difference in AUC between each TF track and their corresponding CRE (or RNA) track 
    # and use this value as the lag between the two tracks. Taking the mean lag across cases is in fact
    # mathematically the same as taking the mean across all TF tracks and then subtracting the mean of CRE
    # (or RNA) tracks.
    
    # Note that this AUC method mostly works when your curve all go up in pseudotime
    
    # E.g, mean(TFAUC - CREAUC) = mean(rowMeans(TFTracks) - rowMeans(CRETracks)) = mean(TFTracks) - mean(CRETrackss)
    
    # Comput lag between TF and CRE
    CRELag <- mean(TFTracks) - mean(CRETracks)
    
    # Comput lag between TF and RNA level
    RNALag <- mean(TFTracks) - mean(RNATracks)
    
    # Compute the distance of the TF to CRE center
    TFBSCenter <- start(resize(TFBSKept, 1, fix = "center"))
    CRECenter <- start(resize(CREs[TFBSCREOv@to], 1, fix = "center"))
    centerDist <- mean(abs(TFBSCenter - CRECenter))
    
    # Combine results
    list(
      TFTracks = data.frame(
        TFScore = colMeans(TFTracks),
        ptime = 1:ncol(TFTracks),
        stderror = TFStde,
        TF = TF
      ),
      CRELag = CRELag,
      RNALag = RNALag,
      centerDist = centerDist
    )
  },
  mc.cores = 16
)

# Remove invalid results (number of sites < 3)
names(TFBSDynamics) <- names(cisBPMotifs)
TFBSDynamics <- TFBSDynamics[!sapply(TFBSDynamics, is.null)]

# Plot aggregate dynamics tracks for specified TFs
TFBSTracks <- data.table::rbindlist(lapply(TFBSDynamics, function(x){x$TFTracks}))

# Also plot aggregate accessibility
CRETracks <- ATACMtx[findOverlaps(TFBSRanges, CREs)@to,]
CRETracks <- t(sapply(1:nrow(CRETracks), function(i){processTrack(CRETracks[i,])}))
CREStde <- t(sapply(1:ncol(CRETracks), function(i){stderror(CRETracks[,i])}))[1,]
TFBSTracks <- rbind(TFBSTracks, 
                    data.frame(
                      TFScore = colMeans(CRETracks), 
                      ptime = 1:length(CREStde),
                      stderror = CREStde,
                      TF = "Accessibility")
)

system("mkdir ../../data/BMMC/plots/TFDynamics/")

plotTFs <- c("GATA1", "GATA2", "TAL1", "NFE2L2", "RXRA", "NR2F1", "Accessibility")
TFBSTracks <- TFBSTracks[TFBSTracks$TF %in% plotTFs,]
pdf("../../data/BMMC/plots/TFDynamics/TFDynamicLinePlot.pdf", width = 4.5, height = 3.5)
ggplot(data = as.data.frame(TFBSTracks)) +
  geom_line(aes(x = ptime, y = TFScore, group = TF, color = TF)) +
  geom_ribbon(aes(x = ptime, ymin = TFScore - stderror, ymax = TFScore + stderror,
                  group = TF), fill = "grey", alpha = 0.3) +
  ylab("Rescaled TF binding score") + xlab("Pseudo-time") +
  theme_classic()
dev.off()

# Compare time lag of each TF with their average distance to CRE centers 
CRELag <- sapply(TFBSDynamics, function(x){x$CRELag})
RNALag <- sapply(TFBSDynamics, function(x){x$RNALag})
centerDist <- sapply(TFBSDynamics, function(x){x$centerDist})
plotData <- data.frame(
  CRELag = CRELag,
  RNALag = RNALag,
  centerDist = centerDist,
  TF = names(TFBSDynamics)
)
labels <- rep("", length(TFBSDynamics))
labelTFs <- c("GATA1", "GATA2", "TAL1", "NFE2L2", "RXRA", "NFE2", "E2F4", "IRF2", "CUX1",
              "JUN", "TCF3", "SMAD9", "ESR2", "KLF1", "RELA", "SP1", "USF2")
labelTFs <- intersect(labelTFs, plotData$TF)
labels[match(labelTFs, plotData$TF)] <- labelTFs
pdf("../../data/BMMC/plots/TFDynamics/TFLag.pdf", width = 6, height = 5)
ggplot(plotData) +
  geom_point(aes(x = CRELag, y = centerDist, color = RNALag)) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  xlab("TF-vs-CRE accessibilty lag") +
  ylab("Average distance of TF\nbinding site to CRE center (bp)") +
  ggrepel::geom_text_repel(x = CRELag, y = centerDist, label = labels,
                           max.overlaps = 1000) +
  theme_classic()
dev.off()

sort(ptimeCor[intersect(names(cisBPMotifs), rownames(ptimeCor)),], decreasing = T)[1:100]
ptimeCor["EGR1",]
plotData[(plotData$centerDist < 70) & (plotData$CRELag > 0.08),]
plotData[(plotData$CRELag > 0.15),]

ptimeCor["GATA1",]

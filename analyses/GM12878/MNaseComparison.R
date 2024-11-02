# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(rtracklayer)
source("../../code/utils.R")

###################
# Load MNase data #
###################

# Load GM12878 MNase bigwig. The dataset is hg19
# https://www.encodeproject.org/experiments/ENCSR000CXP/
MNaseBigWig <- import.bw("../../data/GM12878/MNase/ENCFF000VME.bigWig")

# Lift-over from hg19 
pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)
seqlevelsStyle(MNaseBigWig) <- "UCSC"  # necessary
MNaseBigWig <- unlist(rtracklayer::liftOver(MNaseBigWig, ch))

# Load genomic regions of CRE
regions <- readRDS("../../data/GM12878/regionRanges.rds")

# Only keep data within CRE regions
MNaseBigWigFilt <- subsetByOverlaps(MNaseBigWig, regions)

##################################################
# Compute single bp resolution MNase data tracks #
##################################################

# Get 1bp bins of the CRE regions
regions$regionInd <- 1:length(regions)
region1bpDataFrame <- data.table::rbindlist(
  pbmcapply::pbmclapply(
    1:length(regions),
    function(regionInd){
      region <- regions[regionInd]
      data.frame(
        chr = as.character(seqnames(region)),
        start = start(region):end(region),
        pos = 1:width(region),
        regionInd = regionInd
      )
    },
    mc.cores = 16
  )
)

# Convert 1bp bins to GRanges object
region1bpPos <- GRanges(
  seqnames = region1bpDataFrame$chr,
  ranges = IRanges(start = region1bpDataFrame$start,
                   end = region1bpDataFrame$start)
)
region1bpPos$pos <- region1bpDataFrame$pos
region1bpPos$regionInd <- region1bpDataFrame$regionInd

# Generate region-by-position matrix of MNase data
MNaseOv <- findOverlaps(MNaseBigWigFilt, region1bpPos)
MNaseTracks <- Matrix::sparseMatrix(
  i = region1bpPos[MNaseOv@to]$regionInd,
  j = region1bpPos[MNaseOv@to]$pos,
  x = MNaseBigWigFilt[MNaseOv@from]$score,
  dims = c(length(regions), unique(width(regions)))
)

###############################
# Load ATAC footprinting data #
###############################

# Load nucleosome footprints
chunkDir <- "../../data/GM12878/chunkedFootprintResults/50/"
chunkFiles <- gtools::mixedsort(list.files(chunkDir))
nucFootprints <- Reduce(rbind, pbmcapply::pbmclapply(
  chunkFiles,
  function(chunkFile){
    chunkFootprints <- readRDS(paste0(chunkDir, chunkFile))
    chunkFootprints <- t(sapply(
      chunkFootprints,
      function(x){
        x$aggregateScores
      }
    ))
  },
  mc.cores = 16
))

# Load Tn5 insertion counts
ATACTracks <- Reduce(rbind, pbmcapply::pbmclapply(
  chunkFiles,
  function(chunkFile){
    chunkFootprints <- readRDS(paste0(chunkDir, chunkFile))
    chunkFootprints <- t(sapply(
      chunkFootprints,
      function(x){
        x$aggregateATAC
      }
    ))
  },
  mc.cores = 16
))

# Remove 100 bp at each edge to prevent edge effect
nucFootprints <- nucFootprints[, 101:900]
MNaseTracks <- MNaseTracks[, 101:900]

##################################
# Compare ATAC and MNase results #
##################################

regionDepth <- rowSums(ATACTracks)
regionFilter <- which((regionDepth > 1000) & (rowMaxs(MNaseTracks) > 2))
corr <- pbmcapply::pbmcmapply(
  function(regionInd){
    cor(nucFootprints[regionInd, ], MNaseTracks[regionInd, ])
  },
  regionFilter,
  mc.cores = 16
)

pdf(paste0("../../data/GM12878/plots/MNase/correlationHist.pdf"), width = 5, height = 5)
data.frame(corr) %>%
  ggplot(aes(x = corr)) +
  geom_histogram(bins=50, fill = '#9C5D41') +
  xlab("Correlation between MNase signal\nand ATAC nucleosome footprints") + 
  ylab("Number of examples") +
  ggtitle(paste0("Median correlation = ", round(median(corr), 2))) +
  theme_classic()
dev.off()

library(patchwork)
positions <- -399:400
regionDepth <- rowSums(ATACTracks)
regionInd <- sample(regionFilter[which(abs(corr - 0.5) < 0.01)], 1)
regionCor <- cor(nucFootprints[regionInd, ], MNaseTracks[regionInd, ], method = "pearson")

# Plot Tn5 insertion counts
plotData <- data.frame(x1 = positions - 2, 
                       x2 = positions + 2, 
                       baseline = rep(0, 800), 
                       Tn5Insertion = ATACTracks[regionInd, 101:900],
                       ATACFootprints = nucFootprints[regionInd, ],
                       MNase = MNaseTracks[regionInd, ])
p1 <- ggplot(plotData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = baseline, ymax = Tn5Insertion), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Tn5\ninsertion") +
  theme_classic() + 
  ggtitle(paste(resize(regions[regionInd], 800, fix = "center"),
                "\ncorrelation = ", round(regionCor, 2)))

p2 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "x1", ymin = 0, ymax = "ATACFootprints"), 
              fill = "#69ACD5") + xlab("") + ylab("ATAC\nnucleosme\nfootprints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

p3 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "x1", ymin = 0, ymax = "MNase"), 
              fill = "#69ACD5") + xlab("") + ylab("MNase\nnucleosome\nfootprints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

pdf(paste0("../../data/GM12878/plots/MNase/region_", regionInd, ".pdf"), width = 5, height = 5)
p1 / p2 / p3 
dev.off()

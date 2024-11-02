# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(DECIPHER)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(png)
library(stats)
source("../../code/utils.R")
source("../../code/getFootprints.R")

#############################################
# Simulate digestion and get fragment sizes #
#############################################

# Load restriction enzyme data
data(RESTRICTION_ENZYMES)
site <- RESTRICTION_ENZYMES[c("SbfI")]

# Get human chromosome DNA sequences
chrList <- paste0("chr", c(1:22, "X", "Y"))
gDNA <- getSeq(BSgenome.Hsapiens.UCSC.hg38, chrList)
names(gDNA) <- chrList

# Simulate digestion and get digested sequences
digestion <- unlist(DigestDNA(site, gDNA))
widths <- width(digestion)
widthsFilt <- widths[widths < 1e5]
hist(log10(widthsFilt), xlab = "log10(fragment size)", ylab = "Number of fragments",
     main = "") 

# Generate the genomic ranges of the digested fragments
cutSites <- unlist(DigestDNA(site, gDNA, type = "positions"))

# For each enzyme cutting site, we will get a top and a bottom coordinate that are off by 3 bp
# To determine the fragments let's only keep the top cut site
cutSites <- cutSites[stringr::str_detect(names(cutSites), "top")]
nCuts <- length(cutSites)
starts <- cutSites[1:(nCuts - 1)]
ends <- cutSites[2:nCuts]
filter <- starts < ends 
fragRanges <- GRanges(seqnames = stringr::str_split_fixed(names(cutSites), "\\.", 2)[1:(nCuts - 1),][filter,1],
                      ranges = IRanges(start = starts[filter], end = ends[filter]))
rtracklayer::export.bed(fragRanges, "../../data/humanGenomicDNA/predictedFrags.bed")

fragSizes <- seq(1, 4.5, 0.2)
CGContent <- sapply(
  fragSizes,
  function(binStart){
    sizeBin <- digestion[(widths > 10 ^ binStart) & (widths < 10 ^ (binStart + 0.2))]
    CGContent <- (vcountPattern("C", sizeBin[1]) + vcountPattern("G", sizeBin[1])) / width(sizeBin[1])
  }
)
qplot(10 ^ fragSizes, CGContent) + scale_x_continuous(trans = "log10") +
  xlab("log10(fragment size)") + ylab("CG content")

#############################################
# Simulate digestion and get fragment sizes #
#############################################

ladderImg <- readPNG("../../data/humanGenomicDNA/2logladder.png")
width <- dim(ladderImg)[2]

findSummits <- function(x){
  sapply
}

ladderIntensity <- conv(ladderImg[,as.integer(width / 2),1], 5)/10
ladderPosition <- 1:dim(ladderImg)[1]
ggplotly(plt(ladderIntensity))
bandPosition <- findSummits(ladderIntensity, threshold = 0.25, r = 15)
bandSize <- c(1e4, 8e3, 6e3, 5e3, 4e3, 3e3, 2e3, 1500, 1200, 1000, 
              900, 800, 700, 600, 500, 400, 300, 200, 100)
bandInfo <- data.frame(
  position = ladderPosition,
  intensity = ladderIntensity
)

# Fit model that predicts log10 frag size using position
reg <- lm(log10(bandSize) ~ bandPosition)

nSlices <- 100
sliceThickness <- round(1000 / nSlices)
results <- as.data.frame(t(sapply(
  1:nSlices,
  function(sliceInd){
    slicePosition <- data.frame(bandPosition = sliceThickness * (sliceInd + 0:1))
    sliceFragSizes <- 10 ^ predict(reg, slicePosition)
    covered <- sum(widths[(widths < sliceFragSizes[1]) & (widths > sliceFragSizes[2])]) / 2
    Reduce(c, list(slicePosition$bandPosition, round(sliceFragSizes), covered / 1e6))
  }
)))

colnames(results) <- c("startPos", "endPos", "startSize", "endSize", "coverage")

ggplot() +
  geom_line(data = bandInfo, mapping = aes(x = position, y = intensity)) +
  xlab("Position") + ylab("Intensity") +
  ggrepel::geom_text_repel(data = bandInfo[bandInfo$position %in% bandPosition,],
                           aes(x = position, y = intensity), label = bandSize)

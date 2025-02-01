# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(GenomicRanges)

# Data is downloaded from https://zenodo.org/record/3905306#.Y5iwNezMI8M
# For details, see https://www.vierstra.org/resources/dgf

# Load footprint-by-sample footprint score matrix
footprintMat <- data.table::fread(
  "../../data/DNaseIFootprint/vierstra/consensus_index_matrix_full_hg38.txt.gz",
  showProgress = T,
  nThread = 16
)[, -1] # First column is not scores. Remove it
footprintMat <- as.matrix(footprintMat)

# Load sample metadata 
sampleMeta <- read.csv(
  "../../data/DNaseIFootprint/vierstra/sampleMeta.tsv",
  header = T, sep = "\t")
stopifnot(dim(sampleMeta)[1] == dim(footprintMat)[2])

# Samples for K562: K562-DS15363, K562-DS16924, h.K562-DS52908
dataset <- "K562-DS15363"
footprintScores <- footprintMat[, match(dataset, sampleMeta$Identifier)]

# Load row metadata for the above matrix (genomic ranges of footprints)
metaData <- data.table::fread(
  "../../data/DNaseIFootprint/vierstra/consensus_footprints_and_motifs_hg38.bed.gz",
  showProgress = T,
  nThread = 16
)

# Convert to GRanges object
footprintRanges <- GRanges(seqnames = metaData$V1,
                           ranges = IRanges(start = metaData$V2,
                                            end = metaData$V3))

# Load TF motif matches across the genome
TFMotifRanges <- readRDS("../../data/BMMC/TFMotifRanges.rds")

# Retrieve footprints for each TF based on overlap with motifs
DNaseIFootprints <- pbmcapply::pbmclapply(
  names(TFMotifRanges),
  function(TF){
    ov <- findOverlaps(footprintRanges, TFMotifRanges[[TF]])
    TFFootprintRanges <- footprintRanges[ov@from]
    TFFootprintRanges$score <- footprintScores[ov@from]
    TFFootprintRanges
  },
  mc.cores = 16
)
names(DNaseIFootprints) <- names(TFMotifRanges)

# Save results
saveRDS(DNaseIFootprints, paste0("../../data/DNaseIFootprint/vierstra/", dataset, "FootprintRanges.rds"))

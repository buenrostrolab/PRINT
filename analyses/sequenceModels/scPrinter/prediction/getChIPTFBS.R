# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(ggplot2)
library(GenomicRanges)

##################################
# Determine which dataset to use #
##################################

dataset <- "GM12878"

#################
# Get ChIP data #
#################

# Get the list of TFs with ChIP data
ChIPFolder <- paste0("/data/PRINT/multiScaleFootprinting/data/", dataset, "/ENCODEChIP/")
ChIPMetadata <- read.csv(paste0("/data/PRINT/multiScaleFootprinting/data/",  dataset, "/metadata.tsv"), sep = "\t")
ChIPMetadata <- ChIPMetadata[ChIPMetadata$File.format == "bed narrowPeak",]
ChIPMetadata <- ChIPMetadata[ChIPMetadata$File.assembly == "GRCh38",]
ChIPTFs <- unname(sapply(
  ChIPMetadata$Experiment.target,
  function(target){
    strsplit(target, "-")[[1]][1]
  }
))

# Retrieve the genomic ranges of the TFs
TFChIPRanges <- pbmcapply::pbmclapply(
  sort(unique(ChIPTFs)),
  function(TF){
    ChIPRangeList <- lapply(
      which(ChIPTFs %in% TF),
      function(entry){
        ChIPBed <- read.csv(paste0(ChIPFolder, ChIPMetadata$File.accession[entry], ".bed.gz"), 
                            sep = "\t", header = F)
        ChIPRanges <- GRanges(
          seqnames = ChIPBed[,1],
          ranges = IRanges(start = ChIPBed[, 2], end = ChIPBed[,3])
        )
        ChIPRanges$signal <- ChIPBed[, 9]
        ChIPRanges
      }
    )
    # For TFs with more than one ChIP files, take the intersection of regions in all files
    ChIPRanges <- Reduce(subsetByOverlaps, ChIPRangeList)
    ChIPRanges
  },
  mc.cores = 16
)
names(TFChIPRanges) <- sort(unique(ChIPTFs))
saveRDS(TFChIPRanges, paste0("/data/PRINT/multiScaleFootprinting/data/", dataset, "/ENCODEChIPRanges.rds"))

# Also convert the integrated list to bed format
TFChIPBed <- lapply(
  sort(unique(ChIPTFs)),
  function(TF){
    data.frame(
      chr = as.character(seqnames(TFChIPRanges[[TF]])),
      start = start(TFChIPRanges[[TF]]),
      end = end(TFChIPRanges[[TF]]),
      TF = TF
    )
  }
)
TFChIPBed <- as.data.frame(data.table::rbindlist(TFChIPBed))
write.table(TFChIPBed, paste0("/data/PRINT/multiScaleFootprinting/data/", dataset, "/ENCODEChIPBed.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = F)

#########################
# Functions we will use #
#########################

# Calculate entropy
entropy <- function(x, nbins = 100){
  freq <- hist(x, breaks = nbins, plot = F)$counts
  prob <- freq / sum(freq)
  entropy <- -sum((prob * log(prob))[prob > 0])
}

# Max entropy thresholding
autothreshold <- function(x, nThresholds = 100){
  thresholds <- seq(min(x), max(x), length.out = nThresholds + 2)
  thresholds <- thresholds[2:(length(thresholds) - 1)]
  entropyList <- sapply(
    thresholds,
    function(threshold){
      entropy(x[x < threshold]) + entropy(x[x >= threshold])
    }
  )
  optimalThreshold <- thresholds[entropyList == max(entropyList)]
  return(optimalThreshold)
}

###################################################
# For each TF, find high confidence binding sites #
###################################################

# Load ENCODE ChIP-seq data
ENCODEChIP <- readRDS(paste0("/data/PRINT/multiScaleFootprinting/data/", dataset, "/ENCODEChIPRanges.rds"))

# Load motif PWMs
cisBPMotifs <- readRDS("/data/PRINT/multiScaleFootprinting/data/shared/cisBP_human_pwms_2021.rds")

# Load ATAC-peak regions
regions <- readRDS(paste0("/data/PRINT/multiScaleFootprinting/data/", dataset, "/regionRanges.rds"))

keptTFs <- intersect(names(ENCODEChIP), names(cisBPMotifs))
TFBS <- data.table::rbindlist(pbmcapply::pbmclapply(
  keptTFs,
  function(TF){
    
    # Retreive ChIP ranges
    ChIPRanges <- ENCODEChIP[[TF]]
    ChIPCenters <- resize(ChIPRanges, 1, fix = "center")
    
    # Get motif match positions
    motifSites <- motifmatchr::matchMotifs(
      pwms = cisBPMotifs[TF],
      subject = regions,
      out = "position",
      genome = "hg38",
      p.cutoff = 5e-5
    )[[1]]
    motifCenters <- resize(motifSites, 1, fix = "center")
    
    # For each motif site, find the overlapping ChIP-site
    # Calculate relative position of motif to ChIP peak center
    ov <- findOverlaps(motifCenters, ChIPRanges)
    if(length(ov) <= 200){return(NULL)}
    relativePos <- start(motifCenters[ov@from]) - start(ChIPCenters[ov@to])
    relativeDistance <- abs(relativePos)
    
    # Detect optimal threshold for the distance between motif and ChIP-peak center
    selectedSites <- motifSites[ov@from][relativeDistance <= autothreshold(relativeDistance)]
    
    if(length(selectedSites) > 0){
      data.frame(
        chr = as.character(seqnames(selectedSites)),
        start = start(selectedSites) - 1,
        end = end(selectedSites),
        strand = strand(selectedSites),
        TF = TF
      )
    }else{
      NULL
    }
  },
  mc.cores = 8
))

write.table(TFBS, paste0("/data/PRINT/multiScaleFootprinting/data/", dataset, "/filteredTFBS.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = F)

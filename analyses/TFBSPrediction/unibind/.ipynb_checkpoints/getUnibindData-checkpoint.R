# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
myPath <- .libPaths()
myPath <- c(myPath,'/packages')
.libPaths(myPath)
library(GenomicRanges)
source("../../../code/utils.R")

dataset <- "K562"

# Get dataset metadata
ChIPDir <- "../../../data/shared/unibind/damo_hg38_all_TFBS/"
ChIPSubDirs <- list.files(ChIPDir)
ChIPMeta <- as.data.frame(t(sapply(
  ChIPSubDirs, 
  function(ChIPSubDir){
    strsplit(ChIPSubDir, "\\.")[[1]]
  }
)))
rownames(ChIPMeta) <- NULL
colnames(ChIPMeta) <- c("ID", "dataset", "TF")
ChIPMeta$Dir <- ChIPSubDirs

# Select cell types
datasetName <- list("K562" = "K562_myelogenous_leukemia",
                    "HepG2" = "HepG2_hepatoblastoma",
                    "GM12878" = "GM12878_female_B-cells_lymphoblastoid_cell_line",
                    "A549" = "A549_lung_carcinoma")
ChIPMeta <- ChIPMeta[ChIPMeta$dataset == datasetName[[dataset]], ]
ChIPTFs <- ChIPMeta$TF

# If we have multiple files for the same TF, decide whether to keep the intersection or union of sites
mode <- "union"

# Extract TFBS ChIP ranges
unibindTFBS <- pbmcapply::pbmclapply(
  sort(unique(ChIPTFs)),
  function(TF){
    
    # Retrieve the list of ChIP Ranges
    ChIPRangeList <- lapply(
      which(ChIPTFs %in% TF),
      function(entry){
        ChIPSubDir <- ChIPMeta$Dir[entry]
        ChIPFiles <- list.files(paste0(ChIPDir, ChIPSubDir))
        ChIPRanges <- lapply(
          ChIPFiles,
          function(ChIPFile){
            ChIPBed <- read.table(paste0(ChIPDir, ChIPSubDir, "/", ChIPFile))
            ChIPRange <- GRanges(seqnames = ChIPBed$V1, ranges = IRanges(start = ChIPBed$V2, end = ChIPBed$V3))
            ChIPRange$score <- ChIPBed$V5
            ChIPRange
          }
        )
        if(mode == "intersect"){
          ChIPRanges <- Reduce(subsetByOverlaps, ChIPRanges)
        }else if (mode == "union"){
          ChIPRanges <- mergeRegions(Reduce(c, ChIPRanges))
        }
        
        ChIPRanges
      }
    )
    
    # For TFs with more than one ChIP files, take the intersection or union of regions in all files
    if(mode == "intersect"){
      ChIPRanges <- Reduce(subsetByOverlaps, ChIPRangeList)
    }else if (mode == "union"){
      ChIPRanges <- mergeRegions(Reduce(c, ChIPRangeList))
    }
    ChIPRanges
  },
  mc.cores = 12
)
names(unibindTFBS) <- sort(unique(ChIPTFs))

# Save results to file
if(mode == "intersect"){
  saveRDS(unibindTFBS, paste0("../../../data/shared/unibind/", dataset, "ChIPRanges.rds"))
}else if (mode == "union"){
  saveRDS(unibindTFBS, paste0("../../../data/shared/unibind/", dataset, "UnionChIPRanges.rds"))
}

# Also save to bed file
unibindBed <- data.table::rbindlist(
  pbmcapply::pbmclapply(
    names(unibindTFBS),
    function(TF){
      sites <- unibindTFBS[[TF]]
      data.frame(
        chr = as.character(seqnames(sites)),
        start = start(sites),
        end = end(sites),
        TF = TF
      )
    },
    mc.cores = 12
  )
)
if(mode == "intersect"){
  path <- paste0("../../../data/shared/unibind/", dataset, "ChIPRanges.bed")
}else if (mode == "union"){
  path <- paste0("../../../data/shared/unibind/", dataset, "UnionChIPRanges.bed")
}
write.table(
  unibindBed, path,
  sep = "\t", quote = F, col.names = F, row.names = F
)


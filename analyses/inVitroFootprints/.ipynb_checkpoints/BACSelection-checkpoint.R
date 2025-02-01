# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

TFChIPDir <- "../../data/inVitroFootprints/ChIP/"
TFChIPFiles <- list.files(TFChIPDir)
ChIPTFs <- stringr::str_split_fixed(TFChIPFiles, "\\.", 7)[, 3]
TFs <- unique(ChIPTFs)
TFChIPRanges <- lapply(
  TFs,
  function(TF){
    mergeRegions(Reduce(c,
      lapply(
        which(ChIPTFs ==TF), 
        function(x){
          rtracklayer::import.bed(paste0(TFChIPDir, TFChIPFiles[x]))
        })
    ))
  }
)
names(TFChIPRanges) <- TFs

# Load genomic ranges of BACs
BACInfo <- read.csv("/n/home09/yanhu/cell_dynamics/SHARE_footprinting/data/BAC/RP11.GCF_000001405.38.118.unique_concordant.gff",
                    sep = "\t", comment.char = "#", header = F)

# Extract info of each clone
cloneIDs <- unname(sapply(BACInfo$V9, function(x){stringr::str_match(x, "Name=(.+?);")[2]}))
chr <- unname(sapply(BACInfo$V1, function(x){
  chrInd <- as.integer(stringr::str_match(x, "NC_(.+?)\\.")[2])
  if(is.na(chrInd)){chrInd <- "NA"}
  if(chrInd == 23){chrInd <- "X"}
  if(chrInd == 24){chrInd <- "Y"}
  paste("chr", chrInd, sep = "")
}))

# Convert to genomic ranges
BACRanges <- GRanges(
  seqnames = chr,
  ranges = IRanges(start = BACInfo$V4, end = BACInfo$V5)
)
names(BACRanges) <- cloneIDs

# Remove entries that are not in standard chromosomes
BACRanges <- BACRanges[!as.character(seqnames(BACRanges)) == "chrNA"]

# Filter out BACs with abnormal sizes
BACRanges <- BACRanges[width(BACRanges) < 5e5]

# Count the number of ChIP peaks in each BAC region
results <- data.table::rbindlist(lapply(
  TFs,
  function(TF){
    BACInds <- sort(table(findOverlaps(BACRanges, TFChIPRanges[[TF]])@from), decreasing = T)[1:10]
    BACIDs <- names(BACRanges)[as.integer(names(BACInds))]
    data.frame(
      BACIDs = BACIDs,
      BACRange = as.character(BACRanges[as.integer(names(BACInds))]),
      numOverlap = unname(as.integer(BACInds)),
      TF = TF)
  }
))

write.table(results, "../../data/inVitroFootprints/selectedBACs.txt", 
            row.names = F, col.names = T, quote = F, sep = "\t")

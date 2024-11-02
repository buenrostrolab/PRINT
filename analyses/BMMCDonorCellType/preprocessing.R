# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(SummarizedExperiment)

#################################
# 1. Get genomic ranges of CREs #
#################################

# Load single cell ATAC
scATAC <- readRDS("../../data/BMMCCellType/atac.se.rds")

# Get region ranges
regions <- rowRanges(scATAC)
regions <- IRanges::resize(regions, width = 1000, fix = "center")
saveRDS(regions, "../../data/BMMCCellType/regionRanges.rds")

############################
# 2. Generate pseudo-bulks #
############################

# Reformat cell barcodes
scBarcodes <- colnames(scATAC)

# Get cell type labels
cellTypeLabels <- scATAC$cistopic.assign.l2.rank
cellTypes <- sort(unique(cellTypeLabels))

# Get donor labels
donorLabels <- scATAC$donor
donors <- sort(unique(donorLabels))

# Pseudobulk each (donor, cell type) combination
barcodeGroups <- list()
groupInfo <- list()
for(donor in donors){
  for(cellType in cellTypes){
    groupID <- paste0(donor, "_",  cellType)
    barcodeGroups[[groupID]] <-
      data.frame(
        barcode = scBarcodes[(donorLabels == donor) & (cellTypeLabels == cellType)],
        group = groupID
      )
    groupInfo[[groupID]] <-
      data.frame(
        group = groupID,
        donor = donor,
        cellType = cellType
      )
  }
}
barcodeGroups <- Reduce(rbind, barcodeGroups)
system("mkdir ../../data/BMMCDonorCellType")
write.table(barcodeGroups, "../../data/BMMCDonorCellType/barcodeGrouping.txt",
            row.names = F, sep = "\t", quote = F)

# Also save metadata for each pseudobulk
groupInfo <- Reduce(rbind, groupInfo)
write.table(groupInfo, "../../data/BMMCDonorCellType/groupInfo.txt",
            row.names = F, sep = "\t", quote = F)

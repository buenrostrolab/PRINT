# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

myPath <- .libPaths()
myPath <- c(myPath,'/packages')
.libPaths(myPath)

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getGroupData.R")

###############################
# Load pseudobulk information #
###############################

# Group pseudo-bulk clusters into sub-populations
pbulkClusters <- read.table("../../data/mHSCAging10xMultiome/pbulkClusters.txt")$V1
subpopAnno <- list(
  "Old_1" = "Old Mk-biased", 
  "Old_2" = "Old intermediate",
  "Old_3" = "Old Mk-biased",
  "Old_4" = "Old multi-lineage",
  "Young_1" = "Young multi-lineage",
  "Young_2" = "Young multi-lineage",
  "Young_3" = "Young Mk-biased")
subpopLabels <- unname(sapply(pbulkClusters, function(x){subpopAnno[[x]]}))

####################################################
# Compute peak-by-pseudobulk matrix of ATAC counts #
####################################################

# Initialize a footprintingProject object
projectName <- "mHSCAging10xMultiome/"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "mm10")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName)
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regions <- readRDS(paste0(projectMainDir, "data/mHSCAging10xMultiome/regionRanges.rds"))

# Set the regionRanges slot
regionRanges(project) <- regions

# Load barcodes for each pseudo-bulk
pathToFragGrouping <- paste0(projectDataDir, "barcodeGrouping.txt")
barcodeGroups <- read.table(pathToFragGrouping, header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- gtools::mixedsort(unique(barcodeGroups$group))

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- "../../data/mHSCAging10xMultiome/all.frags.filt.tsv.gz"
if(!dir.exists("../../data/mHSCAging10xMultiome/chunkedCountTensor/")){
  getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = F)
}

# Get region-by-pseudobulk ATAC matrix 
pseudobulkATACPath <- paste0(projectDataDir, "pseudobulkATAC.rds")
if(!file.exists(pseudobulkATACPath)){
  pbulkCounts <- getGroupATAC(project)
  rownames(pbulkCounts) <- as.character(regions)
  saveRDS(pbulkCounts, pseudobulkATACPath)
}else{
  pbulkCounts <- readRDS(pseudobulkATACPath)
}

###############################################################################
# Calculate differential ATAC peaks between young and old for each cell state #
###############################################################################

for(subpop in c("Mk-biased", "multi-lineage", "all")){
  
  # Get metadata for each pseudobulk
  metadata <- as.data.frame(stringr::str_split_fixed(subpopLabels, " ", 2))
  colnames(metadata) <- c("age", "subpop")
  
  # Generate DDS object
  if(subpop == "all"){
      filter <- 1:dim(pbulkCounts)[2]
  }else{
      filter <- metadata$subpop == subpop
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pbulkCounts[, filter],
    colData = metadata[filter, ],
    design = ~ age)
  
  # Model fitting
  dds <- DESeq2::DESeq(dds)
  
  # Retrieve differential analysis results
  diffATAC <- DESeq2::results(dds, contrast = c("age", "Old", "Young"))
  diffATAC <- diffATAC[!is.na(diffATAC$padj),]
  diffATAC <- diffATAC[order(diffATAC$pvalue),] 
  
  write.table(diffATAC, paste0("../../data/mHSCAging10xMultiome/diffATAC_", subpop, ".tsv"),
              sep = "\t", quote = F, row.names = T)
}




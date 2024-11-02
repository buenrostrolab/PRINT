# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getGroupData.R")

###################
# Load input data #
###################

# Load footprinting project
project <- readRDS("../../data/mHSCAging10xMultiome/project.rds")

# Load single cell RNA data
scRNA <- readRDS("../../data/mHSCAging10xMultiome/scRNA.rds")

# Load barcodes for each pseudo-bulk
barcodeGroups <- read.table("../../data/mHSCAging10xMultiome/barcodeGrouping.txt", header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- gtools::mixedsort(unique(barcodeGroups$group))

#############################
# Differential RNA analysis #
#############################

# Get pseudobulk RNA count matrix
RNAMatrix <- scRNA@assays$RNA@data
pseudobulkRNA <- getGroupRNA(RNAMatrix, barcodeGroups)

# Get pseudobulk metadata
metadata <- t(sapply(colnames(pseudobulkRNA), 
                     function(x){
                       strsplit(x, "_")[[1]][c(1,2)]
                     }))
metadata <- as.data.frame(metadata)
colnames(metadata) <- c("Age", "pbulk_ind")

# Generate DDS object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = pseudobulkRNA,
  colData = metadata,
  design = ~ Age)

# Model fitting
dds <- DESeq2::DESeq(dds)

# Retrieve differential analysis results
diffRNA <- DESeq2::results(dds, contrast = c("Age","Old","Young"))
diffRNA <- diffRNA[!is.na(diffRNA$padj),]
diffRNA <- diffRNA[order(diffRNA$pvalue),]
write.table(diffRNA, "../../data/mHSCAging10xMultiome/diffRNA.tsv",
            sep = "\t", quote = F)

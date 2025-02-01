# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(Seurat)
source("../../code/utils.R")

###################
# Load input data #
###################

# Load pathway gene sets from hypeR
library(hypeR)
c5GO <- msigdb_gsets(species = "Homo sapiens","C5","BP",clean = TRUE)$genesets
names(c5GO) <- stringr::str_replace_all(names(c5GO), " ",  "_")

# Load CRISPRa data
scMatrix <- Matrix::readMM("../../data/K562/perturbSeq/CRISPRa/GSE133344_filtered_matrix.mtx.gz")
cellMeta <- read.csv("../../data/K562/perturbSeq/CRISPRa/GSE133344_filtered_cell_identities.csv.gz")
barcodes <- read.table("../../data/K562/perturbSeq/CRISPRa/GSE133344_filtered_barcodes.tsv.gz", sep = "\t")
genes <- read.table("../../data/K562/perturbSeq/CRISPRa/GSE133344_filtered_genes.tsv.gz", sep = "\t")
colnames(scMatrix) <- barcodes$V1
rownames(scMatrix) <- genes$V2
rownames(cellMeta) <- cellMeta$cell_barcode

# Convert dgTMatrix format to dgCMatrix
scMatrix <- as(scMatrix, "dgCMatrix")

# Only keep cells with assigned guide identities
scBarcodes <- intersect(colnames(scMatrix), cellMeta$cell_barcode)
scMatrix <- scMatrix[, scBarcodes]
cellMeta <- cellMeta[scBarcodes, ]

# Only keep cells with good coverage
coverageFilter <- cellMeta$good_coverage == "True"
scMatrix <- scMatrix[, coverageFilter]
cellMeta <- cellMeta[coverageFilter, ]

# Parse guide info
guideInfo <- as.data.frame(stringr::str_split_fixed(cellMeta$guide_identity, "_", 5)[, 1:2])

# Load motif database
cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")

##########################################
# Get non-targeted negative control data #
##########################################

numNegCtrl <- stringr::str_count(cellMeta$guide_identity, pattern = "NegCtrl") / 2
ctrlInds <- which(numNegCtrl == 2)

# Pseudobulk background cells (cells with non-targeting guides)
nPbulks <- 10
set.seed(42)
pbulkInds <- sample(rep(1:nPbulks, length.out = length(ctrlInds)))
ctrlPbulkCounts <- pbmcapply::pbmcmapply(
  function(i){
    rowSums(scMatrix[, ctrlInds[pbulkInds == i]])
  },
  1:nPbulks,
  mc.cores = 16
)

########################################
# Find pathways enriched in TF targets #
########################################

selectedTFs <- intersect(c(guideInfo$V1, guideInfo$V2), names(cisBPMotifs))

# Run DESeq to compare control-vs-TF perturbation
perturbDiffRNA <- pbmcapply::pbmclapply(
  selectedTFs,
  function(guideTarget){
    
    # Find cells targeted by guide for this gene (1 targeting guide and 1 control guide)
    fgInds <- which(((guideInfo$V1 == guideTarget) | (guideInfo$V2 == guideTarget)) & (numNegCtrl == 1))
    
    # Pseudobulk cells 
    nPbulks <- 5
    set.seed(42)
    pbulkInds <- sample(rep(1:nPbulks, length.out = length(fgInds)))
    fgPbulkCounts <- sapply(
      1:nPbulks,
      function(i){
        rowSums(scMatrix[, fgInds[pbulkInds == i]])
      }
    )
    
    # Get pseudobulk metadata
    metadata <- data.frame(group = c(rep("fg", dim(fgPbulkCounts)[2]),
                                     rep("bg", dim(ctrlPbulkCounts)[2])))
    
    # Generate DDS object
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = cbind(fgPbulkCounts, ctrlPbulkCounts),
      colData = metadata,
      design = ~ group)
    
    # Model fitting
    dds <- DESeq2::DESeq(dds)
    
    # Retrieve differential analysis results
    diffRNA <- DESeq2::results(dds, contrast = c("group","fg","bg"))
    diffRNA <- diffRNA[!is.na(diffRNA$padj),]
    diffRNA <- diffRNA[order(sign(diffRNA$log2FoldChange) * -log10(diffRNA$padj),
                             decreasing = T),]
    
    diffRNA
  },
  mc.cores = 16
)
names(perturbDiffRNA) <- selectedTFs
saveRDS(perturbDiffRNA, "../../data/K562/perturbSeq/CRISPRa/perturbDiffRNA.rds")

# For each TF, get a list of target genes
TFTargetGenes <- lapply(
  selectedTFs,
  function(TF){
    data <- perturbDiffRNA[[TF]]
    rownames(data[(data$padj < 1e-3) & (data$log2FoldChange > 1), ])
  }
)
names(TFTargetGenes) <- selectedTFs

saveRDS(TFTargetGenes, "../../data/K562/perturbSeq/CRISPRa/TFTargetGenes.rds")

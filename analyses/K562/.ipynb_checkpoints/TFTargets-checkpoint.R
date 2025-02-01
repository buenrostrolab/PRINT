# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

myPath <- .libPaths()
myPath <- c(myPath,'/packages')
.libPaths(myPath)

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getTFBS.R")
source("../../code/getGroupData.R")
source("../../code/getTargetPathway.R")
library(ComplexHeatmap)
library(SummarizedExperiment)

###################
# Load input data #
###################

# Load footprinting project object
project <- readRDS("../../data/K562/project.rds")

# Load ATAC-seq peak regions
regions <- readRDS("/data/PRINT/multiScaleFootprinting/data/K562/regionRanges.rds")

# Load TFBS scores
bw <- rtracklayer::import.bw("/data/PRINT/multiScaleFootprinting/data/K562/K562_TFBS.bigwig")

# Load CRE-by-pseudobulk accessibility matrix
path <- "../../data/K562/pseudobulkATAC.rds"
if(!file.exists(path)){
  pseudobulkATAC <- getGroupATAC(project)
  saveRDS(pseudobulkATAC, path)
}else{
  pseudobulkATAC <- readRDS(path)
}
CRESE <- SummarizedExperiment(assay = pseudobulkATAC,
                              rowRanges = regionRanges(project))

# Generate a vector of TF binding scores
TFBindingSE <- SummarizedExperiment(
  assays = list(TFBS = as.matrix(bw$score)),
  rowRanges = bw
)

# Load TF motif matches
TFMotifRanges <- readRDS("../../data/K562/TFMotifRanges.rds")

# Load CRE-gene mapping
CREGeneCorr <- readRDS("../../data/BMMC/CREGeneCorr.rds")
BMMCRegions <- readRDS("../../data/BMMC/regionRanges.rds")
CREGeneCorr$Region <- as.character(BMMCRegions)[CREGeneCorr$CREInd]
CREGeneCorr$fdr <- p.adjust(CREGeneCorr$pvalZ, method = "fdr")
CREGeneCorr <- CREGeneCorr[CREGeneCorr$fdr < 0.2, ]

# Load pathway gene sets from hypeR
library(hypeR)
c5GO <- msigdb_gsets(species = "Homo sapiens","C5","BP",clean = TRUE)$genesets
names(c5GO) <- stringr::str_replace_all(names(c5GO), " ",  "_")

# Create a folder to save plots
system("mkdir -p ../../data/K562/plots/TFTargets/")

##########################################
# Map a single TF to its target pathways #
##########################################

# Decide whether to use PRINT to score TF binding
mode <- "PRINT" # "PRINT" or "motif"
if(mode == "PRINT"){
  siteSE <- TFBindingSE
}else{
  siteSE <- CRESE
}

TF <- "TFEB" #XBP1, EGR1, HLF
enrichment <- mapTargets(
  TFMotifSites = TFMotifRanges[TF],  
  siteSE = siteSE,
  regionGeneCorr = CREGeneCorr,
  geneSets = c5GO,
  threshold = 0.75,
  genome = "hg38")

# Visualize results
plotData <- data.frame(
  Name = enrichment[[TF]]$pathway[1:10],
  logFDR = -log10(enrichment[[TF]]$fdr[1:10])
)
plotData <- plotData[nrow(plotData):1, ]
plotData$Name <- stringr::str_replace_all(plotData$Name, "_", " ")
plotData$Name <- sapply(plotData$Name, function(s){gsub('(.{1,30})(\\s|$)', '\\1\n', s)}) # Add newline to long strings
plotData$Name <- factor(plotData$Name, levels = plotData$Name) # This keeps the entries in the original order when plotting
pdf(paste0("../../data/K562/plots/TFTargets/", TF, "_", mode, "_targets_seq_model.pdf"))
ggplot(plotData) +
  geom_bar(aes(x = Name, y = logFDR), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Top hits")  + 
  ylab("Log10(FDR)") +
  coord_flip() +
  ggtitle(TF) +
  theme_classic() +
  theme(axis.text = element_text(size = 12))  
dev.off()

#############################################################
# Map TFs to target genes and calculate pathway enrichment  #
#############################################################

TFs <- names(TFMotifRanges)
enrichmentList <- mapTargets(
  TFMotifSites = TFMotifRanges,  
  siteSE = siteSE,
  regionGeneCorr = CREGeneCorr,
  geneSets = c5GO,
  threshold = 0.75,
  genome = "hg38")
saveRDS(enrichmentList, "../../data/K562/TFTargetEnrichemnt_seq_model.rds")

# Select top pathways for each TF
selectedPathways <- unique(Reduce(c, lapply(
  TFs,
  function(TF){
    TFFilter <- enrichmentList$TF == TF 
    dt <- TFTargetEnrichment[TFFilter, ]
    dt <- dt[order(dt$fdr),][1:5, ]
    dt$pathway
  }
)))

TFTargetEnrichment <- TFTargetEnrichment[TFTargetEnrichment$pathway %in% selectedPathways, ]

targetPathwayMatrix <- as.matrix(Matrix::sparseMatrix(i = match(TFTargetEnrichment$pathway, selectedPathways),
                                                      j = match(TFTargetEnrichment$TF, TFs),
                                                      x = -log10(TFTargetEnrichment$fdr)))
rownames(targetPathwayMatrix) <- selectedPathways
colnames(targetPathwayMatrix) <- TFs
targetPathwayMatrix <- targetPathwayMatrix[, colMaxs(targetPathwayMatrix) > 5]

pdf("../../data/K562/plots/TFTargets/TFTargetEnrichemnt.pdf", width = 50, height = 100)
Heatmap(pmin(targetPathwayMatrix, 10),
        col = BuenColors::jdb_palette("solar_extra"),
        row_names_gp = grid::gpar(fontsize = 8))
dev.off()

write.table(t(targetPathwayMatrix), "../../data/K562/TFTargetEnrichemntMatrix.tsv",
            quote = F, row.names = T, col.names = T, sep = "\t")

sort(rownames(targetPathwayMatrix))

sort(targetPathwayMatrix[, "TFEB"], decreasing = T)[1:20]

###############################
# Pathway enrichment analysis #
###############################

ID <- "Electron_Transport_Chain"
if(ID %in% names(cisBPMotifs)){
  # If trying to find target pathways of a selected TF
  topHits <- sort(targetPathwayMatrix[, ID], decreasing = T)[10:1]
}else{
  # If trying to find regulators of a selected pathway
  topHits <- sort(targetPathwayMatrix[ID, ], decreasing = T)[10:1]
}

# Visualize results
plotData <- data.frame(
  Name = names(topHits),
  logFDR = unname(topHits)
)
plotData$Name <- stringr::str_replace_all(plotData$Name, "_", " ")
plotData$Name <- sapply(plotData$Name, function(s){gsub('(.{1,25})(\\s|$)', '\\1\n', s)}) # Add newline to long strings
plotData$Name <- factor(plotData$Name, levels = plotData$Name) # This keeps the entries in the original order when plotting
ggplot(plotData) +
  geom_bar(aes(x = Name, y = logFDR), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Top hits")  + 
  ylab("Log10(FDR)") +
  coord_flip() +
  ggtitle(ID) +
  theme_classic() +
  theme(axis.text = element_text(size = 12))  

write.table(enrichment, "test.txt", quote = F, sep = "\t")
write.table(targetGenes, "test.txt", quote = F, row.names = F, col.names = F)

###############################
# Validate using CRISPRa data #
###############################

# Load TF target genes identified by CRISPRa
validatedTargets <- readRDS("../../data/K562/perturbSeq/CRISPRa/TFTargetGenes.rds")

# Keep TFs with > 50 targets
validatedTargets <- validatedTargets[sapply(validatedTargets, length) > 50]

# Map each TF to these validated target lists
TFs <- intersect(names(validatedTargets), names(TFMotifRanges))
targetEnrichment <- mapTargets(
  TFMotifSites = TFMotifRanges[TFs],  
  siteSE = siteSE,
  regionGeneCorr = CREGeneCorr,
  geneSets = validatedTargets,
  threshold = 0.75,
  genome = "hg38",
  proximalOnly = T)

# Reformat into a TF (CRISPR) by TF (footprinted) matrix
enrichmentMatrix <- sapply(
  TFs,
  function(TF){
    data <- targetEnrichment[[TF]][TFs, ]$fdr
    data[is.na(data)] <- 1
    -log10(data)
  }
)
dimnames(enrichmentMatrix) <- list(paste0(TFs, " CRISPRa"), TFs)
while(min(colMaxs(enrichmentMatrix)) < -log10(0.25)){
  filtTFs <- colMaxs(enrichmentMatrix) > -log10(0.25)
  enrichmentMatrix <- enrichmentMatrix[filtTFs, filtTFs]
}

# Visualize on heatmap
pltMtx <- pmin(enrichmentMatrix, 10)
rowOrder <- hclust(dist(t(pltMtx)))$order
pdf("../../data/K562/plots/TFTargets/CRISPRaValidation.pdf")
ComplexHeatmap::Heatmap(pltMtx[rowOrder, rowOrder], 
        name = "z-score",
        border = T,
        col = BuenColors::jdb_palette("solar_rojos"),
        cluster_rows = F, 
        cluster_columns = F)
dev.off()

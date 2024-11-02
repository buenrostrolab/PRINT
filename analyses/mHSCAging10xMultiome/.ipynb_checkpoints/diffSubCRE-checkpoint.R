# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getFootprints.R")
source("../../code/getSubstructures.R")
source("../../code/getTFBS.R")
source("../../code/getTargetPathway.R")
source("../../code/visualization.R")
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

###################
# Load input data #
###################

# Load footprinting project
project <- readRDS("../../data/mHSCAging10xV3/project.rds")
regions <- regionRanges(project)

# Get subCRE-by-pseudobulk SummarizedExperiment object
subCREPath <- "../../data/mHSCAging10xV3/subCRESE.rds"
if(file.exists(subCREPath)){
  subCRESE <- readRDS("../../data/mHSCAging10xV3/subCRESE.rds")
}else{
  subCRESE <- getSubstructureSE(project)
  saveRDS(subCRESE, subCREPath)
}
subCRERanges <- rowRanges(subCRESE)
subCREMat <- assay(subCRESE)

# Load results from differential RNA testing
diffRNA <- read.table("../../data/mHSCAging10xMultiome/diffRNA.tsv")

# Load pseudobulk clustering
pbulkClusters <- read.table("../../data/mHSCAging10xV3/pbulkClusters.txt")
cellStates <- rep("", length(pbulkClusters$V1))
cellStates[pbulkClusters$V1 %in% c("Old_1", "Old_3")] <- "Old Mk-biased"
cellStates[pbulkClusters$V1 %in% c("Old_2")] <- "Young-like Old"
cellStates[pbulkClusters$V1 %in% c("Old_4")] <- "Old Multi-lineage"
cellStates[pbulkClusters$V1 %in% c("Young_3")] <- "Young Mk-biased"
cellStates[pbulkClusters$V1 %in% c("Young_1", "Young_2")] <- "Young Multi-lineage"
groupCellType(project) <- cellStates

########################
# Differential testing #
########################

# Filter out sites with low signal
subCREFilter <- rowMaxs(subCREMat) > 0.3
subCRERanges <- subCRERanges[subCREFilter]
subCREMat <- subCREMat[subCREFilter,]

# Differential testing for CREs and candidate TF binding sites
CREMat <- groupATAC(project)
featureMatList <- list("CRE" = CREMat, "subCRE" = subCREMat)
diffResults <- list("CRE" = list(), "subCRE" = list())
for(feature in c("subCRE", "CRE")){
  
  featureMat <- featureMatList[[feature]]
  
  pbulkAge <- unname(sapply(colnames(featureMat), function(x){strsplit(x, "_")[[1]][1]}))
  
  # Calculate mean feature scores per site
  featureMeans <- rowMeans(featureMat)
  
  # Calculate the difference between young and old
  if(feature == "CRE"){
    featureDiff <- log2(rowMeans(featureMat[, pbulkAge == "Old"] + 0.01) / rowMeans(featureMat[, pbulkAge == "Young"] + 0.01))
  }else if(feature == "subCRE"){
    featureDiff <- rowMeans(featureMat[, pbulkAge == "Old"]) - rowMeans(featureMat[, pbulkAge == "Young"])
  }
  
  featurePvals <- twoSampTTest(featureMat[, pbulkAge == "Old"], featureMat[, pbulkAge == "Young"])
  featureFDRs <- p.adjust(featurePvals, method = "fdr")
  featureFDRs[is.na(featureFDRs)] <- 1
  diffResults[[feature]][["pvals"]] <- featurePvals
  diffResults[[feature]][["diff"]] <- featureDiff
  diffResults[[feature]][["FDRs"]] <- featureFDRs
  diffResults[[feature]][["average"]] <- featureMeans
  
}

# Save subCRE differential results
HSCAgingDiffRanges <- subCRERanges
HSCAgingDiffRanges$pvals <- diffResults$subCRE$pvals
HSCAgingDiffRanges$diffMean <- diffResults$subCRE$diff
HSCAgingDiffRanges$FDRs <- diffResults$subCRE$FDRs
saveRDS(HSCAgingDiffRanges, "../../data/mHSCAging10xV3/HSCAgingDiffRanges.rds")

# Match results for subCRE and the corresponding CREs
overlap <- findOverlaps(regions, subCRERanges)
diffCompare <- data.frame(
  CREDiff = -log10(diffResults[["CRE"]][["FDRs"]][overlap@from]) * sign(diffResults[["CRE"]][["diff"]][overlap@from]),
  CREInd = overlap@from,
  subCREDiff = -log10(diffResults[["subCRE"]][["FDRs"]][overlap@to]) * sign(diffResults[["subCRE"]][["diff"]][overlap@to]),
  subCREInd = overlap@to
)

# For each CRE, keep the most differential subCRE
diffCompare <- as.data.frame(diffCompare %>% group_by(CREInd) %>% filter(abs(subCREDiff) == max(abs(subCREDiff))))
write.table(diffCompare, "../../data/mHSCAging10xV3/diff_CRE_subCRE_compare.tsv",
            sep = "\t", row.names = F, quote = F)

###################################################################
# Compare differential target gene RNA with differential CRE/TFBS #
###################################################################

# Find TSSs of differential genes
TSS <- BuenRTools::mm10TSSRanges
TSS$gene_name <- as.character(TSS$gene_name)
diffTSS <- TSS[TSS$gene_name %in% rownames(diffRNA)[diffRNA$padj < 0.1]]
diffTSS <- resize(diffTSS, 1000, fix = "center")

# Keep differential testing results around these TSS sites
diffGeneData <- diffCompare[diffCompare$CREInd %in% findOverlaps(diffTSS, regions)@to, ]

# Retrieve the corresponding gene name and signed log10(fdr)
diffGeneDataRNA <- data.table::rbindlist(pbmcapply::pbmclapply(
  diffGeneData$CREInd,
  function(CREInd){
    diffGene <- as.character(subsetByOverlaps(diffTSS, regions[CREInd])$gene_name)
    RNADiff <- -log10(diffRNA[diffGene, ]$padj) * sign(diffRNA[diffGene, ]$log2FoldChange)
    filter <- abs(RNADiff) == max(abs(RNADiff))
    data.frame(gene = diffGene[filter], RNADiff = RNADiff[filter])
  },
  mc.cores = 16
))
diffGeneData <- as.data.frame(cbind(diffGeneData, diffGeneDataRNA))

# Visualize the comparison
plotX <- diffGeneData$CREDiff
plotY <- diffGeneData$subCREDiff
plotX[!is.finite(plotX)] <- 0
plotY[!is.finite(plotY)] <- 0
density <- get_density(plotX, plotY)
plotData <- data.frame(CREDiff = plotX, subCREDiff = plotY, density = density,
                       RNADiff = diffGeneData$RNADiff)
pdf("../../data/mHSCAging10xV3/plots/diffCompareCRESubstructure.pdf",
    width = 5.5, height = 5)
ggplot(plotData) + 
  geom_point(aes(x = CREDiff, y = subCREDiff, color = density ^ 0.8), size = 0.5) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  xlab("Differential CRE signed log10(FDR)") + ylab("Differential subCRE signed log10(FDR)") +
  ylim(-15, 15) + xlim(-15, 15) +
  geom_vline(xintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_vline(xintercept = -1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = -1, linetype ="dashed", 
             color = "black", size = 0.5) +
  theme_classic()
dev.off()

# Load pathway gene sets from hypeR
library(hypeR)
c5GO <- msigdb_gsets(species = "Mus musculus","C5","BP",clean = TRUE)$genesets
names(c5GO) <- stringr::str_replace_all(names(c5GO), " ",  "_")
classI <- (diffGeneData$CREDiff > 1) & (diffGeneData$subCREDiff > 1) 
classII <- (abs(diffGeneData$CREDiff) < 1) & (diffGeneData$subCREDiff > 1) 
enrichment <- pathwayEnrichment(
  fgGenes = diffGeneData$gene[classII], 
  bgGenes = diffGeneData$gene[classII | classI],
  geneSets = c5GO,
  pvalThrshold = 0.05)
as.data.frame(enrichment[, c("pathway", "pvals", "fdrs")])
write.table(enrichment, "../../data/mHSCAging10xV3/classIIPathways.tsv",
            sep = "\t", quote = F)

# Visualize results
plotData <- as.data.frame(enrichment[10:1, ])
plotData$pathway <- stringr::str_replace_all(plotData$pathway, "_", " ")
plotData$pathway <- sapply(plotData$pathway, function(s){gsub('(.{1,25})(\\s|$)', '\\1\n', s)}) # Add newline to long strings
plotData$pathway <- factor(plotData$pathway, levels = plotData$pathway) # This keeps the entries in the original order when plotting
plotData$logP <- -log10(plotData$pval)
pdf("../../data/mHSCAging10xV3/plots/classIIPathways.pdf",
    width = 5.5, height = 5)
ggplot(plotData) +
  geom_bar(aes(x = pathway, y = logP), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Enriched pathways")  + 
  ylab("Log10(p-value)") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 12))    
dev.off()

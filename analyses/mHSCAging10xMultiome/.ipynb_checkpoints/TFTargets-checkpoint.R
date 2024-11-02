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
source("../../code/getTargetPathway.R")
library("ComplexHeatmap")

###################
# Load input data #
###################

# Load region ranges
regions <- readRDS("/data/PRINT/multiScaleFootprinting/data/mHSCAging10xV3/regionRanges.rds")

# Load TF binding scores and reformat into a vector
bw <- rtracklayer::import.bw("/data/PRINT/multiScaleFootprinting/data/mHSCAging10xV3/mouse_HSC_TFBS.bigwig")

# Generate a position-by-sample matrix of TF bindings scores
TFBindingSE <- SummarizedExperiment(
  assays = list(TFBS = as.matrix(bw$score)),
  rowRanges = bw
)

# Load pathway gene sets from hypeR
library(hypeR)
c5GO <- msigdb_gsets(species = "Mus musculus","C5","BP",clean = TRUE)$genesets
names(c5GO) <- stringr::str_replace_all(names(c5GO), " ",  "_")

# Get CRE to gene mapping
TSS <- FigR::mm10TSSRanges
TSSCREOv <- findOverlaps(resize(TSS, 1000, fix = "center"), regions)
CREGeneMapping <- data.frame(
  Gene = as.character(TSS$gene_name)[TSSCREOv@from],
  CREInd = TSSCREOv@to,
  CRE = as.character(regions)[TSSCREOv@to],
  Region = as.character(regions[TSSCREOv@to])
)

##########################
# Get TF motif positions #
##########################

# Load PWM data
cisBPMotifs <- readRDS("/data/PRINT/multiScaleFootprinting/data/shared/cisBP_mouse_pwms_2021.rds")

# Find motif matches for all TFs
motifPath <- "/data/PRINT/multiScaleFootprinting/data/mHSCAging10xV3/TFMotifRanges.rds"
if(!file.exists(motifPath)){
  TFMotifRanges <- pbmcapply::pbmclapply(
    names(cisBPMotifs),
    function(TF){
      motifmatchr::matchMotifs(cisBPMotifs[TF], 
                               regions, 
                               genome = "mm10",
                               out = "positions",
                               p.cutoff = 1e-4)[[1]]
    },
    mc.cores = 16
  )
  names(TFMotifRanges) <- names(cisBPMotifs)
  saveRDS(TFMotifRanges, motifPath)
}else{
  TFMotifRanges <- readRDS(motifPath)
}

#########################################
# Find binding sites of a particular TF #
#########################################

TF <- "E2f1"
enrichment <- mapTargets(
  TFMotifSites = TFMotifRanges[TF], 
  siteSE = TFBindingSE,
  regionGeneCorr = CREGeneMapping,
  geneSets = c5GO,
  threshold = 0.75,
  genome = "mm10")

enrichment[[TF]][1:20, -7]

# Visualize results
plotData <- as.data.frame(enrichment[[TF]][10:1, ])
plotData$pathway <- stringr::str_replace_all(plotData$pathway, "_", " ")
plotData$pathway <- sapply(plotData$pathway, function(s){gsub('(.{1,30})(\\s|$)', '\\1\n', s)}) # Add newline to long strings
plotData$pathway <- factor(plotData$pathway, levels = plotData$pathway) # This keeps the entries in the original order when plotting
plotData$logP <- -log10(plotData$pval)
system("mkdir ../../data/mHSCAging10xV3/plots/TFPathwayMapping")
pdf(paste0("../../data/mHSCAging10xV3/plots/TFPathwayMapping/", TF, "_target_pathways.pdf"),
    width = 8, height = 6)
ggplot(plotData) +
  geom_bar(aes(x = pathway, y = logP), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Enriched pathways")  + 
  ylab("Log10(p-value)") +
  coord_flip() +
  theme_classic() +
  ggtitle(TF) +
  theme(axis.text = element_text(size = 10))  
dev.off()

#######################################
# Get TF regulators of the gene lists #
#######################################

TFs <- names(cisBPMotifs)
enrichmentList <- mapTargets(
  TFMotifSites = TFMotifRanges, 
  siteSE = TFBindingSE,
  regionGeneCorr = CREGeneMapping,
  geneSets = c5GO,
  threshold = 0.9,
  genome = "mm10")
saveRDS(enrichmentList, "../../data/mHSCAging10xV3/regulatorEnrichment.rds")
gc()

####################################
# Get pathway-to-TF mapping matrix #
####################################

enrichmentList <- readRDS("../../data/mHSCAging10xV3/regulatorEnrichment.rds")

# Generate pathway-by-TF regulation matrix
pathways <- names(c5GO)
enrichmentMat <- pbmcapply::pbmcmapply(
  function(TF){
    enrichment <- enrichmentList[[TF]]
    enrichVec <- rep(0, length(pathways))
    names(enrichVec) <- pathways
    enrichVec[enrichment$pathway] <- -log10(enrichment$fdrs)
    enrichVec
  },
  names(enrichmentList),
  mc.cores = 16
)
enrichmentMatFilt <- enrichmentMat[rowMaxs(enrichmentMat) > 5, colMaxs(enrichmentMat) > -log10(0.25)]

############################################
# Get top regulators of a specific pathway #
############################################

pathway <- "Oxidative_Phosphorylation"
enrichmentMatScaled <- t(t(enrichmentMat) / colMaxs(enrichmentMat))
sortedHits <- sort(enrichmentMatScaled[pathway,], decreasing = T)
plotData <- data.frame(
  TF = names(sortedHits),
  score = sortedHits
)

# Visualize results
plotData <- plotData[10:1, ]
plotData$TF <- factor(plotData$TF, levels = plotData$TF) # This keeps the entries in the original order when plotting
system("mkdir ../../data/mHSCAging10xV3/plots/TFPathwayMapping")
pdf(paste0("../../data/mHSCAging10xV3/plots/TFPathwayMapping/", pathway, "_regulator_TFs.pdf"),
    width = 8, height = 6)
ggplot(plotData) +
  geom_bar(aes(x = TF, y = score), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Top inferred regulator TFs")  + 
  ylab("Normalized enrichment") +
  coord_flip() +
  theme_classic() +
  ggtitle(TF) +
  theme(axis.text = element_text(size = 10))  
dev.off()

######################
# Identify aging TFs #
######################

# Selected age-associated TFs identified by Seq2PRINT 
upTFs <- sort(Reduce(c, lapply(
    readLines("../../data//mHSCAging10xV3/seqTF_foot_stat_aging_TF_grouping_up.txt"),
    function(x){
        stringr::str_split(x, ",")[[1]]
    }
)))

downTFs <- sort(Reduce(c, lapply(
    readLines("../../data//mHSCAging10xV3/seqTF_foot_stat_aging_TF_grouping_down.txt"),
    function(x){
        stringr::str_split(x, ",")[[1]]
    }
)))

# Filter based on age-associated differential RNA
diffRNA <- read.table("../../data/mHSCAging10xV3/diffRNA.tsv", header = T, sep = "\t")
diffRNATFs <- intersect(
  colnames(enrichmentMatFilt), 
  rownames(diffRNA)[diffRNA$padj < 1e-1]
)
upTFs <- intersect(upTFs, diffRNATFs[diffRNA[diffRNATFs, ]$log2FoldChange > 0])
downTFs <- intersect(downTFs, diffRNATFs[diffRNA[diffRNATFs, ]$log2FoldChange < 0])

#####################
# Visualize results #
#####################

# The complete TF-to-pathway mapping
pltMtx <- enrichmentMatFilt
pltMtx <- pmin(pltMtx, 10)
pdf("../../data/mHSCAging10xV3/plots/TFPathwayMapping/TFPathwayMapping.pdf",
    height = 200, width = 200)
Heatmap(
  pltMtx,
  col = BuenColors::jdb_palette("solar_rojos"),
  row_names_max_width = unit(20, "cm"),
  column_names_max_height = unit(20, "cm")
)
dev.off()

# Visualize aging TFs
pdf("../../data/mHSCAging10xV3/plots/TFPathwayMapping/seq2PRINT_aging_TF_to_pathway_full.pdf",
    width = 150, height = 200)
plotMtx <- enrichmentMat[rowMaxs(enrichmentMat) > 2, c(upTFs, downTFs)]
plotMtx <- t(t(plotMtx) / colMaxs(plotMtx))
plotMtx <- plotMtx[, !is.na(colMaxs(plotMtx))]
plotMtx <- plotMtx[rowSums(plotMtx > 0.2) > 1, ]
diffSign <- c("Down", "Up")[(colnames(plotMtx) %in% upTFs) + 1]
colors <- circlize::colorRamp2(seq(0, quantile(plotMtx, 0.99),length.out=9),
                               colors = BuenColors::jdb_palette("solar_rojos"))
Heatmap(
  plotMtx,
  col = colors,
  column_split = diffSign,
  row_names_max_width = unit(20, "cm"),
  column_names_max_height = unit(20, "cm")
)
dev.off()

##########################
# Plot selected pathways #
##########################

# Plot a subset of pathways regulated by aging-up-regulated TFs
selectedPathways <- c(
  "Macroautophagy", "Regulation_Of_Macroautophagy", "Regulation_Of_Autophagy", 
  "Endoplasmic_Reticulum_Unfolded_Protein_Response", "Cellular_Response_To_Topologically_Incorrect_Protein",
  "Response_To_Endoplasmic_Reticulum_Stress",
  "Proteasomal_Protein_Catabolic_Process", "Proteolysis", "Protein_Catabolic_Process",
  "Mitochondrion_Organization", "Mitochondrial_Translation", "Mitochondrial_Transport",
  "Respiratory_Electron_Transport_Chain", "Mitotic_Cell_Cycle",
  "Regulation_Of_Cell_Cycle", "Ribonucleoprotein_Complex_Biogenesis", "Ribosome_Biogenesis", "Rrna_Metabolic_Process",
  "Double_Strand_Break_Repair", "Recombinational_Repair", 
  "Response_To_Virus", "Response_To_Type_I_Interferon", "Innate_Immune_Response", "Response_To_Biotic_Stimulus",
  "Microtubule_Based_Process", "Microtubule_Cytoskeleton_Organization", "Cilium_Movement", "Cell_Projection_Assembly",
  "Regulation_Of_Lymphocyte_Mediated_Immunity", "Regulation_Of_Cd4_Positive_Alpha_Beta_T_Cell_Differentiation",
  "Alpha_Beta_T_Cell_Activation", "Antigen_Processing_And_Presentation_Of_Peptide_Antigen_Via_Mhc_Class_I",
  "Telomere_Organization", "Positive_Regulation_Of_Cytokine_Production"
)
selectedPathways <- intersect(selectedPathways, rownames(enrichmentMat))

pdf("../../data/mHSCAging10xV3/plots/TFPathwayMapping/seq2PRINT_aging_TF_to_pathway.pdf",
    width = 24, height = 12)
plotMtx <- enrichmentMat[selectedPathways, c(upTFs, downTFs)]
plotMtx <- plotMtx[, colMaxs(plotMtx) > -log10(0.2)]
plotMtx <- t(t(plotMtx) / colMaxs(plotMtx))
plotMtx <- plotMtx[, !is.na(colMaxs(plotMtx))]
rownames(plotMtx) <- sapply(rownames(plotMtx), function(x){stringr::str_replace_all(x, "_", " ")})
colors <- circlize::colorRamp2(seq(0, quantile(plotMtx, 0.99),length.out=9),
                         colors = BuenColors::jdb_palette("solar_rojos"))
diffSign <- c("Down", "Up")[(colnames(plotMtx) %in% upTFs) + 1]
Heatmap(
  plotMtx,
  col = colors,
  column_split = diffSign,
  name = "Rescaled\n enrichment",
  cluster_rows = F,
  rect_gp = gpar(col = "grey", lwd = 2),
  column_names_max_height = unit(10, "cm"),
  row_names_max_width = unit(10, "cm")
)
dev.off()
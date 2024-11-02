# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(patchwork)
source("../../../code/utils.R")

###################
# Load input data #
###################

# Load MPRA data
MPRAData <- read.table("../../../data/sequenceModels/variants/Kircher_et_al_2019/GRCh38_ALL.tsv", header = 1)

# Correct typo from the original dataset
MPRAData$Element <- stringr::str_replace_all(MPRAData$Element, "GP1BA", "GP1BB")

# Write variants to a vcf file (1-based indexing)
vcf <- data.frame(
  chr = paste0("chr", MPRAData$Chromosome), 
  pos = MPRAData$Position,
  ID = Reduce(function(x, y){paste(x, y, sep = "-")}, 
              list(paste0("chr", MPRAData$Chromosome), MPRAData$Position,
                   MPRAData$Ref, MPRAData$Alt)),
  ref = MPRAData$Ref,
  alt = MPRAData$Alt,
  qual = ".", filt = ".")
write.table(
  vcf, "../../../data/sequenceModels/variants/Kircher_et_al_2019/MPRA.vcf",
  sep = "\t", quote = F, row.names = F, col.names = F)

# Load MPRA metadata for each target locus
metadata <- read.csv("../../../data/sequenceModels/variants/Kircher_et_al_2019/MPRA_Info.txt", sep = "\t")

########################
# Load Enformer scores #
########################

scores <- list()

# Match cell types
cellTypeMatching <- setNames(
  c("HepG2", "HepG2", "HepG2", "K562", "K562", "K562", "K562", "K562"), # 
  c("F9", "LDLR", "SORT1", "PKLR", "BCL11A+58", "GP1BB", "HBB", "HBG1") #  
)
genes <- names(cellTypeMatching)

# Load Enformer scores for MPRA variants
EnformerResults <- data.table::fread(
  "../../../data/sequenceModels/variants/Kircher_et_al_2019/MPRA_scores.tsv", 
  sep = "\t", nThread = 16, showProgress = T)
EnformerMat <- as.data.frame(EnformerResults[, 7:ncol(EnformerResults)])
EnformerRanges <- GRanges(
  seqnames = EnformerResults$chrom,
  ranges = IRanges(
    start = EnformerResults$pos,
    end = EnformerResults$pos))

# For each target gene with a matched cell type, retrieve 
# Enformer scores for each MPRA varianty
features <- colnames(EnformerMat)
scores[["Enformer"]] <- pbmcapply::pbmclapply(
  genes,
  function(gene){
    
    # Find cell type-matched features
    cellType <- unname(cellTypeMatching[gene])
    matchedFeatures <- features[
      (stringr::str_detect(features, "DNASE") | stringr::str_detect(features, "CAGE")) & 
        (stringr::str_detect(features, cellType))]
    filtMat <- EnformerMat[, matchedFeatures]
    
    # Use first PC as summary statistic
    #pca <- prcomp(t(scale(filtMat)), retx = T)
    #PC1 <- pca$rotation[, 1]
    
    # Only keep results for the current locus
    locus <- metadata[metadata$Name == gene, ]$Genomic.coordinates..GRCh38.
    locus <- GRanges(stringr::str_replace_all(locus, ",", ""))
    variantInd <- findOverlaps(EnformerRanges, locus)@from
    
    # Compile results
    results <- EnformerResults[variantInd, 2:6]
    results$scores <- rowMeans(EnformerMat[variantInd, matchedFeatures, drop = F])
    #results$scores <- PC1[variantInd]
    colnames(results) <- c("Chromosome", "Position", "ID", "Ref", "Alt", "Score")
    
    # Remove duplicate entries
    # (Duplicates might arise because in the original MPRA dataset the same gene has multiple experiments)
    results <- as.data.frame(results[!duplicated(results$ID), ])
    rownames(results) <- results$ID
    results
  },
  mc.cores = 12
)
names(scores[["Enformer"]]) <- genes

##########################
# Load chromBPNet scores #
##########################

#####################
# Load PRINT scores #
#####################

PRINTScores <- read.table(paste0("../../../data/sequenceModels/variants/Kircher_et_al_2019/MPRA_res.tsv"), sep = ",", header = T)

geneRemap <- list(
  "MYC (rs6983267)" = "MYCrs6983267",
  "HNF4A (P2)" = "HNF4A",
  "BCL11A+58" = "BCL11A",
  "PKLR" = "PKLR-24h",
  "TERT" = "TERT-HEK")

genes <- names(cellTypeMatching)
scores[["PRINT"]] <- pbmcapply::pbmclapply(
  genes,
  function(gene){
    
    # Find the corresponding locus
    locus <- metadata[metadata$Name == gene, ]$Genomic.coordinates..GRCh38.
    locus <- GRanges(stringr::str_replace_all(locus, ",", ""))
    
    # Find the best-matched cell type
    cellType <- cellTypeMatching[gene]
    
    # Only keep results for the current locus
    df <- as.data.frame(PRINTScores %>% 
                          filter(chrom == as.character(seqnames(locus))) %>%
                          filter(pos >= start(locus)) %>%
                          filter(pos <= end(locus))) 
    
    # Retrive sequence attribution score in the corresponding cell type
    # Model was trained with 5-fold validation so we average the scores across all folds
    df$Score <- rowMeans(df[, paste0(cellType, "_fold", 0:4, "_tf_footprint")])
    df <- df[, c("chrom", "pos", "REF", "ALT", "Score")]
    colnames(df) <- c("Chromosome", "Position", "Ref", "Alt", "Score")
    
    # Assign an ID to each variant
    df$ID <- Reduce(
      function(x, y){paste(x, y, sep = "-")}, 
      list(df$Chromosome, df$Position, df$Ref, df$Alt))
    
    # The same locus could be included in multiple experiments, for example SORT/SORT1-flip/SORT1.2
    # This results in some duplicate entries in the vcf. We here only keep unique entries
    df <- df[match(unique(df$ID), df$ID), ]
    rownames(df) <- df$ID
    
    df
  },
  mc.cores = 12
)
names(scores[["PRINT"]]) <- genes

########################
# Retrieve MPRA scores #
########################

scores[["MPRA"]] <- pbmcapply::pbmclapply(
  genes,
  function(gene){
    
    if(gene %in% names(geneRemap)){gene <- geneRemap[[gene]]}
    
    # Retrieve MPRA data for the current target gene
    targetMPRA <- MPRAData[MPRAData$Element == gene, ]
    
    # Find the genomic range of the gene
    targetStart <- min(targetMPRA$Position)
    targetEnd <- max(targetMPRA$Position)
    targetRange <- GRanges(paste0("chr", targetMPRA$Chromosome[1], ":", targetStart, "-", targetEnd))
    
    # Filter by pval and number of tags. Also keep start and end position
    filter1 <- (targetMPRA$Tags > 10) & (targetMPRA$P.Value < 0.05) 
    filter2 <- targetMPRA$Position %in% c(targetStart, targetEnd) 
    targetMPRA <- targetMPRA[filter1 | filter2, ]
    
    # Rename columns and add ID
    colnames(targetMPRA) <- c(
      "Chromosome", "Position", "Ref", "Alt", "Tags", 
      "DNA", "RNA", "Score", "P.Value", "Element")
    targetMPRA$ID <- Reduce(
      function(x, y){paste(x, y, sep = "-")}, 
      list(paste0("chr",targetMPRA$Chromosome), targetMPRA$Position, targetMPRA$Ref, targetMPRA$Alt))
    rownames(targetMPRA) <- targetMPRA$ID
    targetMPRA
  },
  mc.cores = 16
) 
names(scores[["MPRA"]]) <- genes

##################################################
# Benchmark across methods and visualize results #
##################################################

methods <- c("Enformer", "MPRA", "PRINT")
gene <- "HBB"
plots <- lapply(
  methods,
  function(method){
    results <- scores[[method]][[gene]]
    results$conversion <- paste0(results$Ref, "-", results$Alt)
    ggplot(results, aes(x=Position, y=Score)) +
      geom_point(aes(shape = Ref, fill = conversion, color = conversion)) + 
      geom_segment(aes(x=Position, xend=Position, y=0, yend=Score, color = conversion)) +
      theme_classic() +
      ggtitle(paste(method, gene))
  }
)
names(plots) <- methods
plots$Enformer / plots$MPRA / plots$PRINT

# Benchmark model prediction on MPRA as ground truth
method <- "PRINT"
cors <- sapply(
  genes,
  function(gene){
    ovIDs <- intersect(scores[[method]][[gene]]$ID, scores$MPRA[[gene]]$ID)
    predScores <- scores[[method]][[gene]][ovIDs,]$Score
    MPRAScores <- scores$MPRA[[gene]][ovIDs,]$Score
    cor(predScores, MPRAScores)
  }
)
data.frame(gene = names(cors), correlation = unname(cors))
mean(cors)

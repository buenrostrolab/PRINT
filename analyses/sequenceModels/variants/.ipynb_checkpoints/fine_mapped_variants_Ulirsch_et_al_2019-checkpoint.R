# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(GenomicRanges)
library(rtracklayer)
library(PRROC)
library(Matrix)
set.seed(42)

###################
# Load input data #
###################

# Load fine-mapped variant data
variantInfo <- read.table(
  "../../../data/sequenceModels/variants/Ulirsch_et_al_2019/fine_mapped_variants.txt", 
  sep = "\t", header = 1)

# Convert to GRanges object
variants <- GRanges(
  seqnames = variantInfo$chr, 
  ranges = IRanges(start = variantInfo$start, end = variantInfo$start))
variants$ref <- variantInfo$ref
variants$alt <- variantInfo$alt
variants$PP <- variantInfo$PP
variants$trait <- variantInfo$trait
variants$ID <- paste(GRanges(variants), paste(variants$ref, variants$alt, sep = "-"), sep = "-")
variants <- sort(variants)

# Load chain file for liftOver
pathToChain <- "../../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Lift over from hg19 to hg38
seqlevelsStyle(variants) <- "UCSC"  # necessary
variants <- unlist(rtracklayer::liftOver(variants, ch))

# Also save to bed file for visualization
rtracklayer::export(
  variants[variants$PP > 0.9], format = "bed",
  con = "../../../data/sequenceModels/variants/Ulirsch_et_al_2019/variants.bed")

# Define a list of positive (functional) and negative (non-functional) variants 
posVariants <- variants[variants$PP > 0.5]
negVariants <- variants[variants$PP < 0.005]

# Only keep variants overlapping with cCRE regions
regions <- readRDS("../../../data/BMMC/regionRanges.rds")
resizedRegions <- resize(regions, 2114, fix = "center")
posVariants <- subsetByOverlaps(posVariants, regions)
negVariants <- subsetByOverlaps(negVariants, regions)

# Only keep unique variants
posVariants <- posVariants[match(unique(posVariants$ID), posVariants$ID)]

# Down sample negative variants to the same number as positive variants
negVariantID <- sample(unique(negVariants$ID), length(posVariants))
negVariants <- negVariants[match(negVariantID, negVariants$ID)]

################################
# Score variants with Enformer #
################################

# Write VCF file 
variantList <- list("pos" = posVariants, "neg" = negVariants)
for(type in names(variantList)){
  var <- variantList[[type]]
  vcf <- data.frame(
    chr = as.character(seqnames(var)),
    pos = start(var),
    ID = as.character(var),
    ref = var$ref, alt = var$alt,
    qual = ".", filt = "."
  )
  write.table(
    vcf, 
    paste0("../../../data/sequenceModels/variants/Ulirsch_et_al_2019/", type, "_variants.vcf"),
    row.names = F, col.names = F, quote = F, sep = "\t")
}

# Next run enformer_variant_scoring_Ulirsch_et_al.ipynb

# Load from VCF file (this makes sure the pos/neg varianst we use are the same that's scored by enformer)
posVariants <- read.table("../../../data/sequenceModels/variants/Ulirsch_et_al_2019/pos_variants.vcf")
negVariants <- read.table("../../../data/sequenceModels/variants/Ulirsch_et_al_2019/neg_variants.vcf")
posVariants <- GRanges(posVariants$V3)
negVariants <- GRanges(negVariants$V3)

# Load the Enformer scores for each variant
posScores <- read.csv("../../../data/sequenceModels/variants/Ulirsch_et_al_2019/pos_variants_scores.tsv", sep = "\t")
negScores <- read.csv("../../../data/sequenceModels/variants/Ulirsch_et_al_2019/neg_variants_scores.tsv", sep = "\t")

# Get variant-by-feature prediction score matrix
enformerScores <- rbind(posScores, negScores)
scoreMat <- enformerScores[, 7:ncol(enformerScores)]

# Convert variants to GRanges objects
posVariants <- GRanges(posScores$id)
negVariants <- GRanges(negScores$id)
enformerScores <- GRanges(enformerScores$id)

# K-means clustering of all Enformer output features
set.seed(42)
print("K-means clustering of motifs")
kmClustering <- kmeans(t(scale(scoreMat)), 10, iter.max = 100, nstart = 5)
clusterLabels <- kmClustering$cluster

# Manually annotate each cluster
# Cluster 1: HEK; Cluster 2: H3K27me3; Cluster 4: MCF-7; Cluster 5: HepG2 ;Cluster 6: K562
# Cluster 7: H3K4me3; Cluster 8: Hematopoietic cells; Cluster 9: H3K9me3; Cluster 10: CAGE
colnames(scoreMat)[clusterLabels == 8]
sapply(
  sort(unique(clusterLabels)),
  function(cluster){
    mean(stringr::str_detect(colnames(scoreMat), "CAGE")[clusterLabels == cluster])
  }
)

# Keep relevant features and get a final score for each variant
enformerScores$score <- matrixStats::rowMaxs(as.matrix(abs(scoreMat[, clusterLabels == 8])))

###############################################################################
# Evaluate ability of model to distinguish functional/non-functional variants #
###############################################################################

scores <- list()
scores[["chromBPNet"]] <- import.bw("../../../data/sequenceModels/chromBPNet/BMMC/chrombp_BMMC_early.profile_scores.bw")
scores[["Enformer"]] <- enformerScores

method <- "Enformer"

# Retrieve model scores for functional and non-functional variants
posOv <- findOverlaps(scores[[method]], posVariants)
posScores <- abs(scores[[method]][posOv@from]$score)
negOv <- findOverlaps(scores[[method]], negVariants)
negScores <- abs(scores[[method]][negOv@from]$score)

# Calculate AUPRC
pr <- pr.curve(scores.class0 = posScores, scores.class1 = negScores, curve = T)
pr$auc.integral
plot(pr)

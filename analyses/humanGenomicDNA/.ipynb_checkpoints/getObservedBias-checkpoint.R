# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(GenomicRanges)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(patchwork)

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/visualization.R")
source("../../code/getBias.R")
source("../../code/getAggregateFootprint.R")

########################
# Get insertion counts #
########################

# Initialize a footprintingProject object
projectName <- "humanGenomicDNA"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Select sample to run
sample <- "gDNA_3_4_rep1"
sizeFilter <- list(
  "gDNA_2_25" = c(2000, 2500), # 2-2.5 kb fragments
  "gDNA_25_3" = c(2500, 3000), # 2.5-3 kb fragments
  "gDNA_3_4_rep1" = c(3000, 4000), # 3-4 kb fragments, replicate 1
  "gDNA_3_4_rep2" = c(3000, 4000), # 3-4 kb fragments, replicate 2
  "gDNA_4_5" = c(4000, 5000) # 4-5 kb fragments
)

# Select the fragments of the right lengths
regionSize <- sizeFilter[[sample]][1]
fragRanges <- rtracklayer::import.bed("../../data/humanGenomicDNA/predictedFrags.bed")
fragWidth <- width(fragRanges)
widthFilter <- (width(fragRanges) > sizeFilter[[sample]][1]) & (width(fragRanges) < sizeFilter[[sample]][2])
fragRanges <- fragRanges[widthFilter]
fragWidthFilt <- fragWidth[widthFilter]

# Set the regionRanges slot
regionRanges(project) <- resize(fragRanges, regionSize, fix = "center")

# Load barcode grouping
barcodeGroups <- data.frame(barcode = sample, group = 1)

# Load fragments
frags <- data.table::fread(paste0("../../data/humanGenomicDNA/fragments/", sample, ".frags.gz"),
                           showProgress = T, nThread = 16)
nFrags <- dim(frags)[1]

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
biasPath <- paste0(projectDataDir, sample, "_predBias.rds")
if(file.exists(biasPath)){
  regionBias(project) <- readRDS(biasPath)
}else{
  project <- getPrecomputedBias(project, nCores = 24)
  saveRDS(regionBias(project), biasPath)
}

##############################
# Down-sample insertion data #
##############################

# Down-sample fragments data
fragDir <- paste0("../../data/humanGenomicDNA/fragments/downSampledFragments/", sample, "/")
if(!dir.exists(fragDir)){
  system(paste("mkdir -p", fragDir))
}
for(downSampleRate in c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){
  print(paste0("Downsampling rate: ", downSampleRate))
  downSampleInd <- sample(1:nFrags, as.integer(nFrags * downSampleRate))
  downSampledFrags <- frags[downSampleInd, ]
  gz <- gzfile(paste0(fragDir, "/fragsDownsample", downSampleRate, ".tsv.gz"), "w")
  write.table(downSampledFrags, gz, quote = F, row.names = F, col.names = F, sep = "\t")
  close(gz)
}

# Get Tn5 insertion count tensor for fragments with and without down-sampling
countsPath <- paste0("../../data/humanGenomicDNA/", sample, "_countTensors.rds")
if(!file.exists(countsPath)){
  counts <- list()
  pathToFrags <- paste0("../../data/humanGenomicDNA/fragments/", sample, ".frags.gz")
  if(dir.exists("../../data/humanGenomicDNA/chunkedCountTensor")){
    system("rm -r ../../data/humanGenomicDNA/chunkedCountTensor")
  }
  counts[["all"]] <- countTensor(getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T))
  for(downSampleRate in c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){
    system("rm -r ../../data/humanGenomicDNA/chunkedCountTensor")
    pathToFrags <- paste0(fragDir, "fragsDownsample", downSampleRate, ".tsv.gz")
    counts[[as.character(downSampleRate)]] <- countTensor(getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T))
  }
  system("rm -r ../../data/humanGenomicDNA/chunkedCountTensor")
  saveRDS(counts, paste0("../../data/humanGenomicDNA/", sample, "_countTensors.rds"))
}else{
  counts <- readRDS(paste0("../../data/humanGenomicDNA/", sample, "_countTensors.rds"))
}

#################
# Load BAC data #
#################

# Also load BAC data
BACRanges <- readRDS("../../data/BAC/tileRanges.rds")
BACBias <- readRDS("../../data/BAC/predBias.rds")
BACCounts <- readRDS("../../data/BAC/tileCounts.rds")$all

###############################################
# Compare model prediction with observed bias #
###############################################

# Before this step, first make sure we run the above code to get counts for each sample
samples <- c(names(sizeFilter), paste0("BAC_rep", 1:5))
corrList <- lapply(
  samples,
  function(sample){
    
    print(paste0("Processing sample ", sample))
    
    # Get count tensor of the current sample
    if(stringr::str_detect(sample, "BAC")){
      
      # Load BAC counts for the current replicate
      repInd <- as.integer(stringr::str_match(sample, "BAC_rep(\\d)")[2])
      counts <- lapply(BACCounts, function(x){x[x$group == repInd, ]})
      
      # Retrieve predicted bias
      predBiasMtx <- BACBias
      regionSize <- dim(BACBias)[2]
      
    }else{
      
      # Load human gDNA counst for the current sample
      counts <- readRDS(paste0("../../data/humanGenomicDNA/", sample, "_countTensors.rds"))$all
      
      # Retrieve predicted bias
      biasPath <- paste0(projectDataDir, sample, "_predBias.rds")
      predBiasMtx <- readRDS(biasPath)
      regionSize <- sizeFilter[[sample]][1]
      
    }
    
    # Calculate coverage for each region
    coverage <- sapply(counts, function(x){sum(x$count)}) / regionSize
    
    # Compare observed Tn5 insertion with predicted bias
    corrs <- pbmcapply::pbmcmapply(
      function(regionInd){
        
        # Calculate observed Tn5 bias
        obsBias <- getObsBias(counts, regionInd, regionSize)$obsBias
        
        # Calculate correlation with predicted value
        predBias <- predBiasMtx[regionInd, ]
        corr <- cor(obsBias, predBias)
        
      },
      1:length(coverage),
      mc.cores = 16
    )
    list(coverage = coverage, corrs = corrs, sample = rep(sample, length(corrs)))
  }
)

# Combine results
coverage <- Reduce(c, lapply(corrList,function(x){x$coverage}))
corrs <- Reduce(c, lapply(corrList,function(x){x$corrs}))
sampleLabel <- Reduce(c, lapply(corrList,function(x){x$sample}))

# Remove NA values (NAs result from cases where localAvg is 0 at certain positions)
filter <- !is.na(corrs)
coverage <- coverage[filter]
corrs <- corrs[filter]
sampleLabel <- sampleLabel[filter]

plotData <- data.frame(
  coverage = coverage,
  sample = sampleLabel,
  corrs = corrs
)
isBAC <- stringr::str_detect(sampleLabel, "BAC")

pdf(paste0("../../data/humanGenomicDNA/plots/coverage_correlation_scatter_gDNA.pdf"), 
    width = 6.5, height = 5)
ggplot(plotData[!isBAC, ][sample(1:sum(!isBAC)),]) +
  geom_point(aes(x = coverage, y = corrs, color = sample), size = 0.01) + 
  xlab("Coverage (insertion per bp)") + ylab("Observed-predicted correlation") +
  geom_hline(yintercept=0.92, linetype="dashed") + 
  ggtitle("Human genomic DNA") +
  theme_classic()
dev.off()

pdf(paste0("../../data/humanGenomicDNA/plots/coverage_correlation_scatter_BAC.pdf"), 
    width = 6.5, height = 5)
ggplot(plotData[isBAC, ][sample(1:sum(isBAC)),]) +
  geom_point(aes(x = coverage, y = corrs, color = sample), size = 0.01) + 
  xlab("Coverage (insertion per bp)") + ylab("Observed-predicted correlation") +
  geom_hline(yintercept=0.92, linetype="dashed") + 
  ggtitle("BAC DNA") +
  theme_classic()
dev.off()

##################################
# Visualize an individual region #
##################################

#counts <- readRDS(paste0("../../data/humanGenomicDNA/gDNA_3_4_rep2_countTensors.rds"))

regionInd <- sample(which(coverage > 15), 1)
predBias <- regionBias(project)[regionInd, ]
Tn5Insertion <- getObsBias(counts, regionInd)$insertion
obsBias <- getObsBias(counts, regionInd)$obsBias
region <- regionRanges(project)[regionInd]

position = start(region):end(region)
plotData <- data.frame(
  x1 = position - 2,
  x2 = position + 2,
  predBias = predBias,
  Tn5Insertion = Tn5Insertion,
  obsBias = obsBias
)
plotData <- plotData[1:1000,]
corr <- cor(plotData$predBias, plotData$obsBias)
plotRegion <- resize(region,1000,fix="center")
p1 <- ggplot(plotData) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = 0, ymax = Tn5Insertion)) +
  ylab("Observed\nTn5 insertion") +
  ggtitle(paste0(as.character(plotRegion),", Pearson correlation:" , round(corr, 2))) +
  theme_classic()
p2 <- ggplot(plotData) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = 0, ymax = obsBias)) +
  ylab("Observed\nTn5 bias") +
  theme_classic()
p3 <- ggplot(plotData) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = 0, ymax = predBias)) +
  ylab("Predicted\nTn5 bias") +
  theme_classic()
pdf(paste0("../../data/humanGenomicDNA/plots/", sample, "_", as.character(plotRegion), ".pdf"),
    width = 6, height = 4)
p1 / p2 / p3
dev.off()

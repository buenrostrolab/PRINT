# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../../code/utils.R")
source("../../../code/getCounts.R")
source("../../../code/getFootprints.R")
source("../../../code/getBias.R")
library(GenomicRanges)

###################
# Load input data #
###################

# Select sample to run
sizeFilter <- list(
  "gDNA_2_25" = c(2000, 2500), # 2-2.5 kb fragments
  "gDNA_25_3" = c(2500, 3000), # 2.5-3 kb fragments
  "gDNA_3_4_rep1" = c(3000, 4000), # 3-4 kb fragments, replicate 1
  "gDNA_3_4_rep2" = c(3000, 4000), # 3-4 kb fragments, replicate 2
  "gDNA_4_5" = c(4000, 5000) # 4-5 kb fragments
)
samples <- names(sizeFilter)

# Load CRE ranges
dataset <- "K562"
regions <- readRDS(paste0("../../../data/", dataset, "/regionRanges.rds"))

# Load chromBPNet predicted bias
chromBPNetBias <- rtracklayer::import.bw(paste0(
  "../../../data/sequenceModels/chromBPNet/", dataset, "/fold_0_predicted_liftoverpeak_bias.bw"))

# Load PRINT predicted bias
PRINTBias <- readRDS(paste0("../../../data/", dataset, "/predBias.rds"))

######################################
# Get Tn5 bias GenomicRanges objects #
######################################

# Get data.frame of PRINT predicted bias
PRINTBiasDF <- data.table::rbindlist(pbmcapply::pbmclapply(
  1:length(regions),
  function(regionInd){
    region <- regions[regionInd]
    data.frame(
      chr = as.character(seqnames(region)),
      pos = start(region):end(region),
      bias = PRINTBias[regionInd, ]
    )
  },
  mc.cores = 16
))

# Convert to GRanges object
PRINTBias <- GRanges(
  seqnames = PRINTBiasDF$chr,
  ranges = IRanges(start = PRINTBiasDF$pos, end = PRINTBiasDF$pos))
PRINTBias$bias <- PRINTBiasDF$bias

# Get observed Tn5 bias
obsBias <- list()
for(sample in samples){
  
  print(paste("Processing", sample))
  
  # Select the fragments of the right lengths
  regionSize <- sizeFilter[[sample]][1]
  fragRanges <- rtracklayer::import.bed("../../../data/humanGenomicDNA/predictedFrags.bed")
  fragWidth <- width(fragRanges)
  widthFilter <- (width(fragRanges) > sizeFilter[[sample]][1]) & (width(fragRanges) < sizeFilter[[sample]][2])
  fragRanges <- fragRanges[widthFilter]

  # Set the regionRanges slot
  gDNARanges <- resize(fragRanges, regionSize, fix = "center")
  
  # Load gDNA Tn5 insertion counts
  gDNACounts <- readRDS(paste0("../../../data/humanGenomicDNA/", sample, "_countTensors.rds"))$all

  # Get data.frame of observed bias
  locaRadius <- 50
  obsBiasDF <- data.table::rbindlist(pbmcapply::pbmclapply(
    1:length(gDNARanges),
    function(regionInd){
      
      # Compute observed bias
      obsBias <- getObsBias(gDNACounts, regionInd, width(gDNARanges[regionInd]))
      
      # Get local coverage
      coverage <- conv(obsBias$insertion, locaRadius) / (locaRadius * 2)
      
      # Return as data.frame
      region <- gDNARanges[regionInd]
      data.frame(
        chr = as.character(seqnames(region)),
        pos = start(region):end(region),
        coverage = coverage,
        obsBias = obsBias$obsBias
      )
    },
    mc.cores = 16
  ))
  obsBiasDF <- obsBiasDF[obsBiasDF$coverage > 10,]
  
  # Convert to GRanges object
  obsBias[[sample]] <- GRanges(
    seqnames = obsBiasDF$chr,
    ranges = IRanges(start = obsBiasDF$pos, end = obsBiasDF$pos))
  obsBias[[sample]]$bias <- obsBiasDF$obsBias
  obsBias[[sample]] <- obsBias[[sample]][!is.na(obsBias[[sample]]$bias)]
}

#######################
# Compare performance #
#######################

plotData <- data.table::rbindlist(pbmcapply::pbmclapply(
  names(sizeFilter),
  function(sample){

    ov <- findOverlaps(obsBias[[sample]], chromBPNetBias)
    chromBPNetCor <- cor(obsBias[[sample]]$bias[ov@from], chromBPNetBias$score[ov@to])
    
    ov <- findOverlaps(obsBias[[sample]], PRINTBias)
    PRINTCor <- cor(obsBias[[sample]]$bias[ov@from], PRINTBias$bias[ov@to])
    
    data.frame(
      method = c("chromBPNet", "PRINT"),
      cor = c(chromBPNetCor, PRINTCor),
      sample = sample
    )
  }
))

pdf(paste0("../../../data/sequenceModels/plots/", dataset, "_gDNA_Tn5_bias_benchmark.pdf"), 
    width = 8, height = 7)
ggplot(plotData, aes(x = sample, y = cor, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Correlation between observed\nand predicted Tn5 bias") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

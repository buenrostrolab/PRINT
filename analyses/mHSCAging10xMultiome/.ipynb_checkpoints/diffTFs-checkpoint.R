# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

myPath <- .libPaths()
myPath <- c(myPath,'/packages')
.libPaths(myPath)

source("../../code/utils.R")

###################
# Load input data #
###################

# Load ATAC peak ranges
regions <- readRDS("../../data/mHSCAging10xV3/regionRanges.rds")

# Load differential RNA results
diffRNA <- read.table("../../data/mHSCAging10xV3/diffRNA.tsv", sep = "\t")

# Load PWM data
cisBPMotifs <- readRDS("/data/PRINT/multiScaleFootprinting/data/shared/cisBP_mouse_pwms_2021.rds")

# Get differential TFs in aging
diffRNATFs <- intersect(names(cisBPMotifs), rownames(diffRNA[(diffRNA$padj < 0.01), ]))

##########################
# Get TF motif positions #
##########################

# Find motif matches for all TFs
bedPath <- "/data/PRINT/multiScaleFootprinting/data/mHSCAging10xV3/TFMotifRanges.bed"
if(!file.exists(bedPath)){
  TFMotifBed <- data.table::rbindlist(pbmcapply::pbmclapply(
    diffRNATFs,
    function(TF){
      sites <- motifmatchr::matchMotifs(cisBPMotifs[TF], 
                               regions, 
                               genome = "mm10",
                               out = "positions")[[1]]
      data.frame(
        chr = as.character(seqnames(sites)),
        start = start(sites) - 1,
        end = end(sites),
        strand = strand(sites),
        TF = TF
      )
    },
    mc.cores = 16
  ))
  write.table(TFMotifBed, bedPath, quote = F, col.names = T, row.names = F, sep = "\t")
}
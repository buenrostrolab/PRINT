# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
library(data.table)
system("mkdir ../../data/BMMCTutorial")

# Read in all fragments
fragPath <- "../../data/BMMC/merged.fragments.tsv.gz"
frags <- data.table::fread(fragPath, sep = "\t", showProgress = TRUE, nrows = Inf)
colnames(frags) <- c("chr", "start", "end", "barcode", "count")

# Read CRE ranges
regions <- readRDS("../../data/BMMC/regionRanges.rds")

# Select example ranges and subset fragments falling into this region
regionsFilt <- regions[178629:178638]
fragsFilt <- frags[frags$chr == as.character(unique(seqnames(regionsFilt))), ]
fragsFilt <- fragsFilt[fragsFilt$start >= min(start(regionsFilt)),]
fragsFilt <- fragsFilt[fragsFilt$start <= max(end(regionsFilt)),]

# Get bed file of selected regions
regionsFiltBed <- data.frame(
  chr = as.character(seqnames(regionsFilt)),
  start = start(regionsFilt),
  end = end(regionsFilt)
)

# Filter barcode grouping
barcodeGroups <- data.table::fread("../../data/BMMC/barcodeGrouping.txt", sep = "\t", showProgress = TRUE, nrows = Inf)
barcodeGroupsFilt <- barcodeGroups[barcodeGroups$barcode %in% fragsFilt$barcode, ]

# Save subsetted data 
fwrite(fragsFilt, file = "../../data/BMMCTutorial/BMMCTutorialFragments.tsv", quote = F, 
       col.names = T, row.names = F, sep = "\t")
fwrite(barcodeGroupsFilt, file = "../../data/BMMCTutorial/barcodeGrouping.txt", quote = F, 
       col.names = T, row.names = F, sep = "\t")
write.table(regionsFiltBed, "../../data/BMMCTutorial/BMMCTutorialRegions.bed", quote = F,
            col.names = T, row.names = F, sep = "\t")

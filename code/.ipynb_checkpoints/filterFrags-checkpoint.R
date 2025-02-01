library("data.table")
library(GenomicRanges)
library(dplyr)
library(SummarizedExperiment)

getCountsFromFrags <- function(fragFile, peaks, barcodeList = NULL, maxFragLength = NULL, 
    addColData = TRUE) 
{
    start_time <- Sys.time()
    if (class(fragFile) == "character") {
        GA <- BuenRTools::fragsToRanges(fragFile, barcodeList = barcodeList, 
            startsAre0based = TRUE)
    }
    else if (class(fragFile) == "GRanges") {
        GA <- fragFile
        if (!is.null(barcodeList)) {
            cat("Retaining only select barcodes specified within list ..\\n")
            GA <- GA[mcols(GA)[, 1] %in% barcodeList]
        }
    }
    if (ncol(mcols(GA)) > 1) {
        if (inherits(mcols(GA)[2][, 1], "character")) {
            colnames(mcols(GA)) <- c("barcodeID", "readID")
        }
        else if (inherits(mcols(GA)[2][, 1], "integer")) {
            colnames(mcols(GA)) <- c("barcodeID", "pcrDup")
        }
    }
    else {
        colnames(mcols(GA)) <- "barcodeID"
    }
    if (!is.null(maxFragLength)) {
        cat("Removing frags with length > ", maxFragLength, " bp ..\n")
        GA <- GA[width(GA) <= maxFragLength]
        if (length(GA) == 0) 
            stop("Fragment filtering resulting in 0 aligned fragments. Please check / change the provided filter size ..\n")
    }
    barcodes <- as.character(GA$barcodeID)
    denom <- table(barcodes)
    uniqueBarcodes <- names(denom)
    id <- factor(barcodes, levels = uniqueBarcodes)
    cat("Finding overlap between peaks and fragments in data ..\n")
    ovPEAKStarts <- findOverlaps(query = peaks, subject = resize(GA, 
        width = 1, fix = "start"))
    ovPEAKEnds <- findOverlaps(query = peaks, subject = resize(GA, 
        width = 1, fix = "end"))
    cat("Filtering for valid fragment-peak overlaps based on cut site start/end coordinates ..\n")
    validHits <- unique.data.frame(rbind(as.data.frame(ovPEAKStarts), 
        as.data.frame(ovPEAKEnds)))
    require(dplyr)
    cat("Generating matrix of counts ..\n")
    countdf <- data.frame(peaks = validHits$queryHits, sample = as.numeric(id)[validHits$subjectHits]) %>% 
        dplyr::group_by(peaks, sample) %>% dplyr::summarise(count = n()) %>% 
        data.matrix()
    m <- Matrix::sparseMatrix(i = c(countdf[, 1], length(peaks)), 
        j = c(countdf[, 2], length(uniqueBarcodes)), x = c(countdf[, 
            3], 0))
    colnames(m) <- uniqueBarcodes
    if (addColData) {
        cat("Computing sample read depth and FRIP ..\n")
        colData <- data.frame(sample = uniqueBarcodes, depth = as.numeric(denom), 
            FRIP = Matrix::colSums(m)/as.numeric(denom), stringsAsFactors = FALSE)
        stopifnot(all.equal(colData$sample, colnames(m)))
        if (any(colData$FRIP > 1)) 
            warning("One or more barcodes ended up with FRIP score > 1 .. check your fragment file as it may contain some abnormally large fragments that should be removed ..\n")
        cat("Generating SummarizedExperiment object ..\n")
        SE <- SummarizedExperiment(rowRanges = peaks, assays = list(counts = m), 
            colData = colData)
    }
    else {
        cat("Generating SummarizedExperiment object ..\n")
        SE <- SummarizedExperiment(rowRanges = peaks, assays = list(counts = m))
    }
    cat("Done!\n")
    end_time <- Sys.time()
    cat("Time elapsed: ", end_time - start_time, units(end_time - 
        start_time), " \n\n")
    return(SE)
}

# Get ATAC peak ranges
peakBed <- read.table("peaks.bed", sep = "\t", header = F)
peaks <- GRanges(seqnames = peakBed$V1,
                 ranges = IRanges(start = peakBed$V2, end = peakBed$V3))

# Get peak-by-cell count matrix
scATAC <- getCountsFromFrags("all.frags.tsv.gz", peaks)

# Filter by FRIP and depth
FRIPThreshold <- 0.3
depthThreshold <- 100
scATACFilt <- scATAC[, (scATAC$FRIP > FRIPThreshold) & (scATAC$FRIP < 1) & (scATAC$depth > depthThreshold)]

# Save scATAC data after filtering
saveRDS(scATACFilt, "scATAC.rds")
write.table(colnames(scATACFilt), "scBarcodes.txt", sep = "\t", quote = F,
            row.names = F, col.names = F)

# Filter the fragments file
frags <- data.table::fread("all.frags.tsv.gz", sep = "\t", showProgress = TRUE) %>% data.frame() 
fragsFilt <- frags %>% filter(V4 %in% colnames(scATACFilt))

# Save results to file
fragsGz <- gzfile("all.frags.filt.tsv.gz", "w")
write.table(fragsFilt, fragsGz, quote = F, sep = "\t", col.names = F, row.names = F)
close(fragsGz)


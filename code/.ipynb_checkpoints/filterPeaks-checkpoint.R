library("data.table")
library(GenomicRanges)

# Get reference genome from command line arguments
args <-  commandArgs(trailingOnly=TRUE)
genome <- args[1]
print(genome)

# Path to summit file(s)
peakCallingFiles <- list.files("peakCalling")
summitFile <- peakCallingFiles[stringr::str_detect(peakCallingFiles, "summits")]

# Path to blacklist bed file that specifies regions to be filtered out.(for particular reference genome) 
# Also get chromosome size file
if(genome == "mm10"){
  blacklist <- "../shared/BlacklistFiles/mm10.full.blacklist.bed"
  chrSizes <- "../shared/mm10.chrom.sizes"
}else if(genome == "hg19"){
  blacklist <- "../shared/BlacklistFiles/hg19.full.blacklist.bed"
  chrSizes <- "../shared/hg19.chrom.sizes"
}else if(genome == "hg38"){
  blacklist <- "../shared/BlacklistFiles/hg38-blacklist.v2.bed"
  chrSizes <- "../shared/hg38.chrom.sizes"
}

summitsToCleanPeaks <- function(summitFiles, peakWidth = 300, blackList, chromSizes, 
    topNPeaks = NULL, useQuantileRanks = FALSE, FDRcutoff = 0.01, 
    outDir, expName, resizeTo = NULL) 
{
    "%ni%" <- Negate("%in%")
    makePeaksDF <- function(summit_files, peak_width, blacklist, 
        chrom_sizes, n = 999999999, fdr_threshold = 0.01, useQuantileScores = FALSE) {
        pad <- round(as.numeric(peak_width)/2)
        cat("Reading in peak summit file(s): \n")
        cat(summit_files, sep = "\n")
        cat(" ....\n\n")
        cat("NOTE: Assuming all start coordinates are 0-based ..\n")
        cat("Padding peak summits by: ", pad, " bp on either side..\n")
        cat("Removing peaks overlapping with blacklisted regions and out of bound peaks based on chromosome sizes ..\n")
        sizedf <- read.table(chrom_sizes, stringsAsFactors = FALSE)
        names(sizedf) <- c("seqnames", "end")
        sizedf$start <- 0
        chrranges <- makeGRangesFromDataFrame(sizedf)
        chrs <- as.character(sizedf[, 1])
        chrs <- chrs[chrs %ni% c("chrY", "chrM", "MT")]
        bdf <- data.frame(fread(paste0(input = blacklist), header = FALSE))
        bg <- makeGRangesFromDataFrame(setNames(data.frame(bdf[, 
            1], bdf[, 2], bdf[, 3]), c("seqnames", "start", "end")))
        l <- lapply(summit_files, function(file) {
            if (tools::file_ext(file) == "gz") {
                dt <- data.table::fread(paste0("zcat < ", file), 
                  stringsAsFactors = TRUE)
            }
            else {
                dt <- data.table::fread(paste0(file), stringsAsFactors = TRUE)
            }
            dt <- as.data.frame(dt)
            dt <- dt[dt$V1 %in% chrs, ]
            dt <- dt[dt$V5 > -1 * log10(fdr_threshold), ]
            if (useQuantileScores) {
                dt$V6 <- trunc(rank(dt$V5))/nrow(dt)
                scoreColIndex <- 6
            }
            else {
                scoreColIndex <- 5
            }
            peaksR <- makeGRangesFromDataFrame(setNames(data.frame(dt[, 
                1], dt[, 2], dt[, 3], dt[, eval(scoreColIndex)]), 
                c("seqnames", "start", "end", "score")), keep.extra.columns = TRUE, 
                starts.in.df.are.0based = TRUE)
            start(peaksR) <- start(peaksR) - pad
            end(peaksR) <- end(peaksR) + pad
            peaksR <- peaksR[!(1:length(peaksR) %in% data.frame(findOverlaps(peaksR, 
                bg))$queryHits)]
            peaksR <- subsetByOverlaps(peaksR, chrranges, type = "within")
            peaksR
        })
        peaks <- do.call(c, l)
        peaks <- sortSeqlevels(peaks)
        peaks <- sort(peaks)
        if (useQuantileScores) {
            cat("Filtering overlapping peaks based on normalized quantiles of peak summit scores ..\n")
        }
        else {
            cat("Filtering overlapping peaks based on peak summit score ..\n")
        }
        keep_peaks <- 1:length(peaks)
        while (!(isDisjoint(peaks[keep_peaks]))) {
            chr_names <- as.character(seqnames(peaks[keep_peaks]))
            starts <- start(peaks[keep_peaks])
            ends <- end(peaks[keep_peaks])
            scores <- mcols(peaks)$score
            overlap_next <- intersect(which(chr_names[1:(length(keep_peaks) - 
                1)] == chr_names[2:(length(keep_peaks))]), which(ends[1:(length(keep_peaks) - 
                1)] >= starts[2:(length(keep_peaks))]))
            overlap_previous <- overlap_next + 1
            overlap_comparison <- scores[keep_peaks[overlap_previous]] > 
                scores[keep_peaks[overlap_next]]
            discard <- keep_peaks[c(overlap_previous[!overlap_comparison], 
                overlap_next[overlap_comparison])]
            keep_peaks <- keep_peaks[keep_peaks %ni% discard]
        }
        cat("Resorting filtered peaks by chromosome  and coordinate..\n")
        peaks <- sortSeqlevels(peaks)
        peaks <- sort(peaks)
        fP <- data.frame(peaks[keep_peaks], rank = 1:length(keep_peaks))
        if (!is.null(n)) {
            if (as.numeric(n) < dim(fP)[1]) {
                nout <- as.numeric(n)
                cat("Returning only the top ", n, " peaks post-filtering ..\n")
            }
            else {
                nout <- dim(fP)[1]
                cat("The number of peaks specified to be returned (topNPeaks) is greater than the number of peaks left after filtering.. no extra filtering was done\n")
            }
        }
        else {
            nout <- dim(fP)[1]
        }
        odf <- head(fP[order(fP$score, decreasing = TRUE), ], 
            nout)
        cat("FINISHED!\n\n")
        cat("Saving merged peak file to output (NOTE: START COORDINATES ARE NOW 1-BASED)..\n")
        return(odf)
    }
    odf <- makePeaksDF(summit_files = summitFiles, peak_width = peakWidth, 
        blacklist = blackList, chrom_sizes = chromSizes, n = topNPeaks, 
        fdr_threshold = FDRcutoff, useQuantileScores = useQuantileRanks)
    outDir <- gsub(x = outDir, pattern = "/$", replacement = "")
    outFile <- paste0(outDir, "/", expName, ".fixedwidthpeaks_", 
        peakWidth, "bp.bed")
    write.table(odf[sort(odf$rank, decreasing = FALSE, index.return = TRUE)$ix, 
        c(1, 2, 3)], file = outFile, col.names = FALSE, row.names = FALSE, 
        sep = "\t", quote = FALSE)
    cat("\nDone!\n")
    if (!is.null(resizeTo)) 
        BuenRTools::reSizePeaks(outFile, targetSize = resizeTo)
}

# Call summitsToCleanPeaks function for peak clean-up
summitsToCleanPeaks(summitFiles = paste0("peakCalling/", summitFile),
                    peakWidth = 800, # Window to use for peak summit padding
                    blackList = blacklist,
                    chromSizes = chrSizes,
                    topNPeaks = NULL, # Filter top N peaks after?
                    useQuantileRanks = FALSE, # Use normalized ranks instead of raw MACS score?
                    FDRcutoff = 0.01, # MACS peak FDR cut-off used for all peaks
                    outDir = "peakCalling/", # Make sure this directory exists first
                    expName= "filtered", # Name to give as prefix to resulting filtered peak bed file
                    resizeTo = 1000 # Re-size peaks after filtering?
)

system("mv peakCalling/filtered.fixedwidthpeaks_800bp_1000reSized.bed peaks.bed")
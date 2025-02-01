library(GenomicRanges)
library(dplyr)
library(reticulate)
library(SummarizedExperiment)
library(ggplot2)

setClassUnion("groupsClass", c("integer", "character"))
setClass("footprintingProject", 
         slots = c(
           projectName = "character", # Name of the project
           dataDir = "character", # Where we store input and output data
           mainDir = "character", # The main folder containing the code/analyses/data subfolders
           fragFile = "character", # Path to the ATAC fragments file
           refGenome = "character", # Reference genome. Example "hg38", "mm10"
           countTensor = "list", # Sample-by-region-by-position 3D Tensor of ATAC insertions
           regionRanges = "GRanges", # Regions within which we wish to perform footprinting
           barcodeGrouping = "data.frame", # Two-column data.frame (barcode and group)
           groups = "groupsClass", # Names of all samples/pseudobulks. Can be integer or character
           regionChunkSize = "integer", # We sometimes chunk the region list into batches. This is batch size
           regionWidth = "integer", # Width of each region. Should be an integer vector
           footprints = "list", # Footprinting results
           regionBias = "matrix", # Region-by-position predicted Tn5 bias matrix. All regions should be the same size
           dispModel = "list", # List of dispersion models for each window size
           groupATAC = "matrix", # Region-by-sample matrix of ATAC data
           groupRNA = "matrix", # Gene-by-sample matrix of RNA data
           groupCellType = "character", # Cell type of each sample
           groupUMAP = "matrix", # UMAP coordinate of each sample
           groupPseudoTime = "numeric", # Pseudotime assigned to each sample
           ATACTracks = "matrix", # Region-by-position matrix of ATAC data
           TFBindingModel = "list", # Model that uses multi-scale footprints to predict TF binding
           predictionModel = "list" # Slot for a general model that uses multi-scale footprints to make predictions
         ))

# Constructor fucntion for our main class
footprintingProject <- function(refGenome, # Reference genome. Example "hg38", "mm10"
                                projectName, # Character. Name of the project
                                mainDir = "../",
                                dataDir = "../") {
  
  new("footprintingProject", 
      projectName = projectName,
      refGenome = refGenome,
      mainDir = mainDir,
      dataDir = paste0(dataDir, projectName, "/"),
      regionChunkSize = 2000L)
  
}

# Method to get region width from the footprintingProject object
setGeneric("regionWidth", function(x) standardGeneric("regionWidth"))

setMethod("regionWidth", "footprintingProject", function(x) x@regionWidth)

# Method to set region width in the footprintingProject object
setGeneric("regionWidth<-", function(x, value) standardGeneric("regionWidth<-"))

setMethod("regionWidth<-", "footprintingProject", function(x, value) {
  x@regionWidth <- value
  x
})

# Method to get region ranges from the footprintingProject object
setGeneric("regionRanges", function(x) standardGeneric("regionRanges"))

setMethod("regionRanges", "footprintingProject", function(x) x@regionRanges)

# Method to set region ranges in the footprintingProject object
setGeneric("regionRanges<-", function(x, value) standardGeneric("regionRanges<-"))

setMethod("regionRanges<-", "footprintingProject", function(x, value) {
  x@regionRanges <- value
  x@regionWidth <- width(x@regionRanges)
  x
})

# Method to get barcode grouping from the footprintingProject object
setGeneric("barcodeGrouping", function(x) standardGeneric("barcodeGrouping"))

setMethod("barcodeGrouping", "footprintingProject", function(x) x@barcodeGrouping)

# Method to set barcode grouping in the footprintingProject object
setGeneric("barcodeGrouping<-", function(x, value) standardGeneric("barcodeGrouping<-"))

setMethod("barcodeGrouping<-", "footprintingProject", function(x, value) {
  x@barcodeGrouping <- value
  x
})

# Method to get region chunk size from the footprintingProject object
setGeneric("regionChunkSize", function(x) standardGeneric("regionChunkSize"))

setMethod("regionChunkSize", "footprintingProject", function(x) x@regionChunkSize)

# Method to set region chunk size in the footprintingProject object
setGeneric("regionChunkSize<-", function(x, value) standardGeneric("regionChunkSize<-"))

setMethod("regionChunkSize<-", "footprintingProject", function(x, value) {
  x@regionChunkSize <- value
  x
})

# Method to get group names from the footprintingProject object
setGeneric("groups", function(x) standardGeneric("groups"))

setMethod("groups", "footprintingProject", function(x) x@groups)

# Method to set group names in the footprintingProject object
setGeneric("groups<-", function(x, value) standardGeneric("groups<-"))

setMethod("groups<-", "footprintingProject", function(x, value) {
  x@groups <- value
  x
})

# Method to get the path to the fragment file from the footprintingProject object
setGeneric("fragFile", function(x) standardGeneric("fragFile"))

setMethod("fragFile", "footprintingProject", function(x) x@fragFile)

# Method to set the path to the fragment file in the footprintingProject object
setGeneric("fragFile<-", function(x, value) standardGeneric("fragFile<-"))

setMethod("fragFile<-", "footprintingProject", function(x, value) {
  x@fragFile <- value
  x
})

# Method to get the directory for storing data from the footprintingProject object
setGeneric("dataDir", function(x) standardGeneric("dataDir"))

setMethod("dataDir", "footprintingProject", function(x) x@dataDir)

# Method to set the directory for storing data for the footprintingProject object
setGeneric("dataDir<-", function(x, value) standardGeneric("dataDir<-"))

setMethod("dataDir<-", "footprintingProject", function(x, value) {
  x@dataDir <- value
  if(!dir.exists(value)) system(paste("mkdir", value))
  x
})

# Method to get the main directory
setGeneric("mainDir", function(x) standardGeneric("mainDir"))

setMethod("mainDir", "footprintingProject", function(x) x@mainDir)

# Method to set the main directory
setGeneric("mainDir<-", function(x, value) standardGeneric("mainDir<-"))

setMethod("mainDir<-", "footprintingProject", function(x, value) {
  x@mainDir <- value
  x
})

# Method to get the reference genome name from the footprintingProject object
setGeneric("refGenome", function(x) standardGeneric("refGenome"))

setMethod("refGenome", "footprintingProject", function(x) x@refGenome)

# Method to set the reference genome name in the footprintingProject object
setGeneric("refGenome<-", function(x, value) standardGeneric("refGenome<-"))

setMethod("refGenome<-", "footprintingProject", function(x, value) {
  x@refGenome <- value
  x
})

# Method to get the Tn5 bias stored in the footprintingProject object
setGeneric("regionBias", function(x) standardGeneric("regionBias"))

setMethod("regionBias", "footprintingProject", function(x) x@regionBias)

# Method to set the Tn5 bias in the footprintingProject object
setGeneric("regionBias<-", function(x, value) standardGeneric("regionBias<-"))

setMethod("regionBias<-", "footprintingProject", function(x, value) {
  if(length(regionRanges(project)) == 0){
    stop("Need to set the regionRanges slot first!")
  }else if(dim(value)[1] != length(regionRanges(project))){
    stop("Number of regions in region ranges and region bias don't match")
  }
  x@regionBias <- value
  x
})

# Method to get the ATAC data for each cell group
setGeneric("groupATAC", function(x) standardGeneric("groupATAC"))

setMethod("groupATAC", "footprintingProject", function(x) x@groupATAC)

# Method to set the the ATAC data for each cell group
setGeneric("groupATAC<-", function(x, value) standardGeneric("groupATAC<-"))

setMethod("groupATAC<-", "footprintingProject", function(x, value) {
  if(length(regionRanges(project)) == 0){
    stop("Need to set the regionRanges slot first!")
  }else if(dim(value)[1] != length(regionRanges(project))){
    stop("Number of regions in region ranges and  bias don't match")
  }
  x@groupATAC <- value
  x
})

# Method to get the RNA data for each cell group
setGeneric("groupRNA", function(x) standardGeneric("groupRNA"))

setMethod("groupRNA", "footprintingProject", function(x) x@groupRNA)

# Method to set the the RNA data for each cell group
setGeneric("groupRNA<-", function(x, value) standardGeneric("groupRNA<-"))

setMethod("groupRNA<-", "footprintingProject", function(x, value) {
  x@groupRNA <- value
  x
})

# Method to get the cell type label for each cell group
setGeneric("groupCellType", function(x) standardGeneric("groupCellType"))

setMethod("groupCellType", "footprintingProject", function(x) x@groupCellType)

# Method to set the the cell type label for each cell group
setGeneric("groupCellType<-", function(x, value) standardGeneric("groupCellType<-"))

setMethod("groupCellType<-", "footprintingProject", function(x, value) {
  x@groupCellType <- value
  x
})

# Method to get the UMAP coordinates for each cell group
setGeneric("groupUMAP", function(x) standardGeneric("groupUMAP"))

setMethod("groupUMAP", "footprintingProject", function(x) x@groupUMAP)

# Method to set the the UMAP coordinates for each cell group
setGeneric("groupUMAP<-", function(x, value) standardGeneric("groupUMAP<-"))

setMethod("groupUMAP<-", "footprintingProject", function(x, value) {
  x@groupUMAP <- value
  x
})

# Method to get the pseudotime for each cell group
setGeneric("groupPseudoTime", function(x) standardGeneric("groupPseudoTime"))

setMethod("groupPseudoTime", "footprintingProject", function(x) x@groupPseudoTime)

# Method to set the pseudotime for each cell group
setGeneric("groupPseudoTime<-", function(x, value) standardGeneric("groupPseudoTime<-"))

setMethod("groupPseudoTime<-", "footprintingProject", function(x, value) {
  x@groupPseudoTime <- value
  x
})

# Method to get the ATAC Tracks for each region
setGeneric("ATACTracks", function(x) standardGeneric("ATACTracks"))

setMethod("ATACTracks", "footprintingProject", function(x) x@ATACTracks)

# Method to set the ATAC Tracks for each 
setGeneric("ATACTracks<-", function(x, value) standardGeneric("ATACTracks<-"))

setMethod("ATACTracks<-", "footprintingProject", function(x, value) {
  x@ATACTracks <- value
  x
})

# Prepares cluster for parallel computing using foreach
prep_cluster <- function(len, # Number of elemnts in the iterable list
                         n_cores = 12 # Number of cores to use
){
  library(doParallel)
  opts <- list()
  pb <- txtProgressBar(min = 0, max = len, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(n_cores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  list("opts" = opts, "cl" = cl)
}

# Chunk a vector/list x into chunks. Return starts and ends of chunks
getChunkInterval <- function(x, # Vector or list
                             chunkSize = 2000 # Size of a single chunk
){
  
  chunkSize <- min(length(x), chunkSize)
  nData <- length(x)
  starts <- seq(1, nData, chunkSize)
  ends <- starts + chunkSize - 1
  ends[length(ends)] <- nData
  
  list("starts" = starts,
       "ends" = ends)
}

# Calculate cosine distance among observations
distCosine <- function(data # observation-by-feature matrix
){
  Matrix <- as.matrix(data)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  dist <- as.dist(1 - sim)
  dist[is.na(dist)] <- 1
  dist
}

# Plotting function used for debugging purposes
plt <- function(data){
  ggplot(data = data.frame(x = 1:length(data), y = data)) +
    geom_line(aes(x = x, y = y))
}

# Performs column-normalization on matrix-like data
# After normalization each column will sum up to the same number
centerCounts <- function(obj, # Input data
                         doInChunks = TRUE, # Whether chunk data when processing 
                         chunkSize = 1000 
){
  if (!class(obj) %in% c("SummarizedExperiment", "RangedSummarizedExperiment", 
                         "dgCMatrix", "dgeMatrix", "Matrix", "matrix")) 
    stop("Supplied object must be either of class SummarizedExperiment, matrix or sparse Matrix ..\n")
  if (ncol(obj) > 10000) 
    doInChunks <- TRUE
  if (doInChunks) {
    cat("Centering counts for cells sequentially in groups of size ", 
        chunkSize, " ..\n\n")
    starts <- seq(1, ncol(obj), chunkSize)
  }
  else {
    starts <- 1
  }
  counts.l <- list()
  for (i in 1:length(starts)) {
    beginning <- starts[i]
    if (i == length(starts)) {
      ending <- ncol(obj)
    }
    else {
      ending <- starts[i] + chunkSize - 1
    }
    cat("Computing centered counts for cells: ", beginning, 
        " to ", ending, "..\n")
    if (class(obj) == "RangedSummarizedExperiment" | class(obj) == 
        "SummarizedExperiment") {
      m <- SummarizedExperiment::assay(obj[, beginning:ending])
    }
    else {
      m <- obj[, beginning:ending]
    }
    obsMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    cCounts <- Matrix::t(Matrix::t(m)/obsMeans)
    counts.l[[i]] <- cCounts
    gc()
  }
  cat("Merging results..\n")
  centered.counts <- do.call("cbind", counts.l)
  cat("Done!\n")
  if (class(obj) == "RangedSummarizedExperiment" | class(obj) == 
      "SummarizedExperiment") {
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  }
  else {
    return(centered.counts)
  }
}

# Merge overlapping regions in a GRanges object
# The input should have $score attribute so we can merge based on this score
mergeRegions <- function(regions # The GRanges object to be merged. It must have a $score attribute
){
  
  "%ni%" <- Negate("%in%")
  
  # Sort
  regions <- sortSeqlevels(regions)
  regions <- sort(regions)
  
  # Filter regions based on summit score/quantile normalized summit score
  keptRegions <- 1:length(regions)
  while (!(isDisjoint(regions[keptRegions]))) {
    
    # Fast variable access
    chrNames <- as.character(seqnames(regions[keptRegions]))
    starts <- start(regions[keptRegions])
    ends <- end(regions[keptRegions])
    scores <- regions$score
    
    # See if consecutive regions are overlapping
    overlapNext <- intersect(
      which(chrNames[1:(length(keptRegions) - 1)] == chrNames[2:(length(keptRegions))]),
      which(ends[1:(length(keptRegions) - 1)] >= starts[2:(length(keptRegions))] )
    )
    
    # Compare consectuive regions
    overlapPrevious <- overlapNext + 1
    overlapComparison <- scores[keptRegions[overlapPrevious]] > scores[keptRegions[overlapNext]]
    discard <- keptRegions[c(overlapPrevious[!overlapComparison], overlapNext[overlapComparison])]
    keptRegions <- keptRegions[keptRegions %ni% discard]
  }
  regions <- sortSeqlevels(regions); regions <- sort(regions)
  
  regions[keptRegions]
  
}

# Load SHARE-seq RNA data from h5 files
ReadCB_h5 <-
  function(filename,
           use.names = TRUE,
           unique.features = TRUE) {
    
    library(hdf5r)
    
    if (!requireNamespace('hdf5r', quietly = TRUE)) {
      stop("Please install hdf5r to read HDF5 files")
    }
    if (!file.exists(filename)) {
      stop("File not found")
    }
    infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
    genomes <- names(x = infile)
    output <- list()
    if (hdf5r::existsGroup(infile, 'matrix')) {
      # cellranger version 3
      message('CellRanger version 3+ format H5')
      if (use.names) {
        feature_slot <- 'features/name'
      } else {
        feature_slot <- 'features/id'
      }
    } else {
      message('CellRanger version 2 format H5')
      if (use.names) {
        feature_slot <- 'gene_names'
      } else {
        feature_slot <- 'genes'
      }
    }
    for (genome in genomes) {
      counts <- infile[[paste0(genome, '/data')]]
      indices <- infile[[paste0(genome, '/indices')]]
      indptr <- infile[[paste0(genome, '/indptr')]]
      shp <- infile[[paste0(genome, '/shape')]]
      features <- infile[[paste0(genome, '/', feature_slot)]][]
      barcodes <- infile[[paste0(genome, '/barcodes')]]
      sparse.mat <- sparseMatrix(
        i = indices[] + 1,
        p = indptr[],
        x = as.numeric(x = counts[]),
        dims = shp[],
        repr = "C"
        # giveCsparse = FALSE
      )
      if (unique.features) {
        features <- make.unique(names = features)
      }
      rownames(x = sparse.mat) <- features
      colnames(x = sparse.mat) <- barcodes[]
      sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
      # Split v3 multimodal
      if (infile$exists(name = paste0(genome, '/features'))) {
        types <- infile[[paste0(genome, '/features/feature_type')]][]
        types.unique <- unique(x = types)
        if (length(x = types.unique) > 1) {
          message(
            "Genome ",
            genome,
            " has multiple modalities, returning a list of matrices for this genome"
          )
          sparse.mat <- sapply(
            X = types.unique,
            FUN = function(x) {
              return(sparse.mat[which(x = types == x),])
            },
            simplify = FALSE,
            USE.NAMES = TRUE
          )
        }
      }
      output[[genome]] <- sparse.mat
    }
    infile$close_all()
    if (length(x = output) == 1) {
      return(output[[genome]])
    } else{
      return(output)
    }
  }

# Cluster TF motifs based on Jaccard distance of their genome-wide motif matching patterns
clusterMotifs <- function(motifs, # PWMatrixList object. 
                          regions, # GRanges object. Specifies the regions for motif matching
                          genome, # Ref genome. Example: "hg38
                          nClusters = 100 # Number of clusters
){
  
  # Find matches of the remaining motifs in regions of interest
  print("Finding matches of the remaining motifs in footprint regions")
  motifMatches <- t(pbmcapply::pbmcmapply(
    function(TF){
      as.matrix(assay(motifmatchr::matchMotifs(subject = regions, 
                                               pwms = motifs[TF], 
                                               genome = genome)))
    },
    names(motifs),
    mc.cores = 16
  ))
  motifMatches <- t(as.matrix(motifMatches * 1)) # Convert to a matrix with 1s and 0s
  
  # Calculating simialrity among motifs
  # Similarity = Jaccard index between motif matches of two TFs
  print("Calculating simialrity among motifs")
  nTFs <- dim(motifMatches)[2]
  intersection <- t(motifMatches) %*% motifMatches
  unmatched <- 1 - motifMatches
  union <- array(dim(unmatched)[1], dim = c(nTFs, nTFs)) - t(unmatched) %*% unmatched
  seqSimilarity <- intersection / union
  
  # K-means clustering
  set.seed(42)
  print("K-means clustering of motifs")
  kmClustering <- kmeans(seqSimilarity, nClusters, iter.max = 100, nstart = 5)
  clusterLabels <- kmClustering$cluster
  
  # Use hierarchical clustering to reorder the K-mean clusters
  clusterCenters <- t(sapply(
    sort(unique(clusterLabels)),
    function(cluster){
      colMeans(seqSimilarity[clusterLabels %in% cluster, , drop = F])
    }
  ))
  hclustTree <- hclust(dist(clusterCenters))
  clusterOrder <- hclustTree$order
  reOrderMap <- 1:nClusters
  names(reOrderMap) <- clusterOrder
  clusterLabels <- reOrderMap[as.character(clusterLabels)]
  motifOrder <- order(clusterLabels)
  
  # Summarize motif clustering results into a data.frame
  motifClustering <- as.data.frame(cbind(names(motifs)[motifOrder], clusterLabels[motifOrder]))
  colnames(motifClustering) <- c("TF", "cluster")
  rownames(motifClustering) <- motifClustering$TF
  
  # Find representative member for each cluster
  representatives <- sapply(
    1:nClusters,
    function(clusterInd){
      clusterTFs <- motifClustering[motifClustering$cluster == clusterInd,]$TF
      names(sort(rowMeans(seqSimilarity[clusterTFs, clusterTFs]), decreasing = T))[1]
    }
  )
  motifClustering$representative <- F
  motifClustering[motifClustering$TF %in% representatives,]$representative <- T
  
  motifClustering
}

# Perform t-test for each feature in the two feature-by-observation matrices x and y
# Here we are using unequal variance Welch test
# For instructions see https://www.statsdirect.co.uk/help/parametric_methods/utt.htm
twoSampTTest <- function(x, # Feature-by-observation matrices x 
                         y, # Feature-by-observation matrices y
                         return_stats=FALSE # Whether to return t-statistics too
){
  if(is.null(dim(x)) & is.null(dim(y))){
    n1 <- length(x)
    n2 <- length(y)
    s1 <- sum((x - mean(x)) ^ 2) / (n1 - 1)
    s2 <- sum((y - mean(y)) ^ 2) / (n2 - 1)
    statistic <- (mean(x) - mean(y)) /  sqrt(s1 / n1 + s2 / n2)
  }else{
    n1 <- dim(x)[2]
    n2 <- dim(y)[2]
    s1 <- rowSums((x - rowMeans(x)) ^ 2) / (n1 - 1)
    s2 <- rowSums((y - rowMeans(y)) ^ 2) / (n2 - 1)
    statistic <- (rowMeans(x) - rowMeans(y)) /  sqrt(s1 / n1 + s2 / n2)
  }
  df <- (s1 / n1 + s2 / n2) ^ 2 / ((s1 / n1) ^ 2 / (n1 - 1) + (s2 / n2) ^ 2 / (n2 - 1))
  statistic[is.na(statistic)] <- 0
  
  pvals <- pt(statistic, df = df)
  pvals[is.na(pvals)] <- 0.5
  pvals[pvals > 0.5] <- 1 - pvals[pvals > 0.5]
  pvals <- pvals * 2 # Two tailed test

  if(return_stats){
      list("pvals"=pvals, "stats"=statistic)
  }else{
      pvals
  }
  
}

# Helper function for making density plot
get_density <- function(x, y, n = 100) {
  library(ggpointdensity)
  dens <- MASS::kde2d(x = x, y = y, n = n, 
                      h = c((max(x) - min(x)) / 50,
                            (max(y) - min(y)) / 50))
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

extractTFNames <- function(motifIDs) 
{
    if (all(grepl("_", motifIDs, fixed = TRUE))) {
        sapply(strsplit(sapply(strsplit(motifIDs, "_LINE.", fixed = FALSE), 
            "[[", 2), "_", fixed = FALSE), "[[", 2)
    }
    else {
        message("One or more provided motif IDs do not contain any '_' characters .. returning IDs as is")
        motifIDs
    }
}

# Calculate motif enrichment in a selected set of region (genomic regions)
motifRegionZTest <-   function(regionSet, ### region_set: vector of region ids you wish to test for motif enrichment
                             bgRegions, ### bg_regions: matrix of background region selection iterations by chromVAR
                             tfMat ### tf: binary matrix of region by motif
) {
  
  if(nrow(tfMat)!=nrow(bgRegions))
    stop("Reference region set used for TF and background regions matrix must match..\n")
  
  if(!all(regionSet %in% 1:nrow(bgRegions)))
    stop("One or more of the provided region indices are out of the background region set range ..\n")
  
  
  # get selected region motif frequencies
  cat("Getting selected region motif frequencies ..\n")
  
  # get frequency of motifs in test set (observed)
  p.tab <- Matrix::colMeans(tfMat[regionSet, ])
  
  # get the background frequency in region sets of the same size
  cat("Getting background region motif frequencies ..\n")
  # extract relevant rows (i.e. regionset being tested) from background region matrix
  bg.f <- as.matrix(bgRegions[regionSet, ])
  
  # calculate (background) motif frequencies in each iteration of background regions corresponding to regionset
  bg.tab <- apply(bg.f[, c(1:ncol(bgRegions))], 2, function(bg_iter) {
    
    b.i <- Matrix::colMeans(tfMat[bg_iter, ])
    return(b.i)
    
  })
  
  cat("Calculating empirical p values and z score p values ..\n")
  
  # loop over each motif and generate enrichment statistics compared to background
  m.p <- dplyr::bind_rows(lapply(names(p.tab), function(motif) {
    
    # calculate sd and mean frequencies for bg and selected regions
    s <- sd(bg.tab[motif, ])
    bg_freq <- mean(bg.tab[motif, ])
    
    z_score <- (p.tab[motif] - bg_freq) / s
    
    # generate data.frame object of relevant statistics
    d <- data.frame(
      motifID = motif,
      gene = extractTFNames(motif),
      motif_obs_freq = p.tab[motif],
      motif_bg_freq = mean(bg.tab[motif, ]),
      motif_counts = p.tab[motif] * length(regionSet),
      emp_pval = 1 - (sum(bg.tab[motif, ] < p.tab[motif]) / ncol(bg.tab)),
      z_test = z_score,
      pval.z = 2 * pnorm(-abs(z_score)),
      signed.log10p = -log10(2 * pnorm(-abs(z_score))) * sign(z_score)
    )
    return(d)
  }))
  # sort by enrichment pval, motif observed frequency
  m.p <- dplyr::arrange(m.p,pval.z, motif_obs_freq)
  # return df of enrichment scores
  return(m.p)
}

# Function for calculating TSS enrichment. We modified code from the ArchR package (https://www.archrproject.com/)
# For details, see https://rdrr.io/github/GreenleafLab/ArchR/src/R/QualityControl.R
getTSSEnrichment <- function(
    frags, # ATAC fragments. Should be a GRanges object of all fragments.
    genome, # One of "hg38", "hg19", "mm10"
    norm = 100,
    flank = 2000
){
  
  library(ArchR)
  rleSumsStranded <- function(rleList, grList, width, as_integer) {
    .Call('_ArchR_rleSumsStranded', PACKAGE = 'ArchR', rleList, grList, width, as_integer)
  }
  
  if(genome == "hg38"){
    TSS <- sort(sortSeqlevels(FigR::hg38TSSRanges))
  }else if(genome == "hg19"){
    TSS <- sort(sortSeqlevels(FigR::hg19TSSRanges))
  }else if(genome == "mm10"){
    TSS <- sort(sortSeqlevels(FigR::mm10TSSRanges))
  }
  
  chr <- paste0(seqnames(TSS))
  chr <- gtools::mixedsort(intersect(chr, paste0(seqnames(TSS))))
  splitTSS <- split(GenomicRanges::resize(TSS,1,"start"), seqnames(TSS))[chr]
  window <- 2 * flank + 1
  
  sumTSS <- pbmcapply::pbmcmapply(
    function(k){
      
      #TSS for Chr
      TSSi <- splitTSS[[chr[k]]]
      
      #Set TSS To be a dummy chr1
      TSSi <- GRanges(seqnames = rep("chr1",length(TSSi)), ranges = ranges(TSSi), strand = strand(TSSi))
      
      #Extract Fragments
      covi <- frags[frags$V1 == chr[k]]
      covi <- GRanges(seqnames = rep("chr1", length(covi$V1)),
                      ranges = IRanges(start = covi$V2, end = covi$V3))
      
      #Get Insertions
      covi <- sort(c(start(covi), end(covi)))
      
      #IRanges
      covi <- IRanges(start = covi, width = 1)
      
      #Coverage
      covi <- IRanges::coverage(covi)
      
      #Compute Sum
      sumTSSi <- rleSumsStranded(list(chr1=covi), list(chr1=TSSi), window, as.integer)
      
    },
    seq_along(chr),
    mc.cores = 2
  )
  
  sumTSS <- rowSums(sumTSS)
  
  normBy <- mean(sumTSS[c(1:norm,(flank*2-norm+1):(flank*2+1))])
  
  TSSEnrichment <- data.frame(
    position = seq_along(sumTSS) - flank - 1, 
    value = sumTSS, 
    normValue = sumTSS / normBy
  )
  
}

# Function for smooything and rescaling tracks
processTrack <- function(x, smoothRadius = 5){
  x <- conv(x, smoothRadius)/(2*smoothRadius)
  x <- x[smoothRadius:(length(x) - smoothRadius)]
  x <- (x - min(x)) / (max(x) - min(x))
  x
}

# Fast running window sum. This runs a window with radius = r across the input vector x and calculates running sum
# For any position i, this sums x[(i - r + 1) : (i + r)]
conv <- function(x, # Input vector x
                 r # Window radius
                 ){
  smoothKernel <- rep(1, 2 * r)
  xConv <- cladoRcpp::rcpp_convolve(x, smoothKernel)
  xConv[(r + 1):(length(x) + r)]
}

fragsToRanges <- function(fragFile, barcodeList = NULL, startsAre0based = TRUE, 
    ...) 
{
    cat("Reading in fragment file ..\n")
    if (is.null(barcodeList)) {
        cat("Reading all barcodes found within file ..\n")
        GA <- data.table::fread(fragFile, sep = "\t", showProgress = TRUE, 
            ...) %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
            start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, 
            starts.in.df.are.0based = startsAre0based)
    }
    else {
        cat("Reading only select barcodes specified within list from file ..\n")
        GA <- data.table::fread(fragFile, sep = "\t", showProgress = TRUE, 
            ...) %>% as.data.frame() %>% dplyr::filter(V4 %in% 
            barcodeList) %>% GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
            start.field = "V2", end.field = "V3", keep.extra.columns = TRUE, 
            starts.in.df.are.0based = startsAre0based)
    }
    return(GA)
}

getCountsFromFrags <- function(fragFile, peaks, barcodeList = NULL, maxFragLength = NULL, 
    addColData = TRUE) 
{
    start_time <- Sys.time()
    if (class(fragFile) == "character") {
        GA <- fragsToRanges(fragFile, barcodeList = barcodeList, 
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

getGeneScoresFromPeaks <- function(SE, geneList = NULL, genome = c("hg19", "hg38", "mm10"), 
    TSSwindow = 10000, getWeightsOnly = FALSE) 
{
    if (length(genome) > 1) 
        stop("Must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
    if (!genome %in% c("hg19", "hg38", "mm10")) 
        stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
    switch(genome, hg19 = {
        TSSg <- FigR::hg19TSSRanges
    }, hg38 = {
        TSSg <- FigR::hg38TSSRanges
    }, mm10 = {
        TSSg <- FigR::mm10TSSRanges
    })
    TSSg$gene_name <- as.character(TSSg$gene_name)
    if (!is.null(geneList)) {
        if (!(all(geneList %in% as.character(TSSg$gene_name)))) 
            stop("One or more of the gene names supplied is not present in the annotation provided..\n")
        cat("Running gene-peak mapping for genes:", geneList, 
            sep = "\n")
        cat("........\n")
        TSSg <- TSSg[TSSg$gene_name %in% geneList, ]
    }
    else {
        cat("Running gene-peak mapping for all genes in annotation! (n = ", 
            length(TSSg), ") This is bound to take more time than querying specific markers ..\n", 
            sep = "")
    }
    cat("Using window of: ", TSSwindow, " bp (total) around TSS per gene ..\n")
    TSSflank <- GenomicRanges::flank(TSSg, width = TSSwindow/2, 
        both = TRUE)
    if (!all(start(TSSflank) > 0)) {
        cat("WARNING: ", sum(start(TSSflank) < 0), " flanked TSS window(s) found to extend beyond the chromosomal boundary .. Resetting start coordinates will result in uneven bins for these TSSs..\n")
        start(TSSflank)[start(TSSflank) < 0] <- 0
    }
    peakg <- GenomicRanges::granges(SE)
    peakSummits <- GenomicRanges::start(peakg) + GenomicRanges::width(peakg)/2
    GenomicRanges::start(peakg) <- GenomicRanges::end(peakg) <- peakSummits
    TSSPeakOverlap <- GenomicRanges::findOverlaps(TSSflank, peakg)
    cat("Determining peak weights based on exponential inverse distance to TSS ..\n")
    distToTSS <- abs(start(peakg)[subjectHits(TSSPeakOverlap)] - 
        start(TSSg)[queryHits(TSSPeakOverlap)])
    weights <- exp(-1 * (distToTSS/1000))
    cat("Assembling Peaks x Genes weights matrix ..\n")
    m <- Matrix::sparseMatrix(i = c(subjectHits(TSSPeakOverlap), 
        length(peakg)), j = c(queryHits(TSSPeakOverlap), length(TSSflank)), 
        x = c(weights, 0))
    colnames(m) <- TSSg$gene_name
    if (getWeightsOnly) 
        return(m)
    weights.m <- m[, Matrix::colSums(m) != 0]
    cat("Assembling Gene x Cells scores matrix ..\n")
    geneScoresMat <- Matrix::t(weights.m) %*% SummarizedExperiment::assay(SE)
    cat("Done!\n\n")
    return(geneScoresMat)
}

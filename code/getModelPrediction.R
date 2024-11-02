library(hdf5r)

# Use multi-scale footprints as input for model prediction
modelPredict <- function(inputData, # Model input data
                         model, # Weights of the neural net model
                         mode = "nucleosome" # Can be either "nucleosome" or "TF"
){
  
  # First layer
  x <- inputData %*% model[[1]] # Linear transform
  x <- t(t(x) + as.numeric(model[[2]])) # Bias
  x[x < 0] <- 0 # ReLU
  
  # Second layer
  x <- x %*% model[[3]] # Linear transform
  x <- t(t(x) + as.numeric(model[[4]])) # Bias
  x[x < 0] <- 0 # ReLU
  
  # Third layer
  x <- x %*% model[[5]] # Linear transform
  x <- t(t(x) + as.numeric(model[[6]])) # Bias
  
  if(mode == "TF"){
    x <- 1 / (1 + exp(-x)) # Sigmoid
  }
  
  x
}

# Method to get the prediction model from the footprintingProject object
setGeneric("predictionModel", function(x) standardGeneric("predictionModel"))

setMethod("predictionModel", "footprintingProject", function(x) x@predictionModel)

# Method to set the prediction model in the footprintingProject object
setGeneric("predictionModel<-", function(x, value) standardGeneric("predictionModel<-"))

setMethod("predictionModel<-", "footprintingProject", function(x, value) {
  x@predictionModel <- value
  x
})

# Get predicted TF binding scores for a specific region
getRegionPred <- function(regionATAC, # Position-by-pseudobulk matrix of ATAC data for the current region
                          Tn5Bias, # Numeric vector of predicted Tn5 bias for the current region
                          region, # GRanges object for the current region
                          dispModels, # Background dispersion model for center-vs-(center + flank) ratio insertion ratio
                          model, # Model for predicting TF binding at any motif site
                          sites = NULL, # GRanges object. Genomic ranges of motif matches in the current region
                          # Must have $score attribute which indicates how well the motif is matched (between 0 and 1)
                          # If NULL, the region will be divided into equally sized tiles for TF binding scoring
                          tileSize = 10, # Size of tiles if sites is NULL
                          #Whether to return model predicted value or confidence in the prediction
                          contextRadius = 100 # Local radius of model input (in bp)
){
  
  if(!is.null(sites)){
    # Only keep motif matches within the region
    # Although we specified that the input sites should already be matches in the current region
    # this step prevents erros in case the user for got to take the overlap
    sites <- subsetByOverlaps(sites, region, type = "within")
  }else{
    sites <- GenomicRanges::slidingWindows(region, width = tileSize, step = tileSize)[[1]]
    sites$score <- 1
    sites$TF <- ""
  }
  
  width <- length(Tn5Bias)
  
  # These are the footprint radii used by the model. 
  # For our first version, we use 10bp, 20bp, 30bp, 50bp, 80bp, 100bp footprints
  scales <- as.character(model$scales)
  
  # Get the pseudobulk-by-position matrix of footprint scores for the current region for each scale
  multiScaleFootprints <- lapply(
    scales,
    function(scale){
      
      smoothRadius <- as.integer(as.integer(scale) / 2)
      
      # Get the pseudobulk-by-position matrix of footprint scores for the current region with the current scale (footprint size)
      footprintScores <- t(sapply(
        1:dim(regionATAC)[2],
        function(groupInd){
          
          # Get Tn5 insertion track for one pseudobulk
          Tn5Insertion <- regionATAC[, groupInd]
          
          # Calculate footprint pval for each position
          pvals <- footprintScoring(
            Tn5Insertion = Tn5Insertion,
            Tn5Bias = Tn5Bias,
            dispersionModel = dispModels[[scale]],
            footprintRadius = as.integer(scale),
            flankRadius = as.integer(scale)
          )
          
          # Convert pvals to scores
          scores <- -log10(pvals)      
          scores[!is.finite(scores)] <- 0
          scores <- caTools::runmax(scores, 2 * smoothRadius)
          scores <- conv(scores, smoothRadius) / (2 * smoothRadius)
          scores
        }
      ))
      
      footprintScores
    }
  )
  names(multiScaleFootprints) <- scales
  
  # Calculate positions of TF sites relative to the start of the CRE region
  relativePos <- start(resize(sites, 1, fix = "center")) - start(region) + 1
  
  # Only keep sites with distance to CRE edge >= contextRadius
  siteFilter <- (relativePos > contextRadius) & (relativePos <= (width - contextRadius))
  sites <- sites[siteFilter]
  relativePos <- relativePos[siteFilter]
  
  # Go through each site and calculate predicted TF binding score for each pseudobulk
  if(length(sites) > 0){
    predScores <- t(sapply(
      1:length(sites),
      function(siteInd){
        siteFootprints <- do.call("rbind", lapply(
          scales,
          function(scale){
            contextInds <- (relativePos[siteInd] - contextRadius) : (relativePos[siteInd] + contextRadius)
            scaleFootprints <- collapse::fsubset(t(multiScaleFootprints[[scale]]), 1:width %in% contextInds)
            #Flip the pattern if the site is on the negative strand
            if(as.character(strand(sites[siteInd])) == "-"){
              scaleFootprints <- scaleFootprints[rev(1:dim(scaleFootprints)[1]), , drop = F]
            }
            scaleFootprints
          }
        ))
        siteFootprints <- t((siteFootprints - model$footprintMean) / model$footprintSd)
        
        # Calculate model prediction scores
        predScore <- modelPredict(siteFootprints, model)
        
      }
    ))
    
    # If there is only 1 pseudobulk, flip the axis to make sure prediction scores is site-by-pseudobulk
    if(dim(regionATAC)[2] == 1){
      predScores <- t(predScores)
    }
    
    colnames(predScores) <- colnames(regionATAC)
    
  }else{
    predScores = NULL
  }
  
  list("position" = relativePos,
       "region" = region,
       "sites" = sites,
       "predScores" = predScores)
  
}

# Calculate prediction scores for all candidate sites in all regions of interest
setGeneric("getPrediction",
           function(project, # footprintingProject object
                    motifs = NULL, # PWMatrixList object of TF motifs. 
                    # If NULL, the region will be divided into equally sized tiles for TF binding scoring
                    contextRadius = 100, # Local radius of model input (in bp)
                    chunkSize = 2000, # Number of regions per chunk
                    chunkInds = NULL, # Whether to only run the model on a specific subset of chunks. If so, provide chunk indices here
                    innerChunkSize = 10, # Number of regions per inner chunk
                    nCores = 16 # Number of cores to use
           ) standardGeneric("getPrediction"))

setMethod("getPrediction", 
          signature = c(project = "footprintingProject"),
          function(project,
                   motifs = NULL,
                   contextRadius = 100,
                   chunkSize = 2000,
                   chunkInds = NULL,
                   innerChunkSize = 10,
                   nCores = 16) {
            
            # Directory for storing intermediate results
            tmpDir <- dataDir(project)
            
            # Determine chunk size
            if(is.null(chunkSize)){
              chunkSize <- regionChunkSize(project)
            }
            print(paste0("Using chunk size = ", chunkSize))
            
            if(length(countTensor(project)) != 0){
              chunkSize = min(chunkSize, length(countTensor(project)))
            }
            
            # Create a folder for saving intermediate results
            chunkTmpDir <- paste(tmpDir, "chunkedPredResults/", sep = "")
            if(!dir.exists(chunkTmpDir)){
              system(paste("mkdir -p", chunkTmpDir))
            }
            
            # Get region data we need to use later
            width <- regionWidth(project)
            seqBias <- regionBias(project)
            regions <- regionRanges(project)
            model <- predictionModel(project) 
            scales <- as.character(model$scales)
            dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
            names(dispModels) <- scales
            
            # Find motif matches across all regions
            if(!is.null(motifs)){
              motifPositions <- getMotifPositions(project, motifs)
            }else{
              motifPositions <- NULL
            }
            
            # To reduce memory usage, we chunk the region list in to smaller chunks
            cat("Chunking data ..\n")
            groupIDs <- mixedsort(groups(project))
            chunkIntervals <- getChunkInterval(regionRanges(project), chunkSize = chunkSize)
            starts <- chunkIntervals[["starts"]]
            ends <- chunkIntervals[["ends"]]
            
            # See whether the user has specified that we only run on selected chunks 
            if(is.null(chunkInds)){chunkInds <- 1:length(starts)}
            
            # Process each chunk
            for(i in chunkInds){
              
              gc()
              
              # Select regions in the current chunk
              print(paste0("Processing region chunk ", i, " out of ", length(starts), " chunks"))
              print(Sys.time())
              
              # Get ATAC insertion data for the current chunk
              chunkRegions <- starts[i]:ends[i]
              if(length(countTensor(project)) == 0){
                chunkTensorDir <- paste0(tmpDir, "chunkedCountTensor/")
                chunkCountTensor <- readRDS(paste(chunkTensorDir, "chunk_",i, ".rds", sep = ""))
              }else{
                chunkCountTensor <- countTensor(project)[chunkRegions]
              }
              names(chunkCountTensor) <- chunkRegions
              
              # Skip current chunk if result already exists
              if(file.exists(paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))){
                next
              }
              
              print(Sys.time(), "\n")
              
              # The outer for loop iterates through each chunk. Within each iteration,
              # we split the data into even smaller chunks so we can process them in parallel
              # Why do we need inner chunks? Main reason: parallel computation speeds 
              # things up when computing time for each worker is significantly larger than the time
              # required for transferring the data to each worker. We want each worker to do more
              # computation in each iteration and fewer iterations in total
              if(length(chunkRegions) <= 1){stop("Must have more than one regions")}
              if(innerChunkSize >= length(chunkRegions)){innerChunkSize <- as.integer(length(chunkRegions) / 2)}
              innerChunkIntervals <- getChunkInterval(chunkRegions, chunkSize = innerChunkSize)
              innerChunkStarts <- chunkRegions[innerChunkIntervals[["starts"]]]
              innerChunkEnds <- chunkRegions[innerChunkIntervals[["ends"]]]
              
              chunkFootprintResults <- Reduce(c,pbmcapply::pbmclapply(
                1:length(innerChunkStarts),
                function(innerChunkInd){
                  lapply(
                    innerChunkStarts[innerChunkInd] : innerChunkEnds[innerChunkInd],
                    function(regionInd){
                      
                      # Get pseudobulk-by-position matrix of ATAC data for the current region
                      regionATAC <- getRegionATAC(chunkCountTensor, as.character(regionInd), groupIDs, width[regionInd])
                      
                      # Get sites for TF binding scoring
                      if(!is.null(motifs)){
                        sites <- subsetByOverlaps(motifPositions, regions[regionInd])
                      }else{
                        sites <- NULL
                      }
                      
                      # Calculate footprint scores for each region of the current batch
                      getRegionPred(regionATAC = regionATAC,
                                    Tn5Bias = seqBias[regionInd,],
                                    region = regions[regionInd],
                                    dispModels = dispModels,
                                    model = model,
                                    sites = sites)
                      
                    }
                  )
                },
                mc.cores = nCores
              ))
              
              # Save results
              saveRDS(chunkFootprintResults,
                      paste(chunkTmpDir, "chunk_",i, ".rds", sep = ""))
            }
            
            project
          })

# Get the SummarizedExperiment object of prediction scores
getPredSE <- function(project, # footprintingProject object
                      selectedTFs = NULL, # Whether to only retrieve data for selected TFs
                      nCores = 16 # Number of cores to use
){
  
  predDir <- paste0(dataDir(project), "chunkedPredResults/")
  predChunkFiles <- gtools::mixedsort(list.files(predDir))
  predRanges <- NULL
  predScores <- NULL
  regions <- regionRanges(project)
  if(length(regions) == 0){stop("Must load region ranges first!")}
  
  # Get the genomic ranges of prediction sites
  predRangesDf <- data.table::rbindlist(
    pbmcapply::pbmclapply(
      predChunkFiles,
      function(file){
        
        # Load a new chunk of data
        predChunkData <- readRDS(paste0(predDir, file))
        
        # Get TF binding site genomic ranges of the current chunk
        data.table::rbindlist(lapply(
          predChunkData, 
          function(x){
            sites <- x$sites
            if(length(sites) > 0){
              df <- data.frame(chr = as.character(seqnames(sites)),
                               start = start(sites),
                               end = end(sites),
                               TF = sites$TF,
                               regionInd = findOverlaps(x$region, regions, type = "equal")@to)
              if(!is.null(selectedTFs)){
                df <- df[df$TF %in% selectedTFs,]
              }
            }else{
              df <- NULL
            }
            df
          }))
        
      },
      mc.cores = 16
    )
  )
  predRanges <- GRanges(seqnames = as.character(predRangesDf$chr),
                        ranges = IRanges(start = predRangesDf$start,
                                         end = predRangesDf$end))
  predRanges$regionInd <- predRangesDf$regionInd
  predRanges$TF <- predRangesDf$TF
  
  # Get pred-by-pseudobulk matrix of TF binding scores
  predScores <- data.table::rbindlist(
    pbmcapply::pbmclapply(
      predChunkFiles,
      function(file){
        
        # Load a new chunk of data
        predChunkData <- readRDS(paste0(predDir, file))
        
        # Get TF binding site genomic ranges of the current chunk
        data.table::rbindlist(lapply(
          predChunkData, 
          function(x){
            sites <- x$sites
            if(length(sites) > 0){
              df <- as.data.frame(x$predScores)
              if(!is.null(selectedTFs)){
                df <- df[sites$TF %in% selectedTFs,]
              }
            }else{
              df <- NULL
            }
            df
          }))
        
      },
      mc.cores = 16
    )
  )
  
  predScores <- as.matrix(predScores)
  TFBindingSE <- SummarizedExperiment(assay = list(predScores), rowRanges = predRanges)
  
  if(length(groups(project)) == dim(TFBindingSE)[2]){
    colnames(TFBindingSE) <- groups(project)
  }
  
  TFBindingSE
}

# Get the model prediction score matrix for a specific region
getRegionPredMatrix <- function(project, # footprintingProject object
                                regionInd, # Index of the region
                                #Whether to return model predicted value or confidence in the prediction
                                groupIDs = NULL # Whether to use only selected pseudobulks. If so, provide pseudobulk indices here
){
  
  if(is.null(groupIDs)){
    groupIDs <- groups(project)
    if(length(groupIDs) == 0){
      stop("groups(project) slot needs to be set first!")
    }
  }
  
  ret <- getCountData(project, regionInd)
  countData <- ret[["countData"]]
  adjustedRegionInd <- ret[["regionInd"]]
  
  # Get position-by-pseudobulk ATAC insertion matrix
  regionATAC <- getRegionATAC(countData, adjustedRegionInd, groupIDs, regionWidth(project)[regionInd])
  
  model <- predictionModel(project)
  scales <- model$scales # Scales are footprint radii we use for multi-scale footprinting
  dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
  names(dispModels) <- scales
  region <- regionRanges(project)[regionInd]
  
  # Calculate footprint scores for each region of the current batch
  regionPred <- getRegionPred(regionATAC = regionATAC,
                              Tn5Bias = regionBias(project)[regionInd,],
                              region = region,
                              dispModels = dispModels,
                              model = predictionModel(project),
                              tileSize = 1)
  
  featureMatrix <- array(-Inf, dim = dim(regionATAC))
  rownames(featureMatrix) <- rep("", dim(featureMatrix)[1])
  for(ind in 1:length(regionPred$position)){
    siteStart <- start(regionPred$sites[ind]) - start(region) + 1
    siteEnd <- end(regionPred$sites[ind]) - start(region) + 1
    for(pos in siteStart:siteEnd){
      featureMatrix[pos, ] <- pmax(featureMatrix[pos, ],regionPred$predScores[ind,])
    }
  }
  
  featureMatrix
}

# Load TF binding / habitation prediction model
loadModel <- function(
    h5Path
){
  
  library(tensorflow)
  library(keras)
  
  # Load feature scaling factors
  h5file <- H5File$new(h5Path, mode="r")
  footprintMean <- h5file[["footprint_mean"]]
  footprintMean <- footprintMean[1:footprintMean$dims]
  footprintSd <- h5file[["footprint_sd"]]
  footprintSd <- footprintSd[1:footprintSd$dims]
  scales <- h5file[["scales"]]
  scales <- scales[1:scales$dims]
  h5file$close_all()
  
  # Create an empty model as placeholder
  model <- keras_model_sequential() %>%
    layer_dense(64, activation = 'relu', input_shape = shape(201 * length(scales))) %>%
    layer_dense(16, activation = 'relu') %>%
    layer_dense(1)
  
  # Load pre-trained model weights from h5 file
  load_model_weights_hdf5(
    model,
    filepath = h5Path
  )
  model <- keras::get_weights(model)
  
  model$footprintMean <- footprintMean
  model$footprintSd <- footprintSd
  model$scales <- c(10,20,30,50,80)
  
  # Return the retrieved values
  model
  
}

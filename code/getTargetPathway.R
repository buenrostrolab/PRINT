# Calculate pathway enrichment
pathwayEnrichment <- function(fgGenes, # Foreground gene symbols
                              bgGenes, # Background gene symbols, must contain all genes in the foreground
                              geneSets, # Gene sets. A named list of genes in each pathway
                              pvalThrshold = 0.01, # Threshold to filter results
                              nGeneThreshold = 10 # Remove pathways with foreground genes fewer than this number
                              
){
  
  # Remove repetitive entries from gene lists
  fgGenes <- unique(fgGenes)
  bgGenes <- unique(bgGenes)
  
  # Get gene-to-pathway/annotation mapping 
  geneSetMatrix <- sapply(
    geneSets,
    function(geneSet){
      bgGenes %in% geneSet
    }
  )
  rownames(geneSetMatrix) <- bgGenes
  colnames(geneSetMatrix) <- names(geneSets)
  
  # Subset the gene-to-pathway matrix by keeping only background genes
  geneSetMatrix <- geneSetMatrix[bgGenes, ]
  
  isForeground <- (bgGenes %in% fgGenes)
  contingencyMat <- Reduce(rbind, list(
    isForeground %*% geneSetMatrix,
    isForeground %*% (!geneSetMatrix),
    (!isForeground) %*% geneSetMatrix,
    (!isForeground) %*% (!geneSetMatrix)
  ))
  
  # Calculate pvalues of enrichment using Fisher's exact test
  pvals <- apply(
    contingencyMat, 2,
    function(x){
      fisher.test(array(x , dim = c(2,2)))$p.value 
    }
  )
  
  # Calculate odds ratio
  oddsRatio <- (contingencyMat[1,] / contingencyMat[2, ]) / (contingencyMat[3, ] / contingencyMat[4, ])
  
  # Find target genes
  mappedTargets <- apply(
    geneSetMatrix, 2,
    function(x){
      paste(sort(bgGenes[x & isForeground]), collapse = ",")
    }
  )
  
  # Combine results
  enrichment <- data.frame(
    pathway = colnames(geneSetMatrix),
    pvals = pvals,
    fdrs = p.adjust(pvals, method = "fdr"),
    logOR = log2(oddsRatio),
    nDiff = contingencyMat[1, ],
    nExpected = colSums(geneSetMatrix) * sum(isForeground) / dim(geneSetMatrix)[1],
    mappedTargets = mappedTargets
  )
  enrichment <- enrichment[(enrichment$logOR > 0) & (enrichment$nDiff > nGeneThreshold), ]
  enrichment <- enrichment[(enrichment$pval < pvalThrshold), ]
  enrichment <- enrichment[order(enrichment$pval),]
  enrichment
}

# Mapping TFs to target pathways
mapTargets <- function(
    TFMotifSites, # A named list of motif sites. Each list element is a GRanges object.
    siteSE, # A RangedSummarizedExperiment object where rows are genomic regions and 
    # columns are samples (cells, pseudobulks, etc). Values are any type of signal (TF binding score, accessibility, etc.)
    geneSets, # A named list. E.g., list("pathway1" = c("gene1", "gene2"), "pathway2" = ...)
    genome, # Reference genome to use
    nTargets = NULL, # Number of target sites to keep. Overrides the percentile argument.
    threshold = 0.9, # If nTargets is not set, keep top binding sites passing this percentile (0.9 means 90% percentile)
    nMotifThreshold = 100, # Only keep TFs with motif matches larger than this number
    mode = "max",  # How to calculate signal of a region. One of "mean", "min", "max"
    proximalOnly = T, # Whether to include only promoter regions
    regionGeneCorr = NULL, # Must have the following columns: Gene and Region (e.g., "chr1:100-200")
    # If proximalOnly is set to FALSE, this argument must be provided
    returnTargets = F, # Whether to return mapped target genes for each TF
    nCores = 4, # Number of cores to use
    rankByMotifScores = F # Whether to rank motif sites by motif match scores instead of scores from siteSE
){
  
  # Remove TFs with low number of motif sites
  nMotifs <- sapply(TFMotifSites, length)
  TFMotifSites <- TFMotifSites[nMotifs > nMotifThreshold]
  TFs <- names(TFMotifSites)
  
  # Retrieve signal (or binding score) from siteSE
  siteRangesAll <- rowRanges(siteSE)
  siteScoresAll <- assay(siteSE)
  
  # For each genomic location, obtain a single score by summarizing across samples (columns of siteSE)
  if(mode == "mean"){
    siteSignal <- rowMeans(siteScoresAll)
  }else if(mode == "min"){
    siteSignal <- rowMins(siteScoresAll)
  }else if(mode == "max"){
    siteSignal <- rowMaxs(siteScoresAll)
  }
  
  # Get TSS coordinates
  if(genome == "hg38"){
    TSS <- FigR::hg38TSSRanges
  }else if(genome == "mm10"){
    TSS <- FigR::mm10TSSRanges
  }
  allGenes <- as.character(TSS$gene_name)
  
  # For each TF, calculate pathway enrichment among target genes
  enrichmentList <- pbmcapply::pbmclapply(
    names(TFMotifSites),
    function(TF){
      
      #################### Find binding sites of a particular TF ####################
      
      # Find sites with a matched motif
      sites <- TFMotifSites[[TF]]
      
      # Randomly shuffle sites. This is for the case where there are a lot of tied values in binding scores
      set.seed(42)
      sites <- sites[sample(1:length(sites))]
      
      # For each TF motif site, calculate its binding score
      ov <- findOverlaps(siteRangesAll, sites)
      ovMat <- Matrix::sparseMatrix(
        i = ov@from,
        j = ov@to,
        x = 1, dims = c(length(siteRangesAll), length(sites))
      )
      siteScores <- t(siteScoresAll) %*% ovMat
      siteScores <- siteScores[1, ]
      
      if(rankByMotifScores){
        siteScores <- sites$score
      }
      
      if(is.null(nTargets)){
        targetSites <- sites[siteScores >= quantile(siteScores, threshold)]
      }else{
        targetSites <- sites[order(-siteScores)][1:min(nTargets, length(sites))]
      }
      
      #################### Map TF bound regions to target genes ####################
      
      # Map the above bound regions to target genes by overlapping with promoters
      names(TSS) <- as.character(TSS$gene_name)
      TSS <- resize(TSS, 1000, fix = "center")
      proximalTargets <- sort(unique(names(subsetByOverlaps(TSS, targetSites))))
      
      # Merge the above two target gene lists
      if(proximalOnly){
        targetGenes <- sort(unique(proximalTargets))
      }else{
        
        # If we want to include distal targets (TF-to-enhancer-to-gene mapping), we need to make sure
        # the cCRE-gene correlation table is provided for enhancer-to-gene mapping
        if(is.null(regionGeneCorr)){
          stop("If proximalOnly is set to FALSE, then regionGeneCorr must be provided")
        }
        stopifnot(c("Gene", "Region") %in% colnames(regionGeneCorr))
        
        # Map the above bound regions to target genes by overlapping with enhancers
        ov <- findOverlaps(GRanges(regionGeneCorr$Region), targetSites)
        distalTargets <- unique(regionGeneCorr[ov@from, ]$Gene)
        targetGenes <- sort(unique(c(proximalTargets, distalTargets)))
      }
      
      # Also obtain background gene list
      bgGenes <- sort(unique(names(TSS)))
      
      # Calculate enrichment of targets in pathways
      enrichment <- pathwayEnrichment(targetGenes, bgGenes, geneSets, 
                                      pvalThrshold = 1, nGeneThreshold = 0)
      if(nrow(enrichment) > 0){enrichment$TF <- TF}
      
      # Append to list
      if(returnTargets){
        list(enrichment = enrichment, targets = targetGenes)
      }else{
        enrichment
      }
      
    },
    mc.cores = nCores
  )
  
  names(enrichmentList) <- TFs
  enrichmentList
  
}

# Mapping TFs to target pathways
mapTargetsGSEA <- function(
    TFMotifSites, # A named list of motif sites. Each list element is a GRanges object.
    siteSE, # A RangedSummarizedExperiment object where rows are genomic ranges of candidate binding sites and 
    # columns are samples (cells, pseudobulks, etc). Values are any type of signal (TF binding score, accessibility, etc.)
    regionGeneCorr, # Must have the following columns: Gene and Region (e.g., "chr1:100-200")
    geneSets, # A named list. E.g., list("pathway1" = c("gene1", "gene2"), "pathway2" = ...)
    nMotifThreshold = 100, # Only keep TFs with motif matches larger than this number
    # If False, assign a numeric score to each gene and compare genes in and outside of a gene set.
    genome = "hg38",
    nPermute = 500, # Number of permutations
    proximalOnly = F # Whether to include only promoter sites
){
  
  stopifnot(c("Gene", "Region") %in% colnames(regionGeneCorr))
  
  nMotifs <- sapply(TFMotifSites, length)
  TFMotifSites <- TFMotifSites[nMotifs > nMotifThreshold]
  TFs <- names(TFMotifSites)
  
  # Get genomic ranges and scores of candidate binding sites
  siteRangesAll <- rowRanges(siteSE)
  siteScoresAll <- assay(siteSE)
  
  # Get all regions linked to genes
  regions <- GRanges(regionGeneCorr$Region)
  
  # Assign a score to each candidate site
  siteSignal <- rowMaxs(siteScoresAll)
  
  # Get TSSs
  if(genome == "hg38"){
    TSS <- FigR::hg38TSSRanges
  }else if(genome == "mm10"){
    TSS <- FigR::mm10TSSRanges
  }
  TSS <- resize(TSS, 1000, fix = "center")
  
  # Get a list of all genes
  allGenes <- as.character(TSS$gene_name)
  
  # Construct a gene-to-site matrix
  proximalOv <- findOverlaps(TSS, siteRangesAll)
  proximalMapping <- Matrix::sparseMatrix(
    i = proximalOv@from,
    j = proximalOv@to,
    x = 1, dims = c(length(TSS), length(siteRangesAll))
  )
  
  if(!proximalOnly){
    
    # First get a distal gene-to-region matrix
    distalGeneRegionMapping <- Matrix::sparseMatrix(
      i = match(regionGeneCorr$Gene, allGenes),
      j = match(regionGeneCorr$Region, as.character(regions)),
      x = 1, dims = c(length(allGenes), length(regions))
    )
    
    # Then get a distal region-to-site matrix
    distalOv <- findOverlaps(regions, siteRangesAll) 
    distalRegionSiteMapping <- Matrix::sparseMatrix(
      i = distalOv@from,
      j = distalOv@to,
      x = 1, dims = c(length(regions), length(siteRangesAll))
    )
    
    # Combine the two matrices to get a gene-to-site matrix
    distalMapping <- distalGeneRegionMapping %*% distalRegionSiteMapping
    
    # Combine proximal and distal mapping
    combinedMapping <- proximalMapping + distalMapping
    
  }else{
    
    combinedMapping <- proximalMapping
  }
  rownames(combinedMapping) <- allGenes
  
  # Binarize the mapping
  combinedMapping <- 1 * (combinedMapping > 0)
  
  # Assign siteSignal to each column (site)
  combinedMapping <- t(t(combinedMapping) * siteSignal)
  
  # Construct a geneSet-to-gene matrix
  geneSetNames <- names(geneSets)
  geneSetMatrixIJV <- data.table::rbindlist(
    pbmcapply::pbmclapply(
      geneSetNames,
      function(geneSet){
        data.frame(i = match(geneSet, geneSetNames), j = match(geneSets[[geneSet]], allGenes), x = 1)
      },
      mc.cores = 16
    )
  )
  geneSetMatrixIJV <- geneSetMatrixIJV[!is.na(geneSetMatrixIJV$j),]
  geneSetMatrix <- Matrix::sparseMatrix(
    i = geneSetMatrixIJV$i, 
    j = geneSetMatrixIJV$j,
    x = 1, dims = c(length(geneSets), length(allGenes))
  )
  
  # Generate random permutations of the above matrix but keeping the number of genes per geneSet the same
  permutedGeneSetMatrices <- lapply(
    1:nPermute,
    function(i){
      geneSetMatrix[, sample(1:ncol(geneSetMatrix))]
    }
  )
  
  # For each TF, calculate pathway enrichment among target genes
  resultList <- pbmcapply::pbmclapply(
    names(TFMotifSites),
    function(TF){
      
      # Find sites with a matched motif
      motifSites <- TFMotifSites[[TF]]
      motifSiteOv <- findOverlaps(motifSites, siteRangesAll)
      
      # For each gene, summarize the scores of all associated sites to a single score
      # This results in a score vector of length number_of_genes
      motifSiteScoreMatrix <- combinedMapping[, motifSiteOv@to]
      geneScores <- as.matrix(rowMaxs(motifSiteScoreMatrix))
      names(geneScores) <- allGenes
      
      # Calculate scores for the permuted background
      geneSetScores <- geneSetMatrix %*% geneScores
      permutedGeneSetScores <- sapply(
        permutedGeneSetMatrices,
        function(permutation){
          as.matrix(permutation %*% geneScores)
        }
      )
      
      # Get z-scores 
      zScores <- as.numeric((geneSetScores - rowMeans(permutedGeneSetScores)) / rowSds(permutedGeneSetScores))
      names(zScores) <- names(geneSets)
      
      # Convert z-scores to p-values
      pvals <- pnorm(zScores, lower.tail = F)
      
      # Create data frame for results
      results <- data.frame(
        pathway = names(geneSets),
        pvals = pvals,
        zscores = zScores,
        TF = TF
      )
      rownames(results) <- geneSetNames
      
      # Get FDRs
      results$fdr <- p.adjust(results$pvals, method = "fdr")
      
      results <- results[order(results$pvals),]
      results
    },
    mc.cores = 16
  )
  
  names(resultList) <- TFs
  resultList
  
}
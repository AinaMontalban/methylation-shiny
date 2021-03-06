#directory <- "/mnt/ElRaid/amontalban/PROJECTS/methylation/aina/shiny"
#source(paste0(directory, "/", "shinyMethylSet.R"))

setGeneric("shinySummarizepr",
           function(object, ...) standardGeneric("shinySummarizepr"))

setGeneric("shinySummarizeNorm",
           function(object, ...) standardGeneric("shinySummarizeNorm"))
setMethod("shinySummarizeNorm", signature(object = "GenomicRatioSet"),
          function(object) {
            object <- updateObject(object)
            .createIndices <- function(object, betaMatrix) {
              ann  <- getAnnotation(object)
              chr  <- ann$chr
              probeType <- paste0(ann$Type,ann$Color)
              probeNames <- rownames(ann)
              autosomal <- (chr %in% paste0("chr", 1:22))
              indices <- list(IGrn = probeNames[probeType == "IGrn" & autosomal],
                              IRed = probeNames[probeType == "IRed" & autosomal],
                              II = probeNames[probeType == "II" & autosomal],
                              X = probeNames[chr == "chrX"],
                              Y = probeNames[chr == "chrY"])
              for (i in 1:length(indices)){
                indices[[i]] <- which(rownames(betaMatrix) %in% indices[[i]])
              }
              indices
            }
            
            
            cat("[shinySummarizeNorm] Computing methylation values beta-values \n")
            betaMatrix <- minfi::getBeta(object)
            mMatrix    <- minfi::getM(object)
            cnMatrix   <- minfi::getCN(object)
            
            probe.indices <- .createIndices(object, betaMatrix)
            autosomal <- unlist(probe.indices[1:3])
            probs <- seq(from = 0, to = 1, length.out = 500)
            
            betaQuantiles <- vector("list", length(probe.indices)) 
            names(betaQuantiles) <- names(probe.indices)
            mQuantiles <- cnQuantiles <-  betaQuantiles 
            
            cat("[shinySummarizeNorm] Computing quantiles \n")   
            ## To compute the quantiles: 
            for (i in 1:length(probe.indices)){
              betaQuantiles[[i]] <- t(colQuantiles(betaMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              mQuantiles[[i]] <- t(colQuantiles(mMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              cnQuantiles[[i]] <- t(colQuantiles(cnMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              names(betaQuantiles[[i]]) <- names(mQuantiles[[i]]) <- names(cnQuantiles[[i]]) <- colnames(object)  
            }
            
            cat("[shinySummarizeNorm] Computing principal components beta-values \n")
            ## To compute the principal components:
            numPositions = 20000
            autMatrix <- betaMatrix[autosomal,]
            rm(betaMatrix)
            o <- order(-rowVars(autMatrix))[1:numPositions]
            pca <- prcomp(t(autMatrix[o,]))
            pca <- list(scores = pca$x, percs = (pca$sdev^2)/(sum(pca$sdev^2))*100)
            names(pca$percs) <- colnames(object)
            
            cat("[shinySummarizeNorm] Computing principal components m-values \n")
            ## To compute the principal components:
            numPositions = 20000
            autMatrix <- mMatrix[autosomal,]
            # rm(mMatrix)
            gc(verbose=FALSE)
            o <- order(-rowVars(autMatrix))[1:numPositions]
            pca_m <- prcomp(t(autMatrix[o,]))
            pca_m = list(scores = pca_m$x, percs = (pca_m$sdev^2)/(sum(pca_m$sdev^2))*100)
            names(pca_m$percs) <- colnames(object)
            GRSet <- object
            
            object <- shinyMethylSet(sampleNames = colnames(object),
                                     phenotype   = as.data.frame(minfi::pData(object)),
                                     mQuantiles  = mQuantiles,
                                     betaQuantiles = betaQuantiles,
                                     betaMatrix = minfi::getBeta(object),
                                     mMatrix = minfi::getM(object),
                                     
                                     methQuantiles = list(NULL),
                                     unmethQuantiles = list(NULL),
                                     cnQuantiles = cnQuantiles,
                                     greenControls = list(NULL),
                                     redControls = list(NULL),
                                     pca = pca,
                                     pca_m = pca_m,
                                     
                                     detP = matrix(NA),
                                     originObject = "GenomicRatioSet",
                                     array = object@annotation[["array"]],
                                     snps = matrix(NA),
                                     GRSet = GRSet
            )
            object
          })
setMethod("shinySummarizepr", signature(object = "RGChannelSet"),
          function(object) {
            .createIndices <- function(object, betaMatrix) {
              ann  <- getAnnotation(object)
              chr  <- ann$chr
              probeType <- paste0(ann$Type,ann$Color)
              probeNames <- rownames(ann)
              autosomal <- (chr %in% paste0("chr", 1:22))
              indices <- list(IGrn = probeNames[probeType == "IGrn" & autosomal],
                              IRed = probeNames[probeType == "IRed" & autosomal],
                              II = probeNames[probeType == "II" & autosomal],
                              X = probeNames[chr == "chrX"],
                              Y = probeNames[chr == "chrY"])
              for (i in 1:length(indices)){
                indices[[i]] <- which(rownames(betaMatrix) %in% indices[[i]])
              }
              indices
            }
            
            controlType <- c("BISULFITE CONVERSION I",
                             "BISULFITE CONVERSION II",
                             "EXTENSION",
                             "HYBRIDIZATION",
                             "NEGATIVE",
                             "NON-POLYMORPHIC",
                             "NORM_A",
                             "NORM_C",
                             "NORM_G",
                             "NORM_T", 
                             "SPECIFICITY I",
                             "SPECIFICITY II",
                             "TARGET REMOVAL",
                             "STAINING")
            
            object <- updateObject(object)
            cat("[shinySummarize] Extracting Red and Green channels \n")
            r <- minfi::getRed(object)
            g <- minfi::getGreen(object)
            
            ## Extraction of the controls
            greenControls <- redControls <- vector("list",length(controlType))
            names(greenControls) <- names(redControls) <- controlType
            
            for (i in 1:length(controlType)){
              ctrlAddress <- minfi::getControlAddress(
                object, controlType = controlType[i])
              
              redControls[[i]]=r[ctrlAddress,]
              greenControls[[i]]=g[ctrlAddress,]
            }
            
            rm(r)
            rm(g)
            gc(verbose=FALSE)
            
            cat("[shinySummarize] Raw preprocessing \n")
            mset <- minfi::preprocessRaw(object)
            methMatrix <- minfi::getMeth(mset) 
            unmethMatrix <- minfi::getUnmeth(mset)
            rm(mset)
            betaMatrix <- methMatrix/(methMatrix+unmethMatrix+100)
            mMatrix    <- log2((methMatrix+1)/(unmethMatrix+1))
            cnMatrix   <- log2(methMatrix+unmethMatrix)
            
            cat("[shinySummarize] Calculating detection p-values \n")
            targets <- pData(object)
            snames <- targets$Sample_Name
            detP <- detectionP(object)
            colnames(detP) <- snames
            
            cat("[shinySummarize] Mapping to genome \n")
            gmSet <- minfi::mapToGenome(object)
            
            probe.indices <- .createIndices(gmSet, betaMatrix)
            autosomal <- unlist(probe.indices[1:3])
            
            rm(gmSet)
            gc(verbose=FALSE)
            probs <- seq(from = 0, to = 1, length.out = 500)
            betaQuantiles <- vector("list", length(probe.indices)) 
            names(betaQuantiles) <- names(probe.indices)
            mQuantiles <- methQuantiles <- unmethQuantiles <- cnQuantiles <- betaQuantiles
            
            cat("[shinySummarize] Computing quantiles \n")   
            ## To compute the quantiles: 
            for (i in 1:length(probe.indices)){
              betaQuantiles[[i]] <- t(colQuantiles(betaMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              mQuantiles[[i]] <- t(colQuantiles(mMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              methQuantiles[[i]] <- t(colQuantiles(methMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              unmethQuantiles[[i]] <- t(colQuantiles(unmethMatrix[probe.indices[[i]],], probs=probs,na.rm=TRUE))
              cnQuantiles[[i]] <- t(colQuantiles(cnMatrix[probe.indices[[i]],], probs=probs, na.rm=TRUE))
              
              colnames(betaQuantiles[[i]]) <- colnames(mQuantiles[[i]])      <- 
                colnames(methQuantiles[[i]])   <-
                colnames(unmethQuantiles[[i]]) <-
                colnames(cnQuantiles[[i]]) <- colnames(object)
            }
            rm(methMatrix)
            rm(unmethMatrix)
           # rm(mMatrix)
            rm(cnMatrix)
            gc(verbose=FALSE)
            
            cat("[shinySummarize] Computing the SNPs beta values \n")
            ## To compute the principal components:
            snps <- getSnpBeta(object)
            
            
            cat("[shinySummarize] Computing principal components b-values \n")
            ## To compute the principal components:
            numPositions = 20000
            autMatrix <- betaMatrix[autosomal,]
            rm(betaMatrix)
            gc(verbose=FALSE)
            o <- order(-rowVars(autMatrix))[1:numPositions]
            pca <- prcomp(t(autMatrix[o,]))
            pca = list(scores = pca$x, percs = (pca$sdev^2)/(sum(pca$sdev^2))*100)
            names(pca$percs) <- colnames(object)
            
            cat("[shinySummarize] Computing principal components m-values \n")
            ## To compute the principal components:
            numPositions = 20000
            autMatrix <- mMatrix[autosomal,]
           # rm(mMatrix)
            gc(verbose=FALSE)
            o <- order(-rowVars(autMatrix))[1:numPositions]
            pca_m <- prcomp(t(autMatrix[o,]))
            pca_m = list(scores = pca_m$x, percs = (pca_m$sdev^2)/(sum(pca_m$sdev^2))*100)
            names(pca_m$percs) <- colnames(object)
            
         
            cat("[shinySummarize] Normalizing with ssNoob \n")
            # To compute the normalization:
            MSet.noob <- preprocessNoob(object)
            mSetSq <- ratioConvert(MSet.noob)
           #Convert to GenomicRatioSet object.
            cat("[shinySummarize] ------> Convert to GenomicRatioSet object. \n")
            mSetSq <- mapToGenome(mSetSq)
            cat("[shinySummarize] ------> Convert to ShinyMethylSet object. \n")
            norm <- shinySummarizeNorm(mSetSq)
            
            RGSET <- object
            
            object <-  shinyMethylSet(sampleNames = colnames(object),
                                      phenotype   = as.data.frame(minfi::pData(object)),
                                      mQuantiles  = mQuantiles,
                                      betaQuantiles = betaQuantiles,
                                      betaMatrix = minfi::getBeta(object),
                                      mMatrix = mMatrix,
                                      methQuantiles = methQuantiles,
                                      unmethQuantiles = unmethQuantiles,
                                      cnQuantiles = cnQuantiles ,
                                      greenControls = greenControls ,
                                      redControls = redControls ,
                                      pca = pca,
                                      pca_m = pca_m,
                                      detP = detP,
                                      originObject = "RGChannelSet",
                                      array = object@annotation[["array"]],
                                      RGSET = RGSET,
                                      snps = snps,
                                      norm = norm
                                      
            )
            object
          })



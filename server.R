server.methylation <- function(shinyMethylSet1, shinyMethylSet2=NULL){
  function(input, output, session) { 
    betaQuantiles   <-  getBeta(shinyMethylSet1)
    mQuantiles      <-  getM(shinyMethylSet1)
    methQuantiles   <-  getMeth(shinyMethylSet1)
    unmethQuantiles <-  getUnmeth(shinyMethylSet1)
    cnQuantiles     <-  getCN(shinyMethylSet1)
    bMatrix <- shinyMethylSet1@betaMatrix
    mMatrix <- shinyMethylSet1@mMatrix
    greenControls   <-  getGreenControls(shinyMethylSet1)
    redControls     <-  getRedControls(shinyMethylSet1)
    covariates      <<- pData(shinyMethylSet1)
    pca             <-  getPCA(shinyMethylSet1)$scores
    pca_m           <-  shinyMethylSet1@pca_m$scores
    detP            <- shinyMethylSet1@detP
    sampleNames     <-  sampleNames(shinyMethylSet1)
    slideNames      <- shinyMethylSet1@phenotype$Slide
    arrayNames      <- shinyMethylSet1@phenotype$Array
    plateNames      <- shinyMethylSet1@phenotype$Sample_Plate
    groupNames      <- shinyMethylSet1@phenotype$Sample_Group
    controlNames    <-  names(greenControls)
    targets <- shinyMethylSet1@phenotype
    snps <- shinyMethylSet1@snps
    colnames(snps) <- targets$Sample_Name
    method <- shinyMethylSet1@originObject
    RGSET <- shinyMethylSet1@RGSET
    normalized <- shinyMethylSet1@norm
    otherNorm <- FALSE
    ann <- getAnnotation(RGSET)
    
    GRSet <- normalized@GRSet
    
    
    ## In the case covariates is empty:
    if (ncol(covariates)==0){
      covariates <- data.frame(slide = slideNames, plate = plateNames)
      rownames(covariates) <- sampleNames
      covariates <<- covariates
    }
    
    ## Global variables:
    ## that will be accessed by the knitr report
    mouse.click.indices <<- c()
    colorSet     <<- "Set1" # Default color set
    sampleColors <<- as.numeric(as.factor(plateNames)) ## Default sample colors
    genderCutoff <<- -0.4 # Default cutoff for gender prediction
    current.control.type <<- "BISULFITE CONVERSION I"
    current.probe.type   <<- "I Green"
    current.density.type <<- "Beta-value"
    
    #########################################################
    
    ###########        Colors 
    
    #########################################################
    
    
    # To change the global color scheme:
    setColor <- reactive({
      colorSet <<- input$colorChoice
    })
    
    set.palette <- function(n, name){
      ## The name of the colors are part of the RColorBrewer package
      if (name != "Pault" & name != "Rainbow"){
        palette(brewer.pal(n = n, name = name ))
        
        ## Custom palette
      } else if (name == "Pault"){
        colors <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
                    "#DDCC77","#661100","#CC6677","#882255","#AA4499")
        colors <- c("#4477AA","#CC6677","#DDCC77","#117733","#88CCEE",
                    "#AA4499","#44AA99","#999933","#882255","#661100",
                    "#6699CC", "#AA4466")
        palette(colors)
      } else if (name == "Rainbow"){
        colors <- c("#781C81", "#3F56A7","#4B91C0","#5FAA9F","#91BD61",
                    "#D8AF3D","#E77C30","#D92120")
        palette(colors)
      }
    }
    ## To choose the colors according to the phenotype:
    sampleColors <- reactive({
      return(as.numeric(as.factor(covariates[,match(input$phenotype,
                                                    colnames(covariates))])))
    })
    sampleColorsPCA <- reactive({
      return(as.numeric(as.factor(covariates[,match(input$phenoID,
                                                    colnames(covariates))])))
    })
    #########################################################
    ###########        Computation of the densities
    #########################################################
    
    ## Return x and y matrices of the densities for raw data
    returnDensityMatrix <- reactive({
      index <- match(input$probeType,c("I Green","I Red","II","X","Y"))
      if (input$mOrBeta=="Beta-value"){
        bw <-1
        quantiles <- betaQuantiles[[index]]
      } else {
        bw <- 1
        quantiles <- mQuantiles[[index]]
      }
      
      matrix <- apply(quantiles, 2, function(x){
        d <- density(x, bw = bw, n =512)
        c(d$x,d$y)
      })
      matrix.x <- matrix[1:512,]
      matrix.y <- matrix[513:(2*512),]
      return(list(matrix.x = matrix.x, matrix.y = matrix.y))
    })
    
    returnDensityMatrixNorm <- reactive({
        index <- match(input$probeType,
                       c("I Green","I Red","II","X","Y"))
      if (input$mOrBeta=="Beta-value"){
        bw <- 1
        if (input$normButton){
          quantiles <- norm()@betaQuantiles[[index]]
        } else {
          quantiles <- normalized@betaQuantiles[[index]]
        }
      } else {
        bw <- 1
        if (input$normButton){
        quantiles <- norm()@mQuantiles[[index]]
        } else {
          quantiles <- normalized@mQuantiles[[index]]
          
        }
      }
      matrix <- apply(quantiles, 2, function(x){
        d <- density(x, bw = bw, n =512)
        c(d$x,d$y)
      })
      matrix.x <- matrix[1:512,]
      matrix.y <- matrix[513:(2*512),]
      return(list(matrix.x = matrix.x, matrix.y = matrix.y))	
    })
    
    #########################################################
    ###########         Densities plots
    #########################################################
    
    ## Densities plot for raw data  
    output$rawDensities <- renderPlot({
      
      set.palette(n=8, name=setColor())
      colors <- sampleColors()
      lwd <- as.numeric(input$lwd)
      lty <- as.numeric(input$lty)
      index = match(input$probeType,c("I Green","I Red","II","X","Y"))
      
      ## Plot specifications                  
      if (input$mOrBeta=="Beta-value"){
        xlim <- c(-0.2,1.2)
        if (input$probeType == "II"){
          ylim <- c(0,1)
        } else {
          ylim <- c(0,8)
        }
        from = -4; to = 4;
        main = "BETA-VALUE DENSITIES"
        xlab = "Beta-values"
        bw <- 1
        #quantiles <- betaQuantiles[[index]]
        minfi::densityPlot(bMatrix, sampGroups = groupNames)
        
      } else {
        xlim <- c(-8,8)
        ylim <- c(0, 0.35)
        from = -10; to = 10;
        main = "M-VALUE DENSITIES" 
        xlab = "M-values"
        quantiles <- mQuantiles[[index]]
        bw <- 1
        densitiesPlot(matrix.x = returnDensityMatrix()[[1]],
                      matrix.y = returnDensityMatrix()[[2]],
                      quantiles = quantiles,
                      sampleNames = sampleNames,
                      main = main, xlab = xlab,
                      xlim = xlim, ylim = ylim, col = colors,
                      mean = input$mean, 
                      lwd = lwd, lty = lty,
                      from = from, to = to, bw = bw)
        
      }})

    # Density plot for normalized data:
    output$normDensities <- renderPlot({


        set.palette(n=8, name=setColor())
        colors <- sampleColors()
        lwd <- as.numeric(input$lwd)
        lty <- as.numeric(input$lty)
        index <- match(input$probeType,
                       c("I Green","I Red","II","X","Y"))

        ## Plot specifications
        if (input$mOrBeta=="Beta-value"){
          xlim <- c(-0.2,1.2)
          if (input$probeType == "II"){
            ylim <- c(0,6)
          } else {
            ylim <- c(0,10)
          }
          from = -4; to = 4;
          main = "Normalized data (Beta values)"
          xlab = "Beta-values"
          bw <- 1
          if (input$normButton){
            minfi::densityPlot(norm()@betaMatrix, sampGroups = groupNames)
          } else {
            minfi::densityPlot(normalized@betaMatrix, sampGroups = groupNames)
          }
          #quantiles =  shinyMethylSet2@betaQuantiles[[index]]
        } else {
          xlim <- c(-8,8)
          ylim <- c(0, 0.35)
          from = -10; to = 10;
          main = "Normalized data (M values)"
          xlab = "M-values"
          otherNorm
          if (input$normButton){
          quantiles =  norm()@mQuantiles[[index]]
          } else {
            quantiles =  normalized@mQuantiles[[index]]
          }
          bw <- 1


          densitiesPlot(matrix.x = returnDensityMatrixNorm()[[1]],
                        matrix.y = returnDensityMatrixNorm()[[2]],
                        quantiles = quantiles,
                        sampleNames = sampleNames,
                        main = main, xlab = xlab,
                        xlim = xlim, ylim = ylim, col = colors,
                        mean = input$mean,
                        lwd = lwd, lty = lty,
                        from = from, to = to, bw=bw)}

      })
    
    output$click_info <- renderPrint({
      cat("str(input$click_action):")
      str(input$click_action)
    })
    
    
    output$betamatrixDownload <- downloadHandler(
      
      filename = function() {
        paste("raw_bvalues", ".csv", sep = "")
      },
      content = function(file) {
        bMatrix <- shinyMethylSet1@betaMatrix
        write.csv(bMatrix, file)
      }
    )
    
    output$betamatrixNormDownload <- downloadHandler(
      
      filename = function() {
        paste("norm_bvalues", ".csv", sep = "")
      },
      content = function(file) {
        if (input$normButton){
        bMatrix <- norm()@betaMatrix
        } else {
          bMatrix <- normalized@betaMatrix
        }
        write.csv(bMatrix, file)
      }
    )

    output$mMatrixDownload <- downloadHandler(
      
      filename = function() {
        paste("raw_bvalues", ".csv", sep = "")
      },
      content = function(file) {
        mMatrix <- shinyMethylSet1@mMatrix
        write.csv(mMatrix, file)
      }
    )
    output$mMatrixNormDownload <- downloadHandler(
      
      filename = function() {
        paste("norm_bvalues", ".csv", sep = "")
      },
      content = function(file) {
        if (input$normButton){
          mMatrix <- norm()@mMatrix
        } else {
          mMatrix <- normalized@mMatrix
        }
        write.csv(mMatrix, file)
      }
    )
    
    # output$DrawButton <- downloadHandler({
    #   filename =  "rawDensityPlot_bvals.png"
    #   content = function(file){
    #     png(file)
    #     minfi::densityPlot(bMatrix, sampGroups = groupNames)
    #     dev.off()
    #   }
    #   
    # }
    # )

    output$controlTypePlotGreen <- renderPlot({
      cnt <- input$controlType
      
      if (input$controlType %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I", 
                                   "SPECIFICITY II", "TARGET REMOVAL")){
        threshold <- 1
      } else if (input$controlType %in% c("EXTENSION", "STAINING", "NON-POLYMORPHIC")){
        threshold <- 5 # you can increase the threshold
      } else {threshold <- 0}
      
      if (length(sampleNames) >= 50){
        arr <- input$arrayID
        column_names <- colnames(greenControls[[cnt]])
        idx_col <- grep(arr, column_names)
        subset <- greenControls[[cnt]][, idx_col]
        title <- paste("-", arr)
        # array_names <- substr(colnames(subset), nchar(slideNames[1]) + 2, nchar(slideNames[1]) + 7)
        # colnames(subset) <- array_names
      }  else {
        subset <- greenControls[[cnt]]
        title <- ""
      }
      
      log2_subset_GC <- log2(subset)
      df_subset_GC <- melt(log2_subset_GC)
      ggplot(data=as.data.frame(df_subset_GC), aes(x=Var2, y=value)) + 
        geom_point(color="darkgreen", size=1.5) + scale_y_continuous(limits = c(-1, 20)) + 
        theme(axis.text.x = element_text(hjust = 1, angle=45)) +
        geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") + 
        scale_x_discrete(labels=groupNames) + xlab("Samples") + ggtitle(paste("Green Channel", title))
    }
    )
    
    output$controlTypePlotRed <- renderPlot({
      cnt <- input$controlType
      
      if (input$controlType %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I", 
                                   "SPECIFICITY II", "TARGET REMOVAL")){
        threshold <- 1
      } else if (input$controlType %in% c("EXTENSION", "STAINING", "NON-POLYMORPHIC")){
        threshold <- 5 # you can increase the threshold
      } else {threshold <- 0}
      
      if (length(sampleNames) >= 50){
        arr <- input$arrayID
        column_names <- colnames(redControls[[cnt]])
        idx_col <- grep(arr, column_names)
        subset <- redControls[[cnt]][, idx_col]
        title <- paste("-", arr)
        # array_names <- substr(colnames(subset), nchar(slideNames[1]) + 2, nchar(slideNames[1]) + 7)
        # colnames(subset) <- array_names
      }  else {
        subset <- redControls[[cnt]]
        title <- ""
      }
      
      log2_subset_GC <- log2(subset)
      df_subset_GC <- melt(log2_subset_GC)
      ggplot(data=as.data.frame(df_subset_GC), aes(x=Var2, y=value)) + 
        geom_point(color="red", size=1.5) + scale_y_continuous(limits = c(-1, 20)) + 
        theme(axis.text.x = element_text(hjust = 1, angle=45)) +
        geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") + 
        scale_x_discrete(labels=groupNames) + xlab("Samples") + ggtitle(paste("Red Channel", title))
    }
    )
    
    #########################################################
    
    ###########        PCA plot
    
    #########################################################
    
    output$pcaPlot <- renderPlot({
      set.palette(n=8, name=setColor())
      ##palette(brewer.pal(8,setColor()))
      colors <- sampleColorsPCA()
      if (method!="Merging"){
        plotPCA(pca = pca, 
                pc1 = input$pc1,
                pc2 = input$pc2,
                col = colors,
                covariates = covariates,
                selectedCov = input$phenotype,
                title = "Beta-values (PCA)"
        )
      }
    })
    
    output$pca_mPlot <- renderPlot({
      set.palette(n=8, name=setColor())
      ##palette(brewer.pal(8,setColor()))
      colors <- sampleColorsPCA()
      if (method!="Merging"){
        plotPCA(pca = pca_m, 
                pc1 = input$pc1,
                pc2 = input$pc2,
                col = colors,
                covariates = covariates,
                selectedCov = input$phenotype,
                title = "M-values (PCA)"
                
        )
      }
    })
    
    
    ## Print the summary of the PCA regression: 
    output$modelPrint <- renderPrint({
      if (method!="Merging"){
        y <- shinyMethylSet1@pca$scores[,as.numeric(input$pcToExplore)]
        cov <- covariates[match(rownames(shinyMethylSet1@pca$scores),
                                rownames(covariates)),]
        x <- (as.factor(cov[,match(input$covToRegress,colnames(cov))]))
        model <- lm(y~ x)
        return(summary(model))
      }
      else {
        return("Merging was used to create the shinyMethylSet1; no PCA was performed. ")
      }
    })
    
    output$barplotPval <- renderPlot({
      #pal <- brewer.pal(8,"Paired")
      pval_means <- colMeans(detP)
      df_pval_means <- as.data.frame(pval_means)
      colnames(df_pval_means) <- "pvals"
      ggplot(df_pval_means, aes(x=rownames(df_pval_means), y=pvals)) +
        geom_col(show.legend = FALSE, color="darkgrey") +
        theme_classic() + 
        scale_y_continuous(limits=c(0,0.08)) + 
        geom_hline(yintercept = 0.05, color="red") + 
        theme(axis.text.x = element_text(hjust = 1, angle=45)) +
        labs(y="P-values", x="")+
        geom_hline(yintercept = 0.01, color="green")
    })
    
    output$probesFailedPlot <- renderPlot({
      plotFailedPropProbes(detP = detP, targets$Sample_Name)
    })
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.pdf",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(directory, "qc.Rmd")
        #file.copy("report1.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(data = shinyMethylSet1)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )


    
    output$beanPlot <-renderPlot({
      # Function to subset the CpGs
      if (input$mOrBeta=="Beta-value"){
        numberOfCpGs = 10000
        pal = brewer.pal(8, "Dark2")
        idx_cg <- sample(nrow(bMatrix),numberOfCpGs)
        b_subset <- as.matrix(bMatrix[idx_cg, ])
        x <- melt(b_subset, varnames=c("cpg","sample"))
        o <- order(colnames(bMatrix))
        #ggplot(x, aes(x=sample, y=value)) + geom_violin() + coord_flip()
        minfi::densityBeanPlot(bMatrix, sampGroups = covariates$Sample_Group, sampNames = covariates$Sample_Name)
      }
    }
    )
    
    
    output$SampleSheetSubset <- renderTable(
      covariates[, 1:10]    
    )
    
    output$SampleSheetInfo <- renderText({
      print(paste("The experiment contains", length(covariates$Sample_Name), "samples.", sep = " "))
      print(paste("The experiment contains", length(unique(covariates$Sample_Group)), "groups", sep = " "))
    }
    )
    
    
    norm <- eventReactive(input$normButton, {
      otherNorm <- TRUE
      norm_method <- input$normID
      if (norm_method == "Quantile"){
        mSetSq <- preprocessQuantile(RGSET)
        shinySummarizeNorm(mSetSq)
      } else if (norm_method == "Funnorm") {
        mSetSq <- preprocessFunnorm(RGSET)
        shinySummarizeNorm(mSetSq)
      } else if (norm_method == "SWAN") {
        mSetSq <- preprocessSWAN(RGSET)
        # Convert to RatioSet object.
        mSetSq <- ratioConvert(mSetSq)
        # Convert to GenomicRatioSet object.
        mSetSq <- mapToGenome(mSetSq)
        shinySummarizeNorm(mSetSq)
      } else if (norm_method == "ssNoob") {
        normalized
      } else if (norm_method == "Illumina") {
        mSetSq <- preprocessIllumina(RGSET, bg.correct = TRUE, normalize = "controls",
                                     reference = 1)
        # Convert to RatioSet object.
        mSetSq <- ratioConvert(mSetSq)
        # Convert to GenomicRatioSet object.
        mSetSq <- mapToGenome(mSetSq)
        shinySummarizeNorm(mSetSq)
      } else {
        stop("[ERROR] normalization method not correctly specified.")
      }
      
      #print(paste("shinyMethylSet created!"))
    })
    
    output$SNPsPlot <- renderPlot({
    colnames(snps) <- covariates[,"Sample_Name"]
      #heatmap.2(snps, col=rev(heat.colors(16)), trace = "none")
      heatmap.2(snps, trace = "none", cexCol = 0.5, col=colorRampPalette(c("yellow","orange","red"))(100))
    }
    )
    
  
    output$design <- renderPlot({
      color <- covariates[,match(input$phenotype,colnames(covariates))]
      plotDesign450k(as.character(sampleNames), covariates = color ,
                     legendTitle = input$phenotype)
    }
    )
    
    ###############################################################
    #    Remove poor quality samples
    ###############################################################
    # 
    # output$poorQualitySample <- renderPrint({
    #   remove_samples()[1]
    #   remove_samples()[2]
    #   remove_samples()[3]
    #   
    # })

    remove_samples <- eventReactive(input$removeSamples, {
      flag <- TRUE
      keep <- colMeans(detP) < input$pvalThreshold
      RGSET.filtered <- RGSET[, keep]
      targets.filtered <- covariates[keep,]
      # remove poor quality samples from detection p-value table
      detP.filtered <- detP[,keep]
      list(RGSET.filtered, targets.filtered, detP.filtered)      
      # print(dim(detP))
      # print(paste("Poor quality samples removed"))
    })
    
    ###############################################################
    #    Remove failed probes
    ###############################################################
    
    output$failedProbes <- renderPrint({
       remove_probes()
    })
    
    remove_probes <- eventReactive(input$filtering, {
      
      if ("failed" %in% input$whatRemove){
        detP <- detP[match(featureNames(GRSet), rownames(detP)),]
        # remove any probes that have failed in one or more samples
        keep_nonfailedProbes <- rowSums(detP < input$pvalThresholdProbes) == ncol(GRSet) 
        GRSet <- GRSet[keep_nonfailedProbes, ]
      } 
      
      if ("sex" %in% input$whatRemove){
        # if your data includes males and females, remove probes on the sex chromosomes
        keep_nonsexchr <- !(featureNames(GRSet) %in% ann$Name[ann$chr %in% 
                                                              c("chrX","chrY")])
        GRSet <- GRSet[keep_nonsexchr, ]
      }
     
      
      if ("reactive" %in% input$whatRemove){
        xReactiveProbes <- read.csv(file=paste("/mnt/ElRaid/amontalban/PROJECTS/methylation/illumina450k_filtering-master",
                                               "48639-non-specific-probes-Illumina450k.csv",
                                               sep="/"), stringsAsFactors=FALSE)
        keep_reactiveProbes <- !(featureNames(GRSet) %in% xReactiveProbes$TargetID)
        GRSet <- GRSet[keep_reactiveProbes, ]
      }
      
      
      if ("SNP" %in% input$whatRemove){
        GRSet <- dropLociWithSnps(GRSet)
      } 
      GRSet
    })
    
    output$factor <- renderPrint({
      factorOfInterest <- factor(groupNames)
      print((factorOfInterest[1]))
    })
    
    run_analysis <- eventReactive(
      input$analysis, {
        GRSet.filtered <- remove_probes()
        mVals <- getM(GRSet.filtered)
        factorOfInterest <- factor(groupNames)
        design <- model.matrix(~0+factorOfInterest, data = covariates)
        colnames(design) <- levels(factorOfInterest)
        fit <- lmFit(mVals, design)
        print(unique(factorOfInterest)[[1]])
        contr <- paste(as.character(unique(factorOfInterest)[2]), as.character(unique(factorOfInterest)[1]), sep = "-")
        print(contr)
        contMatrix <- makeContrasts(contrasts = contr, levels = design)
        fit2 <- contrasts.fit(fit, contMatrix)
        fit2 <- eBayes(fit2)
        annSub <- ann[match(rownames(mVals), ann$Name), c(1:4, 22:24)]
        DMPs <- topTable(fit2, num=Inf, coef = 1, genelist = annSub)
        DMPs
    })
    
    output$resultsCpG <- renderTable(
      head(run_analysis())
    )
    
    
    output$topCpGs <- renderPlot({
      GRSet.filtered <- remove_probes()
      bVals <- getBeta(GRSet.filtered)
      par(mfrow=c(1,2))
      sapply(rownames(run_analysis())[1:2], function(cpg){
        plotCpg(bVals, cpg=cpg, pheno=groupNames, ylab = "Beta values")
      })
    })
    
  }}

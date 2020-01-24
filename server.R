server.methylation <- function(shinyMethylSet1, shinyMethylSet2=NULL){
  function(input, output, session) { 
    betaQuantiles   <-  getBeta(shinyMethylSet1)
    mQuantiles      <-  getM(shinyMethylSet1)
    methQuantiles   <-  getMeth(shinyMethylSet1)
    unmethQuantiles <-  getUnmeth(shinyMethylSet1)
    cnQuantiles     <-  getCN(shinyMethylSet1)
    bMatrix <- shinyMethylSet1@betaMatrix
    if (!is.null(shinyMethylSet2)){
      bMatrix2 <- shinyMethylSet2@betaMatrix
      mMatrix2 <- shinyMethylSet2@mMatrix
      
    }
    #mMatrix <- shinyMethylSet1@mMatrix
    greenControls   <-  getGreenControls(shinyMethylSet1)
    redControls     <-  getRedControls(shinyMethylSet1)
    covariates      <<- pData(shinyMethylSet1)
    pca             <-  getPCA(shinyMethylSet1)$scores
    detP            <- shinyMethylSet1@detP
    sampleNames     <-  sampleNames(shinyMethylSet1)
    slideNames      <- substr(sampleNames,1,10)
    arrayNames      <- substr(sampleNames,14,19)
    plateNames      <- substr(sampleNames,21,30)
    groupNames      <- substr(sampleNames, 32, 43)
    controlNames    <-  names(greenControls)

    RGSET           <- shinyMethylSet1@RGSET
    
    
    
    
    
    method <- shinyMethylSet1@originObject
    
    
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
      if (is.null(shinyMethylSet2)){
        return(NULL)
      } else {
        index <- match(input$probeType,
                       c("I Green","I Red","II","X","Y"))
      }
      if (input$mOrBeta=="Beta-value"){
        bw <- 1
        quantiles <- shinyMethylSet2@betaQuantiles[[index]]
      } else {
        bw <- 1
        quantiles <- shinyMethylSet2@mQuantiles[[index]]
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
    
    
    ## Density plot for normalized data:
    output$normDensities <- renderPlot({
      if (is.null(shinyMethylSet2)){
        return(NULL)
      }	else {
        
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
          quantiles =  shinyMethylSet2@betaQuantiles[[index]]
          minfi::densityPlot(bMatrix2, sampGroups = groupNames)
        } else {
          xlim <- c(-8,8)
          ylim <- c(0, 0.35)
          from = -10; to = 10;
          main = "Normalized data (M values)" 
          xlab = "M-values"
          quantiles =  shinyMethylSet2@mQuantiles[[index]]
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
        
      }})
    
    output$click_info <- renderPrint({
      cat("str(input$click_action):")
      str(input$click_action)
    })
    
    
    output$betamatrixDownload <- downloadHandler(
      
      filename = function() {
        paste("b_values", ".csv", sep = "")
      },
      content = function(file) {
        write.csv(bMatrix2, file)
      }
    )
    output$plotDownload <- downloadHandler(
      
      filename = function() {
        paste("qcMB_report", ".png", sep = "")
      },
      content = function(file) {
        png(file)
        minfi::densityPlot(bMatrix, sampGroups = covariates$Sample_Group)
        minfi::densityPlot(bMatrix2, sampGroups = covariates$Sample_Group)
        dev.off()
      }
    )  
    
    output$controlTypePlotGreen <- renderPlot({
      
      if (input$controlType %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I", 
                                   "SPECIFICITY II", "TARGET REMOVAL")){
        threshold <- 1
      } else if (input$controlType %in% c("EXTENSION", "STAINING", "NON-POLYMORPHIC")){
        threshold <- 5 # you can increase the threshold
      } else {threshold <- 0}
      
      cnt <- input$controlType
      log2_subset_GC <- log2(greenControls[[cnt]])
      df_subset_GC <- melt(log2_subset_GC)
      ggplot(data=as.data.frame(df_subset_GC), aes(x=Var2, y=value)) + 
        geom_point(color="darkgreen", size=1.5) + scale_y_continuous(limits = c(-1, 20)) + 
        theme(axis.text.x = element_text(hjust = 1, angle=45)) +
        geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") +
        scale_x_discrete(labels=plateNames) + xlab("Samples") + ggtitle("Green Channel")
      
    }
    )
    
    output$controlTypePlotRed <- renderPlot({
      
      if (input$controlType %in% c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I", 
                                   "SPECIFICITY II", "TARGET REMOVAL")){
        threshold <- 1
      } else if (input$controlType %in% c("EXTENSION", "STAINING", "NON-POLYMORPHIC")){
        threshold <- 5 # you can increase the threshold
      } else {threshold <- 0}
      
      cnt <- input$controlType
      log2_subset_GC <- log2(redControls[[cnt]])
      df_subset_GC <- melt(log2_subset_GC)
      ggplot(data=as.data.frame(df_subset_GC), aes(x=Var2, y=value)) + 
        geom_point(color="red", size=1.5) + scale_y_continuous(limits = c(-1, 20)) + 
        theme(axis.text.x = element_text(hjust = 1, angle=45)) +
        geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") + 
        scale_x_discrete(labels=plateNames) + xlab("Samples") + ggtitle("Red Channel") 
    }
    )
    
    #########################################################
    
    ###########        PCA plot
    
    #########################################################
    
    output$pcaPlot <- renderPlotly({
      set.palette(n=8, name=setColor())
      ##palette(brewer.pal(8,setColor()))
      if (method!="Merging"){
        plotPCA(pca = pca, 
                pc1 = input$pc1,
                pc2 = input$pc2,
                col = sampleColors(),
                covariates = covariates,
                selectedCov = input$phenotype
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
      pal <- brewer.pal(8,"Paired")
      pval_means <- colMeans(detP)
      df_pval_means <- as.data.frame(pval_means)
      colnames(df_pval_means) <- "pvals"
      ggplot(df_pval_means, aes(x=rownames(df_pval_means), y=pvals, fill=pal[factor(covariates$Sample_Name)])) +
        geom_col(show.legend = FALSE) +
        scale_y_continuous(limits=c(0,0.08)) + 
        geom_hline(yintercept = 0.05, color="red") + 
        theme(axis.text.x = element_text(hjust = 1, angle=45)) +
        geom_hline(yintercept = 0.01, color="green")
    })
    
    output$probesFailedPlot <- renderPlot({
      plotFailedPropProbes(detP = detP, covariates$Sample_Name)
    })
    
    output$reportDownload <- downloadHandler(
      filename = function(){
        paste("report", ".pdf", sep = "")
      },
      content = function(file){
        pdf(file)
        par(mfrow=c(2,2))
        minfi::densityPlot(bMatrix, sampGroups = covariates$Sample_Group)
        minfi::densityPlot(bMatrix2, sampGroups = covariates$Sample_Group)   
        minfi::densityPlot(mMatrix2, sampGroups = covariates$Sample_Group)   
        plotFailedPropProbes(detP = detP, covariates$Sample_Name)
        dev.off()
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

    output$SampleSheetInfo <- renderPrint({
      print(paste("The experiment contains", length(covariates$Sample_Name), "samples.", sep = " "))
      print(paste("The experiment contains", length(unique(covariates$Sample_Group)), "groups", sep = " "))
    }
    )
    
    
    ntext <- eventReactive(input$normButton, {
      norm_method <- input$normID
      if (norm_method == "Quantile"){
        mSetSq <- preprocessQuantile(RGSET)
      } else if (norm_method == "Funnorm") {
        mSetSq <- preprocessFunnorm(RGSET)
      } else if (norm_method == "SWAN") {
        mSetSq <- preprocessSWAN(RGSET)
      } else if (norm_method == "ssNoob") {
        mSetSq <- preprocessNoob(RGSET, dyeMethod = "single")
        # Convert to RatioSet object.
        mSetSq <- ratioConvert(mSetSq)
        # Convert to GenomicRatioSet object.
        mSetSq <- mapToGenome(mSetSq)
      } else if (norm_method == "Illumina") {
        mSetSq <- preprocessIllumina(RGSET, bg.correct = TRUE, normalize = "controls",
                                     reference = 1)
        # Convert to RatioSet object.
        mSetSq <- ratioConvert(mSetSq)
        # Convert to GenomicRatioSet object.
        mSetSq <- mapToGenome(mSetSq)
      } else {
        stop("[ERROR] normalization method not correctly specified.")
      }
      prov <- shinySummarizeNorm(mSetSq)
      #print(paste("shinyMethylSet created!"))
    })
    
    
    
    output$nText <- renderTable({
      ntext()@phenotype
    })
    
  }}

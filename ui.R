ui.methylation <- function(shinyMethylSet1, shinyMethylSet2 = NULL){
  
  betaQuantiles   <- getBeta(shinyMethylSet1) 
  mQuantiles      <- getM(shinyMethylSet1) 
  methQuantiles   <- getMeth(shinyMethylSet1)
  unmethQuantiles <- getUnmeth(shinyMethylSet1)
  cnQuantiles     <- getCN(shinyMethylSet1)
  bMatrix <- shinyMethylSet1@betaMatrix
  if (!is.null(shinyMethylSet2)){
    bMatrix2 <- shinyMethylSet2@betaMatrix
    mMatrix2 <- shinyMethylSet2@mMatrix
  }
  #mMatrix <- shinyMethylSet1@mMatrix
  greenControls   <- getGreenControls(shinyMethylSet1)
  redControls     <- getRedControls(shinyMethylSet1)
  covariates      <<-pData(shinyMethylSet1)
  pca             <- getPCA(shinyMethylSet1)$scores
  detP            <- shinyMethylSet1@detP 
  sampleNames     <- sampleNames(shinyMethylSet1)
  slideNames      <- substr(sampleNames,1,10)
  arrayNames      <- substr(sampleNames,14,19)
  plateNames      <- substr(sampleNames,21,30)
  groupNames      <- substr(sampleNames, 32, 43)
  targets <- shinyMethylSet1@phenotype
  controlNames    <- names(greenControls)
  slides <- unique(slideNames)
  method <- shinyMethylSet1@originObject
  sampleColors <<- as.numeric(as.factor(plateNames))
  
  ## In the case covariates is empty:
  if (ncol(covariates)==0){
    covariates <- data.frame(slide = slideNames, plate = plateNames)
    rownames(covariates) <- sampleNames
    covariates <<- covariates
  }
  fluidPage(
    ###########################  ---  Header ------------ ##############
    mainPanel(
      navbarPage(
        ###########################  ---  Home   ------------ 
        tabPanel("Home"),
        tabPanel("Parameters",
                 
                 HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Color choice:</span></p>"),
                 
                 selectInput("colorChoice", "Color set:",
                             list("Pault","Rainbow","Set1","Set2","Set3",
                                  "Paired","Dark2","Accent"),
                             multiple=FALSE, 
                             selected="Set1"),
                 HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Quality control exploration:</span></p>"),
                 
                 # radioButtons("mOrBeta", "Methylation measure:",
                 #             list("Beta-value","M-value"),
                 #            selected="M-value"),
                 
                 
                 if ("plate" %in% colnames(covariates)){
                   choices <- colnames(covariates)
                   selectInput("phenotype", "Choose a phenotype:",
                               choices,
                               multiple=FALSE, 
                               selected ="plate")
                 } else {
                   choices <- colnames(covariates)
                   selectInput("phenotype", "Choose a phenotype:",
                               choices,
                               multiple=FALSE)},          	
                 
                 checkboxInput("mean","Average density by phenotypic level"),
                 selectInput("lty", "Density line type (lty):",
                             list(1,2,3,4,5,6),
                             multiple=FALSE, 
                             selected=1),
                 selectInput("lwd", "Density line width (lwd):",
                             list(0.5,1,1.5,2,3,4,5,6,7,8,9,10),
                             multiple=FALSE, 
                             selected=1)	
        ),
        tabPanel("Samples", 
                 verticalLayout(
                   plotOutput("barplotPval"),
                   plotOutput("probesFailedPlot")
                 )),
        ###########################  ---  Quality control --- 
        tabPanel("Methylation Measures",
                 sidebarPanel(
                   radioButtons("mOrBeta", "Methylation measure:",
                                list("Beta-value","M-value"),
                                selected="M-value"),
                   selectInput("probeType", "Choose a probe type for the density curves:", 
                               choices = c("I Green","I Red","II","X","Y"),
                               selected="I Green")
                 ),
                 ## Densities plot:
                 HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),
                 HTML('</td><td>'),
                 ## Fast quality control plot:
                 HTML('</td></tr></table>'),
                 plotOutput("rawDensities", click="click_action"),
                 verbatimTextOutput("click_info"),
                 HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Normalized data:</span></p>"),
                 conditionalPanel(condition= "!is.null(shinyMethylSet2)",
                                  plotOutput("normDensities",
                                             click = "normHover")),
                 HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Download beta-values:</span></p>"),
                 downloadLink("betamatrixDownload", "b_values.csv"),
                 HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Download report:</span></p>"),
                 downloadLink("plotDownload", "qcM.png")
                 #verbatimTextOutput(outputId = "cumulativeListPrint"),
                 #downloadLink("selectedSamples","selectedSamples.csv")
                 
        ),
        tabPanel("Control Type",
                 selectInput("controlType", "Choose a control type:", 
                             choices = controlNames,selected=controlNames[1]),
                 HTML("<p><span style=\"color:#336666;font-size:16px\">
			      Step control</span></p>"),
                 verticalLayout(
                   plotOutput("controlTypePlotGreen"),
                   plotOutput("controlTypePlotRed"))
        ),
        tabPanel("PCA",
                 sidebarLayout(
                   plotlyOutput("pcaPlot"),
                   verticalLayout(
                     HTML("
<p style=\"color:#000000;font-size:17px\">A. Choose two principal components to visualize: </span></p>
"),
                     selectInput("pc1", "PC in x:",
                                 seq(1,ncol(betaQuantiles[[1]]),1),
                                 multiple=FALSE, 
                                 selected=1),
                     selectInput("pc2", "PC in y:",
                                 seq(1,ncol(betaQuantiles[[1]]),1),
                                 multiple=FALSE, 
                                 selected=2),
                     HTML("
<p style=\"color:#000000;font-size:17px\">B. Choose a principal component to explore below:</span></p>
"),   
                     selectInput("pcToExplore", "PC:",
                                 seq(1,ncol(betaQuantiles[[1]]),1),
                                 multiple=FALSE, 
                                 selected=1),
                     HTML("
<p style=\"color:#000000;font-size:17px\">C. Choose a covariate to regress against the chosen PC:</span></p>
"),   
                     if (exists("covariates")){
                       choices <- colnames(covariates)
                       selectInput("covToRegress", "Covariate:",
                                   choices,
                                   multiple= FALSE)
                     } else {
                       selectInput("covToRegress", "Covariate:",
                                   list("Batch"),
                                   multiple= FALSE, 
                                   selected="Batch")})),
                 HTML("
<p style=\"color:#000000;font-size:17px\">Association with the PC:</span></p>
"), 
                 verbatimTextOutput(outputId = "modelPrint")
        ),
        ######################   ----   Type I/TypeII Bias --------  
        tabPanel("Probe Type Diff",
                 selectInput("selectedSampleBias", "Sample:",
                             sampleNames,
                             multiple= FALSE),
                 plotOutput("probeBiasPlot"),
                 conditionalPanel(condition= "!is.null(shinyMethylSet2)", plotOutput("probeBiasPlotNorm"))
        ),
        tabPanel("Downloads", 
                 checkboxGroupInput("selectedPlots", "Select:", choices = c("Raw Beta-values", 
                                                                            "Raw m-values", "Normalized beta-values", "Failed Probes", "PCA")),
                 downloadButton("reportDownload", label = "Download"))
      )
    )
  )
  
}
ui.methylation <- function(shinyMethylSet1, shinyMethylSet2 = NULL){
  
  betaQuantiles   <- getBeta(shinyMethylSet1) 
  mQuantiles      <- getM(shinyMethylSet1) 
  methQuantiles   <- getMeth(shinyMethylSet1)
  unmethQuantiles <- getUnmeth(shinyMethylSet1)
  cnQuantiles     <- getCN(shinyMethylSet1)
  bMatrix <- shinyMethylSet1@betaMatrix
  
  mMatrix <- shinyMethylSet1@mMatrix
  greenControls   <- getGreenControls(shinyMethylSet1)
  redControls     <- getRedControls(shinyMethylSet1)
  covariates      <<-pData(shinyMethylSet1)
  pca             <- getPCA(shinyMethylSet1)$scores
  pca_m <- shinyMethylSet1@pca_m$scores
  detP            <- shinyMethylSet1@detP 
  sampleNames     <- sampleNames(shinyMethylSet1)
  slideNames      <- substr(sampleNames,1,10)
  sampleNames     <-  sampleNames(shinyMethylSet1)
  slideNames      <- shinyMethylSet1@phenotype$Slide
  arrayNames      <- shinyMethylSet1@phenotype$Array
  plateNames      <- shinyMethylSet1@phenotype$Sample_Plate
  groupNames      <- shinyMethylSet1@phenotype$Sample_Group
  RGSET           <- shinyMethylSet1@RGSET
  controlNames    <- names(greenControls)
  slides <- unique(slideNames)
  method <- shinyMethylSet1@originObject
  sampleColors <<- as.numeric(as.factor(plateNames))
  snps <- shinyMethylSet1@snps

  ## In the case covariates is empty:
  if (ncol(covariates)==0){
    covariates <- data.frame(slide = slideNames, plate = plateNames)
    rownames(covariates) <- sampleNames
    covariates <<- covariates
  }
  fluidPage(tags$head(includeCSS("/home/amontalban/Documents/methylation-shiny/www/united.min.css")),
            tags$head(
              tags$style(HTML("
      @import url('//fonts.googleapis.com/css?family=Lobster|Cabin:400,700');
    "))
            ),
            
            headerPanel(
              h1("Shiny Methyl ",
                 style = "font-family: 'Sans-serif', bold;
        font-weight: 500; line-height: 1.1;
        color:   #d66958;")),
            ###########################  ---  Header ------------ ##############
            mainPanel(
              navbarPage(
                ###########################  ---  Home   ------------ 
                title=p(strong("Quality control")),
                tabPanel("Home", 
                         h3("Quality control"), 
                         br(),
                        
                         downloadButton("report", "Generate report"),
                         hr(),
                         HTML("<div id=label>  </div>") ,
                         p(strong("Summary")),
                         p(strong("Tab 1: Experiment")),
                         p(strong("Tab 2: Samples")),
                         p(strong("Tab 3: Methylation Measures")),
                         p(strong("Tab 4: Control Type")),
                         p(strong("Tab 5: PCA")),
                         p(strong("Tab 6: Report Generator")),
                         br(),
                         br(),
                         br(),
                         br(),
                         
                         
                        
                         
                ),
                tabPanel("Experiment",
                         tableOutput("SampleSheetSubset"),
                         textOutput("SampleSheetInfo"),
                         plotOutput("SNPsPlot")
                ),
                tabPanel("Samples", 
                         verticalLayout(
                           HTML("
                       <div id=label> · Examine mean detection p-values across all samples to identify any failed 
samples. </div>"),
                           br(),
                           br(),
                           br(),
                           br(),
                           plotOutput("barplotPval"),
                           HTML("
                       <div id=label> · Examine the proportion of failed probes. </div>"), 
                           br(),
                           br(),
                           br(),
                           plotOutput("probesFailedPlot")
                         )),
                ###########################  ---  Quality control --- 
                tabPanel("Methylation Measures",
                         sidebarPanel(
                           HTML("<p><span style=\"color:#d66958;font-size:16px\">
			      Color choice:</span></p>"),
                           
                           selectInput("colorChoice", "Color set:",
                                       list("Pault","Rainbow","Set1","Set2","Set3",
                                            "Paired","Dark2","Accent"),
                                       multiple=FALSE, 
                                       selected="Set1"),
                           HTML("<p><span style=\"color:#d66958;font-size:16px\">
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
                                       selected=1),
                           radioButtons("mOrBeta", "Methylation measure:",
                                        list("Beta-value","M-value"),
                                        selected="Beta-value"),
                           selectInput("probeType", "Choose a probe type for the density curves:", 
                                       choices = c("I Green","I Red","II","X","Y"),
                                       selected="I Green"),
                           hr(),
                           hr(),
                           hr(),
                           HTML("<p><span style=\"color:#d66958;font-size:16px\">
			      Download beta-values:</span></p>"),
                           downloadLink("betamatrixDownload", "raw_bvalues.csv"),
                           br(),
                           br(),
                           HTML("*Only "),
                           br(), 
                           downloadLink("betamatrixNormDownload", "norm_bvalues.csv"),
                         ),
                         mainPanel(
                           selectInput("normID", "Choose Normalization method:", choices = c("Illumina", "SWAN", "Quantile", "ssNoob", "Funnorm", selected = "ssNoob")),
                           actionButton("normButton", "Normalize"),
                           ## Densities plot:
                           HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),
                           HTML('</td><td>'),
                           ## Fast quality control plot:
                           HTML('</td></tr></table>'),
                           plotOutput("rawDensities", click="click_action"),
                           verbatimTextOutput("click_info"),
                           
                           plotOutput("normDensities"),
                           plotOutput("beanPlot")
                           
                           #verbatimTextOutput(outputId = "cumulativeListPrint"),
                           #downloadLink("selectedSamples","selectedSamples.csv")
                           
                         )),
                tabPanel("Control Type",
                         HTML("<p><span style=\"color:#d66958;font-size:16px\">
			      Step control</span></p>"),
                         selectInput("controlType", "Choose a control type:", 
                                     choices = controlNames,selected=controlNames[1]),
                         conditionalPanel(condition="length(sampleNames) >= 50", selectInput("arrayID", "Select array:", arrayNames)),
                         verticalLayout(
                           plotOutput("controlTypePlotGreen"),
                           plotOutput("controlTypePlotRed"))
                ),
                tabPanel("PCA",
                         sidebarLayout(
                          
                           splitLayout(
                           plotOutput("pcaPlot"),
                           plotOutput("pca_mPlot")
                           ),
                           
                           verticalLayout(
                             HTML("
<p style=\"color:#000000;font-size:17px\">A. Choose two  components to visualize: </span></p>
"),
                             if (exists("covariates")){
                               choices <- colnames(covariates)
                               selectInput("phenoID", "Phenotype:",
                                           choices,
                                           multiple= FALSE)
                             } else {
                               selectInput("phenoID", "Phenotype:",
                                           list("Batch"),
                                           multiple= FALSE, 
                                           selected="Batch")})),
                              
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
                                           selected="Batch")},
                         HTML("
<p style=\"color:#000000;font-size:17px\">Association with the PC:</span></p>
"), 
                         verbatimTextOutput(outputId = "modelPrint")
                ),
                ######################   ----   Type I/TypeII Bias --------  
                
                tabPanel("Report", 
                         checkboxGroupInput("selectedPlots", "Select:", choices = c("Raw Beta-values", 
                                                                                    "Raw m-values", "Normalized beta-values", "Failed Probes", "PCA")),
                         #downloadButton("report", "Generate report")
                ))
            )

  )
  
}
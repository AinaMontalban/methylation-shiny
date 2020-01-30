# Aina Montalban
require(knitr)
require(limma)
require(minfi)
#require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#require(IlluminaHumanMethylation450kmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
require(RColorBrewer)
require(missMethyl)
require(Gviz)
#require(DMRcate)
require(stringr)
require(ggplot2)
require(shiny)
require(base)
require(reshape2)
require(plotly)
library(d3heatmap)
require(gplots)
require(shinycssloaders)
# Set the directories where the data is stored
# Full path of the folder
full_path <- "/mnt/ElRaid/amontalban/PROJECTS/methylation/subset"
# Path of the IDAT files
directory.IDAT<- file.path(full_path, "tmp", fsep = "/")
# Path of the samplesheet file
directory.targets <- file.path(full_path, "data", fsep = "/")
print(paste(directory.targets))
# Directory to save the results
directory.results <- file.path(full_path, "aina/results", fsep = "/")

#===================#
# A) Read targets
#===================#
targets <- read.metharray.sheet(directory.targets, pattern = "sample_sheet_PDX_subset.csv")

# Create a column called Bassename which specifies the location of each individual IDAT file
# in the experiment
targets$Basename <- file.path(directory.IDAT, paste(targets$Slide, targets$Array, sep = "_"))

# Convert Sample_Name as character.
targets$Sample_Name <- as.character(targets$Sample_Name)
# Print how many samples do we have
print(paste("Full sample sheet rows (samples):", nrow(targets)))

# Remove samples without responder info.
targets <- targets[targets$Sample_Group != "", ]
print(paste("Sample sheet rows (samples) after removing samples w/o info:", 
            nrow(targets)))


#===================#
# B) Raw data
#===================#

# Read raw data
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
# Give the samples descriptive names
targets$ID <- paste(targets$Slide,targets$Array, targets$Sample_Name,targets$Sample_Group, sep=".")
# Define the sample names into the rgSet
sampleNames(rgSet) <- targets$ID
# print
rgSet


directory <- "~/Documents/methylation-shiny"
#directory <- "/home/aina/Internship/methylation-shiny"
source(paste0(directory, "/", "ui.R"))
source(paste0(directory, "/", "server.R"))
source(paste0(directory, "/", "shinyMethylSet.R"))
source(paste0(directory, "/", "plotDensities.R"))
source(paste0(directory, "/", "plotPCA.R"))
source(paste0(directory, "/", "plotPropFailedProbes.R"))
source(paste0(directory, "/", "run-methylation.R"))
source(paste0(directory, "/", "plotDesign.R"))



### Example minfi - bioconductor
baseDir <- system.file("extdata", package="minfiData")
targetsEx <- read.metharray.sheet(baseDir)
RGSetEx <- read.metharray.exp(targets = targetsEx)
#summary1 <- shinySummarize(RGSetEx)
summary2 <- shinySummarizepr(RGSetEx)

targetsEx$ID <- paste(targetsEx$Slide,targetsEx$Array, targetsEx$Sample_Name,targetsEx$Sample_Group, sep=".")
mSetSq <- preprocessQuantile(RGSetEx)
summary1.norm <- shinySummarizeNorm(mSetSq)
shinyApp(ui=ui.methylation(summary1), server = server.methylation(summary1))
shinyApp(ui=ui.methylation(summary2), server = server.methylation(summary2))
shinyApp(ui=ui.methylation(summary1,summary1.norm ), server = server.methylation(summary1, summary1.norm))



### Real data
summary <- shinySummarizepr(rgSet)
#GRSet.norm <- preprocessNoob(rgSet)
#mSetSq <- ratioConvert(GRSet.norm)
# Convert to GenomicRatioSet object.
mSetSq <- mapToGenome(mSetSq)
summary.norm <- shinySummarize(mSetSq)
shinyApp(ui=ui.methylation(summary), server = server.methylation(summary))
shinyApp(ui=ui.methylation(summary,summary.norm ), server = server.methylation(summary, summary.norm))
  
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
# list the files
list.files(dataDirectory, recursive = TRUE)
targetsEx1 <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
targetsEx1
rgSetEx1 <- read.metharray.exp(targets=targetsEx1)
# give the samples descriptive names
targetsEx1$ID <- paste(targetsEx1$Sample_Group,targetsEx1$Sample_Name,sep=".")
sampleNames(rgSetEx1) <- targetsEx1$ID
rgSetEx1
# Comprovacio preprocessIllumina
# rgSetEx1Norm <-preprocessIllumina(rgSetEx1)
# densityBeanPlot(rgSetEx1Norm)
# densityPlot(rgSetEx1Norm, sampGroups = pData(rgSetEx1Norm)$Sample_Group)
summary3 <- shinySummarizepr(rgSetEx1)
shinyApp(ui=ui.methylation(summary3), server = server.methylation(summary3))

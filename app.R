#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library("xlsx")
library("ggplot2")
library("vegan")
library("phyloseq")
# library("factoextra") # only a dependency for currently non-implemented functions... does not compile so forget it for now

options(repos = BiocInstaller::biocinstallRepos())

# Define UI
ui <- fluidPage(
   
   # Application title
   titlePanel("CLPP Statistical analysis", windowTitle="CLPP Analysis"),
   
   # Sidebar with a file input
   sidebarLayout(
      sidebarPanel(
        # Excel file input
        fileInput(inputId = "excelFile",
                  "Choose an Excel file containing CLPP data:",
                  accept = c(".xlsx",".xls")),
        actionButton("analyze", "Analyze"),
        # Number of timepoints
        # number_of_sheets=7
        numericInput(inputId = "numSheets",
                     "Number of sheets in Excel workbook (i.e., the number of treatments). The default is 7.",
                     value = 7),
        numericInput(inputId = "timepoint",
                     "Which timepoint is being analyzed? For example, enter 1 for week 1.",
                     value = 1),
        # Location of first well in Excel sheet
        # Default for R1 is A68
        numericInput(inputId = "R1",
                     "Row number within Excel sheets where Replicate 1 readings begin. Default is 68.",
                     value = 68),
        # Default for R2 is A77
        # Default for R3 is A86
        selectInput(inputId = "means",
                    "How should the means be calculated?",
                    c("Across replicate wells of the same chemical" = "chem",
                      "Across replicate wells, then across replicate plates" = "plate",
                      "Do not calculate means" = "no")),
        
        radioButtons(inputId = "ordMethod",
                     "What ordination method should be used?",
                     choices=c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")),
        
        radioButtons(inputId = "myDist",
                     "What distance method should be used?",
                     choices=c("bray","unifrac","jaccard","gower"))
        
      ),
      
      # Show a heatmap of the CLPP data
      mainPanel(
         plotOutput("heatmap"),
         plotOutput("ordination")
      )
   )
)

# Define server logic
server <- function(input, output) {
  
  compoundNames <- as.vector(read.table("./www/compound_names.txt", sep="\n"))
  metaData <- read.table("./www/biolog_metadata.txt", sep="\t", header=T)

  
  loadRawCLPP <- function(excelFile) {

    number_of_sheets=input$numSheets
    timepoint=input$timepoint
    
    myDataSummary <- NULL
    myDataFull <- NULL
    myDataSummaryReps <- NULL

    for (idx in 1:number_of_sheets) {
      
      print(paste0("Numer of sheets: ", number_of_sheets))
      message("Looping through index number:")
      print(idx)
      
      R1rowStart = as.numeric(input$R1)
      R1rowEnd = R1rowStart+7
      spaceBetweenReps = 2
      print(R1rowStart)
      print(R1rowEnd)
      R2rowStart = R1rowEnd+spaceBetweenReps
      R2rowEnd = R2rowStart+7
      R3rowStart = R2rowEnd+spaceBetweenReps
      R3rowEnd = R3rowStart+7
      
      R1 <- as.vector(apply(read.xlsx(excelFile,
                                      sheetIndex = idx,
                                      rowIndex = R1rowStart:R1rowEnd,
                                      colIndex = 1:12,
                                      header = FALSE),
                            2,print))

      R2 <- as.vector(apply(read.xlsx(excelFile,
                                      sheetIndex = idx,
                                      rowIndex = 77:84,
                                      colIndex = 1:12,
                                      header = FALSE),
                            2,print))
      
      R3 <- as.vector(apply(read.xlsx(excelFile,
                                      sheetIndex = idx,
                                      rowIndex = 86:93,
                                      colIndex = 1:12,
                                      header = FALSE),
                            2,print))
      
      R1.summary <- apply(cbind(R1[1:32],R1[33:64],R1[65:96]), 1, mean)
      R2.summary <- apply(cbind(R2[1:32],R2[33:64],R2[65:96]), 1, mean)
      R3.summary <- apply(cbind(R3[1:32],R3[33:64],R3[65:96]), 1, mean)
      
      myDataFull <- cbind(myDataFull, cbind(R1[1:32],R1[33:64],R1[65:96]))
      colnames(myDataFull)[length(colnames(myDataFull))] <- paste0("W",timepoint,"-T",idx,"-R1-ChemRep",seq(1:3))
      myDataFull <- cbind(myDataFull, cbind(R2[1:32],R2[33:64],R2[65:96]))
      colnames(myDataFull)[length(colnames(myDataFull))] <- paste0("W",timepoint,"-T",idx,"-R2-ChemRep",seq(1:3))
      myDataFull <- cbind(myDataFull, cbind(R3[1:32],R3[33:64],R3[65:96]))
      colnames(myDataFull)[length(colnames(myDataFull))] <- paste0("W",timepoint,"-T",idx,"-R3-ChemRep",seq(1:3))

      myDataSummary <- cbind(myDataSummary, R1.summary)
      colnames(myDataSummary)[length(colnames(myDataSummary))] <- paste0("W",timepoint,"-T",idx,"-R1")
      myDataSummary <- cbind(myDataSummary, R2.summary)
      colnames(myDataSummary)[length(colnames(myDataSummary))] <- paste0("W",timepoint,"-T",idx,"-R2")
      myDataSummary <- cbind(myDataSummary, R3.summary)
      colnames(myDataSummary)[length(colnames(myDataSummary))] <- paste0("W",timepoint,"-T",idx,"-R3")
      
      repMean <- apply(cbind(R1.summary, R2.summary, R3.summary), 1, mean)
      myDataSummaryReps <- cbind(myDataSummaryReps, repMean)
      colnames(myDataSummaryReps)[length(colnames(myDataSummaryReps))] <- paste0("W",timepoint,"-T",idx)
      
      message("Finished loop for:")
      print(paste0("W",timepoint,"-T",idx))
      
    }
    
    row.names(myDataFull) <- metaData[,1]
    row.names(myDataSummary) <- metaData[,1]
    row.names(myDataSummaryReps) <- metaData[,1]
    
    my_list <- list(myDataSummary, myDataFull, myDataSummaryReps)
    
    return(my_list)
    
  }
  
  analyzeData <- eventReactive(input$analyze, {
    number_of_sheets=input$numSheets
    timepoint=input$timepoint
    inFile = input$excelFile
    dataList <- loadRawCLPP(inFile$datapath)
    # Slot 1: myDataSummary = means across chemicals
    # Slot 3: myDataSummaryReps = means across chemicals first, then replicates
    # Get data from returned list
    
    compoundTable <- as.matrix(data.frame(row.names=metaData[,1],
                                          Substrate=metaData[,1],
                                          Substrate_class=metaData[,2]))
    
    
    message("Summary across reps:")
    print(dataList[[3]])
    transformedDataSummaryRepsCombined <- t(dataList[[3]])
    print(transformedDataSummaryRepsCombined)
    
    print(dataList[[1]])
    transformedDataSummaryChemsCombined <- t(dataList[[1]])
    print(transformedDataSummaryChemsCombined)
    
    print(dataList[[2]])
    transformedDataSummaryNoMeans <- t(dataList[[2]])
    print(transformedDataSummaryNoMeans)
    
    # Load into Phyloseq
    
    # OTU tables
    OTUreps <- otu_table(transformedDataSummaryRepsCombined, taxa_are_rows = FALSE)
    OTU <- otu_table(transformedDataSummaryChemsCombined, taxa_are_rows = FALSE)
    OTUnomeans <- otu_table(transformedDataSummaryNoMeans, taxa_are_rows = FALSE)
    
    message("OTU objects made...")
    
    # Sample data for means across plates
    sampledataReps <- data.frame(Treatment = rep(as.character(seq(1, number_of_sheets)), times=1),
                                 row.names = colnames(dataList[[3]]),
                                 stringsAsFactors = TRUE)
  
    sampledataReps$Week <- sapply(strsplit(row.names(sampledataReps), "-"), `[`, 1)
    sampledataReps$SampleID <- row.names(sampledataReps)
    
    # Sample data for means across chemicals
    sampledata <- data.frame(Replicate = rep(1:3, times=number_of_sheets, each=1),
                             Treatment = rep(as.character(seq(1, number_of_sheets)), times=1, each=3),
                             row.names = colnames(dataList[[1]]), stringsAsFactors=TRUE)
    sampledata$Week <- sapply(strsplit(row.names(sampledata), "-"), `[`, 1)
    sampledata$SampleID <- row.names(sampledata)
    
    # Sample data for no means
    sampledata_nomeans <- data.frame(Replicate = rep(1:3, times=number_of_sheets*3, each=1),
                             Treatment = rep(as.character(seq(1, number_of_sheets)), times=1, each=3),
                             row.names = colnames(dataList[[2]]), stringsAsFactors=TRUE)
    sampledata_nomeans$Week <- sapply(strsplit(row.names(sampledataReps), "-"), `[`, 1)
    sampledata_nomeans$SampleID <- row.names(sampledata_nomeans)

    # Create phyloseq objects
    physeqReps = phyloseq(otu_table(OTUreps),
                          sample_data(sampledataReps),
                          tax_table(compoundTable))
    
    physeq = phyloseq(otu_table(OTU),
                      sample_data(sampledata),
                      tax_table(compoundTable))
    
    physeqNoMeans = phyloseq(otu_table(OTUnomeans),
                           sample_data(sampledata_nomeans),
                           tax_table(compoundTable))
    
    myOutput <- list(physeq, physeqReps, physeqNoMeans)
    
    return(myOutput)
    # compoundTable <- data.frame(row.names=compoundNames[,1],
    #                             Substrate=compoundNames[,1],
    #                             Substrate_class=gsub('[[:digit:]]+',
    #                             '', compoundNames[,1]))
    # compoundTable <- as.matrix(compoundTable)
    
  })
  
  output$heatmap <- renderPlot({
     message("Rendering plot...")
     # Run data analysis function!
     phyloseqObjectsList <- analyzeData()
     # Use if statement to establish which phyloseq object to use based on user input
     if (input$means=="plate") {
       physeq=phyloseqObjectsList[[2]]
     } 
       else if (input$means=="chem") {
       physeq=phyloseqObjectsList[[1]]
     }
       else if (input$means=="no") {
       physeq=phyloseqObjectsList[[3]]
     }
       
     # Plot heatmap
     plot_heatmap(physeq,
                  "NMDS",
                  "bray",
                  sample.order = sample_names(physeq),
                  taxa.order = metaData$Compound )
   })
  
  output$ordination <- renderPlot({
    # width=input$numSheets
    ordMethod = input$ordMethod
    myDist = input$myDist
    message("Rendering plot...")
    # Run data analysis function!
    phyloseqObjectsList <- analyzeData()
    # Use if statement to establish which phyloseq object to use based on user input
    if (input$means=="plate") {
      physeq=phyloseqObjectsList[[2]]
    } 
    else if (input$means=="chem") {
      physeq=phyloseqObjectsList[[1]]
    }
    else if (input$means=="no") {
      physeq=phyloseqObjectsList[[3]]
    }
    # Plot ordination
    plot_ordination(physeq,
                    ordinate(physeq,
                             ordMethod,
                             myDist),
                    color="Treatment",
                    shape="Week") + 
      geom_point(size=5)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


library(shiny)
library(openxlsx)
library(rstatix)



# Define UI ----
ui <- fluidPage(
  
  titlePanel("RT-qPCR analysis"),
  fluidRow(column(5,wellPanel(
  fileInput("file", h3("File input")),
  textInput("veh", h3("Enter name of reference group"), 
            value = ""),   
  checkboxInput("duplicates", label = "Calculate duplicate means", value = FALSE),
  checkboxInput("outAnalysis", label = "Remove outliers", value = FALSE),
  downloadButton( "download_excel",  "Download Data to Excel"))),
  column(7,tabsetPanel(
    tabPanel("Description", helpText("Simple script for ddCt analysis of qPCR data for pharmacologists. 
It calculates fold change based on raw CT values for reference gene and gene of interest for specific 
treatment groups and returns XLSX file containing summary statistics (including means, standard deviations, 
N, SEM, min, max and normality distributtion) and descriptive statistics from ANOVA and Tukey post-hoc analysis 
of between group differences.",br(),br(),"To use the app:", tags$ol(
  tags$li("format your data accordingly to example below"), 
  tags$li("do not change the column name in first column - treatment"), 
  tags$li("save you file as CSV and upload it (use semicolon as column separator and comma as decimal separator)"),
  tags$li("Type the name of the reference group, eg. VEH"),
  tags$li("Calculate duplicate means or remove outliers if required"),
  tags$li("You have all the results below and graph in the tab in upper-right corner"),
  tags$li("To save the graph just right-click it and choose save as and add image type extensions, for example
          *.jpg or *.png")),
  tags$a(href="https://fritzthescientist.pl/testy/exampleData.csv", "Example data file"),br(),
  tags$a(href="https://fritzthescientist.pl/testy/exampleDuplicates.csv", "Example duplicate data file"))),
  tabPanel("Plot", plotOutput("plot"))
))),

h1(textOutput("error")),
fluidRow(column(12,
tabsetPanel(
  tabPanel("Data", dataTableOutput("dataTable")), 
  tabPanel("Fold change", dataTableOutput("test"),h1(textOutput("errorRef"))), 
  tabPanel("Summary statistics", dataTableOutput("summary"),h1(textOutput("errorRef2"))),
  tabPanel("ANOVA", verbatimTextOutput("aovSummary")),
  tabPanel("Tukey",verbatimTextOutput("tukey"))
)),


  

))
# Define server logic ----
server <- function(input, output) {
   

  #pre-process input data
  dataInput <- reactive({
    req(input$file)
    rawTable <- read.table(input$file$datapath,header=T,sep = ";", dec=",")
    
    #calculate duplicate means
    if(input$duplicates){if (length(rawTable)==5){
      geneOfInterest <- names(rawTable[2])
      geneReference <- names(rawTable[4])
      meanReferenceGene<-data.frame(rawTable[[2]],rawTable[[3]])
      meanReferenceGene<-rowMeans(meanReferenceGene, na.rm=TRUE	)
      meanGeneOfInterest<-data.frame(rawTable[[4]],rawTable[[5]])
      meanGeneOfInterest<-rowMeans(meanGeneOfInterest, na.rm=TRUE)
      newNames <- c("treatment", geneOfInterest, geneReference)
      rawTable<-data.frame(rawTable[[1]], meanReferenceGene, meanGeneOfInterest)
      names(rawTable)<-newNames
      }}
    rawTable
  })
  
  #output errors
  output$error<-renderText({
    req(input$file)
    rawTable<-dataInput()
    if (length(rawTable)!=3){
      x<-"Your data table doesn't have approriate number of columns, please use example templates and don't use whitespace or special characters in headers"
      x
      if (length(rawTable)!=5){x<-"Your data table doesn't have approriate number of columns, please use example templates and don't use whitespace or special characters in headers"
      x
      }}
  })
  output$errorRef<-renderText({
    rawTable<-dataInput()
    req(input$veh=="" || input$veh!=unique(rawTable$treatment))
    if(input$veh==""){
      x<-"Please type the name of the reference group, eg. VEH"
      x}
    else if(any("VEH"!=unique(rawTable$treatment))==FALSE){
      x<-"Please type correct name of the reference group"
      x}
  })
  output$errorRef2<-renderText({
    rawTable<-dataInput()
    req(input$veh=="" || input$veh!=unique(rawTable$treatment))
    if(input$veh==""){
    x<-"Please type the name of the reference group, eg. VEH"
    x}
    else if(any("VEH"!=unique(rawTable$treatment))==FALSE){
      x<-"Please type correct name of the reference group"
      x}
  })
  #output input data table
  output$dataTable <- renderDataTable({ 
    dataInput()
  }, options = list(pageLength = 10))
  
  #Calculate mean fold change
  dataTest <- reactive({
    
    req(input$veh)
    req(input$file)
    rawTable <- dataInput()
    
    #Functions
    meanRefGeneCT<-function(){
      y <- rawTable$treatment==input$veh
      dCTmean<- rawTable[y,]
      dCTmean<- c(dCTmean$dCT)
      dCTmean <- mean(dCTmean)
      return (dCTmean)
    }
    
    ddCTcalc<-function(){
      x<-c(rawTable$dCT)
      ddCT<-x-dCTmean
      return(ddCT)
    }
    
    foldChangeCalc<-function(){
      y<-c(rawTable$ddCT)
      foldChange<-2^(-y)
      return(foldChange)
    }

    #Asign target and reference gene
    referenceGene <- c(rawTable[[2]])
    targetGene <- c(rawTable[[3]])
    #Calculate dCT
    dCT=targetGene-referenceGene
    rawTable=cbind(rawTable,dCT)
    #Mean dCT for reference gene
    dCTmean<-meanRefGeneCT()
  
    #calculate ddCT
    rawTable$ddCT<-ddCTcalc()
    
    #calculate fold change
    rawTable$foldChange<-foldChangeCalc()
    
    #Calculate outliers
    groups<-unique(rawTable$treatment)
    outlier=c()
    outlierPosition=c()
    test=c()
    analysedTable<-rawTable
    tableLength<-1:length(rawTable$foldChange)
    for (name in groups){
      y <- rawTable$treatment==name
      x <- rawTable[y,]
      out <- boxplot(x$foldChange, plot=FALSE)$out
      out_ind <- which(rawTable$foldChange %in% out)
      outlierPosition<-append(outlierPosition, out_ind)
    }
    for (i in tableLength){
      if (i %in% outlierPosition){
        outlier <- append(outlier, TRUE)
      } else {
        outlier <- append(outlier, FALSE)
      }
    } 
    
    analysedTable<-cbind(analysedTable, outlier)
    #OPTIONAL remove outliers
    if (input$outAnalysis){
      outName="_outlier"
      rawTable<-subset(analysedTable, outlier==FALSE)
      dCTmean<-meanRefGeneCT()
      rawTable$ddCT<-ddCTcalc()
      rawTable$foldChange<-foldChangeCalc()
      rawTable
    } else{
      outName=""
      analysedTable
    }
    
  })
  
  #output fold change table
  output$test<-renderDataTable({
    dataTest()
  }, options = list(pageLength = 10))
  

  #calculate summary statistics 
  summaryTable<-reactive({
    rawTable<-dataTest()
    groups<-unique(rawTable$treatment)
    mean_Fold_Change=c()
    sd_Fold_Change=c()
    N=c()
    SEM=c()
    min=c()
    max=c()
    normality_P=c()
    for (name in groups){
      y <- rawTable$treatment==name
      x <- rawTable[y,]
      mean_Fold_Change<-append(mean_Fold_Change, mean(x$foldChange))
      sd_Fold_Change<-append(sd_Fold_Change, sd(x$foldChange))
      N<-append(N, length(x$foldChange))
      SEM<-append(SEM, sd(x$foldChange)/sqrt(length(x$foldChange)))
      min<-append(min, min(x$foldChange))
      max<-append(max,max(x$foldChange))
      shapiro=shapiro.test(x$foldChange)
      normality_P<-append(normality_P, shapiro$p.value)
    }
    DT<-data.frame(groups,mean_Fold_Change,SEM, N, min, max, sd_Fold_Change, normality_P)
    
  })
  
  #output summary statistics
  output$summary<-renderDataTable({
    summaryTable()
  })  
  
  #output descriptive statistics and create variables for xlsx file
  output$aovSummary = reactivePrint(function() {
    summary(aov(foldChange ~ treatment, data = dataTest()))
  })
  
  output$tukey<-reactivePrint(function(){
    foldChangeStat<-aov(foldChange ~ treatment, data = dataTest())
    tukey<-TukeyHSD(foldChangeStat)
    tukey
  })
  ANOVA<-reactive({
    foldChangeStat<-aov(foldChange ~ treatment, data = dataTest())
    anova<-anova_summary(foldChangeStat)
    anova<-anova[1:5]
    anova
  })
  tukeyTest<-reactive({
    foldChangeStat<-aov(foldChange ~ treatment, data = dataTest())
    tukey<-TukeyHSD(foldChangeStat)
    tukey<-data.frame(tukey$treatment)
    tukey<-cbind(comparison = rownames(tukey), tukey)
  })
  
  tableLength<-reactive({
    x<-dataInput()
    x<-length(x[[1]])
    x
  })
  
  #render plot
  output$plot<-renderPlot({
    rawTable<-dataTest()
    boxplot(rawTable$foldChange~rawTable$treatment, ylab="Fold Change",  xlab="")
  })
  
  #output excel file
  output$download_excel <- downloadHandler(
    filename = function() {
      "summary.xlsx"
    },    content = function(file) {
      my_workbook <- createWorkbook()
      addWorksheet(
        wb = my_workbook,
        sheetName = "Fold Change"
      )
      addWorksheet(
        wb = my_workbook,
        sheetName = "Summary Statistics"
      )
      addWorksheet(
        wb = my_workbook,
        sheetName = "Tukey"
      )
      
      writeData(
        my_workbook,
        sheet = 1,
        dataTest()
  )
      writeData(
        wb=my_workbook,
        sheet=2,
        summaryTable()
      )
      writeData(
        wb=my_workbook,
        sheet=3,
        "ANOVA",
        startCol = 1,
        startRow = 1,
      )
      writeData(
        wb=my_workbook,
        sheet=3,
        ANOVA(),
        startCol = 1,
        startRow = 2,
      )
      writeData(
        wb=my_workbook,
        sheet=3,
        "TUKEY",
        startCol = 1,
        startRow = 4,
      )
      writeData(
        wb=my_workbook,
        sheet=3,
        tukeyTest(),
        startCol = 1,
        startRow = 5,
      )
      posStyle <- createStyle(fontColour = "#c20000")
      max<-tableLength()
      conditionalFormatting(
        wb=my_workbook,
        sheet=3,
        cols=1:5,
        rows=1:max,
        rule = "$E1<=0.05",
        style = posStyle
      )
      conditionalFormatting(
        wb=my_workbook,
        sheet=2,
        cols=1:8,
        rows=1:max,
        rule = "$H1<=0.05",
        style = posStyle
      )
      conditionalFormatting(
        wb=my_workbook,
        sheet=1,
        cols=1:7,
        rows=1:max,
        rule = "$G1==TRUE",
        style = posStyle
      )
      saveWorkbook(my_workbook, file)      
})
}

# Run the app ----
shinyApp(ui = ui, server = server)

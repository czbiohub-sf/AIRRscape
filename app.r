### THIS IS THE MAIN SHINY FILE - IN R STUDIO OR SHINYAPPS.IO IT RUNS THE SHOW

#install.packages("shiny")
#install.packages("ggplot2")
#install.packages("tidyverse")
library(shiny)
library(ggplot2)
library(tidyverse)
library(DT)

## for plotting
#library(Biostrings)
library(seqinr)
library(phangorn)
library(ape)

## function here - maybe make it reactive?? - can avoid using alternate conversion code...

# writeFasta<-function(data, filename){
#   fastaLines = c()
#   for (rowNum in 1:nrow(data)){
#     fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"SEQUENCE_ID"], sep = "")))
#     fastaLines = c(fastaLines,as.character(data[rowNum,"JUNCTIONAA"]))
#   }
#   fileConn<-file(filename)
#   writeLines(fastaLines, fileConn)
#   close(fileConn)
# }


### idea about making custom functions - comments by Hadley Wickham:
# Next you could pull out the repeated code into a function:
#   
#   load_path <- function(path) {
#     req(path)
#     
#     ext <- tools::file_ext(path)
#     
#     if (ext == "csv") {
#       read.csv(path, header = TRUE)
#     } else if (ext == "xls" || ext == "xlsx") {
#       read_excel(path)
#     } else {
#       stop("Unknown extension: '.", ext, "'")
#     }
#   }
# file1 <- reactive(load_path(input$selection$datapath[[1]]))
# file2 <- reactive(load_path(input$selection$datapath[[2]]))

##########

## idea of including dataset loading here:
toshiny.few <- read.delim("toshiny_few.tab")
toshiny.few1 <- read.delim("toshiny_few1.tab")
toshiny.few1.h <- read.delim("toshiny_few1_h.tab")
# toshiny.few2 <- read.delim("toshiny_few2.tab")
# toshiny.few3 <- read.delim("toshiny_few3.tab")
# toshiny.few3b <- read.delim("toshiny_few3b.tab")
toshiny.fewhcbg <- read.delim("toshiny_fewhcbg.tab")
toshiny.fewhcbgc <- read.delim("toshiny_fewhcbgc.tab")
# toshiny.natalia4 <- read.delim("toshiny_natalia.tab")
toshiny.cov2sarsmers <- read.delim("toshiny_cov2sarsmers.tab")
toshiny.cov2sarsmersc <- read.delim("toshiny_cov2sarsmersc.tab")

toshiny.hivcov2sarsmers <- read.delim("toshiny_hivcov2sarsmers.tab")
toshiny.hivcov2sarsmersc <- read.delim("toshiny_hivcov2sarsmersc.tab")

# toshiny.augvsnoaug <- read.delim("toshiny_augvsnoaug.tab")
# toshiny.augvsnoaugc <- read.delim("toshiny_augvsnoaugc.tab")

# toshiny.fewhcbghiv <- read.delim("toshiny_fewhcbghiv.tab")
# toshiny.fewhcbghivc <- read.delim("toshiny_fewhcbghivc.tab")
# toshiny.fewhivnih45 <- read.delim("toshiny_fewhivnih45.tab")



toshiny.fewhivmabsnih45rd1214 <- read.delim("toshiny_mabsnih45rd1214.tab")
toshiny.fewhivmabsnih45rd1214c <- read.delim("toshiny_mabsnih45rd1214c.tab")

# toshiny.dengueall <- read.delim("toshiny_dengueall.tab")
# toshiny.dengueallc <- read.delim("toshiny_dengueallc.tab")

# toshiny.dengueallandmice <- read.delim("toshiny_dengueallandmice.tab")
# toshiny.dengueallandmicec <- read.delim("toshiny_dengueallandmicec.tab")

toshiny.fluhcmabs <- read.delim("toshiny_fluhcmabs.tab")
toshiny.fluhcmabsc <- read.delim("toshiny_fluhcmabsc.tab")

## were mutate.if but now as.character but only after loading datasets
#toshiny.fewer %>% mutate_if(is.factor, as.character)
## rename in few2 id to sequence_id - also relabel IDs in HC & comet...
toshiny.few$SEQUENCE_ID <- as.character(toshiny.few$SEQUENCE_ID)
toshiny.few$JUNCTIONAA <- as.character(toshiny.few$JUNCTIONAA)
toshiny.few1$SEQUENCE_ID <- as.character(toshiny.few1$SEQUENCE_ID)
toshiny.few1$JUNCTIONAA <- as.character(toshiny.few1$JUNCTIONAA)
toshiny.few1.h$SEQUENCE_ID <- as.character(toshiny.few1.h$SEQUENCE_ID)
toshiny.few1.h$JUNCTIONAA <- as.character(toshiny.few1.h$JUNCTIONAA)

# toshiny.few2$SEQUENCE_ID <- as.character(toshiny.few2$SEQUENCE_ID)
# toshiny.few2$JUNCTIONAA <- as.character(toshiny.few2$JUNCTIONAA)
# toshiny.few3$SEQUENCE_ID <- as.character(toshiny.few3$SEQUENCE_ID)
# toshiny.few3$JUNCTIONAA <- as.character(toshiny.few3$JUNCTIONAA)
# toshiny.few3b$SEQUENCE_ID <- as.character(toshiny.few3b$SEQUENCE_ID)
# toshiny.few3b$JUNCTIONAA <- as.character(toshiny.few3b$JUNCTIONAA)

toshiny.fewhcbg$SEQUENCE_ID <- as.character(toshiny.fewhcbg$SEQUENCE_ID)
toshiny.fewhcbg$JUNCTIONAA <- as.character(toshiny.fewhcbg$JUNCTIONAA)
toshiny.fewhcbgc$SEQUENCE_ID <- as.character(toshiny.fewhcbgc$SEQUENCE_ID)
toshiny.fewhcbgc$JUNCTIONAA <- as.character(toshiny.fewhcbgc$JUNCTIONAA)

# toshiny.natalia4$SEQUENCE_ID <- as.character(toshiny.natalia4$SEQUENCE_ID)
# toshiny.natalia4$JUNCTIONAA <- as.character(toshiny.natalia4$JUNCTIONAA)

toshiny.cov2sarsmers$SEQUENCE_ID <- as.character(toshiny.cov2sarsmers$SEQUENCE_ID)
toshiny.cov2sarsmers$JUNCTIONAA <- as.character(toshiny.cov2sarsmers$JUNCTIONAA)
toshiny.cov2sarsmersc$SEQUENCE_ID <- as.character(toshiny.cov2sarsmersc$SEQUENCE_ID)
toshiny.cov2sarsmersc$JUNCTIONAA <- as.character(toshiny.cov2sarsmersc$JUNCTIONAA)

toshiny.hivcov2sarsmers$SEQUENCE_ID <- as.character(toshiny.hivcov2sarsmers$SEQUENCE_ID)
toshiny.hivcov2sarsmers$JUNCTIONAA <- as.character(toshiny.hivcov2sarsmers$JUNCTIONAA)
toshiny.hivcov2sarsmersc$SEQUENCE_ID <- as.character(toshiny.hivcov2sarsmersc$SEQUENCE_ID)
toshiny.hivcov2sarsmersc$JUNCTIONAA <- as.character(toshiny.hivcov2sarsmersc$JUNCTIONAA)

# toshiny.augvsnoaug$SEQUENCE_ID <- as.character(toshiny.augvsnoaug$SEQUENCE_ID)
# toshiny.augvsnoaug$JUNCTIONAA <- as.character(toshiny.augvsnoaug$JUNCTIONAA)
# toshiny.augvsnoaugc$SEQUENCE_ID <- as.character(toshiny.augvsnoaugc$SEQUENCE_ID)
# toshiny.augvsnoaugc$JUNCTIONAA <- as.character(toshiny.augvsnoaugc$JUNCTIONAA)

# toshiny.fewhcbghiv$SEQUENCE_ID <- as.character(toshiny.fewhcbghiv$SEQUENCE_ID)
# toshiny.fewhcbghiv$JUNCTIONAA <- as.character(toshiny.fewhcbghiv$JUNCTIONAA)
# toshiny.fewhcbghivc$SEQUENCE_ID <- as.character(toshiny.fewhcbghivc$SEQUENCE_ID)
# toshiny.fewhcbghivc$JUNCTIONAA <- as.character(toshiny.fewhcbghivc$JUNCTIONAA)
# 
# toshiny.fewhivnih45$SEQUENCE_ID <- as.character(toshiny.fewhivnih45$SEQUENCE_ID)
# toshiny.fewhivnih45$JUNCTIONAA <- as.character(toshiny.fewhivnih45$JUNCTIONAA)

toshiny.fewhivmabsnih45rd1214$SEQUENCE_ID <- as.character(toshiny.fewhivmabsnih45rd1214$SEQUENCE_ID)
toshiny.fewhivmabsnih45rd1214$JUNCTIONAA <- as.character(toshiny.fewhivmabsnih45rd1214$JUNCTIONAA)
toshiny.fewhivmabsnih45rd1214c$SEQUENCE_ID <- as.character(toshiny.fewhivmabsnih45rd1214c$SEQUENCE_ID)
toshiny.fewhivmabsnih45rd1214c$JUNCTIONAA <- as.character(toshiny.fewhivmabsnih45rd1214c$JUNCTIONAA)

# toshiny.dengueall$SEQUENCE_ID <- as.character(toshiny.dengueall$SEQUENCE_ID)
# toshiny.dengueall$JUNCTIONAA <- as.character(toshiny.dengueall$JUNCTIONAA)
# toshiny.dengueallc$SEQUENCE_ID <- as.character(toshiny.dengueallc$SEQUENCE_ID)
# toshiny.dengueallc$JUNCTIONAA <- as.character(toshiny.dengueallc$JUNCTIONAA)

# toshiny.dengueallandmicec$SEQUENCE_ID <- as.character(toshiny.dengueallandmicec$SEQUENCE_ID)
# toshiny.dengueallandmicec$JUNCTIONAA <- as.character(toshiny.dengueallandmicec$JUNCTIONAA)
# toshiny.dengueallandmice$SEQUENCE_ID <- as.character(toshiny.dengueallandmice$SEQUENCE_ID)
# toshiny.dengueallandmice$JUNCTIONAA <- as.character(toshiny.dengueallandmice$JUNCTIONAA)

toshiny.fluhcmabs$SEQUENCE_ID <- as.character(toshiny.fluhcmabs$SEQUENCE_ID)
toshiny.fluhcmabs$JUNCTIONAA <- as.character(toshiny.fluhcmabs$JUNCTIONAA)
toshiny.fluhcmabsc$SEQUENCE_ID <- as.character(toshiny.fluhcmabsc$SEQUENCE_ID)
toshiny.fluhcmabsc$JUNCTIONAA <- as.character(toshiny.fluhcmabsc$JUNCTIONAA)

## to re-order any rows for plotting re-run here
toshiny.hivcov2sarsmers$id <- factor(toshiny.hivcov2sarsmers$id, levels = c("anti-CoV2 mAbs", "anti-SARS mAbs", "anti-MERS mAbs", "anti-HIV mAbs"))
toshiny.cov2sarsmers$id <- factor(toshiny.cov2sarsmers$id, levels = c("anti-CoV2 mAbs", "anti-SARS mAbs", "anti-MERS mAbs"))
toshiny.fewhivmabsnih45rd1214$id <- factor(toshiny.fewhivmabsnih45rd1214$id, levels = c("HIV+ patient AD1214", "HIV+ patient NIH45", "anti-HIV mAbs"))
# toshiny.dengueall$id <- factor(toshiny.dengueall$id, levels = c("DEN plasmablasts", "DEN patient 13", "DEN patient 20", "DEN patients OAS data"))
# toshiny.dengueallandmice$id <- factor(toshiny.dengueallandmice$id, levels = c("DEN plasmablasts", "DEN patient 13", "DEN patient 20", "DEN patients OAS data","DEN vaccinated humanized mice"))

ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      # 2 ways to define depending on reactive definitions below
      selectInput("dataset", "Dataset:",
                  c("COVID mAbs",
                    "COVID mAbs - IgH by binding",
                    "COVID mAbs - IgH by neutralization",
                    "Healthy Control vs. COVID mAbs - IgH",
                    "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH",
                    "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH combined",
                    "COVID mAbs vs. SARS vs. MERS mAbs - IgH",
                    "COVID mAbs vs. SARS vs. MERS mAbs - IgH combined",
                    "COVID mAbs vs. SARS vs. MERS vs. HIV mAbs - IgH",
                    "COVID mAbs vs. SARS vs. MERS vs. HIV mAbs - IgH combined",
                    "HIV mAbs vs. HIV+ patient NIH45 vs. HIV+ patient AD1214 - IgH",
                    "HIV mAbs vs. HIV+ patient NIH45 vs. HIV+ patient AD1214 - IgH combined",
                    # "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS - IgH",
                    # "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS - IgH combined",
                    # "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS vs. 2 humanized mice - IgH",
                    # "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS vs. 2 humanized mice - IgH combined",
                    "anti-flu mAbs vs. 3 healthy vaccinated controls - IgH",
                    "anti-flu mAbs vs. 3 healthy vaccinated controls - IgH combined"), selectize = FALSE), 
      selectInput("plotcolors", "Plot Colors:",
                  c("Average SHM",
                    "Percentage of total reads"), selectize = FALSE), 
      p("When plots appear, click on a bin to get a list of mAbs in the lower table. Hovering over a bin will also show some basic stats."),
      br(),
      p("Alternately if you want to see more than a bin you can create a box and all mAbs within will appear in the top table."),
      br(),
      p("From the lower table you can download all or selected mAbs in the chosen bin, download the distance matrix of all mAbs, or create phylogenies of selected mAbs. The final phylogeny options are to find the nearest 50 or 250 sequences of a single selected mAb. Finally make sure to check all mAbs in the table have the same CDR3 length or the distance matrix/tree building will fail."),
      width = 2
      # c("COVID mAbs" = "BX.clusters.mabs",
      #   "Comet Data" = "biohub.all2",
      #   "Healthy Control" = "BX.clusters.hc",
      #   "Healthy Control vs. Comet Data vs. COVID mAbs" = "toshiny.few"), selectize = FALSE),
      
    ),

    mainPanel(
      plotOutput("ggplot1", height = "800px", hover = "plot_hover", click = "plot_click", brush = "plot_brush"),
      uiOutput("hover_info"),
      # plotOutput("ggplot2", height = "800px", click = "plot_click"),
      # uiOutput("hover2_info"),
#      verbatimTextOutput("click_info") - instead trying datatable, next line
## new datatable if double-clicked
      DT::dataTableOutput("brush_info"),
      DT::dataTableOutput("click_info"),
      downloadButton("downloadfilter","Download all data in the clicked bin"),
      downloadButton("downloadfilter2","Download only selected rows"),
      downloadButton("downloadfilter3","Download distance matrix of all data in the clicked bin"),
# actionButton("goButton", "Make simple NJ network of selected CDR3 AA motifs (maximum 25)"),
# actionButton("goButton2", "Make simple parsimony network of selected CDR3 AA motifs (maximum 25)"),
# actionButton("reset", "Clear phylogeny"),
      actionButton("go", "Make phylogeny of selected CDR3 AA motifs (maximum 75)"), 
      selectInput("plottab", "Phylogeny:",
            c("NJ",
              "Parsimony",
              "up to 50 nearest sequences to a single selected mAb - Parsimony",
              "up to 250 nearest sequences to a single selected mAb - Parsimony"), selectize = FALSE), 
      plotOutput("phyloPlot", width = "120%", height = "1700px"),
#      plotOutput("phyloPlot", width = "120%", height = "850px"),
    )
  )
)
  


server <- function(input, output, session) {
  options(width = 180, DT.options = list(pageLength = 10)) # Increase text width for printing table ALSO ADDING DEFAULT NUMBER OF 25 ROWS IN DATATABLE
#  options(DT.options = list(pageLength = 5, language = list(search = 'Filter:')))
  ## adding get to allow selectinput to accept entire dataframe...
    # inputdataset=reactive({get(input$dataset)}) ## then ggplot inputdataset() instead of BX.clusters.cometvsmabsvshc.h.few

    ## early in server now defining facets here
    facetvar1 <- reactive({
      switch(input$dataset, "COVID mAbs" = "CREGION", "COVID mAbs - IgH by binding" = "binding", "COVID mAbs - IgH by neutralization" = "neutralization", "COVID mAbs vs. SARS vs. MERS mAbs - IgH" = "id", "COVID mAbs vs. SARS vs. MERS mAbs - IgH combined" = "CREGION", "COVID mAbs vs. SARS vs. MERS vs. HIV mAbs - IgH" = "id", "COVID mAbs vs. SARS vs. MERS mAbs - IgH combined" = "CREGION", "Healthy Control" = "CREGION", "Healthy Control vs. COVID mAbs - IgH" = "id", "COVID mAbs vs. SARS vs. MERS vs. HIV mAbs - IgH combined" = "CREGION", "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH" = "id", "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH combined" = "CREGION", "HIV mAbs vs. HIV+ patient NIH45 vs. HIV+ patient AD1214 - IgH" = "id", "HIV mAbs vs. HIV+ patient NIH45 vs. HIV+ patient AD1214 - IgH combined" = "CREGION", "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS - IgH" = "id", "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS - IgH combined" = "CREGION", "anti-flu mAbs vs. 3 healthy vaccinated controls - IgH" = "id","Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS vs. 2 humanized mice - IgH" = "id", "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS vs. 2 humanized mice - IgH combined" = "CREGION","anti-flu mAbs vs. 3 healthy vaccinated controls - IgH combined" = "CREGION")
    })
    ## to change x columns vs. leaving fixed (for IgH only datasets) - note if I leave out default is fixed...
    facetvar2 <- reactive({
      switch(input$dataset, "COVID mAbs" = "free_x", "COVID mAbs - IgH by binding" = "fixed", "COVID mAbs - IgH by neutralization" = "fixed", "Healthy Control" = "free_x", "Healthy Control vs. COVID mAbs - IgH" = "fixed", "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAb - IgH" = "fixed")
    })
    
  ## alternate to reactive + get from Shiny advanced example:
  inputdataset <- reactive({
#    switch(input$dataset, "COVID mAbs" = BX.clusters.mabs, "Comet Data" = biohub.all2, "Healthy Control" = BX.clusters.hc, "Healthy Control vs. Comet Data vs. COVID mAbs" = toshiny.few)
## using new names that reduce variables, round
    switch(input$dataset, "COVID mAbs" = toshiny.few1, "COVID mAbs - IgH by binding" = toshiny.few1.h, "COVID mAbs - IgH by neutralization" = toshiny.few1.h, "COVID mAbs vs. SARS vs. MERS mAbs - IgH" = toshiny.cov2sarsmers, "COVID mAbs vs. SARS vs. MERS mAbs - IgH combined" = toshiny.cov2sarsmersc, "COVID mAbs vs. SARS vs. MERS vs. HIV mAbs - IgH" = toshiny.hivcov2sarsmers, "COVID mAbs vs. SARS vs. MERS vs. HIV mAbs - IgH combined" = toshiny.hivcov2sarsmersc, "Healthy Control vs. COVID mAbs - IgH" = toshiny.few, "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH" = toshiny.fewhcbg, "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH combined" = toshiny.fewhcbgc, "HIV mAbs vs. HIV+ patient NIH45 vs. HIV+ patient AD1214 - IgH" = toshiny.fewhivmabsnih45rd1214, "HIV mAbs vs. HIV+ patient NIH45 vs. HIV+ patient AD1214 - IgH combined" = toshiny.fewhivmabsnih45rd1214c, "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS - IgH" = toshiny.dengueall, "Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS - IgH combined" = toshiny.dengueallc,"Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS vs. 2 humanized mice - IgH" = toshiny.dengueallandmice ,"Dengue mAbs vs. DEN patient d13 vs. DEN patient d20 vs. DEN patients via OAS vs. 2 humanized mice - IgH combined" = toshiny.dengueallandmicec,"anti-flu mAbs vs. 3 healthy vaccinated controls - IgH" = toshiny.fluhcmabs, "anti-flu mAbs vs. 3 healthy vaccinated controls - IgH combined" = toshiny.fluhcmabsc)
#    switch(input$dataset, "COVID mAbs" = toshiny.few1, "COVID mAbs - IgH by binding" = toshiny.few1.h, "COVID mAbs - IgH by neutralization" = toshiny.few1.h, "Comet Data" = toshiny.few2, "Healthy Control" = toshiny.few3b, "Healthy Control vs. COVID mAbs - IgH" = toshiny.few, "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH" = toshiny.fewhcbg, "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH combined" = toshiny.fewhcbgc)
    
## if using online version remove sensitive sequences
#    switch(input$dataset, "COVID mAbs" = toshiny.fewer1, "COVID mAbs - IgH by binding" = toshiny.fewer1.h, "COVID mAbs - IgH by neutralization" = toshiny.fewer1.h, "Comet Data" = toshiny.few2, "Healthy Control" = toshiny.few3, "Healthy Control vs. COVID mAbs - IgH" = toshiny.fewer)
    
      })
## trying extra filter for downloading
  # dat <- inputdataset()
  # filteredDS <- nearPoints(dat, input$plot_click, threshold = 3, maxpoints = 99999, addDist = FALSE)
## trying slightly differently - NOVEMBER 2020 CHANGING MAX TO 5000 IN CASE TOO MANY POINTS...
    filteredDS <- reactive({
    nearPoints(inputdataset(), input$plot_click, threshold = 5, maxpoints = 99999, addDist = FALSE)
  })

### 8/28 trying to add code to allow for trees of subsetted data:
    ### notes on how to do this:
    # Thanks to @yihui this is now possible using the DT package and input$tableId_rows_all 
    # where tableID is the id assigned to your table. See the link below for details.
    # 
    # http://rstudio.github.io/DT/shiny.html
    
 #   filteredDS()[input[["click_info_rows_all"]], ]
    ### these 2 commands are okay
    filteredDSpartial <- reactive({
      ids <- input$click_info_rows_selected
      filteredDS()[ids,]
    })
    filteredDSall <- reactive({
      ids <- input$click_info_rows_all
      filteredDS()[ids,]
    })

## then rename sequence_id (note will be same if I rename HC & comet ids beforehand)
#    filteredDSpartial$SEQUENCE_ID <- paste(filteredDSpartial$SEQUENCE_ID,filteredDSpartial$JUNCTIONAA,filteredDSpartial$GENE,sep="_") 
#    filteredDSall$SEQUENCE_ID <- paste(filteredDSall$SEQUENCE_ID,filteredDSall$JUNCTIONAA,filteredDSall$GENE,sep="_") 

    ## these 2 are not working correctly - trying slightly different 8/31 THIS NOW WORKS...
    ## EXCEPT 9/1/2020  the as.character part seems to be catching in the plot function steps below...SOLVED BY RUNNING AS.CHARACTER AT VERY START
    filteredDSpartial2 <- reactive({
      partial2 <- filteredDSpartial() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT)
#      partial2 <- filteredDSpartial() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT) %>% mutate_if(is.factor, as.character)
      partial2$SEQUENCE_ID <- paste(partial2$SEQUENCE_ID,partial2$GENE,partial2$JUNCTIONAA,sep="_")
      # partial2$SEQUENCE_ID <- as.character(partial2$SEQUENCE_ID)
      # partial2$JUNCTIONAA <- as.character(partial2$JUNCTIONAA)
      partial2
    })
    
    filteredDSall2 <- reactive({
      all2 <- filteredDSall() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT)
#      all2 <- filteredDSall() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT) %>% mutate_if(is.factor, as.character)
      all2$SEQUENCE_ID <- paste(all2$SEQUENCE_ID,all2$GENE,all2$JUNCTIONAA,sep="_")
      all2
    })

    filteredDSpartial2id <- reactive({
      partial2id <- filteredDSpartial() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT)
      #      partial2 <- filteredDSpartial() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT) %>% mutate_if(is.factor, as.character)
      partial2id$SEQUENCE_ID <- paste(partial2id$SEQUENCE_ID,partial2id$GENE,partial2id$JUNCTIONAA,sep="_")
      ## error here because I really only want a string from this output
      partial2id$SEQUENCE_ID
      # partial2id$JUNCTIONAA <- NULL
      # partial2id$GENE <- NULL
      # partial2id$GF_JGENE <- NULL
      # partial2id$CDR3LENGTH_IMGT <- NULL
      # partial2id
    })

### ADDING NOV 2020 NEW VARIABLE SEPARATELY CALCULATING DISTANCE MATRIX AS I HAD BEEN BELOW BUT WITHIN TREE-PLOTTING STEP
    ## THIS WOULD BE TO SEPARATELY HAVE AVAILABLE FOR DOWNLOADING...
    matrixDSall2 <- reactive({
    filteredDataall <- filteredDSall2()
    # filteredDataall <- isolate({filteredDSall2()})
    filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
    # filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
    rownames(filteredDataall.y) <- filteredDataall[,1]
    as.AAbin(filteredDataall.y)
    seqsall = as.phyDat(filteredDataall.y, type = "AA")
    dmalldf = dist.hamming(seqsall, ratio = TRUE)
    # dmalldf = dist.ml(seqsall, model="JTT")
    #    dmalldf <- as.data.frame(dmalldf)
    dm.matrixall <- as.matrix(dmalldf) %>% as.data.frame()
    dm.matrixall
    # ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
    # filteredData1ID <- isolate({filteredDSpartial2id()})
    # dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at(filteredData1ID)
    # dmall.matrix <- rownames_to_column(dmall.matrix, var = "SEQUENCE_ID")
    # dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID)
    # dm.matrix <- dmall.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
    #   select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 9999) ## okay if this is more than what is in the set!
    # ## now with new subset make new matrix & tree
    # filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
    # rownames(filteredData.y) <- dm.matrix[,1]
    # as.AAbin(filteredData.y)
    # seqs = as.phyDat(filteredData.y, type = "AA")
    # dm = dist.ml(seqs, model="JTT")
    # dm
    })        

    ## now  need to write as fasta, then reload using ape, then plot tree (tree will have to be in a reactive output)
## think all of the below steps are okay at start of server command, except for final plotting
    
## see up top for custom function and Hadley W commentary
# file2 <- reactive(load_path(input$selection$datapath[[2]]))

    # try writeFasta as a button much like current download button?
  ### MONDAY TRY BELOW COMMAND INSIDE SOMETHING ELSE DOWN AT BOTTOM (NOT DOWNLOAD BUT SIMILAR ACTION BUTTON?)
    ## THEN will somehow need to re-import said fasta file...maybe using Hadley line just above?
    
    
    # filteredDSpartialcommand <- reactive({
    #   writeFasta(filteredDSpartial(), "filteredDSpartial.fasta")    })
    # filteredDSallcommand <- reactive({
    #   writeFasta(filteredDSall(), "filteredDSall.fasta")    })
    
## later commands if I can get fasta files parsed...    see new commands below

    

    
    ## instead of below - easier to manipulate sequence_id name in toshiny.few in HC and comet data (just relabelling TO HC1, HC2, COMET1, COMET2, ETC)
    ## ALSO NEED TO HAVE JUNCTIONAA IN SECOND COLUMN...
    
    # filteredDStoplot <- reactive({
    #   if (input$dataset == "Comet Data"    || input$dataset == "Healthy Control")
    #   {
    #     filteredDStoplot$SEQUENCE_ID <- paste(filteredDStoplot$JUNCTIONAA,filteredDStoplot$GENE,sep="_") 
    #     
    #   } else {
    #     filteredDStoplot$SEQUENCE_ID <- paste(filteredDStoplot$SEQUENCE_ID,filteredDStoplot$JUNCTIONAA,filteredDStoplot$GENE,sep="_") 
    #   }
    # })

    # filteredTable_selected <- reactive({
    #   ids <- input$filteredTable_rows_selected
    #   filteredTable_data()[ids,]
    # })
    
### first plot code    
#     output$ggplot1 <- renderPlot({
#     ## first line for picking among different graphs
#       facet_formula <- as.formula(paste("~", facetvar1()))
# #      facet_wrap(facet_formula, ncol = 2)
#       # WHEN TOP PLOT IS RAW COUNTS IN BINS, USING NEXT LINE
# #   ggplot(inputdataset(), aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
# #    }, height = "auto")
#       # WHEN TOP PLOT IS SHM IN BINS USE THIS LINE
#    ggplot(inputdataset(), aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Mean Somatic Hypermutation as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
#     }, width = 1315, height = "auto")
      
        ## note adding extra expression here from Shiny advanced example
#      dat <- inputdataset()
## if manually choosing input      
  #         ggplot(toshiny.few, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
  # }, height = "auto")
    
### hover data over first plot   
  output$hover_info <- renderUI({
    hover <- input$plot_hover
         dat <- inputdataset()
    ## first line for picking among different graphs- NOTE CHANGE TO 'DAT'
       point <- nearPoints(dat, input$plot_hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    ## if manually choosing input
    # point <- nearPoints(toshiny.few, input$plot_hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px / 2 - 20, "px; top:", top_px / 1.8, "px;")
    
    # actual tooltip created as wellPanel
    ## NOTE I THOUGHT I WAS GETTING COUNT BUT THAT WAS JUST ROW NUMBER - TO ACTUALLY GET COUNTS NEED TO SOMEHOW GET FROM STAT_BIN ANALYSIS...
    wellPanel(
      style = style,
      p(HTML(paste0("<b> V-gene & J-gene: </b>", point$GF_JGENE, "<br/>",
                    "<b> CDR3 Length (aa): </b>", point$CDR3LENGTH_IMGT, "<br/>",
                    "<b> Mean Somatic Hypermutation (%): </b>", point$shm.mean, "<br/>",
                    "<b> Count: </b>", point$ncount, "<br/>")))
    )
  })
  
  # output$hover_info <- renderUI({
  #   req(input$plot_hover) 
  #   verbatimTextOutput("vals")
  # })
  # 
  # 
  # output$vals <- renderPrint({
  #    # nearPoints(BX.clusters.cometvsmabsvshc.h.few, input$plot_hover, threshold = 5, maxpoints = 25, addDist = FALSE)
  #   hover <- input$plot_hover
  #   print(str(hover)) # list
  #   y <- nearPoints(BX.clusters.cometvsmabsvshc.h.few, input$plot_hover)[input$var_y]
  #   req(nrow(y) != 0)
  #   y
  # })
  # 
### second plot code    NOW ONLY PLOT
    output$ggplot1 <- renderPlot({
      ## first line for picking among different graphs
#      ggplot(inputdataset(), aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill= (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ get(input$facet_by), ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
## if manually choosing input  
#      isolate({
      facet_formula <- as.formula(paste("~", facetvar1()))
      dat <- inputdataset()
#      ggplot(dat, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with % of total reads as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      data2 <- if (input$plotcolors == "Average SHM") {
        toplot <- ggplot(dat, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "Mean \nSHM (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Mean Somatic Hypermutation as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      } else {
        toplot <- ggplot(dat, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "% of \nReads  ", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with % of total reads as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      }     
#      })
      toplot
    }, width = 1200, height = "auto")

#      output$click_info <- renderPrint({   # instead line below
    
## NEW SECOND HOVER - doesn't seem to work...
    # ### hover data over first plot   
    # output$hover2_info <- renderUI({
    #   hover2 <- input$plot2_hover
    #   dat <- inputdataset()
    #   ## first line for picking among different graphs- NOTE CHANGE TO 'DAT'
    #   point2 <- nearPoints(dat, input$plot2_hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    #   ## if manually choosing input
    #   # point <- nearPoints(toshiny.few, input$plot_hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    #   if (nrow(point2) == 0) return(NULL)
    #   
    #   # calculate point position INSIDE the image as percent of total dimensions
    #   # from left (horizontal) and from top (vertical)
    #   left_pct2 <- (hover2$x - hover2$domain$left) / (hover2$domain$right - hover2$domain$left)
    #   top_pct2 <- (hover2$domain$top - hover2$y) / (hover2$domain$top - hover2$domain$bottom)
    #   
    #   # calculate distance from left and bottom side of the picture in pixels
    #   left_px2 <- hover2$range$left + left_pct2 * (hover2$range$right - hover2$range$left)
    #   top_px2 <- hover2$range$top + top_pct2 * (hover2$range$bottom - hover2$range$top)
    #   
    #   # create style property fot tooltip
    #   # background color is set so tooltip is a bit transparent
    #   # z-index is set so we are sure are tooltip will be on top
    #   style2 <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
    #                   "left:", left_px2 / 2 - 20, "px; top:", top_px2 / 1.8, "px;")
    #   
    #   # actual tooltip created as wellPanel
    #   ## NOTE I THOUGHT I WAS GETTING COUNT BUT THAT WAS JUST ROW NUMBER - TO ACTUALLY GET COUNTS NEED TO SOMEHOW GET FROM STAT_BIN ANALYSIS...
    #   wellPanel(
    #     style = style2,
    #     p(HTML(paste0("<b> V-gene & J-gene: </b>", point2$GF_JGENE, "<br/>",
    #                   "<b> CDR3 Length (aa): </b>", point2$CDR3LENGTH_IMGT, "<br/>",
    #                   "<b> Mean Somatic Hypermutation (%): </b>", point2$shm.mean, "<br/>",
    #                   "<b> Count: </b>", point2$ncount, "<br/>")))
    #   )
    # })
    # 
    
### datatable under second plot 
      output$click_info <- DT::renderDataTable({
          
        # With ggplot2, no need to tell it what the x and y variables are.
        # threshold: set max distance, in pixels
        # maxpoints: maximum number of rows to return
        # addDist: add column with distance, in pixels
        ## TO DO: SHOW FEWER COLUMNS, BUT ALSO SHOW SUMMARY OF EACH SQUARE???
        ## first line for picking among different graphs - NOTE CHANGE TO 'DAT'
             dat <- inputdataset()
             old <- options(width = 1600); on.exit(options(old))
        nearPoints(dat, input$plot_click, threshold = 5, maxpoints = 99999,
        # nearPoints(toshiny.few, input$plot_click, threshold = 3, maxpoints = 25,
                              
        #      nearPoints(h3.cometvsmabsvshc.h2.stat$count, xvar = BX.clusters.cometvsmabsvshc.h$GF_JGENE, yvar = BX.clusters.cometvsmabsvshc.h$CDR3LENGTH_IMGT, input$plot_click, threshold = 5, maxpoints = 25,
                   addDist = FALSE)
      }, width = 1600)
      
## dec 2020 trying to add (separately?) a datatable with all of the points when double-clicking plot_brush brush_info
      output$brush_info <- DT::renderDataTable({
        
        # With ggplot2, no need to tell it what the x and y variables are.
        # threshold: set max distance, in pixels
        # maxpoints: maximum number of rows to return
        # addDist: add column with distance, in pixels
        ## TO DO: SHOW FEWER COLUMNS, BUT ALSO SHOW SUMMARY OF EACH SQUARE???
        ## first line for picking among different graphs - NOTE CHANGE TO 'DAT'
        dat0 <- inputdataset()
        old <- options(width = 1600); on.exit(options(old))
        brushedPoints(dat0, input$plot_brush, allRows = FALSE)
        # nearPoints(dat0, input$plot_dblclick, threshold = 500, maxpoints = 99999, allRows = FALSE,
        #            addDist = FALSE)
      }, width = 1600)
      
      
## this allows user to download all sequences in clicked bin (first filter) or only selected rows (second filter)
      output$downloadfilter <- downloadHandler(
        filename = function() {
          paste('Filtered data-', Sys.Date(), '.csv', sep = '')
        },
        content = function(file){
          write.csv(filteredDS()[input[["click_info_rows_all"]], ],file, row.names = FALSE)
#          write.csv(inputdataset()[input[["click_info_rows_all"]], ],file)
        }
      )
      
      output$downloadfilter2 <- downloadHandler(
        filename = function() {
          paste('Filtered data-', Sys.Date(), '.csv', sep = '')
        },
        content = function(file){
          write.csv(filteredDS()[input[["click_info_rows_selected"]], ],file, row.names = FALSE)
#          write.csv(inputdataset()[input[["click_info_rows_selected"]], ],file)
        }
      )      

### adding nov 2020 third download button for matrix of distances...
      output$downloadfilter3 <- downloadHandler(
        filename = function() {
          paste('Distance matrix-', Sys.Date(), '.csv', sep = '')
        },
        content = function(file){
          write.csv(matrixDSall2(),file, row.names = FALSE)
          # write.csv(matrixDSall2()[input[["click_info_rows_all"]], ],file, row.names = FALSE)
        }
      )          
#    }


      
      ## 9/14 dueling (NJ vs parsimony) phylogeny of CDR3 motifs of selected sequences from datatable
      v <- reactiveValues(doPlot = FALSE)
      
      observeEvent(input$go, {
        # 0 will be coerced to FALSE
        # 1+ will be coerced to TRUE
        v$doPlot <- input$go
      })
      
      observeEvent(input$plottab, {
        v$doPlot <- FALSE
      })  
      
      output$phyloPlot <- renderPlot({
        if (v$doPlot == FALSE) return()
        
        isolate({
          data <- if (input$plottab == "NJ") {
            filteredData <- isolate({filteredDSpartial2()})
            filteredData.y <- t(sapply(strsplit(filteredData[,2],""), tolower))
            rownames(filteredData.y) <- filteredData[,1]
            as.AAbin(filteredData.y)
            seqs = as.phyDat(filteredData.y, type = "AA")
            dm = dist.ml(seqs, model="JTT")
            tree = midpoint(NJ(dm))
            tree$edge.length[tree$edge.length<0] <- 0
            tree$edge.length <- tree$edge.length * 0.15
          } else if (input$plottab == "Parsimony") {
            filteredData <- isolate({filteredDSpartial2()})
            filteredData.y <- t(sapply(strsplit(filteredData[,2],""), tolower))
            rownames(filteredData.y) <- filteredData[,1]
            as.AAbin(filteredData.y)
            seqs = as.phyDat(filteredData.y, type = "AA")
            dm = dist.ml(seqs, model="JTT")
            tree1 = NJ(dm)
            tree1$edge.length[tree1$edge.length<0] <- 0
            tree <- pratchet(seqs, start = tree1, maxit=100,
                             minit=5, k=5, trace=0)
          } else if (input$plottab == "up to 50 nearest sequences to a single selected mAb - Parsimony") {
            ## initial test of tree of all sequences - this works fine
            # filteredData <- isolate({filteredDSall2()})
            # filteredData.y <- t(sapply(strsplit(filteredData[,2],""), tolower))
            # rownames(filteredData.y) <- filteredData[,1]
            # as.AAbin(filteredData.y)
            # seqs = as.phyDat(filteredData.y, type = "AA")
            # dm = dist.ml(seqs, model="JTT")
            # tree = midpoint(NJ(dm))
            # tree$edge.length[tree$edge.length<0] <- 0
            ## first make matrix of all sequences in the bin
            filteredDataall <- isolate({filteredDSall2()})
            filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
            rownames(filteredDataall.y) <- filteredDataall[,1]
            as.AAbin(filteredDataall.y)
            seqsall = as.phyDat(filteredDataall.y, type = "AA")
            dmall = dist.ml(seqsall, model="JTT")
            ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
            filteredData1ID <- isolate({filteredDSpartial2id()})
            dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at(filteredData1ID)
            dmall.matrix <- rownames_to_column(dmall.matrix, var = "SEQUENCE_ID")
            ### for distance thresholding
            ## MAR 2021 INSTEAD TRY THESE LINES - THESE WORK:
            dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID, filteredData1ID)
            colnames(dmall.matrix)[2] <- "DIST"
            # dmall.matrix <- rename.vars(dmall.matrix, c("SEQUENCE_ID","filteredData1ID"), c("SEQUENCE_ID","DIST"))
            dmall.matrix <- dmall.matrix %>% filter(DIST < 0.505)
            #### end distance thresholding
            dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID)
            dm.matrix <- dmall.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
              select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 50) ## okay if this is more than what is in the set!
            ## now with new subset make new matrix & tree
            filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
            rownames(filteredData.y) <- dm.matrix[,1]
            as.AAbin(filteredData.y)
            seqs = as.phyDat(filteredData.y, type = "AA")
            dm = dist.ml(seqs, model="JTT")
            tree1 = NJ(dm)
            tree1$edge.length[tree1$edge.length<0] <- 0
            tree <- pratchet(seqs, start = tree1, maxit=100,
                             minit=5, k=5, trace=0)
            filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
          } else {
            ## initial test of tree of all sequences - this works fine
            # filteredData <- isolate({filteredDSall2()})
            # filteredData.y <- t(sapply(strsplit(filteredData[,2],""), tolower))
            # rownames(filteredData.y) <- filteredData[,1]
            # as.AAbin(filteredData.y)
            # seqs = as.phyDat(filteredData.y, type = "AA")
            # dm = dist.ml(seqs, model="JTT")
            # tree = midpoint(NJ(dm))
            # tree$edge.length[tree$edge.length<0] <- 0
            ## first make matrix of all sequences in the bin
            filteredDataall <- isolate({filteredDSall2()})
            filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
            rownames(filteredDataall.y) <- filteredDataall[,1]
            as.AAbin(filteredDataall.y)
            seqsall = as.phyDat(filteredDataall.y, type = "AA")
            dmall = dist.ml(seqsall, model="JTT")
            ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
            filteredData1ID <- isolate({filteredDSpartial2id()})
            dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at(filteredData1ID)
            dmall.matrix <- rownames_to_column(dmall.matrix, var = "SEQUENCE_ID")
            ### for distance thresholding
            ## MAR 2021 INSTEAD TRY THESE LINES - THESE WORK:
            dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID, filteredData1ID)
            colnames(dmall.matrix)[2] <- "DIST"
            # dmall.matrix <- rename.vars(dmall.matrix, c("SEQUENCE_ID","filteredData1ID"), c("SEQUENCE_ID","DIST"))
            dmall.matrix <- dmall.matrix %>% filter(DIST < 0.66)
            #### end distance thresholding
            dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID)
            dm.matrix <- dmall.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
              select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 260) ## okay if this is more than what is in the set!
            ## now with new subset make new matrix & tree
            filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
            rownames(filteredData.y) <- dm.matrix[,1]
            as.AAbin(filteredData.y)
            seqs = as.phyDat(filteredData.y, type = "AA")
            dm = dist.ml(seqs, model="JTT")
            tree1 = NJ(dm)
            tree1$edge.length[tree1$edge.length<0] <- 0
            tree <- pratchet(seqs, start = tree1, maxit=100,
                             minit=5, k=5, trace=0)
            filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
          }
          
          # par(family = "mono", mar=c(2, 1, 2, 0) + 0.5)
          par(family = "mono", mar=c(2, 0.0, 2, 0) + 0.0)
          ### this changes the colors of all 'controls' to gray but leaves covid mabs black
          tipcolors <- def(tree$tip.label, "HC" = "gray", "Boyd" = "black", "Galson" = "black", "clonotype" = "black", "N152" = "thistle", "IAVI84" = "plum", "AD1214" = "orchid", "nih45" = "pink", "-au" = "brown", "d13" = "orange", "d20" = "blue", "OAS" = "gold", "mab" = "orchid", default = "orange", regexp = TRUE)
          # plot(bird.orders, tip.color = co2)
          plot(tree, lab4ut="axial",
               edge.width=2, label.offset = 0, cex = 0.70, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors)  ## x.lim breaks parsimony and some NJ
          #               edge.width=2, label.offset = 0, cex = 0.95, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
          gjcdr3.title <- filteredData
          gjcdr3.title$CDR3LENGTH_IMGT <- as.numeric(gjcdr3.title$CDR3LENGTH_IMGT)
          # gjcdr3.title$G_J_CDR3 <- paste(gjcdr3.title$GF_JGENE,gjcdr3.title$CDR3LENGTH_IMGT,sep="_")
          gjcdr3.title$G_J_CDR3 <- paste("Phylogeny of selected CDR3 motifs:", gjcdr3.title$GF_JGENE, gjcdr3.title$CDR3LENGTH_IMGT, "aa", sep=" ")
          gjcdr3.title <- gjcdr3.title %>%
            select(G_J_CDR3)
          gjcdr3.titleID <- gjcdr3.title$G_J_CDR3[1]
          title(gjcdr3.titleID, family = "sans")
          # title("Phylogeny of CDR3 motifs", family = "sans")  "Significance level is" ~
          # title(main = expression("Phylogeny of selected CDR3 motifs -" ~ gjcdr3.titleID), family = "sans")

        })
      })
      
}

### 9/22 thoughts on selecting single sequence and sorting all sequences in a bin by distance to select sequence:

# will need to use filteredDSall2
# then will need to make dist.ml on this full dataset, then sort the matrix and somehow extract largest numbers...

##You can coerce it into a matrix, and then sort it
##Then as @Roland suggested, convert it back to dist
#T.mat <- as.matrix(T)[ordering, ordering]
#T  <- as.dist(T.mat)

### use of a regexp (so we need to quote it) to colour all orders
### with names starting with "C" (and change the default):
# co2 <- def(bird.orders$tip.label, "^C" = "gold", default = "grey", regexp = TRUE)
# plot(bird.orders, tip.color = co2)

## plan B - take the saved csv from above and use as import for tree vis app (either separate or add on to this)
## note need to change the toshiny.few columns so first 2 are id and junction_AA
## the idea is to not make phylogenies of the entire V genes but just cluster the CDR3 aa motifs...

## note John had idea to query CDR3 motifs and look for similar motifs outside of a given bin (this seems really hard)



## inputdataset() not saving the bin, but selected from original dataset. here is explanation:
# The problem arises because when downloading, you are not referencing your filtered table, 
# but the original table, and you apply the filtered line numbers on the original table. 
# When downloading, you need to reference your filtered table, that is possible if you store 
# that in a reactive value and use that reactive in building the datatable AND the download:

## trying to add extra variable:
#filteredDS <- nearPoints(dat, input$plot_click, threshold = 5, maxpoints = 99999, addDist = FALSE)
#filteredDS <- nearPoints(inputdataset(), input$plot_click, threshold = 5, maxpoints = 99999, addDist = FALSE)

## store the currently filtered DS in a reactive
# filteredDS <- reactive({
#   if (!"ALL" %in% input$product){
#     return(DS[DS$PRODUCT %in% input$product,])
#   }else{
#     return(DS)
#   }
# })
# 
# # display the currently filtered DS
# output$ds <- DT::renderDataTable({
#   filteredDS()
# },
# rownames = T,
# server = F)
# 
# # download filtered rows
# output$downloadFiltered <- downloadHandler(
#   filename = "filteredData.csv",
#   content = function(file){
#     s = input$ds_rows_all
#     write.csv(filteredDS()[s, , drop = F], file, row.names = T)
#   })
# 
# # download selected rows
# output$downloadSelected <- downloadHandler(
#   filename = "selectedData.csv",
#   content = function(file){
#     s = input$ds_rows_selected
#     write.csv(filteredDS()[s, , drop = F], file, row.names = T)
#   }
# )
# }

## width isn't working here - very strange possible work-around - note instead went wiith DT renderDataTable
# output$out2 <- renderText({
#   old <- options(width = 20); on.exit(options(old))
#   paste(capture.output(print(1:100)), collapse = "\n")
# })

# xvar <- reactive({
#   if ((input$plot_type == "base"    && input$plot_scaletype   == "x_factor") ||
#       (input$plot_type == "ggplot2" && input$ggplot_scaletype == "x_factor"))
#   {
#     switch(input$dataset, mtcars = "cyl", diamonds = "cut", grid = "xf")
#     
#   } else if ((input$plot_type == "base"    && input$plot_scaletype   == "datetime") ||
#              (input$plot_type == "ggplot2" && input$ggplot_scaletype == "datetime")) {
#     "datetime"
#   } else {
#     switch(input$dataset, mtcars = "wt", diamonds = "carat", grid = "x")
#   }
# })


shinyApp(ui = ui, server = server)


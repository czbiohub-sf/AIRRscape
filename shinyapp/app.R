#install.packages("shiny")
#install.packages("ggplot2")
#install.packages("alakazam")
#install.packages("tidyverse")
#install.packages("DT")
library(shiny)
library(ggplot2)
#library(alakazam)
library(tidyverse)
library(DT)

## for plotting
#install.packages("seqinr")
#install.packages("phangorn")
#install.packages("ape")
#install.packages("shinyscreenshot")
library(seqinr)
library(phangorn)
library(ape)
library(shinyscreenshot)

setwd("~/data_carpentry/AIRRscape/shinyapp")
options(shiny.maxRequestSize=200*1024^2)
##########

## datasets to load
toshiny.cov2.abdab <- read_tsv("toshiny_cov2_abdab.tab")
toshiny.cov2.abdab.h <- read_tsv("toshiny_cov2_abdab_h.tab")

toshiny.cov2hiv <- read_tsv("toshiny_cov2hiv.tab")
toshiny.cov2hivc <- read_tsv("toshiny_cov2hivc.tab")

toshiny.cov2.all <- read_tsv("toshiny_cov2_all.tab")
toshiny.cov2.allc <- read_tsv("toshiny_cov2_allc.tab")

toshiny.hiv.all <- read_tsv("toshiny_hiv_all.tab")
toshiny.hiv.allc <- read_tsv("toshiny_hiv_allc.tab")

toshiny.den.all <- read_tsv("toshiny_den_all.tab")
toshiny.den.allc <- read_tsv("toshiny_den_allc.tab")

#toshiny.cov2hivden.allc <- read.delim("toshiny_cov2hivden_allc.tab")

## then for each dataframe change sequence_id & cdr3_aa_imgt columns to character
## if using read_tsv do not need, but also need to change to non-tibble df for treebuilding
toshiny.cov2.abdab <- toshiny.cov2.abdab %>% as.data.frame()
toshiny.cov2.abdab.h <- toshiny.cov2.abdab.h %>% as.data.frame()

toshiny.cov2hiv <- toshiny.cov2hiv %>% as.data.frame()
toshiny.cov2hivc <- toshiny.cov2hivc %>% as.data.frame()

toshiny.cov2.all <- toshiny.cov2.all %>% as.data.frame()
toshiny.cov2.allc <- toshiny.cov2.allc %>% as.data.frame()

toshiny.hiv.all <- toshiny.hiv.all %>% as.data.frame()
toshiny.hiv.allc <- toshiny.hiv.allc %>% as.data.frame()

toshiny.den.all <- toshiny.den.all %>% as.data.frame()
toshiny.den.allc <- toshiny.den.allc %>% as.data.frame()

# toshiny.cov2hivden.allc$sequence_id <- as.character(toshiny.cov2hivden.allc$sequence_id)
# toshiny.cov2hivden.allc$cdr3_aa_imgt <- as.character(toshiny.cov2hivden.allc$cdr3_aa_imgt)

## to re-order any rows for plotting re-run here
toshiny.cov2.abdab.h$binding <- factor(toshiny.cov2.abdab.h$binding, levels = c("RBD", "non-RBD"))
toshiny.cov2.all$id <- factor(toshiny.cov2.all$id, levels = c("SARS-CoV2 mAbs", "DEN patient 13 bulk repertoire", "COVID-19 patient Binder p11 bulk repertoire", "COVID-19 patient Galson p1 bulk repertoire","COVID-19 patient Kuri-Cervantes m5 bulk repertoire", "COVID-19 patient Nielsen p7450 bulk repertoire", "Healthy control bulk repertoire"))
toshiny.cov2hiv$id <- factor(toshiny.cov2hiv$id, levels = c("SARS-CoV2 mAbs", "HIV mAbs"))
toshiny.den.all$id <- factor(toshiny.den.all$id, levels = c("Dengue plasmablasts", "Dengue patient d13 bulk repertoire", "Dengue Parameswaran 2013 patient bulk repertoires"))


ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      selectInput("dataset", "Dataset:",
                  c("SARS-CoV2 mAbs - heavy chains & light chains",
                    "SARS-CoV2 mAbs - IgH by binding",
                    "SARS-CoV2 mAbs - IgH by neutralization",
                    "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH",
                    "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH combined",
                    "SARS-CoV2 mAbs vs. HIV mAbs - IgH",
                    "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined",
                    "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire vs. HIV Setliff 2018 patient bulk repertoires - IgH",
                    "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire vs. HIV Setliff 2018 patient bulk repertoires - IgH combined",
                    "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH",
                    "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH combined",
                    "Your datasets - IgH",
                    "Your datasets - IgH combined"), selectize = FALSE),
      selectInput("plotcolors", "Plot Colors:",
                  c("Average SHM",
                    "Maximum SHM",
                    "Percentage of total antibody sequences"), selectize = FALSE), 
      fileInput(
        inputId = "calfile1", 
        label = "Select converted & combined tsv/tab files - first the separate datasets", 
        multiple = TRUE,
        accept = c(".tsv",".tab")),
      fileInput(
        inputId = "calfile2", 
        label = "Select converted & combined tsv/tab files - then the combined datasets", 
        multiple = TRUE,
        accept = c(".tsv",".tab")),
      p("When plots appear, click on a bin to get a list of antibodies in the lower table. Hovering over a bin will show some basic stats."),
      br(),
      p("Alternately if you want to see more than a bin you can create a box and all antibodies within will appear in the top table."),
      br(),
      p("From the lower table you can download all or selected antibodies in the chosen bin, download the distance matrix of all antibodies, or create topologies of selected antibodies. The last topology options are to find the nearest sequences (up to 500) of a single selected antibody, with four possible distance thresholds. Note that you can change the window size of the topology using the slider."),
      br(),
      p("Finally make sure to check all antibodies in the table have the same CDR3 length or the topology calculation will fail."),
      width = 2
      
    ),
    
    mainPanel(
      plotOutput("ggplot1", width = "120%", height = "800px", hover = hoverOpts(id = "plot_hover", delay = 300, delayType = c("debounce", "throttle")), click = "plot_click", brush = "plot_brush"),
      uiOutput("hover_info"),
      # uiOutput("hover_info", style = "pointer-events: none"),
      DT::dataTableOutput("brush_info"),
      DT::dataTableOutput("click_info"),
      actionButton("go", "Make topology of selected CDR3 AA motifs"), 
      downloadButton("downloadfilter","Download all data in the clicked bin"),
      downloadButton("downloadfilter2","Download only selected rows"),
      downloadButton("downloadfilter3","Download distance matrix of all data in the clicked bin"),
      selectInput("plottab", "Topology:",
                  c("NJ",
                    "Parsimony",
                    "up to 500 nearest sequences to a single selected mAb - Parsimony; 100% CDR3 identity (Briney 2019)",
                    "up to 500 nearest sequences to a single selected mAb - Parsimony; 80% CDR3 identity (Soto 2019)",
                    "up to 500 nearest sequences to a single selected mAb - Parsimony; 70% CDR3 identity (Setliff 2018)",
                    "up to 500 nearest sequences to a single selected mAb - Parsimony; 50% CDR3 identity"), selectize = FALSE),
      div(style="display: inline-block; width: 300px;",
          sliderInput("height", "Topology height", min = 200, max = 4200, value = 1000)),
      div(style="display: inline-block; width: 300px;",
          sliderInput("width2", "Topology width", min = 0, max = 45, value = 5)),
      # div(style="display: inline-block; width: 300px;",
      #     sliderInput("width", "Topology width2", min = 100, max = 3000, value = 1600)),
      div(HTML("<br>")),br(),
      plotOutput("phyloPlot", inline = TRUE),
      ## new code ot take screenshots from https://deanattali.com/blog/shinyscreenshot-release/
      actionButton("screensht", "Take a screenshot")
    )
  )
)


server <- function(input, output, session) {
  options(width = 180, DT.options = list(pageLength = 10)) # Increase text width for printing table ALSO ADDING DEFAULT NUMBER OF 10 ROWS IN DATATABLE
  ## early in server now defining facets here - this is where all of the sidebar dataset options are tied to a dataset
  
## this only works if you include the if isn't null argument...also another tricky part is below where you specify these if selected (the inputdataset <- reactive({ part...), you need to add parentheses, so uploads1() and uploads2() !!!
  uploads1 <- reactive({
    if (!is.null(input$calfile1)) {
    upload1 <- read.delim(input$calfile1$datapath)
    # upload1$sequence_id <- as.character(upload1$sequence_id)
    # upload1$cdr3_aa_imgt <- as.character(upload1$cdr3_aa_imgt)
    upload1
    }
  })

  uploads2 <- reactive({
    if (!is.null(input$calfile2)) {
      upload2 <- read.delim(input$calfile2$datapath)
      # upload2$sequence_id <- as.character(upload2$sequence_id)
      # upload2$cdr3_aa_imgt <- as.character(upload2$cdr3_aa_imgt)
      upload2
    }
  })
  
## JW idea: add an else if, running the convert function if needed...? but also need to combine datasets, not just convert...
  # convert1 <- NULL
  # convert2 <- NULL
  convert1 <- reactive({
    uploaded_dataset1 <- uploads1()
    if (is.null(uploaded_dataset1$gf_jgene)) {
      upload1afterconv <- shinyprocess(uploaded_dataset1)
      } else {
        upload1afterconv <- uploaded_dataset1
      }
    upload1afterconv$sequence_id <- as.character(upload1afterconv$sequence_id)
    upload1afterconv$cdr3_aa_imgt <- as.character(upload1afterconv$cdr3_aa_imgt)
    upload1afterconv
    })
  convert2 <- reactive({
    uploaded_dataset2 <- uploads2()
    if (is.null(uploaded_dataset2$gf_jgene)) {
      upload2afterconv <- shinyprocess(uploaded_dataset2)
    } else {
      upload2afterconv <- uploaded_dataset2
    }
    upload2afterconv$sequence_id <- as.character(upload2afterconv$sequence_id)
    upload2afterconv$cdr3_aa_imgt <- as.character(upload2afterconv$cdr3_aa_imgt)
    upload2afterconv
  })

  facetvar1 <- reactive({
    switch(input$dataset, "SARS-CoV2 mAbs - heavy chains & light chains" = "cregion", "SARS-CoV2 mAbs - IgH by binding" = "binding", "SARS-CoV2 mAbs - IgH by neutralization" = "neutralization", "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH" = "id", "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH combined" = "cregion", "SARS-CoV2 mAbs vs. HIV mAbs - IgH" = "id", "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined" = "cregion", "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire vs. HIV Setliff 2018 patient bulk repertoires - IgH" = "id", "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire vs. HIV Setliff 2018 patient bulk repertoires - IgH combined" = "cregion", "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH" = "id", "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH combined" = "cregion", "Your datasets - IgH" = "id", "Your datasets - IgH combined" = "cregion")
  })
  ## to change x-axis columns vs. leaving fixed (for IgH only datasets) - note if you leave out default is fixed...
  facetvar2 <- reactive({
    switch(input$dataset, "SARS-CoV2 mAbs - heavy chains & light chains" = "free_x", "SARS-CoV2 mAbs - IgH by binding" = "fixed", "SARS-CoV2 mAbs - IgH by neutralization" = "fixed")
  })
  
  inputdataset <- reactive({
    ## using new names that reduce variables, round
    switch(input$dataset, "SARS-CoV2 mAbs - heavy chains & light chains" = toshiny.cov2.abdab, "SARS-CoV2 mAbs - IgH by binding" = toshiny.cov2.abdab.h, "SARS-CoV2 mAbs - IgH by neutralization" = toshiny.cov2.abdab.h, "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH" = toshiny.cov2.all, "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH combined" = toshiny.cov2.allc, "SARS-CoV2 mAbs vs. HIV mAbs - IgH" = toshiny.cov2hiv, "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined" = toshiny.cov2hivc, "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire vs. HIV Setliff 2018 patient bulk repertoires - IgH" = toshiny.hiv.all, "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire vs. HIV Setliff 2018 patient bulk repertoires - IgH combined" = toshiny.hiv.allc, "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH" = toshiny.den.all, "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH combined" = toshiny.den.allc, "Your datasets - IgH" = convert1(), "Your datasets - IgH combined" = convert2())
  })
  ## extra filter for downloading
  filteredDS <- reactive({
    nearPoints(inputdataset(), input$plot_click, threshold = 5, maxpoints = 99999, addDist = FALSE)
  })
  
  ## isolating the hover x & y for highlighting - note works but makes hover pop up window very momentary
  ## tried to below within renderPlot, doesn't change (trying isolate - same issue, then observe & isolate - breaks)
  highlightedpointsDSx <- reactive({
    dat <- inputdataset()
    highlightedpoint <- nearPoints(dat, input$plot_hover, threshold = 5, maxpoints = 1, addDist = FALSE)
    highlightedpoint$gf_jgene
  })
  highlightedpointsDSy <- reactive({
    dat <- inputdataset()
    highlightedpoint <- nearPoints(dat, input$plot_hover, threshold = 5, maxpoints = 1, addDist = FALSE)
    highlightedpoint$cdr3length_imgt
  })
  
  ### added code to allow for trees of subsetted data:
  # Thanks to @yihui this is now possible using the DT package and input$tableId_rows_all 
  # where tableID is the id assigned to your table. See the link for details. http://rstudio.github.io/DT/shiny.html
  filteredDSpartial <- reactive({
    ids <- input$click_info_rows_selected
    filteredDS()[ids,]
  })
  filteredDSall <- reactive({
    ids <- input$click_info_rows_all
    filteredDS()[ids,]
  })
  
  ## then rename sequence_id
  filteredDSpartial2 <- reactive({
    partial2 <- filteredDSpartial() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt)
    partial2$sequence_id <- paste(partial2$sequence_id,partial2$gene,partial2$cdr3_aa_imgt,sep="_")
    partial2
  })
  
  filteredDSall2 <- reactive({
    all2 <- filteredDSall() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt)
    all2$sequence_id <- paste(all2$sequence_id,all2$gene,all2$cdr3_aa_imgt,sep="_")
    all2
  })
  
  filteredDSpartial2id <- reactive({
    partial2id <- filteredDSpartial() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt)
    partial2id$sequence_id <- paste(partial2id$sequence_id,partial2id$gene,partial2id$cdr3_aa_imgt,sep="_")
    partial2id$sequence_id
  })
  
  ### ADDING NEW VARIABLE SEPARATELY CALCULATING DISTANCE MATRIX AS BELOW BUT WITHIN TREE-PLOTTING STEP
  ## THIS SHOULD SEPARATELY BE AVAILABLE FOR DOWNLOADING...
  matrixDSall2 <- reactive({
    filteredDataall <- filteredDSall2()
    filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
    rownames(filteredDataall.y) <- filteredDataall[,1]
    as.AAbin(filteredDataall.y)
    seqsall = as.phyDat(filteredDataall.y, type = "AA")
    dmalldf = dist.hamming(seqsall, ratio = TRUE)
    dm.matrixall <- as.matrix(dmalldf) %>% as.data.frame()
    dm.matrixall
  })        
  
  ### Hover output: popup window
  output$hover_info <- renderUI({
    hover <- input$plot_hover
    dat <- inputdataset()
    ## first line for picking among different graphs- NOTE CHANGE TO 'DAT'
    point <- nearPoints(dat, input$plot_hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    # tags$style( '#ggplot1 { cursor: pointer; }')
    ## new bit added to change cursor when hovering: (doesn't seem change the cursor though)
    # if (nrow(point) == 0){
    #   css_string <- '
    #   #hover_info {
    #   cursor: default !important; }'
    # } else {
    #   css_string <- '
    #   #hover_info {
    #   cursor: crosshair !important; }'
    # }
    # tags$style(HTML(css_string))
    
    # # calculate point position INSIDE the image as percent of total dimensions
    # # from left (horizontal) and from top (vertical)
    ## simplifiying this per https://gitlab.com/-/snippets/16220
    left_px <- hover$coords_css$x
    top_px <- hover$coords_css$y
    
    # create style property for tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px, "px; top:", top_px + 20, "px;")
    # actual tooltip created as wellPanel - TO ACTUALLY GET COUNTS NEED TO GET FROM STAT_BIN ANALYSIS...
    wellPanel(
      style = style,
      p(HTML(paste0("<b> V-gene & J-gene: </b>", point$gf_jgene, "<br/>",
                    "<b> CDR3 Length (aa): </b>", point$cdr3length_imgt, "<br/>",
                    "<b> Mean Somatic Hypermutation (%): </b>", point$shm_mean, "<br/>",
                    "<b> Max Somatic Hypermutation (%): </b>", point$shm_max, "<br/>",
                    "<b> Count: </b>", point$ncount, "<br/>")))
    )
  })
  
  ### the heatmap (now the only plot)
  output$ggplot1 <- renderPlot({
    ## first line for picking among different graphs
    facet_formula <- as.formula(paste("~", facetvar1()))
    dat <- inputdataset()
    #    note at end of the 3 plot commands below has a final geom_point plotting step to add a green square under the hover, currently it disappears after a moment - tried adding an isolate command, didn't solve
    data2 <- if (input$plotcolors == "Average SHM") {
      toplot <- ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 16) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "Mean \nSHM (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Mean Somatic Hypermutation as fill color") + theme(plot.title = element_text(size = 20, face = "bold")) #+ geom_point(data=subset(dat, gf_jgene ==  highlightedpointsDSx() & cdr3length_imgt == highlightedpointsDSy()), color = "green", shape = "square", size = 4)
    } else if (input$plotcolors == "Maximum SHM") {
      toplot <- ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_max)) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 16) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "Max \nSHM (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Maximum Somatic Hypermutation as fill color") + theme(plot.title = element_text(size = 20, face = "bold")) #+ geom_point(data=subset(dat, gf_jgene ==  highlightedpointsDSx() & cdr3length_imgt == highlightedpointsDSy()), color = "green", shape = "square", size = 4)
    } else {
      toplot <- ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(bins = 40, aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 16) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "% of \nReads  ", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Percentage of total antibody sequences as fill color") + theme(plot.title = element_text(size = 20, face = "bold")) #+ geom_point(data=subset(dat, gf_jgene ==  highlightedpointsDSx() & cdr3length_imgt == highlightedpointsDSy()), color = "green", shape = "square", size = 4)
    }     
    toplot
  }, width = 1200, height = "auto")
  
  
  ### datatable under plot 
  output$click_info <- DT::renderDataTable({
    dat <- inputdataset()
    old <- options(width = 1600); on.exit(options(old))
    nearPoints(dat, input$plot_click, threshold = 5, maxpoints = 99999,
               addDist = FALSE)
  }, width = 1600)
  
  ## adding a datatable with all of the points when double-clicking plot_brush brush_info
  output$brush_info <- DT::renderDataTable({
    dat0 <- inputdataset()
    old <- options(width = 1600); on.exit(options(old))
    brushedPoints(dat0, input$plot_brush, allRows = FALSE)
  }, width = 1600)
  
  
  ## this allows user to download all sequences in clicked bin (first filter) or only selected rows (second filter)
  output$downloadfilter <- downloadHandler(
    filename = function() {
      paste('Filtered data-', Sys.Date(), '.csv', sep = '')
    },
    content = function(file){
      write.csv(filteredDS()[input[["click_info_rows_all"]], ],file, row.names = FALSE)
    }
  )
  
  output$downloadfilter2 <- downloadHandler(
    filename = function() {
      paste('Filtered data-', Sys.Date(), '.csv', sep = '')
    },
    content = function(file){
      write.csv(filteredDS()[input[["click_info_rows_selected"]], ],file, row.names = FALSE)
    }
  )      
  
  ### adding third download button for matrix of distances...
  output$downloadfilter3 <- downloadHandler(
    filename = function() {
      paste('Distance matrix-', Sys.Date(), '.csv', sep = '')
    },
    content = function(file){
      write.csv(matrixDSall2(),file, row.names = FALSE)
    }
  )          
  
  ## alternate phylogenies of CDR3 motifs of selected sequences from datatable
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$go, {
    v$doPlot <- input$go
  })
  
  observeEvent(input$plottab, {
    v$doPlot <- FALSE
  })  
  
  output$phyloPlot <- renderPlot(
    width = 1600,
    # width = function() input$width, ## can comment this width parameter out, only use width2
    height = function() input$height, {
      par(family = "mono", mar=c(0.3, 0, 0.3, 0) + 2.2)
      if (v$doPlot == FALSE) return()
      
      
      isolate({
        data <- if (input$plottab == "NJ") {
          filteredData <- isolate({filteredDSpartial2()})
          filteredData.y <- t(sapply(strsplit(filteredData[,2],""), tolower))
          rownames(filteredData.y) <- filteredData[,1]
          as.AAbin(filteredData.y)
          seqs = as.phyDat(filteredData.y, type = "AA")
          dm = dist.ml(seqs, model="WAG")
          tree = midpoint(NJ(dm))
          tree$edge.length[tree$edge.length<0] <- 0
          tree$edge.length <- tree$edge.length * 0.15
        } else if (input$plottab == "Parsimony") {
          filteredData <- isolate({filteredDSpartial2()})
          filteredData.y <- t(sapply(strsplit(filteredData[,2],""), tolower))
          rownames(filteredData.y) <- filteredData[,1]
          as.AAbin(filteredData.y)
          seqs = as.phyDat(filteredData.y, type = "AA")
          dm = dist.ml(seqs, model="WAG")
          tree1 = NJ(dm)
          tree1$edge.length[tree1$edge.length<0] <- 0
          tree <- pratchet(seqs, start = tree1, maxit=500,
                           minit=10, k=10, trace=0)
          tree <- acctran(tree1, seqs) # added
        } else if (input$plottab == "up to 500 nearest sequences to a single selected mAb - Parsimony; 100% CDR3 identity (Briney 2019)") {
          filteredDataall <- isolate({filteredDSall2()})
          filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
          rownames(filteredDataall.y) <- filteredDataall[,1]
          as.AAbin(filteredDataall.y)
          seqsall = as.phyDat(filteredDataall.y, type = "AA")
          dmall = dist.ml(seqsall, model="WAG")
          ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
          filteredData1ID <- isolate({filteredDSpartial2id()})
          dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(all_of(filteredData1ID)))   ## was arrange_at(1) but this is superceded by arrange(across(1))
          dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
          ### for distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))
          colnames(dmall.matrix)[2] <- "DIST"
          dmall.matrix <- dmall.matrix %>% filter(DIST < 0.02) # now this will change in each else if
          #### end distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id)
          dm.matrix <- dmall.matrix %>% separate(sequence_id, into = c("SEQ", "gene", "cdr3_aa_imgt"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
            select(sequence_id, cdr3_aa_imgt) %>% slice_head(n = 500) ## okay if this is more than what is in the set!
          ## now with new subset make new matrix & tree
          filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
          rownames(filteredData.y) <- dm.matrix[,1]
          as.AAbin(filteredData.y)
          seqs = as.phyDat(filteredData.y, type = "AA")
          dm = dist.ml(seqs, model="WAG")
          tree1 = NJ(dm)
          tree1$edge.length[tree1$edge.length<0] <- 0
          tree <- pratchet(seqs, start = tree1, maxit=500,
                           minit=10, k=10, trace=0)
          tree <- acctran(tree1, seqs) # added
          filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
        } else if (input$plottab == "up to 500 nearest sequences to a single selected mAb - Parsimony; 80% CDR3 identity (Soto 2019)") {
          filteredDataall <- isolate({filteredDSall2()})
          filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
          rownames(filteredDataall.y) <- filteredDataall[,1]
          as.AAbin(filteredDataall.y)
          seqsall = as.phyDat(filteredDataall.y, type = "AA")
          dmall = dist.ml(seqsall, model="WAG")
          ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
          filteredData1ID <- isolate({filteredDSpartial2id()})
          dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(all_of(filteredData1ID)))
          dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
          ### for distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))
          colnames(dmall.matrix)[2] <- "DIST"
          dmall.matrix <- dmall.matrix %>% filter(DIST < 0.301) ## now this will change in each else if
          #### end distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id)
          dm.matrix <- dmall.matrix %>% separate(sequence_id, into = c("SEQ", "gene", "cdr3_aa_imgt"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
            select(sequence_id, cdr3_aa_imgt) %>% slice_head(n = 500) ## okay if this is more than what is in the set!
          ## now with new subset make new matrix & tree
          filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
          rownames(filteredData.y) <- dm.matrix[,1]
          as.AAbin(filteredData.y)
          seqs = as.phyDat(filteredData.y, type = "AA")
          dm = dist.ml(seqs, model="WAG")
          tree1 = NJ(dm)
          tree1$edge.length[tree1$edge.length<0] <- 0
          tree <- pratchet(seqs, start = tree1, maxit=500,
                           minit=10, k=10, trace=0)
          tree <- acctran(tree1, seqs) # added
          filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
        } else if (input$plottab == "up to 500 nearest sequences to a single selected mAb - Parsimony; 70% CDR3 identity (Setliff 2018)") {
          filteredDataall <- isolate({filteredDSall2()})
          filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
          rownames(filteredDataall.y) <- filteredDataall[,1]
          as.AAbin(filteredDataall.y)
          seqsall = as.phyDat(filteredDataall.y, type = "AA")
          dmall = dist.ml(seqsall, model="WAG")
          ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
          filteredData1ID <- isolate({filteredDSpartial2id()})
          dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(all_of(filteredData1ID)))
          dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
          ### for distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))  ## getting this error message Use `all_of(filteredData1ID)` instead of `filteredData1ID` to silence this message.
          colnames(dmall.matrix)[2] <- "DIST"
          dmall.matrix <- dmall.matrix %>% filter(DIST < 0.451) ## now this will change in each else if
          #### end distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id)
          dm.matrix <- dmall.matrix %>% separate(sequence_id, into = c("SEQ", "gene", "cdr3_aa_imgt"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
            select(sequence_id, cdr3_aa_imgt) %>% slice_head(n = 500) ## okay if this is more than what is in the set!
          ## now with new subset make new matrix & tree
          filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
          rownames(filteredData.y) <- dm.matrix[,1]
          as.AAbin(filteredData.y)
          seqs = as.phyDat(filteredData.y, type = "AA")
          dm = dist.ml(seqs, model="WAG")
          tree1 = NJ(dm)
          tree1$edge.length[tree1$edge.length<0] <- 0
          tree <- pratchet(seqs, start = tree1, maxit=500,
                           minit=10, k=10, trace=0)
          tree <- acctran(tree1, seqs) # added
          filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
        } else {
          filteredDataall <- isolate({filteredDSall2()})
          filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
          rownames(filteredDataall.y) <- filteredDataall[,1]
          as.AAbin(filteredDataall.y)
          seqsall = as.phyDat(filteredDataall.y, type = "AA")
          dmall = dist.ml(seqsall, model="WAG")
          ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
          filteredData1ID <- isolate({filteredDSpartial2id()})
          dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(all_of(filteredData1ID)))
          dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
          ### for distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))
          colnames(dmall.matrix)[2] <- "DIST"
          dmall.matrix <- dmall.matrix %>% filter(DIST < 0.751)  ## now this will change in each else if
          #### end distance thresholding
          dmall.matrix <- dmall.matrix %>% select(sequence_id)
          dm.matrix <- dmall.matrix %>% separate(sequence_id, into = c("SEQ", "gene", "cdr3_aa_imgt"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
            select(sequence_id, cdr3_aa_imgt) %>% slice_head(n = 500) ## okay if this is more than what is in the set!
          ## now with new subset make new matrix & tree
          filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
          rownames(filteredData.y) <- dm.matrix[,1]
          as.AAbin(filteredData.y)
          seqs = as.phyDat(filteredData.y, type = "AA")
          dm = dist.ml(seqs, model="WAG")
          tree1 = NJ(dm)
          tree1$edge.length[tree1$edge.length<0] <- 0
          tree <- pratchet(seqs, start = tree1, maxit=500,
                           minit=10, k=10, trace=0)
          tree <- acctran(tree1, seqs) # added, adds edge lengths
          filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
        }
        ### this changes the colors based on source - note these are unique to the datasets...
        tipcolors <- def(tree$tip.label, "hc" = "gray20", "nielsen" = "coral", "galson" = "indianred", "binder" = "orange", "kc" = "sienna", "mt1214" = "limegreen", "nih45" = "seagreen", "bulk-cap" = "darkgreen", "d13" = "blue", "Parameswaran" = "goldenrod", "SARS-CoV2-mAb" = "orchid", "plasmablasts" = "orchid", "denmab" = "orchid", "HIV-IEDBmAb" = "orchid", "HIV-CATNAPmAb" = "orchid", "HIV-Yacoob-mAb" = "orchid",default = "black", regexp = TRUE)
        ### below changing to midpoint of tree
        plot(midpoint(tree), lab4ut="axial",
             edge.width=2, label.offset = 0, cex = 1.2, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors, x.lim = input$width2)  ## x.lim breaks parsimony and some NJ
        add.scale.bar(cex = 1.2, font = 4, lwd = 2)  ## was x = 1, y = -0.1, not sure if this will always work or be meaningful  , x.lim = 15 can also try x.lim = 50 - seems to work only with nearest 50/500 not with NJ or parsimony
        gjcdr3.title <- filteredData
        gjcdr3.title$cdr3length_imgt <- as.numeric(gjcdr3.title$cdr3length_imgt)
        gjcdr3.title$G_J_CDR3 <- paste("Topology of selected CDR3 motifs:", gjcdr3.title$gf_jgene, gjcdr3.title$cdr3length_imgt, "aa", sep=" ")
        gjcdr3.title <- gjcdr3.title %>%
          select(G_J_CDR3)
        gjcdr3.titleID <- gjcdr3.title$G_J_CDR3[1]
        title(gjcdr3.titleID, family = "sans")
      })
    })
  
  observeEvent(input$screensht, {
    screenshot()
  })  
  
}

shinyApp(ui = ui, server = server)

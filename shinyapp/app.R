#install.packages("shiny")
#install.packages("ggplot2")
#install.packages("tidyverse")
library(shiny)
library(ggplot2)
library(tidyverse)
library(DT)

## for plotting
library(seqinr)
library(phangorn)
library(ape)

##########

## dataset loading:
toshiny.cov2.abdab <- read.delim("toshiny_cov2_abdab.tab")
toshiny.cov2.abdab.h <- read.delim("toshiny_cov2_abdab_h.tab")

toshiny.cov2hiv <- read.delim("toshiny_cov2hiv.tab")
toshiny.cov2hivc <- read.delim("toshiny_cov2hivc.tab")

toshiny.cov2.all <- read.delim("toshiny_cov2_all.tab")
toshiny.cov2.allc <- read.delim("toshiny_cov2_allc.tab")

toshiny.hiv.all <- read.delim("toshiny_hiv_all.tab")
toshiny.hiv.allc <- read.delim("toshiny_hiv_allc.tab")

toshiny.den.all <- read.delim("toshiny_den_all.tab")
toshiny.den.allc <- read.delim("toshiny_den_allc.tab")

toshiny.cov2.abdab$sequence_id <- as.character(toshiny.cov2.abdab$sequence_id)
toshiny.cov2.abdab$cdr3_aa_imgt <- as.character(toshiny.cov2.abdab$cdr3_aa_imgt)
toshiny.cov2.abdab.h$sequence_id <- as.character(toshiny.cov2.abdab.h$sequence_id)
toshiny.cov2.abdab.h$cdr3_aa_imgt <- as.character(toshiny.cov2.abdab.h$cdr3_aa_imgt)

toshiny.cov2hiv$sequence_id <- as.character(toshiny.cov2hiv$sequence_id)
toshiny.cov2hiv$cdr3_aa_imgt <- as.character(toshiny.cov2hiv$cdr3_aa_imgt)
toshiny.cov2hivc$sequence_id <- as.character(toshiny.cov2hivc$sequence_id)
toshiny.cov2hivc$cdr3_aa_imgt <- as.character(toshiny.cov2hivc$cdr3_aa_imgt)

toshiny.cov2.all$sequence_id <- as.character(toshiny.cov2.all$sequence_id)
toshiny.cov2.all$cdr3_aa_imgt <- as.character(toshiny.cov2.all$cdr3_aa_imgt)
toshiny.cov2.allc$sequence_id <- as.character(toshiny.cov2.allc$sequence_id)
toshiny.cov2.allc$cdr3_aa_imgt <- as.character(toshiny.cov2.allc$cdr3_aa_imgt)

toshiny.hiv.all$sequence_id <- as.character(toshiny.hiv.all$sequence_id)
toshiny.hiv.all$cdr3_aa_imgt <- as.character(toshiny.hiv.all$cdr3_aa_imgt)
toshiny.hiv.allc$sequence_id <- as.character(toshiny.hiv.allc$sequence_id)
toshiny.hiv.allc$cdr3_aa_imgt <- as.character(toshiny.hiv.allc$cdr3_aa_imgt)

toshiny.den.all$sequence_id <- as.character(toshiny.den.all$sequence_id)
toshiny.den.all$cdr3_aa_imgt <- as.character(toshiny.den.all$cdr3_aa_imgt)
toshiny.den.allc$sequence_id <- as.character(toshiny.den.allc$sequence_id)
toshiny.den.allc$cdr3_aa_imgt <- as.character(toshiny.den.allc$cdr3_aa_imgt)


# toshiny.hcsle <- read.delim("toshiny_hcsle.tsv")
# toshiny.hcsle$sequence_id <- as.character(toshiny.hcsle$sequence_id)
# toshiny.hcsle$cdr3_aa_imgt <- as.character(toshiny.hcsle$cdr3_aa_imgt)
# 

## to re-order any rows for plotting re-run here
toshiny.cov2.abdab.h$binding <- factor(toshiny.cov2.abdab.h$binding, levels = c("RBD", "non-RBD"))
toshiny.cov2.all$id <- factor(toshiny.cov2.all$id, levels = c("SARS-CoV2 mAbs", "DEN patient 13 bulk repertoire", "COVID-19 patient Binder p11 bulk repertoire", "COVID-19 patient Galson p1 bulk repertoire","COVID-19 patient Kuri-Cervantes m5 bulk repertoire", "COVID-19 patient Nielsen p7450 bulk repertoire", "Healthy control bulk repertoire"))
toshiny.cov2hiv$id <- factor(toshiny.cov2hiv$id, levels = c("SARS-CoV2 mAbs", "HIV mAbs"))
toshiny.den.all$id <- factor(toshiny.den.all$id, levels = c("Dengue plasmablasts", "Dengue patient d13 bulk repertoire", "Dengue Parameswaran 2013 patient bulk repertoires"))


## 4 COVID-19 patients vs. Healthy control                                       4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire
## HIV mAbs vs. HIV patient MT1214 vs. HIV patient NIH45                         HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire
## Dengue mAbs vs. Dengue patient d13 vs. Dengue Parameswaran 2013 patients      Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires
ui <- fluidPage(
  
  sidebarLayout(
    
    sidebarPanel(
      # 2 ways to define depending on reactive definitions below
      selectInput("dataset", "Dataset:",
                  c("SARS-CoV2 mAbs - heavy chains & light chains",
                    "SARS-CoV2 mAbs - IgH by binding",
                    "SARS-CoV2 mAbs - IgH by neutralization",
                    "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH",
                    "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH combined",
                    "SARS-CoV2 mAbs vs. HIV mAbs - IgH",
                    "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined",
                    "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire - IgH",
                    "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire - IgH combined",
                    "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH",
                    "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH combined"), selectize = FALSE),
      selectInput("plotcolors", "Plot Colors:",
                  c("Average SHM",
                    "Maximum SHM",
                    "Percentage of total reads"), selectize = FALSE), 
      p("When plots appear, click on a bin to get a list of antibodies in the lower table. Hovering over a bin will also show some basic stats."),
      br(),
      p("Alternately if you want to see more than a bin you can create a box and all antibodies within will appear in the top table."),
      br(),
      p("From the lower table you can download all or selected antibodies in the chosen bin, download the distance matrix of all antibodies, or create networks of selected antibodies. The last network options are to find the nearest sequences (up to 500) of a single selected antibody, with 4 possible distance thresholds. Note that you can change the window size of the network using the slider."),
      br(),
      p("Finally make sure to check all antibodies in the table have the same CDR3 length or the distance matrix/tree building will fail."),
      width = 2
      
    ),

    mainPanel(
      plotOutput("ggplot1", width = "120%", height = "800px", hover = "plot_hover", click = "plot_click", brush = "plot_brush"),
      uiOutput("hover_info"),
      # plotOutput("ggplot2", height = "800px", click = "plot_click"),
      # uiOutput("hover2_info"),
#      verbatimTextOutput("click_info") - instead trying datatable, next line
## new datatable if double-clicked
      DT::dataTableOutput("brush_info"),
      DT::dataTableOutput("click_info"),
      actionButton("go", "Make network of selected CDR3 AA motifs"), 
      downloadButton("downloadfilter","Download all data in the clicked bin"),
      downloadButton("downloadfilter2","Download only selected rows"),
      downloadButton("downloadfilter3","Download distance matrix of all data in the clicked bin"),
      selectInput("plottab", "Network:",
            c("NJ",
              "Parsimony",
              "up to 500 nearest sequences to a single selected mAb - Parsimony; 100% CDR3 identity (Briney 2019)",
              "up to 500 nearest sequences to a single selected mAb - Parsimony; 80% CDR3 identity (Soto 2019)",
              "up to 500 nearest sequences to a single selected mAb - Parsimony; 75% CDR3 identity",
              "up to 500 nearest sequences to a single selected mAb - Parsimony; 50% CDR3 identity"), selectize = FALSE),
      div(style="display: inline-block; width: 300px;",
          sliderInput("height", "Network height", min = 200, max = 4200, value = 1000)),
      div(style="display: inline-block; width: 300px;",
          sliderInput("width", "Network width", min = 1000, max = 3000, value = 1600)),
      div(HTML("<br>")),br(),
      plotOutput("phyloPlot", inline = TRUE),
#plotOutput("phyloPlot", width = "120%", height = "1700px", inline = FALSE),
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
      switch(input$dataset, "SARS-CoV2 mAbs - heavy chains & light chains" = "cregion", "SARS-CoV2 mAbs - IgH by binding" = "binding", "SARS-CoV2 mAbs - IgH by neutralization" = "neutralization", "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH" = "id", "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH combined" = "cregion", "SARS-CoV2 mAbs vs. HIV mAbs - IgH" = "id", "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined" = "cregion", "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire - IgH" = "id", "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire - IgH combined" = "cregion", "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH" = "id", "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH combined" = "cregion")
    })
    ## to change x columns vs. leaving fixed (for IgH only datasets) - note if I leave out default is fixed...
    facetvar2 <- reactive({
      switch(input$dataset, "SARS-CoV2 mAbs - heavy chains & light chains" = "free_x", "SARS-CoV2 mAbs - IgH by binding" = "fixed", "SARS-CoV2 mAbs - IgH by neutralization" = "fixed")
    })
    
  ## alternate to reactive + get from Shiny advanced example:
  inputdataset <- reactive({
#    switch(input$dataset, "COVID mAbs" = BX.clusters.mabs, "Comet Data" = biohub.all2, "Healthy Control" = BX.clusters.hc, "Healthy Control vs. Comet Data vs. COVID mAbs" = toshiny.few)
## using new names that reduce variables, round
    switch(input$dataset, "SARS-CoV2 mAbs - heavy chains & light chains" = toshiny.cov2.abdab, "SARS-CoV2 mAbs - IgH by binding" = toshiny.cov2.abdab.h, "SARS-CoV2 mAbs - IgH by neutralization" = toshiny.cov2.abdab.h, "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH" = toshiny.cov2.all, "SARS-CoV2 mAbs vs. 4 COVID-19 patient bulk repertoires vs. Healthy control bulk repertoire - IgH combined" = toshiny.cov2.allc, "SARS-CoV2 mAbs vs. HIV mAbs - IgH" = toshiny.cov2hiv, "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined" = toshiny.cov2hivc, "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire - IgH" = toshiny.hiv.all, "HIV mAbs vs. HIV patient MT1214 bulk repertoire vs. HIV patient NIH45 bulk repertoire - IgH combined" = toshiny.hiv.allc, "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH" = toshiny.den.all, "Dengue mAbs vs. Dengue patient d13 bulk repertoire vs. Dengue Parameswaran 2013 patient bulk repertoires - IgH combined" = toshiny.den.allc)
#    switch(input$dataset, "COVID mAbs" = toshiny.few1, "COVID mAbs - IgH by binding" = toshiny.few1.h, "COVID mAbs - IgH by neutralization" = toshiny.few1.h, "Comet Data" = toshiny.few2, "Healthy Control" = toshiny.few3b, "Healthy Control vs. COVID mAbs - IgH" = toshiny.few, "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH" = toshiny.fewhcbg, "HC vs. COVID Boyd7450 vs. COVID Galson1 vs. COVID mAbs - IgH combined" = toshiny.fewhcbgc)
    
    
    
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
#    filteredDSpartial$sequence_id <- paste(filteredDSpartial$sequence_id,filteredDSpartial$cdr3_aa_imgt,filteredDSpartial$gene,sep="_") 
#    filteredDSall$sequence_id <- paste(filteredDSall$sequence_id,filteredDSall$cdr3_aa_imgt,filteredDSall$gene,sep="_") 

    ## these 2 are not working correctly - trying slightly different 8/31 THIS NOW WORKS...
    ## EXCEPT 9/1/2020  the as.character part seems to be catching in the plot function steps below...SOLVED BY RUNNING AS.CHARACTER AT VERY START
    filteredDSpartial2 <- reactive({
      partial2 <- filteredDSpartial() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt)
#      partial2 <- filteredDSpartial() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt) %>% mutate_if(is.factor, as.character)
      partial2$sequence_id <- paste(partial2$sequence_id,partial2$gene,partial2$cdr3_aa_imgt,sep="_")
      # partial2$sequence_id <- as.character(partial2$sequence_id)
      # partial2$cdr3_aa_imgt <- as.character(partial2$cdr3_aa_imgt)
      partial2
    })
    
    filteredDSall2 <- reactive({
      all2 <- filteredDSall() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt)
#      all2 <- filteredDSall() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt) %>% mutate_if(is.factor, as.character)
      all2$sequence_id <- paste(all2$sequence_id,all2$gene,all2$cdr3_aa_imgt,sep="_")
      all2
    })

    filteredDSpartial2id <- reactive({
      partial2id <- filteredDSpartial() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt)
      #      partial2 <- filteredDSpartial() %>% select(sequence_id,cdr3_aa_imgt,gene,gf_jgene,cdr3length_imgt) %>% mutate_if(is.factor, as.character)
      partial2id$sequence_id <- paste(partial2id$sequence_id,partial2id$gene,partial2id$cdr3_aa_imgt,sep="_")
      ## error here because I really only want a string from this output
      partial2id$sequence_id
      # partial2id$cdr3_aa_imgt <- NULL
      # partial2id$gene <- NULL
      # partial2id$gf_jgene <- NULL
      # partial2id$cdr3length_imgt <- NULL
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
    })        

## adding distance treshold as separate reactive variable
    distthreshold <- reactive({
      distthresholdvalue <- input$distthreshold
      distthresholdvalue
    })        
    
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
      p(HTML(paste0("<b> V-gene & J-gene: </b>", point$gf_jgene, "<br/>",
                    "<b> CDR3 Length (aa): </b>", point$cdr3length_imgt, "<br/>",
                    "<b> Mean Somatic Hypermutation (%): </b>", point$shm_mean, "<br/>",
                    "<b> Max Somatic Hypermutation (%): </b>", point$shm_max, "<br/>",
                    "<b> Count: </b>", point$ncount, "<br/>")))
    )
  })
  

### second plot code    NOW ONLY PLOT
    output$ggplot1 <- renderPlot({
      ## first line for picking among different graphs
#      ggplot(inputdataset(), aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill= (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ get(input$facet_by), ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
## if manually choosing input  
#      isolate({
      facet_formula <- as.formula(paste("~", facetvar1()))
      dat <- inputdataset()
#      ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(4, 35)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with % of total reads as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      data2 <- if (input$plotcolors == "Average SHM") {
        toplot <- ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "Mean \nSHM (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Mean Somatic Hypermutation as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      } else if (input$plotcolors == "Maximum SHM") {
        toplot <- ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_max)) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "Max \nSHM (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Maximum Somatic Hypermutation as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      } else {
        toplot <- ggplot(dat, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(bins = 40, aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(facet_formula, ncol = 1, scales = facetvar2()) + scale_fill_viridis_c(name = "% of \nReads  ", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8)) + ggtitle("Bins of V+J gene families vs. CDR3 length, with Percentage of Total Reads as fill color") + theme(plot.title = element_text(size = 16, face = "bold"))
      }     
#      })
      toplot
    }, width = 1200, height = "auto")

#      output$click_info <- renderPrint({   # instead line below

    
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
                              
        #      nearPoints(h3.cometvsmabsvshc.h2.stat$count, xvar = BX.clusters.cometvsmabsvshc.h$gf_jgene, yvar = BX.clusters.cometvsmabsvshc.h$cdr3length_imgt, input$plot_click, threshold = 5, maxpoints = 25,
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
      
      output$phyloPlot <- renderPlot(
        height = function() input$height,
        width = function() input$width, {
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
            tree <- acctran(tree1, seqs) # adding Aug 2021
          } else if (input$plottab == "up to 500 nearest sequences to a single selected mAb - Parsimony; 100% CDR3 identity (Briney 2019)") {
            filteredDataall <- isolate({filteredDSall2()})
            filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
            rownames(filteredDataall.y) <- filteredDataall[,1]
            as.AAbin(filteredDataall.y)
            seqsall = as.phyDat(filteredDataall.y, type = "AA")
            dmall = dist.ml(seqsall, model="WAG")
            ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
            filteredData1ID <- isolate({filteredDSpartial2id()})
            dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(filteredData1ID))   ## was arrange_at(1) but this is superceded by arrange(across(1))
            dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
            ### for distance thresholding
            ## MAR 2021 INSTEAD TRY THESE LINES - THESE WORK:
            dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))
            colnames(dmall.matrix)[2] <- "DIST"
            # dmall.matrix <- rename.vars(dmall.matrix, c("sequence_id","filteredData1ID"), c("sequence_id","DIST"))
            dmall.matrix <- dmall.matrix %>% filter(DIST < 0.02) ## was .505, also tried 0.8 - finally 0.75 is about 50% of aa; trying distthresholdvalue but doesn't seem to work
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
            tree <- acctran(tree1, seqs) # adding Aug 2021
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
            dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(filteredData1ID))
            dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
            ### for distance thresholding
            ## MAR 2021 INSTEAD TRY THESE LINES - THESE WORK:
            dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))
            colnames(dmall.matrix)[2] <- "DIST"
            # dmall.matrix <- rename.vars(dmall.matrix, c("sequence_id","filteredData1ID"), c("sequence_id","DIST"))
            dmall.matrix <- dmall.matrix %>% filter(DIST < 0.301) ## was .505, also tried 0.8 - finally 0.75 is about 50% of aa; trying distthresholdvalue but doesn't seem to work
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
            tree <- acctran(tree1, seqs) # adding Aug 2021
            filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
          } else if (input$plottab == "up to 500 nearest sequences to a single selected mAb - Parsimony; 75% CDR3 identity") {
            filteredDataall <- isolate({filteredDSall2()})
            filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
            rownames(filteredDataall.y) <- filteredDataall[,1]
            as.AAbin(filteredDataall.y)
            seqsall = as.phyDat(filteredDataall.y, type = "AA")
            dmall = dist.ml(seqsall, model="WAG")
            ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
            filteredData1ID <- isolate({filteredDSpartial2id()})
            dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(filteredData1ID))
            dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
            ### for distance thresholding
            ## MAR 2021 INSTEAD TRY THESE LINES - THESE WORK:
            dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))  ## getting this error message Use `all_of(filteredData1ID)` instead of `filteredData1ID` to silence this message.
            colnames(dmall.matrix)[2] <- "DIST"
            # dmall.matrix <- rename.vars(dmall.matrix, c("sequence_id","filteredData1ID"), c("sequence_id","DIST"))
            dmall.matrix <- dmall.matrix %>% filter(DIST < 0.376) ## was .505, also tried 0.8 - finally 0.75 is about 50% of aa; trying distthresholdvalue but doesn't seem to work
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
            tree <- acctran(tree1, seqs) # adding Aug 2021
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
            dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange(across(filteredData1ID))
            dmall.matrix <- rownames_to_column(dmall.matrix, var = "sequence_id")
            ### for distance thresholding
            ## MAR 2021 INSTEAD TRY THESE LINES - THESE WORK:
            dmall.matrix <- dmall.matrix %>% select(sequence_id, all_of(filteredData1ID))
            colnames(dmall.matrix)[2] <- "DIST"
            # dmall.matrix <- rename.vars(dmall.matrix, c("sequence_id","filteredData1ID"), c("sequence_id","DIST"))
            dmall.matrix <- dmall.matrix %>% filter(DIST < 0.751)  ## was 0.66, also tried 0.8  #also sept21 tried 0.5 ..go with 0.75? trying distthresholdvalue but doesn't seem to work
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
            tree <- acctran(tree1, seqs) # adding Aug 2021, adds edge lengths
            filteredData <- isolate({filteredDSpartial2()})  ### added because not in this option but now need for changing title below
          }
          # par(family = "mono", mar=c(2, 1, 2, 0) + 0.5) most recently mar=c(2, 0.0, 2, 0) + 0.8) ## think I need to remove the + bit
          # par(family = "mono", mai=c(1.02, 0.82, 0.82, 0.82))  ## moving to start of command
          ### this changes the colors of all 'controls' to gray but leaves covid mabs black
          tipcolors <- def(tree$tip.label, "hc" = "gray70", "nielsen" = "coral", "galson" = "indianred", "binder" = "orange", "kc" = "sienna", "mt1214" = "plum", "nih45" = "pink", "d13" = "blue", "Parameswaran" = "gold", "SARS-CoV2-mAb" = "orchid", "plasmablasts" = "orchid", "denmab" = "orchid","HIV-mAb" = "orchid", default = "black", regexp = TRUE)
          # plot(bird.orders, tip.color = co2)
          ### below changing to midpoint of tree Aug 2021
          plot(midpoint(tree), lab4ut="axial",
               edge.width=2, label.offset = 0, cex = 0.95, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors)  ## x.lim breaks parsimony and some NJ
          add.scale.bar(cex = 1, font = 4, lwd = 2)  ## was x = 1, y = -0.1, not sure if this will always work or be meaningful  , x.lim = 15 can also try x.lim = 50 - seems to work only with nearest 50/500 not with NJ or parsimony
          #               edge.width=2, label.offset = 0, cex = 0.95, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
          gjcdr3.title <- filteredData
          gjcdr3.title$cdr3length_imgt <- as.numeric(gjcdr3.title$cdr3length_imgt)
          # gjcdr3.title$G_J_CDR3 <- paste(gjcdr3.title$gf_jgene,gjcdr3.title$cdr3length_imgt,sep="_")
          gjcdr3.title$G_J_CDR3 <- paste("Network of selected CDR3 motifs:", gjcdr3.title$gf_jgene, gjcdr3.title$cdr3length_imgt, "aa", sep=" ")
          gjcdr3.title <- gjcdr3.title %>%
            select(G_J_CDR3)
          gjcdr3.titleID <- gjcdr3.title$G_J_CDR3[1]
          title(gjcdr3.titleID, family = "sans")
        })
      })
      
}

shinyApp(ui = ui, server = server)

library(alakazam)
library(igraph)
library(dplyr)
library(RColorBrewer)
library(hexbin)
library(scales)
library(grid)
library(lattice)
library(gdata)
library(gridExtra)
library(ape)
library(shazam)
library(limma)
library(reshape2)
library(DT)

library(ggplot2)
library(seqinr)
library(phangorn)
library(shiny)
library(tidyverse)

## 2 functions
#parse an attribute from tree tip labels based on their position after the 
#last "sep" (default = "_") in the tip name
parseAttribute = function(graph,label,position=0,unknown="unknown",
                          germline="GERM",maskgermline=FALSE,sep="_"){
  if(class(graph) == "list"){
    graphs = lapply(graph,function(x)
      parseAttribute(x,label,position,germline))
    return(graphs)
  }else{
    tips = V(graph)$name
    st = strsplit(tips,split=sep)
    state = unlist(lapply(st,function(x)
      if(length(x) > 1){
        return(x[length(x)-position])
      }else{
        return(unknown)
      }))
    if(maskgermline){
      state[grepl(germline,tips)] = unknown
    }
    vertex_attr(graph, label) = state
    return(graph)
  }
}

#combine multiple graphs into a single igraph object
#changes all node names
combineGraphs = function(graphs,directed=FALSE){
  gu = data.frame()
  attributes = data.frame()
  for(i in 1:length(graphs)){
    g = graphs[[i]]
    E(g)$weight = rep(1,length(E(g)$weight))
    V(g)$name = paste0(i,"-",V(g)$name)
    attributes = 
      rbind(attributes,data.frame(vertex_attr(g),stringsAsFactors=FALSE))
    gu = rbind(gu,igraph::as_data_frame(g))
  }
  graph = graph_from_data_frame(gu,directed=directed)
  oa = attributes[order(match(attributes$name,V(graph)$name)),]
  for(n in names(oa)){
    vertex_attr(graph, n) = oa[,n]
  }
  return(graph)
}


#################################
### COMMANDS FOR UPDATING COVID ABDAB DATABASE BY JUST GRABBING RAW DATA...


#### NEED TO RELOAD ALL TAB FILES BEFORE UPDATING

toshiny.few <- read.delim("toshiny_few.tab")
toshiny.few1 <- read.delim("toshiny_few1.tab")
toshiny.few1.h <- read.delim("toshiny_few1_h.tab")
## realoading few2 - will want to combine COMET for internal comparisons (moving to app_test in _2 folder)
toshiny.few2 <- read.delim("toshiny_few2b.tab")

toshiny.few3 <- read.delim("toshiny_few3.tab")
# toshiny.few3b <- read.delim("toshiny_few3b.tab")
toshiny.fewhcbg <- read.delim("toshiny_fewhcbg.tab")
toshiny.fewhcbgc <- read.delim("toshiny_fewhcbgc.tab")
# toshiny.natalia4 <- read.delim("toshiny_natalia.tab")
toshiny.cov2sarsmers <- read.delim("toshiny_cov2sarsmers.tab")
toshiny.cov2sarsmersc <- read.delim("toshiny_cov2sarsmersc.tab")

toshiny.hivcov2sarsmers <- read.delim("toshiny_hivcov2sarsmers.tab")
toshiny.hivcov2sarsmersc <- read.delim("toshiny_hivcov2sarsmersc.tab")

toshiny.fewhivmabsnih45rd1214 <- read.delim("toshiny_mabsnih45rd1214.tab")
toshiny.fewhivmabsnih45rd1214c <- read.delim("toshiny_mabsnih45rd1214c.tab")

toshiny.dengueall <- read.delim("toshiny_dengueall.tab")
toshiny.dengueallc <- read.delim("toshiny_dengueallc.tab")

toshiny.fluhcmabs <- read.delim("toshiny_fluhcmabs.tab")
toshiny.fluhcmabsc <- read.delim("toshiny_fluhcmabsc.tab")

toshiny.fewg1 <- read.delim("toshiny_fewg1.tab")
toshiny.fewb1 <- read.delim("toshiny_fewb1.tab")
toshiny.few3.h <- subset(toshiny.few3, CREGION %in% c("IgH"))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
###################################################################################################

## FIRST TAKE DATA AND ADD TO GOOGLE DRIVE SPREADSHEET
## THEN RUN A JOIN COMMAND - better yet an rbind if all columns are fixed

### latest version of this is from end of October - starting in early November adding directly to toshiny

#BX.clusters.mabs.toadd.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_111020add.tab")
#BX.clusters.mabs.toadd.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_251120add.tab")
BX.clusters.mabs.toadd.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_060121add.tab")
BX.clusters.mabs.toadd.l <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_060121add.tab")


BX.clusters.mabs.toadd.h$binding <- gsub("S; ","",BX.clusters.mabs.toadd.h$binding)
BX.clusters.mabs.toadd.h$neutralization <- gsub("non-neutralizingunknown","non-neutralizing or unknown",BX.clusters.mabs.toadd.h$neutralization)

BX.clusters.mabs.toadd.h$SEQUENCE_IDLC <- NULL
BX.clusters.mabs.toadd.h$GENELC <- NULL
BX.clusters.mabs.toadd.h$JGENELC <- NULL

BX.clusters.mabs.toadd.h$JUNCTIONAALC <- NULL
BX.clusters.mabs.toadd.h$V_IDENTITYLC <- NULL

BX.clusters.mabs.toadd.h$GENE <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.h$GENE)
BX.clusters.mabs.toadd.h$JGENE <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.h$JGENE)



BX.clusters.mabs.toadd.l$binding <- gsub("S; ","",BX.clusters.mabs.toadd.l$binding)
BX.clusters.mabs.toadd.h$neutralization <- gsub("non-neutralizingunknown","non-neutralizing or unknown",BX.clusters.mabs.toadd.h$neutralization)

BX.clusters.mabs.toadd.l$SEQUENCE_ID <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, SEQUENCE_ID = SEQUENCE_IDLC)


BX.clusters.mabs.toadd.l$GENE <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, GENE = GENELC)
BX.clusters.mabs.toadd.l$JGENE <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, JGENE = JGENELC)

BX.clusters.mabs.toadd.l$JUNCTIONAA <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, JUNCTIONAA = JUNCTIONAALC)
BX.clusters.mabs.toadd.l$V_IDENTITY <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, V_IDENTITY = V_IDENTITYLC)

BX.clusters.mabs.toadd.l$GENE <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.l$GENE)
BX.clusters.mabs.toadd.l$JGENE <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.l$JGENE)

BX.clusters.mabs.toadd <- rbind(BX.clusters.mabs.toadd.h, BX.clusters.mabs.toadd.l)

## THIS NEXT COMMAND MIGHT NOT BE THE CASE FOR ALL - SHOULD CHANGE IN .TAB FILE
#BX.clusters.mabs.toadd$neutralization <- "neutralizing"

##########

BX.clusters.mabs.toadd$JUNCTIONAA <- as.character(BX.clusters.mabs.toadd$JUNCTIONAA)
BX.clusters.mabs.toadd$SEQUENCE_ID <- as.character(BX.clusters.mabs.toadd$SEQUENCE_ID)

BX.clusters.mabs.toadd$CDR3LENGTH_IMGT <- nchar(BX.clusters.mabs.toadd$JUNCTIONAA)

BX.clusters.mabs.toadd$JUNCTION_LENGTH <- (3 * (BX.clusters.mabs.toadd$CDR3LENGTH_IMGT)) + 6
BX.clusters.mabs.toadd$JGF <- BX.clusters.mabs.toadd$JGENE
BX.clusters.mabs.toadd$GF <- str_sub(BX.clusters.mabs.toadd$GENE, end=5)
### note updated str_sub for substring (think both work)
# BX.clusters.mabs.toadd$GF <- substring(BX.clusters.mabs.toadd$GENE, 1,5)
# BX.clusters.mabs.toadd$JGF <- substring(BX.clusters.mabs.toadd$JGENE, 1,5)
BX.clusters.mabs.toadd$VJ_JUNCTION_PATTERN <- paste(BX.clusters.mabs.toadd$GF,BX.clusters.mabs.toadd$JGF,BX.clusters.mabs.toadd$JUNCTION_LENGTH,sep="_") 

### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
BX.clusters.mabs.toadd$SHM <- (100 - (BX.clusters.mabs.toadd$V_IDENTITY))
# BX.clusters.mabs.toadd$SHM <- (100 - (BX.clusters.mabs.toadd$V_IDENTITY * 100))
BX.clusters.mabs.toadd <- BX.clusters.mabs.toadd %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

BX.clusters.mabs.toadd$JUNCTION_LENGTH <- NULL
BX.clusters.mabs.toadd$V_IDENTITY <- NULL

## ADDING CREGION BY WAY OF GF...
BX.clusters.mabs.toadd$CREGION <- str_sub(BX.clusters.mabs.toadd$GF, end=3)
BX.clusters.mabs.toadd$CREGION <- gsub("IGK","Kappa",BX.clusters.mabs.toadd$CREGION)
BX.clusters.mabs.toadd$CREGION <- gsub("IGL","Lambda",BX.clusters.mabs.toadd$CREGION)
BX.clusters.mabs.toadd$CREGION <- gsub("IGH","IgH",BX.clusters.mabs.toadd$CREGION)



### NOV 2020 - GOING TO REMOVE AUGMENTA FROM ALL DATASETS FROM NOW ON...
### manually removing from nov12 version and reloading  - for duplicates jan21 need to manually remove and reload..
toshiny.few1 <- read.delim("toshiny_few1.tab")

BX.clusters.mabs.alladded <- full_join(toshiny.few1, BX.clusters.mabs.toadd)
### now this can replace toshiny.few1 
rm(BX.clusters.mabs.toadd)
rm(BX.clusters.mabs.toadd.h)
rm(BX.clusters.mabs.toadd.l)

#BX.clusters.fewermabs <- BX.clusters.mabs %>% filter(source != "Augmenta")

toshiny.few1$binding <- gsub("Unk","unknown",toshiny.few1$binding)

toshiny.few1 <- BX.clusters.mabs.alladded %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

## NOW NEED TO RECALCULATE NCOUNT AND SHM.MEAN
toshiny.few1$ncount <- NULL
toshiny.few1$shm.mean <- NULL

toshiny.few1 <- toshiny.few1 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.few1.h <- subset(toshiny.few1, CREGION %in% c("IgH"))
# toshiny.few1.h$binding <- gsub("Unk","unknown",toshiny.few1.h$binding)

#write.table(toshiny.few1, "toshiny_few1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.few1.h, "toshiny_few1_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

rm(BX.clusters.mabs.alladded)

## NOTE THAT THIS WILL NOT HAVE SHM DATA (UNLESS FOUND IN A PAPER) - WOULD HAVE TO SEPARATELY GO BACK AND FIND RAW NUCLEOTIDE SEQUENCE IN NCBI TO CALCULATE SHM...

## removing all Augmenta from loaded datasets...
# BX.clusters.mabs <- BX.clusters.mabs %>% filter(source != "Augmenta")
# BX.clusters.mabs.h <- BX.clusters.mabs.h %>% filter(source != "Augmenta")
# BX.clusters.mabsvshc <- BX.clusters.mabsvshc %>% filter(source != "Augmenta")
# BX.clusters.mabsvshc.h <- BX.clusters.mabsvshc.h %>% filter(source != "Augmenta")


#######################################
### NEXT: re-calculating all of the combinations...

### 1 toshiny.few
## note HC only is few3 aka .hc - few3b is HC but all reads not just clones
toshiny.few3.h <- subset(toshiny.few3, CREGION %in% c("IgH"))

toshiny.few <- bind_rows(toshiny.few1.h, toshiny.few3.h, .id = "id")

toshiny.few$id <- gsub("2","Healthy control bulk repertoire",toshiny.few$id)
toshiny.few$id <- gsub("1","anti-CoV2 mAbs",toshiny.few$id)

#write.table(toshiny.few, "toshiny_few.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### 2

BX.clusters.hcvsboydvsgalsonmabs <- bind_rows(toshiny.few3.h, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, .id = "id")
#toshiny.fewhcbg <- subset(BX.clusters.hcvsboydvsgalsonmabs, CREGION %in% c("IgH"))
toshiny.fewhcbg <- BX.clusters.hcvsboydvsgalsonmabs[ grep("Kappa|Lambda", BX.clusters.hcvsboydvsgalsonmabs$CREGION, invert = TRUE) , ] %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
#rm(BX.clusters.hcvsboydvsgalsonmabs)

toshiny.fewhcbg$id <- gsub("1","Healthy control",toshiny.fewhcbg$id)
## be careful with 1 in name here....
toshiny.fewhcbg$id <- gsub("4","anti-CoVII mAbs",toshiny.fewhcbg$id)
toshiny.fewhcbg$id <- gsub("2","Boyd7450",toshiny.fewhcbg$id)
toshiny.fewhcbg$id <- gsub("3","Galson1",toshiny.fewhcbg$id)
toshiny.fewhcbg$id <- gsub("anti-CoVII mAbs","anti-CoV2 mAbs",toshiny.fewhcbg$id)

## for combined need to recalculate ncount and shm.mean
toshiny.fewhcbgc <- toshiny.fewhcbg
toshiny.fewhcbgc$ncount <- NULL
toshiny.fewhcbgc$shm.mean <- NULL

toshiny.fewhcbgc <- toshiny.fewhcbgc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.fewhcbg, "toshiny_fewhcbg.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhcbgc, "toshiny_fewhcbgc.tab", sep = "\t", row.names = FALSE, quote = FALSE)



### 3

BX.clusters.hcvsboydvsgalsonmabshiv <- bind_rows(toshiny.few3.h, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, toshiny.fewhiv, .id = "id")
#toshiny.fewhcbghiv <- subset(BX.clusters.hcvsboydvsgalsonmabshiv, CREGION %in% c("IgH"))
toshiny.fewhcbghiv <- BX.clusters.hcvsboydvsgalsonmabshiv[ grep("Kappa|Lambda", BX.clusters.hcvsboydvsgalsonmabshiv$CREGION, invert = TRUE) , ] %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
#rm(BX.clusters.hcvsboydvsgalsonmabshiv)

toshiny.fewhcbghiv$id <- gsub("1","Healthy control",toshiny.fewhcbghiv$id)
## be careful with 1 in name here....
toshiny.fewhcbghiv$id <- gsub("4","anti-CoVII mAbs",toshiny.fewhcbghiv$id)
toshiny.fewhcbghiv$id <- gsub("2","Boyd7450",toshiny.fewhcbghiv$id)
toshiny.fewhcbghiv$id <- gsub("3","Galson1",toshiny.fewhcbghiv$id)
toshiny.fewhcbghiv$id <- gsub("anti-CoVII mAbs","anti-CoV2 mAbs",toshiny.fewhcbghiv$id)
toshiny.fewhcbghiv$id <- gsub("Boyd7450","Boyd74s0",toshiny.fewhcbghiv$id)
toshiny.fewhcbghiv$id <- gsub("5","HIV+ patient",toshiny.fewhcbghiv$id)
toshiny.fewhcbghiv$id <- gsub("Boyd74s0","Boyd7450",toshiny.fewhcbghiv$id)

## for combined need to recalculate ncount and shm.mean
toshiny.fewhcbghivc <- toshiny.fewhcbghiv
toshiny.fewhcbghivc$ncount <- NULL
toshiny.fewhcbghivc$shm.mean <- NULL

toshiny.fewhcbghivc <- toshiny.fewhcbghivc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

## updaating HIV+ patient to HIV+ patient AD1214
toshiny.fewhcbghiv$id <- gsub("HIV+ patient","HIV+ patient AD1214", toshiny.fewhcbghiv$id)
toshiny.fewhcbghivc$id <- gsub("HIV+ patient","HIV+ patient AD1214", toshiny.fewhcbghivc$id)


#write.table(toshiny.fewhcbghiv, "toshiny_fewhcbghiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhcbghivc, "toshiny_fewhcbghivc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

rm(toshiny.fewhcbghiv)
rm(toshiny.fewhcbghivc)
rm(toshiny.fewhcbg)
rm(toshiny.fewhcbgc)
### 4

toshiny.cov2sarsmers <- bind_rows(toshiny.few1.h, toshiny.sars.h, toshiny.mers.h, .id = "id")

toshiny.cov2sarsmers$id <- gsub("2","anti-SARS mAbs",toshiny.cov2sarsmers$id)
toshiny.cov2sarsmers$id <- gsub("3","anti-MERS mAbs",toshiny.cov2sarsmers$id)
## because has 2 needs to be last
toshiny.cov2sarsmers$id <- gsub("1","anti-CoV2 mAbs",toshiny.cov2sarsmers$id)

toshiny.cov2sarsmers$id <- factor(toshiny.cov2sarsmers$id, levels = c("anti-CoV2 mAbs", "anti-SARS mAbs", "anti-MERS mAbs"))
toshiny.cov2sarsmers$source <- NULL

## for combined need to recalculate ncount and shm.mean
toshiny.cov2sarsmersc <- toshiny.cov2sarsmers
toshiny.cov2sarsmersc$ncount <- NULL
toshiny.cov2sarsmersc$shm.mean <- NULL

toshiny.cov2sarsmersc <- toshiny.cov2sarsmersc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.cov2sarsmers, "toshiny_cov2sarsmers.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.cov2sarsmersc, "toshiny_cov2sarsmersc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### 5

toshiny.hivcov2sarsmers <- bind_rows(toshiny.few1.h, toshiny.sars.h, toshiny.mers.h, toshiny.fewhivmabs.h, .id = "id")

toshiny.hivcov2sarsmers$id <- gsub("2","anti-SARS mAbs",toshiny.hivcov2sarsmers$id)
toshiny.hivcov2sarsmers$id <- gsub("3","anti-MERS mAbs",toshiny.hivcov2sarsmers$id)
toshiny.hivcov2sarsmers$id <- gsub("4","anti-HIV mAbs",toshiny.hivcov2sarsmers$id)
## because has 2 needs to be last
toshiny.hivcov2sarsmers$id <- gsub("1","anti-CoV2 mAbs",toshiny.hivcov2sarsmers$id)

toshiny.hivcov2sarsmers$id <- factor(toshiny.hivcov2sarsmers$id, levels = c("anti-CoV2 mAbs", "anti-SARS mAbs", "anti-MERS mAbs", "anti-HIV mAbs"))
toshiny.hivcov2sarsmers$source <- NULL

## for combined need to recalculate ncount and shm.mean
toshiny.hivcov2sarsmersc <- toshiny.hivcov2sarsmers
toshiny.hivcov2sarsmersc$ncount <- NULL
toshiny.hivcov2sarsmersc$shm.mean <- NULL

toshiny.hivcov2sarsmersc <- toshiny.hivcov2sarsmersc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.hivcov2sarsmers, "toshiny_hivcov2sarsmers.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.hivcov2sarsmersc, "toshiny_hivcov2sarsmersc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## then putting into app-6

#install.packages("rsconnect")
library(rsconnect)
rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-6')


#### jan 2021 ideas:
## change plot title to name of V+J+CDR3length?? key is here:
#filteredData <- isolate({filteredDSpartial2()}) - note this line is in 2 of 3 options, would need to add to third
## then near bottom of plot code something similar to this line to add to actual title:
# gjcdr3.title <- filteredData %>% separate(SEQUENCE_ID, into = c("GF_JGENE", "CDR3LENGTH_IMGT"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
#   select(SEQUENCE_ID)

filteredData$CDR3LENGTH_IMGT <- as.numeric(filteredData$CDR3LENGTH_IMGT)
filteredData$G_J_CDR3 <- paste(filteredData$GF_JGENE,filteredData$CDR3LENGTH_IMGT,sep="_")
gjcdr3.title <- filteredData %>%
  select(G_J_CDR3)
gjcdr3.titleID <- gjcdr3.title$G_J_CDR3[1]
## then use gjcdr3.titleID to replace current title in " "

#################################
### END COMMANDS FOR UPDATING COVID ABDAB DATABASE BY JUST GRABBING RAW DATA...

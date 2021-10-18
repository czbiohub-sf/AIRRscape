## THIS IS A LARGE FILE WITH MANY R COMMANDS - PRIMARILY INTENDED TO PRE-PROCESS MAB DATASETS TO THEN RUN IN THE SHINY TOOL...
## BEFORE GOING PUBLIC WOULD BE GOOD TO CLEAN THIS UP...

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

# pkgs = c("alakazam", "igraph", "dplyr","RColorBrewer", "hexbin", "scales","grid", "lattice", "gdata","gridExtra", "ape", "shazam","limma", "reshape2", "DT","ggplot2", "seqinr", "phangorn","shiny", "tidyverse") # package names
pkgs = c("alakazam", "igraph", "dplyr","RColorBrewer", "hexbin", "scales","grid", "lattice", "gdata","gridExtra", "ape", "shazam","reshape2", "DT","ggplot2", "seqinr", "phangorn","shiny", "tidyverse") # package names
install.packages(pkgs)
  inst = lapply(pkgs, library, character.only = TRUE) # load them

  ## Warning in install.packages :
  # package ‘grid’ is a base package, and should not be updated
  # Warning in install.packages :
  #   package ‘limma’ is not available for this version of R
  
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
########
## end of 2 functions




##################################################################################################################
### COMMANDS FOR READING IN AIRR-COMPLIANT DATASETS FOR SHINY APP
### note that for non-AIRR (e.g. older Immcantation data or databases like AbDab) datasets, column names will be different
###### particularly note that 1) JUNCTION MAY NOT NEED TO BE TRIMMED 1 AA ON EACH END TO BE KABAT CDR3
######                        2) V_IDENTITY MAY BE BETWEEN 0-1 NOT 0-100, AND THUS SHM FORMULA NEEDS TO CHANGE
##################################################################################################################

### for heavy + light datasets, better to seprately load hc & lc data, then combine later if necessary
## reason - see below about removing ends of airr-defined junction but only in heavy datasets!!!

BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/hiv_iavi84and152hc.tab")

## this removes first and last aa from junction (because ireceptor includes these) BE CAREFUL NOT TO EVER RUN THIS ON LC DATA...WRONG  SEE BELOW
# str_sub(BX.clusters.hiv$JUNCTIONAA, 1, 1) <- ""
# str_sub(BX.clusters.hiv$JUNCTIONAA, -1, -1) <- ""

## these are headers in latest Immcantation germ-pass.tsv files NOTE THEY ARE NOW TSV!!
## sequence_id	sequence	rev_comp	productive	v_call	d_call	j_call	sequence_alignment	germline_alignment	junction	
## junction_aa	v_cigar	d_cigar	j_cigar	stop_codon	vj_in_frame	locus	junction_length	np1_length	np2_length
## v_sequence_start	v_sequence_end	v_germline_start	v_germline_end	d_sequence_start	d_sequence_end	d_germline_start	d_germline_end	j_sequence_start	j_sequence_end
## j_germline_start	j_germline_end	v_score	v_identity	v_support	d_score	d_identity	d_support	j_score	j_identity
## j_support	fwr1	fwr2	fwr3	fwr4	cdr1	cdr2	cdr3	consensus_count	duplicate_count
## clone_id	prcons	seqorient	cregion	ganda_subtype	germline_alignment_d_mask	germline_v_call	germline_d_call	germline_j_call

BX.clusters.test <- read.delim("bcells_1M_germ-pass.tsv")

## NOTE THAT THIS ALREADY HAS junction_aa ALSO THIS COLUMN IS SAME IN BOTH HEAVY AND LIGHT, JUST NEEDS FIRST AND LAST RESIDUES REMOVED

### AIRR-COMPLIANT TSV FILES NEED THE FOLLOWING MODIFICATIONS
# CHANGING RELEVANT HEADERS FROM LOWERCASE TO UPPER CASE
# REMOVING FIRST AND LAST RESIDUES FROM JUNCTION
# CONVERT V_IDENTITY TO SHM - FOR SHINY WANT 0-100


## EXAMPLE TO MODIFY...
BX.clusters.test <- read.delim("bcells_1M_germ-pass.tsv")

## BUT ALSO WILL WANT TO CHANGE THE NAME OF EACH SEQUENCE!!! DONE HERE...
BX.clusters.test$sequence_id <- paste0("test",1:nrow(BX.clusters.test))

str_sub(BX.clusters.test$junction_aa, 1, 1) <- ""
str_sub(BX.clusters.test$junction_aa, -1, -1) <- ""

BX.clusters.test$JUNCTIONAA <- as.character(BX.clusters.test$junction_aa)
BX.clusters.test$CDR3LENGTH_IMGT <- nchar(BX.clusters.test$junction_aa)

BX.clusters.test$GENE <- getGene(BX.clusters.test$germline_v_call, first=TRUE, strip_d=TRUE)
BX.clusters.test$GF <- substring(BX.clusters.test$GENE, 1,5)
BX.clusters.test$JGENE <- getGene(BX.clusters.test$germline_j_call, first=TRUE, strip_d=TRUE)
BX.clusters.test$JGF <- substring(BX.clusters.test$JGENE, 1,5)

BX.clusters.test <- BX.clusters.test %>% add_count(clone_id) %>%
  rename(reads_per_clone = n)

### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
#BX.clusters.mabs.toadd$SHM <- (100 - (BX.clusters.mabs.toadd$V_IDENTITY))
BX.clusters.test$SHM <- (100 - (BX.clusters.test$v_identity * 100))
BX.clusters.test <- BX.clusters.test %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

## MAY NOT ALWAYS WANT TO DO THIS - also note IgD in here (could also remove??)
BX.clusters.test$cregion <- gsub("IgA","IgH",BX.clusters.test$cregion)
BX.clusters.test$cregion <- gsub("IgG","IgH",BX.clusters.test$cregion)
BX.clusters.test$cregion <- gsub("IgM","IgH",BX.clusters.test$cregion)
BX.clusters.test$cregion <- gsub("IgD","IgH",BX.clusters.test$cregion)

BX.clusters.test <- BX.clusters.test %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2))) %>%
  rename(SEQUENCE_ID = sequence_id) %>%
  rename(CREGION = cregion)

BX.clusters.test <- BX.clusters.test %>% filter(CDR3LENGTH_IMGT > 3.8)  
##


toshiny.test <- BX.clusters.test %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

rm(BX.clusters.test)

toshiny.test.h <- subset(toshiny.test, CREGION %in% c("IgH"))

#write.table(toshiny.test, "toshiny_test.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.test.h, "toshiny_test_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# ggplot(toshiny.test.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.test.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## columns used in app:
# SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT, CREGION, shm.mean, SHM


## NEXT GENERIC CODE TO COMBINE 2-6 DATASETS INTO A SINGLE FILE FOR SHINY VISUALIZATION (BOTH SEPARATE & COMBINED)

toshiny.test.many <- bind_rows(toshiny.test.h, toshiny.test2.h, toshiny.test3.h, toshiny.test4.h, .id = "id")

toshiny.test.many$id <- gsub("1","test1",toshiny.fewhcbg$id)
toshiny.test.many$id <- gsub("2","test2",toshiny.fewhcbg$id)
toshiny.test.many$id <- gsub("3","test3",toshiny.fewhcbg$id)
toshiny.test.many$id <- gsub("4","test4",toshiny.fewhcbg$id)

## for combined need to recalculate ncount and shm.mean
toshiny.test.manyc <- toshiny.test.many
toshiny.test.manyc$ncount <- NULL
toshiny.test.manyc$shm.mean <- NULL

toshiny.test.manyc <- toshiny.fewhcbgc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.test.many, "toshiny_testmany.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.test.manyc, "toshiny_testmanyc.tab", sep = "\t", row.names = FALSE, quote = FALSE)


##################################################################################################################
### END COMMANDS FOR READING IN AIRR-COMPLIANT DATASETS FOR SHINY APP


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
### 3/3/2021 reloading comet data & combining with HC & covid datasets...
## note comet data has no shm...
toshiny.few2$SHM <- NA
toshiny.few2$JGENE <- toshiny.few2$JGF
toshiny.few2$JUNCTION_LENGTH <- (3 * (toshiny.few2$CDR3LENGTH_IMGT)) + 6
toshiny.few2$VJ_JUNCTION_PATTERN <- paste(toshiny.few2$GF_JGENE,toshiny.few2$JUNCTION_LENGTH,sep="_") 
toshiny.few2$JUNCTION_LENGTH <- NULL

comet.pos <- subset(toshiny.few2, id %in% c("COMET COVID-positive"))
comet.neg <- subset(toshiny.few2, id %in% c("COMET COVID-negative"))
comet.pos$id <- NULL
comet.neg$id <- NULL
#comet.pos <- filter(BX.clusters.hc3, JGENE == "IGLJ1" & GF == "IGHV3")
# toshiny.few3 <- read.delim("toshiny_few3.tab")
# toshiny.few3.h <- subset(toshiny.few3, CREGION %in% c("IgH"))

BX.clusters.hcvsboydvsgalsonmabscomet <- bind_rows(toshiny.few3.h, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, comet.pos, comet.neg, .id = "id")
#toshiny.fewhcbg <- subset(BX.clusters.hcvsboydvsgalsonmabs, CREGION %in% c("IgH"))
toshiny.fewhcbgcomet <- BX.clusters.hcvsboydvsgalsonmabscomet[ grep("Kappa|Lambda", BX.clusters.hcvsboydvsgalsonmabscomet$CREGION, invert = TRUE) , ] %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.fewhcbgcomet$id <- gsub("1","Healthy control",toshiny.fewhcbgcomet$id)
## be careful with 1 in name here....
toshiny.fewhcbgcomet$id <- gsub("4","anti-CoVII mAbs",toshiny.fewhcbgcomet$id)
toshiny.fewhcbgcomet$id <- gsub("5","COMET COVID-positive",toshiny.fewhcbgcomet$id)
toshiny.fewhcbgcomet$id <- gsub("6","COMET COVID-negative",toshiny.fewhcbgcomet$id)


toshiny.fewhcbgcomet$id <- gsub("2","Boyd7450",toshiny.fewhcbgcomet$id)
toshiny.fewhcbgcomet$id <- gsub("3","Galson1",toshiny.fewhcbgcomet$id)
toshiny.fewhcbgcomet$id <- gsub("anti-CoVII mAbs","anti-CoV2 mAbs",toshiny.fewhcbgcomet$id)
#write.table(toshiny.fewhcbgcomet, "toshiny_fewhcbgcomet.tab", sep = "\t", row.names = FALSE, quote = FALSE)

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


### for updating HIV mabs with Filip's newly found mab:

toshiny.fewhivmabsnih45rd1214 <- bind_rows(toshiny.fewhivmabs.h, toshiny.fewhivnih45, toshiny.fewhiv, .id = "id")

toshiny.fewhivmabsnih45rd1214$id <- gsub("1","anti-HIV mAbs",toshiny.fewhivmabsnih45rd1214$id)
toshiny.fewhivmabsnih45rd1214$id <- gsub("2","HIV+ patient NIH45",toshiny.fewhivmabsnih45rd1214$id)
toshiny.fewhivmabsnih45rd1214$id <- gsub("3","HIV+ patient AD1214",toshiny.fewhivmabsnih45rd1214$id)

toshiny.fewhivmabsnih45rd1214$id <- factor(toshiny.fewhivmabsnih45rd1214$id, levels = c("HIV+ patient AD1214", "HIV+ patient NIH45", "anti-HIV mAbs"))

## for combined need to recalculate ncount and shm.mean
toshiny.fewhivmabsnih45rd1214c <- toshiny.fewhivmabsnih45rd1214
toshiny.fewhivmabsnih45rd1214c$ncount <- NULL
toshiny.fewhivmabsnih45rd1214c$shm.mean <- NULL

toshiny.fewhivmabsnih45rd1214c <- toshiny.fewhivmabsnih45rd1214c %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.fewhivmabsnih45rd1214, "toshiny_mabsnih45rd1214.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhivmabsnih45rd1214c, "toshiny_mabsnih45rd1214c.tab", sep = "\t", row.names = FALSE, quote = FALSE)



##################################################################################################################
##################################################################################################################
##################################################################################################################

########
## reading in previous healthy germ-pass samples to get clonalclusters
healthy <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/BXmay10m_stim_germ-pass.tab")
BX.full0 <- healthy

healthy2 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/20181125_d147km_10m_stimulated_germ-pass.tab")
BX.full0 <- healthy2

healthy3 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/20181029r_d682_10m_stimulated_germ-pass.tab")
BX.full0 <- healthy3


BX.full0$V_CALL <- as.character(BX.full0$V_CALL)
BX.full0$J_CALL <- as.character(BX.full0$J_CALL)
BX.full0$V_CALL <- as.character(BX.full0$V_CALL)

#BX.full0 <- BX.full0 %>% add_column(CONSCOUNT = 1, DUPCOUNT = 1)
BX.full0 <- BX.full0 %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

BX.clones.hc <- collapseClones(BX.full0, method="mostCommon", cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT", germlineColumn = "GERMLINE_IMGT_D_MASK",
                               includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#BX.clusters.hc <- BX.clusters2.hc
BX.clusters.hc3 <- groupGenes(BX.clones.hc, v_call = "V_CALL", j_call = "J_CALL",
                              junc_len = "JUNCTION_LENGTH", cell_id = NULL, locus = NULL, only_igh = TRUE,
                              first = TRUE)
rm(BX.clones.hc)
BX.clusters.hc3$JUNCTIONAA <- translateDNA(BX.clusters.hc3$JUNCTION, trim=TRUE)
BX.clusters.hc3$GENE <- getGene(BX.clusters.hc3$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hc3$GF <- substring(BX.clusters.hc3$GENE, 1,5)
BX.clusters.hc3$JGENE <- getGene(BX.clusters.hc3$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hc3$JGF <- substring(BX.clusters.hc3$JGENE, 1,5)
BX.clusters.hc3$VJ_JUNCTION_PATTERN <- paste(BX.clusters.hc3$GF,BX.clusters.hc3$JGF,BX.clusters.hc3$JUNCTION_LENGTH,sep="_") 

write.table(BX.clusters.hc, "Galson_healthycontrol_clonalclusters.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(BX.clusters.hc2, "Galson_healthycontrol2_clonalclusters.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(BX.clusters.hc3, "Galson_healthycontrol3_clonalclusters.tab", sep = "\t", row.names = FALSE, quote = FALSE)



rm(BX.full0)
rm(healthy3)


### late Nov 2020 re-loading hc2 & hc3, to combine and examine with anti-Flu mabs (MAKE THIS FROM IEDB SET!)

BX.clusters.hc2 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/Galson_healthycontrol2_clonalclusters.tab")
BX.clusters.hc3 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/Galson_healthycontrol3_clonalclusters.tab")


BX.clusters.hc2$JUNCTIONAA <- as.character(BX.clusters.hc2$JUNCTIONAA)
BX.clusters.hc2$CDR3LENGTH_IMGT <- nchar(BX.clusters.hc2$JUNCTIONAA)
### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
#BX.clusters.mabs.toadd$SHM <- (100 - (BX.clusters.mabs.toadd$V_IDENTITY))
BX.clusters.hc2$SHM <- (100 - (BX.clusters.hc2$V_IDENTITY * 100))
BX.clusters.hc2 <- BX.clusters.hc2 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

BX.clusters.hc2$CREGION <- gsub("IgA","IgH",BX.clusters.hc2$CREGION)
BX.clusters.hc2$CREGION <- gsub("IgG","IgH",BX.clusters.hc2$CREGION)
BX.clusters.hc2$CREGION <- gsub("IgM","IgH",BX.clusters.hc2$CREGION)

BX.clusters.hc2 <- BX.clusters.hc2 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))
##
BX.clusters.hc3$JUNCTIONAA <- as.character(BX.clusters.hc3$JUNCTIONAA)
BX.clusters.hc3$CDR3LENGTH_IMGT <- nchar(BX.clusters.hc3$JUNCTIONAA)
### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
#BX.clusters.mabs.toadd$SHM <- (100 - (BX.clusters.mabs.toadd$V_IDENTITY))
BX.clusters.hc3$SHM <- (100 - (BX.clusters.hc3$V_IDENTITY * 100))
BX.clusters.hc3 <- BX.clusters.hc3 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

BX.clusters.hc3$CREGION <- gsub("IgA","IgH",BX.clusters.hc3$CREGION)
BX.clusters.hc3$CREGION <- gsub("IgG","IgH",BX.clusters.hc3$CREGION)
BX.clusters.hc3$CREGION <- gsub("IgM","IgH",BX.clusters.hc3$CREGION)

BX.clusters.hc3 <- BX.clusters.hc3 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))


## NOTE HC1 IS CALLED toshiny.few3
toshiny.hc2 <- BX.clusters.hc2 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.hc3 <- BX.clusters.hc3 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

rm(BX.clusters.hc2)
rm(BX.clusters.hc3)

toshiny.hc2.h <- subset(toshiny.hc2, CREGION %in% c("IgH"))
toshiny.hc3.h <- subset(toshiny.hc3, CREGION %in% c("IgH"))

rm(toshiny.hc2)
rm(toshiny.hc3)

## save, then rename id's, then reload

#write.table(toshiny.hc2.h, "toshiny_hc2h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.hc3.h, "toshiny_hc3h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
toshiny.hc2.h <- read.delim("toshiny_hc2h.tab")
toshiny.hc3.h <- read.delim("toshiny_hc3h.tab")

### flu mabs - ADDING IN A SIMILAR FASHION TO THE ITERATIVELY ADDED ABDAB DATASET...
toshiny.flumabs.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/flu_mablist.tab")
toshiny.flumabs.l <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/flu_mablist.tab")

toshiny.flumabs.h$SEQUENCE_IDLC <- NULL
toshiny.flumabs.h$GENELC <- NULL
toshiny.flumabs.h$JGENELC <- NULL

toshiny.flumabs.h$JUNCTIONAALC <- NULL
toshiny.flumabs.h$V_IDENTITYLC <- NULL

toshiny.flumabs.l$SEQUENCE_ID <- NULL
toshiny.flumabs.l <- rename(toshiny.flumabs.l, SEQUENCE_ID = SEQUENCE_IDLC)

toshiny.flumabs.l$GENE <- NULL
toshiny.flumabs.l <- rename(toshiny.flumabs.l, GENE = GENELC)
toshiny.flumabs.l$JGENE <- NULL
toshiny.flumabs.l <- rename(toshiny.flumabs.l, JGENE = JGENELC)

toshiny.flumabs.l$JUNCTIONAA <- NULL
toshiny.flumabs.l <- rename(toshiny.flumabs.l, JUNCTIONAA = JUNCTIONAALC)
toshiny.flumabs.l$V_IDENTITY <- NULL
toshiny.flumabs.l <- rename(toshiny.flumabs.l, V_IDENTITY = V_IDENTITYLC)


toshiny.flumabs.h <- rename(toshiny.flumabs.h, V_CALL = GENE)
toshiny.flumabs.h <- rename(toshiny.flumabs.h, J_CALL = JGENE)
toshiny.flumabs.l <- rename(toshiny.flumabs.l, V_CALL = GENE)
toshiny.flumabs.l <- rename(toshiny.flumabs.l, J_CALL = JGENE)

toshiny.flumabs.h$V_CALL <- as.character(toshiny.flumabs.h$V_CALL)
toshiny.flumabs.h$J_CALL <- as.character(toshiny.flumabs.h$J_CALL)
toshiny.flumabs.h$GENE <- getGene(toshiny.flumabs.h$V_CALL, first=TRUE, strip_d=TRUE)
toshiny.flumabs.h$JGENE <- getGene(toshiny.flumabs.h$J_CALL, first=TRUE, strip_d=TRUE)
toshiny.flumabs.l$V_CALL <- as.character(toshiny.flumabs.l$V_CALL)
toshiny.flumabs.l$J_CALL <- as.character(toshiny.flumabs.l$J_CALL)
toshiny.flumabs.l$GENE <- getGene(toshiny.flumabs.l$V_CALL, first=TRUE, strip_d=TRUE)
toshiny.flumabs.l$JGENE <- getGene(toshiny.flumabs.l$J_CALL, first=TRUE, strip_d=TRUE)


toshiny.flumabs.h$Chain.1.Nucleotide <- NULL
toshiny.flumabs.h$Chain.2.Nucleotide <- NULL
toshiny.flumabs.h$Chain.1.Full.Sequence <- NULL
toshiny.flumabs.h$Chain.2.Full.Sequence <- NULL

toshiny.flumabs.l$Chain.1.Nucleotide <- NULL
toshiny.flumabs.l$Chain.2.Nucleotide <- NULL
toshiny.flumabs.l$Chain.1.Full.Sequence <- NULL
toshiny.flumabs.l$Chain.2.Full.Sequence <- NULL

toshiny.flumabs <- rbind(toshiny.flumabs.h, toshiny.flumabs.l)

toshiny.flumabs$JUNCTIONAA <- as.character(toshiny.flumabs$JUNCTIONAA)
toshiny.flumabs$SEQUENCE_ID <- as.character(toshiny.flumabs$SEQUENCE_ID)

toshiny.flumabs$CDR3LENGTH_IMGT <- nchar(toshiny.flumabs$JUNCTIONAA)

toshiny.flumabs$JUNCTION_LENGTH <- (3 * (toshiny.flumabs$CDR3LENGTH_IMGT)) + 6
toshiny.flumabs$JGF <- toshiny.flumabs$JGENE
toshiny.flumabs$GF <- str_sub(toshiny.flumabs$GENE, end=5)
### note updated str_sub for substring (think both work)
# toshiny.flumabs$GF <- substring(toshiny.flumabs$GENE, 1,5)
# toshiny.flumabs$JGF <- substring(toshiny.flumabs$JGENE, 1,5)
toshiny.flumabs$VJ_JUNCTION_PATTERN <- paste(toshiny.flumabs$GF,toshiny.flumabs$JGF,toshiny.flumabs$JUNCTION_LENGTH,sep="_") 

### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
# toshiny.flumabs$SHM <- toshiny.flumabs$V_IDENTITY * 100
# toshiny.flumabs$SHM <- (100 - (toshiny.flumabs$V_IDENTITY))
toshiny.flumabs$SHM <- (100 - (toshiny.flumabs$V_IDENTITY * 100))
toshiny.flumabs <- toshiny.flumabs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

toshiny.flumabs$JUNCTION_LENGTH <- NULL
toshiny.flumabs$V_IDENTITY <- NULL

toshiny.flumabs$V_CALL <- NULL
toshiny.flumabs$J_CALL <- NULL

## ADDING CREGION BY WAY OF GF...
toshiny.flumabs$CREGION <- str_sub(toshiny.flumabs$GF, end=3)
toshiny.flumabs$CREGION <- gsub("IGK","Kappa",toshiny.flumabs$CREGION)
toshiny.flumabs$CREGION <- gsub("IGL","Lambda",toshiny.flumabs$CREGION)
toshiny.flumabs$CREGION <- gsub("IGH","IgH",toshiny.flumabs$CREGION)

toshiny.flumabs <- toshiny.flumabs %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

rm(toshiny.flumabs.h)
rm(toshiny.flumabs.l)
toshiny.flumabs.h <- subset(toshiny.flumabs, CREGION %in% c("IgH"))

#write.table(toshiny.flumabs, "toshiny_flumabs.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.flumabs.h, "toshiny_flumabs_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## then combine hc hc2 hc3 (also flu mabs?)

toshiny.fluhcmabs <- bind_rows(toshiny.flumabs.h, toshiny.few3.h, toshiny.hc2.h, toshiny.hc3.h, .id = "id")

toshiny.fluhcmabs$id <- gsub("1","anti-flu mAbs",toshiny.fluhcmabs$id)
## be careful with 1 in name here....
toshiny.fluhcmabs$id <- gsub("4","Healthy vaccinated control III",toshiny.fluhcmabs$id)
toshiny.fluhcmabs$id <- gsub("2","Healthy vaccinated control",toshiny.fluhcmabs$id)
toshiny.fluhcmabs$id <- gsub("3","Healthy vaccinated control II",toshiny.fluhcmabs$id)

## for combined need to recalculate ncount and shm.mean
toshiny.fluhcmabsc <- toshiny.fluhcmabs
toshiny.fluhcmabsc$ncount <- NULL
toshiny.fluhcmabsc$shm.mean <- NULL

toshiny.fluhcmabsc <- toshiny.fluhcmabsc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.fluhcmabs, "toshiny_fluhcmabs.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fluhcmabsc, "toshiny_fluhcmabsc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

######################################
################# dec 2020 adding some paired OAS data (2 HIV donors, one is NIH45 & 6 healthy donors)

## nih45
toshiny.pairedhivnih45.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/paired_OAS/libraseq_nih45list.tab")
toshiny.pairedhivnih45.l <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/paired_OAS/libraseq_nih45list.tab")

toshiny.pairedhivnih45.h$SEQUENCE_IDLC <- NULL
toshiny.pairedhivnih45.h$GENELC <- NULL
toshiny.pairedhivnih45.h$JGENELC <- NULL

toshiny.pairedhivnih45.h$JUNCTIONAALC <- NULL
toshiny.pairedhivnih45.h$V_IDENTITYLC <- NULL

toshiny.pairedhivnih45.l$SEQUENCE_ID <- NULL
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, SEQUENCE_ID = SEQUENCE_IDLC)

toshiny.pairedhivnih45.l$GENE <- NULL
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, GENE = GENELC)
toshiny.pairedhivnih45.l$JGENE <- NULL
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, JGENE = JGENELC)

toshiny.pairedhivnih45.l$JUNCTIONAA <- NULL
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, JUNCTIONAA = JUNCTIONAALC)
toshiny.pairedhivnih45.l$V_IDENTITY <- NULL
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, V_IDENTITY = V_IDENTITYLC)

## changed for OAS
toshiny.pairedhivnih45.h <- rename(toshiny.pairedhivnih45.h, V_CALL = GENE)
toshiny.pairedhivnih45.h <- rename(toshiny.pairedhivnih45.h, J_CALL = JGENE)
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, V_CALL = GENE)
toshiny.pairedhivnih45.l <- rename(toshiny.pairedhivnih45.l, J_CALL = JGENE)


toshiny.pairedhivnih45.h$Chain.1.Nucleotide <- NULL
toshiny.pairedhivnih45.h$Chain.2.Nucleotide <- NULL

toshiny.pairedhivnih45.l$Chain.1.Nucleotide <- NULL
toshiny.pairedhivnih45.l$Chain.2.Nucleotide <- NULL
toshiny.pairedhivnih45.l$SUBTYPE <- NULL

toshiny.pairedhivnih45 <- full_join(toshiny.pairedhivnih45.h, toshiny.pairedhivnih45.l)

toshiny.pairedhivnih45$JUNCTIONAA <- as.character(toshiny.pairedhivnih45$JUNCTIONAA)
toshiny.pairedhivnih45$SEQUENCE_ID <- as.character(toshiny.pairedhivnih45$SEQUENCE_ID)

### FOR OAS NEED TO REMOVE FIRST AND LAST AA FROM ALL JUNCTIONS (BOTH HEAVY & LIGHT)
str_sub(toshiny.pairedhivnih45$JUNCTIONAA, 1, 1) <- ""
str_sub(toshiny.pairedhivnih45$JUNCTIONAA, -1, -1) <- ""


toshiny.pairedhivnih45$CDR3LENGTH_IMGT <- nchar(toshiny.pairedhivnih45$JUNCTIONAA)

## CHECK THIS FOR THESE OAS PAIRED DATASETS!!!
toshiny.pairedhivnih45$JUNCTION_LENGTH <- (3 * (toshiny.pairedhivnih45$CDR3LENGTH_IMGT)) + 6

toshiny.pairedhivnih45$GENE <- getGene(toshiny.pairedhivnih45$V_CALL, first=TRUE, strip_d=TRUE)
toshiny.pairedhivnih45$GF <- substring(toshiny.pairedhivnih45$GENE, 1,5)
toshiny.pairedhivnih45$JGENE <- getGene(toshiny.pairedhivnih45$J_CALL, first=TRUE, strip_d=TRUE)
toshiny.pairedhivnih45$JGF <- substring(toshiny.pairedhivnih45$JGENE, 1,5)

# toshiny.pairedhivnih45$JGF <- toshiny.pairedhivnih45$JGENE
# toshiny.pairedhivnih45$GF <- str_sub(toshiny.pairedhivnih45$GENE, end=5)
toshiny.pairedhivnih45$VJ_JUNCTION_PATTERN <- paste(toshiny.pairedhivnih45$GF,toshiny.pairedhivnih45$JGF,toshiny.pairedhivnih45$JUNCTION_LENGTH,sep="_") 

### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
# toshiny.pairedhivnih45$SHM <- toshiny.pairedhivnih45$V_IDENTITY * 100
toshiny.pairedhivnih45$SHM <- (100 - (toshiny.pairedhivnih45$V_IDENTITY))
# toshiny.pairedhivnih45$SHM <- (100 - (toshiny.pairedhivnih45$V_IDENTITY * 100))
toshiny.pairedhivnih45 <- toshiny.pairedhivnih45 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

toshiny.pairedhivnih45$JUNCTION_LENGTH <- NULL
toshiny.pairedhivnih45$V_IDENTITY <- NULL

toshiny.pairedhivnih45$V_CALL <- NULL
toshiny.pairedhivnih45$J_CALL <- NULL

### NOTE OAS PAIRED HAS ISOTYPE & SUBTYPE!!!!!
## ADDING CREGION BY WAY OF GF...
toshiny.pairedhivnih45$CREGION <- str_sub(toshiny.pairedhivnih45$GF, end=3)
toshiny.pairedhivnih45$CREGION <- gsub("IGK","Kappa",toshiny.pairedhivnih45$CREGION)
toshiny.pairedhivnih45$CREGION <- gsub("IGL","Lambda",toshiny.pairedhivnih45$CREGION)
toshiny.pairedhivnih45$CREGION <- gsub("IGH","IgH",toshiny.pairedhivnih45$CREGION)

toshiny.pairedhivnih45 <- toshiny.pairedhivnih45 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

rm(toshiny.pairedhivnih45.h)
rm(toshiny.pairedhivnih45.l)
toshiny.pairedhivnih45.h <- subset(toshiny.pairedhivnih45, CREGION %in% c("IgH"))

#write.table(toshiny.pairedhivnih45, "toshiny_pairedhivnih45.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.pairedhivnih45.h, "toshiny_pairedhivnih45_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# ggplot(toshiny.pairedhivnih45.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ SUBTYPE, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.pairedhivnih45.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.pairedhivnih45.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## n90
toshiny.pairedhivn90.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/paired_OAS/libraseq_n90list.tab")
toshiny.pairedhivn90.l <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/paired_OAS/libraseq_n90list.tab")

toshiny.pairedhivn90.h$SEQUENCE_IDLC <- NULL
toshiny.pairedhivn90.h$GENELC <- NULL
toshiny.pairedhivn90.h$JGENELC <- NULL

toshiny.pairedhivn90.h$JUNCTIONAALC <- NULL
toshiny.pairedhivn90.h$V_IDENTITYLC <- NULL

toshiny.pairedhivn90.l$SEQUENCE_ID <- NULL
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, SEQUENCE_ID = SEQUENCE_IDLC)

toshiny.pairedhivn90.l$GENE <- NULL
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, GENE = GENELC)
toshiny.pairedhivn90.l$JGENE <- NULL
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, JGENE = JGENELC)

toshiny.pairedhivn90.l$JUNCTIONAA <- NULL
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, JUNCTIONAA = JUNCTIONAALC)
toshiny.pairedhivn90.l$V_IDENTITY <- NULL
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, V_IDENTITY = V_IDENTITYLC)

## changed for OAS
toshiny.pairedhivn90.h <- rename(toshiny.pairedhivn90.h, V_CALL = GENE)
toshiny.pairedhivn90.h <- rename(toshiny.pairedhivn90.h, J_CALL = JGENE)
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, V_CALL = GENE)
toshiny.pairedhivn90.l <- rename(toshiny.pairedhivn90.l, J_CALL = JGENE)

toshiny.pairedhivn90.h$Chain.1.Nucleotide <- NULL
toshiny.pairedhivn90.h$Chain.2.Nucleotide <- NULL

toshiny.pairedhivn90.l$Chain.1.Nucleotide <- NULL
toshiny.pairedhivn90.l$Chain.2.Nucleotide <- NULL
toshiny.pairedhivn90.l$SUBTYPE <- NULL

toshiny.pairedhivn90 <- full_join(toshiny.pairedhivn90.h, toshiny.pairedhivn90.l)

toshiny.pairedhivn90$JUNCTIONAA <- as.character(toshiny.pairedhivn90$JUNCTIONAA)
toshiny.pairedhivn90$SEQUENCE_ID <- as.character(toshiny.pairedhivn90$SEQUENCE_ID)

### FOR OAS NEED TO REMOVE FIRST AND LAST AA FROM ALL JUNCTIONS (BOTH HEAVY & LIGHT)
str_sub(toshiny.pairedhivn90$JUNCTIONAA, 1, 1) <- ""
str_sub(toshiny.pairedhivn90$JUNCTIONAA, -1, -1) <- ""


toshiny.pairedhivn90$CDR3LENGTH_IMGT <- nchar(toshiny.pairedhivn90$JUNCTIONAA)

## CHECK THIS FOR THESE OAS PAIRED DATASETS!!!
toshiny.pairedhivn90$JUNCTION_LENGTH <- (3 * (toshiny.pairedhivn90$CDR3LENGTH_IMGT)) + 6

toshiny.pairedhivn90$GENE <- getGene(toshiny.pairedhivn90$V_CALL, first=TRUE, strip_d=TRUE)
toshiny.pairedhivn90$GF <- substring(toshiny.pairedhivn90$GENE, 1,5)
toshiny.pairedhivn90$JGENE <- getGene(toshiny.pairedhivn90$J_CALL, first=TRUE, strip_d=TRUE)
toshiny.pairedhivn90$JGF <- substring(toshiny.pairedhivn90$JGENE, 1,5)

# toshiny.pairedhivn90$JGF <- toshiny.pairedhivn90$JGENE
# toshiny.pairedhivn90$GF <- str_sub(toshiny.pairedhivn90$GENE, end=5)
toshiny.pairedhivn90$VJ_JUNCTION_PATTERN <- paste(toshiny.pairedhivn90$GF,toshiny.pairedhivn90$JGF,toshiny.pairedhivn90$JUNCTION_LENGTH,sep="_") 

### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
# toshiny.pairedhivn90$SHM <- toshiny.pairedhivn90$V_IDENTITY * 100
toshiny.pairedhivn90$SHM <- (100 - (toshiny.pairedhivn90$V_IDENTITY))
# toshiny.pairedhivn90$SHM <- (100 - (toshiny.pairedhivn90$V_IDENTITY * 100))
toshiny.pairedhivn90 <- toshiny.pairedhivn90 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

toshiny.pairedhivn90$JUNCTION_LENGTH <- NULL
toshiny.pairedhivn90$V_IDENTITY <- NULL

toshiny.pairedhivn90$V_CALL <- NULL
toshiny.pairedhivn90$J_CALL <- NULL

### NOTE OAS PAIRED HAS ISOTYPE & SUBTYPE!!!!!
## ADDING CREGION BY WAY OF GF...
toshiny.pairedhivn90$CREGION <- str_sub(toshiny.pairedhivn90$GF, end=3)
toshiny.pairedhivn90$CREGION <- gsub("IGK","Kappa",toshiny.pairedhivn90$CREGION)
toshiny.pairedhivn90$CREGION <- gsub("IGL","Lambda",toshiny.pairedhivn90$CREGION)
toshiny.pairedhivn90$CREGION <- gsub("IGH","IgH",toshiny.pairedhivn90$CREGION)

toshiny.pairedhivn90 <- toshiny.pairedhivn90 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

rm(toshiny.pairedhivn90.h)
rm(toshiny.pairedhivn90.l)
toshiny.pairedhivn90.h <- subset(toshiny.pairedhivn90, CREGION %in% c("IgH"))

#write.table(toshiny.pairedhivn90, "toshiny_pairedhivn90.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.pairedhivn90.h, "toshiny_pairedhivn90_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# ggplot(toshiny.pairedhivn90.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ SUBTYPE, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.pairedhivn90.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.pairedhivn90.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


#### 6 healthy juveniles

toshiny.paired6juvs.h <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/paired_OAS/pairedseq_6juveniles.tab")
toshiny.paired6juvs.l <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/paired_OAS/pairedseq_6juveniles.tab")

toshiny.paired6juvs.h$SEQUENCE_IDLC <- NULL
toshiny.paired6juvs.h$GENELC <- NULL
toshiny.paired6juvs.h$JGENELC <- NULL

toshiny.paired6juvs.h$JUNCTIONAALC <- NULL
toshiny.paired6juvs.h$V_IDENTITYLC <- NULL

toshiny.paired6juvs.l$SEQUENCE_ID <- NULL
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, SEQUENCE_ID = SEQUENCE_IDLC)

toshiny.paired6juvs.l$GENE <- NULL
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, GENE = GENELC)
toshiny.paired6juvs.l$JGENE <- NULL
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, JGENE = JGENELC)

toshiny.paired6juvs.l$JUNCTIONAA <- NULL
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, JUNCTIONAA = JUNCTIONAALC)
toshiny.paired6juvs.l$V_IDENTITY <- NULL
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, V_IDENTITY = V_IDENTITYLC)

## changed for OAS
toshiny.paired6juvs.h <- rename(toshiny.paired6juvs.h, V_CALL = GENE)
toshiny.paired6juvs.h <- rename(toshiny.paired6juvs.h, J_CALL = JGENE)
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, V_CALL = GENE)
toshiny.paired6juvs.l <- rename(toshiny.paired6juvs.l, J_CALL = JGENE)

toshiny.paired6juvs.h$Chain.1.Nucleotide <- NULL
toshiny.paired6juvs.h$Chain.2.Nucleotide <- NULL

toshiny.paired6juvs.l$Chain.1.Nucleotide <- NULL
toshiny.paired6juvs.l$Chain.2.Nucleotide <- NULL
toshiny.paired6juvs.l$SUBTYPE <- NULL

toshiny.paired6juvs <- full_join(toshiny.paired6juvs.h, toshiny.paired6juvs.l)

toshiny.paired6juvs$JUNCTIONAA <- as.character(toshiny.paired6juvs$JUNCTIONAA)
toshiny.paired6juvs$SEQUENCE_ID <- as.character(toshiny.paired6juvs$SEQUENCE_ID)

### FOR OAS NEED TO REMOVE FIRST AND LAST AA FROM ALL JUNCTIONS (BOTH HEAVY & LIGHT)
str_sub(toshiny.paired6juvs$JUNCTIONAA, 1, 1) <- ""
str_sub(toshiny.paired6juvs$JUNCTIONAA, -1, -1) <- ""


toshiny.paired6juvs$CDR3LENGTH_IMGT <- nchar(toshiny.paired6juvs$JUNCTIONAA)

## CHECK THIS FOR THESE OAS PAIRED DATASETS!!!
toshiny.paired6juvs$JUNCTION_LENGTH <- (3 * (toshiny.paired6juvs$CDR3LENGTH_IMGT)) + 6

toshiny.paired6juvs$GENE <- getGene(toshiny.paired6juvs$V_CALL, first=TRUE, strip_d=TRUE)
toshiny.paired6juvs$GF <- substring(toshiny.paired6juvs$GENE, 1,5)
toshiny.paired6juvs$JGENE <- getGene(toshiny.paired6juvs$J_CALL, first=TRUE, strip_d=TRUE)
toshiny.paired6juvs$JGF <- substring(toshiny.paired6juvs$JGENE, 1,5)

# toshiny.paired6juvs$JGF <- toshiny.paired6juvs$JGENE
# toshiny.paired6juvs$GF <- str_sub(toshiny.paired6juvs$GENE, end=5)
toshiny.paired6juvs$VJ_JUNCTION_PATTERN <- paste(toshiny.paired6juvs$GF,toshiny.paired6juvs$JGF,toshiny.paired6juvs$JUNCTION_LENGTH,sep="_") 

### check V_IDENTITY IF 0-100 OR 0-1 - WE WANT SHM TO BE BETWEEN 0-100
# toshiny.paired6juvs$SHM <- toshiny.paired6juvs$V_IDENTITY * 100
toshiny.paired6juvs$SHM <- (100 - (toshiny.paired6juvs$V_IDENTITY))
# toshiny.paired6juvs$SHM <- (100 - (toshiny.paired6juvs$V_IDENTITY * 100))
toshiny.paired6juvs <- toshiny.paired6juvs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

toshiny.paired6juvs$JUNCTION_LENGTH <- NULL
toshiny.paired6juvs$V_IDENTITY <- NULL

toshiny.paired6juvs$V_CALL <- NULL
toshiny.paired6juvs$J_CALL <- NULL

### NOTE OAS PAIRED HAS ISOTYPE & SUBTYPE!!!!!
## ADDING CREGION BY WAY OF GF...
toshiny.paired6juvs$CREGION <- str_sub(toshiny.paired6juvs$GF, end=3)
toshiny.paired6juvs$CREGION <- gsub("IGK","Kappa",toshiny.paired6juvs$CREGION)
toshiny.paired6juvs$CREGION <- gsub("IGL","Lambda",toshiny.paired6juvs$CREGION)
toshiny.paired6juvs$CREGION <- gsub("IGH","IgH",toshiny.paired6juvs$CREGION)

toshiny.paired6juvs <- toshiny.paired6juvs %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

rm(toshiny.paired6juvs.h)
rm(toshiny.paired6juvs.l)
toshiny.paired6juvs.h <- subset(toshiny.paired6juvs, CREGION %in% c("IgH"))

#write.table(toshiny.paired6juvs, "toshiny_paired6juvs.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.paired6juvs.h, "toshiny_paired6juvs_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# ggplot(toshiny.paired6juvs.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ SUBTYPE, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.paired6juvs.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.paired6juvs.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

# ggplot(toshiny.paired6juvs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

######################################
## next combine toshiny.fewhivnih45 & toshiny.pairedhivnih45

toshiny.hivnih45 <- full_join(toshiny.fewhivnih45, toshiny.pairedhivnih45)

## for combined need to recalculate ncount and shm.mean
toshiny.hivnih45$ncount <- NULL
toshiny.hivnih45$shm.mean <- NULL

toshiny.hivnih45 <- toshiny.hivnih45 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.hivnih45.h <- subset(toshiny.hivnih45, CREGION %in% c("IgH"))

#write.table(toshiny.hivnih45, "toshiny_hivnih45all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.hivnih45.h, "toshiny_hivnih45all_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# ggplot(toshiny.hivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ SUBTYPE, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.hivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.hivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


# ggplot(toshiny.hivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.pairedhivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


######################################
## combining new paired HIV data with existing HIV data...

## but first need to load N152 & IAVI84 data, then combine RD1214, newly combined nih45 & hiv mabs...
## but note IAVI84 is only HV1 & HV4 - N152 has HV3 & HV4  also note 2 HIV datasets (N152 & IAVI84) are included in toshiny.fewerhcbghiv2

## so instead just re-running HIV combination steps but with original & paired NIH45 data - toshiny.hivnih45

toshiny.fewhivmabsnih45rd1214 <- bind_rows(toshiny.fewhivmabs.h, toshiny.hivnih45.h, toshiny.fewhiv, .id = "id")

toshiny.fewhivmabsnih45rd1214$id <- gsub("1","anti-HIV mAbs",toshiny.fewhivmabsnih45rd1214$id)
toshiny.fewhivmabsnih45rd1214$id <- gsub("2","HIV+ patient NIH45",toshiny.fewhivmabsnih45rd1214$id)
toshiny.fewhivmabsnih45rd1214$id <- gsub("3","HIV+ patient AD1214",toshiny.fewhivmabsnih45rd1214$id)

toshiny.fewhivmabsnih45rd1214$id <- factor(toshiny.fewhivmabsnih45rd1214$id, levels = c("HIV+ patient AD1214", "HIV+ patient NIH45", "anti-HIV mAbs"))

## for combined need to recalculate ncount and shm.mean
toshiny.fewhivmabsnih45rd1214c <- toshiny.fewhivmabsnih45rd1214
toshiny.fewhivmabsnih45rd1214c$ncount <- NULL
toshiny.fewhivmabsnih45rd1214c$shm.mean <- NULL

toshiny.fewhivmabsnih45rd1214c <- toshiny.fewhivmabsnih45rd1214c %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.fewhivmabsnih45rd1214, "toshiny_mabsnih45rd1214.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhivmabsnih45rd1214c, "toshiny_mabsnih45rd1214c.tab", sep = "\t", row.names = FALSE, quote = FALSE)
ggplot(toshiny.fewhivmabsnih45rd1214, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.fewhivmabsnih45rd1214, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
######## to get comparison plots

### processing of small set recently run through immcantation:
# BX.clusters.toprocess <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_toprocess.tab")
# 
# BX.clusters.toprocess$V_CALL <- as.character(BX.clusters.toprocess$V_CALL)
# BX.clusters.toprocess$J_CALL <- as.character(BX.clusters.toprocess$J_CALL)
# 
# BX.clusters.toprocess$JUNCTIONAA <- translateDNA(BX.clusters.toprocess$JUNCTION, trim=TRUE)
# BX.clusters.toprocess$GENE <- getGene(BX.clusters.toprocess$V_CALL, first=TRUE, strip_d=TRUE)
# BX.clusters.toprocess$GF <- substring(BX.clusters.toprocess$GENE, 1,5)
# BX.clusters.toprocess$JGENE <- getGene(BX.clusters.toprocess$J_CALL, first=TRUE, strip_d=TRUE)
# BX.clusters.toprocess$JGF <- substring(BX.clusters.toprocess$JGENE, 1,5)
# BX.clusters.toprocess$VJ_JUNCTION_PATTERN <- paste(BX.clusters.toprocess$GF,BX.clusters.toprocess$JGF,BX.clusters.toprocess$JUNCTION_LENGTH,sep="_") 
# write.table(BX.clusters.toprocess, "COVID_mablist_toprocess2.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# 


### 950 anti-CoV2 mabs

# BX.clusters.mabs <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_clonalclusters.tab")
# BX.clusters.mabs$CDR3KABAT <- ((BX.clusters.mabs$JUNCTION_LENGTH) / 3)
# BX.clusters.mabs$CDR3LENGTH_IMGT <- ((BX.clusters.mabs$JUNCTION_LENGTH) / 3) - 2
# #BX.clusters.mabs <- BX.clusters.mabs %>% mutate( CDR3_bins = cut( CDR3LENGTH_IMGT, breaks = c(0,8,10,12,16,20,24,28) )) %>% unite(GF_CDR3, GF, CDR3_bins, sep = "_", remove = FALSE, na.rm = FALSE)
# BX.clusters.mabs$JGENE[BX.clusters.mabs$JGENE==""]<-NA
# BX.clusters.mabs <- BX.clusters.mabs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)
# BX.clusters.mabs$CREGION <- gsub("IGK","Kappa",BX.clusters.mabs$CREGION)
# BX.clusters.mabs$CREGION <- gsub("IGL","Lambda",BX.clusters.mabs$CREGION)
# BX.clusters.mabs$CREGION <- gsub("IGH","IgH",BX.clusters.mabs$CREGION)
# #BX.clusters.mabs$reads_per_clone <- 1
# BX.clusters.mabs$SHM <- (100 - (BX.clusters.mabs$V_IDENTITY * 100))
# 
# 
# BX.clusters.mabs$neutralization <- gsub("weakneutralizing","neutralizing",BX.clusters.mabs$neutralization)
# BX.clusters.mabs$neutralization <- gsub("neutralizingbroad","neutralizing",BX.clusters.mabs$neutralization)
# 
# ###
# ## NOW ADDTHE COUNTS OF COMBINED GF_JGENE,CDR3LENGTH_IMGT
# ## BUT HAVE TO DO IT BEFORE COMBINING DATASETS...
# BX.clusters.mabs <- BX.clusters.mabs %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   rename(ncount = n)
# 
# ## THIS ADDS MEAN SHM FOR EACH COMBINED GF_JGENE,CDR3LENGTH_IMGT 
# BX.clusters.mabs.h <- BX.clusters.mabs.h %>%
#   group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   mutate(shm.mean = mean(SHM, na.rm = TRUE))

## cleaned-up commands for single dataset input

BX.clusters.mabs <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_clonalclusters.tab")
#BX.clusters.mabs$CDR3KABAT <- ((BX.clusters.mabs$JUNCTION_LENGTH) / 3)
BX.clusters.mabs$CDR3LENGTH_IMGT <- ((BX.clusters.mabs$JUNCTION_LENGTH) / 3) - 2
#BX.clusters.mabs$JGENE[BX.clusters.mabs$JGENE==""]<-NA
# toshiny.few$SEQUENCE_ID <- as.character(toshiny.few$SEQUENCE_ID)
# toshiny.few$JUNCTIONAA <- as.character(toshiny.few$JUNCTIONAA)

BX.clusters.mabs$SHM <- (100 - (BX.clusters.mabs$V_IDENTITY * 100))
BX.clusters.mabs <- BX.clusters.mabs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
#BX.allclusters.b1 <- BX.allclusters.b1 %>% filter(is.wholenumber(CDR3LENGTH_IMGT))




## specific to dataset
#BX.clusters.mabs$reads_per_clone <- 1
BX.clusters.mabs$neutralization <- gsub("weakneutralizing","neutralizing",BX.clusters.mabs$neutralization)
BX.clusters.mabs$neutralization <- gsub("neutralizingbroad","neutralizing",BX.clusters.mabs$neutralization)
BX.clusters.mabs$CREGION <- gsub("IGK","Kappa",BX.clusters.mabs$CREGION)
BX.clusters.mabs$CREGION <- gsub("IGL","Lambda",BX.clusters.mabs$CREGION)
BX.clusters.mabs$CREGION <- gsub("IGH","IgH",BX.clusters.mabs$CREGION)
#BX.clusters.mabs.h$neutralization <- str_replace(BX.clusters.mabs.h$neutralization,"","non-neutralizing")
BX.clusters.mabs$neutralization <- sub("^$", "non-neutralizing", BX.clusters.mabs$neutralization)
BX.clusters.mabs.h$neutralization <- sub("^$", "non-neutralizing", BX.clusters.mabs.h$neutralization)


BX.clusters.mabs$neutralization <- sub("non-neutralizing", "non-neutralizing or unknown", BX.clusters.mabs$neutralization)
BX.clusters.mabs.h$neutralization <- sub("non-neutralizing", "non-neutralizing or unknown", BX.clusters.mabs.h$neutralization)

## rerunning after changing name of Augmenta for easier removal

BX.clusters.mabs$CREGION <- gsub("IGK","Kappa",BX.clusters.mabs$CREGION)
BX.clusters.mabs$CREGION <- gsub("IGL","Lambda",BX.clusters.mabs$CREGION)
BX.clusters.mabs$CREGION <- gsub("IGH","IgH",BX.clusters.mabs$CREGION)



BX.clusters.fewermabs <- BX.clusters.mabs %>% filter(source != "Augmenta")




BX.clusters.mabs.h <- subset(BX.clusters.mabs, CREGION %in% c("IgH"))
BX.clusters.fewermabs.h <- subset(BX.clusters.fewermabs, CREGION %in% c("IgH"))

# BX.clusters.mabs.h <- BX.clusters.mabs.h %>%
#   group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   mutate(shm.mean = mean(SHM, na.rm = TRUE))
toshiny.few1 <- BX.clusters.mabs %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.few1.h <- BX.clusters.mabs.h %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.fewer1 <- BX.clusters.fewermabs %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.fewer1.h <- BX.clusters.fewermabs.h %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))


######################################
### 2 COVID positive bulk datasets for comparisons 9/15

## cleaned-up commands for single dataset input

BX.clusters.g1 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/Galson1_clonalclusters.tab")
BX.clusters.1 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/Boyd_results/Boyd_7450ab10_clonalclusters.tab")

#BX.clusters.g1$CDR3KABAT <- ((BX.clusters.g1$JUNCTION_LENGTH) / 3)
BX.clusters.g1$CDR3LENGTH_IMGT <- ((BX.clusters.g1$JUNCTION_LENGTH) / 3) - 2
#BX.clusters.g1$JGENE[BX.clusters.g1$JGENE==""]<-NA
BX.clusters.g1$SHM <- (100 - (BX.clusters.g1$V_IDENTITY * 100))
BX.clusters.g1 <- BX.clusters.g1 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

## specific to dataset
#BX.clusters.g1$reads_per_clone <- 1
BX.clusters.g1$CREGION <- "IgH"

#BX.clusters.1$CDR3KABAT <- ((BX.clusters.1$JUNCTION_LENGTH) / 3)
BX.clusters.1$CDR3LENGTH_IMGT <- ((BX.clusters.1$JUNCTION_LENGTH) / 3) - 2
#BX.clusters.1$JGENE[BX.clusters.1$JGENE==""]<-NA
BX.clusters.1$SHM <- (100 - (BX.clusters.1$V_IDENTITY * 100))
BX.clusters.1 <- BX.clusters.1 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

## specific to dataset
#BX.clusters.1$reads_per_clone <- 1
BX.clusters.1$CREGION <- "IgH"

## first reduce dimensions, then combine toshiny.few tables??

toshiny.fewg1 <- BX.clusters.g1 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.fewb1 <- BX.clusters.1 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

### combining Boyd, Galson & HC1 & mabs

#BX.clusters.hcvsboydvsgalson <- bind_rows(BX.clusters.hc, BX.clusters.1, BX.clusters.g1, .id = "id")
BX.clusters.hcvsboydvsgalson <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, .id = "id")
BX.clusters.hcvsboydvsgalsonmabs <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, .id = "id")
toshiny.fewhcbg <- subset(BX.clusters.hcvsboydvsgalsonmabs, CREGION %in% c("IgH"))

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

  #write.table(toshiny.fewhcbgc, "toshiny_fewhcbgc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## for combining full set
## use toshiny.few3b as well...


####################################
### HIV positive bulk datasets for comparisons 9/28
## cleaned-up commands for single dataset input

## new combination because IAVI84 is only HV1 & HV4 - N152 has HV3 & HV4 (so still no HV2 or HV5-7)
BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/hiv_iavi84and152hc.tab")

## this removes first and last aa from junction (because ireceptor includes these) BE CAREFUL NOT TO EVER RUN THIS ON LC DATA...
# str_sub(BX.clusters.hiv$JUNCTIONAA, 1, 1) <- ""
# str_sub(BX.clusters.hiv$JUNCTIONAA, -1, -1) <- ""

BX.clusters.hiv$GENE <- getGene(BX.clusters.hiv$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$GF <- substring(BX.clusters.hiv$GENE, 1,5)
BX.clusters.hiv$JGENE <- getGene(BX.clusters.hiv$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$JGF <- substring(BX.clusters.hiv$JGENE, 1,5)
BX.clusters.hiv$VJ_JUNCTION_PATTERN <- paste(BX.clusters.hiv$GF,BX.clusters.hiv$JGF,BX.clusters.hiv$JUNCTION_LENGTH,sep="_") 

BX.clusters.hiv$CDR3LENGTH_IMGT <- ((BX.clusters.hiv$JUNCTION_LENGTH) / 3) - 2
#BX.clusters.hiv$JGENE[BX.clusters.hiv$JGENE==""]<-NA
BX.clusters.hiv$SHM <- (100 - (BX.clusters.hiv$V_IDENTITY * 100))
BX.clusters.hiv <- BX.clusters.hiv %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

BX.clusters.hiv <- BX.clusters.hiv %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
BX.clusters.hiv$CREGION <- "IgH"
toshiny.fewhiv <- BX.clusters.hiv %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))


####################################
### combining Boyd, Galson & HC1 & mabs & hiv

## need to rename b1 & g1 before combining
#write.table(toshiny.fewg1, "toshiny_fewg1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewb1, "toshiny_fewb1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
toshiny.fewg1 <- read.delim("toshiny_fewg1.tab")
toshiny.fewb1 <- read.delim("toshiny_fewb1.tab")


#BX.clusters.hcvsboydvsgalson <- bind_rows(BX.clusters.hc, BX.clusters.1, BX.clusters.g1, .id = "id")
#BX.clusters.hcvsboydvsgalson <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, .id = "id")
BX.clusters.hcvsboydvsgalsonmabshiv <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, toshiny.fewhiv, .id = "id")
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

#write.table(toshiny.fewhcbghiv, "toshiny_fewhcbghiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhcbghivc, "toshiny_fewhcbghivc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## removing some entries missing Jgene

#toshiny.fewhiv <- toshiny.fewhiv[ grep("__", toshiny.fewhiv$VJ_JUNCTION_PATTERN, invert = TRUE) , ]
#toshiny.fewhcbghivc <- toshiny.fewhcbghivc[ grep("__", toshiny.fewhcbghivc$VJ_JUNCTION_PATTERN, invert = TRUE) , ]
#toshiny.fewhcbghiv <- toshiny.fewhcbghiv[ grep("__", toshiny.fewhcbghiv$VJ_JUNCTION_PATTERN, invert = TRUE) , ]


d858.spike.clones <- d858.spike %>%
  arrange(desc(DUPCOUNTmax),(-DUPCOUNT)) %>%
  group_by(CLONE) %>% slice(1)
d858.spike.clones <- as.data.frame(d858.spike.clones)
## THIS ISN'T WORKING CORRECTLY, SAYING ALL COUNTS ARE 1???
d858.spike.clones <- d858.spike.clones %>%
  add_count(CREGION)

## but this works
starwars2 <- starwars %>% add_count(species)
## and this works?
BX.clusters.nat <- BX.clusters.nat %>% add_count(VJ_JUNCTION_PATTERN)
## this explains why... not a pure dataframe..
# > class(BX.clusters.nat)
# [1] "data.frame"
# > class(d858.spike.clones)
# [1] "grouped_df" "tbl_df"     "tbl"        "data.frame"



####################################
### reading Natalia's 4 samples 9/28 (testing input of raw tab files...)
d858.spike.clones <- as.data.frame(d858.spike.clones)
d858.spike.clones$V_CALL <- as.character(d858.spike.clones$V_CALL)
d858.spike.clones$J_CALL <- as.character(d858.spike.clones$J_CALL)

d858.spike.clones$GENE <- getGene(d858.spike.clones$V_CALL, first=TRUE, strip_d=TRUE)
d858.spike.clones$GF <- substring(d858.spike.clones$GENE, 1,5)
d858.spike.clones$JGENE <- getGene(d858.spike.clones$J_CALL, first=TRUE, strip_d=TRUE)
d858.spike.clones$JGF <- substring(d858.spike.clones$JGENE, 1,5)
d858.spike.clones$VJ_JUNCTION_PATTERN <- paste(d858.spike.clones$GF,d858.spike.clones$JGF,d858.spike.clones$JUNCTION_LENGTH,sep="_") 
d858.spike.clones$CDR3LENGTH_IMGT <- ((d858.spike.clones$JUNCTION_LENGTH) / 3) - 2
d858.spike.clones$SHM <- (100 - (d858.spike.clones$V_IDENTITY * 100))
d858.spike.clones$JUNCTIONAA <- translateDNA(d858.spike.clones$JUNCTION, trim=TRUE)
d858.spike.clones <- d858.spike.clones %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
toshiny.d858 <- d858.spike.clones %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

####
d903.spike.clones <- as.data.frame(d903.spike.clones)
d903.spike.clones$V_CALL <- as.character(d903.spike.clones$V_CALL)
d903.spike.clones$J_CALL <- as.character(d903.spike.clones$J_CALL)

d903.spike.clones$GENE <- getGene(d903.spike.clones$V_CALL, first=TRUE, strip_d=TRUE)
d903.spike.clones$GF <- substring(d903.spike.clones$GENE, 1,5)
d903.spike.clones$JGENE <- getGene(d903.spike.clones$J_CALL, first=TRUE, strip_d=TRUE)
d903.spike.clones$JGF <- substring(d903.spike.clones$JGENE, 1,5)
d903.spike.clones$VJ_JUNCTION_PATTERN <- paste(d903.spike.clones$GF,d903.spike.clones$JGF,d903.spike.clones$JUNCTION_LENGTH,sep="_") 
d903.spike.clones$CDR3LENGTH_IMGT <- ((d903.spike.clones$JUNCTION_LENGTH) / 3) - 2
d903.spike.clones$SHM <- (100 - (d903.spike.clones$V_IDENTITY * 100))
d903.spike.clones$JUNCTIONAA <- translateDNA(d903.spike.clones$JUNCTION, trim=TRUE)
d903.spike.clones <- d903.spike.clones %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
toshiny.d903 <- d903.spike.clones %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

###
d968.spike.clones <- as.data.frame(d968.spike.clones)
d968.spike.clones$V_CALL <- as.character(d968.spike.clones$V_CALL)
d968.spike.clones$J_CALL <- as.character(d968.spike.clones$J_CALL)

d968.spike.clones$GENE <- getGene(d968.spike.clones$V_CALL, first=TRUE, strip_d=TRUE)
d968.spike.clones$GF <- substring(d968.spike.clones$GENE, 1,5)
d968.spike.clones$JGENE <- getGene(d968.spike.clones$J_CALL, first=TRUE, strip_d=TRUE)
d968.spike.clones$JGF <- substring(d968.spike.clones$JGENE, 1,5)
d968.spike.clones$VJ_JUNCTION_PATTERN <- paste(d968.spike.clones$GF,d968.spike.clones$JGF,d968.spike.clones$JUNCTION_LENGTH,sep="_") 
d968.spike.clones$CDR3LENGTH_IMGT <- ((d968.spike.clones$JUNCTION_LENGTH) / 3) - 2
d968.spike.clones$SHM <- (100 - (d968.spike.clones$V_IDENTITY * 100))
d968.spike.clones$JUNCTIONAA <- translateDNA(d968.spike.clones$JUNCTION, trim=TRUE)
d968.spike.clones <- d968.spike.clones %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
toshiny.d968 <- d968.spike.clones %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))


###
d969.spike.clones <- as.data.frame(d969.spike.clones)
d969.spike.clones$V_CALL <- as.character(d969.spike.clones$V_CALL)
d969.spike.clones$J_CALL <- as.character(d969.spike.clones$J_CALL)

d969.spike.clones$GENE <- getGene(d969.spike.clones$V_CALL, first=TRUE, strip_d=TRUE)
d969.spike.clones$GF <- substring(d969.spike.clones$GENE, 1,5)
d969.spike.clones$JGENE <- getGene(d969.spike.clones$J_CALL, first=TRUE, strip_d=TRUE)
d969.spike.clones$JGF <- substring(d969.spike.clones$JGENE, 1,5)
d969.spike.clones$VJ_JUNCTION_PATTERN <- paste(d969.spike.clones$GF,d969.spike.clones$JGF,d969.spike.clones$JUNCTION_LENGTH,sep="_") 
d969.spike.clones$CDR3LENGTH_IMGT <- ((d969.spike.clones$JUNCTION_LENGTH) / 3) - 2
d969.spike.clones$SHM <- (100 - (d969.spike.clones$V_IDENTITY * 100))
d969.spike.clones$JUNCTIONAA <- translateDNA(d969.spike.clones$JUNCTION, trim=TRUE)
d969.spike.clones <- d969.spike.clones %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
toshiny.d969 <- d969.spike.clones %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

## combining 4 natalia samples
toshiny.natalia4 <- bind_rows(toshiny.few1.h, toshiny.d858, toshiny.d903, toshiny.d968, toshiny.d969, .id = "id")

## be careful with 3,2 & 1 in name here....
toshiny.natalia4$id <- gsub("4","d968",toshiny.natalia4$id)
toshiny.natalia4$id <- gsub("5","d969",toshiny.natalia4$id)
toshiny.natalia4$id <- gsub("3","d903",toshiny.natalia4$id)
toshiny.natalia4$id <- gsub("2","d858",toshiny.natalia4$id)
toshiny.natalia4$id <- gsub("1","anti-CoV2 mAbs",toshiny.natalia4$id)

#write.table(toshiny.natalia4, "toshiny_natalia.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## for combined need to recalculate ncount and shm.mean



#### to reload all correct tab files per app-2 and app-5
toshiny.few <- read.delim("toshiny_few.tab")
toshiny.few1 <- read.delim("toshiny_few1.tab")
toshiny.few1.h <- read.delim("toshiny_few1_h.tab")
toshiny.few2 <- read.delim("toshiny_few2.tab")
toshiny.few3 <- read.delim("toshiny_few3.tab")
toshiny.few3b <- read.delim("toshiny_few3b.tab")
toshiny.fewhcbg <- read.delim("toshiny_fewhcbg.tab")
toshiny.fewhcbgc <- read.delim("toshiny_fewhcbgc.tab")
toshiny.fewhcbghiv <- read.delim("toshiny_fewhcbghiv.tab")
toshiny.fewhcbghivc <- read.delim("toshiny_fewhcbghivc.tab")
toshiny.natalia4 <- read.delim("toshiny_natalia.tab")

toshiny.fewer <- read.delim("toshiny_fewer.tab")
toshiny.fewer1 <- read.delim("toshiny_fewer1.tab")
toshiny.fewer1.h <- read.delim("toshiny_fewer1_h.tab")
# toshiny.few2 <- read.delim("toshiny_few2.tab")
# toshiny.few3 <- read.delim("toshiny_few3.tab")
toshiny.fewerhcbg <- read.delim("toshiny_fewerhcbg.tab")
toshiny.fewerhcbgc <- read.delim("toshiny_fewerhcbgc.tab")


########################################## RD1214 HIV dataset

BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/binder_data/AD1214HC_germ-pass_clones.tab")

## unique to this dataset
BX.clusters.hiv <- BX.clusters.hiv %>% filter(JUNCTION_LENGTH > 4)

BX.clusters.hiv$V_CALL <- as.character(BX.clusters.hiv$V_CALL)
BX.clusters.hiv$J_CALL <- as.character(BX.clusters.hiv$J_CALL)

BX.clusters.hiv$GENE <- getGene(BX.clusters.hiv$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$GF <- substring(BX.clusters.hiv$GENE, 1,5)
BX.clusters.hiv$JGENE <- getGene(BX.clusters.hiv$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$JGF <- substring(BX.clusters.hiv$JGENE, 1,5)
BX.clusters.hiv$VJ_JUNCTION_PATTERN <- paste(BX.clusters.hiv$GF,BX.clusters.hiv$JGF,BX.clusters.hiv$JUNCTION_LENGTH,sep="_") 

BX.clusters.hiv$CDR3LENGTH_IMGT <- ((BX.clusters.hiv$JUNCTION_LENGTH) / 3) - 2
#BX.clusters.hiv$JGENE[BX.clusters.hiv$JGENE==""]<-NA
BX.clusters.hiv$SHM <- (100 - (BX.clusters.hiv$V_IDENTITY * 100))
BX.clusters.hiv$JUNCTIONAA <- translateDNA(BX.clusters.hiv$JUNCTION, trim=TRUE)

BX.clusters.hiv <- BX.clusters.hiv %>% filter(CDR3LENGTH_IMGT > 4.8)
#write.table(BX.clusters.hiv, "AD1214HC_germ-pass_clones2.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#BX.clusters.hiv <- BX.clusters.hiv[ grep("*", BX.clusters.hiv$JUNCTIONAA, invert = TRUE) , ]
#BX.clusters.hiv$JGENE[BX.clusters.hiv$JGENE==""]<-NA
#BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/binder_data/AD1214HC_germ-pass_clones2.tab")

BX.clusters.hiv <- BX.clusters.hiv %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

BX.clusters.hiv <- BX.clusters.hiv %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
BX.clusters.hiv$CREGION <- "IgH"


toshiny.fewhiv <- BX.clusters.hiv %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.fewhiv <- toshiny.fewhiv[ grep("Kappa|Lambda|IGKJ|__", toshiny.fewhiv$JGENE, invert = TRUE) , ]
#write.table(toshiny.fewhiv, "toshiny_fewhiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)


## combining this hiv instead of previous hiv run...actually add both!
toshiny.fewhiv1 <- read.delim("toshiny_fewhiv1.tab")
toshiny.fewhiv1$ncount <- NULL
toshiny.fewhiv1$shm.mean <- NULL

toshiny.fewhiv1 <- toshiny.fewhiv1 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))


BX.clusters.hcvsboydvsgalsonmabshiv <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, toshiny.fewhiv1, .id = "id")
#toshiny.fewhcbghiv <- subset(BX.clusters.hcvsboydvsgalsonmabshiv, CREGION %in% c("IgH"))
toshiny.fewhcbghiv <- BX.clusters.hcvsboydvsgalsonmabshiv[ grep("Kappa|Lambda|IGKJ|__", BX.clusters.hcvsboydvsgalsonmabshiv$CREGION, invert = TRUE) , ]
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

#write.table(toshiny.fewhcbghiv, "toshiny_fewhcbghiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhcbghivc, "toshiny_fewhcbghivc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## note that first 2 HIV datasets (N152 & IAVI84) are in toshiny.fewhcbghiv2, while only AD1214 is in toshiny.fewhcbghiv3
## above toshiny.fewhcbghiv includes all 3 combined...


## fewer
### 1 hiv dataset (AD1214)
BX.clusters.hcvsboydvsgalsonmabshivfew <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, toshiny.fewer1.h, toshiny.fewhiv, .id = "id")
### 3 hiv datasets
BX.clusters.hcvsboydvsgalsonmabshivfew <- bind_rows(toshiny.few3, toshiny.fewb1, toshiny.fewg1, toshiny.fewer1.h, toshiny.fewhiv1, .id = "id")
#toshiny.fewhcbghiv <- subset(BX.clusters.hcvsboydvsgalsonmabshiv, CREGION %in% c("IgH"))
toshiny.fewerhcbghiv <- BX.clusters.hcvsboydvsgalsonmabshivfew[ grep("Kappa|Lambda|IGKJ|__", BX.clusters.hcvsboydvsgalsonmabshivfew$JGENE, invert = TRUE) , ]
toshiny.fewerhcbghiv <- BX.clusters.hcvsboydvsgalsonmabshivfew[ grep("Kappa|Lambda|IGKJ|__", BX.clusters.hcvsboydvsgalsonmabshivfew$CREGION, invert = TRUE) , ]
#rm(BX.clusters.hcvsboydvsgalsonmabshiv)

toshiny.fewerhcbghiv$id <- gsub("1","Healthy control",toshiny.fewerhcbghiv$id)
## be careful with 1 in name here....
toshiny.fewerhcbghiv$id <- gsub("4","anti-CoVII mAbs",toshiny.fewerhcbghiv$id)
toshiny.fewerhcbghiv$id <- gsub("2","Boyd7450",toshiny.fewerhcbghiv$id)
toshiny.fewerhcbghiv$id <- gsub("3","Galson1",toshiny.fewerhcbghiv$id)
toshiny.fewerhcbghiv$id <- gsub("anti-CoVII mAbs","anti-CoV2 mAbs",toshiny.fewerhcbghiv$id)
toshiny.fewerhcbghiv$id <- gsub("Boyd7450","Boyd74s0",toshiny.fewerhcbghiv$id)
toshiny.fewerhcbghiv$id <- gsub("5","HIV+ patient",toshiny.fewerhcbghiv$id)
toshiny.fewerhcbghiv$id <- gsub("Boyd74s0","Boyd7450",toshiny.fewerhcbghiv$id)

## for combined need to recalculate ncount and shm.mean
toshiny.fewerhcbghivc <- toshiny.fewerhcbghiv
toshiny.fewerhcbghivc$ncount <- NULL
toshiny.fewerhcbghivc$shm.mean <- NULL

toshiny.fewerhcbghivc <- toshiny.fewerhcbghivc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.fewerhcbghiv, "toshiny_fewerhcbghiv3.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewerhcbghivc, "toshiny_fewerhcbghiv3c.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(toshiny.fewerhcbghiv, "toshiny_fewerhcbghiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewerhcbghivc, "toshiny_fewerhcbghivc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## removing some entries missing Jgene

toshiny.fewhiv <- toshiny.fewhiv[ grep("__", toshiny.fewhiv$VJ_JUNCTION_PATTERN, invert = TRUE) , ]
toshiny.fewhiv1 <- toshiny.fewhiv1[ grep("__", toshiny.fewhiv1$VJ_JUNCTION_PATTERN, invert = TRUE) , ]
toshiny.fewhcbghivc <- toshiny.fewhcbghivc[ grep("__", toshiny.fewhcbghivc$VJ_JUNCTION_PATTERN, invert = TRUE) , ]
toshiny.fewhcbghiv <- toshiny.fewhcbghiv[ grep("__", toshiny.fewhcbghiv$VJ_JUNCTION_PATTERN, invert = TRUE) , ]


## need to recalculated combined for just HC & CoV2 mabs

## for combined need to recalculate ncount and shm.mean
toshiny.fewerc <- toshiny.fewer
toshiny.fewerc$ncount <- NULL
toshiny.fewerc$shm.mean <- NULL

toshiny.fewerc <- toshiny.fewerc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.fewerc$CREGION <- "IgH"

#write.table(toshiny.fewerc, "toshiny_fewerc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

################################################################################################################
################################################################################################################

### NIH45 HIV dataset (donor with VRC01)
## NOTE THIS DONOR HAS ONLY HV1 & HV4 IN LARGE NUMBERS...
BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/hiv-nih45/airr-nih45.tsv")

BX.clusters.hiva <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/hiv-nih45/airr-nih45a.tsv")
BX.clusters.hivat <- BX.clusters.hiva %>% filter(junction != "") %>% 
  filter(productive == "TRUE")

BX.clusters.hiva <- BX.clusters.hiva %>% filter(junction != "") %>% 
  filter(productive == "TRUE") %>% 
  sample_n(30000)


BX.clusters.hivb <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/hiv-nih45/airr-nih45b.tsv")

BX.clusters.hivb <- BX.clusters.hivb %>% filter(junction != "") %>% 
  filter(productive == "TRUE") %>% 
  sample_n(30000)

BX.clusters.hivc <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/hiv-nih45/airr-nih45c.tsv")

BX.clusters.hivc <- BX.clusters.hivc %>% filter(junction != "") %>% 
  filter(productive == "TRUE") %>% 
  sample_n(30000)

BX.clusters.hiv <- bind_rows(BX.clusters.hiva, BX.clusters.hivb, BX.clusters.hivc)

#write.table(BX.clusters.hiv, "airr-nih45thin90k.tab", sep = "\t", row.names = FALSE, quote = FALSE)
## adding manually ~100 VRC01 lineage:

BX.clusters.hivp <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/hiv-nih45/airr-nih45thin90k.tab")
#BX.clusters.hivp <- BX.clusters.hiv
BX.clusters.hivvrc01 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/hiv-nih45/airr-nih45vrc01.tsv")

BX.clusters.hiv <- bind_rows(BX.clusters.hivp, BX.clusters.hivvrc01)

#BX.clusters.hiva <- BX.clusters.hiva %>% filter(productive == "F")
## next thin, then repeat for a & b, then add all together
#sample_n(mtcars, 10)



## because ireceptor need to capitalize column names
#(SEQUENCE_ID,CREGION,JUNCTIONAA,,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM

## remove any igk, igl, blank junctions, also trim junctionaa
#BX.clusters.hiv <- BX.clusters.hiv[ grep("IGL|IGK", BX.clusters.hiv$v_call, invert = TRUE) , ]

BX.clusters.hiv <- BX.clusters.hiv %>% filter(junction != "")


## this removes first and last aa from junction (because ireceptor includes these) BE CAREFUL NOT TO EVER RUN THIS ON LC DATA...
str_sub(BX.clusters.hiv$junction_aa, 1, 1) <- ""
str_sub(BX.clusters.hiv$junction_aa, -1, -1) <- ""

# rename(iris, petal_length = Petal.Length) new.name = old.name

BX.clusters.hiv <- BX.clusters.hiv %>%
  rename(SEQUENCE_ID = sequence_id) %>%
  rename(V_CALL = v_call) %>%
  rename(J_CALL = j_call) %>%
  rename(JUNCTIONAA = junction_aa) %>%
  rename(V_IDENTITY = v_identity) %>%
  rename(JUNCTION = junction) %>%
  rename(JUNCTION_LENGTH = junction_aa_length)

## unique to this dataset
#BX.clusters.hiv <- BX.clusters.hiv %>% filter(JUNCTION_LENGTH > 4)

BX.clusters.hiv$V_CALL <- as.character(BX.clusters.hiv$V_CALL)
BX.clusters.hiv$J_CALL <- as.character(BX.clusters.hiv$J_CALL)

BX.clusters.hiv$GENE <- getGene(BX.clusters.hiv$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$GF <- substring(BX.clusters.hiv$GENE, 1,5)
BX.clusters.hiv$JGENE <- getGene(BX.clusters.hiv$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$JGF <- substring(BX.clusters.hiv$JGENE, 1,5)
BX.clusters.hiv$VJ_JUNCTION_PATTERN <- paste(BX.clusters.hiv$GF,BX.clusters.hiv$JGF,BX.clusters.hiv$JUNCTION_LENGTH,sep="_") 

### UNIQUE TO THIS DATASET DO NOT USE UNLESS YOU ARE SURE OF JUNCTION DEFINITIONS
BX.clusters.hiv$junction_length <- NULL
BX.clusters.hiv$CDR3LENGTH_IMGT <- (BX.clusters.hiv$JUNCTION_LENGTH) - 2
#BX.clusters.hiv$JGENE[BX.clusters.hiv$JGENE==""]<-NA
## also already multiplied by 100 here
BX.clusters.hiv$SHM <- (100 - BX.clusters.hiv$V_IDENTITY)
#BX.clusters.hiv$JUNCTIONAA <- translateDNA(BX.clusters.hiv$JUNCTION, trim=TRUE)

BX.clusters.hiv <- BX.clusters.hiv %>% filter(CDR3LENGTH_IMGT > 4.8)
#write.table(BX.clusters.hiv, "AD1214HC_germ-pass_clones2.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#BX.clusters.hiv <- BX.clusters.hiv[ grep("*", BX.clusters.hiv$JUNCTIONAA, invert = TRUE) , ]
#BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/binder_data/AD1214HC_germ-pass_clones2.tab")

BX.clusters.hiv <- BX.clusters.hiv %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

BX.clusters.hiv <- BX.clusters.hiv %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
BX.clusters.hiv$CREGION <- "IgH"

## NOTE THIS DONOR HAS ONLY HV1 & HV4 IN LARGE NUMBERS...
toshiny.fewhivnih45 <- BX.clusters.hiv %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

ggplot(toshiny.fewhivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.fewhivnih45, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

#toshiny.fewhiv <- toshiny.fewhiv[ grep("Kappa|Lambda|IGKJ|__", toshiny.fewhiv$JGENE, invert = TRUE) , ]
#write.table(toshiny.fewhivnih45, "toshiny_fewhivnih45.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## in Excel renaming id's, then reload
toshiny.fewhivnih45 <- read.delim("toshiny_fewhivnih45.tab")


########################################
### first list of anti-HIV bnAbs..
BX.clusters.hivmabs <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/HIV_mablist_clonalclusters.tab")

### just to collapse multiple rows & get median IC50s...
gm_mean = function(a){prod(a)^(1/length(a))}
or
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
# BX.clusters.hivmabsforIC50 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/HIV_mablist_forIC50s.tab")
# BX.clusters.hivmabsforIC50b <- BX.clusters.hivmabsforIC50 %>%
#   group_by(Antibody) %>%
#   # summarize(median_IC50 = median(IC50))
#   summarize(median_IC50 = gm_mean(IC50))

#write.table(BX.clusters.hivmabsforIC50b, "HIV_mablist_forIC50s_medians.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## unique to this dataset

BX.clusters.hivmabs$JUNCTIONAA <- as.character(BX.clusters.hivmabs$JUNCTIONAA)
BX.clusters.hivmabs$V_CALL <- as.character(BX.clusters.hivmabs$V_CALL)
BX.clusters.hivmabs$J_CALL <- as.character(BX.clusters.hivmabs$J_CALL)

BX.clusters.hivmabs$CDR3LENGTH_IMGT <- nchar(BX.clusters.hivmabs$JUNCTIONAA)

BX.clusters.hivmabs$JUNCTION_LENGTH <- (BX.clusters.hivmabs$CDR3LENGTH_IMGT) + 2


BX.clusters.hivmabs$GENE <- getGene(BX.clusters.hivmabs$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hivmabs$GF <- substring(BX.clusters.hivmabs$GENE, 1,5)
BX.clusters.hivmabs$JGENE <- getGene(BX.clusters.hivmabs$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hivmabs$JGF <- substring(BX.clusters.hivmabs$JGENE, 1,5)
BX.clusters.hivmabs$VJ_JUNCTION_PATTERN <- paste(BX.clusters.hivmabs$GF,BX.clusters.hivmabs$JGF,BX.clusters.hivmabs$JUNCTION_LENGTH,sep="_") 

### UNIQUE TO THIS DATASET DO NOT USE UNLESS YOU ARE SURE OF JUNCTION DEFINITIONS
## also already multiplied by 100 here
BX.clusters.hivmabs$SHM <- (100 - BX.clusters.hivmabs$V_IDENTITY)
#BX.clusters.hivmabs$JUNCTIONAA <- translateDNA(BX.clusters.hivmabs$JUNCTION, trim=TRUE)


BX.clusters.hivmabs <- BX.clusters.hivmabs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

BX.clusters.hivmabs <- BX.clusters.hivmabs %>% filter(is.wholenumber(CDR3LENGTH_IMGT))


BX.clusters.hivmabs$CREGION <- gsub("IGK","Kappa",BX.clusters.hivmabs$CREGION)
BX.clusters.hivmabs$CREGION <- gsub("IGL","Lambda",BX.clusters.hivmabs$CREGION)
BX.clusters.hivmabs$CREGION <- gsub("IGH","IgH",BX.clusters.hivmabs$CREGION)

toshiny.fewhivmabs <- BX.clusters.hivmabs %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

ggplot(toshiny.fewhivmabs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.fewhivmabs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#write.table(toshiny.fewhivmabs, "toshiny_hivmabs.tab", sep = "\t", row.names = FALSE, quote = FALSE)

toshiny.fewhivmabs.h <- subset(toshiny.fewhivmabs, CREGION %in% c("IgH"))


########################################################
### combining NIH45, RD1214, list of anti-HIV bnAbs..

toshiny.fewhiv <- toshiny.fewhiv[ grep("Kappa|Lambda|IGKJ|__", toshiny.fewhiv$JGENE, invert = TRUE) , ]

toshiny.fewhivmabsnih45rd1214 <- bind_rows(toshiny.fewhivmabs.h, toshiny.fewhivnih45, toshiny.fewhiv, .id = "id")

toshiny.fewhivmabsnih45rd1214$id <- gsub("1","anti-HIV mAbs",toshiny.fewhivmabsnih45rd1214$id)
toshiny.fewhivmabsnih45rd1214$id <- gsub("2","HIV+ patient NIH45",toshiny.fewhivmabsnih45rd1214$id)
toshiny.fewhivmabsnih45rd1214$id <- gsub("3","HIV+ patient AD1214",toshiny.fewhivmabsnih45rd1214$id)

toshiny.fewhivmabsnih45rd1214$id <- factor(toshiny.fewhivmabsnih45rd1214$id, levels = c("HIV+ patient AD1214", "HIV+ patient NIH45", "anti-HIV mAbs"))

## for combined need to recalculate ncount and shm.mean
toshiny.fewhivmabsnih45rd1214c <- toshiny.fewhivmabsnih45rd1214
toshiny.fewhivmabsnih45rd1214c$ncount <- NULL
toshiny.fewhivmabsnih45rd1214c$shm.mean <- NULL

toshiny.fewhivmabsnih45rd1214c <- toshiny.fewhivmabsnih45rd1214c %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.fewhivmabsnih45rd1214, "toshiny_mabsnih45rd1214.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhivmabsnih45rd1214c, "toshiny_mabsnih45rd1214c.tab", sep = "\t", row.names = FALSE, quote = FALSE)


ggplot(toshiny.fewhivmabsnih45rd1214, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.fewhivmabsnih45rd1214, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))



########################################################
## and anti-HIV vs. anti-COV2 vs. anti-SARSMERS

toshiny.fewhivmabs$CREGION <- gsub("IGK","Kappa",toshiny.fewhivmabs$CREGION)
toshiny.fewhivmabs$CREGION <- gsub("IGL","Lambda",toshiny.fewhivmabs$CREGION)
toshiny.fewhivmabs$CREGION <- gsub("IGH","IgH",toshiny.fewhivmabs$CREGION)
toshiny.fewhivmabs.h <- subset(toshiny.fewhivmabs, CREGION %in% c("IgH"))

toshiny.sarsmers.h <- read.delim("toshiny_sarsmers_h.tab")

list2env(x = split(x = toshiny.sarsmers.h,
                   f = toshiny.sarsmers.h$source),
         envir = globalenv())

toshiny.sars.h <- SARS
toshiny.mers.h <- MERS
rm(SARS)
rm(MERS)



toshiny.hivcov2sarsmers <- bind_rows(toshiny.few1.h, toshiny.sars.h, toshiny.mers.h, toshiny.fewhivmabs.h, .id = "id")

toshiny.hivcov2sarsmers$id <- gsub("2","anti-SARS mAbs",toshiny.hivcov2sarsmers$id)
toshiny.hivcov2sarsmers$id <- gsub("3","anti-MERS mAbs",toshiny.hivcov2sarsmers$id)
toshiny.hivcov2sarsmers$id <- gsub("4","anti-HIV mAbs",toshiny.hivcov2sarsmers$id)
## because has 2 needs to be last
toshiny.hivcov2sarsmers$id <- gsub("1","anti-CoV2 mAbs",toshiny.hivcov2sarsmers$id)

toshiny.hivcov2sarsmers$id <- factor(toshiny.hivcov2sarsmers$id, levels = c("anti-CoV2 mAbs", "anti-SARS mAbs", "anti-MERS mAbs", "anti-HIV mAbs"))

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


ggplot(toshiny.hivcov2sarsmers, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.hivcov2sarsmers, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


######################################################################################################
######################################################################################################

########
# Dengue datasets

## dengue
#   363214 OAS_dengue_germ-pass.tab
# plasmablast list ALSO ADD ANY GENERIC MABS????  DENmAbs.tab
OAS_dengue_germ-pass2.tab

d13_2_enriched_UMI2high_germ-pass2.tab
d13_2_iggiga_UMI2high_germ-pass2.tab
d13_3_iggiga_UMI2high_germ-pass2.tab
d13_3_enriched_UMI2high_germ-pass2.tab
d20_2_enriched_UMI2high_germ-pass2.tab
d20_2_iggiga_UMI2high_germ-pass2.tab
d20_3_iggiga_UMI2high_germ-pass2.tab
d20_3_enriched_UMI2high_germ-pass2.tab

## these above are combined in jan 2021 - but as of mar 2021 we want to NOT INCLUDE the conv _3 datasest
## we also want to INCLUDE the stim0 datasets - BUT MAYBE FIRST REMOVING THE OFFENDING IDENTICAL CONVERGENCES?
d13_2_stimulated0ss_germ-pass22.tab
d20_2_stimulated0ss_germ-pass22.tab
## REMOVE THE DC38 SEQUENCES TOO

## JAN 2021 NOTE - WHEN I MADE THESE PASS2.TAB FILES, I REMOVED THE FULL SEQUENCE!! NEED TO RECOVER - now pass22.tab
### CAN MAKE EQUALLY SMALL BY INSTEAD REMOVING FWR1_IMGT FWR2_IMGT FWR3_IMGT FWR4_IGMT CDR1_IMGT CDR2_IMGT
### ALSO NOTE BELOW I RECODE ALL CREGION TO IGH, WANT TO NOT DO THAT OR MAYBE ALSO ADD NEW COLUMN CREGION0 WITH ORIGINAL DATA...

## MARCH UPDATE NO LONGER ADDING THE CONV _3 DATASETS, ADDING INITIAL STIM DATASETS
BX.clusters.dengue.d13.2 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d13_2_enriched_UMI2high_germ-pass23.tab")
BX.clusters.dengue.d13.2iggiga <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d13_2_iggiga_UMI2high_germ-pass23.tab")
# BX.clusters.dengue.d13.3 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d13_3_enriched_UMI2high_germ-pass23.tab")
# BX.clusters.dengue.d13.3iggiga <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d13_3_iggiga_UMI2high_germ-pass23.tab")

BX.clusters.dengue.d20.2 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d20_2_enriched_UMI2high_germ-pass23.tab")
BX.clusters.dengue.d20.2iggiga <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d20_2_iggiga_UMI2high_germ-pass23.tab")
# BX.clusters.dengue.d20.3 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d20_3_enriched_UMI2high_germ-pass23.tab")
# BX.clusters.dengue.d20.3iggiga <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/d20_3_iggiga_UMI2high_germ-pass23.tab")

BX.clusters.dengue.d13.2stim <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/feb2021_newssdensearches/d13_2_stimulated0ss_germ-pass23.tab")
BX.clusters.dengue.d20.2stim <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/feb2021_newssdensearches/d20_2_stimulated0ss_germ-pass23.tab")
### OAS COMMANDS 200 LINES AHEAD

# BX.clusters.dengue.mabs <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/DENmAbs.tab")
# BX.clusters.dengue.oas <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/OAS_dengue_germ-pass2.tab")

### LATE NOV 2020 - NEED TO COLLAPSE IDENTICAL JUNCTIONAA FIRST (ALSO FOR OAS), OTHERWISE TOO MANY SEQUENCES...
# distinct_at() will be superseded by use of across() inside distinct() from version 1.0.0 (dplyr.tidyverse.org/news/index.html#across). 
# The equivalent pattern for your answer would be dat %>% distinct(across(-z)) but distinct_at() will still be available for several years. 
## to make sure it keeps largest CONSCOUNT though - group_by(CONSCOUNT) %>%
BX.clusters.dengue.d13.2$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d13.2$JUNCTION, trim=TRUE)
# BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2 %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

BX.clusters.dengue.d13.2iggiga$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d13.2iggiga$JUNCTION, trim=TRUE)
# BX.clusters.dengue.d13.2iggiga <- BX.clusters.dengue.d13.2iggiga %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

BX.clusters.dengue.d20.2$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d20.2$JUNCTION, trim=TRUE)
# BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2 %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

BX.clusters.dengue.d20.2iggiga$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d20.2iggiga$JUNCTION, trim=TRUE)
# BX.clusters.dengue.d20.2iggiga <- BX.clusters.dengue.d20.2iggiga %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

BX.clusters.dengue.d13.2stim$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d13.2stim$JUNCTION, trim=TRUE)
BX.clusters.dengue.d20.2stim$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d20.2stim$JUNCTION, trim=TRUE)


BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2 %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
BX.clusters.dengue.d13.2iggiga <- BX.clusters.dengue.d13.2iggiga %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2 %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
BX.clusters.dengue.d20.2iggiga <- BX.clusters.dengue.d20.2iggiga %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
BX.clusters.dengue.d13.2stim <- BX.clusters.dengue.d13.2stim %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

## not using mar 2021
# BX.clusters.dengue.d13.3$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d13.3$JUNCTION, trim=TRUE)
# # BX.clusters.dengue.d13.3 <- BX.clusters.dengue.d13.3 %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
# BX.clusters.dengue.d13.3iggiga$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d13.3iggiga$JUNCTION, trim=TRUE)
# # BX.clusters.dengue.d13.3iggiga <- BX.clusters.dengue.d13.3iggiga %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
# 
# BX.clusters.dengue.d20.3$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d20.3$JUNCTION, trim=TRUE)
# # BX.clusters.dengue.d20.3 <- BX.clusters.dengue.d20.3 %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
# BX.clusters.dengue.d20.3iggiga$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d20.3iggiga$JUNCTION, trim=TRUE)
# # BX.clusters.dengue.d20.3iggiga <- BX.clusters.dengue.d20.3iggiga %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

### OR OLDER COMMANDS TO COLLAPSE - BUT GROUP BY WILL CHANGE SO JUNCTIONAA IS NOW PRIMARY COLUMN...
### creating lists per clone not per read
# clonestatsh1 <- BX_hobs %>%
#   group_by(CLONE) %>%
#   summarize_at(c("PRCONS2","GENE","V_CALL","D_CALL","J_CALL","JUNCTION_LENGTH","JUNCTION_LENGTH2","CDRH3KABAT_LENGTH","FAMILY"), first)
# clonestatsh2 <- BX_hobs %>%
#   group_by(CLONE) %>%
#   summarize_if(is.numeric, mean)
# clonestatsh3 <- BX_hobs %>%
#   group_by(CLONE) %>%
#   summarize_at(c("PRCONS2"), n_distinct)
# clonestatsh5 <- BX_hobs %>%
#   group_by(CLONE) %>%
#   summarize_at(c("PRCONS2"), Mode)

## MAR 2021 - BEFORE COMBINING BELOW, LOOK AT THE DATASETS SEPARATELY...

BX.clusters.dengue.d13.2iggiga <- BX.clusters.dengue.d13.2iggiga[ grep("IgD|Kappa|Lambda", BX.clusters.dengue.d13.2iggiga$CREGION, invert = TRUE) , ]
# BX.clusters.dengue.d13.2iggiga <- BX.clusters.dengue.d13.2iggiga[ grep("IGLV|IGKV", BX.clusters.dengue.d13.2iggiga$V_CALL, invert = TRUE) , ]

# toshiny.dengue.d13.2iggiga <- toshiny.dengue.d13.2iggiga[ grep("IGHV7", toshiny.dengue.d13.2iggiga$GF_JGENE, invert = TRUE) , ]
# toshiny.dengue.d13.2stim <- toshiny.dengue.d13.2stim[ grep("IGHV7", toshiny.dengue.d13.2stim$GF_JGENE, invert = TRUE) , ]
# toshiny.dengue.d13.2iggiga <- toshiny.dengue.d13.2iggiga[ grep("ATGDVRGDAIVRKRSPARSTGWPRRPEEDH", toshiny.dengue.d13.2iggiga$JUNCTIONAA, invert = TRUE) , ]

# BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2[ grep("IgD|Kappa|Lambda", BX.clusters.dengue.d20.2$CREGION, invert = TRUE) , ]
# BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2[ grep("IGLV|IGKV", BX.clusters.dengue.d20.2$GF_JGENE, invert = TRUE) , ]

## NEXT FIND THE CONVERGENT SEQUENCES AND REMOVE...
BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2[ grep("ARALFGLVAVASPFDN", BX.clusters.dengue.d13.2$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2[ grep("ITPPLYLMVGGVSRAMAV", BX.clusters.dengue.d13.2$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2[ grep("ARQDRNWFDT", BX.clusters.dengue.d13.2$JUNCTIONAA, invert = TRUE) , ]

BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2[ grep("ARGPGGTSTSCYHCWFDP", BX.clusters.dengue.d20.2$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2[ grep("AKNYGSGTLNWFDS", BX.clusters.dengue.d20.2$JUNCTIONAA, invert = TRUE) , ]


BX.clusters.dengue.d13.2stim <- BX.clusters.dengue.d13.2stim[ grep("DC38", BX.clusters.dengue.d13.2stim$SEQUENCE_ID, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("DC38", BX.clusters.dengue.d20.2stim$SEQUENCE_ID, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARADEMATVQGFYAFDI", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARADEMATIEGFYAFDI", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARADEMATIEGFYAFGI", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]

BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARGPGGTTTSCYHCWFDP", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARGPGGTSSSCYQCWFDP", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("AKNYGSGTLNWFDS", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]

BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARGFATTQWQGHNWFDP", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("AKDVGECSGGNCFSGYFYYMDA", BX.clusters.dengue.d20.2stim$JUNCTIONAA, invert = TRUE) , ]

######################

## need to run this for each of the 6 datasets separately for now
BX.clusters.dengue.d13.2stim$V_CALL <- as.character(BX.clusters.dengue.d13.2stim$V_CALL)
BX.clusters.dengue.d13.2stim$J_CALL <- as.character(BX.clusters.dengue.d13.2stim$J_CALL)
BX.clusters.dengue.d13.2stim$V_CALL <- as.character(BX.clusters.dengue.d13.2stim$V_CALL)
BX.clusters.dengue.d13.2stim <- BX.clusters.dengue.d13.2stim %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

BX.clusters.dengue.d13.2stim$GENE <- getGene(BX.clusters.dengue.d13.2stim$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.d13.2stim$GF <- substring(BX.clusters.dengue.d13.2stim$GENE, 1,5)
BX.clusters.dengue.d13.2stim$JGENE <- getGene(BX.clusters.dengue.d13.2stim$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.d13.2stim$JGF <- substring(BX.clusters.dengue.d13.2stim$JGENE, 1,5)
BX.clusters.dengue.d13.2stim$VJ_JUNCTION_PATTERN <- paste(BX.clusters.dengue.d13.2stim$GF,BX.clusters.dengue.d13.2stim$JGF,BX.clusters.dengue.d13.2stim$JUNCTION_LENGTH,sep="_") 
BX.clusters.dengue.d13.2stim$SHM <- (100 - (BX.clusters.dengue.d13.2stim$V_IDENTITY * 100))
## NOTE UNIQUE TO MY IMMCANTATION RUNS...
BX.clusters.dengue.d13.2stim$CDR3LENGTH_IMGT <- ((BX.clusters.dengue.d13.2stim$JUNCTION_LENGTH) / 3) - 2

BX.clusters.dengue.d13.2stim <- BX.clusters.dengue.d13.2stim %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
toshiny.dengue.d13.2stim <- BX.clusters.dengue.d13.2stim %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,SEQUENCE_INPUT) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))


### NEW COMBINATIONS
toshiny.dengue.d13a <- full_join(toshiny.dengue.d13.2, toshiny.dengue.d13.2iggiga)
toshiny.dengue.d13 <- full_join(toshiny.dengue.d13a, toshiny.dengue.d13.2stim)
toshiny.dengue.d20a <- full_join(toshiny.dengue.d20.2, toshiny.dengue.d20.2iggiga)
toshiny.dengue.d20 <- full_join(toshiny.dengue.d20a, toshiny.dengue.d20.2stim)
rm(toshiny.dengue.d13a)
rm(toshiny.dengue.d20a)


## NEED TO RECALCULATE SHM.MEAN & NCOUNT

toshiny.dengue.d13$ncount <- NULL
toshiny.dengue.d13$shm.mean <- NULL
toshiny.dengue.d13 <- toshiny.dengue.d13 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.dengue.d20$ncount <- NULL
toshiny.dengue.d20$shm.mean <- NULL
toshiny.dengue.d20 <- toshiny.dengue.d20 %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

## KEEPING CREGION0 AS INITIAL DATA BUT CHANGING ALL TO IGH FOR CREGION
toshiny.dengue.d13$CREGION0 <- toshiny.dengue.d13$CREGION
toshiny.dengue.d20$CREGION0 <- toshiny.dengue.d20$CREGION

toshiny.dengue.d13.h <- subset(toshiny.dengue.d13, CREGION %in% c("IgM","IgG","IgA"))
toshiny.dengue.d20.h <- subset(toshiny.dengue.d20, CREGION %in% c("IgM","IgG","IgA"))


toshiny.dengue.d13.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.d13.h$CREGION)
toshiny.dengue.d13.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.d13.h$CREGION)
toshiny.dengue.d13.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.d13.h$CREGION)

toshiny.dengue.d20.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.d20.h$CREGION)
toshiny.dengue.d20.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.d20.h$CREGION)
toshiny.dengue.d20.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.d20.h$CREGION)

#write.table(toshiny.dengue.d13, "toshiny_dengue_d13.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.d20, "toshiny_dengue_d20.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#write.table(toshiny.dengue.d13.h, "toshiny_dengue_d13_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.d20.h, "toshiny_dengue_d20_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)


## NOW COMBINING MABS & OAS
toshiny.dengue.oas.h <- read.delim("toshiny_dengue_oas.tab")


toshiny.dengueall <- bind_rows(toshiny.dengue.mabs.h, toshiny.dengue.oas.h, toshiny.dengue.d13.h, toshiny.dengue.d20.h, .id = "id")

toshiny.dengueall$id <- gsub("2","DEN patients OAS data",toshiny.dengueall$id)
toshiny.dengueall$id <- gsub("3","DEN patient thirteen",toshiny.dengueall$id)
toshiny.dengueall$id <- gsub("4","DEN patient twenty",toshiny.dengueall$id)
## because has 2 needs to be last
toshiny.dengueall$id <- gsub("1","DEN plasmablasts",toshiny.dengueall$id)

toshiny.dengueall$id <- gsub("DEN patient thirteen","DEN patient 13",toshiny.dengueall$id)
toshiny.dengueall$id <- gsub("DEN patient twenty","DEN patient 20",toshiny.dengueall$id)

toshiny.dengueall$id <- factor(toshiny.dengueall$id, levels = c("DEN plasmablasts", "DEN patient 13", "DEN patient 20", "DEN patients OAS data"))

## for combined need to recalculate ncount and shm.mean
toshiny.dengueallc <- toshiny.dengueall
toshiny.dengueallc$ncount <- NULL
toshiny.dengueallc$shm.mean <- NULL

toshiny.dengueallc <- toshiny.dengueallc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.dengueall, "toshiny_dengueall.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengueallc, "toshiny_dengueallc.tab", sep = "\t", row.names = FALSE, quote = FALSE)
## note these are now 352k, was 370k, but much bigger filesize due to full sequences in there...

ggplot(toshiny.dengueall, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengueall, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


######################################################################################
### d13.2 iggiga remove IGKV IGLV

toshiny.dengue.d13.2iggiga <- toshiny.dengue.d13.2iggiga[ grep("IgD|Kappa|Lambda", toshiny.dengue.d13.2iggiga$CREGION, invert = TRUE) , ]
toshiny.dengue.d13.2iggiga <- toshiny.dengue.d13.2iggiga[ grep("IGLV|IGKV", toshiny.dengue.d13.2iggiga$GF_JGENE, invert = TRUE) , ]

# toshiny.dengue.d13.2iggiga <- toshiny.dengue.d13.2iggiga[ grep("IGHV7", toshiny.dengue.d13.2iggiga$GF_JGENE, invert = TRUE) , ]
# toshiny.dengue.d13.2stim <- toshiny.dengue.d13.2stim[ grep("IGHV7", toshiny.dengue.d13.2stim$GF_JGENE, invert = TRUE) , ]
# toshiny.dengue.d13.2iggiga <- toshiny.dengue.d13.2iggiga[ grep("ATGDVRGDAIVRKRSPARSTGWPRRPEEDH", toshiny.dengue.d13.2iggiga$JUNCTIONAA, invert = TRUE) , ]


toshiny.dengue.d20 <- toshiny.dengue.d20[ grep("IgD|Kappa|Lambda", toshiny.dengue.d20$CREGION, invert = TRUE) , ]
toshiny.dengue.d20 <- toshiny.dengue.d20[ grep("IGLV|IGKV", toshiny.dengue.d20$GF_JGENE, invert = TRUE) , ]


toshiny.dengue6separate <- bind_rows(toshiny.dengue.d13.2, toshiny.dengue.d13.2iggiga, toshiny.dengue.d13.2stim, toshiny.dengue.d20.2, toshiny.dengue.d20.2iggiga, toshiny.dengue.d20.2stim, .id = "id")

toshiny.dengue6separate$id <- gsub("1","DEN patient thirteen non-stim",toshiny.dengue6separate$id)
toshiny.dengue6separate$id <- gsub("2","DEN patient thirteen IgGIgA",toshiny.dengue6separate$id)
toshiny.dengue6separate$id <- gsub("3","DEN patient thirteen stim",toshiny.dengue6separate$id)
toshiny.dengue6separate$id <- gsub("4","DEN patient twenty non-stim",toshiny.dengue6separate$id)
toshiny.dengue6separate$id <- gsub("5","DEN patient twenty IgGIgA",toshiny.dengue6separate$id)
toshiny.dengue6separate$id <- gsub("6","DEN patient twenty stim",toshiny.dengue6separate$id)

toshiny.dengue6separate$id <- gsub("DEN patient thirteen","DEN patient 13",toshiny.dengue6separate$id)
toshiny.dengue6separate$id <- gsub("DEN patient twenty","DEN patient 20",toshiny.dengue6separate$id)

toshiny.dengue6separate$id <- factor(toshiny.dengue6separate$id, levels = c("DEN patient 13 non-stim", "DEN patient 13 IgGIgA", "DEN patient 13 stim", "DEN patient 20 non-stim", "DEN patient 20 IgGIgA", "DEN patient 20 stim"))


#write.table(toshiny.dengue6separate, "toshiny_dengue6separate.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue6separate, "toshiny_dengue6separate.tab", sep = "\t", row.names = FALSE, quote = FALSE)
ggplot(toshiny.dengue6separate, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue6separate, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


##################################################################################################################
##################################################################################################################
##################################################################################################################
################3
## then combine all d13 into 1 file, all d20 into 1 file - OLD COMBINATION
BX.clusters.dengue.d13a <- full_join(BX.clusters.dengue.d13.2, BX.clusters.dengue.d13.2iggiga)
BX.clusters.dengue.d13b <- full_join(BX.clusters.dengue.d13.3,BX.clusters.dengue.d13.3iggiga)
BX.clusters.dengue.d20a <- full_join(BX.clusters.dengue.d20.2, BX.clusters.dengue.d20.2iggiga)
BX.clusters.dengue.d20b <- full_join(BX.clusters.dengue.d20.3, BX.clusters.dengue.d20.3iggiga)
BX.clusters.dengue.d13 <- full_join(BX.clusters.dengue.d13a, BX.clusters.dengue.d13b)
BX.clusters.dengue.d20 <- full_join(BX.clusters.dengue.d20a, BX.clusters.dengue.d20b)
rm(BX.clusters.dengue.d13a)
rm(BX.clusters.dengue.d13b)
rm(BX.clusters.dengue.d20a)
rm(BX.clusters.dengue.d20b)

rm(BX.clusters.dengue.d13.2)
rm(BX.clusters.dengue.d13.2iggiga)
rm(BX.clusters.dengue.d13.3)
rm(BX.clusters.dengue.d13.3iggiga)
rm(BX.clusters.dengue.d20.2)
rm(BX.clusters.dengue.d20.2iggiga)
rm(BX.clusters.dengue.d20.3)
rm(BX.clusters.dengue.d20.3iggiga)

###


# BX.clones.hc <- collapseClones(BX.full0, method="mostCommon", cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT", germlineColumn = "GERMLINE_IMGT_D_MASK",
#                                includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
#BX.clusters.hc <- BX.clusters2.hc
# BX.clusters.hc3 <- groupGenes(BX.clones.hc, v_call = "V_CALL", j_call = "J_CALL",
#                               junc_len = "JUNCTION_LENGTH", cell_id = NULL, locus = NULL, only_igh = TRUE,
#                               first = TRUE)

BX.clusters.dengue.d13$V_CALL <- as.character(BX.clusters.dengue.d13$V_CALL)
BX.clusters.dengue.d13$J_CALL <- as.character(BX.clusters.dengue.d13$J_CALL)
BX.clusters.dengue.d13$V_CALL <- as.character(BX.clusters.dengue.d13$V_CALL)
BX.clusters.dengue.d13 <- BX.clusters.dengue.d13 %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

#BX.clusters.dengue.d13$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d13$JUNCTION, trim=TRUE)
BX.clusters.dengue.d13$GENE <- getGene(BX.clusters.dengue.d13$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.d13$GF <- substring(BX.clusters.dengue.d13$GENE, 1,5)
BX.clusters.dengue.d13$JGENE <- getGene(BX.clusters.dengue.d13$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.d13$JGF <- substring(BX.clusters.dengue.d13$JGENE, 1,5)
BX.clusters.dengue.d13$VJ_JUNCTION_PATTERN <- paste(BX.clusters.dengue.d13$GF,BX.clusters.dengue.d13$JGF,BX.clusters.dengue.d13$JUNCTION_LENGTH,sep="_") 
BX.clusters.dengue.d13$SHM <- (100 - (BX.clusters.dengue.d13$V_IDENTITY * 100))
## NOTE UNIQUE TO MY IMMCANTATION RUNS...
BX.clusters.dengue.d13$CDR3LENGTH_IMGT <- ((BX.clusters.dengue.d13$JUNCTION_LENGTH) / 3) - 2


BX.clusters.dengue.d13 <- BX.clusters.dengue.d13 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))


BX.clusters.dengue.d20$V_CALL <- as.character(BX.clusters.dengue.d20$V_CALL)
BX.clusters.dengue.d20$J_CALL <- as.character(BX.clusters.dengue.d20$J_CALL)
BX.clusters.dengue.d20$V_CALL <- as.character(BX.clusters.dengue.d20$V_CALL)
BX.clusters.dengue.d20 <- BX.clusters.dengue.d20 %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

#BX.clusters.dengue.d20$JUNCTIONAA <- translateDNA(BX.clusters.dengue.d20$JUNCTION, trim=TRUE)
BX.clusters.dengue.d20$GENE <- getGene(BX.clusters.dengue.d20$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.d20$GF <- substring(BX.clusters.dengue.d20$GENE, 1,5)
BX.clusters.dengue.d20$JGENE <- getGene(BX.clusters.dengue.d20$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.d20$JGF <- substring(BX.clusters.dengue.d20$JGENE, 1,5)
BX.clusters.dengue.d20$VJ_JUNCTION_PATTERN <- paste(BX.clusters.dengue.d20$GF,BX.clusters.dengue.d20$JGF,BX.clusters.dengue.d20$JUNCTION_LENGTH,sep="_") 
BX.clusters.dengue.d20$SHM <- (100 - (BX.clusters.dengue.d20$V_IDENTITY * 100))
BX.clusters.dengue.d20$CDR3LENGTH_IMGT <- ((BX.clusters.dengue.d20$JUNCTION_LENGTH) / 3) - 2

BX.clusters.dengue.d20 <- BX.clusters.dengue.d20 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

# BX.clusters.hiv <- BX.clusters.hiv %>% filter(CDR3LENGTH_IMGT > 4.8)
#write.table(BX.clusters.hiv, "AD1214HC_germ-pass_clones2.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#BX.clusters.hiv <- BX.clusters.hiv[ grep("*", BX.clusters.hiv$JUNCTIONAA, invert = TRUE) , ]

# BX.clusters.mabs.h <- BX.clusters.mabs.h %>%
#   group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   mutate(shm.mean = mean(SHM, na.rm = TRUE))

#write.table(BX.clusters.dengue.d13, "dengue_d13all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(BX.clusters.dengue.d20, "dengue_d20all.tab", sep = "\t", row.names = FALSE, quote = FALSE)

toshiny.dengue.d13 <- BX.clusters.dengue.d13 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,SEQUENCE_INPUT) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

## in d13, remove IgD, Kappa, Lambda

toshiny.dengue.d20 <- BX.clusters.dengue.d20 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,SEQUENCE_INPUT) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

## JAN 2021 NEW SELECT TO KEEP SEQUENCE_INPUT



toshiny.dengue.d13 <- toshiny.dengue.d13[ grep("IgD|Kappa|Lambda", toshiny.dengue.d13$CREGION, invert = TRUE) , ]
toshiny.dengue.d20 <- toshiny.dengue.d20[ grep("IgD|Kappa|Lambda", toshiny.dengue.d20$CREGION, invert = TRUE) , ]
toshiny.dengue.d13 <- toshiny.dengue.d13[ grep("IGLV|IGKV", toshiny.dengue.d13$GF_JGENE, invert = TRUE) , ]
toshiny.dengue.d20 <- toshiny.dengue.d20[ grep("IGLV|IGKV", toshiny.dengue.d20$GF_JGENE, invert = TRUE) , ]

### create cregion0??
toshiny.dengue.d13$CREGION0 <- toshiny.dengue.d13$CREGION
toshiny.dengue.d20$CREGION0 <- toshiny.dengue.d20$CREGION

toshiny.dengue.d13.h <- subset(toshiny.dengue.d13, CREGION %in% c("IgM","IgG","IgA"))
toshiny.dengue.d20.h <- subset(toshiny.dengue.d20, CREGION %in% c("IgM","IgG","IgA"))


toshiny.dengue.d13.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.d13.h$CREGION)
toshiny.dengue.d13.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.d13.h$CREGION)
toshiny.dengue.d13.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.d13.h$CREGION)

toshiny.dengue.d20.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.d20.h$CREGION)
toshiny.dengue.d20.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.d20.h$CREGION)
toshiny.dengue.d20.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.d20.h$CREGION)

ggplot(toshiny.dengue.d13.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue.d13.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

ggplot(toshiny.dengue.d20.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue.d20.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

ggplot(toshiny.dengue.d20, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

#write.table(toshiny.dengue.d20, "toshiny_dengue_d20full.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.d13, "toshiny_dengue_d13full.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.d20.h, "toshiny_dengue_d20.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.d13.h, "toshiny_dengue_d13.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# plasmablast list ALSO ADD ANY GENERIC MABS????  DENmAbs.tab

BX.clusters.dengue.mabs <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/DENmAbs.tab")

BX.clusters.dengue.mabs$V_CALL <- as.character(BX.clusters.dengue.mabs$V_CALL)
BX.clusters.dengue.mabs$J_CALL <- as.character(BX.clusters.dengue.mabs$J_CALL)
BX.clusters.dengue.mabs$V_CALL <- as.character(BX.clusters.dengue.mabs$V_CALL)
BX.clusters.dengue.mabs <- BX.clusters.dengue.mabs %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

BX.clusters.dengue.mabs$JUNCTIONAA <- translateDNA(BX.clusters.dengue.mabs$JUNCTION, trim=TRUE)
BX.clusters.dengue.mabs$GENE <- getGene(BX.clusters.dengue.mabs$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.mabs$GF <- substring(BX.clusters.dengue.mabs$GENE, 1,5)
BX.clusters.dengue.mabs$JGENE <- getGene(BX.clusters.dengue.mabs$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.mabs$JGF <- substring(BX.clusters.dengue.mabs$JGENE, 1,5)
BX.clusters.dengue.mabs$VJ_JUNCTION_PATTERN <- paste(BX.clusters.dengue.mabs$GF,BX.clusters.dengue.mabs$JGF,BX.clusters.dengue.mabs$JUNCTION_LENGTH,sep="_") 
BX.clusters.dengue.mabs$SHM <- (100 - (BX.clusters.dengue.mabs$V_IDENTITY * 100))
BX.clusters.dengue.mabs$CDR3LENGTH_IMGT <- ((BX.clusters.dengue.mabs$JUNCTION_LENGTH) / 3) - 2

BX.clusters.dengue.mabs <- BX.clusters.dengue.mabs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))



#toshiny.dengue.mabs <- BX.clusters.dengue.mabs %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

ggplot(toshiny.dengue.mabs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue.mabs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

#### then OAS data

BX.clusters.dengue.oas <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/OAS_dengue_germ-pass23.tab")

BX.clusters.dengue.oas$JUNCTIONAA <- translateDNA(BX.clusters.dengue.oas$JUNCTION, trim=TRUE)
#BX.clusters.dengue.oas <- BX.clusters.dengue.oas %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

BX.clusters.dengue.oas$V_CALL <- as.character(BX.clusters.dengue.oas$V_CALL)
BX.clusters.dengue.oas$J_CALL <- as.character(BX.clusters.dengue.oas$J_CALL)
BX.clusters.dengue.oas$V_CALL <- as.character(BX.clusters.dengue.oas$V_CALL)
BX.clusters.dengue.oas <- BX.clusters.dengue.oas %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

#BX.clusters.dengue.oas$JUNCTIONAA <- translateDNA(BX.clusters.dengue.oas$JUNCTION, trim=TRUE)
BX.clusters.dengue.oas$GENE <- getGene(BX.clusters.dengue.oas$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.oas$GF <- substring(BX.clusters.dengue.oas$GENE, 1,5)
BX.clusters.dengue.oas$JGENE <- getGene(BX.clusters.dengue.oas$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.oas$JGF <- substring(BX.clusters.dengue.oas$JGENE, 1,5)
BX.clusters.dengue.oas$VJ_JUNCTION_PATTERN <- paste(BX.clusters.dengue.oas$GF,BX.clusters.dengue.oas$JGF,BX.clusters.dengue.oas$JUNCTION_LENGTH,sep="_") 
BX.clusters.dengue.oas$SHM <- (100 - (BX.clusters.dengue.oas$V_IDENTITY * 100))
BX.clusters.dengue.oas$CDR3LENGTH_IMGT <- ((BX.clusters.dengue.oas$JUNCTION_LENGTH) / 3) - 2

### need to remove some bad sequences
BX.clusters.dengue.oas2 <- BX.clusters.dengue.oas %>% filter(CDR3LENGTH_IMGT > 4.8)
BX.clusters.dengue.oas2 <- BX.clusters.dengue.oas2 %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)

#BX.clusters.dengue.oas2 <- BX.clusters.dengue.oas2[ grep("\\*", BX.clusters.dengue.oas2$JUNCTIONAA, invert = TRUE) , ]
#BX.clusters.dengue.oas2$JGENE[BX.clusters.dengue.oas2$JGENE==""]<-NA


BX.clusters.dengue.oas2 <- BX.clusters.dengue.oas2 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))



toshiny.dengue.oas <- BX.clusters.dengue.oas2 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,SEQUENCE_INPUT) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

rm(BX.clusters.dengue.oas)

ggplot(toshiny.dengue.oas, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue.oas, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

### then combining the dengue datasets...ALSO THE 38 PLASMABLAST SEQUENCES...

#toshiny.dengue.oas.h <- subset(toshiny.dengue.oas, CREGION %in% c("IgH"))

toshiny.dengue.oas$CREGION0 <- toshiny.dengue.oas$CREGION

toshiny.dengue.oas.h <- subset(toshiny.dengue.oas, CREGION %in% c("IgG"))
toshiny.dengue.oas.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.oas.h$CREGION)

rm(toshiny.dengue.oas)
#write.table(toshiny.dengue.oas.h, "toshiny_dengue_oas.tab", sep = "\t", row.names = FALSE, quote = FALSE)

toshiny.dengue.oas.h <- read.delim("toshiny_dengue_oas.tab")


toshiny.dengue.mabs <- BX.clusters.dengue.mabs %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,SEQUENCE_INPUT) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.dengue.mabs$CREGION0 <- toshiny.dengue.mabs$CREGION
toshiny.dengue.mabs.h <- subset(toshiny.dengue.mabs, CREGION %in% c("IgG"))
toshiny.dengue.mabs.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.mabs.h$CREGION)


toshiny.dengue.d13.h <- toshiny.dengue.d13.h %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
toshiny.dengue.d20.h <- toshiny.dengue.d20.h %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
#write.table(toshiny.dengue.d20.h, "toshiny_dengue_d20.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.d13.h, "toshiny_dengue_d13.tab", sep = "\t", row.names = FALSE, quote = FALSE)


toshiny.dengueall <- bind_rows(toshiny.dengue.mabs.h, toshiny.dengue.oas.h, toshiny.dengue.d13.h, toshiny.dengue.d20.h, .id = "id")

toshiny.dengueall$id <- gsub("2","DEN patients OAS data",toshiny.dengueall$id)
toshiny.dengueall$id <- gsub("3","DEN patient thirteen",toshiny.dengueall$id)
toshiny.dengueall$id <- gsub("4","DEN patient twenty",toshiny.dengueall$id)
## because has 2 needs to be last
toshiny.dengueall$id <- gsub("1","DEN plasmablasts",toshiny.dengueall$id)

toshiny.dengueall$id <- gsub("DEN patient thirteen","DEN patient 13",toshiny.dengueall$id)
toshiny.dengueall$id <- gsub("DEN patient twenty","DEN patient 20",toshiny.dengueall$id)

toshiny.dengueall$id <- factor(toshiny.dengueall$id, levels = c("DEN plasmablasts", "DEN patient 13", "DEN patient 20", "DEN patients OAS data"))

## for combined need to recalculate ncount and shm.mean
toshiny.dengueallc <- toshiny.dengueall
toshiny.dengueallc$ncount <- NULL
toshiny.dengueallc$shm.mean <- NULL

toshiny.dengueallc <- toshiny.dengueallc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.dengueall, "toshiny_dengueall.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengueallc, "toshiny_dengueallc.tab", sep = "\t", row.names = FALSE, quote = FALSE)
## note these are now 352k, was 370k, but much bigger filesize due to full sequences in there...

ggplot(toshiny.dengueall, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengueall, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## too big to keep around!!
rm(toshiny.dengueall)
rm(toshiny.dengueallc)


toshiny.dengueall0 <- read.delim("toshiny_dengueall0.tab")

###############################################3
## some dengue network pre-processing with joins:

## also CONSCOUNT
toshiny.dengue.mabsmore <- BX.clusters.dengue.mabs %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,CONSCOUNT,SEQUENCE_IMGT,GERMLINE_IMGT_D_MASK,V_CALL,J_CALL,JUNCTION_LENGTH) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.dengue.mabsmore$CREGION0 <- toshiny.dengue.mabsmore$CREGION
toshiny.dengue.mabsmore.h <- subset(toshiny.dengue.mabsmore, CREGION %in% c("IgG"))
toshiny.dengue.mabsmore.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.mabsmore.h$CREGION)

toshiny.dengue.oasmore <- BX.clusters.dengue.oas2 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,CONSCOUNT,SEQUENCE_IMGT,GERMLINE_IMGT_D_MASK,V_CALL,J_CALL,JUNCTION_LENGTH) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.dengue.oasmore$CREGION0 <- toshiny.dengue.oasmore$CREGION
toshiny.dengue.oasmore.h <- subset(toshiny.dengue.oasmore, CREGION %in% c("IgG"))
toshiny.dengue.oasmore.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.oasmore.h$CREGION)

#write.table(toshiny.dengue.oasmore.h, "toshiny_dengue_oas23.tab", sep = "\t", row.names = FALSE, quote = FALSE)

rm(toshiny.dengue.oasmore)



toshiny.dengue.d13more <- BX.clusters.dengue.d13 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,CONSCOUNT,SEQUENCE_IMGT,GERMLINE_IMGT_D_MASK,V_CALL,J_CALL,JUNCTION_LENGTH) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

## in d13, remove IgD, Kappa, Lambda

toshiny.dengue.d20more <- BX.clusters.dengue.d20 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,CONSCOUNT,SEQUENCE_IMGT,GERMLINE_IMGT_D_MASK,V_CALL,J_CALL,JUNCTION_LENGTH) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

## JAN 2021 NEW SELECT TO KEEP SEQUENCE_INPUT



toshiny.dengue.d13more <- toshiny.dengue.d13more[ grep("IgD|Kappa|Lambda", toshiny.dengue.d13more$CREGION, invert = TRUE) , ]
toshiny.dengue.d20more <- toshiny.dengue.d20more[ grep("IgD|Kappa|Lambda", toshiny.dengue.d20more$CREGION, invert = TRUE) , ]
toshiny.dengue.d13more <- toshiny.dengue.d13more[ grep("IGLV|IGKV", toshiny.dengue.d13more$GF_JGENE, invert = TRUE) , ]
toshiny.dengue.d20more <- toshiny.dengue.d20more[ grep("IGLV|IGKV", toshiny.dengue.d20more$GF_JGENE, invert = TRUE) , ]

### create cregion0??
toshiny.dengue.d13more$CREGION0 <- toshiny.dengue.d13more$CREGION
toshiny.dengue.d20more$CREGION0 <- toshiny.dengue.d20more$CREGION

toshiny.dengue.d13more.h <- subset(toshiny.dengue.d13more, CREGION %in% c("IgM","IgG","IgA"))
toshiny.dengue.d20more.h <- subset(toshiny.dengue.d20more, CREGION %in% c("IgM","IgG","IgA"))


toshiny.dengue.d13more.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.d13more.h$CREGION)
toshiny.dengue.d13more.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.d13more.h$CREGION)
toshiny.dengue.d13more.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.d13more.h$CREGION)

toshiny.dengue.d20more.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.d20more.h$CREGION)
toshiny.dengue.d20more.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.d20more.h$CREGION)
toshiny.dengue.d20more.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.d20more.h$CREGION)
toshiny.dengue.d13more.h <- toshiny.dengue.d13more.h %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)
toshiny.dengue.d20more.h <- toshiny.dengue.d20more.h %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)


####
toshiny.dengueallmore <- bind_rows(toshiny.dengue.mabsmore.h, toshiny.dengue.oasmore.h, toshiny.dengue.d13more.h, toshiny.dengue.d20more.h, .id = "id")

toshiny.dengueallmore$id <- gsub("2","DEN patients OAS data",toshiny.dengueallmore$id)
toshiny.dengueallmore$id <- gsub("3","DEN patient thirteen",toshiny.dengueallmore$id)
toshiny.dengueallmore$id <- gsub("4","DEN patient twenty",toshiny.dengueallmore$id)
## because has 2 needs to be last
toshiny.dengueallmore$id <- gsub("1","DEN plasmablasts",toshiny.dengueallmore$id)

toshiny.dengueallmore$id <- gsub("DEN patient thirteen","DEN patient 13",toshiny.dengueallmore$id)
toshiny.dengueallmore$id <- gsub("DEN patient twenty","DEN patient 20",toshiny.dengueallmore$id)

toshiny.dengueallmore$id <- factor(toshiny.dengueallmore$id, levels = c("DEN plasmablasts", "DEN patient 13", "DEN patient 20", "DEN patients OAS data"))

#write.table(toshiny.dengueallmore, "toshiny_dengue_allwith2sequencecols.tab", sep = "\t", row.names = FALSE, quote = FALSE)

rm(toshiny.dengue.oasmore.h)
rm(toshiny.dengue.oas.h)
rm(toshiny.dengue.d13)
rm(toshiny.dengue.d13.h)
rm(toshiny.dengue.d20)
rm(toshiny.dengue.d20.h)

rm(toshiny.dengue.d13more)
rm(toshiny.dengue.d20more)


rm(toshiny.dengueallmore)
rm(toshiny.dengueallc)
rm(toshiny.dengueall)

## getting CF10 & CF7 (also CF4) extractions

#toshiny.dengue.cf10 <- subset(toshiny.dengue.mabs, CREGION %in% c("IgG"))
toshiny.dengue.cf10 <- toshiny.dengueallmore %>% filter(GF_JGENE == "IGHV4_IGHJ5" & CDR3LENGTH_IMGT == "10")
toshiny.dengue.cf7 <- toshiny.dengueallmore %>% filter(GF_JGENE == "IGHV1_IGHJ5" & CDR3LENGTH_IMGT == "16")

#write.table(toshiny.dengue.cf10, "toshiny_dengue_cf10.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.cf7, "toshiny_dengue_cf4andcf7.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#starwars <- starwars %>% filter(hair_color == "none" & eye_color == "black")




## NEED TO THIN TO SMALLER SUBSET THOUGH...inner join with subset from Geneious
##  dengue_CF10_list.tsv
CF10mergedwithOAS.list <- read.delim("dengue_CF10_list.tsv")

## next run some sort of join...semi-join semi_join() return all rows from x with a match in y.

CF10mergedwithOASfewer <- semi_join(toshiny.dengue.cf10, CF10mergedwithOAS.list)

CF4mergedwithOAS.list <- read.delim("dengue_CF4_list.tsv")

CF4mergedwithOASfewer <- semi_join(toshiny.dengue.cf7, CF4mergedwithOAS.list)

### networks of dengue cf10 & cf 7 - see far below 

######################################################################################################
######################################################################################################
### march 2021 - rechecking DEN vaccinated alloy mice
######################################################################################################

BXss.clusters.dengue.m2 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/mousedenv1234m2ss_germ-pass23.tab")

BXss.clusters.dengue.m2ss <- filter(BXss.clusters.dengue.m2, grepl('DC38', SEQUENCE_ID)) %>% select(CLONE) %>% unique()
BXss.clusters.dengue.m2shm <- BXss.clusters.dengue.m2 %>% add_count(CLONE) %>% rename(reads_per_clone = n) %>% group_by(CLONE) %>% summarize_at(c("V_IDENTITY", "reads_per_clone"), mean) %>% rename(V_IDENTITYmean = V_IDENTITY)
BXss.clusters.dengue.m2maxc <- BXss.clusters.dengue.m2 %>% group_by(CLONE) %>% summarize_at(c("CONSCOUNT"), max) %>% rename(CONSCOUNTmax = CONSCOUNT)
BXss.clusters.dengue.m2ssden0 <- inner_join(BXss.clusters.dengue.m2, BXss.clusters.dengue.m2ss)
BXss.clusters.dengue.m2ssden1 <- inner_join(BXss.clusters.dengue.m2ssden0, BXss.clusters.dengue.m2shm)
BXss.clusters.dengue.m2ssden <- inner_join(BXss.clusters.dengue.m2ssden1, BXss.clusters.dengue.m2maxc)
BXss.clusters.dengue.m2ssden$JUNCTIONAA <- translateDNA(BXss.clusters.dengue.m2ssden$JUNCTION, trim=TRUE)


BXss.clusters.dengue.m1 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/mousedenv1234m1ss_germ-pass23.tab")

BXss.clusters.dengue.m1ss <- filter(BXss.clusters.dengue.m1, grepl('DC38', SEQUENCE_ID)) %>% select(CLONE) %>% unique()
BXss.clusters.dengue.m1shm <- BXss.clusters.dengue.m1 %>% add_count(CLONE) %>% rename(reads_per_clone = n) %>% group_by(CLONE) %>% summarize_at(c("V_IDENTITY", "reads_per_clone"), mean) %>% rename(V_IDENTITYmean = V_IDENTITY)
BXss.clusters.dengue.m1maxc <- BXss.clusters.dengue.m1 %>% group_by(CLONE) %>% summarize_at(c("CONSCOUNT"), max) %>% rename(CONSCOUNTmax = CONSCOUNT)
BXss.clusters.dengue.m1ssden0 <- inner_join(BXss.clusters.dengue.m1, BXss.clusters.dengue.m1ss)
BXss.clusters.dengue.m1ssden1 <- inner_join(BXss.clusters.dengue.m1ssden0, BXss.clusters.dengue.m1shm)
BXss.clusters.dengue.m1ssden <- inner_join(BXss.clusters.dengue.m1ssden1, BXss.clusters.dengue.m1maxc)
BXss.clusters.dengue.m1ssden$JUNCTIONAA <- translateDNA(BXss.clusters.dengue.m1ssden$JUNCTION, trim=TRUE)

# write.table(BXss.clusters.dengue.m2ssden, "BXss.clusters.dengue.m2ssden.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(BXss.clusters.dengue.m1ssden, "BXss.clusters.dengue.m1ssden.tab", sep = "\t", row.names = FALSE, quote = FALSE)

rm(BXss.clusters.dengue.m2)
rm(BXss.clusters.dengue.m2maxc)
rm(BXss.clusters.dengue.m2shm)
rm(BXss.clusters.dengue.m2ssden0)
rm(BXss.clusters.dengue.m2ssden1)

rm(BXss.clusters.dengue.m1)
rm(BXss.clusters.dengue.m1maxc)
rm(BXss.clusters.dengue.m1shm)
rm(BXss.clusters.dengue.m1ssden0)
rm(BXss.clusters.dengue.m1ssden1)

#####################
### running below commands as both m1 and m2, then combine into new allplusmousedengue file...

BX.clusters.dengue.m1 <- read.delim("/users/eric.waltari/immcantation_pipeline/dengue/mousedenv1234m1ss_germ-pass23.tab")
BX.clusters.dengue.m1 <- BX.clusters.dengue.m1[ grep("DC38", BX.clusters.dengue.m1$SEQUENCE_ID, invert = TRUE) , ]


BX.clusters.dengue.m1$JUNCTIONAA <- translateDNA(BX.clusters.dengue.m1$JUNCTION, trim=TRUE)
BX.clusters.dengue.m1 <- BX.clusters.dengue.m1 %>% arrange((desc(CONSCOUNT))) %>% distinct(across(JUNCTIONAA), .keep_all = TRUE)


BX.clusters.dengue.m1$V_CALL <- as.character(BX.clusters.dengue.m1$V_CALL)
BX.clusters.dengue.m1$J_CALL <- as.character(BX.clusters.dengue.m1$J_CALL)
BX.clusters.dengue.m1$V_CALL <- as.character(BX.clusters.dengue.m1$V_CALL)
BX.clusters.dengue.m1 <- BX.clusters.dengue.m1 %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

#BX.clusters.dengue.m1$JUNCTIONAA <- translateDNA(BX.clusters.dengue.m1$JUNCTION, trim=TRUE)
BX.clusters.dengue.m1$GENE <- getGene(BX.clusters.dengue.m1$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.m1$GF <- substring(BX.clusters.dengue.m1$GENE, 1,5)
BX.clusters.dengue.m1$JGENE <- getGene(BX.clusters.dengue.m1$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.dengue.m1$JGF <- substring(BX.clusters.dengue.m1$JGENE, 1,5)
BX.clusters.dengue.m1$VJ_JUNCTION_PATTERN <- paste(BX.clusters.dengue.m1$GF,BX.clusters.dengue.m1$JGF,BX.clusters.dengue.m1$JUNCTION_LENGTH,sep="_") 
BX.clusters.dengue.m1$SHM <- (100 - (BX.clusters.dengue.m1$V_IDENTITY * 100))
BX.clusters.dengue.m1$CDR3LENGTH_IMGT <- ((BX.clusters.dengue.m1$JUNCTION_LENGTH) / 3) - 2

BX.clusters.dengue.m1 <- BX.clusters.dengue.m1 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))


toshiny.dengue.m1 <- BX.clusters.dengue.m1 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone,SEQUENCE_INPUT) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.dengue.m1 <- toshiny.dengue.m1[ grep("IgD|Kappa|Lambda", toshiny.dengue.m1$CREGION, invert = TRUE) , ]
toshiny.dengue.m1 <- toshiny.dengue.m1[ grep("IGLV|IGKV|IGHV8", toshiny.dengue.m1$GF_JGENE, invert = TRUE) , ]

### create cregion0??
toshiny.dengue.m1$CREGION0 <- toshiny.dengue.m1$CREGION

toshiny.dengue.m1.h <- subset(toshiny.dengue.m1, CREGION %in% c("IgM","IgG","IgA"))

toshiny.dengue.m1.h$CREGION <- gsub("IgA","IgH",toshiny.dengue.m1.h$CREGION)
toshiny.dengue.m1.h$CREGION <- gsub("IgG","IgH",toshiny.dengue.m1.h$CREGION)
toshiny.dengue.m1.h$CREGION <- gsub("IgM","IgH",toshiny.dengue.m1.h$CREGION)

ggplot(toshiny.dengue.m1.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue.m1.h, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengue.m1, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

#write.table(toshiny.dengue.m1, "toshiny_dengue_m1full.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.m1.h, "toshiny_dengue_m1.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#######################################
## next combine m1 & m2
toshiny.dengue.m1and2.h <- full_join(toshiny.dengue.m1.h, toshiny.dengue.m2.h)

## first need to recalculate ncount and shm.mean
toshiny.dengue.m1and2.h$ncount <- NULL
toshiny.dengue.m1and2.h$shm.mean <- NULL

toshiny.dengue.m1and2.h <- toshiny.dengue.m1and2.h %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))


toshiny.dengueallandmice <- bind_rows(toshiny.dengue.mabs.h, toshiny.dengue.oas.h, toshiny.dengue.d13.h, toshiny.dengue.d20.h, toshiny.dengue.m1and2.h, .id = "id")

toshiny.dengueallandmice$id <- gsub("2","DEN patients OAS data",toshiny.dengueallandmice$id)
toshiny.dengueallandmice$id <- gsub("3","DEN patient thirteen",toshiny.dengueallandmice$id)
toshiny.dengueallandmice$id <- gsub("4","DEN patient twenty",toshiny.dengueallandmice$id)
toshiny.dengueallandmice$id <- gsub("5","DEN vaccinated humanized mice",toshiny.dengueallandmice$id)
## because has 2 needs to be last
toshiny.dengueallandmice$id <- gsub("1","DEN plasmablasts",toshiny.dengueallandmice$id)

toshiny.dengueallandmice$id <- gsub("DEN patient thirteen","DEN patient 13",toshiny.dengueallandmice$id)
toshiny.dengueallandmice$id <- gsub("DEN patient twenty","DEN patient 20",toshiny.dengueallandmice$id)

toshiny.dengueallandmice$id <- factor(toshiny.dengueallandmice$id, levels = c("DEN plasmablasts", "DEN patient 13", "DEN patient 20", "DEN patients OAS data","DEN vaccinated humanized mice"))

## for combined need to recalculate ncount and shm.mean
toshiny.dengueallandmicec <- toshiny.dengueallandmice
toshiny.dengueallandmicec$ncount <- NULL
toshiny.dengueallandmicec$shm.mean <- NULL

toshiny.dengueallandmicec <- toshiny.dengueallandmicec %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.dengueallandmice, "toshiny_dengueallandmice.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengueallandmicec, "toshiny_dengueallandmicec.tab", sep = "\t", row.names = FALSE, quote = FALSE)
## note these are now 352k, was 370k, but much bigger filesize due to full sequences in there...

ggplot(toshiny.dengueallandmice, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.dengueallandmice, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## too big to keep around!!
rm(toshiny.dengueallandmice)
rm(toshiny.dengueallandmicec)

toshiny.dengue.d20
######################################################################################################

########################################################
### Binder data

### no longer needed - I split up outside R (mainly because for ireceptor data need to manually check & convert V score to V idenity)
# ## first need to split by unique values in column
# 
# ## need to clean up - remove * & _ in junctionaa, also any blanks in 
# BX.clusters.binderall <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/binder_data/airr-binder.tsv")
# 
# ## unique to this dataset
# BX.clusters.binderall <- BX.clusters.binderall %>% split(.$repertoire_id)
# 
# #toshiny.natalia4 <- toshiny.natalia4 %>% split(.$id)
# 
# split(toshiny.natalia4, toshiny.natalia4$id)
# 
# ### this works (splits into 4 separate dataframes)
# list2env(x = split(x = toshiny.natalia4,
#                    f = toshiny.natalia4$id),
#          envir = globalenv())
# 
# 
# list2env(x = split(x = BX.clusters.binderall,
#                    f = BX.clusters.binderall$repertoire_id),
#          envir = globalenv())
# 

## binder data for patients 2,3,6-14,16 (12 total)


## this removes first and last aa from junction (because ireceptor includes these) BE CAREFUL NOT TO EVER RUN THIS ON LC DATA...
# str_sub(BX.clusters.hiv$JUNCTIONAA, 1, 1) <- ""
# str_sub(BX.clusters.hiv$JUNCTIONAA, -1, -1) <- ""



################################

## cleaned-up commands for single dataset input

### MERS/SARS mab list

BX.clusters.sarsmers <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/SARSMERS_mablist_clonalclusters.tab")
BX.clusters.sarsmers$JUNCTIONAA <- as.character(BX.clusters.sarsmers$JUNCTIONAA)


BX.clusters.sarsmers$CDR3LENGTH_IMGT <- nchar(BX.clusters.sarsmers$JUNCTIONAA)
BX.clusters.sarsmers$SHM <- (100 - (BX.clusters.sarsmers$V_IDENTITY * 100))
#BX.clusters.sarsmers$SHM <- 0

BX.clusters.sarsmers <- BX.clusters.sarsmers %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
#BX.allclusters.b1 <- BX.allclusters.b1 %>% filter(is.wholenumber(CDR3LENGTH_IMGT))


## specific to dataset
#BX.clusters.mabs$reads_per_clone <- 1
BX.clusters.sarsmers$neutralization <- gsub("weakneutralizing","neutralizing",BX.clusters.sarsmers$neutralization)
BX.clusters.sarsmers$neutralization <- gsub("neutralizingbroad","neutralizing",BX.clusters.sarsmers$neutralization)
BX.clusters.sarsmers$CREGION <- gsub("IGK","Kappa",BX.clusters.sarsmers$CREGION)
BX.clusters.sarsmers$CREGION <- gsub("IGL","Lambda",BX.clusters.sarsmers$CREGION)
BX.clusters.sarsmers$CREGION <- gsub("IGH","IgH",BX.clusters.sarsmers$CREGION)


BX.clusters.sarsmers.h <- subset(BX.clusters.sarsmers, CREGION %in% c("IgH"))
# BX.clusters.mabs.h <- BX.clusters.mabs.h %>%
#   group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   mutate(shm.mean = mean(SHM, na.rm = TRUE))

toshiny.sarsmers <- BX.clusters.sarsmers %>% select(SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.sarsmers.h <- BX.clusters.sarsmers.h %>% select(SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
#write.table(toshiny.sarsmers, "toshiny_sarsmers.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.sarsmers.h, "toshiny_sarsmers_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## splitting up SARS and MERS datasets before merging with CoV2 mabs...

list2env(x = split(x = BX.clusters.sarsmers,
                   f = BX.clusters.sarsmers$source),
         envir = globalenv())

BX.clusters.sars <- SARS
BX.clusters.mers <- MERS
rm(SARS)
rm(MERS)

list2env(x = split(x = BX.clusters.sarsmers.h,
                   f = BX.clusters.sarsmers.h$source),
         envir = globalenv())

BX.clusters.sars.h <- SARS
BX.clusters.mers.h <- MERS
rm(SARS)
rm(MERS)

####################################
### combining CoV2, MERS & SARS datasets for vis tool
toshiny.sars.h <- BX.clusters.sars.h %>% select(SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.mers.h <- BX.clusters.mers.h %>% select(SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))


toshiny.cov2sarsmers <- bind_rows(toshiny.few1.h, toshiny.sars.h, toshiny.mers.h, .id = "id")

toshiny.cov2sarsmers$id <- gsub("2","anti-SARS mAbs",toshiny.cov2sarsmers$id)
toshiny.cov2sarsmers$id <- gsub("3","anti-MERS mAbs",toshiny.cov2sarsmers$id)
## because has 2 needs to be last
toshiny.cov2sarsmers$id <- gsub("1","anti-CoV2 mAbs",toshiny.cov2sarsmers$id)

toshiny.cov2sarsmers$id <- factor(toshiny.cov2sarsmers$id, levels = c("anti-CoV2 mAbs", "anti-SARS mAbs", "anti-MERS mAbs"))

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


####################################
## new combination because IAVI84 is only HV1 & HV4 - N152 has HV3 & HV4 (so still no HV2 or HV5-7)
BX.clusters.hiv <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/hiv_iavi84and152hc.tab")

## this removes first and last aa from junction (because ireceptor includes these) BE CAREFUL NOT TO EVER RUN THIS ON LC DATA...
# str_sub(BX.clusters.hiv$JUNCTIONAA, 1, 1) <- ""
# str_sub(BX.clusters.hiv$JUNCTIONAA, -1, -1) <- ""

BX.clusters.hiv$GENE <- getGene(BX.clusters.hiv$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$GF <- substring(BX.clusters.hiv$GENE, 1,5)
BX.clusters.hiv$JGENE <- getGene(BX.clusters.hiv$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$JGF <- substring(BX.clusters.hiv$JGENE, 1,5)
BX.clusters.hiv$VJ_JUNCTION_PATTERN <- paste(BX.clusters.hiv$GF,BX.clusters.hiv$JGF,BX.clusters.hiv$JUNCTION_LENGTH,sep="_") 

BX.clusters.hiv$CDR3LENGTH_IMGT <- ((BX.clusters.hiv$JUNCTION_LENGTH) / 3) - 2
#BX.clusters.hiv$JGENE[BX.clusters.hiv$JGENE==""]<-NA
BX.clusters.hiv$SHM <- (100 - (BX.clusters.hiv$V_IDENTITY * 100))
BX.clusters.hiv <- BX.clusters.hiv %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

BX.clusters.hiv <- BX.clusters.hiv %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
BX.clusters.hiv$CREGION <- "IgH"
toshiny.fewhiv <- BX.clusters.hiv %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,VJ_JUNCTION_PATTERN,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))



##############################################################################################################################
### for full set of reads not clones:

BX.allclusters.b1 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/Boyd_results/Boyd_7450ab10_germ-pass.tab")
BX.allclusters.g1 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/Galson1_germ-pass.tab")

BX.allclusters.g1$JUNCTIONAA <- translateDNA(BX.allclusters.g1$JUNCTION, trim=TRUE)
BX.allclusters.g1$GENE <- getGene(BX.allclusters.g1$V_CALL, first=TRUE, strip_d=TRUE)
BX.allclusters.g1$GF <- substring(BX.allclusters.g1$GENE, 1,5)
BX.allclusters.g1$JGENE <- getGene(BX.allclusters.g1$J_CALL, first=TRUE, strip_d=TRUE)
BX.allclusters.g1$JGF <- substring(BX.allclusters.g1$JGENE, 1,5)
BX.allclusters.g1$VJ_JUNCTION_PATTERN <- paste(BX.allclusters.g1$GF,BX.allclusters.g1$JGF,BX.allclusters.g1$JUNCTION_LENGTH,sep="_") 

BX.allclusters.g1$CDR3LENGTH_IMGT <- ((BX.allclusters.g1$JUNCTION_LENGTH) / 3) - 2
BX.allclusters.g1$SHM <- (100 - (BX.allclusters.g1$V_IDENTITY * 100))
BX.allclusters.g1 <- BX.allclusters.g1 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
BX.allclusters.g1 <- BX.allclusters.g1 %>% filter(is.wholenumber(CDR3LENGTH_IMGT))


## specific to dataset
BX.allclusters.g1$CREGION <- gsub("IgA","IgH",BX.allclusters.g1$CREGION)
BX.allclusters.g1$CREGION <- gsub("IgD","IgH",BX.allclusters.g1$CREGION)
BX.allclusters.g1$CREGION <- gsub("IgM","IgH",BX.allclusters.g1$CREGION)
BX.allclusters.g1$CREGION <- gsub("IgG","IgH",BX.allclusters.g1$CREGION)


#####
BX.allclusters.b1$JUNCTIONAA <- translateDNA(BX.allclusters.b1$JUNCTION, trim=TRUE)
BX.allclusters.b1$GENE <- getGene(BX.allclusters.b1$V_CALL, first=TRUE, strip_d=TRUE)
BX.allclusters.b1$GF <- substring(BX.allclusters.b1$GENE, 1,5)
BX.allclusters.b1$JGENE <- getGene(BX.allclusters.b1$J_CALL, first=TRUE, strip_d=TRUE)
BX.allclusters.b1$JGF <- substring(BX.allclusters.b1$JGENE, 1,5)
BX.allclusters.b1$VJ_JUNCTION_PATTERN <- paste(BX.allclusters.b1$GF,BX.allclusters.b1$JGF,BX.allclusters.b1$JUNCTION_LENGTH,sep="_") 

BX.allclusters.b1$CDR3LENGTH_IMGT <- ((BX.allclusters.b1$JUNCTION_LENGTH) / 3) - 2
BX.allclusters.b1$SHM <- (100 - (BX.allclusters.b1$V_IDENTITY * 100))
BX.allclusters.b1 <- BX.allclusters.b1 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))
BX.allclusters.b1 <- BX.allclusters.b1 %>% filter(is.wholenumber(CDR3LENGTH_IMGT))


## specific to dataset
BX.allclusters.b1$CREGION <- gsub("IgA","IgH",BX.allclusters.b1$CREGION)
BX.allclusters.b1$CREGION <- gsub("IgD","IgH",BX.allclusters.b1$CREGION)
BX.allclusters.b1$CREGION <- gsub("IgM","IgH",BX.allclusters.b1$CREGION)
BX.allclusters.b1$CREGION <- gsub("IgG","IgH",BX.allclusters.b1$CREGION)


toshiny.fewg1b <- BX.clusters.g1 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.fewb1b <- BX.clusters.1 %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))


##################################################
### FULL SET OF COMMANDS IF NOT CLUSTERED USING PREVIOUS COMMANDS - USING FOR FULL HEALTHY DATASET (NOT COLLAPSED BY CLONES)

#healthy <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/BXmay10m_stim_germ-pass.tab")
BX.allclusters.hc <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/BXmay10m_stim_germ-pass.tab")

BX.allclusters.hc$V_CALL <- as.character(BX.allclusters.hc$V_CALL)
BX.allclusters.hc$J_CALL <- as.character(BX.allclusters.hc$J_CALL)
BX.allclusters.hc$V_CALL <- as.character(BX.allclusters.hc$V_CALL)

#BX.full0 <- BX.full0 %>% add_column(CONSCOUNT = 1, DUPCOUNT = 1)
BX.allclusters.hc <- BX.allclusters.hc %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)

# BX.clones.hc <- collapseClones(BX.full0, method="mostCommon", cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT", germlineColumn = "GERMLINE_IMGT_D_MASK",
#                                includeAmbiguous=FALSE, breakTiesStochastic=FALSE)
# BX.clusters.hc <- groupGenes(BX.clones.hc, v_call = "V_CALL", j_call = "J_CALL",
#                               junc_len = "JUNCTION_LENGTH", cell_id = NULL, locus = NULL, only_igh = TRUE,
#                               first = TRUE)

BX.allclusters.hc$JUNCTIONAA <- translateDNA(BX.allclusters.hc$JUNCTION, trim=TRUE)
BX.allclusters.hc$GENE <- getGene(BX.allclusters.hc$V_CALL, first=TRUE, strip_d=TRUE)
BX.allclusters.hc$GF <- substring(BX.allclusters.hc$GENE, 1,5)
BX.allclusters.hc$JGENE <- getGene(BX.allclusters.hc$J_CALL, first=TRUE, strip_d=TRUE)
BX.allclusters.hc$JGF <- substring(BX.allclusters.hc$JGENE, 1,5)
BX.allclusters.hc$VJ_JUNCTION_PATTERN <- paste(BX.allclusters.hc$GF,BX.allclusters.hc$JGF,BX.allclusters.hc$JUNCTION_LENGTH,sep="_") 

#BX.allclusters.hc$CDR3KABAT <- ((BX.allclusters.hc$JUNCTION_LENGTH) / 3)
BX.allclusters.hc$CDR3LENGTH_IMGT <- ((BX.allclusters.hc$JUNCTION_LENGTH) / 3) - 2
BX.allclusters.hc$SHM <- (100 - (BX.allclusters.hc$V_IDENTITY * 100))
BX.allclusters.hc <- BX.allclusters.hc %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

BX.allclusters.hc <- BX.allclusters.hc %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
## specific to dataset
BX.allclusters.hc$CREGION <- gsub("IgA","IgH",BX.allclusters.hc$CREGION)
BX.allclusters.hc$CREGION <- gsub("IgD","IgH",BX.allclusters.hc$CREGION)
BX.allclusters.hc$CREGION <- gsub("IgM","IgH",BX.allclusters.hc$CREGION)
BX.allclusters.hc$CREGION <- gsub("IgG","IgH",BX.allclusters.hc$CREGION)


## then subsets


####################################
## healthy controls

BX.clusters.hc <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/Galson_healthycontrol_clonalclusters.tab")
#BX.clusters.hc$CDR3KABAT <- ((BX.clusters.hc$JUNCTION_LENGTH) / 3)
BX.clusters.hc$CDR3LENGTH_IMGT <- ((BX.clusters.hc$JUNCTION_LENGTH) / 3) - 2
BX.clusters.hc$SHM <- (100 - (BX.clusters.hc$V_IDENTITY * 100))
BX.clusters.hc <- BX.clusters.hc %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))


## specific to dataset
BX.clusters.hc$CREGION <- gsub("IgA","IgH",BX.clusters.hc$CREGION)
BX.clusters.hc$CREGION <- gsub("IgD","IgH",BX.clusters.hc$CREGION)
BX.clusters.hc$CREGION <- gsub("IgM","IgH",BX.clusters.hc$CREGION)
BX.clusters.hc$CREGION <- gsub("IgG","IgH",BX.clusters.hc$CREGION)
## because of misassigned genes, removing
#BX.clusters.hc$VJ_JUNCTION_PATTERN <- gsub("_IGKJ1_","_",BX.clusters.hc$VJ_JUNCTION_PATTERN)
#BX.clusters.hc$VJ_JUNCTION_PATTERN <- gsub("_IGLJ7_","_",BX.clusters.hc$VJ_JUNCTION_PATTERN)



####
## repeat for hc2 and hc3
#BX.clusters.hc2$CDR3KABAT <- ((BX.clusters.hc2$JUNCTION_LENGTH) / 3)
BX.clusters.hc2$CDR3LENGTH_IMGT <- ((BX.clusters.hc2$JUNCTION_LENGTH) / 3) - 2
BX.clusters.hc2$CREGION <- gsub("IgA","IgH",BX.clusters.hc2$CREGION)
BX.clusters.hc2$CREGION <- gsub("IgD","IgH",BX.clusters.hc2$CREGION)
BX.clusters.hc2$CREGION <- gsub("IgM","IgH",BX.clusters.hc2$CREGION)
BX.clusters.hc2$CREGION <- gsub("IgG","IgH",BX.clusters.hc2$CREGION)
#BX.clusters.hc2$JGENE[BX.clusters.hc2$JGENE==""]<-NA
BX.clusters.hc2 <- BX.clusters.hc2 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

#BX.clusters.hc3$CDR3KABAT <- ((BX.clusters.hc3$JUNCTION_LENGTH) / 3)
BX.clusters.hc3$CDR3LENGTH_IMGT <- ((BX.clusters.hc3$JUNCTION_LENGTH) / 3) - 2
BX.clusters.hc3$CREGION <- gsub("IgA","IgH",BX.clusters.hc3$CREGION)
BX.clusters.hc3$CREGION <- gsub("IgD","IgH",BX.clusters.hc3$CREGION)
BX.clusters.hc3$CREGION <- gsub("IgM","IgH",BX.clusters.hc3$CREGION)
BX.clusters.hc3$CREGION <- gsub("IgG","IgH",BX.clusters.hc3$CREGION)
#BX.clusters.hc3$JGENE[BX.clusters.hc3$JGENE==""]<-NA
BX.clusters.hc3 <- BX.clusters.hc3 %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

#BX.clusters.hc2$VJ_JUNCTION_PATTERN <- gsub("_IGLJ7_","_",BX.clusters.hc2$VJ_JUNCTION_PATTERN)
#BX.H <- BX.H[ grep("IGKV", BX.H$V_CALL, invert = TRUE) , ]

BX.clusters.hc2r <- filter(BX.clusters.hc2, CREGION == "Lambda" & GF == "IGHV3")
BX.clusters.hc2 <- anti_join(BX.clusters.hc2, BX.clusters.hc2r)

BX.clusters.hc3r <- filter(BX.clusters.hc3, JGENE == "IGLJ1" & GF == "IGHV3")
BX.clusters.hc3 <- anti_join(BX.clusters.hc3, BX.clusters.hc3r)

####################################
### comet mabs

biohub.all2$CREGION <- "IgH"
biohub.all2$GF <- substring(biohub.all2$v_gene, 1,5)
#biohub.all2$JUNCTION_LENGTH <- nchar(biohub.all2$JUNCTIONAA)
#biohub.all2$CDR3KABAT <- ((biohub.all2$JUNCTION_LENGTH) / 3)
#biohub.all2$CDR3LENGTH_IMGT <- ((biohub.all2$JUNCTION_LENGTH) / 3) - 2
biohub.all2$JUNCTION_LENGTH <- NULL
biohub.all2$CDR3LENGTH_IMGT <- nchar(biohub.all2$JUNCTIONAA)
biohub.all2 <- biohub.all2 %>% unite(GF_JGENE, GF, j_gene, sep = "_", remove = FALSE, na.rm = TRUE)


#test <- 1850 %/% 500
#biohub.all2$c_umi <- NULL
biohub.all2$reads_per_clone <- biohub.all2$n_umi %/% 100

comet.pos <- subset(biohub.all2, COVID_status %in% c("Positive"))
comet.neg <- subset(biohub.all2, COVID_status %in% c("Negative"))
#comet.pos <- filter(BX.clusters.hc3, JGENE == "IGLJ1" & GF == "IGHV3")

## note comet data has no shm...



### healthy controls 2 & 3

hc2 <- ggplot(BX.clusters.hc2, aes(GF_JGENE,CDR3LENGTH_IMGT))
hc2 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

hc3 <- ggplot(BX.clusters.hc3, aes(GF_JGENE,CDR3LENGTH_IMGT))
hc3 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(BX.clusters.hc3, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


## combining all 3 healthy controls
BX.clusters.hc$CLONE <- as.character(BX.clusters.hc$CLONE)

BX.clusters.allhc <- bind_rows(BX.clusters.hc, BX.clusters.hc2, BX.clusters.hc3, .id = "id")

BX.clusters.allhc$id <- gsub("3","Healthy control 3",BX.clusters.allhc$id)
BX.clusters.allhc$id <- gsub("2","Healthy control 2",BX.clusters.allhc$id)
BX.clusters.allhc$id <- gsub("1","Healthy control 1",BX.clusters.allhc$id)

BX.clusters.allhc.h <- subset(BX.clusters.allhc, CREGION %in% c("IgH"))
BX.clusters.allhc.k <- subset(BX.clusters.allhc, CREGION %in% c("Kappa"))
BX.clusters.allhc.l <- subset(BX.clusters.allhc, CREGION %in% c("Lambda"))

h3.allhc.h <- ggplot(BX.clusters.allhc.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.allhc.h + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.allhc.gfbyjgenebycdr3.viridis.h <-  h3.allhc.h + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.allhc.gfbyjgenebycdr3_h.pdf", BX.clusters.allhc.gfbyjgenebycdr3.viridis.h, width = 16, height = 12, units = "in")

h3.allhc.k <- ggplot(BX.clusters.allhc.k, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.allhc.k + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.allhc.gfbyjgenebycdr3.viridis.k <-  h3.allhc.k + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.allhc.gfbyjgenebycdr3_k.pdf", BX.clusters.allhc.gfbyjgenebycdr3.viridis.k, width = 16, height = 12, units = "in")

h3.allhc.l <- ggplot(BX.clusters.allhc.l, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.allhc.l + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.allhc.gfbyjgenebycdr3.viridis.l <-  h3.allhc.l + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.allhc.gfbyjgenebycdr3_l.pdf", BX.clusters.allhc.gfbyjgenebycdr3.viridis.l, width = 16, height = 12, units = "in")

## to uncount (will be huge)
#BX.clusters.allhc.h2 <- uncount(BX.clusters.allhc.h, BX.clusters.mabsvshc.h$reads_per_clone)




############################################################################################################
### combining datasets


#############
### combining mabs & HC

BX.clusters.mabsvshc <- bind_rows(BX.clusters.mabs, BX.clusters.hc, .id = "id")

BX.clusters.mabsvshc$id <- gsub("2","Healthy control bulk repertoire",BX.clusters.mabsvshc$id)
BX.clusters.mabsvshc$id <- gsub("1","anti-CoV2 mAbs",BX.clusters.mabsvshc$id)

h3.mabsvshc <- ggplot(BX.clusters.mabsvshc, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.mabsvshc + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION + id, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nClones",  option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


## but need to figure out how to plot by id as well...
#h3.mabsvshc <- ggplot(BX.clusters.mabsvshc, aes(GF_JGENE,CDR3LENGTH_IMGT))
#h3.mabsvshc + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_grid(rows = vars(CREGION), cols = vars(id), scales = "free") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

BX.clusters.mabsvshc.h <- subset(BX.clusters.mabsvshc, CREGION %in% c("IgH"))
BX.clusters.mabsvshc.k <- subset(BX.clusters.mabsvshc, CREGION %in% c("Kappa"))
BX.clusters.mabsvshc.l <- subset(BX.clusters.mabsvshc, CREGION %in% c("Lambda"))

h3.mabsvshc.h <- ggplot(BX.clusters.mabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
## counts
h3.mabsvshc.h + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## percent of each group
h3.mabsvshc.h + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nClones",  option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#test - we want to use this to get actual percentages
h3.mabsvshc.h + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nClones",  option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


## for shm
h3shm.mabsvshc.h <- ggplot(BX.clusters.mabsvshc.h, aes(GF_JGENE,SHM))
## percent of each group
h3shm.mabsvshc.h + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + theme_bw() + ylab("Somatic Hypermutation (%)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nClones",  option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# raw counts
h3shm.mabsvshc.h + geom_bin2d(aes(fill= log10(..count..))) + theme_bw() + ylab("Somatic Hypermutation (%)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))



BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.h <-  h3.mabsvshc.h + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabsvshc.gfbyjgenebycdr3_h.pdf", BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.h, width = 16, height = 12, units = "in")

h3.mabsvshc.k <- ggplot(BX.clusters.mabsvshc.k, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.mabsvshc.k + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.k <-  h3.mabsvshc.k + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabsvshc.gfbyjgenebycdr3_k.pdf", BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.k, width = 16, height = 12, units = "in")

h3.mabsvshc.l <- ggplot(BX.clusters.mabsvshc.l, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.mabsvshc.l + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.l <-  h3.mabsvshc.l + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nClones", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabsvshc.gfbyjgenebycdr3_l.pdf", BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.l, width = 16, height = 12, units = "in")

## to cluster by reads not by clonal families...but all mabs are 1 or 100, not the sum
## first uncount by reads_per_clone
BX.clusters.mabsvshc.h2 <- uncount(BX.clusters.mabsvshc.h, BX.clusters.mabsvshc.h$reads_per_clone)


h3.mabsvshc.h2 <- ggplot(BX.clusters.mabsvshc.h2, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.mabsvshc.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.h2 <-  h3.mabsvshc.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabsvshc.gfbyjgenebycdr3_h2.pdf", BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.h2, width = 16, height = 12, units = "in")

BX.clusters.mabsvshc.k2 <- uncount(BX.clusters.mabsvshc.k, BX.clusters.mabsvshc.k$reads_per_clone)
h3.mabsvshc.k2 <- ggplot(BX.clusters.mabsvshc.k2, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.mabsvshc.k2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous() + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.k2 <-  h3.mabsvshc.k2 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabsvshc.gfbyjgenebycdr3_k2.pdf", BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.k2, width = 16, height = 12, units = "in")

BX.clusters.mabsvshc.l2 <- uncount(BX.clusters.mabsvshc.l, BX.clusters.mabsvshc.l$reads_per_clone)
h3.mabsvshc.l2 <- ggplot(BX.clusters.mabsvshc.l2, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.mabsvshc.l2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous() + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.l2 <-  h3.mabsvshc.l2 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRL3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabsvshc.gfbyjgenebycdr3_l2.pdf", BX.clusters.mabsvshc.gfbyjgenebycdr3.viridis.l2, width = 16, height = 12, units = "in")


####################################
### HC vs. COMET vs. mab list

BX.clusters.cometvsmabsvshc <- bind_rows(BX.clusters.hc, comet.pos, comet.neg, BX.clusters.mabs, .id = "id")
BX.clusters.cometvsmabsvshc$id <- gsub("1","Healthy control bulk repertoire",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("2","COMET COVID-positive",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("3","COMET COVID-negative",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("4","anti-CoV2 mAbs",BX.clusters.cometvsmabsvshc$id)

BX.clusters.cometvsmabsvshc.h <- subset(BX.clusters.cometvsmabsvshc, CREGION %in% c("IgH"))

BX.clusters.cometvsmabsvshc.h$id <- factor(BX.clusters.cometvsmabsvshc.h$id, levels = c("Healthy control bulk repertoire", "COMET COVID-positive", "COMET COVID-negative", "anti-CoV2 mAbs"))

## uncount both HC & COMET by reads_per_clone
BX.clusters.cometvsmabsvshc.h2 <- uncount(BX.clusters.cometvsmabsvshc.h, BX.clusters.cometvsmabsvshc.h$reads_per_clone)

## full set of HC reads not clones & UMIs
h3.cometvsmabsvshc.h22 <- ggplot(BX.clusters.cometvsmabsvshc.h2, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.cometvsmabsvshc.h22 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
h3.cometvsmabsvshc.h22 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3.viridis.h22 <- h3.cometvsmabsvshc.h22 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_h_fullUMIs_bypercentage.pdf", BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3.viridis.h22, width = 16, height = 12, units = "in")


####################################
## Natalia random 500 cells
random500m200 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/Natalia_sequencing/random_500_germ-pass_min200.tab")

BX.clusters.nat <- groupGenes(random500m200, v_call = "V_CALL", j_call = "J_CALL",
                              junc_len = "JUNCTION_LENGTH", cell_id = NULL, locus = NULL, only_igh = TRUE,
                              first = TRUE)
BX.clusters.nat$JUNCTIONAA <- translateDNA(BX.clusters.nat$JUNCTION, trim=TRUE)
BX.clusters.nat$GENE <- getGene(BX.clusters.nat$V_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.nat$GF <- substring(BX.clusters.nat$GENE, 1,5)
BX.clusters.nat$JGENE <- getGene(BX.clusters.nat$J_CALL, first=TRUE, strip_d=TRUE)
BX.clusters.nat$JGF <- substring(BX.clusters.nat$JGENE, 1,5)
BX.clusters.nat$VJ_JUNCTION_PATTERN <- paste(BX.clusters.nat$GF,BX.clusters.nat$JGF,BX.clusters.nat$JUNCTION_LENGTH,sep="_") 

#BX.clusters.nat$CDR3KABAT <- ((BX.clusters.nat$JUNCTION_LENGTH) / 3)
BX.clusters.nat$CDR3LENGTH_IMGT <- ((BX.clusters.nat$JUNCTION_LENGTH) / 3) - 2
BX.clusters.nat <- BX.clusters.nat %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

h3.natalia <- ggplot(BX.clusters.nat, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.natalia + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

### NATALIA + ANTI-COV2 + HC
BX.clusters.hc$CLONE <- as.integer(BX.clusters.hc$CLONE)
BX.clusters.hc.h <- subset(BX.clusters.hc, CREGION %in% c("IgH"))
BX.clusters.mabs.h <- subset(BX.clusters.mabs, CREGION %in% c("IgH"))


BX.clusters.mabsvshcvsnat <- bind_rows(BX.clusters.mabs.h, BX.clusters.hc.h, BX.clusters.nat, .id = "id")
BX.clusters.mabsvshcvsnat$id <- gsub("3","Natalia 500 randomly sorted cells",BX.clusters.mabsvshcvsnat$id)
BX.clusters.mabsvshcvsnat$id <- gsub("2","Healthy control bulk repertoire",BX.clusters.mabsvshcvsnat$id)
BX.clusters.mabsvshcvsnat$id <- gsub("1","anti-CoV2 mAbs",BX.clusters.mabsvshcvsnat$id)

h3.cometvsmabsvsnat.h2 <- ggplot(BX.clusters.mabsvshcvsnat, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.cometvsmabsvsnat.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
h3.cometvsmabsvsnat.h2 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h <-  h3.cometvsmabsvsnat.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("mabsvshcvsnatalia.gfbyjgenebycdr3_h.pdf", cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h, width = 16, height = 12, units = "in")
cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h2 <-  h3.cometvsmabsvsnat.h2 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("mabsvshcvsnatalia.gfbyjgenebycdr3_h2.pdf", cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h2, width = 16, height = 12, units = "in")


############ latest commands used

### ALSO HAVE TO ADD NCOUNT FOR EACH DATASET - thus really hard to automate
## BUT HAVE TO DO IT BEFORE COMBINING DATASETS...
# BX.clusters.cometvsmabsvshc.h.few <- BX.clusters.cometvsmabsvshc.h.few %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   rename(ncount = n)

BX.clusters.hc <- BX.clusters.hc %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
comet.pos <- comet.pos %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
comet.neg <- comet.neg %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
BX.clusters.mabs <- BX.clusters.mabs %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)

BX.clusters.cometvsmabsvshc <- bind_rows(BX.clusters.hc, comet.pos, comet.neg, BX.clusters.mabs, .id = "id")
BX.clusters.cometvsmabsvshc$id <- gsub("1","Healthy control bulk repertoire",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("2","COMET COVID-positive",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("3","COMET COVID-negative",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("4","anti-CoV2 mAbs",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- factor(BX.clusters.cometvsmabsvshc$id, levels = c("Healthy control bulk repertoire", "COMET COVID-positive", "COMET COVID-negative", "anti-CoV2 mAbs"))

BX.clusters.cometvsmabsvshc.h <- subset(BX.clusters.cometvsmabsvshc, CREGION %in% c("IgH"))

## if combining any dataset need to relabel id's, also factor (so might be hard to automate this)
## also would have to ensure CREGION names are consistent...

BX.clusters.mabsvshc <- bind_rows(BX.clusters.mabs, BX.clusters.hc, .id = "id")

BX.clusters.mabsvshc$id <- gsub("2","Healthy control bulk repertoire",BX.clusters.mabsvshc$id)
BX.clusters.mabsvshc$id <- gsub("1","anti-CoV2 mAbs",BX.clusters.mabsvshc$id)

BX.clusters.mabsvshc$id <- factor(BX.clusters.mabsvshc$id, levels = c("Healthy control bulk repertoire", "anti-CoV2 mAbs"))

BX.clusters.mabsvshc.h <- subset(BX.clusters.mabsvshc, CREGION %in% c("IgH"))
#BX.clusters.mabsvshc.k <- subset(BX.clusters.mabsvshc, CREGION %in% c("Kappa"))
#BX.clusters.mabsvshc.l <- subset(BX.clusters.mabsvshc, CREGION %in% c("Lambda"))

## only make .few after converting to toshiny.few


#################################
#################################
##### AUGUST 2020 CODE...

## healthy control trying by SHM
BX.clusters.hc$SHM <- (100 - (BX.clusters.hc$V_IDENTITY * 100))
hc.byshm <- ggplot(BX.clusters.hc, aes(GF_JGENE,SHM))
hc.byshm + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("Somatic Hypermutation (%)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
hc.byshm + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + theme_bw() + ylab("Somatic Hypermutation (%)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

BX.clusters.mabs$SHM <- (100 - (BX.clusters.mabs$V_IDENTITY * 100))
mabs.byshm <- ggplot(BX.clusters.mabs, aes(GF_JGENE,SHM))
mabs.byshm + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("Somatic Hypermutation (%)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
mabs.byshm + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + theme_bw() + ylab("Somatic Hypermutation (%)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## mabs by binding & by neutralization (first subset to HC)

BX.clusters.mabs.h <- subset(BX.clusters.mabs, CREGION %in% c("IgH"))
mabs.h <- ggplot(BX.clusters.mabs.h, aes(GF_JGENE,CDR3LENGTH_IMGT))


## by neutralizing vs. not
mabs.h + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ neutralization, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## by binding
mabs.h + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ binding, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


ggplot(BX.clusters.mabs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


## trying plot where SHM is color not counts...
#h3.mabsvshc.h <- ggplot(BX.clusters.mabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
##this doesn't work because of fractions of CDR3 lengths - after the extra command to get to toshiny.few it does work though...
#h3.mabsvshc.h + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Somatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## THIS WORKS - COULD BE USED AS ANOTHER SHINY PLOT
ggplot(toshiny.few, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## need to update comet data

BX.clusters.comet <- bind_rows(comet.pos, comet.neg, .id = "id")
BX.clusters.comet$id <- gsub("1","COMET COVID-positive",BX.clusters.comet$id)
BX.clusters.comet$id <- gsub("2","COMET COVID-negative",BX.clusters.comet$id)

### code to replace dummy shiny dataset name with particular dataset

## need to remove Augmenta for online version...also maybe not have Comet data up?

##### THEN REPLACE BELOW ↓↓↓↓
toshiny.few <- BX.clusters.mabsvshc.h %>% select(id,SEQUENCE_ID,binding,neutralization,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.few$initial <- NULL
toshiny.few$final <- NULL
toshiny.few$JUNCTION_LENGTH0 <- NULL

#switch(input$dataset, "COVID mAbs" = BX.clusters.mabs, "Comet Data" = biohub.all2, "Healthy Control" = BX.clusters.hc, "Healthy Control vs. Comet Data vs. COVID mAbs" = toshiny.few)

toshiny.few1 <- BX.clusters.mabs %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.few1$initial <- NULL
toshiny.few1$final <- NULL
toshiny.few1$JUNCTION_LENGTH0 <- NULL

toshiny.few1.h <- BX.clusters.mabs.h %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.few1.h$initial <- NULL
toshiny.few1.h$final <- NULL
toshiny.few1.h$JUNCTION_LENGTH0 <- NULL


## no shm so slightly different
toshiny.few2 <- BX.clusters.comet %>% select(id,clonotype_id,CREGION,chain:c_gene,JUNCTIONAA:ncount) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
toshiny.few3 <- BX.clusters.hc %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))

#toshiny.few <- BX.clusters.mabsvshc.k %>% select(id,SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
#toshiny.few <- BX.clusters.mabsvshc.l %>% select(id,SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
## full HC set
toshiny.few3b <- BX.allclusters.hc %>% select(SEQUENCE_ID,CREGION,JUNCTIONAA:VJ_JUNCTION_PATTERN,SHM,ncount,shm.mean,CDR3LENGTH_IMGT,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.few3b$SEQUENCE_ID <- as.character(toshiny.few3b$SEQUENCE_ID)
toshiny.few3b$JUNCTIONAA <- as.character(toshiny.few3b$JUNCTIONAA)
# write.table(toshiny.few3b, "toshiny_few3b.tab", sep = "\t", row.names = FALSE, quote = FALSE)


####################################
### removing Augmenta for 2 of these files (.fewer)
BX.clusters.mabsnoaugmenta <- filter(BX.clusters.mabs, source != "Augmenta")
BX.clusters.mabsnoaugmenta.h <- subset(BX.clusters.mabsnoaugmenta, CREGION %in% c("IgH"))

toshiny.fewer1 <- BX.clusters.mabsnoaugmenta %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.fewer1$initial <- NULL
toshiny.fewer1$final <- NULL
toshiny.fewer1$JUNCTION_LENGTH0 <- NULL

toshiny.fewer1.h <- BX.clusters.mabsnoaugmenta.h %>% select(SEQUENCE_ID,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.fewer1.h$initial <- NULL
toshiny.fewer1.h$final <- NULL
toshiny.fewer1.h$JUNCTION_LENGTH0 <- NULL

BX.clusters.mabsvshcnoaugmenta <- bind_rows(BX.clusters.mabsnoaugmenta, BX.clusters.hc, .id = "id")
BX.clusters.mabsvshcnoaugmenta$id <- gsub("2","Healthy control bulk repertoire",BX.clusters.mabsvshcnoaugmenta$id)
BX.clusters.mabsvshcnoaugmenta$id <- gsub("1","anti-CoV2 mAbs",BX.clusters.mabsvshcnoaugmenta$id)
BX.clusters.mabsvshcnoaugmenta$id <- factor(BX.clusters.mabsvshcnoaugmenta$id, levels = c("Healthy control bulk repertoire", "anti-CoV2 mAbs"))
BX.clusters.mabsvshcnoaugmenta.h <- subset(BX.clusters.mabsvshcnoaugmenta, CREGION %in% c("IgH"))
toshiny.fewer <- BX.clusters.mabsvshcnoaugmenta.h %>% select(id,SEQUENCE_ID,binding,neutralization,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
toshiny.fewer$initial <- NULL
toshiny.fewer$final <- NULL
toshiny.fewer$JUNCTION_LENGTH0 <- NULL


#write.table(toshiny.fewer1, "toshiny_fewer1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewer, "toshiny_fewer.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## to round CDR3 values - note added to above command
### new function
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
toshiny.few <- toshiny.few %>% filter(is.wholenumber(CDR3LENGTH_IMGT))

## plots used in Shiny app (looking at full healthy control dataset by CREGION here)
ggplot(BX.allclusters.hc, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_tile(aes(fill = shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol = 1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
ggplot(BX.allclusters.hc, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol = 1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
ggplot(BX.clusters.hc, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 12) + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol = 1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=8))

####################################
### to compare Augmenta vs. non-Augmenta


BX.clusters.mabsnoaugmenta <- filter(BX.clusters.mabs, source != "Augmenta")
BX.clusters.mabsnoaugmenta.h <- subset(BX.clusters.mabsnoaugmenta, CREGION %in% c("IgH"))

## need to reset shm.mean and ncount
BX.clusters.mabsnoaugmenta$ncount <- NULL
BX.clusters.mabsnoaugmenta$shm.mean <- NULL

BX.clusters.mabsnoaugmenta <- BX.clusters.mabsnoaugmenta %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

toshiny.fewer1 <- BX.clusters.mabsnoaugmenta %>% select(SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
# toshiny.fewer1$initial <- NULL
# toshiny.fewer1$final <- NULL
# toshiny.fewer1$JUNCTION_LENGTH0 <- NULL

toshiny.fewer1.h <- BX.clusters.mabsnoaugmenta.h %>% select(SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
# toshiny.fewer1.h$initial <- NULL
# toshiny.fewer1.h$final <- NULL
# toshiny.fewer1.h$JUNCTION_LENGTH0 <- NULL


BX.clusters.mabsaugmenta <- filter(BX.clusters.mabs, source == "Augmenta")
BX.clusters.mabsaugmenta.h <- subset(BX.clusters.mabsaugmenta, CREGION %in% c("IgH"))

BX.clusters.mabsaugmenta$ncount <- NULL
BX.clusters.mabsaugmenta$shm.mean <- NULL

BX.clusters.mabsaugmenta <- BX.clusters.mabsaugmenta %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

BX.clusters.augmentavsnoaugmenta <- bind_rows(BX.clusters.mabsnoaugmenta, BX.clusters.mabsaugmenta, .id = "id")
BX.clusters.augmentavsnoaugmenta$id <- gsub("2","collaborative mAbs",BX.clusters.augmentavsnoaugmenta$id)
BX.clusters.augmentavsnoaugmenta$id <- gsub("1","anti-CoV2 mAbs",BX.clusters.augmentavsnoaugmenta$id)
#BX.clusters.augmentavsnoaugmenta$id <- factor(BX.clusters.augmentavsnoaugmenta$id, levels = c("Healthy control bulk repertoire", "anti-CoV2 mAbs"))
BX.clusters.augmentavsnoaugmenta <- subset(BX.clusters.augmentavsnoaugmenta, CREGION %in% c("IgH"))

toshiny.augvsnoaug <- BX.clusters.augmentavsnoaugmenta %>% select(id,SEQUENCE_ID,source,binding,neutralization,CREGION,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))
# toshiny.augvsnoaug$initial <- NULL
# toshiny.augvsnoaug$final <- NULL
# toshiny.augvsnoaug$JUNCTION_LENGTH0 <- NULL


## for combined need to recalculate ncount and shm.mean
toshiny.augvsnoaugc <- toshiny.augvsnoaug
toshiny.augvsnoaugc$ncount <- NULL
toshiny.augvsnoaugc$shm.mean <- NULL

toshiny.augvsnoaugc <- toshiny.augvsnoaugc %>%
  add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n) %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE)) %>%
  mutate_at(vars(shm.mean), funs(round(., 2)))

#write.table(toshiny.augvsnoaug, "toshiny_augvsnoaug.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.augvsnoaugc, "toshiny_augvsnoaugc.tab", sep = "\t", row.names = FALSE, quote = FALSE)


## for John all-in-one
# write.table(toshiny.few, "toshiny_few.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(toshiny.few1, "toshiny_few1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(toshiny.few1.h, "toshiny_few1_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(toshiny.few2, "toshiny_few2.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(toshiny.few3, "toshiny_few3.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# 
# write.table(toshiny.fewer1, "toshiny_fewer1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(toshiny.fewer, "toshiny_fewer.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(toshiny.fewer1.h, "toshiny_fewer1_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## to re-arrange columns so name & junction_aa are first 2
#toshiny.few <- toshiny.few[c(1,3,2,4:16)]


#toshiny.few <- read.delim("toshiny_few.tab")

toshiny.few <- read.delim("toshiny_few.tab")
toshiny.few1 <- read.delim("toshiny_few1.tab")
toshiny.few1.h <- read.delim("toshiny_few1_h.tab")
toshiny.few2 <- read.delim("toshiny_few2.tab")
toshiny.few3 <- read.delim("toshiny_few3.tab")
toshiny.few3b <- read.delim("toshiny_few3b.tab")

## were mutate.if but now as.character but only after loading datasets
#toshiny.fewer %>% mutate_if(is.factor, as.character)
## rename in few2 id to sequence_id - also relabel IDs in HC & comet...
toshiny.few$SEQUENCE_ID <- as.character(toshiny.few$SEQUENCE_ID)
toshiny.few$JUNCTIONAA <- as.character(toshiny.few$JUNCTIONAA)
toshiny.few1$SEQUENCE_ID <- as.character(toshiny.few1$SEQUENCE_ID)
toshiny.few1$JUNCTIONAA <- as.character(toshiny.few1$JUNCTIONAA)
toshiny.few1.h$SEQUENCE_ID <- as.character(toshiny.few1.h$SEQUENCE_ID)
toshiny.few1.h$JUNCTIONAA <- as.character(toshiny.few1.h$JUNCTIONAA)

toshiny.few2$SEQUENCE_ID <- as.character(toshiny.few2$SEQUENCE_ID)
toshiny.few2$JUNCTIONAA <- as.character(toshiny.few2$JUNCTIONAA)
toshiny.few3$SEQUENCE_ID <- as.character(toshiny.few3$SEQUENCE_ID)
toshiny.few3$JUNCTIONAA <- as.character(toshiny.few3$JUNCTIONAA)
toshiny.few3b$SEQUENCE_ID <- as.character(toshiny.few3b$SEQUENCE_ID)
toshiny.few3b$JUNCTIONAA <- as.character(toshiny.few3b$JUNCTIONAA)


## to deploy app on ShinyApps.io  at: ewaltari.shinyapps.io
## MAKE SURE THAT THE TAB INPUT FILES ARE IN THE APP-2 FOLDER
library(rsconnect)
#rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-2')
#rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-3')

#rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-4')

#rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-5')

rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-6')

### FULL SET OF STEPS NEEDED FOR VIS TOOL:
# 1) RUN DATASET THROUGH IMMCANTATION (OR CUSTOM BUILD TO HAVE SIMILAR COLUMNS)
# 2) IMPORT DATASET (OR MULTIPLE)
# 3) RUN SOME ADDITIONAL COMMANDS TO GET KEY COLUMNS
#  -relabelling constant regions to be consistent
#  -defining cdr3 kabat & shm columns
#  -creating combined GF/JGENE column
#  -adding count of each GF/JGENE + CDR3length combination bin
#  -similarly calculating mean SHM for each of these bins
#   
# 4) IF COMBINING DATASETS BIND_ROWS COMMAND HERE, INCLUDES RENAMING COLUMN & ADDING FACTORS FOR INTENDED ORDER
#  - FOR FULL SET OF VIS OPTIONS WOULD NEED EVERY POSSIBLE COMBINATION OF MULTIPLE datasets
#  - ALSO WOULD BE GOOD TO HAVE .H .K .L SUBSETS OF ALL DATASETS IF COMBINING
#  - in addition might want to uncount by reads_per clone for bulk datasets if each row is a clonal family
# 5) RELABEL 'TOSHINY.FEW' TO BE WANTED DATASET
# 6) THEN CAN RUN CURRENT APP.R SCRIPT

## EVENTUALLY WOULD LIKE BOTH DROPDOWN MENU OPTIONS + HOVER/CLICK CAPABILITY THAT IS CURRENTLY LACKING
## ADDITIONALLY WOULD LIKE TO BE ABLE TO CREATE A NETWORK FOR SELECTED CLONAL FAMILY

BX.clusters.mabs <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_clonalclusters.tab")
#BX.clusters.mabs$CDR3KABAT <- ((BX.clusters.mabs$JUNCTION_LENGTH) / 3)
BX.clusters.mabs$CDR3LENGTH_IMGT <- ((BX.clusters.mabs$JUNCTION_LENGTH) / 3) - 2
BX.clusters.mabs$JGENE[BX.clusters.mabs$JGENE==""]<-NA
BX.clusters.mabs <- BX.clusters.mabs %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)
BX.clusters.mabs$CREGION <- gsub("IGK","Kappa",BX.clusters.mabs$CREGION)
BX.clusters.mabs$CREGION <- gsub("IGL","Lambda",BX.clusters.mabs$CREGION)
BX.clusters.mabs$CREGION <- gsub("IGH","IgH",BX.clusters.mabs$CREGION)
BX.clusters.mabs$SHM <- (100 - (BX.clusters.mabs$V_IDENTITY * 100))
## NOW ADDTHE COUNTS OF COMBINED GF_JGENE,CDR3LENGTH_IMGT
## BUT HAVE TO DO IT BEFORE COMBINING DATASETS...
BX.clusters.mabs <- BX.clusters.mabs %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)

## THIS ADDS MEAN SHM FOR EACH COMBINED GF_JGENE,CDR3LENGTH_IMGT 
BX.clusters.mabs.h <- subset(BX.clusters.mabs, CREGION %in% c("IgH"))
BX.clusters.mabs.h <- BX.clusters.mabs.h %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

  
#h3.mabsvshc.h <- ggplot(BX.clusters.mabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
#h3.mabsvshc.k <- ggplot(BX.clusters.mabsvshc.k, aes(GF_JGENE,CDR3LENGTH_IMGT))
#h3.mabsvshc.l <- ggplot(BX.clusters.mabsvshc.l, aes(GF_JGENE,CDR3LENGTH_IMGT))


### fewer columns for testing in Shiny 8/4/2020 - now incorporating only when calling it toshiny.few
BX.clusters.cometvsmabsvshc.h.few <- BX.clusters.cometvsmabsvshc.h %>% select(id,SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean)
BX.clusters.cometvsmabsvshc.h.few <- BX.clusters.cometvsmabsvshc.h.few %>% filter(!is.na(CDR3LENGTH_IMGT))
# !is.na()
#filter(starwars, mass > 1000)




############################################################################################################
### EARLY SUMMER CODE THAT MAKES INDIVIDUAL PLOTS - FOCUS IS NOW (FALL 2020) ON PREPARING DATASETS FOR SHINY
############################################################################################################
### plotting just one dataset

###########################
## mabs plots
## first tests
#h2 <- ggplot(BX.clusters.mabs, aes(GF_CDR3,JGENE))
#h2 + geom_bin2d() + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_gradient(low = "gray", high = "magenta") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## better comparisons
h3 <- ggplot(BX.clusters.mabs, aes(GF_JGENE,CDR3LENGTH_IMGT))
## gray to magenta
h3 + geom_bin2d() + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_gradient(low = "gray", high = "magenta") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#BX.clusters.mabs.gfbyjgenebycdr3 <- h3 + geom_bin2d() + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_gradient(low = "gray", high = "magenta") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabs.gfbyjgenebycdr3.pdf", BX.clusters.mabs.gfbyjgenebycdr3, width = 16, height = 12, units = "in")
h3 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_gradient(low = "gray", high = "magenta", name = "# of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## good
h3 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nReads",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## testing different viridis options
h3 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDR3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_gradient2(low = "gray", mid = "seagreen", high = "gold", midpoint = 1, name = "# of \nReads", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## D is default, A bad, B ok, C better, D default also good, E bad
BX.clusters.mabs.gfbyjgenebycdr3.viridis <- h3 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nReads", option = "D",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabs.gfbyjgenebycdr3.pdf", BX.clusters.mabs.gfbyjgenebycdr3.viridis, width = 16, height = 12, units = "in")

## all CoV2 mabs heavy & light
h3 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(BX.clusters.mabs, aes(GF_JGENE,CDR3LENGTH_IMGT)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


BX.clusters.mabs.gfbyjgenebycdr3.viridis <- h3 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.mabs.gfbyjgenebycdr3.pdf", BX.clusters.mabs.gfbyjgenebycdr3.viridis, width = 16, height = 12, units = "in")


## other plot color options
# replace
#scale_fill_viridis_c(name = "# of \nReads",  option = "D"
# with
#scale_fill_gradient2(name = "# of \nReads", low = "gray", mid = "seagreen", high = "gold"


## dot plot test
#h3 + geom_count() + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nClones") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


#geom_hex(aes(fill=log10(..count..))) + scale_fill_gradient(low = "cyan", high = "magenta", name = "Number of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000))

###########################
## HC plots
h3.hc <- ggplot(BX.clusters.hc, aes(GF_JGENE,CDR3LENGTH_IMGT))
#BX.clusters.hc.gfbyjgenebycdr3 <- h3.hc + geom_bin2d() + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_gradient(low = "gray", high = "magenta") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.hc.gfbyjgenebycdr3.pdf", BX.clusters.hc.gfbyjgenebycdr3, width = 16, height = 12, units = "in")
#h3.hc + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_gradient(low = "gray", high = "magenta", name = "# of Reads",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
h3.hc + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
h3.hc + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "% of \nClones",  option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

BX.clusters.hc.gfbyjgenebycdr3.bycount <- h3.hc + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("HC1_heavyandlight_bycount.pdf", BX.clusters.hc.gfbyjgenebycdr3.bycount, width = 16, height = 12, units = "in")
BX.clusters.hc.gfbyjgenebycdr3.bypercentage <- h3.hc + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "% of \nClones",  option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("HC1_heavyandlight_bypercentage.pdf", BX.clusters.hc.gfbyjgenebycdr3.bypercentage, width = 16, height = 12, units = "in")


#h3.hc + geom_bin2d(aes(fill=reads_per_clone)) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nClones, "option = "D") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## to cluster by reads not by clonal families...SEE BELOW HAVE TO UNCOUNT FIRST
#h3.hc2 <- ggplot(BX.clusters.hc, aes(GF_JGENE,CDR3LENGTH_IMGT, fill=reads_per_clone, group = reads_per_clone))
#h3.hc2 + geom_bin2d(aes(fill=log10(reads_per_clone))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nClones",  option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


BX.clusters.hc.gfbyjgenebycdr3.viridis <- h3.hc + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ CREGION, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nClones", option = "D",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.hc.gfbyjgenebycdr3.pdf", BX.clusters.hc.gfbyjgenebycdr3.viridis, width = 16, height = 12, units = "in")
## h3 histographs (V+J by CDR3), but now combinining 2 datasets into 1 graph...

##################################################

## test to plot not heat map of counts but dotplot with shm as color...

#BX.clusters.mabs.h <- subset(BX.clusters.mabs, CREGION %in% c("IgH"))
#mabs.h <- ggplot(BX.clusters.mabs.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
mabs.h.geompoint <- ggplot(BX.clusters.mabs.h, aes(x=GF_JGENE, y=CDR3LENGTH_IMGT, fill=shm.mean)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
  ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") +
  scale_y_continuous(limits = c(3, 42)) +
#  scale_fill_viridis_c(name = "Somatic \nHypermutation", option = "C") +
  geom_point(shape = 21, color = "white") +
  facet_wrap(~ binding, ncol=1, scales = "fixed")
mabs.h.geompoint + scale_fill_viridis_c(name = "Somatic \nHypermutation", option = "")
## this works but only single point for many reads?
## maybe instead flip between counts & shm
## below doesn't work with multiple SHM values

#mabs.h + geom_bin2d(aes(fill=SHM)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ binding, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Somatic \nHypermutation", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## need to separately calculate SHM for each gene/cdr3 combination in a separate calculation...

#iris %>% group_by(Species) %>% mutate(…)
BX.clusters.mabs.h <- BX.clusters.mabs.h %>%
  group_by(GF_JGENE,CDR3LENGTH_IMGT) %>%
  mutate(shm.mean = mean(SHM, na.rm = TRUE))

## too many na's given the dataset looks pretty good? geom_bin2d(aes(fill=log10(..count..)))
mabs.h + geom_bin2d(aes(fill=shm.mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ binding, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Somatic \nHypermutation", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
## maybe use regular plot and have shm.mean values on the bar?
mabs.h + geom_bin2d(geom = "text",aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ binding, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


###############
## NOW ADDING 8/7 THE COUNTS OF COMBINED GF_JGENE,CDR3LENGTH_IMGT
## BUT HAVE TO DO IT BEFORE COMBINING DATASETS...
# BX.clusters.cometvsmabsvshc.h.few <- BX.clusters.cometvsmabsvshc.h.few %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
#   rename(ncount = n)

BX.clusters.hc <- BX.clusters.hc %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
comet.pos <- comet.pos %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
comet.neg <- comet.neg %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
BX.clusters.mabs <- BX.clusters.mabs %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)
BX.clusters.cometvsmabsvshc <- bind_rows(BX.clusters.hc, comet.pos, comet.neg, BX.clusters.mabs, .id = "id")
BX.clusters.cometvsmabsvshc$id <- gsub("1","Healthy control bulk repertoire",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("2","COMET COVID-positive",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("3","COMET COVID-negative",BX.clusters.cometvsmabsvshc$id)
BX.clusters.cometvsmabsvshc$id <- gsub("4","anti-CoV2 mAbs",BX.clusters.cometvsmabsvshc$id)

BX.clusters.cometvsmabsvshc.h <- subset(BX.clusters.cometvsmabsvshc, CREGION %in% c("IgH"))

BX.clusters.cometvsmabsvshc.h$id <- factor(BX.clusters.cometvsmabsvshc.h$id, levels = c("Healthy control bulk repertoire", "COMET COVID-positive", "COMET COVID-negative", "anti-CoV2 mAbs"))

BX.clusters.cometvsmabsvshc.h.few <- BX.clusters.cometvsmabsvshc.h %>% select(id,SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean)
### new function
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
BX.clusters.cometvsmabsvshc.h.few <- BX.clusters.cometvsmabsvshc.h.few %>% filter(is.wholenumber(CDR3LENGTH_IMGT))
##############

h3.cometvsmabsvshc.h2 <- ggplot(BX.clusters.cometvsmabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw(base_size = 10) + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

#h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill=..count../sum(..count..), group=id)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill=stat(..count../sum(..count..), group=id))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill=..count.., group=id)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#geom_histogram(aes(y=0.01*..density..), binwidth = 0.01)  aes(y = stat(density*width))
#to get around use of facets use geom_bar(aes(y = (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) instead. Each facet should sum to 100% 
## THIS WORKS - GETS PERCENT FOR EACH FACET...
h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3.viridis.h <-  h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3.viridis.h2 <-  h3.cometvsmabsvshc.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#ggsave("BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_h.pdf", BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3.viridis.h2, width = 16, height = 12, units = "in")
#ggsave("BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_h_bypercentage.pdf", BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3.viridis.h, width = 16, height = 12, units = "in")



#####################
### Briney data

BX.clusters.briney <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/Briney_control1_clonalclusters.tab")
#BX.clusters.briney$CDR3KABAT <- ((BX.clusters.briney$JUNCTION_LENGTH) / 3)
BX.clusters.briney$CDR3LENGTH_IMGT <- ((BX.clusters.briney$JUNCTION_LENGTH) / 3) - 2
BX.clusters.briney <- BX.clusters.briney %>% unite(GF_JGENE, GF, JGENE, sep = "_", remove = FALSE, na.rm = TRUE)

## to uncount
BX.clusters.briney.2 <- uncount(BX.clusters.briney, BX.clusters.briney$reads_per_clone)


h3.briney.h2 <- ggplot(BX.clusters.briney.2, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.briney.h2 + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "# of \nReads", option = "D", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#h3.briney.h2 + geom_bin2d(aes(fill= 10*(..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "‱ of \nReads", option = "C", trans = "log") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


###########################
## to separately calculate rather than just plot:
stat_bin2d(): geom_bin2d()

## also trying stat_summary_2d
d <- ggplot(diamonds, aes(carat, depth, z = price))
d + stat_summary_2d()

# Specifying function
d + stat_summary_2d(fun = function(x) sum(x^2))
d + stat_summary_2d(fun = ~ sum(.x^2))

#ggplot(BX.clusters.mabsvshcvsnat, aes(GF_JGENE,CDR3LENGTH_IMGT, z = ..count..)) + stat_summary_2d()

h3.cometvsmabsvsnat.h2 <- ggplot(BX.clusters.mabsvshcvsnat, aes(GF_JGENE,CDR3LENGTH_IMGT))
cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h2 <-  h3.cometvsmabsvsnat.h2 + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ id, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "% of \nReads", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
#write.table(cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h2$data, "BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_htest.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h2$coordinates, "BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_htestcoords.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#ggsave("mabsvshcvsnatalia.gfbyjgenebycdr3_h2.pdf", cometvsmabsvsnatalia.gfbyjgenebycdr3.viridis.h2, width = 16, height = 12, units = "in")

d <- ggplot(diamonds, aes(x, y)) + xlim(4, 10) + ylim(4, 10)
d.stat <- ggplot(diamonds, aes(x, y)) + xlim(4, 10) + ylim(4, 10) + stat_bin2d()

d.stat2 <- ggplot_build(d.stat)$data

p <- ggplot(mpg, aes(x=factor(cyl), y=..count..))
p + geom_bar(stat="bin")
p + geom_histogram()   

## this gets the count ( is 'n') ### THIS WORRRKSSSS!!!!!
t <- BX.clusters.mabs.h %>% group_by(GF_JGENE,CDR3LENGTH_IMGT) %>% tally
## next try a join to get the ncount into the original dataset????
## or try this:
BX.clusters.mabs.h <- BX.clusters.mabs.h %>% add_count(GF_JGENE,CDR3LENGTH_IMGT) %>%
  rename(ncount = n)


p <- ggplot(diamonds,aes(carat,price)) + stat_summary2d(fun=sum,aes(z=depth))
#Use function ggplot_build() to get data used for plot. Coordinates of rectangles are in columns xmin, xmax, ymin and ymax and sum are in column value.

df <- ggplot_build(p)$data[[1]]


# legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10)

## testing
#write.table(h3.cometvsmabsvshc.h2.stat$data, "BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#BX.clusters.mabs.gfbyjgenebycdr3.viridis <- h3 + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + facet_wrap(~ CREGION, ncol=1, scales = "free") + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nReads", option = "D",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## this works and ACTUALLY CAN GET PER FACET SEE BELOW
#h3.cometvsmabsvshc.h2 <- ggplot(BX.clusters.cometvsmabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT))
h3.cometvsmabsvshc.h2.stat <- h3.cometvsmabsvshc.h2 + stat_bin2d(aes(fill=..count..)) + facet_wrap(~ id)
h3.cometvsmabsvshc.h2.stattable <- ggplot_build(h3.cometvsmabsvshc.h2.stat)
write.table(h3.cometvsmabsvshc.h2.stattable$data, "BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_htable.tab", sep = "\t", row.names = FALSE, quote = FALSE)

h3.cometvsmabsvshc.h2.stattable2 <- ggplot_build(h3.cometvsmabsvshc.h2.stat)$data
write.table(h3.cometvsmabsvshc.h2.stat, "BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

newdat <- ggplot_build(h3.cometvsmabsvshc.h2.stat)$data[[1]]

#h3.cometvsmabsvshc.h2.stat <- h3.cometvsmabsvshc.h2 + stat_bin2d(aes(fill=..count.., group = (..count..)[..PANEL..])) + facet_wrap(~ id)
## next 2 don't work
#h3.cometvsmabsvshc.h2.stat <- h3.cometvsmabsvshc.h2 + stat_bin2d(aes(fill=..count.., group = (..count..)*(reads_per_clone)[..PANEL..])) + facet_wrap(~ id)
#h3.cometvsmabsvshc.h2.stat <- h3.cometvsmabsvshc.h2 + stat_bin2d(aes(fill=..count.., group = (..count..)*(reads_per_clone))) + facet_wrap(~ id)
#h3.cometvsmabsvshc.h2b <- ggplot(BX.clusters.cometvsmabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT, group = (..count..)*(reads_per_clone)))

#h3.cometvsmabsvshc.h2b <- ggplot(BX.clusters.cometvsmabsvshc.h, aes(GF_JGENE,CDR3LENGTH_IMGT, fill = reads_per_clone))
#h3.cometvsmabsvshc.h2.stat <- h3.cometvsmabsvshc.h2b + stat_bin2d(aes(fill=..count..)) + facet_wrap(~ id)
#h3.cometvsmabsvshc.h2.stat <- ggplot_build(h3.cometvsmabsvshc.h2.stat)$data
write.table(h3.cometvsmabsvshc.h2.stat, "BX.clusters.cometvsmabsvshc.gfbyjgenebycdr3_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
## trying to count reads_per_clone, but this doesn't look different than before - will likely have to uncount still

## to visualize x & y from  BX.clusters.cometvsmabsvshc.h.few
## but z is h3.cometvsmabsvshc.h2.stat$count


  # Specifying function
#d <- ggplot(diamonds, aes(carat, depth, z = price))
#d + stat_summary_2d(fun = function(x) sum(x^2))
#If you want the dataframe that has positions and other information used in assembling geoms, use
#ggplot_build(p)$data
### also plots of top clonal clusters



############################################################################################################
############################################################################################################

## simple phylogenies test commands

# library(treeio)
# library(tidytree) 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree")
#library(ggtree)

library(seqinr)
library(phangorn)
library(ape)

### commmands for plotting tree of CDR3 motifs after extracting dataset...
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"SEQUENCE_ID"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"JUNCTIONAA"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
# writeFasta2<-function(data, filename){
#   fastaLines = c()
#   for (rowNum in 1:nrow(data)){
#     fastaLines = c(fastaLines, as.character(paste(">", data["SEQUENCE_ID"], sep = "")))
#     fastaLines = c(fastaLines,as.character(data["JUNCTIONAA"]))
#   }
#   fileConn<-file(filename)
#   writeLines(fastaLines, fileConn)
#   close(fileConn)
# }


exampleData <- read.csv("filtereddata.csv")

exampleData <- read.csv("filtereddata.csv") %>% mutate_all(as.character)

## relabeling tips for tree (note if using HC, do not need initial sequence id at all!)
exampleData$SEQUENCE_ID <- paste(exampleData$SEQUENCE_ID,exampleData$GENE,exampleData$JUNCTIONAA,sep="_") 


#writeFasta(exampleData, "example.fasta")

# ### test to convert rather than save and have to reload - doesn't work
# convertFasta<-function(data){
#   fastaLines = c()
#   for (rowNum in 1:nrow(data)){
#     fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"SEQUENCE_ID"], sep = "")))
#     fastaLines = c(fastaLines,as.character(data[rowNum,"JUNCTIONAA"]))
#   }
#   return(fastaLines)
# }
# test <- convertFasta(exampleData)
test <- as.data.frame(test)
# dat.test = read.fasta(test, type = "AA")
as.AAbin(test)
# dat.test = read.phyDat(test, format = "fasta", type = "AA")
dat.test = as.phyDat(test, format = "fasta", type = "AA")
dat.test = as.phyDat(test, type = "AA")
dat.test = as.AAbin(test, format = "fasta")

## this works, trying to replicate above without saving as .fasta
dat = read.phyDat("example.fasta", format = "fasta", type = "AA")

## read.phyDat is from phangorn package
## as.AAbin is from Ape
## as.AAbin.phyDat is from phangorn

## trying to skip writing to .fasta and the re-reading - but doesn't work...
# x <- structure(c("55548", "43297", "35309", "34468", "AATTCAATGCTCGGGAAGCAAGGAAAGCTGGGGACCAACTTCTCTTGGAGACATGAGCTTAGTGCAGTTAGATCGGAAGAGCA", "AATTCCTAAAACACCAATCAAGTTGGTGTTGCTAATTTCAACACCAACTTGTTGATCTTCACGTTCACAACCGTCTTCACGTT", "AATTCACCACCACCACTAGCATACCATCCACCTCCATCACCACCACCGGTTAAGATCGGAAGAGCACACTCTGAACTCCAGTC", "AATTCTATTGGTCATCACAATGGTGGTCCGTGGCTCACGTGCGTTCCTTGTGCAGGTCAACAGGTCAAGTTAAGATCGGAAGA"), .Dim = c(4L, 2L))
# y <- t(sapply(strsplit(x[,2],""), tolower))
# rownames(y) <- x[,1]
# dat.test <- as.DNAbin(y) ## but not a list!
# 
# 
# exampleData.y <- t(sapply(strsplit(exampleData[,2],""), tolower))
# rownames(exampleData.y) <- exampleData[,1]
# as.DNAbin(y)
#Then as.DNAbin(y) shows
dat.test = as.phyDat(y, type = "DNA")
dat.test = as.phyDat(y, type = "AA")


writeFasta(exampleData, "example.fasta")


### THIS WORKS AT CONVERTING!!

## 2021 start here for testing
toshiny.fewhivnih45 <- read.delim("toshiny_fewhivnih45.tab")
exampleData <- read.delim("toshiny_dengue_cf10.tab")

exampleData$SEQUENCE_ID <- as.character(exampleData$SEQUENCE_ID)
exampleData$JUNCTIONAA <- as.character(exampleData$JUNCTIONAA)

exampleData <- exampleData %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT)

exampleData$SEQUENCE_ID <- paste(exampleData$SEQUENCE_ID,exampleData$GENE,exampleData$JUNCTIONAA,sep="_")

exampleData.y <- t(sapply(strsplit(exampleData[,2],""), tolower))
rownames(exampleData.y) <- exampleData[,1]
as.AAbin(exampleData.y)
# dat.test = as.phyDat(exampleData.y, type = "AA")
# seqsall = as.phyDat(filteredDataall.y, type = "AA")
seqsall = as.phyDat(exampleData.y, type = "AA")
dmall = dist.ml(seqsall, model="JTT")

## need to now sort by single read?
# dm.matrix <- as.matrix(dm) %>% as.data.frame() %>% arrange_at("C003_IGHV3-53_ARDYGDFYFDY")

## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
# filteredData1ID <- isolate({filteredDSpartial2id()})
dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at("mab503-B10")
dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at("mab503-B10_IGHV4-39_ARQDRNWFDS")
dmall.matrix <- rownames_to_column(dmall.matrix, var = "SEQUENCE_ID")
## in between here need to trim to fewer sequences using distance threshold...

# dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID, filteredData1ID)
# dmall.matrix[2] = lapply(dmall.matrix[2], as.numeric)
# dmall.matrix <- dmall.matrix %>% filter(filteredData1ID < 0.5)

dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID, "mab503-B10_IGHV4-39_ARQDRNWFDS")
# dmall.matrix[2] = lapply(dmall.matrix[2], as.numeric)
## try renaming 2nd column...THIS WORKED!!!! 3/2/2021
colnames(dmall.matrix)[2] <- "DIST"
#dmall.matrix <- rename.vars(dmall.matrix, c("SEQUENCE_ID","mab503-B10"), c("SEQUENCE_ID","DIST"))
dmall.matrix <- dmall.matrix %>% filter(DIST < 0.5)

## this works
#dmall.matrix <- dmall.matrix %>% filter(SEQUENCE_ID == "mab503-B10")
## but this doesn't
#dmall.matrix <- dmall.matrix %>% filter(filteredData1ID < 0.5)

## then back to code to trim to first 50 sequences...
dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID)
dm.matrix <- dmall.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
  select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 80) ## okay if this is more than what is in the set!

## last tree plotting steps

filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
rownames(filteredData.y) <- dm.matrix[,1]
as.AAbin(filteredData.y)
seqs = as.phyDat(filteredData.y, type = "AA")
dm = dist.ml(seqs, model="JTT")
tree1 = NJ(dm)
tree1$edge.length[tree1$edge.length<0] <- 0
tree <- pratchet(seqs, start = tree1, maxit=100,
                 minit=5, k=5, trace=0)
par(family = "mono", mar=c(2, 0.0, 2, 0) + 0.0)
### this changes the colors of all 'controls' to gray but leaves covid mabs black
tipcolors <- def(tree$tip.label, "HC" = "gray", "Boyd" = "black", "Galson" = "black", "clonotype" = "black", "N152" = "thistle", "IAVI84" = "plum", "AD1214" = "orchid", "nih45" = "pink", "-au" = "brown", "d13" = "orange", "d20" = "blue", "OAS" = "gold", "mab" = "orchid", default = "orange", regexp = TRUE)
# plot(bird.orders, tip.color = co2)
plot(tree, lab4ut="axial",
     edge.width=2, y.lim = c(-10, 10), label.offset = 0, cex = 0.70, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors)  ## x.lim breaks parsimony and some NJ
## y.lim Limits on the vertical axis NULL (d), two numeric values x.lim = c(-0.1, 2.2)

## trying to make tips closer together
# split.width = 1 DOESN'T WORK
# edge.width = 1  DOESN'T WORK
# y.lim = c(-10, 10),  WORKS - BUT ONLY FOR TREES WITH VERY FEW TIPS, NOT SURE IF THERE IS A WAY TO AUTOMATE THE SPACING...
## trying to understand best JTT threshold:
## plotting a full distance matrix
CF6_matrix <- read.csv("CF6_matrix.csv")
CF6_matrix <- read.csv("CF1_matrix.csv")
# melt(CF6_matrix) -> CF6_matrix_melted
CF6_matrix_melted <- melt(replace(CF6_matrix, lower.tri(CF6_matrix, TRUE), NA), na.rm = TRUE)
c <- ggplot(CF6_matrix_melted, aes(value))
c + geom_histogram(bins = 80) + scale_y_log10()



## using column names if possible - is not, need to reorder in shiny dataset first
# exampleData.y <- t(sapply(strsplit(exampleData[,JUNCTIONAA],""), tolower))
# rownames(exampleData.y) <- exampleData[,SEQUENCE_ID]
# as.AAbin(exampleData.y)
# dat.test = as.phyDat(exampleData.y, type = "AA")

#### jan 2021 ideas:
## change plot title to name of V+J+CDR3length?? key is here:
#filteredData <- isolate({filteredDSpartial2()}) - note this line is in 2 of 3 options, would need to add to third
## then near bottom of plot code something similar to this line to add to actual title:
# gjcdr3.title <- filteredData %>% separate(SEQUENCE_ID, into = c("GF_JGENE", "CDR3LENGTH_IMGT"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
#   select(SEQUENCE_ID)
exampleData$CDR3KABAT <- as.numeric(exampleData$CDR3KABAT)
exampleData$CDR3LENGTH_IMGT <- exampleData$CDR3KABAT - 2

## this should do it...
exampleData$CDR3LENGTH_IMGT <- as.numeric(exampleData$CDR3LENGTH_IMGT)
exampleData$G_J_CDR3 <- paste("Phylogeny of selected CDR3 motifs:", exampleData$GF_JGENE,exampleData$CDR3LENGTH_IMGT,"aa",sep=" ")

gjcdr3.title <- exampleData %>%
  select(G_J_CDR3)
gjcdr3.titleID <- gjcdr3.title$G_J_CDR3[1]

## also need to rename sequence_id in .hc and comet, then recalculate all toshiny files...
toshiny.fewer <- BX.clusters.mabsvshcnoaugmenta.h %>% select(id,SEQUENCE_ID,binding,neutralization,JUNCTIONAA,GENE,GF_JGENE,GF,JGENE,JGF,CDR3LENGTH_IMGT,SHM,ncount,shm.mean,reads_per_clone) %>% filter(!is.na(CDR3LENGTH_IMGT)) %>% filter(is.wholenumber(CDR3LENGTH_IMGT)) %>% mutate_at(vars(SHM), funs(round(., 2))) %>% mutate_at(vars(shm.mean), funs(round(., 2)))



## alternate function
# TabularToFasta <- function (filename){
#   file <- read.csv (file=filename, header=TRUE)
#   
#   file = as.data.frame(file)
#   #delete if any existing file 
#   
#   unlink(c("dna_fasta.fasta"), force=TRUE)
#   
#   #give output filename
#   
#   sink("dna_fasta.fasta")
#   
#   for (i in 1:nrow(file)){
#     name = paste0(">",file[i,1])
#     sequence = paste(file[i,2])
#     cat(name,sep="\n")
#     cat(sequence,sep="\n")
#   }
#   #this is sink all the console output to the file 
#   sink()
# }


##########

### 9/22 thoughts on selecting single sequence and sorting all sequences in a bin by distance to select sequence:

# will need to use filteredDSall2
# then will need to make dist.ml on this full dataset, then sort the matrix and somehow extract SMALLEST numbers...

##You can coerce it into a matrix, and then sort it
##Then as @Roland suggested, convert it back to dist

#T.mat <- as.matrix(T)[ordering, ordering]

#T  <- as.dist(T.mat)

dm.matrix <- as.matrix(dm)

## think arrange(matrix, selected_sequence) will arrange from smallest to largest the distances...

## then will want to extract just first column of first xx rows (let's say 25 but could change)

## this new set will be a new dataset to run the commands on - could call it filteredDSnearest2
## then tree making commands would go under a third option after NJ & parsimony...



##dm.matrix <- as.data.frame(dm) can't do this
#dm.matrix <- as.matrix(dm) %>% arrange(filteredDSpartial2id())

#dm.matrix <- as.matrix(dm) %>% as.data.frame()

#dm.matrix <- as.matrix(dm) %>% as.data.frame() %>% arrange("C003_IGHV3-53_ARDYGDFYFDY")

dm.matrix <- as.matrix(dm) %>% as.data.frame() %>% arrange_at("C003_IGHV3-53_ARDYGDFYFDY")



### THURSDAY ADD FOLLOWING TO GENERIC CODE BELOW...
## first make the rownames their own column again
dm.matrix <- rownames_to_column(dm.matrix, var = "SEQUENCE_ID")
#dm.matrix$SEQUENCE_ID <- as.character(dm.matrix$SEQUENCE_ID)
#dm.matrix$more <- "filler"
dm.matrix <- dm.matrix %>% select(SEQUENCE_ID)

## now split rownames into 3 components but keep full name (NEED TO HAVE NO UNDERSCORES IN REGULAR NAMES NOW)
# then select just first column and fourth column, also just first xxx rows (6 in this example for name issue above)
dm.matrix <- dm.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>% 
  select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 5) ## okay if this is more than what is in the set!

#mtcars %>% slice_head(n = 5)

## now to tree making test...
exampleData.y2 <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
rownames(exampleData.y2) <- dm.matrix[,1]
as.AAbin(exampleData.y2)
dat.test2 = as.phyDat(exampleData.y2, type = "AA")
dm2 = dist.ml(dat.test2, model="JTT")
tree2 = midpoint(NJ(dm2))
tree2$edge.length[tree2$edge.length<0] <- 0
par(family = "mono", mar=c(4, 1, 4, 0) + 0.5)
tipcolors <- def(tree2$tip.label, "-66" = "orange", default = "black", regexp = TRUE)
# plot(bird.orders, tip.color = co2)
plot(midpoint(tree2), lab4ut="axial",
     edge.width=1, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4, tip.color = tipcolors)
title("NJ phylogeny of CDR3 motifs", family = "sans")

## trying to make tips closer together
split.width = 1
edge.width = 1

# filteredData1ID <- "C003_IGHV3-53_ARDYGDFYFDY"
# dm.matrix2 <- arrange(dm.matrix, filteredData1ID)
# mtcars2 <- arrange(mtcars, qsec)
# mtcars2 <- rownames_to_column(mtcars, var = "cars") %>% select(cars)
# mtcars2 %>% separate(cars, into = c("cars1", "cars2", "cars3"), sep = " ", remove = FALSE, convert = TRUE, extra = "merge", fill = "left")



## would need to use this near top of app code for selected_sequence above:
filteredDSpartial2id <- reactive({
  partial2id <- filteredDSpartial() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT)
  #      partial2 <- filteredDSpartial() %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT) %>% mutate_if(is.factor, as.character)
  partial2id$SEQUENCE_ID <- paste(partial2id$SEQUENCE_ID,partial2id$GENE,partial2id$JUNCTIONAA,sep="_")
  partial2id$JUNCTIONAA <- NULL
  partial2id$GENE <- NULL
  partial2id$GF_JGENE <- NULL
  partial2id$CDR3LENGTH_IMGT <- NULL
  partial2id
})

### this will be 3rd tree function
} else {
  ## first make matrix of all sequences in the bin
  ## IT WOULD BE NICE TO HAVE A DISTANCE THRESHOLD, BUT EXAMINING THE CALCULATED MATRIX SHOWS THIS IS VERY TRICKY
   ## YOU'D HAVE TO SEPARATELY EXTRACT THE DISTANCE OF EACH SEQUENCE AGAINST THE TARGET, NOT CURRENTLY DONE
  filteredDataall <- isolate({filteredDSall2()})
  filteredDataall.y <- t(sapply(strsplit(filteredDataall[,2],""), tolower))
  rownames(filteredDataall.y) <- filteredDataall[,1]
  as.AAbin(filteredDataall.y)
  seqsall = as.phyDat(filteredDataall.y, type = "AA")
  dmall = dist.ml(seqsall, model="JTT")
  ## now sort this matrix by single selected sequence (just the sequence_id, so it is filteredDSpartial2id from above)
  filteredData1ID <- isolate({filteredDSpartial2id()})
  dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at(filteredData1ID)
  # dmall.matrix <- as.matrix(dmall) %>% as.data.frame() %>% arrange_at("Boyd7450-47660_IGHV4-4_ARDRGIAAAADHPFGYYYYYMDV")
  dmall.matrix <- rownames_to_column(dmall.matrix, var = "SEQUENCE_ID")
  # dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID)
  # dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID, "Boyd7450-47660_IGHV4-4_ARDRGIAAAADHPFGYYYYYMDV")
  dmall.matrix$nn <- dmall.matrix$"Boyd7450-47660_IGHV4-4_ARDRGIAAAADHPFGYYYYYMDV"
  # dmall.matrix <- dmall.matrix %>% filter(nn < 0.75)
  dmall.matrix <- dmall.matrix %>% select(SEQUENCE_ID, filteredData1ID)
  dmall.matrix <- dmall.matrix %>% filter(filteredData1ID < 0.75)

    dm.matrix <- dmall.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>% 
    select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 15) ## okay if this is more than what is in the set!
  ## now with new subset make new matrix & tree
  filteredData.y <- t(sapply(strsplit(dm.matrix[,2],""), tolower))
  rownames(filteredData.y) <- filteredData[,1]
  as.AAbin(filteredData.y)
  seqs = as.phyDat(filteredData.y, type = "AA")
  dm = dist.ml(seqs, model="JTT")
  tree = midpoint(NJ(dm))
  tree$edge.length[tree$edge.length<0] <- 0
}


# use embracing when wrapping in a function;
# see ?dplyr_data_masking for more details
# tidy_eval_arrange <- function(.data, var) {
#   .data %>%
#     arrange({{ var }})
# }
# tidy_eval_arrange(mtcars, mpg)

###########

#dat = read.phyDat("example.fasta", format = "fasta", type = "AA")
#dat = read.phyDat("Twist.fasta", format = "fasta", type = "AA")
dm = dist.ml(dat.test, model="JTT")
tree = NJ(dm)
#plot(tree, main = "Neighbor Joining")
tree$edge.length[tree$edge.length < 0] <- 0

# plot(tree, lab4ut="axial",
#      edge.width=2, label.offset = 0.2)

## for app cases where trees are cut off - it is usually because of some superlong branch lengths.
## maybe rescale all of the lengths??
#tree$edge.length <- tree$edge.length * 0.75
# tree$edge.length <- tree$edge.length ^ 0.75
# tree$edge.length <- tree$edge.length ^ 0.25
# tree$edge.length <- tree$edge.length * 0.35
tree$edge.length <- tree$edge.length ^ 0.2

## idea to make trees narrower - use log DOESN'T SEEM TO HELP
tree$edge.length[tree$edge.length < 0] <- 0.01
tree$edge.length[tree$edge.length == 0] <- 0.01

tree$edge.length <- tree$edge.length + 1
tree$edge.length <- log(tree$edge.length)




### version in app
midpoint(tree)
plot(tree, lab4ut="axial",
     edge.width=2, label.offset = 0.2, cex = 0.75, align.tip.label = 0, adj = 1, no.margin = FALSE, x.lim = 8)

par(family = "mono", mar=c(3, 3, 3, 3) + 0.5)
plot(midpoint(tree), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
title("NJ phylogeny of CDR3 motifs", family = "sans")


## using 4,1,4,0 in app script...
par(family = "mono", mar=c(4, 1, 4, 0) + 0.5)
plot(midpoint(tree), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
title("NJ phylogeny of CDR3 motifs", family = "sans")

par(family = "mono", mar=c(4, 4, 4, 4) + 1.5)
plot(midpoint(tree), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
title("NJ phylogeny of CDR3 motifs", family = "sans")


### this changes the colors of some tips
par(family = "mono", mar=c(4, 1, 4, 0) + 0.5)
tipcolors <- def(tree$tip.label, "-66" = "orange", default = "black", regexp = TRUE)
# plot(bird.orders, tip.color = co2)
plot(midpoint(tree), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4, tip.color = tipcolors)
title("NJ phylogeny of CDR3 motifs", family = "sans")


## using 2,1,2,0 in app script...but some trees are still cut off...
par(family = "mono", mar=c(2, 0.0, 2, 0) + 0.0)
plot(midpoint(tree), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
title("NJ phylogeny of CDR3 motifs", family = "sans")


###################
### tests with different dataset on fitting long names + long branch lengths in tree
library(DT)

## for plotting
#library(Biostrings)
library(seqinr)
library(phangorn)
library(ape)

exampleDatab <- read.csv("Filtereddata-all.csv") %>% select(SEQUENCE_ID,JUNCTIONAA,GENE,GF_JGENE,CDR3LENGTH_IMGT) %>% mutate_all(as.character)
## relabeling tips for tree (note if using HC, do not need initial sequence id at all!)
exampleDatab$SEQUENCE_ID <- paste(exampleDatab$SEQUENCE_ID,exampleDatab$GENE,exampleDatab$JUNCTIONAA,sep="_") 

# exampleData$SEQUENCE_ID <- as.character(exampleData$SEQUENCE_ID)
# exampleData$JUNCTIONAA <- as.character(exampleData$JUNCTIONAA)


## separate test in nov 2020 to get just full matrix downloaded separately
#as.AAbin(filteredData.y)
seqsall = as.phyDat(filteredData.y, type = "AA")
dmalldf = dist.ml(seqsall, ratio = TRUE, model="JTT")
dmalldf = dist.hamming(seqsall, ratio = TRUE)
# as.phyDat.data.frame(dmalldf)
# writeDist(dmalldf, file = "dmalldffile")
# dmalldffile <- readDist("dmalldffile")
# 
# dmalldffile <- read.alignment(file = "dmalldffile", format = "phylip")
## this works
dm.matrix <- as.matrix(dmalldf) %>% as.data.frame()

#####


exampleDatab.y <- t(sapply(strsplit(exampleDatab[,2],""), tolower))
rownames(exampleDatab.y) <- exampleDatab[,1]
as.AAbin(exampleDatab.y)
seqsb = as.phyDat(exampleDatab.y, type = "AA")
dmb = dist.ml(seqsb, model="JTT")


dmb.matrix <- as.matrix(dmb) %>% as.data.frame() %>% arrange_at("COV2-2335_IGHV4-39_ARHDGSGEMDTITWGPIYYYMDV")

dmb.matrix <- rownames_to_column(dmb.matrix, var = "SEQUENCE_ID")
dmb.matrix <- dmb.matrix %>% select(SEQUENCE_ID)

dmb.matrix <- dmb.matrix %>% separate(SEQUENCE_ID, into = c("SEQ", "GENE", "JUNCTIONAA"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>% 
  select(SEQUENCE_ID, JUNCTIONAA) %>% slice_head(n = 20) ## okay if this is more than what is in the set!


filteredData.y <- t(sapply(strsplit(dmb.matrix[,2],""), tolower))
rownames(filteredData.y) <- dmb.matrix[,1]
as.AAbin(filteredData.y)
seqsbb = as.phyDat(filteredData.y, type = "AA")
dmbb = dist.ml(seqsbb, model="JTT")
treeb = midpoint(NJ(dmbb))
treeb$edge.length[treeb$edge.length<0] <- 0

par(family = "mono", mar=c(2, 0.0, 2, 0) + 0.0)
### this changes the colors of all 'controls' to gray but leaves covid mabs black
tipcolors <- def(tree$tip.label, "HC" = "gray", "Boyd" = "black", "Galson" = "black", "clonotype" = "black", "N152" = "thistle", "IAVI84" = "plum", "AD1214" = "orchid", "nih45" = "pink", "-au" = "brown", "d13" = "orange", "d20" = "blue", "OAS" = "gold", "mab" = "orchid", default = "orange", regexp = TRUE)
plot(treeb, lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.95, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors)  ## x.lim breaks parsimony and some NJ
title("Phylogeny of CDR3 motifs", family = "sans")


treeb$edge.length <- treeb$edge.length * 0.15

treeb$edge.length <- treeb$edge.length ^ 0.25


par(family = "mono", mar=c(2, 0.0, 2, 0) + 0.0)
### this changes the colors of all 'controls' to gray but leaves covid mabs black
tipcolors <- def(treeb$tip.label, "HC" = "gray", "Boyd" = "black", "Galson" = "black", "clonotype" = "black", default = "orange", regexp = TRUE)
plot(treeb, lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.95, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors)  ## x.lim breaks parsimony and some NJ
title("Phylogeny of CDR3 motifs", family = "sans")



par(family = "mono", mar=c(2, 0.0, 2, 0) + 0.0)
plot(midpoint(treeb), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 10, font = 4)
title("NJ phylogeny of CDR3 motifs", family = "sans")

############################################################################
############################################################################
## lots more testing mostly not used in app
plot(tree, lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 1, align.tip.label = TRUE, adj = 1, no.margin = FALSE, font = 4, tip.color = tipcolors)

library(ape)

# Plotting in ggtree
# BiocManager::install("ggtree")
#library(ggtree)

mytree <- ggtree(tree, layout="equal_angle", size=0.5, linetype=1)
mytree

mytree <- ggtree(dm, layout="daylight")


trees <- getTrees(clones[1:2,])

# simple tree plotting with ggtree R package
plots <- plotTrees(trees, tips="c_call")


## trying parsimony in fangorn
ptree <- pratchet(dat.test, start = tree, maxit=100,
                  minit=5, k=5, trace=0) 
## BELOW UPDATE TO PAR COMMAND IS THE TRICK...
par(family = "mono", mar=c(3, 3, 3, 3) + 0.5)
plot(midpoint(ptree), lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, no.margin = FALSE, x.lim = 24, family = "mono", font = 4)

plot(ptree, lab4ut="axial",
     edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, family = "mono", font = 4)
title("Parsimony phylogeny of CDR3 motifs", family = "sans")

data(Laurasiatherian)
dm2 <- dist.hamming(Laurasiatherian)
tree2 <- NJ(dm2)
parsimony(tree2, Laurasiatherian)
#treeRA <- random.addition(Laurasiatherian)
#treeNNI <- optim.parsimony(tree, Laurasiatherian)
treeRatchet <- pratchet(Laurasiatherian, start=tree2, maxit=100,
                        minit=5, k=5, trace=0)

plot(midpoint(treeRatchet))

### trying tidytree to plot trees in ggplot:
## LOTS OF INTERESTING OPTIONS, BUT CAN'T RIGHT JUSTIFY THE LABELS...
library(tidytree)
tree.tidy <- as_tibble(tree)

ggtree(tree.tidy)$data
#The function, ggtree, was implemented as a short cut to visualize a tree, and it works exactly the same as shown above.

#ggtree takes all the advantages of ggplot2. For example, we can change the color, size and type of the lines as we do with ggplot2 (Figure 4.1B).

#ggplot(tree.tidy, aes(x, y)) + geom_tree() + theme_tree()
ggtree(tree)

ggtree(tree, color="firebrick", size=2, linetype="dotted")
#By default, the tree is viewed in ladderize form, user can set the parameter ladderize = FALSE to disable it (Figure 4.1C).

ggtree(tree, ladderize=FALSE)
#The branch.length is used to scale the edge, user can set the parameter branch.length = "none" to only view the tree topology (cladogram, Figure 4.1D) or other numerical variable to scale the tree (e.g. dN/dS, see also in Chapter 5).
ggtree(tree, branch.length="none") + geom_tiplab()

## BELOW UPDATE TO PAR COMMAND IS THE TRICK...
par(family = "mono", mar=c(3, 3, 3, 3) + 0.5)
#ggtree(tree) + geom_tiplab(align = TRUE, color='black') + theme_tree2(plot.margin=margin(6, 150, 6, 6))
ggtree(ptree) + geom_tiplab(align = T, color='black', hjust = -0.1, offset = 0.1) + xlim(0, 6) + theme_light(base_family = "Source Sans Pro")



g <- plot(midpoint(ptree), lab4ut="axial",
          edge.width=2, label.offset = 0, cex = 0.75, align.tip.label = TRUE, adj = 1, family = "mono", font = 4)

# original data is in str(g$plot$data)

# but it's easier to process the data for rendering
g[["data"]][[1]][["size"]] <- 5
g[["data"]][[1]][["colour"]] <- "red"


# svg(width = 10, height = 10, family = "mono")
#dev.off()
# plot_dendrogram(tree, mode="phylo", colbar = palette(),
#                 edge.color = NULL, use.edge.length = FALSE)

### idea for in Shiny app - note the exampleData is already stored as a reactive variable: filteredDS
# filteredDS <- reactive({
#   nearPoints(inputdataset(), input$plot_click, threshold = 3, maxpoints = 99999, addDist = FALSE)
# })

# filteredDS <- reactive({
#   if (input$dataset == "Comet Data"    || input$dataset == "Healthy Control")
#   {
#     filteredDS$SEQUENCE_ID <- paste(filteredDS$JUNCTIONAA,filteredDS$GENE,sep="_") 
#     
#   } else {
#     filteredDS$SEQUENCE_ID <- paste(filteredDS$SEQUENCE_ID,filteredDS$JUNCTIONAA,filteredDS$GENE,sep="_") 
#   }
# })



###################
## plot in ape
plot(tree, type = "c")
plot(tree, type = "u")
#plot(tree, type = "u", cex = 2)
plot(tree, "unrooted", main="NJ")
plot(tree)
#plot(tree, type = "r")

plot(unroot(tree),type="unrooted",no.margin=TRUE,lab4ut="axial",
     edge.width=2)

plot(unroot(tree),type="fan",no.margin=TRUE,lab4ut="axial",
     edge.width=2)


## NEED TO REMOVE NEGATIVE BRANCH LENGTHS...
tree$edge.length[tree$edge.length<0] <- 0

plot(tree,ftype="i")
plot(tree)
plot(tree, lab4ut="axial",
     edge.width=2)

### to rerecreate tip names and allow for moving around
## these are better plots...
plot(tree, lab4ut="axial",
     edge.width=2, label.offset = 0.2)


plot(as.phylo(tree), type="fan",
     edge.width=2, label.offset = 0.2)



# ## Creating the vector of labels
# my_labels <- tree$tip.label
# ## "Removing" the unwanted labels (e.g. label 1, 2 and 7)
# #my_labels[c(1,2,5)] <- ""
# ## Adding the labels
# plot(tree, lab4ut="axial",
#      edge.width=2, show.tip.label = FALSE)
# tiplabels(my_labels, adj = 0.5, cex = 1)
# 
# data(bird.orders)
# plot(bird.orders, use.edge.length = FALSE, font = 1)
# nodelabels(bs, adj = -0.2, frame = "n", cex = 0.8, font = 2)
# nodelabels(bs2, adj = c(1.2, 1), frame = "n", cex = 0.8, font = 3)
# nodelabels(bs3, adj = c(1.2, -0.2), frame = "n", cex = 0.8)

# Plotting in ggtree
mytree <- ggtree(phylo_nj, layout="equal_angle", size=0.5, linetype=1)
mytree

mytree <- ggtree(tree, layout="daylight")

##### test in seqinr - don't work 
fastaFile <- readAAStringSet("Twist.fasta")

### read in Rel sequences from the file
#mySeqs <- readAAStringSet("relseqs")   # from package Biostrings
length(fastaFile)

### Turn your alignment into a tree
# convert the alignment for the seqinr package
myAln2 <- msaConvert(fastaFile, type="seqinr::alignment")
# this object is a list object with 4 elements

# generate a distance matrix using seqinr package
d <- dist.alignment(myAln2, "identity")
# From the manual for the seqinr package
# This function computes a matrix of pairwise distances from aligned sequences
# using similarity (Fitch matrix, for protein sequences only) 

# have a look at the output
as.matrix(d)

# generate the tree with the ape package
# the nj() function allows neighbor-joining tree estimation
myTree <- nj(d)


###
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"SEQUENCE_ID"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"JUNCTIONAA"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
require(dplyr)

# exampleData = dplyr::data_frame(name = c("seq1", "seq2", "seq3"),
#                                 seq = c("AAGGTTTTGCCAA", "TTTTGCCAAGGAA", "TTTAAGGTGCCAA"), 
#                                 other = c("meta1", "meta2", "meta3"))

writeFasta(exampleData, "example.fasta")

## better to not have to convert to fasta then convert back!

########### Sep 2020 Dowser commands

#Alternatively, to install the source code, first install the build dependencies:

#  install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))

#To install the latest development code via devtools, along with the development version of Alakazam:

library(devtools)
install_bitbucket("kleinstein/alakazam@master")
install_bitbucket("kleinstein/dowser@master")

# Load required packages
library(alakazam)
install.packages("dowser")
library(dowser)

# load example AIRR tsv data
data(ExampleDb)

# Process example data into proper format
clones <- formatClones(ExampleDb, trait="c_call")

# Build maxmimum parsimony trees for first two clones using 
# phangorn package in R
trees <- getTrees(clones[1:2,])

# simple tree plotting with ggtree R package
plots <- plotTrees(trees, tips="c_call")

# plot tree of first clone
plots[[1]]




#########################

## TO PLOT TOP 15 NETWORKS - IGRAPH CODE HERE....

## note that for all of these commands need non-collapsed datasets...going to be much bigger
BX.full0 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/BXmay10m_stim_germ-pass.tab")
#BX.full0 <- healthy
BX.full0$V_CALL <- as.character(BX.full0$V_CALL)
BX.full0$J_CALL <- as.character(BX.full0$J_CALL)
BX.full0$V_CALL <- as.character(BX.full0$V_CALL)

#BX.full0 <- BX.full0 %>% add_column(CONSCOUNT = 1, DUPCOUNT = 1)
BX.full0 <- BX.full0 %>% add_count(CLONE) %>%
  rename(reads_per_clone = n)
BX.full0 <- subset(BX.full0, CREGION %in% c("IgM", "IgG", "IgA"))

BX.clones.hc <- collapseClones(BX.full0, method="mostCommon", cloneColumn = "CLONE", sequenceColumn = "SEQUENCE_IMGT", germlineColumn = "GERMLINE_IMGT_D_MASK",
                               includeAmbiguous=FALSE, breakTiesStochastic=FALSE)

#BX.fullgandcloneswithmut$CONSCOUNTm <- as.integer(BX.fullgandcloneswithmut$CONSCOUNTm)
BX.full0$CONSCOUNTXn <- (BX.full0$reads_per_clone * BX.full0$CONSCOUNT)
#BX.clusters.hctop15a <- subset(BX.clusters.hc, CREGION %in% c("IgH")) %>% arrange(desc(reads_per_clone, -CONSCOUNTXn)) %>% slice(1:15)

top15 <- BX.full0 %>%
  group_by(CLONE, reads_per_clone) %>%
  summarize(top15 = mean(CONSCOUNTXn)) %>%
  arrange(desc(reads_per_clone),(-top15)) %>%
  ungroup() %>%
  slice(1:15)
BX.clusters.hctop15 <- BX.full0 %>%
  inner_join(top15, by='CLONE') %>%
  select(-reads_per_clone.y) %>%
  select(-top15) %>%
  rename(reads_per_clone = reads_per_clone.x)


## already saved?
#BX.clusters.hctop15 <- read.delim("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/BX.fullgandcloneswithmuttop15.healthycontrol.preprocess.tab")



#BX.sub <- BX.fullgandcloneswithmut
BX.sub <- BX.clusters.hctop15

#BX.sub <- BX.sub %>% filter(CONSCOUNT > 1) ## CAN CHANGE THIS!
BX.sub$sequence_id <- as.character(BX.sub$SEQUENCE_ID)

BX.sub$CREGION <- as.character(BX.sub$CREGION)
BX.sub$SEQUENCE_IMGT <- as.character(BX.sub$SEQUENCE_IMGT)
BX.sub$GERMLINE_IMGT_D_MASK <- as.character(BX.sub$GERMLINE_IMGT_D_MASK)

BX.sub$CONSCOUNT2 <- (2 * log(BX.sub$CONSCOUNT)) + 1
BX.sub$CONSCOUNT2 <- round(BX.sub$CONSCOUNT2)

write.table(BX.sub, "BXsub.tab", sep = "\t", row.names = FALSE, quote = FALSE)

###

clones <- BX.sub %>%
  group_by(CLONE) %>%
  do(CHANGEO=makeChangeoClone(., id = "SEQUENCE_ID",
                              seq = "SEQUENCE_IMGT",
                              germ = "GERMLINE_IMGT_D_MASK",
                              v_call = "V_CALL",
                              j_call = "J_CALL",
                              junc_len = "JUNCTION_LENGTH",
                              clone = "CLONE",max_mask = NULL, pad_end = TRUE, text_fields=c("CREGION"), 
                              num_fields="CONSCOUNT2"))
dnapars_exec <- "~/phylip/exe/dnapars"

graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
                 phylip_exec=dnapars_exec, rm_temp=TRUE)
graphs[sapply(graphs, is.null)] <- NULL
#graphs <- graphs[sapply(graphs, vcount) >= 134]  ## CAN CHANGE THIS OR LEAVE OUT!

###
graph = combineGraphs(graphs)

V(graph)$size = rep(3,length(V(graph)$name))
V(graph)$size = 1.5
V(graph)$size[is.na(V(graph)$CREGION)] = 1
V(graph)$size[grepl("Germline", V(graph)$name)] <- 1
V(graph)$size[V(graph)$CONSCOUNT2 > 5] = 2
V(graph)$size[V(graph)$CONSCOUNT2 > 10] = 3
#V(graph)$size = V(graph)$CONSCOUNT2

V(graph)$color[grepl("Germline", V(graph)$name)] <- "black"
V(graph)$color[V(graph)$CREGION == "Kappa"] <- "steelblue"
V(graph)$color[grepl("IgM", V(graph)$CREGION)] <- "blue"
V(graph)$color[grepl("IgA", V(graph)$CREGION)] <- "green"
V(graph)$color[grepl("IgG", V(graph)$CREGION)] <- "red"

#pal <- brewer.pal(length(unique(V(graph)$CREGION)),"Set1")
#layout = layout.fruchterman.reingold(graph, niter = 5000)
layout = layout.auto
#layout = layout_with_kk(graph, maxiter = 500)

#color = pal[as.numeric(as.factor(vertex_attr(graph, "CREGION")))]
par(mar=c(3, 3, 3, 3) + 0.1)
#rm(BX.full0)

#V(graph)$label <- V(graph)$name
plot.igraph(graph,
            layout=layout,vertex.label=NA, edge.label=NA)	

## check 'graphs' to see how many clonal families there are - best to not have more then 25
## can adjust using graphs <- graphs[sapply(graphs, vcount) >= 15] & changing the number here
## then save plot

layoutisfr <- layout_with_fr(graph, niter = 480)

plot.igraph(graph,
            layout=layoutisfr, circular=TRUE, mode="out", vertex.label=NA, edge.label=NA)	

############
BX.sub.30275 <- subset(BX.sub, CLONE == 30275)

BX.sub.30275f <- filter(BX.sub.30275, CONSCOUNT > 2)



clonesm <- makeChangeoClone(BX.sub.30275f, id = "SEQUENCE_ID",
                            seq = "SEQUENCE_IMGT",
                            germ = "GERMLINE_IMGT_D_MASK",
                            v_call = "V_CALL",
                            j_call = "J_CALL",
                            junc_len = "JUNCTION_LENGTH",
                            clone = "CLONE",max_mask = NULL, pad_end = TRUE, text_fields=c("CREGION"), 
                            num_fields="CONSCOUNT2")
# clonesm <- makeChangeoClone(BX.sub.30275, max_mask = NULL, pad_end = TRUE, text_fields=c("DONOR","SAMPLE","CREGION"), 
#                            num_fields="CONSCOUNT2")
clonesm@data[, c("CONSCOUNT2")]

# clone2@data[, c("SAMPLE", "DONOR","CONSCOUNT2")]
# clone3@data[, c("SEQUENCE_ID","CREGION", "CONSCOUNT2","SAMPLE","DONOR")]


dnapars_exec <- "~/phylip/exe/dnapars"
#graph <- buildPhylipLineage(clone2, dnapars_exec, dist_mat=getDNAMatrix(gap=-1), rm_temp=TRUE)
graphsm <- buildPhylipLineage(clonesm, dnapars_exec, rm_temp=TRUE)


V(graphsm)$shape <- "circle"
V(graphsm)$shape[V(graphsm)$CREGION == "IgG"] <- "circle"
V(graphsm)$shape[V(graphsm)$CREGION == "IgA"] <- "square"
#V(graphsm)$shape[V(graphsm)$CREGION == "IgA"] <- "triangle"
V(graphsm)$shape[V(graphsm)$CREGION == "IgG,IgM"] <- "circle"
V(graphsm)$shape[V(graphsm)$CREGION == "IgA,IgG"] <- "circle"
V(graphsm)$shape[V(graphsm)$CREGION == "IgA,IgG,IgM"] <- "circle"

V(graphsm)$color[grepl("Germline", V(graphsm)$name)] <- "black"
V(graphsm)$color[grepl("IgM", V(graphsm)$CREGION)] <- "blue"
V(graphsm)$color[grepl("IgA", V(graphsm)$CREGION)] <- "green"
V(graphsm)$color[grepl("IgG", V(graphsm)$CREGION)] <- "red"


V(graphsm)$color[V(graphsm)$name == "Germline"] <- "black"
V(graphsm)$size <- 2
V(graphsm)$size[V(graphsm)$name == "Germline"] <- 2
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 3] <- 3
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 4.4] <- 4
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 5] <- 5
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 6.4] <- 6
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 7] <- 7
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 8.4] <- 9
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 8.8] <- 11

V(graphsm)$size[grepl("Inferred", V(graphsm)$name)] <- 1
V(graphsm)$size[grepl("INF", V(graphsm)$name)] <- 2
V(graphsm)$color[V(graphsm)$name == "Germline"] <- "black"

### plot with CREGION
V(graphsm)$label <- V(graphsm)$CREGION

clonetree_PCRONSnames <- plot(graphsm, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, vertex.frame.color="grey",
                              vertex.label.color="black", edge.label.color="black")

### plot with names only
V(graphsm)$name[grepl("Inferred", V(graphsm)$name)] <- ""
#V(graphsm)$name[grepl("o", V(graphsm)$name)] <- ""

V(graphsm)$label <- V(graphsm)$name
clonetree_names <- plot(graphsm, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, label.cex=0.1, vertex.frame.color="grey",
                        vertex.label.color="black", edge.label.color="black")

#V(graphsm)$name[grepl("132s", V(graphsm)$name)] <- "s"
#clonetree_names <- plot(graphsm, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, label.cex=0.1, vertex.frame.color="grey",
#                        vertex.label.color="black", edge.label.color="black")


### plot with no labels
V(graphsm)$label <- ""
#V(graphsm)$label <- V(graphsm)$name
clonetree_nolabels <- plot(graphsm, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, vertex.frame.color="grey",
                           vertex.label.color="black", edge.label.color="black")

### bigger trees
V(graphsm)$size <- 25
V(graphsm)$size[V(graphsm)$name == "Germline"] <- 20
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 3] <- 30
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 4.4] <- 40
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 5] <- 50
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 6.4] <- 60
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 7] <- 70
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 8.4] <- 90
V(graphsm)$size[V(graphsm)$CONSCOUNT2 > 8.8] <- 110

onetree_nolabels <- plot(graphsm, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(0,14),xlim=c(-8,14), asp = 1.4)

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### jan 2021 dengue networks

### jan 2021 some quick networks of cf10 and cf7 for latest paper version (moved to bottom of repertoire_vis)
# CF10mergedwithOAS <- read.delim("toshiny_dengue_cf10.tab")
# CF7mergedwithOAS <- read.delim("toshiny_dengue_cf4andcf7.tab")
## but already loaded as
toshiny.dengue.cf10
toshiny.dengue.cf7

## BUT NEED LARGER STARTING FILE WITH GERMLINE MASK & SEQUENCE IMGT

## NEED TO THIN TO SMALLER SUBSET THOUGH...inner join with subset from Geneious
##  dengue_CF10_list.tsv
CF10mergedwithOAS.list <- read.delim("dengue_CF10_list.tsv")

## next run some sort of join...semi-join semi_join() return all rows from x with a match in y.

CF10mergedwithOASfewer <- semi_join(toshiny.dengue.cf10, CF10mergedwithOAS.list)


BX.sub <- CF10mergedwithOASfewer

BX.sub <- CF4mergedwithOASfewer


########################################################################################################################
########################################################################################################################

### WILL NEED TO RENAME SEQUENCE_ID, SEQUENCE_IMGT, AND D_MASK
BX.sub$SEQUENCE_IMGT <- as.character(BX.sub$SEQUENCE_IMGT)
BX.sub$GERMLINE_IMGT_D_MASK <- as.character(BX.sub$GERMLINE_IMGT_D_MASK)
BX.sub$SEQUENCE_ID <- as.character(BX.sub$SEQUENCE_ID)

BX.sub <- rename(BX.sub, sequence_id = SEQUENCE_ID)
BX.sub <- rename(BX.sub, sequence_alignment = SEQUENCE_IMGT)
BX.sub <- rename(BX.sub, germline_alignment = GERMLINE_IMGT_D_MASK)
BX.sub <- rename(BX.sub, v_call = V_CALL)
BX.sub <- rename(BX.sub, j_call = J_CALL)
BX.sub <- rename(BX.sub, junction_length = JUNCTION_LENGTH)

# BX.sub$sequence_id <- BX.sub$SEQUENCE_ID
# BX.sub$germline_alignment <- BX.sub$SEQUENCE_IMGT
# BX.sub$germline_alignment_d_mask <- BX.sub$GERMLINE_IMGT_D_MASK



BX.sub$CREGION0 <- as.character(BX.sub$CREGION0)
# BX.sub$SAMPLE <- as.character(BX.sub$SAMPLE)
# BX.sub$DONOR <- as.character(BX.sub$DONOR)

BX.sub$CONSCOUNT2 <- (2 * log(BX.sub$CONSCOUNT)) + 1
BX.sub$CONSCOUNT2 <- round(BX.sub$CONSCOUNT2)


BX.sub$clone_id <- 4

#clone2 <- makeChangeoClone(BX.sub, max_mask = NULL, pad_end = TRUE, text_fields=c("CREGION0"), 
#                           num_fields="CONSCOUNT2")
clone2 <- makeChangeoClone(BX.sub, max_mask = NULL, pad_end = TRUE, text_fields=c("CREGION0"), 
                           num_fields="CONSCOUNT2")
clone3 <- makeChangeoClone(BX.sub, max_mask = NULL, pad_end = TRUE, text_fields=c("sequence_id","CREGION0"), 
                           num_fields="CONSCOUNT2")

#BX.sub$SEQUENCE_ID <- padSeqEnds(BX.sub$SEQUENCE_ID)
#BX.sub$SEQUENCE_IMGT <- padSeqEnds(BX.sub$SEQUENCE_IMGT)
#BX.sub$GERMLINE_IMGT_D_MASK <- padSeqEnds(BX.sub$GERMLINE_IMGT_D_MASK)

clone2@data[, c("CONSCOUNT2")]
clone3@data[, c("sequence_id","CREGION0", "CONSCOUNT2")]


#clone3@data[, c("SEQUENCE_ID","CREGION0", "CONSCOUNT2")]

dnapars_exec <- "~/phylip/exe/dnapars"
#graph <- buildPhylipLineage(clone2, dnapars_exec, dist_mat=getDNAMatrix(gap=-1), rm_temp=TRUE)
graph <- buildPhylipLineage(clone2, dnapars_exec, rm_temp=TRUE)


## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# The graph has shared annotations for the clone
data.frame(CLONE=graph$clone,
           JUNCTION_LENGTH=graph$junc_len,
           V_GENE=graph$v_gene,
           J_GENE=graph$j_gene)

# The vertices have sequence specific annotations
data.frame(sequence_id=V(graph)$name, 
           CREGION0=V(graph)$CREGION0,
           CONSCOUNT2=V(graph)$CONSCOUNT2)

## ---- eval=TRUE----------------------------------------------------------
# Plot graph with defaults
plot(graph)

########################################


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("FELLA")
# library("FELLA")
# 
# add.vertex.shape(
#   "triangle", clip = shapes("circle")$clip,
#   plot = FELLA:::mytriangle)


########################################
########################################
## USED FOR Feb 2020 plots
########################################

## NOTE SOME METHODS: FOR THESE TREES I TOOK THE 40 CLONAL MEMBERS WITH THE HIGHEST UMI COUNT, AND ALSO ANY SEQUENCES IDENTICAL OR NEARLY IDENTICAL TO MABS (0-1BP)
V(graph)$shape <- "circle"
V(graph)$shape[V(graph)$CREGION0 == "IgA"] <- "square"
V(graph)$shape[V(graph)$CREGION0 == "IgM"] <- "rectangle"
V(graph)$shape[V(graph)$CREGION0 == "IgG,IgM"] <- "rectangle"


#grepl("bipartite|merge|norm|sugiyama|tree",

V(graph)$color[grepl("OAS", V(graph)$name)] <- "yellow"

### 3 colors only
V(graph)$color[grepl("d13", V(graph)$name)] <- "orange"
V(graph)$color[grepl("d20", V(graph)$name)] <- "blue"
V(graph)$color[grepl("J8", V(graph)$name)] <- "gray"
V(graph)$color[grepl("J9", V(graph)$name)] <- "gray"
#V(graph)$color[grepl("xxx", V(graph)$name)] <- "red"
#V(graph)$color[grepl("d13,d20", V(graph)$name)] <- "red"
V(graph)$color[V(graph)$name == "d13,d20"] <- "red"
#V(graph)$shape[V(graph)$CREGION0 == "IgG,IgM"] <- "circle"
V(graph)$color[grepl("mab", V(graph)$name)] <- "purple"


### many colors
V(graph)$color[grepl("d13", V(graph)$name)] <- "orange"
V(graph)$color[grepl("d13_3e", V(graph)$name)] <- "palegreen"
V(graph)$color[grepl("d13_3g", V(graph)$name)] <- "green"
#V(graph)$color[grepl("d13_2yg", V(graph)$name)] <- "green"
V(graph)$color[grepl("d13_3s", V(graph)$name)] <- "darkgreen"
V(graph)$color[grepl("d13_2e", V(graph)$name)] <- "beige"
V(graph)$color[grepl("d13_2g", V(graph)$name)] <- "tan"
V(graph)$color[grepl("d13_2s", V(graph)$name)] <- "orange"
#V(graph)$color[grepl("d13_2xg", V(graph)$name)] <- "orange"
V(graph)$color[grepl("d20", V(graph)$name)] <- "blue"

V(graph)$color[grepl("d20_3e", V(graph)$name)] <- "plum"
V(graph)$color[grepl("d20_3g", V(graph)$name)] <- "mediumorchid2"
V(graph)$color[grepl("d20_3s", V(graph)$name)] <- "magenta4"
V(graph)$color[grepl("d20_2e", V(graph)$name)] <- "lightblue"
V(graph)$color[grepl("d20_2s", V(graph)$name)] <- "dodgerblue"
V(graph)$color[grepl("d20_2g", V(graph)$name)] <- "mediumblue"
V(graph)$color[grepl("J8", V(graph)$name)] <- "gray"
V(graph)$color[grepl("J9", V(graph)$name)] <- "gray"
#V(graph)$color[grepl("xxx", V(graph)$name)] <- "red"
V(graph)$color[V(graph)$name == "d13,d20"] <- "red"
V(graph)$color[grepl("mab", V(graph)$name)] <- "gray"
V(graph)$color[V(graph)$name == "Germline"] <- "black"

#V(graph)$color[V(graph)$name == "d20_2g_ACCAAATATCTGTGTT"] <- "slateblue" changing upper right = B10 so middle gray circle is D8, small left is M1
#grepl("bipartite|merge|norm|sugiyama|tree",

## if using original dimensions
#V(graph)$size[V(graph)$name == "Germline"] <- 3
#V(graph)$size[grepl("Inferred", V(graph)$name)] <- 1
#V(graph)$size[grepl("INF", V(graph)$name)] <- 3 
#V(graph)$color[V(graph)$name == "Germline"] <- "black"

# V(graph)$shape <- "circle"
# V(graph)$shape[V(graph)$CREGION0 == "IgA"] <- "square"


V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$size <- 2
V(graph)$size[grepl("Inferred", V(graph)$name)] <- 1
V(graph)$size[grepl("Germline", V(graph)$name)] <- 3
V(graph)$size[V(graph)$CONSCOUNT2 > 3] <- 3
V(graph)$size[grepl("OAS", V(graph)$name)] <- 3
V(graph)$size[grepl("mab", V(graph)$name)] <- 6
V(graph)$size[V(graph)$CONSCOUNT2 > 4.4] <- 4
V(graph)$size[V(graph)$CONSCOUNT2 > 5] <- 5
V(graph)$size[V(graph)$CONSCOUNT2 > 6.4] <- 6
V(graph)$size[V(graph)$CONSCOUNT2 > 7] <- 7
V(graph)$size[V(graph)$CONSCOUNT2 > 8.4] <- 9
V(graph)$size[V(graph)$CONSCOUNT2 > 8.8] <- 11

## to see names
V(graph)$name
V(graph)$CREGION0
V(graph)$CONSCOUNT2

################################################################################################
################################################################################################
# Remove large default margins


par(mar=c(0, 0, 0, 0) + 0.1)
# Plot graph
#layout <- layout.reingold.tilford(graph, circular=F)

### plot with CREGION0
V(graph)$label <- V(graph)$CREGION0

clonetree_PCRONSnames <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, vertex.frame.color="grey",
                              vertex.label.color="black", edge.label.color="black")

### plot with names only
V(graph)$name[grepl("Inferred", V(graph)$name)] <- ""
#V(graph)$name[grepl("o", V(graph)$name)] <- ""

V(graph)$label <- V(graph)$name
clonetree_names <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, label.cex=0.1, vertex.frame.color="grey",
                        vertex.label.color="black", edge.label.color="black")

#V(graph)$name[grepl("132s", V(graph)$name)] <- "s"
#clonetree_names <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, label.cex=0.1, vertex.frame.color="grey",
#                        vertex.label.color="black", edge.label.color="black")


### plot with no labels
V(graph)$label <- ""
#V(graph)$label <- V(graph)$name
clonetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=1, edge.label.cex=1, vertex.frame.color="grey",
                           vertex.label.color="black", edge.label.color="black")

### 12/3/19 HERE TRY TO MAKE BRANCHES DARKER...
#V(graph)$size[V(graph)$name == "Germline"] <- 5
#V(graph)$size[grepl("o", V(graph)$name)] <- 2
#clonetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=3, edge.color="black", edge.label=NA, vertex.frame.color="gray",
#                           vertex.label.color="black", edge.label.color="black")


### saved as Feb20_j8j9lineage_noOAS_1200x1200_3colors or Feb20_j8j9lineage_noOAS_1200x1200_manycolors

### EXPANDING THE TREE SIZE!!! - NOTE THAT I HAD TO INCREASE ALL SIZES BY A FACTOR OF 10
V(graph)$label <- ""
V(graph)$size[V(graph)$name == "Germline"] <- 50
#V(graph)$size[grepl("o", V(graph)$name)] <- 2
#clonetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="gray",
#                           vertex.label.color="black", edge.label.color="black")

## for trees with OAS
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(0,14),xlim=c(-8,14), asp = 1.4)
legend("bottomright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab","OAS"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray","yellow"))

## for J8J9_v2
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(0,12),xlim=c(4,8), asp = 2.2)
legend("bottomright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab","OAS"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray","yellow"))

## for J8J9_v3
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(0,22),xlim=c(4,8), asp = 1.8)
legend("bottomright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab","OAS"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray","yellow"))

## for trees without OAS
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(0,14),xlim=c(-14,8), asp = 1.4)

legend("bottomright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray"))

## M1 trees with OAS
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(2,28),xlim=c(-20,108), asp = 1.4)
#legend("topright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab","OAS"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray","yellow"))

## M1 trees without OAS
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=3, edge.label.cex=1, label.dist =1, edge.label.dist="1", vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(0,15),xlim=c(-34,42), asp = 0)

legend("bottomright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray"))

#onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=3, edge.label.cex=1, label.dist =1, edge.label.dist="1", vertex.frame.color="gray",
#                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(5,15),xlim=c(-34,34), asp = 0.95)



## M1 trees with OAS _v33
onetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="gray",
                         vertex.label.color="black", edge.label.color="black", rescale = FALSE, ylim=c(2,22),xlim=c(-30,50), asp = 1.4)
#legend("topright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab","OAS"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray","yellow"))

## for old 1200x1200 size graphs
V(graph)$size <- 2
V(graph)$size[V(graph)$name == "Germline"] <- 5
V(graph)$size[V(graph)$name == "Germline"] <- 3
V(graph)$size[grepl("Inferred", V(graph)$name)] <- 1
V(graph)$size[grepl("INF", V(graph)$name)] <- 3 
V(graph)$color[V(graph)$name == "Germline"] <- "black"


V(graph)$size[V(graph)$CONSCOUNT2 > 3] <- 3
V(graph)$size[V(graph)$CONSCOUNT2 > 4.4] <- 4
V(graph)$size[V(graph)$CONSCOUNT2 > 5] <- 5
V(graph)$size[V(graph)$CONSCOUNT2 > 6.4] <- 6
V(graph)$size[V(graph)$CONSCOUNT2 > 7] <- 7
V(graph)$size[V(graph)$CONSCOUNT2 > 8.4] <- 9
V(graph)$size[V(graph)$CONSCOUNT2 > 8.8] <- 11
clonetree_nolabels <- plot(graph, layout=layout_as_tree, mode="out", edge.arrow.mode=0, edge.width=2, edge.label.cex=1, vertex.frame.color="grey",
                           vertex.label.color="black", edge.label.color="black")
legend("bottomright",legend=c("d13_2","d13_2g","d13_2s","d13_3","d13_3g","d13_3s","d20_2","d20_2g","d20_2s","d20_3","d20_3g","d20_3s","d13andd20","mab"),fill=c("beige","tan","orange","palegreen","green","darkgreen","lightblue","dodgerblue","mediumblue","plum","mediumorchid2","magenta4","red","gray"))


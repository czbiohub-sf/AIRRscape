
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
# library(limma)
library(reshape2)
library(DT)

library(ggplot2)
library(seqinr)
library(phangorn)
library(shiny)
# install.packages("plotly")
library(plotly)
library(rlang)
library(tidyverse)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}


#dir.create("~/code")
## this works ONLY IF YOU HAVE CREATED THE DIRECTORY FIRST...
# setwd("/Users/eric.waltari/data_carpentry/AIRRScape")


##################################################################################################################
### COMMANDS FOR READING IN AIRR-COMPLIANT DATASETS FOR SHINY APP
### note that for non-AIRR (e.g. older Immcantation data or databases like AbDab) datasets, column names will be different
######### solution is to use only AIRR datasets, fixed as of Sept 2021
###### particularly note that 1) junction MAY NOT NEED TO BE TRIMMED 1 AA ON EACH END TO BE KABAT CDR3
######                        2) v_identity MAY BE BETWEEN 0-1 NOT 0-100, AND THUS shm FORMULA NEEDS TO CHANGE
### updated function below now addresses both points 1 & 2
##################################################################################################################

### SEPT 2021 DATASETS

toshiny.cov2.abdab <- read_tsv("toshiny_cov2abdab.tab")

cov2.bulk.binder.p11 <- read_tsv("Binder_p11_germ-pass.tsv")   # v_identity between 0.60 and 1, sequences renamed, have clone_id & germline_alignment_d_mask, NOT cdr3_aa
cov2.bulk.nielsen.p7450 <- read_tsv("Nielsen_7450_airr-covid-19.tsv")   # v_identity between 60 and 100, has cdr3_aa
cov2.bulk.galson.p1 <- read_tsv("Galson_p1_germ-pass.tsv")   # v_identity between 0.60 and 1, sequences renamed, have clone_id & germline_alignment_d_mask, NOT cdr3_aa
#cov2.bulk.kc.m5.rep1 <- read_tsv("Kuri-Cervantes_M5-rep1_airr-covid-19.tsv")   # v_identity between 0.60 and 1, NOT cdr3_aa
cov2.bulk.kc.m5.allreps <- read_tsv("Kuri-Cervantes_M5-reps1to4_airr-covid-19.tsv")   # v_identity between 0.60 and 1, NOT cdr3_aa

toshiny.hiv.mabs <- read_tsv("toshiny_hivmabs.tab")
hiv.bulk.nih45 <- read_tsv("hivnih45_vdjserver.tsv")   # v_identity between 
hiv.bulk.mt1214 <- read_tsv("MT1214downloaded.tab")   # v_identity between 

toshiny.den.mabs <- read_tsv("toshiny_denmabs.tab")

den.bulk.OAS <- read_tsv("OAS_sept21_germ-pass2.tsv")
den.bulk.d13enrich <- read_tsv("d13_2enrichHC_germ-pass.tsv")
den.bulk.d13stim <- read_tsv("d13_2stimHC_germ-pass.tsv")
hc.BXmay.10mstim <- read_tsv("BXmay10mstimHC_germ-pass.tsv")

## one-off rennaming of sequence_id in OAS dataset...
den.bulk.OAS$sequence_id <- gsub("\\_","\\-",den.bulk.OAS$sequence_id)
den.bulk.OAS$sequence_id <- gsub("2014","2013",den.bulk.OAS$sequence_id)

head(den.bulk.OAS)
head(den.bulk.d13enrich)
head(den.bulk.d13stim)

head(hc.BXmay.10mstim)

head(hiv.bulk.nih45)
head(hiv.bulk.mt1214)

head(cov2.bulk.binder.p11)
head(cov2.bulk.nielsen.p7450)
head(cov2.bulk.galson.p1)
head(cov2.bulk.kc.m5.rep1)

### AIRR-COMPLIANT TSV FILES NEED THE FOLLOWING MODIFICATIONS
# CHANGING RELEVANT HEADERS FROM LOWERCASE TO UPPER CASE
# REMOVING FIRST AND LAST RESIDUES FROM junction
# CONVERT v_identity TO shm - FOR SHINY WANT 0-100

## mutate_at has been superseded by across()
## also funs() is deprecated, use a list
#df %>% mutate(across(cols, round, 3)) 
## updating this command

## all of these
#mutate_at(vars(shm_mean), funs(round(., 2)))
## are now
#mutate(across(shm_mean, round, 2)) 


## check to see if there is cdr3_aa column, if not trim start and end from junctionaa
### make sure length of cdr3_aa is cdr3length_imgt
## USE THESE NAMES: cdr3_aa_imgt & cdr3length_imgt

### can just go with some universal column creating...i.e. always go with trim of junction_aa?

## also always remove any stop_codon = true, productive = false, any * in junction_aa

## create a function to do all pre-processing all at once???
## also need to rename sequence_id's??
## required columns: v_identity, v_call, j_call

## REMEMBER - FOR SHINY APP PROCESSING THERE CAN BE NO UNDERSCORES IN THE SEQUENCE_ID NAMES!!!
## USE DASHES INSTEAD OF UNDERSCORES!!!

### THIS IS NOW BETTER THAN PREVIOUS
shinyprocess <- function(x, filter_columns = TRUE, renumber_sequences = TRUE, filter_after_counting = TRUE) {
  colname <- substitute(x)
  ## this removes columns with all NAs
  x <- x[!map_lgl(x, ~ all(is.na(.)))]
  ## this calculates SHM but depending on whether v_identity is from 0 to 1 or 0 to 100
  if (mean(x$v_identity) < 1) {
    x$shm <- (100 - (x$v_identity * 100))
  } else {
    x$shm <- (100 - x$v_identity)
  }
  ## this makes new standard cdr3 column (sometimes already exists, but there should always be a junction_aa) by removing both ends of the junction_aa column
  x$cdr3_aa_imgt <- x$junction_aa
  str_sub(x$cdr3_aa_imgt, -1, -1) <- ""
  str_sub(x$cdr3_aa_imgt, 1, 1) <- ""
  ## this calculates the CDR3 length
  x$cdr3length_imgt <- nchar(x$cdr3_aa_imgt)
  ## removing non-productive, out of frame, stop codons, 4Ns
  x <- x %>% filter(productive != "FALSE") %>%
    filter(vj_in_frame != "FALSE") %>%
    filter(productive != "F") %>%
    filter(vj_in_frame != "F")
  x <- x[ grep("\\*", x$junction_aa, invert = TRUE) , ]
  x <- x[ grep("\\X", x$junction_aa, invert = TRUE) , ]  ## USE ONLY FOR KC DATASET - OCT21 USING FOR ALL
  # x <- x[ grep("NNNN", x$sequence, invert = TRUE) , ]   # ...NO LONGER USING
  ### removing all sequences with IMGT CDR3 less than 3
  x <- x %>% filter(cdr3length_imgt > 2.8)  
  ## next lines create V gene family, J gene columns
  x$gene <- getGene(x$v_call, first=TRUE, strip_d=TRUE)
  x$gf <- substring(x$gene, 1,5)
  x$jgene <- getGene(x$j_call, first=TRUE, strip_d=TRUE)
  ## this creates new column gf_jgene which is used in all shiny plots
  x <- x %>% unite(gf_jgene, gf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)
  ## this removes any rows without CDR3, or with junctions that are not 3-mers
  x <- x %>% filter(!is.na(cdr3length_imgt)) %>% 
    filter(is.wholenumber(cdr3length_imgt))
  # if there is a clone_id column this will make a count of reads_per_clone
  if ("clone_id" %in% names(x)) {
    x <- x %>% add_count(clone_id) %>%
      rename(reads_per_clone = n)
  }
  ## if no cregion column, make one  !(x %in% y)
  if (!("cregion" %in% names(x))) {
    x$cregion <- str_sub(x$v_call, end=3)
    x$cregion <- gsub('IG','Ig',x$cregion)
  }
  ## making more important columns used in plotting, also a rounding step
  x <- x %>%
    add_count(gf_jgene,cdr3length_imgt) %>% 
    rename(ncount = n) %>% 
    group_by(gf_jgene,cdr3length_imgt) %>% 
    mutate(shm_mean = mean(shm, na.rm = TRUE)) %>% 
    # NOTE ADDIN MAX SHM AS WELL..
    mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
    mutate(across(shm, round, 2)) %>% 
    mutate(across(shm_max, round, 2)) %>% 
    mutate(across(shm_mean, round, 2))
## this will filter the dataset if filter_columns option is set to true - note the any_of which allows columns to be missing
  vars2 <- c("sequence_id", "binding", "neutralization", "cregion", "cdr3_aa_imgt","gene", "gf_jgene", "gf","jgene", "cdr3length_imgt", "shm", "shm_max", "shm_mean", "ncount", "reads_per_clone")
  if (filter_columns) {
    x <- x %>% select(any_of(vars2))
  }
  ## this will remove all redundant sequences with same gf/gene & cdr3 motif...note we count above so okay to collapse here!!
  if (filter_after_counting) {
    x <- x %>%
      group_by(cdr3_aa_imgt,gf_jgene) %>%
      summarize_all(first) %>%
      rename(ncountfull = ncount) %>% 
      ungroup() %>%
      add_count(gf_jgene,cdr3length_imgt) %>% 
      rename(ncount = n) %>%
      relocate(ncount, .before = shm_mean)
  }
  ## this will make a new sequence_id column with new row names if renumber_sequences option is set to true MOVING LAST TO CHANGE X DEFINITION
  if (renumber_sequences) {
    ## NOTE THIS NOW WORKS, TRICK WAS TO ASSIGN COLNAME VERY EARLY ON BEFORE ANYTHING ELSE...
    ## oct 3 adding change from underscores to dashes...
    x$dataset <- deparse(substitute(colname))
    x$dataset <- gsub('"','',x$dataset)
    x$dataset <- gsub("\\.","\\-",x$dataset)
    x$dataset <- gsub("\\_","\\-",x$dataset)
    x$obs <- 1:nrow(x) 
    x <- x %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
    x <- x %>% relocate(sequence_id, .before = cregion)
    # x$sequence_id <- gsub("\\_","\\-",x$sequence_id) ## moved to always run
  }
  ## need to always check and remove underscores from all names
  x$sequence_id <- gsub("\\_","\\-",x$sequence_id)
  return(x)
}


## NOTES
## latest looks like quo_name but now as_name is the way to do this...  https://community.rstudio.com/t/how-to-extract-the-object-name-when-the-object-is-piped-to-a-function/99554
# variable_name <- rlang::as_name(rlang::enquo(variable))
# x$dataset <- rlang::as_name(rlang::enquo(x))
# these still don't work

## adding second argument with the name??
# nameprinter <- function(x, onename) {
# onename <- deparse(substitute(onename))
# When you have the data-variable in a function argument (i.e. an env-variable that holds a promise2),
# you need to embrace the argument by surrounding it in doubled braces, like filter(df, {{ var }}).
## that doesn't work but maybe enquo()   https://stackoverflow.com/questions/49700912/why-is-enquo-preferable-to-substitute-eval
# UPDATE: rlang 0.4.0 introduced a new operator {{ (pronounced "curly curly"), which is effectively a short hand for !!enquo(). This allows us to simplify the definition of g2 to
#   g2 <- function( myExpr ) {
#     val <- f2( {{myExpr}} )
#     val
#   }
# 
# val <- f2( !!enquo(myExpr) )

#debugonce(shinyprocess)

## works
test.den.bulk.d13enrich <- shinyprocess(den.bulk.d13enrich, renumber_sequences = FALSE, filter_columns = FALSE)
test.den.bulk.d13enrich <- shinyprocess(den.bulk.d13enrich, renumber_sequences = FALSE)

## NOW IT WORKS WITH THE SINGLE ARGUMENT!!
test.den.bulk.d13enrich <- shinyprocess(den.bulk.d13enrich)


#undebug(shinyprocess)

den.bulk.OAS$sequence_id[1]          #name ok
den.bulk.d13enrich$sequence_id[1]
den.bulk.d13stim$sequence_id[1]

hiv.bulk.nih45$sequence_id[1]
hiv.bulk.mt1214$sequence_id[1]       #name ok but change anyway
cov2.bulk.binder.p11$sequence_id[1]       #name ok but change anyway
cov2.bulk.nielsen.p7450$sequence_id[1]
cov2.bulk.galson.p1$sequence_id[1]    #name ok but change anyway
cov2.bulk.kc.m5.rep1$sequence_id[1]
hc.BXmay.10mstim$sequence_id[1]


toshiny.den.bulk.d13enrich <- shinyprocess(den.bulk.d13enrich)
toshiny.den.bulk.d13stim <- shinyprocess(den.bulk.d13stim)
toshiny.den.bulk.OAS <- shinyprocess(den.bulk.OAS, renumber_sequences = FALSE)

toshiny.hiv.bulk.mt1214 <- shinyprocess(hiv.bulk.mt1214)
toshiny.hiv.bulk.nih45 <- shinyprocess(hiv.bulk.nih45)
toshiny.hc.BXmay.10mstim <- shinyprocess(hc.BXmay.10mstim)

toshiny.cov2.bulk.binder.p11 <- shinyprocess(cov2.bulk.binder.p11)
toshiny.cov2.bulk.galson.p1 <- shinyprocess(cov2.bulk.galson.p1)
toshiny.cov2.bulk.nielsen.p7450 <- shinyprocess(cov2.bulk.nielsen.p7450)

toshiny.cov2.bulk.kc.m5.allreps <- shinyprocess(cov2.bulk.kc.m5.allreps)  ## need to run without the NNNN removing...will remove all!! ALSO INSTEAD REMOVE X'S


#######################################
#######################################
## one off commands...
toshiny.hiv.bulk.mt1214$cregion <- toshiny.hiv.bulk.mt1214$cregion %>% replace_na("IgH")
toshiny.hiv.bulk.nih45 <- toshiny.hiv.bulk.nih45[ grep("IgK", toshiny.hiv.bulk.nih45$cregion, invert = TRUE) , ]



# toshiny.den.bulk.OAS$cregion[1]          #caps
# toshiny.den.bulk.d13enrich$cregion[1]
# toshiny.den.bulk.d13stim$cregion[1]
# 
# toshiny.hiv.bulk.nih45$cregion[1]        #caps
# toshiny.hiv.bulk.mt1214$cregion[1]  
# toshiny.cov2.bulk.binder.p11$cregion[1]          #caps
# toshiny.cov2.bulk.nielsen.p7450$cregion[1]        #caps
# toshiny.cov2.bulk.galson.p1$cregion[1]          #caps
# toshiny.cov2.bulk.kc.m5.allreps$cregion[1]        #caps
# toshiny.hc.BXmay.10mstim$cregion[1]

### checking datasets...
unique(toshiny.den.bulk.d13stim$cregion)
unique(toshiny.den.bulk.d13enrich$cregion)
unique(toshiny.den.bulk.OAS$cregion)

unique(toshiny.hiv.bulk.nih45$cregion)  ## ONE IGK!!
unique(toshiny.hiv.bulk.mt1214$cregion) ## NA

unique(toshiny.hc.BXmay.10mstim$cregion)
unique(toshiny.cov2.bulk.binder.p11$cregion)
unique(toshiny.cov2.bulk.nielsen.p7450$cregion)
unique(toshiny.cov2.bulk.galson.p1$cregion)
unique(toshiny.cov2.bulk.kc.m5.allreps$cregion)

unique(toshiny.cov2.abdab$cregion)
unique(toshiny.den.mabs$cregion)
unique(toshiny.hiv.mabs$cregion)

unique(toshiny.den.bulk.d13stim$gf)
unique(toshiny.den.bulk.d13enrich$gf)
unique(toshiny.den.bulk.OAS$gf)
unique(toshiny.hiv.bulk.nih45$gf)  ## only HV1, HV3, HV4
unique(toshiny.hiv.bulk.mt1214$gf) ## NA
unique(toshiny.hc.BXmay.10mstim$gf)
unique(toshiny.cov2.bulk.binder.p11$gf)
unique(toshiny.cov2.bulk.nielsen.p7450$gf)
unique(toshiny.cov2.bulk.galson.p1$gf)
unique(toshiny.cov2.bulk.kc.m5.allreps$gf)


unique(toshiny.cov2.abdab$jgene)
unique(toshiny.den.mabs$jgene)
unique(toshiny.hiv.mabs$jgene)

unique(toshiny.den.bulk.d13stim$jgene)
unique(toshiny.den.bulk.d13enrich$jgene)
unique(toshiny.den.bulk.OAS$jgene)
unique(toshiny.hiv.bulk.nih45$jgene)  ## only HV1, HV3, HV4
unique(toshiny.hiv.bulk.mt1214$jgene) ## NA
unique(toshiny.hc.BXmay.10mstim$jgene)
unique(toshiny.cov2.bulk.binder.p11$jgene)
unique(toshiny.cov2.bulk.nielsen.p7450$jgene)
unique(toshiny.cov2.bulk.galson.p1$jgene)
unique(toshiny.cov2.bulk.kc.m5.allreps$jgene)
unique(toshiny.cov2.abdab$jgene)
unique(toshiny.den.mabs$jgene)
unique(toshiny.hiv.mabs$jgene)




## ALSO MORE ONE-OFF COMMANDS FOR MAB DATASETS
## need to recalculate shm_mean & shm_max (new) for the three mab datasets
# toshiny.den.mabs
# toshiny.hiv.mabs
# toshiny.cov2.abdab

## also need to rename sequence_id for cov2 and hiv mabs
# SARS-CoV2-mAb
# HIV-mAb
## BUT ALSO WILL WANT TO CHANGE THE NAME OF EACH SEQUENCE!!! DONE HERE...
toshiny.hiv.mabs$sequence_id0 <- toshiny.hiv.mabs$sequence_id
toshiny.hiv.mabs$sequence_id <- NULL
toshiny.hiv.mabs$dataset <- "HIV-mAb"
toshiny.hiv.mabs <- toshiny.hiv.mabs %>% unite(sequence_id, dataset, sequence_id0, sep = "-", remove = TRUE, na.rm = TRUE)

toshiny.cov2.abdab$sequence_id0 <- toshiny.cov2.abdab$sequence_id
toshiny.cov2.abdab$sequence_id <- NULL
toshiny.cov2.abdab$dataset <- "SARS-CoV2-mAb"
toshiny.cov2.abdab <- toshiny.cov2.abdab %>% unite(sequence_id, dataset, sequence_id0, sep = "-", remove = TRUE, na.rm = TRUE)

## first part changed very early 
toshiny.cov2.abdab$cdr3_aa_imgt <- toshiny.cov2.abdab$cdr3_aa
toshiny.cov2.abdab$cdr3_aa <- NULL

toshiny.hiv.mabs$cdr3_aa_imgt <- toshiny.hiv.mabs$cdr3_aa
toshiny.hiv.mabs$cdr3_aa <- NULL

toshiny.den.mabs$cdr3_aa_imgt <- toshiny.den.mabs$cdr3_aa
toshiny.den.mabs$cdr3_aa <- NULL

toshiny.hiv.mabs$ncount <- NULL
toshiny.hiv.mabs$shm.mean <- NULL
toshiny.den.mabs$ncount <- NULL
toshiny.den.mabs$shm.mean <- NULL
toshiny.cov2.abdab$ncount <- NULL
toshiny.cov2.abdab$shm.mean <- NULL


toshiny.cov2.abdab <- toshiny.cov2.abdab %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

toshiny.hiv.mabs <- toshiny.hiv.mabs %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

toshiny.den.mabs <- toshiny.den.mabs %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

## anything with no cregion?
## any with LC? if so need to make separate hc only file
## which files need renaming
# not
# OAS, den.mabs, HIV.mabs, cov2.abdab

## too many N's in a row? NNNN
## 5 KV's in d13stim
## any with LC? if so need to make separate hc only file
# only these 2, toshiny.den.mabs is HC only

# distthresholdvalue <- 500
# filter(DIST < distthresholdvalue)
# filter(starwars, mass < distthresholdvalue)

## ALSO NEED TO MAKE SURE SEQUENCE ID IS FIRST FOR MABS...AND CDR3 AFTER CREGION
toshiny.cov2.abdab <- toshiny.cov2.abdab %>% relocate(sequence_id)
toshiny.hiv.mabs <- toshiny.hiv.mabs %>% relocate(sequence_id)

toshiny.cov2.abdab <- toshiny.cov2.abdab %>% relocate(cdr3_aa_imgt, .after = cregion)
toshiny.hiv.mabs <- toshiny.hiv.mabs %>% relocate(cdr3_aa_imgt, .after = cregion)

toshiny.cov2.abdab <- toshiny.cov2.abdab %>% relocate(neutralization, .after = cdr3_aa_imgt)
toshiny.cov2.abdab <- toshiny.cov2.abdab %>% relocate(binding, .after = cdr3_aa_imgt)


toshiny.cov2.abdab.h <- subset(toshiny.cov2.abdab, cregion %in% c("IgH"))
toshiny.hiv.mabs.h <- subset(toshiny.hiv.mabs, cregion %in% c("IgH"))
## END ONE-OFF COMMANDS FOR MAB DATASETS


## next save these individually, then combine them!!!
write.table(toshiny.den.bulk.d13stim, "toshiny_den_bulk_d13stim.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.bulk.d13enrich, "toshiny_den_bulk_d13enrich.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.bulk.OAS, "toshiny_den_bulk_OAS.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.bulk.nih45, "toshiny_hiv_bulk_nih45.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.bulk.mt1214, "toshiny_hiv_bulk_mt1214.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hc.BXmay.10mstim, "toshiny_hc_BXmay_10mstim.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.binder.p11, "toshiny_cov2_bulk_binder_p11.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.nielsen.p7450, "toshiny_cov2_bulk_nielsen_p7450.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.galson.p1, "toshiny_cov2_bulk_galson_p1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.kc.m5.allreps, "toshiny_cov2_bulk_kc.m5_allreps.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### mabs
write.table(toshiny.cov2.abdab, "toshiny_cov2_abdab.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.abdab.h, "toshiny_cov2_abdab_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.mabs, "toshiny_den_mabs.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.mabs, "toshiny_hiv_mabs.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.mabs.h, "toshiny_hiv_mabs_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)


## combine d13 enrich & d 13 stim, also one-off potential contamination removed...
toshiny.den.bulk.d13 <- full_join(toshiny.den.bulk.d13enrich, toshiny.den.bulk.d13stim)
toshiny.den.bulk.d13 <- toshiny.den.bulk.d13[ grep("ARALFGLVAVASPFDN", toshiny.den.bulk.d13$cdr3_aa_imgt, invert = TRUE) , ]
toshiny.den.bulk.d13 <- toshiny.den.bulk.d13[ grep("ITPPLYLMVGGVSRAMAV", toshiny.den.bulk.d13$cdr3_aa_imgt, invert = TRUE) , ]
toshiny.den.bulk.d13 <- toshiny.den.bulk.d13[ grep("ARQDRNWFDT", toshiny.den.bulk.d13$cdr3_aa_imgt, invert = TRUE) , ]

write.table(toshiny.den.bulk.d13, "toshiny_den_bulk_d13.tab", sep = "\t", row.names = FALSE, quote = FALSE)


#######################################
#######################################
## NEXT generic CODE TO COMBINE 2-6 DATASETS INTO A SINGLE FILE FOR SHINY VISUALIZATION (BOTH SEPARATE & COMBINED)


## files to use
# toshiny.cov2.abdab  # by cregion
# toshiny.cov2.abdab.h  # by binding & by neutralization

## then 4 combinations
# toshiny.cov2.all
toshiny.cov2.abdab.h
toshiny.cov2.bulk.binder.p11
toshiny.cov2.bulk.galson.p1
toshiny.cov2.bulk.kc.m5.allreps
toshiny.cov2.bulk.nielsen.p7450
toshiny.hc.BXmay.10mstim

# toshiny.cov2hiv
toshiny.cov2.abdab.h
toshiny.hiv.mabs.h

# all hiv
toshiny.hiv.mabs.h
toshiny.hiv.bulk.nih45
toshiny.hiv.bulk.mt1214

# all den
toshiny.den.mabs
toshiny.den.bulk.OAS
toshiny.den.bulk.d13


# toshiny.cov2.all
toshiny.cov2.all <- bind_rows(toshiny.cov2.abdab.h, toshiny.cov2.bulk.binder.p11, toshiny.cov2.bulk.galson.p1, toshiny.cov2.bulk.kc.m5.allreps, toshiny.cov2.bulk.nielsen.p7450, toshiny.hc.BXmay.10mstim, .id = "id")

toshiny.cov2.all$id <- gsub("1","SARS-CoVtwo mAbs",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("2","COVID-nineteen patient Binder peleven",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("3","COVID-nineteen patient Galson pone",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("4","COVID-nineteen patient Kuri-Cervantes mfive",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("5","COVID-nineteen patient Nielsen psevenfourfivezero",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("6","Healthy control bulk repertoire",toshiny.cov2.all$id)

toshiny.cov2.all$id <- gsub("SARS-CoVtwo mAbs","SARS-CoV2 mAbs",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("COVID-nineteen patient Binder peleven","COVID-19 patient Binder p11 bulk repertoire",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("COVID-nineteen patient Galson pone","COVID-19 patient Galson p1 bulk repertoire",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("COVID-nineteen patient Kuri-Cervantes mfive","COVID-19 patient Kuri-Cervantes m5 bulk repertoire",toshiny.cov2.all$id)
toshiny.cov2.all$id <- gsub("COVID-nineteen patient Nielsen psevenfourfivezero","COVID-19 patient Nielsen p7450 bulk repertoire",toshiny.cov2.all$id)

## check id columns!!
unique(toshiny.cov2.all$id)
head(toshiny.cov2.all$sequence_id)
# 
# "SARS-CoVtwo mAbs"
# "COVID-nineteen patient Binder peleven"
# "COVID-nineteen patient Galson pone"
# "COVID-nineteen patient Kuri-Cervantes mfive"
# "COVID-nineteen patient Nielsen psevenfourfivezero"
# "Healthy control"
# 
# "SARS-CoV2 mAbs"
# "COVID-19 patient Binder p11"
# "COVID-19 patient Galson p1"
# "COVID-19 patient Kuri-Cervantes m5"
# "COVID-19 patient Nielsen p7450"

## for all combined need to recalculate ncount and shm_mean
toshiny.cov2.allc <- toshiny.cov2.all
toshiny.cov2.allc$ncount <- NULL
toshiny.cov2.allc$shm_mean <- NULL
toshiny.cov2.allc$shm_max <- NULL
toshiny.cov2.allc <- toshiny.cov2.allc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

## also every combined - add cregion0 and convert cregion to IgH
toshiny.cov2.allc$cregion0 <- toshiny.cov2.allc$cregion
toshiny.cov2.allc$cregion <- "IgH"


# toshiny.cov2.all <- toshiny.cov.all
# toshiny.cov2.allc <- toshiny.cov.allc
# rm(toshiny.cov.all)
# rm(toshiny.cov.allc)

write.table(toshiny.cov2.all, "toshiny_cov2_all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.allc, "toshiny_cov2_allc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

##############
# toshiny.cov2hiv
#    toshiny.cov2.abdab.h
#     toshiny.hiv.mabs.h

#c("anti-CoV2 mAbs", "anti-HIV mAbs"))
# #"SARS-CoV2 mAbs vs. HIV mAbs - IgH",

toshiny.cov2hiv <- bind_rows(toshiny.cov2.abdab.h, toshiny.hiv.mabs.h, .id = "id")

toshiny.cov2hiv$id <- gsub("2","HIV mAbs",toshiny.cov2hiv$id)
toshiny.cov2hiv$id <- gsub("1","SARS-CoV2 mAbs",toshiny.cov2hiv$id)

## check id columns!!
unique(toshiny.cov2hiv$id)
head(toshiny.cov2hiv$sequence_id)
toshiny.cov2hiv$vj_junction_pattern <- NULL


## for all combined need to recalculate ncount and shm_mean
toshiny.cov2hivc <- toshiny.cov2hiv
toshiny.cov2hivc$ncount <- NULL
toshiny.cov2hivc$shm_mean <- NULL
toshiny.cov2hivc$shm_max <- NULL
toshiny.cov2hivc <- toshiny.cov2hivc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))


## also every combined - add cregion0 and convert cregion to IgH
toshiny.cov2hivc$cregion0 <- toshiny.cov2hivc$cregion
toshiny.cov2hivc$cregion <- "IgH"


write.table(toshiny.cov2hiv, "toshiny_cov2hiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2hivc, "toshiny_cov2hivc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

###################
# toshiny.hiv.all
#     toshiny.hiv.mabs.h
#     toshiny.hiv.bulk.nih45
#     toshiny.hiv.bulk.mt1214

#c("HIV mAbs", "HIV patient MT1214", "HIV patient NIH45"))  #"SARS-CoV2 mAbs vs. HIV mAbs - IgH",

toshiny.hiv.all <- bind_rows(toshiny.hiv.mabs.h, toshiny.hiv.bulk.mt1214, toshiny.hiv.bulk.nih45, .id = "id")

toshiny.hiv.all$id <- gsub("1","HIV mAbs",toshiny.hiv.all$id)
toshiny.hiv.all$id <- gsub("3","HIV patient NIH45 bulk repertoire",toshiny.hiv.all$id)
toshiny.hiv.all$id <- gsub("2","HIV patient MT1214 bulk repertoire",toshiny.hiv.all$id)

## check id columns!!
unique(toshiny.hiv.all$id)
head(toshiny.hiv.all$sequence_id)
toshiny.hiv.all$vj_junction_pattern <- NULL


## for all combined need to recalculate ncount and shm_mean
toshiny.hiv.allc <- toshiny.hiv.all
toshiny.hiv.allc$ncount <- NULL
toshiny.hiv.allc$shm_mean <- NULL
toshiny.hiv.allc$shm_max <- NULL
toshiny.hiv.allc <- toshiny.hiv.allc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

## also every combined - add cregion0 and convert cregion to IgH
toshiny.hiv.allc$cregion0 <- toshiny.hiv.allc$cregion
toshiny.hiv.allc$cregion <- "IgH"

write.table(toshiny.hiv.all, "toshiny_hiv_all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.allc, "toshiny_hiv_allc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#################3
# toshiny.den.all
#    toshiny.den.mabs
#    toshiny.den.bulk.d13
#    toshiny.den.bulk.OAS

#c("Dengue plasmablasts", "Dengue patient d13", "Dengue Parameswaran 2013 patients"))  

toshiny.den.all <- bind_rows(toshiny.den.mabs, toshiny.den.bulk.d13, toshiny.den.bulk.OAS, .id = "id")

toshiny.den.all$id <- gsub("1","Dengue plasmablasts",toshiny.den.all$id)
toshiny.den.all$id <- gsub("2","Dengue patient dthirteen",toshiny.den.all$id)
toshiny.den.all$id <- gsub("3","Dengue Parameswaran 2013 patient bulk repertoires",toshiny.den.all$id)

toshiny.den.all$id <- gsub("Dengue patient dthirteen","Dengue patient d13 bulk repertoire",toshiny.den.all$id)

## check id columns!!
unique(toshiny.den.all$id)
head(toshiny.den.all$sequence_id)


## for all combined need to recalculate ncount and shm_mean
toshiny.den.allc <- toshiny.den.all
toshiny.den.allc$ncount <- NULL
toshiny.den.allc$shm_mean <- NULL
toshiny.den.allc$shm_max <- NULL
toshiny.den.allc <- toshiny.den.allc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # ADD MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

## also every combined - add cregion0 and convert cregion to IgH
toshiny.den.allc$cregion0 <- toshiny.den.allc$cregion
toshiny.den.allc$cregion <- "IgH"

write.table(toshiny.den.all, "toshiny_den_all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.allc, "toshiny_den_allc.tab", sep = "\t", row.names = FALSE, quote = FALSE)



# c("SARS-CoV2 mAbs - heavy chains & light chains",
#   "SARS-CoV2 mAbs - IgH by binding",
#   "SARS-CoV2 mAbs - IgH by neutralization",
#   "SARS-CoV2 mAbs vs. 4 COVID-19 patients vs. Healthy control - IgH",
#   "SARS-CoV2 mAbs vs. 4 COVID-19 patients vs. Healthy control - IgH combined",
#   "SARS-CoV2 mAbs vs. HIV mAbs - IgH",
#   "SARS-CoV2 mAbs vs. HIV mAbs - IgH combined",
#   "HIV mAbs vs. HIV patient MT1214 vs. HIV patient NIH45 - IgH",
#   "HIV mAbs vs. HIV patient MT1214 vs. HIV patient NIH45 - IgH combined",
#   "Dengue mAbs vs. Dengue patient d13 vs. Dengue Parameswaran 2013 patients - IgH",
#   "Dengue mAbs vs. Dengue patient d13 vs. Dengue Parameswaran 2013 patients - IgH combined"), selectize = FALSE), 

toshiny.cov2.abdab
toshiny.cov2.abdab.h
toshiny.cov2.abdab.h

toshiny.cov2hiv
toshiny.cov2hivc

toshiny.cov2.all
toshiny.cov2.allc

toshiny.hiv.all
toshiny.hiv.allc

toshiny.den.all
toshiny.den.allc


"SARS-CoV2 mAbs"
"COVID-19 patient Binder p11"
"COVID-19 patient Galson p1"
"COVID-19 patient Kuri-Cervantes m5"
"COVID-19 patient Nielsen p7450"

#"SARS-CoV2 mAbs vs. HIV mAbs - IgH",

# unique(toshiny.cov2.all$id)
# [1] "SARS-CoV2 mAbs"                     "COVID-19 patient Binder p11"        "COVID-19 patient Galson p1"         "COVID-19 patient Kuri-Cervantes m5"
# [5] "COVID-19 patient Nielsen p7450"     "Healthy control"

# > unique(toshiny.cov2hiv$id)
# [1] "SARS-CoV2 mAbs" "HIV mAbs"      
# > 
#   
# > unique(toshiny.hiv.all$id)
# [1] "HIV mAbs"           "HIV patient MT1214" "HIV patient NIH45" 

# unique(toshiny.den.all$id)
# [1] "Dengue plasmablasts"               "Dengue patient d13"                "Dengue Parameswaran 2013 patients"
# > 

## to test plots 
# ggplot(toshiny.cov2.abdab, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.cov2.abdab, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# 
# ggplot(toshiny.hc.BXmay.10mstim, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "% of \nReads  ", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.hc.BXmay.10mstim, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


## STILL GETTING TREE ERRORS IN MAKING 50 OR 300 CLOSEST...
## BECAUSE SEQUENCE_ID TOO LATE?
# MOVE sequence_id	JUST BEFORE CREGION
# MOVE cdr3_aa_imgt JUST AFTER CREGION
# 

## to deploy app on ShinyApps.io  at: ewaltari.shinyapps.io
## MAKE SURE THAT THE TAB INPUT FILES ARE IN THE APP-2 FOLDER
library(rsconnect)
rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-7')

### fyi got this alert when starting shiny app
# Note: Using an external vector in selections is ambiguous.
# ℹ Use `all_of(filteredData1ID)` instead of `filteredData1ID` to silence this message.
# ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
# This message is displayed once per session.


##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING NEW OAS DENGUE DATA

OAS.clusters.p148a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150126_Heavy_Bulk.csv", skip =1)
OAS.clusters.p148b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150229_Heavy_Bulk.csv", skip =1)
OAS.clusters.p148c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150329_Heavy_Bulk.csv", skip =1)
OAS.clusters.p172a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150457_Heavy_Bulk.csv", skip =1)
OAS.clusters.p172b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150481_Heavy_Bulk.csv", skip =1)
OAS.clusters.p172c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150504_Heavy_Bulk.csv", skip =1)
OAS.clusters.p194a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150549_Heavy_Bulk.csv", skip =1)
OAS.clusters.p194b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150573_Heavy_Bulk.csv", skip =1)
OAS.clusters.p194c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150597_Heavy_Bulk.csv", skip =1)
OAS.clusters.p199a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150643_Heavy_Bulk.csv", skip =1)
OAS.clusters.p199b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150668_Heavy_Bulk.csv", skip =1)
OAS.clusters.p199c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150692_Heavy_Bulk.csv", skip =1)
OAS.clusters.p203a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150715_Heavy_Bulk.csv", skip =1)
OAS.clusters.p203b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150734_Heavy_Bulk.csv", skip =1)
OAS.clusters.p203c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150753_Heavy_Bulk.csv", skip =1)
OAS.clusters.p208a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150802_Heavy_Bulk.csv", skip =1)
OAS.clusters.p208b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150816_Heavy_Bulk.csv", skip =1)
OAS.clusters.p208c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150838_Heavy_Bulk.csv", skip =1)
OAS.clusters.p232a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150935_Heavy_Bulk.csv", skip =1)
OAS.clusters.p232b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150947_Heavy_Bulk.csv", skip =1)
OAS.clusters.p232c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2150972_Heavy_Bulk.csv", skip =1)
OAS.clusters.p237a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151066_Heavy_Bulk.csv", skip =1)
OAS.clusters.p237b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151089_Heavy_Bulk.csv", skip =1)
OAS.clusters.p237c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151105_Heavy_Bulk.csv", skip =1)
OAS.clusters.p238a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151162_Heavy_Bulk.csv", skip =1)
OAS.clusters.p238b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151187_Heavy_Bulk.csv", skip =1)
OAS.clusters.p238c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151211_Heavy_Bulk.csv", skip =1)
OAS.clusters.p240a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151231_Heavy_Bulk.csv", skip =1)
OAS.clusters.p240b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151247_Heavy_Bulk.csv", skip =1)
OAS.clusters.p240c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151270_Heavy_Bulk.csv", skip =1)
OAS.clusters.p249a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151293_Heavy_Bulk.csv", skip =1)
OAS.clusters.p249b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151316_Heavy_Bulk.csv", skip =1)
OAS.clusters.p249c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151330_Heavy_Bulk.csv", skip =1)
OAS.clusters.p252a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151353_Heavy_Bulk.csv", skip =1)
OAS.clusters.p252b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151376_Heavy_Bulk.csv", skip =1)
OAS.clusters.p252c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151395_Heavy_Bulk.csv", skip =1)
OAS.clusters.p255a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151414_Heavy_Bulk.csv", skip =1)
OAS.clusters.p255b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151435_Heavy_Bulk.csv", skip =1)
OAS.clusters.p255c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151450_Heavy_Bulk.csv", skip =1)
OAS.clusters.p265a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151499_Heavy_Bulk.csv", skip =1)
OAS.clusters.p265b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151523_Heavy_Bulk.csv", skip =1)

OAS.clusters.p275a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151538_Heavy_Bulk.csv", skip =1)
OAS.clusters.p275b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151562_Heavy_Bulk.csv", skip =1)
OAS.clusters.p275c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151598_Heavy_Bulk.csv", skip =1)
OAS.clusters.p276a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151713_Heavy_Bulk.csv", skip =1)
OAS.clusters.p276b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151738_Heavy_Bulk.csv", skip =1)
OAS.clusters.p276c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2151761_Heavy_Bulk.csv", skip =1)
OAS.clusters.p287a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153023_Heavy_Bulk.csv", skip =1)
OAS.clusters.p287b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153024_Heavy_Bulk.csv", skip =1)
OAS.clusters.p287c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153025_Heavy_Bulk.csv", skip =1)
OAS.clusters.p289a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153026_Heavy_Bulk.csv", skip =1)
OAS.clusters.p289b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153027_Heavy_Bulk.csv", skip =1)
OAS.clusters.p289c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153028_Heavy_Bulk.csv", skip =1)
OAS.clusters.p299a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153029_Heavy_Bulk.csv", skip =1)
OAS.clusters.p299b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153030_Heavy_Bulk.csv", skip =1)
OAS.clusters.p299c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153031_Heavy_Bulk.csv", skip =1)
OAS.clusters.p301a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153032_Heavy_Bulk.csv", skip =1)
OAS.clusters.p301b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153033_Heavy_Bulk.csv", skip =1)
OAS.clusters.p301c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153034_Heavy_Bulk.csv", skip =1)
OAS.clusters.p307a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153035_Heavy_Bulk.csv", skip =1)
OAS.clusters.p307b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153036_Heavy_Bulk.csv", skip =1)

OAS.clusters.p311a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153037_Heavy_Bulk.csv", skip =1)
OAS.clusters.p311b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153038_Heavy_Bulk.csv", skip =1)

OAS.clusters.p320a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153039_Heavy_Bulk.csv", skip =1)
OAS.clusters.p320b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153040_Heavy_Bulk.csv", skip =1)

OAS.clusters.p346a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153045_Heavy_Bulk.csv", skip =1)

OAS.clusters.p376a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153046_Heavy_Bulk.csv", skip =1)
OAS.clusters.p376b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153047_Heavy_Bulk.csv", skip =1)
OAS.clusters.p376c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153048_Heavy_Bulk.csv", skip =1)
OAS.clusters.p391a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153049_Heavy_Bulk.csv", skip =1)
OAS.clusters.p391b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153050_Heavy_Bulk.csv", skip =1)

OAS.clusters.p422a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153052_Heavy_Bulk.csv", skip =1)
OAS.clusters.p422b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153053_Heavy_Bulk.csv", skip =1)
OAS.clusters.p422c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153054_Heavy_Bulk.csv", skip =1)
OAS.clusters.p444a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153056_Heavy_Bulk.csv", skip =1)
OAS.clusters.p444b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153057_Heavy_Bulk.csv", skip =1)
OAS.clusters.p444c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153058_Heavy_Bulk.csv", skip =1)
OAS.clusters.p455a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153060_Heavy_Bulk.csv", skip =1)
OAS.clusters.p455b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153061_Heavy_Bulk.csv", skip =1)
OAS.clusters.p455c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153062_Heavy_Bulk.csv", skip =1)
OAS.clusters.p479a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153063_Heavy_Bulk.csv", skip =1)
OAS.clusters.p479b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153064_Heavy_Bulk.csv", skip =1)
OAS.clusters.p479c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153065_Heavy_Bulk.csv", skip =1)
OAS.clusters.p481a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153066_Heavy_Bulk.csv", skip =1)
OAS.clusters.p481b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153067_Heavy_Bulk.csv", skip =1)
OAS.clusters.p481c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153068_Heavy_Bulk.csv", skip =1)
OAS.clusters.p489a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153070_Heavy_Bulk.csv", skip =1)
OAS.clusters.p489b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153071_Heavy_Bulk.csv", skip =1)
OAS.clusters.p489c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153072_Heavy_Bulk.csv", skip =1)
OAS.clusters.p500a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153073_Heavy_Bulk.csv", skip =1)
OAS.clusters.p500b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153074_Heavy_Bulk.csv", skip =1)
OAS.clusters.p500c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153075_Heavy_Bulk.csv", skip =1)
OAS.clusters.p514a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153230_Heavy_Bulk.csv", skip =1)
OAS.clusters.p514b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153231_Heavy_Bulk.csv", skip =1)
OAS.clusters.p514c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153232_Heavy_Bulk.csv", skip =1)
OAS.clusters.p515a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153233_Heavy_Bulk.csv", skip =1)
OAS.clusters.p515b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153234_Heavy_Bulk.csv", skip =1)
OAS.clusters.p515c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153235_Heavy_Bulk.csv", skip =1)
OAS.clusters.p517a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153236_Heavy_Bulk.csv", skip =1)
OAS.clusters.p517b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153237_Heavy_Bulk.csv", skip =1)
OAS.clusters.p517c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153238_Heavy_Bulk.csv", skip =1)
OAS.clusters.p520a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153239_Heavy_Bulk.csv", skip =1)
OAS.clusters.p520b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153240_Heavy_Bulk.csv", skip =1)
OAS.clusters.p520c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153241_Heavy_Bulk.csv", skip =1)
OAS.clusters.p524a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153242_Heavy_Bulk.csv", skip =1)
OAS.clusters.p524b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153243_Heavy_Bulk.csv", skip =1)
OAS.clusters.p524c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153244_Heavy_Bulk.csv", skip =1)
OAS.clusters.p529a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153245_Heavy_Bulk.csv", skip =1)
OAS.clusters.p529b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153247_Heavy_Bulk.csv", skip =1)
OAS.clusters.p529c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153248_Heavy_Bulk.csv", skip =1)
OAS.clusters.p543a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153249_Heavy_Bulk.csv", skip =1)
OAS.clusters.p543b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153250_Heavy_Bulk.csv", skip =1)
OAS.clusters.p543c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153251_Heavy_Bulk.csv", skip =1)
OAS.clusters.p551a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153252_Heavy_Bulk.csv", skip =1)
OAS.clusters.p551b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153253_Heavy_Bulk.csv", skip =1)
OAS.clusters.p551c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153254_Heavy_Bulk.csv", skip =1)
OAS.clusters.p555a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153255_Heavy_Bulk.csv", skip =1)
OAS.clusters.p555b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153256_Heavy_Bulk.csv", skip =1)
OAS.clusters.p555c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153258_Heavy_Bulk.csv", skip =1)
OAS.clusters.p558a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153261_Heavy_Bulk.csv", skip =1)
OAS.clusters.p558b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153262_Heavy_Bulk.csv", skip =1)

OAS.clusters.p563a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153263_Heavy_Bulk.csv", skip =1)
OAS.clusters.p563b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153264_Heavy_Bulk.csv", skip =1)
OAS.clusters.p563c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153265_Heavy_Bulk.csv", skip =1)
OAS.clusters.p569a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153266_Heavy_Bulk.csv", skip =1)
OAS.clusters.p569b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153267_Heavy_Bulk.csv", skip =1)
OAS.clusters.p569c <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153268_Heavy_Bulk.csv", skip =1)

OAS.clusters.p346 <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/oas/SRR2153045_Heavy_Bulk.csv", skip =1)

#####


OAS.clusters.p148 <- rbind(OAS.clusters.p148a,OAS.clusters.p148b,OAS.clusters.p148c)
OAS.clusters.p172 <- rbind(OAS.clusters.p172a,OAS.clusters.p172b,OAS.clusters.p172c)
OAS.clusters.p194 <- rbind(OAS.clusters.p194a,OAS.clusters.p194b,OAS.clusters.p194c)
OAS.clusters.p199 <- rbind(OAS.clusters.p199a,OAS.clusters.p199b,OAS.clusters.p199c)
OAS.clusters.p203 <- rbind(OAS.clusters.p203a,OAS.clusters.p203b,OAS.clusters.p203c)
OAS.clusters.p208 <- rbind(OAS.clusters.p208a,OAS.clusters.p208b,OAS.clusters.p208c)
OAS.clusters.p232 <- rbind(OAS.clusters.p232a,OAS.clusters.p232b,OAS.clusters.p232c)
OAS.clusters.p237 <- rbind(OAS.clusters.p237a,OAS.clusters.p237b,OAS.clusters.p237c)
OAS.clusters.p238 <- rbind(OAS.clusters.p238a,OAS.clusters.p238b,OAS.clusters.p238c)
OAS.clusters.p240 <- rbind(OAS.clusters.p240a,OAS.clusters.p240b,OAS.clusters.p240c)
OAS.clusters.p249 <- rbind(OAS.clusters.p249a,OAS.clusters.p249b,OAS.clusters.p249c)
OAS.clusters.p252 <- rbind(OAS.clusters.p252a,OAS.clusters.p252b,OAS.clusters.p252c)
OAS.clusters.p255 <- rbind(OAS.clusters.p255a,OAS.clusters.p255b,OAS.clusters.p255c)
OAS.clusters.p265 <- rbind(OAS.clusters.p265a,OAS.clusters.p265b)
OAS.clusters.p275 <- rbind(OAS.clusters.p275a,OAS.clusters.p275b,OAS.clusters.p275c)
OAS.clusters.p276 <- rbind(OAS.clusters.p276a,OAS.clusters.p276b,OAS.clusters.p276c)
OAS.clusters.p287 <- rbind(OAS.clusters.p287a,OAS.clusters.p287b,OAS.clusters.p287c)
OAS.clusters.p289 <- rbind(OAS.clusters.p289a,OAS.clusters.p289b,OAS.clusters.p289c)
OAS.clusters.p299 <- rbind(OAS.clusters.p299a,OAS.clusters.p299b,OAS.clusters.p299c)
OAS.clusters.p301 <- rbind(OAS.clusters.p301a,OAS.clusters.p301b,OAS.clusters.p301c)
OAS.clusters.p307 <- rbind(OAS.clusters.p307a,OAS.clusters.p307b)
OAS.clusters.p311 <- rbind(OAS.clusters.p311a,OAS.clusters.p311b)
OAS.clusters.p320 <- rbind(OAS.clusters.p320a,OAS.clusters.p320b)

OAS.clusters.p376 <- rbind(OAS.clusters.p376a,OAS.clusters.p376b,OAS.clusters.p376c)
OAS.clusters.p391 <- rbind(OAS.clusters.p391a,OAS.clusters.p391b)
OAS.clusters.p422 <- rbind(OAS.clusters.p422a,OAS.clusters.p422b,OAS.clusters.p422c)
OAS.clusters.p444 <- rbind(OAS.clusters.p444a,OAS.clusters.p444b,OAS.clusters.p444c)
OAS.clusters.p455 <- rbind(OAS.clusters.p455a,OAS.clusters.p455b,OAS.clusters.p455c)
OAS.clusters.p479 <- rbind(OAS.clusters.p479a,OAS.clusters.p479b,OAS.clusters.p479c)
OAS.clusters.p481 <- rbind(OAS.clusters.p481a,OAS.clusters.p481b,OAS.clusters.p481c)
OAS.clusters.p489 <- rbind(OAS.clusters.p489a,OAS.clusters.p489b,OAS.clusters.p489c)
OAS.clusters.p500 <- rbind(OAS.clusters.p500a,OAS.clusters.p500b,OAS.clusters.p500c)
OAS.clusters.p514 <- rbind(OAS.clusters.p514a,OAS.clusters.p514b,OAS.clusters.p514c)
OAS.clusters.p515 <- rbind(OAS.clusters.p515a,OAS.clusters.p515b,OAS.clusters.p515c)
OAS.clusters.p517 <- rbind(OAS.clusters.p517a,OAS.clusters.p517b,OAS.clusters.p517c)
OAS.clusters.p520 <- rbind(OAS.clusters.p520a,OAS.clusters.p520b,OAS.clusters.p520c)
OAS.clusters.p524 <- rbind(OAS.clusters.p524a,OAS.clusters.p524b,OAS.clusters.p524c)
OAS.clusters.p529 <- rbind(OAS.clusters.p529a,OAS.clusters.p529b,OAS.clusters.p529c)
OAS.clusters.p543 <- rbind(OAS.clusters.p543a,OAS.clusters.p543b,OAS.clusters.p543c)
OAS.clusters.p551 <- rbind(OAS.clusters.p551a,OAS.clusters.p551b,OAS.clusters.p551c)
OAS.clusters.p555 <- rbind(OAS.clusters.p555a,OAS.clusters.p555b,OAS.clusters.p555c)
OAS.clusters.p558 <- rbind(OAS.clusters.p558a,OAS.clusters.p558b)
OAS.clusters.p563 <- rbind(OAS.clusters.p563a,OAS.clusters.p563b,OAS.clusters.p563c)
OAS.clusters.p569 <- rbind(OAS.clusters.p569a,OAS.clusters.p569b,OAS.clusters.p569c)

#####
OAS.clusters.p148$dataset <- 'Parameswaran_2013_p148'
OAS.clusters.p172$dataset <- 'Parameswaran_2013_p172'
OAS.clusters.p194$dataset <- 'Parameswaran_2013_p194'
OAS.clusters.p199$dataset <- 'Parameswaran_2013_p199'
OAS.clusters.p203$dataset <- 'Parameswaran_2013_p203'
OAS.clusters.p208$dataset <- 'Parameswaran_2013_p208'
OAS.clusters.p232$dataset <- 'Parameswaran_2013_p232'
OAS.clusters.p237$dataset <- 'Parameswaran_2013_p237'
OAS.clusters.p238$dataset <- 'Parameswaran_2013_p238'
OAS.clusters.p240$dataset <- 'Parameswaran_2013_p240'
OAS.clusters.p249$dataset <- 'Parameswaran_2013_p249'
OAS.clusters.p252$dataset <- 'Parameswaran_2013_p252'
OAS.clusters.p255$dataset <- 'Parameswaran_2013_p255'
OAS.clusters.p265$dataset <- 'Parameswaran_2013_p265'
OAS.clusters.p275$dataset <- 'Parameswaran_2013_p275'
OAS.clusters.p276$dataset <- 'Parameswaran_2013_p276'
OAS.clusters.p287$dataset <- 'Parameswaran_2013_p287'
OAS.clusters.p289$dataset <- 'Parameswaran_2013_p289'
OAS.clusters.p299$dataset <- 'Parameswaran_2013_p299'
OAS.clusters.p301$dataset <- 'Parameswaran_2013_p301'
OAS.clusters.p307$dataset <- 'Parameswaran_2013_p307'
OAS.clusters.p311$dataset <- 'Parameswaran_2013_p311'
OAS.clusters.p320$dataset <- 'Parameswaran_2013_p320'
OAS.clusters.p346$dataset <- 'Parameswaran_2013_p346'
OAS.clusters.p376$dataset <- 'Parameswaran_2013_p376'
OAS.clusters.p391$dataset <- 'Parameswaran_2013_p391'
OAS.clusters.p422$dataset <- 'Parameswaran_2013_p422'
OAS.clusters.p444$dataset <- 'Parameswaran_2013_p444'
OAS.clusters.p455$dataset <- 'Parameswaran_2013_p455'
OAS.clusters.p479$dataset <- 'Parameswaran_2013_p479'
OAS.clusters.p481$dataset <- 'Parameswaran_2013_p481'
OAS.clusters.p489$dataset <- 'Parameswaran_2013_p489'
OAS.clusters.p500$dataset <- 'Parameswaran_2013_p500'
OAS.clusters.p514$dataset <- 'Parameswaran_2013_p514'
OAS.clusters.p515$dataset <- 'Parameswaran_2013_p515'
OAS.clusters.p517$dataset <- 'Parameswaran_2013_p517'
OAS.clusters.p520$dataset <- 'Parameswaran_2013_p520'
OAS.clusters.p524$dataset <- 'Parameswaran_2013_p524'
OAS.clusters.p529$dataset <- 'Parameswaran_2013_p529'
OAS.clusters.p543$dataset <- 'Parameswaran_2013_p543'
OAS.clusters.p551$dataset <- 'Parameswaran_2013_p551'
OAS.clusters.p555$dataset <- 'Parameswaran_2013_p555'
OAS.clusters.p558$dataset <- 'Parameswaran_2013_p558'
OAS.clusters.p563$dataset <- 'Parameswaran_2013_p563'
OAS.clusters.p569$dataset <- 'Parameswaran_2013_p569'

####

OAS.clusters.p148$obs <- 1:nrow(OAS.clusters.p148)
OAS.clusters.p172$obs <- 1:nrow(OAS.clusters.p172)
OAS.clusters.p194$obs <- 1:nrow(OAS.clusters.p194)
OAS.clusters.p199$obs <- 1:nrow(OAS.clusters.p199)
OAS.clusters.p203$obs <- 1:nrow(OAS.clusters.p203)
OAS.clusters.p208$obs <- 1:nrow(OAS.clusters.p208)
OAS.clusters.p232$obs <- 1:nrow(OAS.clusters.p232)
OAS.clusters.p237$obs <- 1:nrow(OAS.clusters.p237)
OAS.clusters.p238$obs <- 1:nrow(OAS.clusters.p238)
OAS.clusters.p240$obs <- 1:nrow(OAS.clusters.p240)
OAS.clusters.p249$obs <- 1:nrow(OAS.clusters.p249)
OAS.clusters.p252$obs <- 1:nrow(OAS.clusters.p252)
OAS.clusters.p255$obs <- 1:nrow(OAS.clusters.p255)
OAS.clusters.p265$obs <- 1:nrow(OAS.clusters.p265)
OAS.clusters.p275$obs <- 1:nrow(OAS.clusters.p275)
OAS.clusters.p276$obs <- 1:nrow(OAS.clusters.p276)
OAS.clusters.p287$obs <- 1:nrow(OAS.clusters.p287)
OAS.clusters.p289$obs <- 1:nrow(OAS.clusters.p289)
OAS.clusters.p299$obs <- 1:nrow(OAS.clusters.p299)
OAS.clusters.p301$obs <- 1:nrow(OAS.clusters.p301)
OAS.clusters.p307$obs <- 1:nrow(OAS.clusters.p307)
OAS.clusters.p311$obs <- 1:nrow(OAS.clusters.p311)
OAS.clusters.p320$obs <- 1:nrow(OAS.clusters.p320)
OAS.clusters.p346$obs <- 1:nrow(OAS.clusters.p346)
OAS.clusters.p376$obs <- 1:nrow(OAS.clusters.p376)
OAS.clusters.p391$obs <- 1:nrow(OAS.clusters.p391)
OAS.clusters.p422$obs <- 1:nrow(OAS.clusters.p422)
OAS.clusters.p444$obs <- 1:nrow(OAS.clusters.p444)
OAS.clusters.p455$obs <- 1:nrow(OAS.clusters.p455)
OAS.clusters.p479$obs <- 1:nrow(OAS.clusters.p479)
OAS.clusters.p481$obs <- 1:nrow(OAS.clusters.p481)
OAS.clusters.p489$obs <- 1:nrow(OAS.clusters.p489)
OAS.clusters.p500$obs <- 1:nrow(OAS.clusters.p500)
OAS.clusters.p514$obs <- 1:nrow(OAS.clusters.p514)
OAS.clusters.p515$obs <- 1:nrow(OAS.clusters.p515)
OAS.clusters.p517$obs <- 1:nrow(OAS.clusters.p517)
OAS.clusters.p520$obs <- 1:nrow(OAS.clusters.p520)
OAS.clusters.p524$obs <- 1:nrow(OAS.clusters.p524)
OAS.clusters.p529$obs <- 1:nrow(OAS.clusters.p529)
OAS.clusters.p543$obs <- 1:nrow(OAS.clusters.p543)
OAS.clusters.p551$obs <- 1:nrow(OAS.clusters.p551)
OAS.clusters.p555$obs <- 1:nrow(OAS.clusters.p555)
OAS.clusters.p558$obs <- 1:nrow(OAS.clusters.p558)
OAS.clusters.p563$obs <- 1:nrow(OAS.clusters.p563)
OAS.clusters.p569$obs <- 1:nrow(OAS.clusters.p569)

####


OAS.clusters.p148 <- OAS.clusters.p148 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p172 <- OAS.clusters.p172 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p194 <- OAS.clusters.p194 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p199 <- OAS.clusters.p199 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p203 <- OAS.clusters.p203 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p208 <- OAS.clusters.p208 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p232 <- OAS.clusters.p232 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p237 <- OAS.clusters.p237 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p238 <- OAS.clusters.p238 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p240 <- OAS.clusters.p240 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p249 <- OAS.clusters.p249 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p252 <- OAS.clusters.p252 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p255 <- OAS.clusters.p255 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p265 <- OAS.clusters.p265 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p275 <- OAS.clusters.p275 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p276 <- OAS.clusters.p276 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p287 <- OAS.clusters.p287 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p289 <- OAS.clusters.p289 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p299 <- OAS.clusters.p299 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p301 <- OAS.clusters.p301 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p307 <- OAS.clusters.p307 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p311 <- OAS.clusters.p311 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p320 <- OAS.clusters.p320 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p346 <- OAS.clusters.p346 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p376 <- OAS.clusters.p376 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p391 <- OAS.clusters.p391 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p422 <- OAS.clusters.p422 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p444 <- OAS.clusters.p444 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p455 <- OAS.clusters.p455 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p479 <- OAS.clusters.p479 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p481 <- OAS.clusters.p481 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p489 <- OAS.clusters.p489 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p500 <- OAS.clusters.p500 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p514 <- OAS.clusters.p514 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p515 <- OAS.clusters.p515 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p517 <- OAS.clusters.p517 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p520 <- OAS.clusters.p520 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p524 <- OAS.clusters.p524 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p529 <- OAS.clusters.p529 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p543 <- OAS.clusters.p543 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p551 <- OAS.clusters.p551 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p555 <- OAS.clusters.p555 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p558 <- OAS.clusters.p558 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p563 <- OAS.clusters.p563 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
OAS.clusters.p569 <- OAS.clusters.p569 %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)


OAS.clusters.all <- rbind(OAS.clusters.p148,OAS.clusters.p172,OAS.clusters.p194,OAS.clusters.p199,OAS.clusters.p203,OAS.clusters.p208,OAS.clusters.p232,OAS.clusters.p237,OAS.clusters.p238,OAS.clusters.p240,OAS.clusters.p249,OAS.clusters.p252,OAS.clusters.p255,OAS.clusters.p265,OAS.clusters.p275,OAS.clusters.p276,OAS.clusters.p287,OAS.clusters.p289,OAS.clusters.p299,OAS.clusters.p301,OAS.clusters.p307,OAS.clusters.p311,OAS.clusters.p320,OAS.clusters.p346,OAS.clusters.p376,OAS.clusters.p391,OAS.clusters.p422,OAS.clusters.p444,OAS.clusters.p455,OAS.clusters.p479,OAS.clusters.p481,OAS.clusters.p489,OAS.clusters.p500,OAS.clusters.p514,OAS.clusters.p515,OAS.clusters.p517,OAS.clusters.p520,OAS.clusters.p524,OAS.clusters.p529,OAS.clusters.p543,OAS.clusters.p551,OAS.clusters.p555,OAS.clusters.p558,OAS.clusters.p563,OAS.clusters.p569)

OAS.clusters.all$ANARCI_numbering <- NULL
OAS.clusters.all$ANARCI_status <- NULL
write.table(OAS.clusters.all, "OAS_sept21_germ-pass2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

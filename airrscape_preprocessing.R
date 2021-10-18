
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

# pkgs = c("alakazam", "igraph", "dplyr","RColorBrewer", "hexbin", "scales","grid", "lattice", "gdata","gridExtra", "ape", "shazam","limma", "reshape2", "DT","ggplot2", "seqinr", "phangorn","shiny", "tidyverse") # package names

pkgs = c("alakazam", "igraph", "dplyr","RColorBrewer", "hexbin", "scales","grid", "lattice", "gdata","gridExtra", "ape", "shazam","reshape2", "DT","ggplot2", "seqinr", "phangorn","shiny", "tidyverse") # package names
# install.packages(pkgs)
#   inst = lapply(pkgs, library, character.only = TRUE) # load them

  ## Warning in install.packages :
  # package ‘grid’ is a base package, and should not be updated
  # Warning in install.packages :
  #   package ‘limma’ is not available for this version of R
  
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}


#dir.create("~/code")
## this works ONLY IF YOU HAVE CREATED THE DIRECTORY FIRST...
setwd("/Users/eric.waltari/data_carpentry/AIRRScape")


##################################################################################################################
### COMMANDS FOR READING IN AIRR-COMPLIANT DATASETS FOR SHINY APP
### note that for non-AIRR (e.g. older Immcantation data or databases like AbDab) datasets, column names will be different
######### solution is to use only AIRR datasets, fixed as of Sept 2021
###### particularly note that 1) junction MAY NOT NEED TO BE TRIMMED 1 AA ON EACH END TO BE KABAT CDR3
######                        2) v_identity MAY BE BETWEEN 0-1 NOT 0-100, AND THUS shm FORMULA NEEDS TO CHANGE
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


sle.tipton.p1 <- read_tsv("Tipton_sle1_germ-pass.tsv")

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
  # x <- x[ grep("NNNN", x$sequence, invert = TRUE) , ]   # REMOVE FOR KC DATASET...NO LONGER USING
  ### removing all sequences with IMGT CDR3 less than 3
  x <- x %>% filter(cdr3length_imgt > 2.8)  
  ## next lines create Vgene V gene family, J gene columns
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
  ## this will remove all redundant sequences with same gf/gene & cdr3 motif...note we count above so fine to collapse here!!
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

## NOW IT WORKS WITH THE 1 ARGUMENT!!
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

toshiny.sle.tipton.p1 <- shinyprocess(sle.tipton.p1)



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


### TO DO OCT 2021
## NEED TO MAKE CREGION0 AND ALSO CHANGE CREGION TO IGH FOR ALL COMBINED .H DATASETS...
## LOOK INTO WHY HIV DATASETS ARE MISSING CDR3_AA (VS CDR3_AA_IMGT
## ISSUE IS TOSHINY.HIV.MABS - NEED TO CHANGE CDR3_AA TO CDR3_AA_IMGT

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


ggplot(toshiny.cov2.abdab, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.cov2.abdab, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

ggplot(toshiny.sle.tipton.p1, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "% of \nReads  ", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
ggplot(toshiny.hc.BXmay.10mstim, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill= (..count..)*100/tapply(..count..,..PANEL..,sum)[..PANEL..])) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "% of \nReads  ", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


ggplot(toshiny.hc.BXmay.10mstim, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))


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

## might not need ANYTHING BELOW...
############################################################################################################################################################
############################################################################################################################################################

### IDEA - TO GET MEAN SHM & NCOUNT, BUT THEN SELECT ONLY ONE PER GROUP...WILL KEEP THE ORIGINAL NCOUNT (of cdr3 length, but think ok for the plot)
# THEN GROUP BY V CALL JCALL AND junction_aa, SELECT ONLY ONE PER GROUP
toshiny.den.bulk.d13stim.few <- toshiny.den.bulk.d13stim %>%
  group_by(cdr3_aa_imgt,gf_jgene) %>%
  summarize_all(first)


## next HC only? any need to filter more for larger datasets?
## THEN GROUP BY V CALL JCALL AND junction_aa, SELECT ONLY ONE PER GROUP
# 
# BX.clusters.rd1214fewest <- BX.clusters.rd1214fewer %>%
#   group_by(junction_aa,gf_jgene) %>%
#   summarize_all(first)

## more filters
#toshiny.fewhiv <- toshiny.fewhiv[ grep("Kappa|Lambda|IGKJ|__", toshiny.fewhiv$jgene, invert = TRUE) , ]

BX.clusters.rd1214$cregion <- str_sub(BX.clusters.rd1214$v_call, end=3)

## ALSO NO *
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer[ grep("\\*", BX.clusters.rd1214fewer$cdr3_aa_imgt, invert = TRUE) , ]

## some IGK & IGL
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer[ grep("IGK|IGL", BX.clusters.rd1214fewer$gf_jgene, invert = TRUE) , ]
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer[ grep("kappa|lambda", BX.clusters.rd1214fewer$PRIMER, invert = TRUE) , ]

BX.clusters.rd1214 <- BX.clusters.rd1214 %>% filter(FUNCTIONAL == "T") %>%
  filter(IN_FRAME == "T") %>%
  filter(STOP == "F")
#filter(source != "Augmenta")
cov2.bulk.kc.m5.rep1b <- cov2.bulk.kc.m5.rep1 %>% filter(productive != "FALSE") %>%
  filter(vj_in_frame != "FALSE") %>%
  filter(productive != "F") %>%
  filter(vj_in_frame != "F")

cov2.bulk.kc.m5.rep1b <- cov2.bulk.kc.m5.rep1b[ grep("\\*", cov2.bulk.kc.m5.rep1b$cdr3_aa_imgt, invert = TRUE) , ]

# productive = col_logical(),
# stop_codon = col_logical(),
# vj_in_frame = col_logical(),

## DO ANY OF THESE DATASETS INCLUDE LIGHT CHAINS? IF SO NEED TO EXTRACT ONLY IGH...

## renumbering
x$dataset <- paste0(x)
x$obs <- 1:nrow(x) 
x <- x %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
file <- paste0(x, ".csv")

BX.clusters.test

test <- paste0(BX.clusters.test)

BX.clusters.test$dataset <- deparse(quote(BX.clusters.test))

BX.clusters.test$dataset <- deparse(substitute(BX.clusters.test))
BX.clusters.test$dataset <- gsub("\\.","\\_",BX.clusters.test$dataset)

BX.clusters.test$obs <- 1:nrow(BX.clusters.test) 
BX.clusters.test <- BX.clusters.test %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
# file <- paste0(x, ".csv")
BX.clusters.test <- BX.clusters.test %>%
  add_count(gf_jgene,cdr3length_imgt) %>% 
  rename(ncount = n) %>% 
  group_by(gf_jgene,cdr3length_imgt) %>% 
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

BX.clusters.test <- BX.clusters.test %>% add_count(clone_id) %>%
  rename(reads_per_clone = n)

clone_id_var <- c("clone_id")
BX.clusters.test <- BX.clusters.test %>% add_count(any_of(clone_id_var))
  

## then need to combine - first all cov2, hiv, den datasets


# str_sub(BX.clusters.test$junction_aa, 1, 1) <- ""
# str_sub(BX.clusters.test$junction_aa, -1, -1) <- ""
# BX.clusters.mabs$cdr3length_imgt <- ((BX.clusters.mabs$junction_length) / 3) - 2

## convert sequence names if not already...

BX.clusters.rd1214e$junction_aatrim <- BX.clusters.rd1214e$junction_aa
str_sub(BX.clusters.rd1214e$junction_aatrim, -1, -1) <- ""
str_sub(BX.clusters.rd1214e$junction_aatrim, 1, 1) <- ""

if (mean(BX.clusters.momseqexample$v_identity) < 1) {
  BX.clusters.momseqexample$shm <- (100 - (BX.clusters.momseqexample$v_identity * 100))
} else {
  BX.clusters.momseqexample$shm <- (100 - BX.clusters.momseqexample$v_identity)
}


######################################################################
BX.clusters.test <- read_tsv("/Users/eric.waltari/data_carpentry/wikipathways/bcells_1M_germ-pass.tsv")

BX.clusters.test$junction_aa <- as.character(BX.clusters.test$junction_aa)
BX.clusters.test$cdr3length_imgt <- nchar(BX.clusters.test$junction_aa)

BX.clusters.test$gene <- getGene(BX.clusters.test$germline_v_call, first=TRUE, strip_d=TRUE)
BX.clusters.test$gf <- substring(BX.clusters.test$gene, 1,5)
BX.clusters.test$jgene <- getGene(BX.clusters.test$germline_j_call, first=TRUE, strip_d=TRUE)
BX.clusters.test$jgf <- substring(BX.clusters.test$jgene, 1,5)

BX.clusters.test <- BX.clusters.test %>% add_count(clone_id) %>%
  rename(reads_per_clone = n)

if (mean(BX.clusters.test$v_identity) < 1) {
  BX.clusters.test$shm <- (100 - (BX.clusters.test$v_identity * 100))
} else {
  BX.clusters.test$shm <- (100 - BX.clusters.test$v_identity)
}

BX.clusters.test <- BX.clusters.test %>% unite(gf_jgene, gf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)

## MAY NOT ALWAYS WANT TO DO THIS - also note IgD in here (could also remove??)
BX.clusters.test$cregion <- gsub("IgA","IgH",BX.clusters.test$cregion)
BX.clusters.test$cregion <- gsub("IgG","IgH",BX.clusters.test$cregion)
BX.clusters.test$cregion <- gsub("IgM","IgH",BX.clusters.test$cregion)
BX.clusters.test$cregion <- gsub("IgD","IgH",BX.clusters.test$cregion)

BX.clusters.test2 <- BX.clusters.test %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2)) %>%
  rename(sequence_id = sequence_id) %>%
  rename(cregion = cregion)


BX.clusters.test2 <- BX.clusters.test %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2)) 


## all of these mutate_at(vars(shm_mean), #funs(round(., 2)))
## are now mutate(across(shm_mean, round, 2)) 


BX.clusters.test <- BX.clusters.test %>% filter(cdr3length_imgt > 3.8)  
##

## mutate_at has been superseded by across()
## also funs() is deprecated, use a list
#df %>% mutate(across(cols, round, 3)) 

toshiny.test <- BX.clusters.test %>% select(sequence_id,cregion,junction_aa,gene,gf_jgene,gf,jgene,jgf,cdr3length_imgt,shm,ncount,shm_mean,reads_per_clone) %>% 
  filter(!is.na(cdr3length_imgt)) %>% 
  filter(is.wholenumber(cdr3length_imgt)) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_mean, round, 2))

rm(BX.clusters.test)

toshiny.test.h <- subset(toshiny.test, cregion %in% c("IgH"))

#write.table(toshiny.test, "toshiny_test.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.test.h, "toshiny_test_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# ggplot(toshiny.test.h, aes(gf_jgene,cdr3length_imgt)) + geom_tile(aes(fill = shm_mean)) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "fixed") + scale_fill_viridis_c(name = "Mean \nSomatic \nHypermutation (%)", option = "C") + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
# ggplot(toshiny.test.h, aes(gf_jgene,cdr3length_imgt)) + geom_bin2d(aes(fill=log10(..count..))) + scale_y_continuous(limits = c(3, 42)) + theme_bw() + ylab("CDRH/L3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free_x") + scale_fill_viridis_c(name = "# of \nReads", option = "C", breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))

## columns used in app:
# sequence_id,junction_aa,gene,gf_jgene,cdr3length_imgt, cregion, shm_mean, shm


## NEXT geneRIC CODE TO COMBINE 2-6 DATASETS INTO A SINGLE FILE FOR SHINY VISUALIZATION (BOTH SEPARATE & COMBINED)

toshiny.test.many <- bind_rows(toshiny.test.h, toshiny.test2.h, toshiny.test3.h, toshiny.test4.h, .id = "id")

toshiny.test.many$id <- gsub("1","test1",toshiny.fewhcbg$id)
toshiny.test.many$id <- gsub("2","test2",toshiny.fewhcbg$id)
toshiny.test.many$id <- gsub("3","test3",toshiny.fewhcbg$id)
toshiny.test.many$id <- gsub("4","test4",toshiny.fewhcbg$id)

## for combined need to recalculate ncount and shm_mean
toshiny.test.manyc <- toshiny.test.many
toshiny.test.manyc$ncount <- NULL
toshiny.test.manyc$shm_mean <- NULL

toshiny.test.manyc <- toshiny.fewhcbgc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))

#write.table(toshiny.test.many, "toshiny_testmany.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.test.manyc, "toshiny_testmanyc.tab", sep = "\t", row.names = FALSE, quote = FALSE)


##################################################################################################################
### END COMMANDS FOR READING IN AIRR-COMPLIANT DATASETS FOR SHINY APP


#################################
### COMMANDS FOR UPDATING COVID ABDAB DATABASE BY JUST GRABBING RAW DATA...


#### NEED TO RELOAD ALL TAB FILES BEFORE UPDATING
toshiny.few1 <- read_tsv("/Users/eric.waltari/data_carpentry/wikipathways/App-6/toshiny_few1.tab")
## old 
toshiny.few1 <- read_tsv("/Users/eric.waltari/data_carpentry/wikipathways/app6_fromjan2021/toshiny_few1.tab")

toshiny.few <- read_tsv("toshiny_few.tab")
toshiny.few1 <- read_tsv("toshiny_few1.tab")
toshiny.few1.h <- read_tsv("toshiny_few1_h.tab")
## realoading few2 - will want to combine COMET for internal comparisons (moving to app_test in _2 folder)
toshiny.few2 <- read_tsv("toshiny_few2b.tab")

toshiny.few3 <- read_tsv("toshiny_few3.tab")
# toshiny.few3b <- read_tsv("toshiny_few3b.tab")
toshiny.fewhcbg <- read_tsv("toshiny_fewhcbg.tab")
toshiny.fewhcbgc <- read_tsv("toshiny_fewhcbgc.tab")
# toshiny.natalia4 <- read_tsv("toshiny_natalia.tab")
toshiny.cov2sarsmers <- read_tsv("toshiny_cov2sarsmers.tab")
toshiny.cov2sarsmersc <- read_tsv("toshiny_cov2sarsmersc.tab")

toshiny.hivcov2sarsmers <- read_tsv("toshiny_hivcov2sarsmers.tab")
toshiny.hivcov2sarsmersc <- read_tsv("toshiny_hivcov2sarsmersc.tab")

toshiny.fewhivmabsnih45rd1214 <- read_tsv("toshiny_mabsnih45rd1214.tab")
toshiny.fewhivmabsnih45rd1214c <- read_tsv("toshiny_mabsnih45rd1214c.tab")

toshiny.dengueall <- read_tsv("toshiny_dengueall.tab")
toshiny.dengueallc <- read_tsv("toshiny_dengueallc.tab")

toshiny.fluhcmabs <- read_tsv("toshiny_fluhcmabs.tab")
toshiny.fluhcmabsc <- read_tsv("toshiny_fluhcmabsc.tab")

toshiny.fewg1 <- read_tsv("toshiny_fewg1.tab")
toshiny.fewb1 <- read_tsv("toshiny_fewb1.tab")
toshiny.few3.h <- subset(toshiny.few3, cregion %in% c("IgH"))



## FIRST TAKE DATA AND ADD TO GOOGLE DRIVE SPREADSHEET
## THEN RUN A JOIN COMMAND - better yet an rbind if all columns are fixed

### latest version of this is from end of October - starting in early November adding directly to toshiny

#BX.clusters.mabs.toadd.h <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_111020add.tab")
#BX.clusters.mabs.toadd.h <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_251120add.tab")
# BX.clusters.mabs.toadd.h <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_060121add.tab")
# BX.clusters.mabs.toadd.l <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_060121add.tab")

BX.clusters.mabs.toadd.h <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_090721add.tab")
BX.clusters.mabs.toadd.l <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/COVID_mablist_090721add.tab")



BX.clusters.mabs.toadd.h$binding <- gsub("S; ","",BX.clusters.mabs.toadd.h$binding)
BX.clusters.mabs.toadd.h$neutralization <- gsub("non-neutralizingunknown","non-neutralizing or unknown",BX.clusters.mabs.toadd.h$neutralization)

BX.clusters.mabs.toadd.h$sequence_idLC <- NULL
BX.clusters.mabs.toadd.h$geneLC <- NULL
BX.clusters.mabs.toadd.h$jgeneLC <- NULL

BX.clusters.mabs.toadd.h$junction_aaLC <- NULL
BX.clusters.mabs.toadd.h$v_identityLC <- NULL

#BX.clusters.mabs.toadd.h$FULLV_LC <- NULL

BX.clusters.mabs.toadd.h$gene <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.h$gene)
BX.clusters.mabs.toadd.h$jgene <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.h$jgene)



BX.clusters.mabs.toadd.l$binding <- gsub("S; ","",BX.clusters.mabs.toadd.l$binding)
BX.clusters.mabs.toadd.h$neutralization <- gsub("non-neutralizingunknown","non-neutralizing or unknown",BX.clusters.mabs.toadd.h$neutralization)

BX.clusters.mabs.toadd.l$sequence_id <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, sequence_id = sequence_idLC)


BX.clusters.mabs.toadd.l$gene <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, gene = geneLC)
BX.clusters.mabs.toadd.l$jgene <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, jgene = jgeneLC)

# BX.clusters.mabs.toadd.l$FULLV <- NULL
# BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, FULLV = FULLV_LC)


BX.clusters.mabs.toadd.l$junction_aa <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, junction_aa = junction_aaLC)
BX.clusters.mabs.toadd.l$v_identity <- NULL
BX.clusters.mabs.toadd.l <- rename(BX.clusters.mabs.toadd.l, v_identity = v_identityLC)

BX.clusters.mabs.toadd.l$gene <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.l$gene)
BX.clusters.mabs.toadd.l$jgene <- gsub(" [(]Human[)]","",BX.clusters.mabs.toadd.l$jgene)

BX.clusters.mabs.toadd <- rbind(BX.clusters.mabs.toadd.h, BX.clusters.mabs.toadd.l)

## THIS NEXT COMMAND MIGHT NOT BE THE CASE FOR ALL - SHOULD CHANGE IN .TAB FILE
#BX.clusters.mabs.toadd$neutralization <- "neutralizing"


#### #### #### #### #### #### #### #### #### #### #### #### 
## IN JULY 2021 ADDING ALL MABS OVER AGAIN, BUT THEN GRABBING shm IF IT EXISTS...THEN WILL GO THROUGH MISSING TO FIND ANY WITH COMPLETE VgeneS
#### THEN SEARCH IN geneIOUS FOR MATCHES TO GET MORE shms


toshiny.few3 <- read_tsv("toshiny_few3.tab")
toshiny.few3.h <- subset(toshiny.few3, cregion %in% c("IgH"))

toshiny.fewg1 <- read_tsv("toshiny_fewg1.tab")
toshiny.fewb1 <- read_tsv("toshiny_fewb1.tab")


toshiny.few1.old <- toshiny.few1 %>% select(sequence_id,shm)

BX.clusters.mabs.toadd2 <- left_join(BX.clusters.mabs.toadd, toshiny.few1.old)
#write.table(BX.clusters.mabs.toadd2, "BX_clusters_mabs_toadd2.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## next in Geneious ran tblastn searches to find 100% matching nucleotide sequences, then internal IgBlasts to get shm

BX.clusters.mabs.toadd3 <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/galson_data/BX_clusters_mabs_toadd3.tab")

BX.clusters.mabs.toadd2b <- left_join(BX.clusters.mabs.toadd2, BX.clusters.mabs.toadd3, by = "sequence_id")
BX.clusters.mabs.toadd2b <- BX.clusters.mabs.toadd2b %>% unite(shm, shm.x, shm.y, na.rm = TRUE, remove = TRUE)
BX.clusters.mabs.toadd2b$shm <- as.numeric(BX.clusters.mabs.toadd2b$shm)
BX.clusters.mabs.toadd2b$FULLV <- NULL
BX.clusters.mabs.toadd2b$v_identity <- NULL

BX.clusters.mabs.toadd <- BX.clusters.mabs.toadd2b
## NOW CONTINUTE WITH CODE BELOW....EXCEPT NOT shm CALCULATION, ALREADY HAVE IT!

rm(BX.clusters.mabs.toadd2)
rm(BX.clusters.mabs.toadd2b)
rm(BX.clusters.mabs.toadd3)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### 

##########

BX.clusters.mabs.toadd$junction_aa <- as.character(BX.clusters.mabs.toadd$junction_aa)
BX.clusters.mabs.toadd$sequence_id <- as.character(BX.clusters.mabs.toadd$sequence_id)

BX.clusters.mabs.toadd$cdr3length_imgt <- nchar(BX.clusters.mabs.toadd$junction_aa)

BX.clusters.mabs.toadd$junction_length <- (3 * (BX.clusters.mabs.toadd$cdr3length_imgt)) + 6
BX.clusters.mabs.toadd$jgf <- BX.clusters.mabs.toadd$jgene
BX.clusters.mabs.toadd$gf <- str_sub(BX.clusters.mabs.toadd$gene, end=5)
### note updated str_sub for substring (think both work)
# BX.clusters.mabs.toadd$gf <- substring(BX.clusters.mabs.toadd$gene, 1,5)
# BX.clusters.mabs.toadd$jgf <- substring(BX.clusters.mabs.toadd$jgene, 1,5)
BX.clusters.mabs.toadd$vj_junction_pattern <- paste(BX.clusters.mabs.toadd$gf,BX.clusters.mabs.toadd$jgf,BX.clusters.mabs.toadd$junction_length,sep="_") 

### check v_identity IF 0-100 OR 0-1 - WE WANT shm TO BE BETWEEN 0-100
BX.clusters.mabs.toadd$shm <- (100 - (BX.clusters.mabs.toadd$v_identity))
# BX.clusters.mabs.toadd$shm <- (100 - (BX.clusters.mabs.toadd$v_identity * 100))
BX.clusters.mabs.toadd <- BX.clusters.mabs.toadd %>% unite(gf_jgene, gf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)

BX.clusters.mabs.toadd$junction_length <- NULL
BX.clusters.mabs.toadd$v_identity <- NULL

## ADDING cregion BY WAY OF gf...
BX.clusters.mabs.toadd$cregion <- str_sub(BX.clusters.mabs.toadd$gf, end=3)
BX.clusters.mabs.toadd$cregion <- gsub("IGK","Kappa",BX.clusters.mabs.toadd$cregion)
BX.clusters.mabs.toadd$cregion <- gsub("IGL","Lambda",BX.clusters.mabs.toadd$cregion)
BX.clusters.mabs.toadd$cregion <- gsub("IGH","IgH",BX.clusters.mabs.toadd$cregion)



### ALSO JUNE 2021 - NEED TO REMOVE BAD shm DATA FROM CVXXX MABS IN toshiny.few1, then recalculate all of the datasets...
### also manually removing from first tab file

toshiny.few1 <- read_tsv("toshiny_few1.tab")

BX.clusters.mabs.alladded <- full_join(toshiny.few1, BX.clusters.mabs.toadd)
### now this can replace toshiny.few1 
rm(BX.clusters.mabs.toadd)
rm(BX.clusters.mabs.toadd.h)
rm(BX.clusters.mabs.toadd.l)

#BX.clusters.fewermabs <- BX.clusters.mabs %>% filter(source != "Augmenta")

toshiny.few1$binding <- gsub("Unk","unknown",toshiny.few1$binding)

toshiny.few1 <- BX.clusters.mabs.alladded %>% select(sequence_id,binding,neutralization,cregion,junction_aa,gene,gf_jgene,gf,jgene,jgf,cdr3length_imgt,shm,ncount,shm_mean,reads_per_clone) %>% filter(!is.na(cdr3length_imgt)) %>% filter(is.wholenumber(cdr3length_imgt)) %>% mutate(across(shm, round, 2)) %>% mutate(across(shm_mean, round, 2))
## ONLY JULY 2021 RE-RUN
#toshiny.few1 <- BX.clusters.mabs.toadd %>% select(sequence_id,binding,neutralization,cregion,junction_aa,gene,gf_jgene,gf,jgene,jgf,cdr3length_imgt,shm) %>% filter(!is.na(cdr3length_imgt)) %>% mutate(across(shm, round, 2))


### finding some errors in the file jul27
## C215 has lost gene family also IG
## 2 mabs have for binding "NTD" and "non-S1" - should be changed to "non-RBD"
#toshiny.few1 <- read_tsv("toshiny_few1.tab")


## NOW NEED TO RECALCULATE NCOUNT AND shm_mean
toshiny.few1$ncount <- NULL
toshiny.few1$shm_mean <- NULL

toshiny.few1 <- toshiny.few1 %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))

toshiny.few1.h <- subset(toshiny.few1, cregion %in% c("IgH"))
# toshiny.few1.h$binding <- gsub("Unk","unknown",toshiny.few1.h$binding)

#write.table(toshiny.few1, "toshiny_few1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.few1.h, "toshiny_few1_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)

rm(BX.clusters.mabs.alladded)

## NOTE THAT THIS WILL NOT HAVE shm DATA (UNLESS FOUND IN A PAPER) - WOULD HAVE TO SEPARATELY GO BACK AND FIND RAW NUCLEOTIDE SEQUENCE IN NCBI TO CALCULATE shm...



#######################################
### NEXT: re-calculating all of the combinations...

### 1 toshiny.few
## note HC only is few3 aka .hc - few3b is HC but all reads not just clones
toshiny.few3.h <- subset(toshiny.few3, cregion %in% c("IgH"))

toshiny.few <- bind_rows(toshiny.few1.h, toshiny.few3.h, .id = "id")

toshiny.few$id <- gsub("2","Healthy control bulk repertoire",toshiny.few$id)
toshiny.few$id <- gsub("1","anti-CoV2 mAbs",toshiny.few$id)

#write.table(toshiny.few, "toshiny_few.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### 2


BX.clusters.hcvsboydvsgalsonmabs <- bind_rows(toshiny.few3.h, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, .id = "id")
#toshiny.fewhcbg <- subset(BX.clusters.hcvsboydvsgalsonmabs, cregion %in% c("IgH"))
toshiny.fewhcbg <- BX.clusters.hcvsboydvsgalsonmabs[ grep("Kappa|Lambda", BX.clusters.hcvsboydvsgalsonmabs$cregion, invert = TRUE) , ] %>% filter(!is.na(cdr3length_imgt)) %>% filter(is.integer(cdr3length_imgt)) %>% mutate(across(shm, round, 2)) %>% mutate(across(shm_mean, round, 2))
#rm(BX.clusters.hcvsboydvsgalsonmabs)

toshiny.fewhcbg$id <- gsub("1","Healthy control",toshiny.fewhcbg$id)
## be careful with 1 in name here....
toshiny.fewhcbg$id <- gsub("4","anti-CoVII mAbs",toshiny.fewhcbg$id)
toshiny.fewhcbg$id <- gsub("2","Boyd7450",toshiny.fewhcbg$id)
toshiny.fewhcbg$id <- gsub("3","Galson1",toshiny.fewhcbg$id)
toshiny.fewhcbg$id <- gsub("anti-CoVII mAbs","anti-CoV2 mAbs",toshiny.fewhcbg$id)

## for combined need to recalculate ncount and shm_mean
toshiny.fewhcbgc <- toshiny.fewhcbg
toshiny.fewhcbgc$ncount <- NULL
toshiny.fewhcbgc$shm_mean <- NULL

toshiny.fewhcbgc <- toshiny.fewhcbgc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))

#write.table(toshiny.fewhcbg, "toshiny_fewhcbg.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.fewhcbgc, "toshiny_fewhcbgc.tab", sep = "\t", row.names = FALSE, quote = FALSE)



### 3

## need to also load these:
#toshiny.sars.h, toshiny.mers.h, toshiny.fewhivmabs.h, toshiny.fewhiv

toshiny.fewhiv <- read_tsv("toshiny_fewhiv.tab")

toshiny.fewhivmabs <- read_tsv("toshiny_fewhivmabs.tab")

toshiny.fewhivmabs.h <- subset(toshiny.fewhivmabs, cregion %in% c("IgH"))
#toshiny.fewhivmabs.h <- read_tsv("toshiny_fewhivmabs_h.tab")
#toshiny.fewhiv <- read_tsv("toshiny_fewhiv.tab")

toshiny.sarsmers.h <- read_tsv("toshiny_sarsmers_h.tab")
list2env(x = split(x = toshiny.sarsmers.h,
                   f = toshiny.sarsmers.h$source),
         envir = globalenv())

toshiny.sars.h <- SARS
toshiny.mers.h <- MERS


BX.clusters.hcvsboydvsgalsonmabshiv <- bind_rows(toshiny.few3.h, toshiny.fewb1, toshiny.fewg1, toshiny.few1.h, toshiny.fewhiv, .id = "id")
#toshiny.fewhcbghiv <- subset(BX.clusters.hcvsboydvsgalsonmabshiv, cregion %in% c("IgH"))
toshiny.fewhcbghiv <- BX.clusters.hcvsboydvsgalsonmabshiv[ grep("Kappa|Lambda", BX.clusters.hcvsboydvsgalsonmabshiv$cregion, invert = TRUE) , ] %>% filter(!is.na(cdr3length_imgt)) %>% filter(is.integer(cdr3length_imgt)) %>% mutate(across(shm, round, 2)) %>% mutate(across(shm_mean, round, 2))
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

## for combined need to recalculate ncount and shm_mean
toshiny.fewhcbghivc <- toshiny.fewhcbghiv
toshiny.fewhcbghivc$ncount <- NULL
toshiny.fewhcbghivc$shm_mean <- NULL

toshiny.fewhcbghivc <- toshiny.fewhcbghivc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))

## updaating HIV+ patient to HIV+ patient MT1214
toshiny.fewhcbghiv$id <- gsub("HIV+ patient","HIV+ patient MT1214", toshiny.fewhcbghiv$id)
toshiny.fewhcbghivc$id <- gsub("HIV+ patient","HIV+ patient MT1214", toshiny.fewhcbghivc$id)


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

## for combined need to recalculate ncount and shm_mean
toshiny.cov2sarsmersc <- toshiny.cov2sarsmers
toshiny.cov2sarsmersc$ncount <- NULL
toshiny.cov2sarsmersc$shm_mean <- NULL

toshiny.cov2sarsmersc <- toshiny.cov2sarsmersc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))

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

## for combined need to recalculate ncount and shm_mean
toshiny.hivcov2sarsmersc <- toshiny.hivcov2sarsmers
toshiny.hivcov2sarsmersc$ncount <- NULL
toshiny.hivcov2sarsmersc$shm_mean <- NULL

toshiny.hivcov2sarsmersc <- toshiny.hivcov2sarsmersc %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))

#write.table(toshiny.hivcov2sarsmers, "toshiny_hivcov2sarsmers.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.hivcov2sarsmersc, "toshiny_hivcov2sarsmersc.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## then putting into app-6

#install.packages("rsconnect")
library(rsconnect)
rsconnect::deployApp('/Users/eric.waltari/data_carpentry/wikipathways/App-6')

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#### jan 2021 ideas:
## change plot title to name of V+J+CDR3length?? key is here:
#filteredData <- isolate({filteredDSpartial2()}) - note this line is in 2 of 3 options, would need to add to third
## then near bottom of plot code something similar to this line to add to actual title:
# gjcdr3.title <- filteredData %>% separate(sequence_id, into = c("gf_jgene", "cdr3length_imgt"), sep = "_", remove = FALSE, convert = TRUE, extra = "merge", fill = "left") %>%
#   select(sequence_id)

filteredData$cdr3length_imgt <- as.numeric(filteredData$cdr3length_imgt)
filteredData$G_J_CDR3 <- paste(filteredData$gf_jgene,filteredData$cdr3length_imgt,sep="_")
gjcdr3.title <- filteredData %>%
  select(G_J_CDR3)
gjcdr3.titleID <- gjcdr3.title$G_J_CDR3[1]
## then use gjcdr3.titleID to replace current title in " "

#################################
### END COMMANDS FOR UPDATING COVID ABDAB DATABASE BY JUST GRABBING RAW DATA...


##################################################################################################################
##################################################################################################################
##################################################################################################################

########
## reading in previous healthy germ-pass samples to get clonalclusters

################################################################################################################################
#######################################################################################################################################



#### to reload all correct tab files per app-2 and app-5
toshiny.few <- read_tsv("toshiny_few.tab")
toshiny.few1 <- read_tsv("toshiny_few1.tab")
toshiny.few1.h <- read_tsv("toshiny_few1_h.tab")
toshiny.few2 <- read_tsv("toshiny_few2.tab")
toshiny.few3 <- read_tsv("toshiny_few3.tab")
toshiny.few3b <- read_tsv("toshiny_few3b.tab")
toshiny.fewhcbg <- read_tsv("toshiny_fewhcbg.tab")
toshiny.fewhcbgc <- read_tsv("toshiny_fewhcbgc.tab")
toshiny.fewhcbghiv <- read_tsv("toshiny_fewhcbghiv.tab")
toshiny.fewhcbghivc <- read_tsv("toshiny_fewhcbghivc.tab")
toshiny.natalia4 <- read_tsv("toshiny_natalia.tab")

toshiny.fewer <- read_tsv("toshiny_fewer.tab")
toshiny.fewer1 <- read_tsv("toshiny_fewer1.tab")
toshiny.fewer1.h <- read_tsv("toshiny_fewer1_h.tab")
# toshiny.few2 <- read_tsv("toshiny_few2.tab")
# toshiny.few3 <- read_tsv("toshiny_few3.tab")
toshiny.fewerhcbg <- read_tsv("toshiny_fewerhcbg.tab")
toshiny.fewerhcbgc <- read_tsv("toshiny_fewerhcbgc.tab")


########################################## RD1214 HIV dataset

BX.clusters.hiv <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/binder_data/MT1214HC_germ-pass_clones.tab")

## unique to this dataset
BX.clusters.hiv <- BX.clusters.hiv %>% filter(junction_length > 4)

BX.clusters.hiv$v_call <- as.character(BX.clusters.hiv$v_call)
BX.clusters.hiv$j_call <- as.character(BX.clusters.hiv$j_call)

BX.clusters.hiv$gene <- getGene(BX.clusters.hiv$v_call, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$gf <- substring(BX.clusters.hiv$gene, 1,5)
BX.clusters.hiv$jgene <- getGene(BX.clusters.hiv$j_call, first=TRUE, strip_d=TRUE)
BX.clusters.hiv$jgf <- substring(BX.clusters.hiv$jgene, 1,5)
BX.clusters.hiv$vj_junction_pattern <- paste(BX.clusters.hiv$gf,BX.clusters.hiv$jgf,BX.clusters.hiv$junction_length,sep="_") 

BX.clusters.hiv$cdr3length_imgt <- ((BX.clusters.hiv$junction_length) / 3) - 2
#BX.clusters.hiv$jgene[BX.clusters.hiv$jgene==""]<-NA
BX.clusters.hiv$shm <- (100 - (BX.clusters.hiv$v_identity * 100))
BX.clusters.hiv$junction_aa <- translateDNA(BX.clusters.hiv$junction, trim=TRUE)

BX.clusters.hiv <- BX.clusters.hiv %>% filter(cdr3length_imgt > 4.8)
#write.table(BX.clusters.hiv, "MT1214HC_germ-pass_clones2.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#BX.clusters.hiv <- BX.clusters.hiv[ grep("*", BX.clusters.hiv$junction_aa, invert = TRUE) , ]
#BX.clusters.hiv$jgene[BX.clusters.hiv$jgene==""]<-NA
#BX.clusters.hiv <- read_tsv("/users/eric.waltari/immcantation_pipeline/COVID_mabs/binder_data/MT1214HC_germ-pass_clones2.tab")

BX.clusters.hiv <- BX.clusters.hiv %>% unite(gf_jgene, gf, jgene, sep = "_", remove = FALSE, na.rm = TRUE) %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE))

BX.clusters.hiv <- BX.clusters.hiv %>% filter(is.wholenumber(cdr3length_imgt))
BX.clusters.hiv$cregion <- "IgH"


toshiny.fewhiv <- BX.clusters.hiv %>% select(sequence_id,cregion,junction_aa,gene,gf_jgene,gf,jgene,jgf,vj_junction_pattern,cdr3length_imgt,shm,ncount,shm_mean) %>% filter(!is.na(cdr3length_imgt)) %>% filter(is.wholenumber(cdr3length_imgt)) %>% mutate(across(shm, round, 2)) %>% mutate(across(shm_mean, round, 2))

toshiny.fewhiv <- toshiny.fewhiv[ grep("Kappa|Lambda|IGKJ|__", toshiny.fewhiv$jgene, invert = TRUE) , ]
#write.table(toshiny.fewhiv, "toshiny_fewhiv.tab", sep = "\t", row.names = FALSE, quote = FALSE)


## combining this hiv instead of previous hiv run...actually add both!
toshiny.fewhiv1 <- read_tsv("toshiny_fewhiv1.tab")
toshiny.fewhiv1$ncount <- NULL
toshiny.fewhiv1$shm_mean <- NULL

toshiny.fewhiv1 <- toshiny.fewhiv1 %>%
  add_count(gf_jgene,cdr3length_imgt) %>%
  rename(ncount = n) %>%
  group_by(gf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  mutate(across(shm_mean, round, 2))


########################################
### sept2021 going back to raw dataset - instead of taking 1 sequence per clone, removing all singleton clones
BX.clusters.rd1214 <- read_tsv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/MT1214HC_germ-pass.tab")

# BX.clusters.rd1214 <- BX.clusters.rd1214 %>% filter(FUNCTIONAL == "T") %>% 
#   filter(IN_FRAME == "T") %>% 
#   filter(STOP == "F")
#toshiny.fewhiv <- toshiny.fewhiv[ grep("Kappa|Lambda|IGKJ|__", toshiny.fewhiv$jgene, invert = TRUE) , ]

BX.clusters.rd1214$cregion <- str_sub(BX.clusters.rd1214$v_call, end=3)
BX.clusters.rd1214 <- subset(BX.clusters.rd1214, cregion %in% c("IgH"))

BX.clusters.rd1214fewer <- BX.clusters.rd1214 %>% add_count(clone_id) %>%
  rename(clonecount = n) %>%
  filter(clonecount > 1)

## THIS IS A ONE-OFF BECAUSE NOT USING LATEST AIRR FORMAT
## but will need similar code for all datasets Sept 2021
BX.clusters.rd1214fewer$junction_aa <- translateDNA(BX.clusters.rd1214fewer$junction, trim=TRUE)
# str_sub(BX.clusters.rd1214fewer$junction_aa, 1, 1) <- ""
# str_sub(BX.clusters.rd1214fewer$junction_aa, -1, -1) <- ""
BX.clusters.rd1214fewer$cdr3length_imgt <- nchar(BX.clusters.rd1214fewer$junction_aa)

BX.clusters.rd1214fewer$junction_length_CHECK <- BX.clusters.rd1214fewer$junction_length / 3
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer %>% filter(is.wholenumber(junction_length_CHECK))

## FILTER junction_aa > 3 & NO NAs...
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer %>% filter(!is.na(junction_aa)) %>% filter(cdr3length_imgt > 2)
## ALSO NO *
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer[ grep("\\*", BX.clusters.rd1214fewer$junction_aa, invert = TRUE) , ]


BX.clusters.rd1214fewer$gene <- getGene(BX.clusters.rd1214fewer$GERMLINE_v_call, first=TRUE, strip_d=TRUE)
BX.clusters.rd1214fewer$gf <- substring(BX.clusters.rd1214fewer$gene, 1,5)
BX.clusters.rd1214fewer$jgene <- getGene(BX.clusters.rd1214fewer$GERMLINE_j_call, first=TRUE, strip_d=TRUE)
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer %>% unite(gf_jgene, gf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)
## some IGK & IGL
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer[ grep("IGK|IGL", BX.clusters.rd1214fewer$gf_jgene, invert = TRUE) , ]
BX.clusters.rd1214fewer <- BX.clusters.rd1214fewer[ grep("kappa|lambda", BX.clusters.rd1214fewer$PRIMER, invert = TRUE) , ]


## THEN GROUP BY V CALL JCALL AND junction_aa, SELECT ONLY ONE PER GROUP

BX.clusters.rd1214fewest <- BX.clusters.rd1214fewer %>%
  group_by(junction_aa,gf_jgene) %>%
  summarize_all(first)
## should be the same! note updating superceded summarize_all per dplyr help...
# BX.clusters.rd1214fewest <- BX.clusters.rd1214fewer %>%
#   group_by(junction_aa,gf_jgene) %>%
#   summarize(first)

h3.hc <- ggplot(BX.clusters.rd1214fewest, aes(gf_jgene,cdr3length_imgt))

BX.clusters.rd1214fewest.viridis <- h3.hc + geom_bin2d(aes(fill=log10(..count..))) + theme_bw() + ylab("CDRH3 Length (aa)") + xlab("V-gene & J-gene") + facet_wrap(~ cregion, ncol=1, scales = "free") + scale_fill_viridis_c(name = "# of \nReads", option = "D",  breaks = c(0, 1, 2, 3, 4), labels = c(1, 10, 100, 1000, 10000)) + theme(axis.text.x = element_text(angle=45, hjust=1, size=5))
BX.clusters.rd1214fewest.viridis
#write.table(BX.clusters.rd1214fewest, "BX_clusters_rd1214fewest.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### I THEN CHANGED HEADERS TO MATCH AIRR FORMAT, ALSO CHANGED NAME TO RD1214HC_germ-pass2.tab

### sept 13 trying to recapitulate RD1214 from OAS download:
# first need to remove first line from every file, then import into R

BX.clusters.rd1214a <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/SRR5811762_Heavy_IGHAn.csv")
BX.clusters.rd1214a$ANARCI_numbering <- NULL
BX.clusters.rd1214a$ANARCI_status <- NULL
BX.clusters.rd1214a <- BX.clusters.rd1214a %>% filter(Redundancy > 1)


BX.clusters.rd1214d <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/SRR5811762_Heavy_IGHDn.csv")
BX.clusters.rd1214e <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/SRR5811762_Heavy_IGHEn.csv")
BX.clusters.rd1214g <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/SRR5811762_Heavy_IGHGn.csv")

BX.clusters.rd1214d$ANARCI_numbering <- NULL
BX.clusters.rd1214d$ANARCI_status <- NULL
BX.clusters.rd1214d <- BX.clusters.rd1214d %>% filter(Redundancy > 1)

BX.clusters.rd1214e$ANARCI_numbering <- NULL
BX.clusters.rd1214e$ANARCI_status <- NULL
BX.clusters.rd1214e <- BX.clusters.rd1214e %>% filter(Redundancy > 1)


BX.clusters.rd1214g$ANARCI_numbering <- NULL
BX.clusters.rd1214g$ANARCI_status <- NULL
BX.clusters.rd1214g <- BX.clusters.rd1214g %>% filter(Redundancy > 1)


BX.clusters.rd1214m <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/SRR5811762_Heavy_IGHMn.csv")
BX.clusters.rd1214m$ANARCI_numbering <- NULL
BX.clusters.rd1214m$ANARCI_status <- NULL
BX.clusters.rd1214m <- BX.clusters.rd1214m %>% filter(Redundancy > 1)

BX.clusters.rd1214b <- read.csv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/sept21/SRR5811762_Heavy_Bulkn.csv")
BX.clusters.rd1214b$ANARCI_numbering <- NULL
BX.clusters.rd1214b$ANARCI_status <- NULL
BX.clusters.rd1214b <- BX.clusters.rd1214b %>% filter(Redundancy > 1)


BX.clusters.rd1214a$cregion <- "IgA"
BX.clusters.rd1214g$cregion <- "IgG"
BX.clusters.rd1214m$cregion <- "IgM"

## combine, then locus IGH
BX.clusters.rd1214agm <- rbind(BX.clusters.rd1214a, BX.clusters.rd1214g, BX.clusters.rd1214m)
BX.clusters.rd1214all <- full_join(BX.clusters.rd1214agm, BX.clusters.rd1214b)
BX.clusters.rd1214all$locus <- "IGH"

## RUN THESE ON ALL DATASETS TO BE USED - ALSO LOOK INTO RENUMBERING WITHIN R
BX.clusters.rd1214all$junction_length_check <- BX.clusters.rd1214all$junction_length / 3
BX.clusters.rd1214all <- BX.clusters.rd1214all %>% filter(is.wholenumber(junction_length_check))

BX.clusters.rd1214all$dataset <- "RD1214"
BX.clusters.rd1214all$obs <- 1:nrow(BX.clusters.rd1214all) 
BX.clusters.rd1214all <- BX.clusters.rd1214all %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
  
#write.table(BX.clusters.rd1214all, "BX_clusters_rd1214downloaded.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### code to thin OAS & BXmay hc datasets:

BX.clusters.OASo <- read_tsv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/OAS_sept21_db-pass.tsv")

## filter productive = FALSE, * in junction_aa
BX.clusters.OAS <- BX.clusters.OASo %>% filter(productive == TRUE) %>% filter(productive == TRUE)
BX.clusters.OAS <- BX.clusters.OAS[ grep("\\*", BX.clusters.OAS$junction_aa, invert = TRUE) , ]
## 700k to 363k 
#write.table(BX.clusters.OAS, "OAS_sept21_db-pass2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
rm(BX.clusters.OASo)


### NEW CODE FOR IMPORTING NEW OAS DENGUE DATA
## SEE AT VERY BOTTOM

## final concatenate command
OAS.clusters.all <- rbind(OAS.clusters.p148,OAS.clusters.p172,OAS.clusters.p194,OAS.clusters.p199,OAS.clusters.p203,OAS.clusters.p208,OAS.clusters.p232,OAS.clusters.p237,OAS.clusters.p238,OAS.clusters.p240,OAS.clusters.p249,OAS.clusters.p252,OAS.clusters.p255,OAS.clusters.p265,OAS.clusters.p275,OAS.clusters.p276,OAS.clusters.p287,OAS.clusters.p289,OAS.clusters.p299,OAS.clusters.p301,OAS.clusters.p307,OAS.clusters.p311,OAS.clusters.p320,OAS.clusters.p346,OAS.clusters.p376,OAS.clusters.p391,OAS.clusters.p422,OAS.clusters.p444,OAS.clusters.p455,OAS.clusters.p479,OAS.clusters.p481,OAS.clusters.p489,OAS.clusters.p500,OAS.clusters.p514,OAS.clusters.p515,OAS.clusters.p517,OAS.clusters.p520,OAS.clusters.p524,OAS.clusters.p529,OAS.clusters.p543,OAS.clusters.p551,OAS.clusters.p555,OAS.clusters.p558,OAS.clusters.p563,OAS.clusters.p569)
OAS.clusters.all$ANARCI_numbering <- NULL
OAS.clusters.all$ANARCI_status <- NULL

#write.table(OAS.clusters.all, "OAS_sept21_germ-pass2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###

BX.clusters.BXmayo <- read_tsv("/Users/eric.waltari/immcantation_pipeline/COVID_mabs/BXmay10mstim_db-passHC.tsv")

BX.clusters.BXmay <- BX.clusters.BXmayo %>% filter(consensus_count > 2)
## 353k to 89k

#write.table(BX.clusters.BXmay, "BXmay10mstim_db-pass2HC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



### end sept2021

# dengue
### LATE NOV 2020 - NEED TO COLLAPSE IDENTICAL junction_aa FIRST (ALSO FOR OAS), OTHERWISE TOO MANY SEQUENCES...
# distinct_at() will be superseded by use of across() inside distinct() from version 1.0.0 (dplyr.tidyverse.org/news/index.html#across). 
# The equivalent pattern for your answer would be dat %>% distinct(across(-z)) but distinct_at() will still be available for several years. 
## to make sure it keeps largest consensus_count though - group_by(consensus_count) %>%
BX.clusters.dengue.d13.2$junction_aa <- translateDNA(BX.clusters.dengue.d13.2$junction, trim=TRUE)
# BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2 %>% distinct(across(junction_aa), .keep_all = TRUE)


## MAR 2021 - BEFORE COMBINING BELOW, LOOK AT THE DATASETS SEPARATELY...

## NEXT FIND THE CONVERGENT SEQUENCES AND REMOVE...
BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2[ grep("ARALFGLVAVASPFDN", BX.clusters.dengue.d13.2$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2[ grep("ITPPLYLMVGGVSRAMAV", BX.clusters.dengue.d13.2$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d13.2 <- BX.clusters.dengue.d13.2[ grep("ARQDRNWFDT", BX.clusters.dengue.d13.2$junction_aa, invert = TRUE) , ]

BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2[ grep("ARGPGGTSTSCYHCWFDP", BX.clusters.dengue.d20.2$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d20.2 <- BX.clusters.dengue.d20.2[ grep("AKNYGSGTLNWFDS", BX.clusters.dengue.d20.2$junction_aa, invert = TRUE) , ]


BX.clusters.dengue.d13.2stim <- BX.clusters.dengue.d13.2stim[ grep("DC38", BX.clusters.dengue.d13.2stim$sequence_id, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("DC38", BX.clusters.dengue.d20.2stim$sequence_id, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARADEMATVQgfYAFDI", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARADEMATIEgfYAFDI", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARADEMATIEgfYAFGI", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]

BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARGPGGTTTSCYHCWFDP", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARGPGGTSSSCYQCWFDP", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("AKNYGSGTLNWFDS", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]

BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("ARgfATTQWQGHNWFDP", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]
BX.clusters.dengue.d20.2stim <- BX.clusters.dengue.d20.2stim[ grep("AKDVGECSGGNCFSGYFYYMDA", BX.clusters.dengue.d20.2stim$junction_aa, invert = TRUE) , ]

######################

## getting CF10 & CF7 (also CF4) extractions

#toshiny.dengue.cf10 <- subset(toshiny.dengue.mabs, cregion %in% c("IgG"))
toshiny.dengue.cf10 <- toshiny.dengueallmore %>% filter(gf_jgene == "IGHV4_IGHJ5" & cdr3length_imgt == "10")
toshiny.dengue.cf7 <- toshiny.dengueallmore %>% filter(gf_jgene == "IGHV1_IGHJ5" & cdr3length_imgt == "16")

#write.table(toshiny.dengue.cf10, "toshiny_dengue_cf10.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.dengue.cf7, "toshiny_dengue_cf4andcf7.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#starwars <- starwars %>% filter(hair_color == "none" & eye_color == "black")




## NEED TO THIN TO SMALLER SUBSET THOUGH...inner join with subset from Geneious
##  dengue_CF10_list.tsv
CF10mergedwithOAS.list <- read_tsv("dengue_CF10_list.tsv")

## next run some sort of join...semi-join semi_join() return all rows from x with a match in y.

CF10mergedwithOASfewer <- semi_join(toshiny.dengue.cf10, CF10mergedwithOAS.list)

CF4mergedwithOAS.list <- read_tsv("dengue_CF4_list.tsv")

CF4mergedwithOASfewer <- semi_join(toshiny.dengue.cf7, CF4mergedwithOAS.list)

### networks of dengue cf10 & cf 7 - see far below 


##############################################################################################################################
##############################################################################################################################

### NEW CODE FOR IMPORTING NEW OAS DENGUE DATA

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

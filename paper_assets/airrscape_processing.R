
## if need to install
# pkgs = c("alakazam", "igraph", "dplyr","RColorBrewer", "hexbin", "scales","grid", "lattice", "gdata","gridExtra", "ape", "shazam","reshape2", "DT","ggplot2", "seqinr", "phangorn","shiny","rlang", "knitr", tidyverse") # package names
# install.packages(pkgs)

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
library(reshape2)
library(DT)

library(ggplot2)
library(seqinr)
library(phangorn)
library(shiny)
library(rlang)
library(knitr)
library(tidyverse)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

## this works ONLY IF YOU HAVE CREATED THE DIRECTORY FIRST...
## adjust the directory structure in the next command to whereever you have downloaded the AIRRscape repo:
#setwd("~/data_carpentry/AIRRscape")


### LOADING AND CONVERTING INTERMEDIATE DATASETS
cov2.bulk.binder.p11 <- read_tsv("paper_assets/intermediate_files/Binder_p11_germ-pass.tsv.gz")   # v_identity between 0.60 and 1, sequences renamed, have clone_id & germline_alignment_d_mask, NOT cdr3_aa
cov2.bulk.nielsen.p7450 <- read_tsv("paper_assets/intermediate_files/Nielsen_7450_airr-covid-19.tsv.gz")   # v_identity between 60 and 100, has cdr3_aa
cov2.bulk.galson.p1 <- read_tsv("paper_assets/intermediate_files/Galson_p1_germ-pass.tsv.gz")   # v_identity between 0.60 and 1, sequences renamed, have clone_id & germline_alignment_d_mask, NOT cdr3_aa
cov2.bulk.kc.m5.allreps <- read_tsv("paper_assets/intermediate_files/Kuri-Cervantes_M5-reps1to4_airr-covid-19.tsv.gz")   # v_identity between 0.60 and 1, NOT cdr3_aa

hiv.bulk.nih45 <- read_tsv("paper_assets/intermediate_files/hivnih45_vdjserver.tsv.gz")
hiv.bulk.mt1214 <- read_tsv("paper_assets/intermediate_files/MT1214downloaded.tab.gz")

den.bulk.OAS <- read_tsv("paper_assets/intermediate_files/OAS_sept21_germ-pass2.tsv.gz")
den.bulk.d13enrich <- read_tsv("paper_assets/intermediate_files/d13_2enrichHC_germ-pass.tsv.gz")
den.bulk.d13stim <- read_tsv("paper_assets/intermediate_files/d13_2stimHC_germ-pass.tsv.gz")
hc.BXmay.10mstim <- read_tsv("paper_assets/intermediate_files/BXmay10mstimHC_germ-pass.tsv.gz")

## Setliff 2018 data is too large to be combined, so 12 separate files
hiv.bulk.cap287.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap287_3y.tab.gz")
hiv.bulk.cap287.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap287_6m.tab.gz")
hiv.bulk.cap301.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap301_3y.tab.gz")
hiv.bulk.cap301.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap301_6m.tab.gz")
hiv.bulk.cap312.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap312_3y.tab.gz")
hiv.bulk.cap312.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap312_6m.tab.gz")
hiv.bulk.cap322.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap322_3y.tab.gz")
hiv.bulk.cap322.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap322_6m.tab.gz")
hiv.bulk.cap335.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap335_3y.tab.gz")
hiv.bulk.cap335.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap335_6m.tab.gz")
hiv.bulk.cap351.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap351_3y.tab.gz")
hiv.bulk.cap351.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap351_6m.tab.gz")


### lists of mabs
toshiny.cov2.abdab <- read_tsv("shinyapp/toshiny_cov2_abdab.tab")
den.mabs <- read_tsv("paper_assets/intermediate_files/toshiny_den_mabs0_germ-pass.tsv")

# for HIV 3 sources: IEDB, CATNAP & Yacoob
toshiny.hiv.mabs.all <- read_tsv("paper_assets/intermediate_files/toshiny_hiv_mabs_all.tab")


##################################################################################################################
### COMMANDS FOR READING IN AIRR-COMPLIANT DATASETS FOR SHINY APP
### note that for non-AIRR (e.g. older Immcantation data or databases like AbDab) datasets, column names will be different
######### solution is to use only AIRR datasets
###### particularly note that 1) junction MAY NOT NEED TO BE TRIMMED 1 AA ON EACH END TO BE KABAT CDR3
######                        2) v_identity MAY BE BETWEEN 0-1 NOT 0-100, AND THUS shm FORMULA NEEDS TO CHANGE
### updated shinyprocess function below now addresses both points 1 & 2

### AIRR-COMPLIANT TSV FILES NEED THE FOLLOWING MODIFICATIONS
# CHANGING RELEVANT HEADERS FROM LOWERCASE TO UPPER CASE
# REMOVING FIRST AND LAST RESIDUES FROM junction
# CONVERT v_identity TO shm - FOR SHINY WANT 0-100

# earlier version removing NNNN was changed to only removing X
## create a function to do all pre-processing all at once
## also need to rename sequence_id's??
## required columns: v_identity, v_call, j_call

## FOR SHINY APP PROCESSING THERE CAN BE NO UNDERSCORES IN THE SEQUENCE_ID NAMES - USE DASHES INSTEAD OF UNDERSCORES!!!
## function has this conversion near end

### single function to do all of this pre-processing:
##################################################################################################################

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
  ## removing non-productive, out of frame, stop codons, any X in CDR3 (was 'nnnn' in sequence)
  x <- x %>% filter(productive != "FALSE") %>%
    filter(vj_in_frame != "FALSE") %>%
    filter(productive != "F") %>%
    filter(vj_in_frame != "F")
  x <- x[ grep("\\*", x$junction_aa, invert = TRUE) , ]
  x <- x[ grep("\\X", x$junction_aa, invert = TRUE) , ]
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
    ## adding change from underscores to dashes...
    x$dataset <- deparse(substitute(colname))
    x$dataset <- gsub('"','',x$dataset)
    x$dataset <- gsub("\\.","\\-",x$dataset)
    x$dataset <- gsub("\\_","\\-",x$dataset)
    x$obs <- 1:nrow(x) 
    x <- x %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
    # x <- x %>% relocate(sequence_id, .before = cregion)  ## changed to default i.e. move to make first column
    x <- x %>% relocate(sequence_id)
    # x$sequence_id <- gsub("\\_","\\-",x$sequence_id) ## moved to always run
  }
  ## need to always check and remove underscores from all names
  x$sequence_id <- gsub("\\_","\\-",x$sequence_id)
  return(x)
}

#################################################################
## checking names
den.bulk.OAS$sequence_id[1]          #name ok
den.bulk.d13enrich$sequence_id[1]
den.bulk.d13stim$sequence_id[1]
hiv.bulk.nih45$sequence_id[1]
hiv.bulk.mt1214$sequence_id[1]
cov2.bulk.binder.p11$sequence_id[1]       #name ok but change anyway
cov2.bulk.nielsen.p7450$sequence_id[1]
cov2.bulk.galson.p1$sequence_id[1]    #name ok but change anyway
cov2.bulk.kc.m5.allreps$sequence_id[1]
hc.BXmay.10mstim$sequence_id[1]


## PROCESSING OF THESE INTERMEDIATE FILES USING SHINYPROCESS
toshiny.den.bulk.d13enrich <- shinyprocess(den.bulk.d13enrich)
toshiny.den.bulk.d13stim <- shinyprocess(den.bulk.d13stim)
toshiny.den.bulk.OAS <- shinyprocess(den.bulk.OAS, renumber_sequences = FALSE)

toshiny.hiv.bulk.mt1214 <- shinyprocess(hiv.bulk.mt1214)
toshiny.hiv.bulk.nih45 <- shinyprocess(hiv.bulk.nih45)

toshiny.hc.BXmay.10mstim <- shinyprocess(hc.BXmay.10mstim)

toshiny.cov2.bulk.binder.p11 <- shinyprocess(cov2.bulk.binder.p11)
toshiny.cov2.bulk.galson.p1 <- shinyprocess(cov2.bulk.galson.p1)
toshiny.cov2.bulk.nielsen.p7450 <- shinyprocess(cov2.bulk.nielsen.p7450)
toshiny.cov2.bulk.kc.m5.allreps <- shinyprocess(cov2.bulk.kc.m5.allreps)

## dengue maps need processing, HIV & CoV2 are processed in airrscape_preprocessing script
toshiny.den.mabs <- shinyprocess(den.mabs, renumber_sequences = FALSE, filter_after_counting = FALSE)

## for Setliff datasets, run Shinyprocess on each timepoint individually
toshiny.hiv.bulk.cap287.3y <- shinyprocess(hiv.bulk.cap287.3y)
toshiny.hiv.bulk.cap287.6m <- shinyprocess(hiv.bulk.cap287.6m)
toshiny.hiv.bulk.cap301.3y <- shinyprocess(hiv.bulk.cap301.3y)
toshiny.hiv.bulk.cap301.6m <- shinyprocess(hiv.bulk.cap301.6m)
toshiny.hiv.bulk.cap312.3y <- shinyprocess(hiv.bulk.cap312.3y)
toshiny.hiv.bulk.cap312.6m <- shinyprocess(hiv.bulk.cap312.6m)
toshiny.hiv.bulk.cap322.3y <- shinyprocess(hiv.bulk.cap322.3y)
toshiny.hiv.bulk.cap322.6m <- shinyprocess(hiv.bulk.cap322.6m)
toshiny.hiv.bulk.cap335.3y <- shinyprocess(hiv.bulk.cap335.3y)
toshiny.hiv.bulk.cap335.6m <- shinyprocess(hiv.bulk.cap335.6m)
toshiny.hiv.bulk.cap351.3y <- shinyprocess(hiv.bulk.cap351.3y)
toshiny.hiv.bulk.cap351.6m <- shinyprocess(hiv.bulk.cap351.6m)


## after running function remove large intermediate files
rm(den.bulk.OAS)
rm(den.bulk.d13enrich)
rm(den.bulk.d13stim)
rm(hiv.bulk.nih45)
rm(hiv.bulk.mt1214)
rm(cov2.bulk.binder.p11)
rm(cov2.bulk.nielsen.p7450)
rm(cov2.bulk.galson.p1)
rm(cov2.bulk.kc.m5.allreps)
rm(hc.BXmay.10mstim)
rm(toshiny.hiv.bulk.cap287.3y)
rm(toshiny.hiv.bulk.cap287.6m)
rm(toshiny.hiv.bulk.cap301.3y)
rm(toshiny.hiv.bulk.cap301.6m)
rm(toshiny.hiv.bulk.cap312.3y)
rm(toshiny.hiv.bulk.cap312.6m)
rm(toshiny.hiv.bulk.cap322.3y)
rm(toshiny.hiv.bulk.cap322.6m)
rm(toshiny.hiv.bulk.cap335.3y)
rm(toshiny.hiv.bulk.cap335.6m)
rm(toshiny.hiv.bulk.cap351.3y)
rm(toshiny.hiv.bulk.cap351.6m)

##############################################################################
##############################################################################

### checking datasets...
unique(toshiny.den.bulk.d13stim$cregion)
unique(toshiny.den.bulk.d13enrich$cregion)
unique(toshiny.den.bulk.OAS$cregion)

unique(toshiny.hiv.bulk.nih45$cregion)  ## ONE IGK!!
unique(toshiny.hiv.bulk.mt1214$cregion) ## NA
unique(toshiny.hiv.bulk.cap$cregion)

unique(toshiny.hc.BXmay.10mstim$cregion)
unique(toshiny.cov2.bulk.binder.p11$cregion)
unique(toshiny.cov2.bulk.nielsen.p7450$cregion)
unique(toshiny.cov2.bulk.galson.p1$cregion)
unique(toshiny.cov2.bulk.kc.m5.allreps$cregion)

unique(toshiny.cov2.abdab$cregion)
unique(toshiny.den.mabs$cregion)
unique(toshiny.hiv.mabs.all$cregion)

unique(toshiny.den.bulk.d13stim$gf)
unique(toshiny.den.bulk.d13enrich$gf)
unique(toshiny.den.bulk.OAS$gf)
unique(toshiny.hiv.bulk.nih45$gf)  ## only HV1, HV3, HV4
unique(toshiny.hiv.bulk.mt1214$gf) ## NA
unique(toshiny.hiv.bulk.mt1214$gf)
unique(toshiny.hc.BXmay.10mstim$gf)
unique(toshiny.cov2.bulk.binder.p11$gf)
unique(toshiny.cov2.bulk.nielsen.p7450$gf)
unique(toshiny.cov2.bulk.galson.p1$gf)
unique(toshiny.cov2.bulk.kc.m5.allreps$gf)

unique(toshiny.cov2.abdab$jgene)
unique(toshiny.den.mabs$jgene)

unique(toshiny.den.bulk.d13stim$jgene)
unique(toshiny.den.bulk.d13enrich$jgene)
unique(toshiny.den.bulk.OAS$jgene)
unique(toshiny.hiv.bulk.nih45$jgene)  ## only HV1, HV3, HV4
unique(toshiny.hiv.bulk.mt1214$jgene) ## NA
unique(toshiny.hiv.bulk.mt1214$jgene)
unique(toshiny.hc.BXmay.10mstim$jgene)
unique(toshiny.cov2.bulk.binder.p11$jgene)
unique(toshiny.cov2.bulk.nielsen.p7450$jgene)
unique(toshiny.cov2.bulk.galson.p1$jgene)
unique(toshiny.cov2.bulk.kc.m5.allreps$jgene)
unique(toshiny.cov2.abdab$jgene)
unique(toshiny.hiv.mabs.all$jgene)



## some one off commands...
toshiny.hiv.bulk.mt1214$cregion <- toshiny.hiv.bulk.mt1214$cregion %>% replace_na("IgH")
toshiny.hiv.bulk.nih45 <- toshiny.hiv.bulk.nih45[ grep("IgK", toshiny.hiv.bulk.nih45$cregion, invert = TRUE) , ]

## for mabs make heavy chain subsets
toshiny.den.mabs.h <- subset(toshiny.den.mabs, cregion %in% c("IgH"))
toshiny.cov2.abdab.h <- subset(toshiny.cov2.abdab, cregion %in% c("IgH"))
toshiny.hiv.mabs.all.h <- subset(toshiny.hiv.mabs.all, cregion %in% c("IgH"))

## combining all Setliff data
toshiny.hiv.bulk.cap <- rbind(toshiny.hiv.bulk.cap287.3y,toshiny.hiv.bulk.cap287.6m,toshiny.hiv.bulk.cap301.3y,toshiny.hiv.bulk.cap301.6m,toshiny.hiv.bulk.cap312.3y,toshiny.hiv.bulk.cap312.6m,toshiny.hiv.bulk.cap322.3y,toshiny.hiv.bulk.cap322.6m,toshiny.hiv.bulk.cap335.3y,toshiny.hiv.bulk.cap335.6m,toshiny.hiv.bulk.cap351.3y,toshiny.hiv.bulk.cap351.6m)
## need to recalculate ncount, shm_mean, shm_max
toshiny.hiv.bulk.cap$ncount <- NULL
toshiny.hiv.bulk.cap$shm_mean <- NULL
toshiny.hiv.bulk.cap$shm_max <- NULL
toshiny.hiv.bulk.cap <- toshiny.hiv.bulk.cap %>%
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


## next you may save these individually - further down are combinaation steps
write.table(toshiny.den.bulk.d13stim, "toshiny_den_bulk_d13stim.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.bulk.d13enrich, "toshiny_den_bulk_d13enrich.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.bulk.OAS, "toshiny_den_bulk_OAS.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.bulk.nih45, "toshiny_hiv_bulk_nih45.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.bulk.mt1214, "toshiny_hiv_bulk_mt1214.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.bulk.cap, "toshiny_hiv_bulk_cap.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hc.BXmay.10mstim, "toshiny_hc_BXmay_10mstim.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.binder.p11, "toshiny_cov2_bulk_binder_p11.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.nielsen.p7450, "toshiny_cov2_bulk_nielsen_p7450.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.galson.p1, "toshiny_cov2_bulk_galson_p1.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.bulk.kc.m5.allreps, "toshiny_cov2_bulk_kc.m5_allreps.tab", sep = "\t", row.names = FALSE, quote = FALSE)

### mabs
write.table(toshiny.cov2.abdab, "toshiny_cov2_abdab.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2.abdab.h, "toshiny_cov2_abdab_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.mabs, "toshiny_den_mabs.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.den.mabs.h, "toshiny_den_mabs_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.mabs.all, "toshiny_hiv_mabs_all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.hiv.mabs.all.h, "toshiny_hiv_mabs_all_h.tab", sep = "\t", row.names = FALSE, quote = FALSE)


## combine d13 enrich & d 13 stim, also one-off potential contamination removed...
toshiny.den.bulk.d13 <- full_join(toshiny.den.bulk.d13enrich, toshiny.den.bulk.d13stim)
toshiny.den.bulk.d13 <- toshiny.den.bulk.d13[ grep("ARALFGLVAVASPFDN", toshiny.den.bulk.d13$cdr3_aa_imgt, invert = TRUE) , ]
toshiny.den.bulk.d13 <- toshiny.den.bulk.d13[ grep("ITPPLYLMVGGVSRAMAV", toshiny.den.bulk.d13$cdr3_aa_imgt, invert = TRUE) , ]
toshiny.den.bulk.d13 <- toshiny.den.bulk.d13[ grep("ARQDRNWFDT", toshiny.den.bulk.d13$cdr3_aa_imgt, invert = TRUE) , ]

write.table(toshiny.den.bulk.d13, "toshiny_den_bulk_d13.tab", sep = "\t", row.names = FALSE, quote = FALSE)


#####################################################################################################################
## CODE TO COMBINE DATASETS INTO A SINGLE FILE FOR SHINY VISUALIZATION (BOTH SEPARATE & COMBINED)
#####################################################################################################################

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
toshiny.hiv.mabs.all.h

# all hiv
toshiny.hiv.mabs.all.h
toshiny.hiv.bulk.nih45
toshiny.hiv.bulk.mt1214
toshiny.hiv.bulk.cap

# all den
toshiny.den.mabs
toshiny.den.bulk.OAS
toshiny.den.bulk.d13

## entire collection combined, will be called toshiny.cov2hivden.all
toshiny.cov2.abdab.h
toshiny.cov2.bulk.binder.p11
toshiny.cov2.bulk.galson.p1
toshiny.cov2.bulk.kc.m5.allreps
toshiny.cov2.bulk.nielsen.p7450
toshiny.hc.BXmay.10mstim
toshiny.hiv.mabs.all.h
toshiny.hiv.bulk.nih45
toshiny.hiv.bulk.mt1214
toshiny.hiv.bulk.cap
toshiny.den.mabs
toshiny.den.bulk.OAS
toshiny.den.bulk.d13


#####################################################
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

#####################################################
# toshiny.cov2hiv
#    toshiny.cov2.abdab.h
#     toshiny.hiv.mabs.h

#c("anti-CoV2 mAbs", "anti-HIV mAbs"))
# #"SARS-CoV2 mAbs vs. HIV mAbs - IgH",

toshiny.cov2hiv <- bind_rows(toshiny.cov2.abdab.h, toshiny.hiv.mabs.all.h, .id = "id")

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

##########################################################
# toshiny.hiv.all
#     toshiny.hiv.mabs.h
#     toshiny.hiv.bulk.nih45
#     toshiny.hiv.bulk.mt1214
#     toshiny.hiv.bulk.cap ## added feb 2022

#c("HIV mAbs", "HIV patient MT1214", "HIV patient NIH45"))  #"SARS-CoV2 mAbs vs. HIV mAbs - IgH",

toshiny.hiv.all <- bind_rows(toshiny.hiv.mabs.all.h, toshiny.hiv.bulk.mt1214, toshiny.hiv.bulk.nih45,toshiny.hiv.bulk.cap, .id = "id")

toshiny.hiv.all$id <- gsub("1","HIV mAbs",toshiny.hiv.all$id)
toshiny.hiv.all$id <- gsub("4","HIV Setliff twozerooneeight patient bulk repertoires",toshiny.hiv.all$id)
toshiny.hiv.all$id <- gsub("3","HIV patient NIH45 bulk repertoire",toshiny.hiv.all$id)
toshiny.hiv.all$id <- gsub("2","HIV patient MT1214 bulk repertoire",toshiny.hiv.all$id)
toshiny.hiv.all$id <- gsub("HIV Setliff twozerooneeight patient bulk repertoires","HIV Setliff 2018 patient bulk repertoires",toshiny.hiv.all$id)

## check id columns!!
unique(toshiny.hiv.all$id)
head(toshiny.hiv.all$sequence_id)
toshiny.hiv.all$vj_junction_pattern <- NULL

unique(toshiny.hiv.all$cregion)
#[1] "IgH" "IgM" "IgG" "IgA" NA    "IgK

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


########################################################
# toshiny.den.all
#    toshiny.den.mabs
#    toshiny.den.bulk.d13
#    toshiny.den.bulk.OAS

#c("Dengue plasmablasts", "Dengue patient d13", "Dengue Parameswaran 2013 patients"))  

toshiny.den.all <- bind_rows(toshiny.den.mabs.h, toshiny.den.bulk.d13, toshiny.den.bulk.OAS, .id = "id")

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

#####################################################
# toshiny.cov2hivden.all
toshiny.cov2.abdab.h
toshiny.cov2.bulk.binder.p11
toshiny.cov2.bulk.galson.p1
toshiny.cov2.bulk.kc.m5.allreps
toshiny.cov2.bulk.nielsen.p7450
toshiny.hc.BXmay.10mstim
toshiny.hiv.mabs.h
toshiny.hiv.bulk.nih45
toshiny.hiv.bulk.mt1214
toshiny.hiv.bulk.cap
toshiny.den.mabs
toshiny.den.bulk.OAS
toshiny.den.bulk.d13

toshiny.cov2hivden.all <- bind_rows(toshiny.cov2.abdab.h, toshiny.cov2.bulk.binder.p11, toshiny.cov2.bulk.galson.p1, toshiny.cov2.bulk.kc.m5.allreps, toshiny.cov2.bulk.nielsen.p7450, toshiny.hc.BXmay.10mstim, toshiny.hiv.mabs.all.h, toshiny.hiv.bulk.mt1214, toshiny.hiv.bulk.nih45, toshiny.hiv.bulk.cap, toshiny.den.mabs.h, toshiny.den.bulk.d13, toshiny.den.bulk.OAS, .id = "id")

toshiny.cov2hivden.all$id <- gsub("10","HIV patient MTonetwoonefour bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("11","Dengue plasmablasts",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("12","Dengue patient dthirteen",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("13","Dengue Parameswaran twozeroonethree patient bulk repertoires",toshiny.cov2hivden.all$id)

toshiny.cov2hivden.all$id <- gsub("1","SARS-CoVtwo mAbs",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("2","COVID-nineteen patient Binder peleven",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("3","COVID-nineteen patient Galson pone",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("4","COVID-nineteen patient Kuri-Cervantes mfive",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("5","COVID-nineteen patient Nielsen psevenfourfivezero",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("6","Healthy control bulk repertoire",toshiny.cov2hivden.all$id)

toshiny.cov2hivden.all$id <- gsub("7","HIV mAbs",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("8","HIV Setliff twozerooneeight patient bulk repertoires",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("9","HIV patient NIHfourfive bulk repertoire",toshiny.cov2hivden.all$id)


toshiny.cov2hivden.all$id <- gsub("SARS-CoVtwo mAbs","SARS-CoV2 mAbs",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("COVID-nineteen patient Binder peleven","COVID-19 patient Binder p11 bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("COVID-nineteen patient Galson pone","COVID-19 patient Galson p1 bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("COVID-nineteen patient Kuri-Cervantes mfive","COVID-19 patient Kuri-Cervantes m5 bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("COVID-nineteen patient Nielsen psevenfourfivezero","COVID-19 patient Nielsen p7450 bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("HIV Setliff twozerooneeight patient bulk repertoires","HIV Setliff 2018 patient bulk repertoires",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("Dengue patient dthirteen","Dengue patient d13 bulk repertoire",toshiny.cov2hivden.all$id)

toshiny.cov2hivden.all$id <- gsub("HIV patient NIHfourfive bulk repertoire","HIV patient NIH45 bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("HIV patient MTonetwoonefour bulk repertoire","HIV patient MT1214 bulk repertoire",toshiny.cov2hivden.all$id)
toshiny.cov2hivden.all$id <- gsub("Dengue Parameswaran twozeroonethree patient bulk repertoires","Dengue Parameswaran 2013 patient bulk repertoires",toshiny.cov2hivden.all$id)


unique(toshiny.cov2hivden.all$id)
head(toshiny.cov2hivden.all$sequence_id)


## for all combined need to recalculate ncount and shm_mean
toshiny.cov2hivden.allc <- toshiny.cov2hivden.all
toshiny.cov2hivden.allc$ncount <- NULL
toshiny.cov2hivden.allc$shm_mean <- NULL
toshiny.cov2hivden.allc$shm_max <- NULL
toshiny.cov2hivden.allc <- toshiny.cov2hivden.allc %>%
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
toshiny.cov2hivden.allc$cregion0 <- toshiny.cov2hivden.allc$cregion
toshiny.cov2hivden.allc$cregion <- "IgH"

write.table(toshiny.cov2hivden.all, "toshiny_cov2hivden_all.tab", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(toshiny.cov2hivden.allc, "toshiny_cov2hivden_allc.tab", sep = "\t", row.names = FALSE, quote = FALSE)


##############################################################################################################################
##############################################################################################################################

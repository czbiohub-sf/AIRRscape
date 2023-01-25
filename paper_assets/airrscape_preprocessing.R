# install.packages("tidyverse")
# install.packages("knitr")
library(knitr)
library(tidyverse)

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING CoV2-Nielsen7450 data from iReceptor

## login to iReceptor, search metadata for 7450 (& Target Substrate = DNA)
## aka https://gateway.ireceptor.org/sequences?query_id=59121
## results in 289023 sequences: saved as Nielsen_7450_airr-covid-19.tsv

## this file is already available in the /paper_assets/intermediate_files folder
# cov2.bulk.nielsen.p7450 <- read_tsv("paper_assets/intermediate_files/Nielsen_7450_airr-covid-19.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING CoV2-Kuri-Cervantes_M5 data from iReceptor

## login to iReceptor, search metadata for SubjectID = M5 & Diagnosis = COVID-19
## aka https://gateway.ireceptor.org/samples?query_id=59127 & https://gateway.ireceptor.org/sequences-download?query_id=59128&n=1136400&page=sequences
## results in 1,136,400 sequences: saved as Kuri-Cervantes_M5-reps1to4_airr-covid-19.tsv

## this file is already available in the /paper_assets/intermediate_files folder
# cov2.bulk.kc.m5.allreps <- read_tsv("paper_assets/intermediate_files/Kuri-Cervantes_M5-reps1to4_airr-covid-19.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING CoV2-Schultheiß-p11 data from iReceptor

## login to iReceptor, search metadata for SubjectID = Pt-11 & Diagnosis = COVID-19 & PCR Target = IGH
## aka https://gateway.ireceptor.org/samples?query_id=59148
## results in 47,936 sequences: Binder-pt11-airr-covid-19.tsv

## BUT THESE ARE MISSING SHM DATA, SO RUNNING THROUGH IMMCANTATION
## converting to fasta, running these commands from the Immcantation docker container
## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
     ## we used (docker pull immcantation/suite:4.1.0)
      # changeo-igblast -s /data/Binder-pt11-airr-covid-19.fasta -n Binder_p11 -o /data
      # changeo-clone -d /data/Binder_p11_db-pass.tsv -n Binder_p11 -o /data -x 0.12
# resulting file is Binder_p11_germ-pass.tsv - 33,852 sequences

## this file is already available in the /paper_assets/intermediate_files folder
# cov2.bulk.binder.p11 <- read_tsv("paper_assets/intermediate_files/Binder_p11_germ-pass.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING CoV2-Galson-p1 data from iReceptor

## login to iReceptor, search metadata for Covid_Barts_01
## aka https://gateway.ireceptor.org/samples?query_id=59155
## results in 64,415 sequences: Galson_Covid_Barts_01_airr-covid-19.tsv

## BUT THESE ARE MISSING SHM DATA, SO RUNNING THROUGH IMMCANTATION
## converting to fasta, running these commands from the Immcantation docker container
## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
     ## we used (docker pull immcantation/suite:4.1.0)
      # changeo-igblast -s /data/galson_01a.fasta -n Galson_p1 -o /data
      # changeo-clone -d /data/Galson_p1_db-pass.tsv -n Galson_p1 -o /data -x 0.12
# resulting file is Galson_p1_germ-pass.tsv - 63,953 sequences

## this file is already available in the /paper_assets/intermediate_files folder
# cov2.bulk.galson.p1 <- read_tsv("paper_assets/intermediate_files/Galson_p1_germ-pass.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING healthy control data from Waltari 2019

## these are only in SRA, so need full processing through Immcantation
## d536M_r1_stim
## https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8647465
# fastq-dump --split-files --orivgfmt --gzip SRR8647465

## first step in converting raw 300x250 MiSeq data uses Immcantation:Presto
## see Waltari 2019 for details https://doi.org/10.3389/fimmu.2019.01452
## result is a fastq file (presto-abseq.sh file changed to require 3 or more consensus reads i.e. MIN_CONSCOUNT=3)
# BXmaykm_10m_stimulated_S2-final_collapse-unique_atleast-3.fastq

## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
## we used (docker pull immcantation/suite:4.1.0)
# changeo-igblast -s /data/BXmaykm_10m_stimulated_S2-final_collapse-unique_atleast-3.fastq -n BXmay10mstim3 -o /data
## FOR THIS FILE - REMOVED LC READS (grep -E -v "Kappa|Lambda" BXmay10mstim3_db-pass.tsv > BXmay10mstim3_db-passHC.tsv)
# changeo-clone -d /data/BXmay10mstim3_db-passHC.tsv -n BXmay10mstimHC -o /data -x 0.12
# resulting file is BXmay10mstimHC_germ-pass.tsv - 77,672 sequences

## this file is already available in the /paper_assets/intermediate_files folder
# hc.BXmay.10mstim <- read_tsv("paper_assets/intermediate_files/BXmay10mstimHC_germ-pass.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING dengue-Durham2019 data

## these are only in SRA, so need full processing through Immcantation
## d13-enriched
## https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR17417727
## d13-stimulated
## https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR17417728
# fastq-dump --split-files --orivgfmt --gzip SRR17417727
# fastq-dump --split-files --orivgfmt --gzip SRR17417728

## first step in converting raw 300x250 MiSeq data uses Immcantation:Presto
## see Waltari 2019 for details https://doi.org/10.3389/fimmu.2019.01452
## result are 2 fastq files
# 20190112_d13_2_enriched-final_collapse-unique_atleast-2.fastq
# 20190112_d13_2_stimulated-final_collapse-unique_atleast-2.fastq

## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
## we used (docker pull immcantation/suite:4.1.0)
# changeo-igblast -s /data/20190112_d13_2_enriched-final_collapse-unique_atleast-2.fastq -n d13_2enrich -o /data
## from resulting db-pass tsv file removed all LC reads
# changeo-clone -d /data/d13_2enrich_db-passHC.tsv -n d13_2enrichHC -o /data -x 0.12
# changeo-igblast -s /data/20190113_d13_2_stimulated-final_collapse-unique_atleast-2.fastq -n d13_2stim -o /data
## from resulting db-pass tsv file removed all LC reads
# changeo-clone -d /data/d13_2stim_db-passHC.tsv -n d13_2stimHC -o /data -x 0.12

# resulting files are d13_2enrichHC_germ-pass.tsv - 19,718 sequences & d13_2stimHC_germ-pass.tsv - 146,529 sequences

## these files are already available in the /paper_assets/intermediate_files folder
# den.bulk.d13enrich <- read_tsv("paper_assets/intermediate_files/d13_2enrichHC_germ-pass.tsv.gz")
# den.bulk.d13stim <- read_tsv("paper_assets/intermediate_files/d13_2stimHC_germ-pass.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING Parameswaran_2013 DENGUE DATA FROM OAS
## these can be found by searching OAS for unpaired dengue sequences:
## result is 127 separate files all from Parameswaran 2013
## wget commands and key combining filenames with datasets found in:
# /AIRRscape/paper_assets/intermediate_files/oas/bulk_download_dengue_paramaswaran.xlsx
## use those wget commands, or just use the files already in the /paper_assets/intermediate_files/oas/ folder

## then run these R commands
OAS.clusters.p148a <- read_csv("paper_assets/intermediate_files/oas/SRR2150126_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p148b <- read_csv("paper_assets/intermediate_files/oas/SRR2150229_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p148c <- read_csv("paper_assets/intermediate_files/oas/SRR2150329_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p172a <- read_csv("paper_assets/intermediate_files/oas/SRR2150457_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p172b <- read_csv("paper_assets/intermediate_files/oas/SRR2150481_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p172c <- read_csv("paper_assets/intermediate_files/oas/SRR2150504_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p194a <- read_csv("paper_assets/intermediate_files/oas/SRR2150549_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p194b <- read_csv("paper_assets/intermediate_files/oas/SRR2150573_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p194c <- read_csv("paper_assets/intermediate_files/oas/SRR2150597_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p199a <- read_csv("paper_assets/intermediate_files/oas/SRR2150643_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p199b <- read_csv("paper_assets/intermediate_files/oas/SRR2150668_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p199c <- read_csv("paper_assets/intermediate_files/oas/SRR2150692_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p203a <- read_csv("paper_assets/intermediate_files/oas/SRR2150715_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p203b <- read_csv("paper_assets/intermediate_files/oas/SRR2150734_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p203c <- read_csv("paper_assets/intermediate_files/oas/SRR2150753_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p208a <- read_csv("paper_assets/intermediate_files/oas/SRR2150802_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p208b <- read_csv("paper_assets/intermediate_files/oas/SRR2150816_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p208c <- read_csv("paper_assets/intermediate_files/oas/SRR2150838_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p232a <- read_csv("paper_assets/intermediate_files/oas/SRR2150935_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p232b <- read_csv("paper_assets/intermediate_files/oas/SRR2150947_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p232c <- read_csv("paper_assets/intermediate_files/oas/SRR2150972_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p237a <- read_csv("paper_assets/intermediate_files/oas/SRR2151066_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p237b <- read_csv("paper_assets/intermediate_files/oas/SRR2151089_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p237c <- read_csv("paper_assets/intermediate_files/oas/SRR2151105_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p238a <- read_csv("paper_assets/intermediate_files/oas/SRR2151162_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p238b <- read_csv("paper_assets/intermediate_files/oas/SRR2151187_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p238c <- read_csv("paper_assets/intermediate_files/oas/SRR2151211_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p240a <- read_csv("paper_assets/intermediate_files/oas/SRR2151231_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p240b <- read_csv("paper_assets/intermediate_files/oas/SRR2151247_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p240c <- read_csv("paper_assets/intermediate_files/oas/SRR2151270_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p249a <- read_csv("paper_assets/intermediate_files/oas/SRR2151293_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p249b <- read_csv("paper_assets/intermediate_files/oas/SRR2151316_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p249c <- read_csv("paper_assets/intermediate_files/oas/SRR2151330_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p252a <- read_csv("paper_assets/intermediate_files/oas/SRR2151353_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p252b <- read_csv("paper_assets/intermediate_files/oas/SRR2151376_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p252c <- read_csv("paper_assets/intermediate_files/oas/SRR2151395_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p255a <- read_csv("paper_assets/intermediate_files/oas/SRR2151414_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p255b <- read_csv("paper_assets/intermediate_files/oas/SRR2151435_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p255c <- read_csv("paper_assets/intermediate_files/oas/SRR2151450_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p265a <- read_csv("paper_assets/intermediate_files/oas/SRR2151499_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p265b <- read_csv("paper_assets/intermediate_files/oas/SRR2151523_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p275a <- read_csv("paper_assets/intermediate_files/oas/SRR2151538_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p275b <- read_csv("paper_assets/intermediate_files/oas/SRR2151562_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p275c <- read_csv("paper_assets/intermediate_files/oas/SRR2151598_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p276a <- read_csv("paper_assets/intermediate_files/oas/SRR2151713_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p276b <- read_csv("paper_assets/intermediate_files/oas/SRR2151738_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p276c <- read_csv("paper_assets/intermediate_files/oas/SRR2151761_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p287a <- read_csv("paper_assets/intermediate_files/oas/SRR2153023_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p287b <- read_csv("paper_assets/intermediate_files/oas/SRR2153024_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p287c <- read_csv("paper_assets/intermediate_files/oas/SRR2153025_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p289a <- read_csv("paper_assets/intermediate_files/oas/SRR2153026_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p289b <- read_csv("paper_assets/intermediate_files/oas/SRR2153027_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p289c <- read_csv("paper_assets/intermediate_files/oas/SRR2153028_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p299a <- read_csv("paper_assets/intermediate_files/oas/SRR2153029_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p299b <- read_csv("paper_assets/intermediate_files/oas/SRR2153030_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p299c <- read_csv("paper_assets/intermediate_files/oas/SRR2153031_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p301a <- read_csv("paper_assets/intermediate_files/oas/SRR2153032_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p301b <- read_csv("paper_assets/intermediate_files/oas/SRR2153033_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p301c <- read_csv("paper_assets/intermediate_files/oas/SRR2153034_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p307a <- read_csv("paper_assets/intermediate_files/oas/SRR2153035_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p307b <- read_csv("paper_assets/intermediate_files/oas/SRR2153036_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p311a <- read_csv("paper_assets/intermediate_files/oas/SRR2153037_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p311b <- read_csv("paper_assets/intermediate_files/oas/SRR2153038_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p320a <- read_csv("paper_assets/intermediate_files/oas/SRR2153039_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p320b <- read_csv("paper_assets/intermediate_files/oas/SRR2153040_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p346a <- read_csv("paper_assets/intermediate_files/oas/SRR2153045_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p376a <- read_csv("paper_assets/intermediate_files/oas/SRR2153046_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p376b <- read_csv("paper_assets/intermediate_files/oas/SRR2153047_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p376c <- read_csv("paper_assets/intermediate_files/oas/SRR2153048_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p391a <- read_csv("paper_assets/intermediate_files/oas/SRR2153049_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p391b <- read_csv("paper_assets/intermediate_files/oas/SRR2153050_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p422a <- read_csv("paper_assets/intermediate_files/oas/SRR2153052_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p422b <- read_csv("paper_assets/intermediate_files/oas/SRR2153053_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p422c <- read_csv("paper_assets/intermediate_files/oas/SRR2153054_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p444a <- read_csv("paper_assets/intermediate_files/oas/SRR2153056_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p444b <- read_csv("paper_assets/intermediate_files/oas/SRR2153057_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p444c <- read_csv("paper_assets/intermediate_files/oas/SRR2153058_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p455a <- read_csv("paper_assets/intermediate_files/oas/SRR2153060_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p455b <- read_csv("paper_assets/intermediate_files/oas/SRR2153061_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p455c <- read_csv("paper_assets/intermediate_files/oas/SRR2153062_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p479a <- read_csv("paper_assets/intermediate_files/oas/SRR2153063_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p479b <- read_csv("paper_assets/intermediate_files/oas/SRR2153064_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p479c <- read_csv("paper_assets/intermediate_files/oas/SRR2153065_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p481a <- read_csv("paper_assets/intermediate_files/oas/SRR2153066_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p481b <- read_csv("paper_assets/intermediate_files/oas/SRR2153067_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p481c <- read_csv("paper_assets/intermediate_files/oas/SRR2153068_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p489a <- read_csv("paper_assets/intermediate_files/oas/SRR2153070_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p489b <- read_csv("paper_assets/intermediate_files/oas/SRR2153071_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p489c <- read_csv("paper_assets/intermediate_files/oas/SRR2153072_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p500a <- read_csv("paper_assets/intermediate_files/oas/SRR2153073_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p500b <- read_csv("paper_assets/intermediate_files/oas/SRR2153074_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p500c <- read_csv("paper_assets/intermediate_files/oas/SRR2153075_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p514a <- read_csv("paper_assets/intermediate_files/oas/SRR2153230_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p514b <- read_csv("paper_assets/intermediate_files/oas/SRR2153231_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p514c <- read_csv("paper_assets/intermediate_files/oas/SRR2153232_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p515a <- read_csv("paper_assets/intermediate_files/oas/SRR2153233_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p515b <- read_csv("paper_assets/intermediate_files/oas/SRR2153234_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p515c <- read_csv("paper_assets/intermediate_files/oas/SRR2153235_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p517a <- read_csv("paper_assets/intermediate_files/oas/SRR2153236_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p517b <- read_csv("paper_assets/intermediate_files/oas/SRR2153237_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p517c <- read_csv("paper_assets/intermediate_files/oas/SRR2153238_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p520a <- read_csv("paper_assets/intermediate_files/oas/SRR2153239_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p520b <- read_csv("paper_assets/intermediate_files/oas/SRR2153240_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p520c <- read_csv("paper_assets/intermediate_files/oas/SRR2153241_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p524a <- read_csv("paper_assets/intermediate_files/oas/SRR2153242_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p524b <- read_csv("paper_assets/intermediate_files/oas/SRR2153243_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p524c <- read_csv("paper_assets/intermediate_files/oas/SRR2153244_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p529a <- read_csv("paper_assets/intermediate_files/oas/SRR2153245_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p529b <- read_csv("paper_assets/intermediate_files/oas/SRR2153247_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p529c <- read_csv("paper_assets/intermediate_files/oas/SRR2153248_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p543a <- read_csv("paper_assets/intermediate_files/oas/SRR2153249_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p543b <- read_csv("paper_assets/intermediate_files/oas/SRR2153250_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p543c <- read_csv("paper_assets/intermediate_files/oas/SRR2153251_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p551a <- read_csv("paper_assets/intermediate_files/oas/SRR2153252_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p551b <- read_csv("paper_assets/intermediate_files/oas/SRR2153253_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p551c <- read_csv("paper_assets/intermediate_files/oas/SRR2153254_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p555a <- read_csv("paper_assets/intermediate_files/oas/SRR2153255_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p555b <- read_csv("paper_assets/intermediate_files/oas/SRR2153256_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p555c <- read_csv("paper_assets/intermediate_files/oas/SRR2153258_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p558a <- read_csv("paper_assets/intermediate_files/oas/SRR2153261_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p558b <- read_csv("paper_assets/intermediate_files/oas/SRR2153262_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p563a <- read_csv("paper_assets/intermediate_files/oas/SRR2153263_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p563b <- read_csv("paper_assets/intermediate_files/oas/SRR2153264_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p563c <- read_csv("paper_assets/intermediate_files/oas/SRR2153265_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p569a <- read_csv("paper_assets/intermediate_files/oas/SRR2153266_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p569b <- read_csv("paper_assets/intermediate_files/oas/SRR2153267_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.p569c <- read_csv("paper_assets/intermediate_files/oas/SRR2153268_Heavy_Bulk.csv.gz", skip =1)

OAS.clusters.p346 <- read_csv("paper_assets/intermediate_files/oas/SRR2153045_Heavy_Bulk.csv.gz", skip =1)
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

den.bulk.OAS <- rbind(OAS.clusters.p148,OAS.clusters.p172,OAS.clusters.p194,OAS.clusters.p199,OAS.clusters.p203,OAS.clusters.p208,OAS.clusters.p232,OAS.clusters.p237,OAS.clusters.p238,OAS.clusters.p240,OAS.clusters.p249,OAS.clusters.p252,OAS.clusters.p255,OAS.clusters.p265,OAS.clusters.p275,OAS.clusters.p276,OAS.clusters.p287,OAS.clusters.p289,OAS.clusters.p299,OAS.clusters.p301,OAS.clusters.p307,OAS.clusters.p311,OAS.clusters.p320,OAS.clusters.p346,OAS.clusters.p376,OAS.clusters.p391,OAS.clusters.p422,OAS.clusters.p444,OAS.clusters.p455,OAS.clusters.p479,OAS.clusters.p481,OAS.clusters.p489,OAS.clusters.p500,OAS.clusters.p514,OAS.clusters.p515,OAS.clusters.p517,OAS.clusters.p520,OAS.clusters.p524,OAS.clusters.p529,OAS.clusters.p543,OAS.clusters.p551,OAS.clusters.p555,OAS.clusters.p558,OAS.clusters.p563,OAS.clusters.p569)

den.bulk.OAS$ANARCI_numbering <- NULL
den.bulk.OAS$ANARCI_status <- NULL
# write.table(den.bulk.OAS, "OAS_sept21_germ-pass2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# rm(list=ls(pattern="OAS."))

## this file "OAS_sept21_germ-pass2.tsv" is already available in the /paper_assets/intermediate_files folder
# den.bulk.OAS <- read_tsv("paper_assets/intermediate_files/OAS_sept21_germ-pass2.tsv.gz")


##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING HIV-NIH45 from iReceptor

## login to iReceptor, search metadata for SRS823746
## aka https://gateway.ireceptor.org/samples?query_id=59112
## results in 186612 sequences: saved as hivnih45_vdjserver.tsv
## this file is already available in the /paper_assets/intermediate_files folder
# hiv.bulk.nih45 <- read_tsv("paper_assets/intermediate_files/hivnih45_vdjserver.tsv.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING HIV-MT1214 DATA from Waltari 2018 (on OAS)

### mt1214 from OAS download:
## can download from OAS, then load removing first line
mt1214e_summaryurl <- paste0("http://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Waltari_2018/csv/SRR5811762_Heavy_IGHE.csv.gz" )
mt1214e_summarydest <- paste0("SRR5811762_Heavy_IGHE.csv.gz" )
download.file(mt1214e_summaryurl, mt1214e_summarydest)
OAS.clusters.mt1214d <- read_csv("SRR5811762_Heavy_IGHE.csv.gz", skip =1)
OAS.clusters.mt1214e <- OAS.clusters.mt1214e %>% filter(Redundancy > 1)  ## none left!
unlink("SRR5811762_Heavy_IGHE.csv.gz")
mt1214d_summaryurl <- paste0("http://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Waltari_2018/csv/SRR5811762_Heavy_IGHD.csv.gz" )
mt1214d_summarydest <- paste0("SRR5811762_Heavy_IGHD.csv.gz" )
download.file(mt1214d_summaryurl, mt1214d_summarydest)
OAS.clusters.mt1214d <- read_csv("SRR5811762_Heavy_IGHD.csv.gz", skip =1)
OAS.clusters.mt1214d <- OAS.clusters.mt1214d %>% filter(Redundancy > 1)  ## none left!
unlink("SRR5811762_Heavy_IGHD.csv.gz")

mt1214a_summaryurl <- paste0("http://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Waltari_2018/csv/SRR5811762_Heavy_IGHA.csv.gz" )
mt1214a_summarydest <- paste0("SRR5811762_Heavy_IGHA.csv.gz" )
download.file(mt1214a_summaryurl, mt1214a_summarydest)
OAS.clusters.mt1214a <- read_csv("SRR5811762_Heavy_IGHA.csv.gz", skip =1)
OAS.clusters.mt1214a <- OAS.clusters.mt1214a %>% filter(Redundancy > 1)
unlink("SRR5811762_Heavy_IGHA.csv.gz")

mt1214g_summaryurl <- paste0("http://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Waltari_2018/csv/SRR5811762_Heavy_IGHG.csv.gz" )
mt1214g_summarydest <- paste0("SRR5811762_Heavy_IGHG.csv.gz" )
download.file(mt1214g_summaryurl, mt1214g_summarydest)
OAS.clusters.mt1214g <- read_csv("SRR5811762_Heavy_IGHG.csv.gz", skip =1)
OAS.clusters.mt1214g <- OAS.clusters.mt1214g %>% filter(Redundancy > 1)
unlink("SRR5811762_Heavy_IGHG.csv.gz")

mt1214m_summaryurl <- paste0("http://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Waltari_2018/csv/SRR5811762_Heavy_IGHM.csv.gz" )
mt1214m_summarydest <- paste0("SRR5811762_Heavy_IGHM.csv.gz" )
download.file(mt1214m_summaryurl, mt1214m_summarydest)
OAS.clusters.mt1214m <- read_csv("SRR5811762_Heavy_IGHM.csv.gz", skip =1)
OAS.clusters.mt1214m <- OAS.clusters.mt1214m %>% filter(Redundancy > 1)
unlink("SRR5811762_Heavy_IGHM.csv.gz")

mt1214b_summaryurl <- paste0("http://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Waltari_2018/csv/SRR5811762_Heavy_Bulk.csv.gz" )
mt1214b_summarydest <- paste0("SRR5811762_Heavy_Bulk.csv.gz" )
download.file(mt1214b_summaryurl, mt1214b_summarydest)
OAS.clusters.mt1214b <- read_csv("SRR5811762_Heavy_Bulk.csv.gz", skip =1)
OAS.clusters.mt1214b <- OAS.clusters.mt1214b %>% filter(Redundancy > 1)
unlink("SRR5811762_Heavy_Bulk.csv.gz")


OAS.clusters.mt1214a$cregion <- "IgA"
OAS.clusters.mt1214g$cregion <- "IgG"
OAS.clusters.mt1214m$cregion <- "IgM"

## combine, then locus IGH
OAS.clusters.mt1214agm <- rbind(OAS.clusters.mt1214a, OAS.clusters.mt1214g, OAS.clusters.mt1214m)
OAS.clusters.mt1214all <- full_join(OAS.clusters.mt1214agm, OAS.clusters.mt1214b)
OAS.clusters.mt1214all$locus <- "IGH"
OAS.clusters.mt1214all$ANARCI_numbering <- NULL
OAS.clusters.mt1214all$ANARCI_status <- NULL

## RUN THESE ON ALL DATASETS TO BE USED - ALSO LOOK INTO RENUMBERING WITHIN R
OAS.clusters.mt1214all$junction_length_check <- OAS.clusters.mt1214all$junction_length / 3
OAS.clusters.mt1214all <- OAS.clusters.mt1214all %>% filter(is.wholenumber(junction_length_check))

OAS.clusters.mt1214all$dataset <- "mt1214"
OAS.clusters.mt1214all$obs <- 1:nrow(OAS.clusters.mt1214all) 
hiv.bulk.mt1214 <- OAS.clusters.mt1214all %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)

#write.table(hiv.bulk.mt1214, "MT1214downloaded.tab", sep = "\t", row.names = FALSE, quote = FALSE)
# rm(list=ls(pattern="OAS."))

## this file "MT1214downloaded.tab" is already available in the /paper_assets/intermediate_files folder
# hiv.bulk.mt1214 <- read_tsv("paper_assets/intermediate_files/MT1214downloaded.tab.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING SETLIFF 2018 HIV DATA FROM OAS AND COMBINING 100 RAW DATASETS INTO 12 DONOR TIMEPOINTS:
## wget commands and key combining filenames with datasets found in:
# /AIRRscape/paper_assets/intermediate_files/bulk_download_hiv_setliff.xlsx
 ## use those wget commands, then run these R commands
hiv.bulk.cap287.3y.b1 <- read_csv("SRR6206380_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap287.3y.b1 <- hiv.bulk.cap287.3y.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap287.3y.b2 <- read_csv("SRR6207009_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap287.3y.b2 <- hiv.bulk.cap287.3y.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap287.3y.A1 <- read_csv("SRR6207009_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap287.3y.A2 <- read_csv("SRR6206380_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap287.3y.D1 <- read_csv("SRR6206380_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap287.3y.D2 <- read_csv("SRR6207009_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap287.3y.G1 <- read_csv("SRR6206380_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap287.3y.M1 <- read_csv("SRR6206380_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap287.6m.b1 <- read_csv("SRR6206379_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap287.6m.b1 <- hiv.bulk.cap287.6m.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap287.6m.b2 <- read_csv("SRR6207016_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap287.6m.b2 <- hiv.bulk.cap287.6m.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap287.6m.A1 <- read_csv("SRR6207016_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap287.6m.A2 <- read_csv("SRR6206379_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap287.6m.D1 <- read_csv("SRR6207016_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap287.6m.D2 <- read_csv("SRR6206379_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap287.6m.G1 <- read_csv("SRR6206379_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap287.6m.M1 <- read_csv("SRR6206379_Heavy_IGHM.csv.gz", skip =1)

hiv.bulk.cap287.3y <- rbind(hiv.bulk.cap287.3y.b1,hiv.bulk.cap287.3y.b2,hiv.bulk.cap287.3y.A1,hiv.bulk.cap287.3y.A2,hiv.bulk.cap287.3y.D1,hiv.bulk.cap287.3y.D2,hiv.bulk.cap287.3y.G1,hiv.bulk.cap287.3y.M1)
hiv.bulk.cap287.3y <- hiv.bulk.cap287.3y  %>% filter(Redundancy > 1)
hiv.bulk.cap287.6m <- rbind(hiv.bulk.cap287.6m.b1,hiv.bulk.cap287.6m.b2,hiv.bulk.cap287.6m.A1,hiv.bulk.cap287.6m.A2,hiv.bulk.cap287.6m.D1,hiv.bulk.cap287.6m.D2,hiv.bulk.cap287.6m.G1,hiv.bulk.cap287.6m.M1)
hiv.bulk.cap287.6m <- hiv.bulk.cap287.6m %>% filter(Redundancy > 1)
#write.table(hiv.bulk.cap287.3y, "hiv_bulk_cap287_3y.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(hiv.bulk.cap287.6m, "hiv_bulk_cap287_6m.tab", sep = "\t", row.names = FALSE, quote = FALSE)


hiv.bulk.cap301.3y.b1 <- read_csv("SRR6207006_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap301.3y.b1 <- hiv.bulk.cap301.3y.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap301.3y.b2 <- read_csv("SRR6206384_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap301.3y.b2 <- hiv.bulk.cap301.3y.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap301.3y.A1 <- read_csv("SRR6206384_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap301.3y.A2 <- read_csv("SRR6207006_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap301.3y.D1 <- read_csv("SRR6207006_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap301.3y.D2 <- read_csv("SRR6206384_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap301.3y.E1 <- read_csv("SRR6207006_Heavy_IGHE.csv.gz", skip =1)
hiv.bulk.cap301.3y.G1 <- read_csv("SRR6206384_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap301.3y.G2 <- read_csv("SRR6207006_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap301.3y.M1 <- read_csv("SRR6207006_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap301.3y.M2 <- read_csv("SRR6206384_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap301.6m.b1 <- read_csv("SRR6207012_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap301.6m.b1 <- hiv.bulk.cap301.6m.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap301.6m.b2 <- read_csv("SRR6206383_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap301.6m.b2 <- hiv.bulk.cap301.6m.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap301.6m.A1 <- read_csv("SRR6206383_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap301.6m.A2 <- read_csv("SRR6207012_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap301.6m.D1 <- read_csv("SRR6206383_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap301.6m.G1 <- read_csv("SRR6207012_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap301.6m.G2 <- read_csv("SRR6206383_Heavy_IGHG.csv.gz", skip =1)

hiv.bulk.cap301.3y <- rbind(hiv.bulk.cap301.3y.b1,hiv.bulk.cap301.3y.b2,hiv.bulk.cap301.3y.A1,hiv.bulk.cap301.3y.A2,hiv.bulk.cap301.3y.D1,hiv.bulk.cap301.3y.D2,hiv.bulk.cap301.3y.E1,hiv.bulk.cap301.3y.G1,hiv.bulk.cap301.3y.G2,hiv.bulk.cap301.3y.M1,hiv.bulk.cap301.3y.M2)
hiv.bulk.cap301.3y <- hiv.bulk.cap301.3y  %>% filter(Redundancy > 1)
hiv.bulk.cap301.6m <- rbind(hiv.bulk.cap301.6m.b1,hiv.bulk.cap301.6m.b2,hiv.bulk.cap301.6m.A1,hiv.bulk.cap301.6m.A2,hiv.bulk.cap301.6m.D1,hiv.bulk.cap301.6m.G1,hiv.bulk.cap301.6m.G2)
hiv.bulk.cap301.6m <- hiv.bulk.cap301.6m  %>% filter(Redundancy > 1)
#write.table(hiv.bulk.cap301.3y, "hiv_bulk_cap301_3y.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(hiv.bulk.cap301.6m, "hiv_bulk_cap301_6m.tab", sep = "\t", row.names = FALSE, quote = FALSE)


hiv.bulk.cap312.3y.b1 <- read_csv("SRR6206375_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap312.3y.b1 <- hiv.bulk.cap312.3y.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap312.3y.b2 <- read_csv("SRR6207010_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap312.3y.b2 <- hiv.bulk.cap312.3y.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap312.3y.A1 <- read_csv("SRR6207010_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap312.3y.A2 <- read_csv("SRR6206375_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap312.3y.D1 <- read_csv("SRR6206375_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap312.3y.D2 <- read_csv("SRR6207010_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap312.3y.G1 <- read_csv("SRR6206375_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap312.3y.G2 <- read_csv("SRR6207010_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap312.6m.b1 <- read_csv("SRR6206369_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap312.6m.b1 <- hiv.bulk.cap312.6m.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap312.6m.b2 <- read_csv("SRR6207014_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap312.6m.b2 <- hiv.bulk.cap312.6m.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap312.6m.A1 <- read_csv("SRR6206369_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap312.6m.A2 <- read_csv("SRR6207014_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap312.6m.D1 <- read_csv("SRR6207014_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap312.6m.D2 <- read_csv("SRR6206369_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap312.6m.G1 <- read_csv("SRR6206369_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap312.6m.G2 <- read_csv("SRR6207014_Heavy_IGHG.csv.gz", skip =1)

hiv.bulk.cap312.3y <- rbind(hiv.bulk.cap312.3y.b1,hiv.bulk.cap312.3y.b2,hiv.bulk.cap312.3y.A1,hiv.bulk.cap312.3y.A2,hiv.bulk.cap312.3y.D1,hiv.bulk.cap312.3y.D2,hiv.bulk.cap312.3y.G1,hiv.bulk.cap312.3y.G2)
hiv.bulk.cap312.3y <- hiv.bulk.cap312.3y  %>% filter(Redundancy > 1)
hiv.bulk.cap312.6m <- rbind(hiv.bulk.cap312.6m.b1,hiv.bulk.cap312.6m.b2,hiv.bulk.cap312.6m.A1,hiv.bulk.cap312.6m.A2,hiv.bulk.cap312.6m.D1,hiv.bulk.cap312.6m.D2,hiv.bulk.cap312.6m.G1,hiv.bulk.cap312.6m.G2)
hiv.bulk.cap312.6m <- hiv.bulk.cap312.6m  %>% filter(Redundancy > 1)
#write.table(hiv.bulk.cap312.3y, "hiv_bulk_cap312_3y.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(hiv.bulk.cap312.6m, "hiv_bulk_cap312_6m.tab", sep = "\t", row.names = FALSE, quote = FALSE)


hiv.bulk.cap322.3y.b1 <- read_csv("SRR6207011_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap322.3y.b1 <- hiv.bulk.cap322.3y.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap322.3y.b2 <- read_csv("SRR6206372_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap322.3y.b2 <- hiv.bulk.cap322.3y.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap322.3y.A1 <- read_csv("SRR6206372_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap322.3y.A2 <- read_csv("SRR6207011_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap322.3y.D1 <- read_csv("SRR6206372_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap322.3y.G1 <- read_csv("SRR6207011_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap322.3y.G2 <- read_csv("SRR6206372_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap322.3y.M1 <- read_csv("SRR6206372_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap322.3y.M2 <- read_csv("SRR6207011_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap322.6m.b1 <- read_csv("SRR6206376_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap322.6m.b1 <- hiv.bulk.cap322.6m.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap322.6m.b2 <- read_csv("SRR6207019_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap322.6m.b2 <- hiv.bulk.cap322.6m.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap322.6m.A1 <- read_csv("SRR6207019_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap322.6m.A2 <- read_csv("SRR6206376_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap322.6m.G1 <- read_csv("SRR6206376_Heavy_IGHG.csv.gz", skip =1)

hiv.bulk.cap322.3y <- rbind(hiv.bulk.cap322.3y.b1,hiv.bulk.cap322.3y.b2,hiv.bulk.cap322.3y.A1,hiv.bulk.cap322.3y.A2,hiv.bulk.cap322.3y.D1,hiv.bulk.cap322.3y.G1,hiv.bulk.cap322.3y.G2,hiv.bulk.cap322.3y.M1,hiv.bulk.cap322.3y.M2)
hiv.bulk.cap322.3y <- hiv.bulk.cap322.3y  %>% filter(Redundancy > 1)
hiv.bulk.cap322.6m <- rbind(hiv.bulk.cap322.6m.b1,hiv.bulk.cap322.6m.b2,hiv.bulk.cap322.6m.A1,hiv.bulk.cap322.6m.A2,hiv.bulk.cap322.6m.G1)
hiv.bulk.cap322.6m <- hiv.bulk.cap322.6m  %>% filter(Redundancy > 1)
#write.table(hiv.bulk.cap322.3y, "hiv_bulk_cap322_3y.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(hiv.bulk.cap322.6m, "hiv_bulk_cap322_6m.tab", sep = "\t", row.names = FALSE, quote = FALSE)

hiv.bulk.cap335.3y.b1 <- read_csv("SRR6207008_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap335.3y.b1 <- hiv.bulk.cap335.3y.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap335.3y.b2 <- read_csv("SRR6206381_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap335.3y.b2 <- hiv.bulk.cap335.3y.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap335.3y.A1 <- read_csv("SRR6207008_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap335.3y.A2 <- read_csv("SRR6206381_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap335.3y.D1 <- read_csv("SRR6206381_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap335.3y.G1 <- read_csv("SRR6206381_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap335.3y.G2 <- read_csv("SRR6207008_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap335.3y.M1 <- read_csv("SRR6207008_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap335.3y.M2 <- read_csv("SRR6206381_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap335.6m.b1 <- read_csv("SRR6206373_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap335.6m.b1 <- hiv.bulk.cap335.6m.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap335.6m.b2 <- read_csv("SRR6207004_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap335.6m.b2 <- hiv.bulk.cap335.6m.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap335.6m.A1 <- read_csv("SRR6206373_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap335.6m.A2 <- read_csv("SRR6207004_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap335.6m.D1 <- read_csv("SRR6207004_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap335.6m.D2 <- read_csv("SRR6206373_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap335.6m.G1 <- read_csv("SRR6206373_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap335.6m.M1 <- read_csv("SRR6207004_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap335.3y.M2 <- read_csv("SRR6206373_Heavy_IGHM.csv.gz", skip =1)

hiv.bulk.cap335.3y <- rbind(hiv.bulk.cap335.3y.b1,hiv.bulk.cap335.3y.b2,hiv.bulk.cap335.3y.A1,hiv.bulk.cap335.3y.A2,hiv.bulk.cap335.3y.D1,hiv.bulk.cap335.3y.G1,hiv.bulk.cap335.3y.G2,hiv.bulk.cap335.3y.M1,hiv.bulk.cap335.3y.M2)
hiv.bulk.cap335.3y <- hiv.bulk.cap335.3y  %>% filter(Redundancy > 1)
hiv.bulk.cap335.6m <- rbind(hiv.bulk.cap335.6m.b1,hiv.bulk.cap335.6m.b2,hiv.bulk.cap335.6m.A1,hiv.bulk.cap335.6m.A2,hiv.bulk.cap335.6m.D1,hiv.bulk.cap335.6m.D2,hiv.bulk.cap335.6m.G1,hiv.bulk.cap335.6m.M1,hiv.bulk.cap335.3y.M2)
hiv.bulk.cap335.6m <- hiv.bulk.cap335.6m  %>% filter(Redundancy > 1)
#write.table(hiv.bulk.cap335.3y, "hiv_bulk_cap335_3y.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(hiv.bulk.cap335.6m, "hiv_bulk_cap335_6m.tab", sep = "\t", row.names = FALSE, quote = FALSE)

hiv.bulk.cap351.3y.b1 <- read_csv("SRR6207005_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap351.3y.b1 <- hiv.bulk.cap351.3y.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap351.3y.b2 <- read_csv("SRR6206370_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap351.3y.b2 <- hiv.bulk.cap351.3y.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap351.3y.A1 <- read_csv("SRR6206370_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap351.3y.A2 <- read_csv("SRR6207005_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap351.3y.D1 <- read_csv("SRR6206370_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap351.3y.G1 <- read_csv("SRR6206370_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap351.3y.G2 <- read_csv("SRR6207005_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap351.3y.M1 <- read_csv("SRR6207005_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap351.3y.M2 <- read_csv("SRR6206370_Heavy_IGHM.csv.gz", skip =1)
hiv.bulk.cap351.6m.b1 <- read_csv("SRR6206386_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap351.6m.b1 <- hiv.bulk.cap351.6m.b1  %>% filter(Redundancy > 1)
hiv.bulk.cap351.6m.b2 <- read_csv("SRR6207007_Heavy_Bulk.csv.gz", skip =1)
hiv.bulk.cap351.6m.b2 <- hiv.bulk.cap351.6m.b2  %>% filter(Redundancy > 1)
hiv.bulk.cap351.6m.A1 <- read_csv("SRR6207007_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap351.6m.A2 <- read_csv("SRR6206386_Heavy_IGHA.csv.gz", skip =1)
hiv.bulk.cap351.6m.D1 <- read_csv("SRR6207007_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap351.6m.D2 <- read_csv("SRR6206386_Heavy_IGHD.csv.gz", skip =1)
hiv.bulk.cap351.6m.G1 <- read_csv("SRR6207007_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap351.6m.G2 <- read_csv("SRR6206386_Heavy_IGHG.csv.gz", skip =1)
hiv.bulk.cap351.6m.M1 <- read_csv("SRR6206386_Heavy_IGHM.csv.gz", skip =1)

hiv.bulk.cap351.3y <- rbind(hiv.bulk.cap351.3y.b1,hiv.bulk.cap351.3y.b2,hiv.bulk.cap351.3y.A1,hiv.bulk.cap351.3y.A2,hiv.bulk.cap351.3y.D1,hiv.bulk.cap351.3y.G1,hiv.bulk.cap351.3y.G2,hiv.bulk.cap351.3y.M2)
hiv.bulk.cap351.3y <- hiv.bulk.cap351.3y %>% filter(Redundancy > 1)
hiv.bulk.cap351.6m <- rbind(hiv.bulk.cap351.6m.b1,hiv.bulk.cap351.6m.b2,hiv.bulk.cap351.6m.A1,hiv.bulk.cap351.6m.A2,hiv.bulk.cap351.6m.D1,hiv.bulk.cap351.6m.D2,hiv.bulk.cap351.6m.G1,hiv.bulk.cap351.6m.G2,hiv.bulk.cap351.6m.M1)
hiv.bulk.cap351.6m <- hiv.bulk.cap351.6m  %>% filter(Redundancy > 1)
#write.table(hiv.bulk.cap351.3y, "hiv_bulk_cap351_3y.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(hiv.bulk.cap351.6m, "hiv_bulk_cap351_6m.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## these 12 saved tab files "hiv.bulk.cap..." are already available in the /paper_assets/intermediate_files folder
# hiv.bulk.cap287.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap287_3y.tab.gz")
# hiv.bulk.cap287.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap287_6m.tab.gz")
# hiv.bulk.cap301.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap301_3y.tab.gz")
# hiv.bulk.cap301.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap301_6m.tab.gz")
# hiv.bulk.cap312.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap312_3y.tab.gz")
# hiv.bulk.cap312.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap312_6m.tab.gz")
# hiv.bulk.cap322.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap322_3y.tab.gz")
# hiv.bulk.cap322.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap322_6m.tab.gz")
# hiv.bulk.cap335.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap335_3y.tab.gz")
# hiv.bulk.cap335.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap335_6m.tab.gz")
# hiv.bulk.cap351.3y <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap351_3y.tab.gz")
# hiv.bulk.cap351.6m <- read_tsv("paper_assets/intermediate_files/hiv_bulk_cap351_6m.tab.gz")

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING HIV ANTIBODY SEQUENCES FROM IEDB, CATNAP, AND YACOOB ET AL.

### first download IEDB data - searching all B cell human data
## extracted all mAbs with organism matching "HIV" or "human immunodeficiency virus", removing duplicates
## code to convert CoV-AbDab columns to AIRRscape format
hiv.iedb <- read_csv("bcell_receptor_table_export_xxyyzz.csv")

hiv.iedb.hiv1 <- hiv.iedb[ grep("HIV", hiv.iedb$Organism) , ]
hiv.iedb.hiv2 <- hiv.iedb[ grep("immunodeficiency virus", hiv.iedb$Organism) , ]
hiv.iedb <- rbind(hiv.iedb.hiv1,hiv.iedb.hiv2)
rm(hiv.iedb.hiv1)
rm(hiv.iedb.hiv2)
##add HIV-IEDBmAb- to names
hiv.iedb$toadd <- "HIV-IEDBmAb"
hiv.iedb <- hiv.iedb %>% unite(sequence_id, toadd, 'Reference Name', sep = "-", remove = FALSE, na.rm = TRUE)

## remove duplicates `Chain 1 CDR3 Curated`
hiv.iedb <- hiv.iedb %>% rename(CDRH3 = 'Chain 1 CDR3 Calculated')

hiv.iedb <- hiv.iedb %>%
   group_by(sequence_id,CDRH3) %>%
   summarize_all(first)

## split HC & LC to separate heavy & light v gene, j gene, and full sequence data
hiv.iedb.h <- hiv.iedb
hiv.iedb.l <- hiv.iedb

hiv.iedb.h <- hiv.iedb.h %>% rename(v_call = `Calculated Chain 1 V Gene`)
hiv.iedb.h <- hiv.iedb.h %>% rename(j_call = `Calculated Chain 1 J Gene`)
hiv.iedb.l <- hiv.iedb.l %>% rename(v_call = `Calculated Chain 2 V Gene`)
hiv.iedb.l <- hiv.iedb.l %>% rename(j_call = `Calculated Chain 2 J Gene`)

hiv.iedb.h <- hiv.iedb.h %>% rename(fullv = `Chain 1 Full Sequence`)
hiv.iedb.l <- hiv.iedb.l %>% rename(fullv = `Chain 2 Full Sequence`)
hiv.iedb.h <- hiv.iedb.h %>% rename(cdr3_aa_imgt = CDRH3)
hiv.iedb.l <- hiv.iedb.l %>% rename(cdr3_aa_imgt = 'Chain 2 CDR3 Calculated')
hiv.iedb.l$CDRH3 <- NULL

## renaming sequence_id, then rejoin hc & lc files
hiv.iedb.l <- hiv.iedb.l %>% rename(sequence_id0 = sequence_id)
hiv.iedb.l$toadd2 <- "LC"
hiv.iedb.l <- hiv.iedb.l %>% unite(sequence_id, sequence_id0, toadd2, sep = "-", remove = FALSE, na.rm = TRUE)
hiv.iedb.l$sequence_id0 <- NULL

hiv.iedb <- full_join(hiv.iedb.h, hiv.iedb.l)
rm(hiv.iedb.h)
rm(hiv.iedb.l)

### want most of commands from shinyprocess function - but can't just rerun because no junction_aa
## subset of those commands:
hiv.iedb$sequence_id <- gsub("\\_","\\-",hiv.iedb$sequence_id)
hiv.iedb$sequence_id <- gsub("\\.","\\-",hiv.iedb$sequence_id)
hiv.iedb$sequence_id <- gsub("\\+","\\-",hiv.iedb$sequence_id)
hiv.iedb$sequence_id <- gsub("\\ ","\\-",hiv.iedb$sequence_id)
hiv.iedb$cdr3length_imgt <- nchar(hiv.iedb$cdr3_aa_imgt)
## next lines create V gene family, J gene columns
hiv.iedb$vgene <- getGene(hiv.iedb$v_call, first=TRUE, strip_d=TRUE)
hiv.iedb$vgf <- substring(hiv.iedb$vgene, 1,5)
hiv.iedb$jgene <- getGene(hiv.iedb$j_call, first=TRUE, strip_d=TRUE)
## this creates new column vgf_jgene which is used in all shiny plots
hiv.iedb <- hiv.iedb %>% unite(vgf_jgene, vgf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)
## this removes any rows without CDR3, or with junctions that are not 3-mers
hiv.iedb <- hiv.iedb %>% filter(!is.na(cdr3length_imgt)) %>% 
   filter(is.wholenumber(cdr3length_imgt))

## this will filter the dataset to AIRRscape-specific columns (plus fullv sequence)
vars2 <- c("sequence_id", "binding", "neutralization", "cregion", "cdr3_aa_imgt","vgene", "vgf_jgene", "vgf","jgene", "cdr3length_imgt", "shm", "shm_max", "shm_mean", "ncount", "fullv")

hiv.iedb.fullv <- hiv.iedb %>% select(any_of(vars2))

### also remove sequences with ND's in cdr3
hiv.iedb.fullv <- hiv.iedb.fullv %>% filter(cdr3_aa_imgt != "NA")

## add cregion
hiv.iedb.fullv$cregion <- str_sub(hiv.iedb.fullv$vgene, end=3)
hiv.iedb.fullv$cregion <- gsub('IG','Ig',hiv.iedb.fullv$cregion)

#write.table(hiv.iedb.fullv, "hiv_iedb_withfullv.tab", sep = "\t", row.names = FALSE, quote = FALSE)

#######
# our resulting file is hiv_iedb_withfullv.tab - after manually removing some duplicates have 96 sequences, all with full v-gene
## for AIRRscape visualization of SHM, we searched Genbank using tblastn (in Geneious) to search for 100% sequence match to nucleotide sequences
## subset with matches are run with commands from the Immcantation docker container
## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
## we used (docker pull immcantation/suite:4.1.0)
# changeo-igblast -s /data/hiv_iedb0.fasta -n hiv_iedb0 -o /data	
# resulting file is hiv_iedb0_germ-pass.tsv - result are 149 with SHM
# toshiny_hiv_iedb.tab - 194 sequences, 149 with SHM

#################################################
## next adding CATNAP - sequences from: https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/heavy_seqs_na.fasta
## and https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/light_seqs_na.fasta
## then combine into one fasta file - CATNAP_seqs_nt.fasta

## running these commands from the Immcantation docker container
## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
## we used (docker pull immcantation/suite:4.1.0)
# changeo-igblast -s /data/CATNAP_seqs_nt.fasta -n CATNAP_seqs -o /data
# changeo-clone -d /data/CATNAP_seqs_db-pass.tsv -n CATNAP_seqs -o /data -x 0.12
## resulting file is CATNAP_seqs_germ-pass.tsv

hiv.mabs.catnap <- read_tsv("CATNAP_seqs_germ-pass.tsv")
toshiny.hiv.mabs.catnap <- AIRRscapeprocess(hiv.mabs.catnap, renumber_sequences = FALSE, filter_after_counting = FALSE)  # shinyprocess function is at bottom of this script
#write.table(toshiny.hiv.mabs.catnap, "toshiny_hiv_mabs_catnap.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## then manually combined IEDB mabs with CATNAP
## resulting file is toshiny_hiv_mabs_all0.tab

#################################################
## lastly added Yacoob et al data - Genbank KX443243–KX443325
## running these commands from the Immcantation docker container
## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
## we used (docker pull immcantation/suite:4.1.0)
# cchangeo-igblast -s /data/yacoob_seqs_nt.fasta -n yacoob_seqs -o /data
# changeo-clone -d /data/yacoob_seqs_db-pass.tsv -n yacoob_seqs -o /data -x 0.12
## resulting file is yacoob_seqs_germ-pass.tsv

hiv.mabs.yacoob <- read_tsv("yacoob_seqs_germ-pass.tsv")
toshiny.hiv.mabs.yacoob <- AIRRscapeprocess(hiv.mabs.yacoob, renumber_sequences = FALSE, filter_after_counting = FALSE) # shinyprocess function is at bottom of this script
#write.table(toshiny.hiv.mabs.yacoob, "toshiny_hiv_mabs_yacoob.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## then combine IEDB, CATNAP & Yacoob mabs
## resulting file is toshiny_hiv_mabs_all.tab - 622 sequences total

## this file is available in the /paper_assets/intermediate_files folder
# toshiny.hiv.mabs.all <- read_tsv("paper_assets/intermediate_files/toshiny_hiv_mabs_all.tab")
## to add newer byvgene stats...
toshiny.hiv.mabs.all$ncount <- NULL
toshiny.hiv.mabs.all$shm_mean <- NULL
toshiny.hiv.mabs.all$shm_max <- NULL
toshiny.hiv.mabs.all <- toshiny.hiv.mabs.all %>%
  add_count(vgf_jgene,cdr3length_imgt) %>% 
  rename(ncount = n) %>% 
  ungroup() %>% 
  add_count(vgene,cdr3length_imgt) %>% 
  rename(ncount_vgene = n) %>% 
  ungroup() %>% 
  group_by(vgene,cdr3length_imgt) %>% 
  mutate(shm_byvgene_mean = mean(shm, na.rm = TRUE)) %>% 
  # NOTE ADDING MAX & MEAN SHM BY V-GENE
  mutate(shm_byvgene_max = max(shm, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(vgf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # NOTE ADDIN MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>%
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>%
  mutate(shm_max = na_if(shm_max, "-Inf")) %>%
  mutate(across(shm, round, 2)) %>%
  mutate(across(shm_max, round, 2)) %>%
  mutate(across(shm_mean, round, 2)) %>%
  mutate(shm_byvgene_mean = na_if(shm_byvgene_mean, "NaN")) %>%
  mutate(shm_byvgene_max = na_if(shm_byvgene_max, "-Inf")) %>%
  mutate(across(shm_byvgene_max, round, 2)) %>%
  mutate(across(shm_byvgene_mean, round, 2))

#write.table(toshiny.hiv.mabs.all, "toshiny_hiv_mabs_all.tab", sep = "\t", row.names = FALSE, quote = FALSE)

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING DENV ANTIBODY SEQUENCES FROM ZANINI PAPER

# data directly from Zanini SOM, then run through Immcantation
# file is available in the /paper_assets/intermediate_files folder - toshiny_den_mabs0.fasta
## converting to fasta, running these commands from the Immcantation docker container
## follow commands on Immcantation website https://immcantation.readthedocs.io/en/stable/docker/intro.html 
## we used (docker pull immcantation/suite:4.1.0)
# changeo-igblast -s /data/toshiny_den_mabs0.fasta -n toshiny_den_mabs0 -o /data	
# note using 20% threshold to match Zanini		
# changeo-clone -d /data/toshiny_den_mabs0_db-pass.tsv -n toshiny_den_mabs0 -o /data -x 0.2	

# resulting file is toshiny_den_mabs0_germ-pass.tsv - 79 sequences
## this file is already available in the /paper_assets/intermediate_files folder
# den.mabs <- read_tsv("paper_assets/intermediate_files/toshiny_den_mabs0_germ-pass.tsv")

## then run shinyprocess - note THIS IS IN THE AIRRSCAPE_PROCESSING SCRIPT
#toshiny.den.mabs <- shinyprocess(den.mabs, renumber_sequences = FALSE, filter_after_counting = FALSE)

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING SARS-COV2 ANTIBODY SEQUENCES FROM CoV-AbDab
## first download full database as CSV (e.g. http://opig.stats.ox.ac.uk/webapps/covabdab/static/downloads/CoV-AbDab_xxyyzz.csv)
### version 1.0 of AIRRscape uses CoV-AbDab from CoV-AbDab_090721 database snapshot

## this list includes all CoV antibodies, and many non-human
## focus is protein sequence, and some entries are not complete (i.e. not full length sequence)
## list is first trimmed to antibodies  binding to SARS-CoV-2, and derived from humans or humanized mice

## code to convert CoV-AbDab columns to AIRRscape format
cov2.abdab <- read_csv("CoV-AbDab_xxyyzz.csv")

cov2.abdab <- cov2.abdab[ grep("SARS-CoV2", cov2.abdab$`Binds to`) , ]
cov2.abdab <- cov2.abdab[ grep("Human", cov2.abdab$`Heavy V Gene`) , ]
cov2.abdab <- subset(cov2.abdab, `Ab or Nb` %in% c("Ab"))
cov2.abdab <- cov2.abdab[ grep("ND", cov2.abdab$`Light V Gene`, invert = TRUE) , ]

cov2.abdab.neut <- cov2.abdab[ grep("SARS-CoV2", cov2.abdab$`Neutralising Vs`) , ]
cov2.abdab.nonneut <- cov2.abdab[ grep("SARS-CoV2", cov2.abdab$`Neutralising Vs`, invert = TRUE) , ]
cov2.abdab.neut$neutralization <- "neutralizing"
cov2.abdab.nonneut$neutralization <- "non-neutralizing or unknown"
cov2.abdab <- bind_rows(cov2.abdab.neut, cov2.abdab.nonneut)
cov2.abdab.rbd <- cov2.abdab %>% filter(`Protein + Epitope` == "S; RBD")
cov2.abdab.nonrbd <- cov2.abdab %>% filter(`Protein + Epitope` != "S; RBD")
cov2.abdab.rbd$binding <- "RBD"
cov2.abdab.nonrbd$binding <- "non-RBD"
cov2.abdab <- bind_rows(cov2.abdab.rbd, cov2.abdab.nonrbd)
rm(cov2.abdab.neut)
rm(cov2.abdab.nonneut)
rm(cov2.abdab.rbd)
rm(cov2.abdab.nonrbd)

cov2.abdab$`Heavy V Gene` <- gsub(" [(]Human[)]","",cov2.abdab$`Heavy V Gene`)
cov2.abdab$`Heavy J Gene` <- gsub(" [(]Human[)]","",cov2.abdab$`Heavy J Gene`)
cov2.abdab$`Light V Gene` <- gsub(" [(]Human[)]","",cov2.abdab$`Light V Gene`)
cov2.abdab$`Light J Gene` <- gsub(" [(]Human[)]","",cov2.abdab$`Light J Gene`)
cov2.abdab$dataset <- "SARS-CoV2-mAb"

## split HC & LC to separate heavy & light v gene, j gene, and full sequence data
cov2.abdab.h <- cov2.abdab
cov2.abdab.l <- cov2.abdab

cov2.abdab.h <- cov2.abdab.h %>% rename(v_call = `Heavy V Gene`)
cov2.abdab.h <- cov2.abdab.h %>% rename(j_call = `Heavy J Gene`)
cov2.abdab.l <- cov2.abdab.l %>% rename(v_call = `Light V Gene`)
cov2.abdab.l <- cov2.abdab.l %>% rename(j_call = `Light J Gene`)

cov2.abdab.h <- cov2.abdab.h %>% rename(fullv = `VH or VHH`)
cov2.abdab.l <- cov2.abdab.l %>% rename(fullv = VL)
cov2.abdab.h <- cov2.abdab.h %>% rename(cdr3_aa_imgt = CDRH3)
cov2.abdab.l <- cov2.abdab.l %>% rename(cdr3_aa_imgt = CDRL3)

## renaming sequence_id, then rejoin hc & lc files
cov2.abdab.l$toadd <- "LC"
cov2.abdab.l <- cov2.abdab.l %>% unite(sequence_id0, Name, toadd, sep = "_", remove = FALSE, na.rm = TRUE)

cov2.abdab.h <- cov2.abdab.h %>% unite(sequence_id, dataset, Name, sep = "_", remove = FALSE, na.rm = TRUE)
cov2.abdab.l <- cov2.abdab.l %>% unite(sequence_id, dataset, sequence_id0, sep = "_", remove = FALSE, na.rm = TRUE)

cov2.abdab <- full_join(cov2.abdab.h, cov2.abdab.l)
rm(cov2.abdab.h)
rm(cov2.abdab.l)

### also remove sequences with ND's in jgene (or manually annotate)
#cov2.abdab <- cov2.abdab %>% filter(jgene != "N")

### want most of commands from shinyprocess function - but can't just rerun because no junction_aa
## subset of those commands:
cov2.abdab$sequence_id <- gsub("\\_","\\-",cov2.abdab$sequence_id)
cov2.abdab$cdr3length_imgt <- nchar(cov2.abdab$cdr3_aa_imgt)
### removing all sequences with IMGT CDR3 less than 3
cov2.abdab <- cov2.abdab %>% filter(cdr3length_imgt > 2.8)  
## next lines create V gene family, J gene columns
cov2.abdab$vgene <- getGene(cov2.abdab$v_call, first=TRUE, strip_d=TRUE)
cov2.abdab$vgf <- substring(cov2.abdab$vgene, 1,5)
cov2.abdab$jgene <- getGene(cov2.abdab$j_call, first=TRUE, strip_d=TRUE)
## this creates new column vgf_jgene which is used in all shiny plots
cov2.abdab <- cov2.abdab %>% unite(vgf_jgene, vgf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)
## this removes any rows without CDR3, or with junctions that are not 3-mers
cov2.abdab <- cov2.abdab %>% filter(!is.na(cdr3length_imgt)) %>% 
   filter(is.wholenumber(cdr3length_imgt))

cov2.abdab <- cov2.abdab %>%
   add_count(vgf_jgene,cdr3length_imgt) %>% 
   rename(ncount = n) %>% 
   group_by(vgf_jgene,cdr3length_imgt)
## this will filter the dataset to AIRRscape-specific columns (plus fullv sequence)
vars2 <- c("sequence_id", "binding", "neutralization", "cregion", "cdr3_aa_imgt","vgene", "vgf_jgene", "vgf","jgene", "cdr3length_imgt", "shm", "shm_max", "shm_mean", "ncount", "fullv")

cov2.abdab.fullv0 <- cov2.abdab %>% select(any_of(vars2))
cov2.abdab.fullv <- cov2.abdab.fullv0 %>% filter(fullv != "ND")
#write.table(cov2.abdab.fullv, "cov2_abdab_withfullv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
   
#######
# our resulting file is cov2_abdab_withfullv.tab - 3,385 of 4,306 total sequences have full vgene sequence
## for AIRRscape visualization of SHM, we searched Genbank using tblastn (in Geneious) to search for 100% sequence match to nucleotide sequences
## subset with matches are run with Immcantation:changeo-igblast
# changeo-igblast -s /data/cov2_abdab0.fasta -n cov2_abdab0 -o /data	
# resulting file is cov2_abdab0_germ-pass.tsv

## these were then checked to remove any with artificially high SHM due to codon optimization (e.g. Seydoux et al., 2020 https://www.biorxiv.org/content/10.1101/2020.05.12.091298v1)
cov2.abdab.withfullv <- read_tsv("paper_assets/intermediate_files/cov2_abdab_withfullv.tab")
cov2.abdab.fullvtoaddshm <- read_tsv("paper_assets/intermediate_files/cov2_abdab0_germ-pass.tsv")

## get 'binding' and 'neutralization' columns to add to cov2.abdab.fullvtoaddshm
cov2.abdab.withfullv.bindneut <- cov2.abdab.withfullv %>% select(sequence_id,binding,neutralization)
cov2.abdab.fullvtoaddshm <- left_join(cov2.abdab.fullvtoaddshm, cov2.abdab.withfullv.bindneut)


## then run shinyprocess - if not already loaded, function is at bottom...
toshiny.cov2.abdab.fullvtoaddshm <- AIRRscapeprocess(cov2.abdab.fullvtoaddshm, renumber_sequences = FALSE)

toshiny.cov2.abdab.fullvtoaddshm <- toshiny.cov2.abdab.fullvtoaddshm %>% relocate(sequence_id)
## need to recover 257 sequences lost during shinyprocess computation
cov2.abdab.257lost <- anti_join(cov2.abdab.fullvtoaddshm, toshiny.cov2.abdab.fullvtoaddshm)

cov2.abdab.257lost <- semi_join(cov2.abdab.withfullv, cov2.abdab.257lost)

## better to do a join that removes all rows in cov2.abdab.withfullv not in cov2.abdab.fullvtoaddshm - anti_join
cov2.abdab.noshm <- anti_join(cov2.abdab.withfullv, cov2.abdab.fullvtoaddshm)
## then add 257 lost sequences
cov2.abdab.noshm <- bind_rows(cov2.abdab.noshm, cov2.abdab.257lost)


toshiny.cov2.abdab.fullv <- full_join(toshiny.cov2.abdab.fullvtoaddshm, cov2.abdab.noshm)
# resulting search adds SHM data to 1,493 sequences - toshiny_cov2_abdab_withfullv.tab

## now some last clean up steps
toshiny.cov2.abdab.fullv$shm_max <- NULL
toshiny.cov2.abdab.fullv$shm_mean <- NULL
toshiny.cov2.abdab.fullv$ncount <- NULL
toshiny.cov2.abdab.fullv$order <- NULL

### change IgK & IgL to Kappa & Lambda
toshiny.cov2.abdab.fullv$cregion <- gsub("IgK","Kappa",toshiny.cov2.abdab.fullv$cregion)
toshiny.cov2.abdab.fullv$cregion <- gsub("IgL","Lambda",toshiny.cov2.abdab.fullv$cregion)

## ALSO NEED TO MAKE SURE SEQUENCE ID IS FIRST FOR MABS...AND CDR3 AFTER CREGION
toshiny.cov2.abdab.fullv <- toshiny.cov2.abdab.fullv %>% relocate(sequence_id)
toshiny.cov2.abdab.fullv <- toshiny.cov2.abdab.fullv %>% relocate(cdr3_aa_imgt, .after = cregion)

toshiny.cov2.abdab.fullv <- toshiny.cov2.abdab.fullv %>% relocate(neutralization, .after = cdr3_aa_imgt)
toshiny.cov2.abdab.fullv <- toshiny.cov2.abdab.fullv %>% relocate(binding, .after = cdr3_aa_imgt)
toshiny.cov2.abdab.fullv <- toshiny.cov2.abdab.fullv %>% relocate(vgf_jgene, .after = gene)

## now add count, shm_mean, shm_max after all sequences are combined
## adding byvgene stats...
toshiny.cov2.abdab <- toshiny.cov2.abdab.fullv %>%
  add_count(vgf_jgene,cdr3length_imgt) %>% 
  rename(ncount = n) %>% 
  ungroup() %>% 
  add_count(vgene,cdr3length_imgt) %>% 
  rename(ncount_vgene = n) %>% 
  ungroup() %>% 
  group_by(vgene,cdr3length_imgt) %>% 
  mutate(shm_byvgene_mean = mean(shm, na.rm = TRUE)) %>% 
  # NOTE ADDING MAX & MEAN SHM BY V-GENE
  mutate(shm_byvgene_max = max(shm, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(vgf_jgene,cdr3length_imgt) %>% 
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>% 
  # NOTE ADDIN MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>% 
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>% 
  mutate(shm_max = na_if(shm_max, "-Inf")) %>% 
  mutate(across(shm, round, 2)) %>% 
  mutate(across(shm_max, round, 2)) %>% 
  mutate(across(shm_mean, round, 2)) %>% 
  mutate(shm_byvgene_mean = na_if(shm_byvgene_mean, "NaN")) %>% 
  mutate(shm_byvgene_max = na_if(shm_byvgene_max, "-Inf")) %>% 
  mutate(across(shm_byvgene_max, round, 2)) %>% 
  mutate(across(shm_byvgene_mean, round, 2))

toshiny.cov2.abdab$fullv <- NULL

#write.table(toshiny.cov2.abdab.fullv, "toshiny_cov2_abdab_fullv.tab", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(toshiny.cov2.abdab, "toshiny_cov2_abdab.tab", sep = "\t", row.names = FALSE, quote = FALSE)

## these files are already available in the /paper_assets/intermediate_files or /shinyapp folder
# cov2.abdab.withfullv <- read_tsv("paper_assets/intermediate_files/cov2_abdab_withfullv.tab")
# toshiny.cov2.abdab.fullv <- read_tsv("paper_assets/intermediate_files/toshiny_cov2_abdab_fullv.tab")
# toshiny.cov2.abdab <- read_tsv("shinyapp/toshiny_cov2_abdab.tab")

## to add newer byvgene stats...
toshiny.cov2.abdab.fullv$ncount <- NULL
toshiny.cov2.abdab.fullv$shm_mean <- NULL
toshiny.cov2.abdab.fullv$shm_max <- NULL
toshiny.cov2.abdab.fullv <- toshiny.cov2.abdab.fullv %>%
  add_count(vgf_jgene,cdr3length_imgt) %>% 
  rename(ncount = n) %>% 
  ungroup() %>% 
  add_count(vgene,cdr3length_imgt) %>% 
  rename(ncount_vgene = n) %>% 
  ungroup() %>% 
  group_by(vgene,cdr3length_imgt) %>% 
  mutate(shm_byvgene_mean = mean(shm, na.rm = TRUE)) %>% 
  # NOTE ADDING MAX & MEAN SHM BY V-GENE
  mutate(shm_byvgene_max = max(shm, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(vgf_jgene,cdr3length_imgt) %>%
  mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
  # NOTE ADDIN MAX SHM AS WELL..
  mutate(shm_max = max(shm, na.rm = TRUE)) %>%
  mutate(shm_mean = na_if(shm_mean, "NaN")) %>%
  mutate(shm_max = na_if(shm_max, "-Inf")) %>%
  mutate(across(shm, round, 2)) %>%
  mutate(across(shm_max, round, 2)) %>%
  mutate(across(shm_mean, round, 2)) %>%
  mutate(shm_byvgene_mean = na_if(shm_byvgene_mean, "NaN")) %>%
  mutate(shm_byvgene_max = na_if(shm_byvgene_max, "-Inf")) %>%
  mutate(across(shm_byvgene_max, round, 2)) %>%
  mutate(across(shm_byvgene_mean, round, 2))

##############################################################################################################################
## shinyprocess function updated to AIRRscapeprocess
##############################################################################################################################

## two functions from https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na
not_all_na <- function(x) any(!is.na(x))
# not_any_na <- function(x) all(!is.na(x))

##### ADDING TCR CAPABILITY...note there are <20 TRBV genes, so need to change code a bit...
## note required columns...and cannot have any missing!!  v_call, j_call, v_identity, junction_aa - but note TR dataset will have NAs for v_identity
AIRRscapeprocess <- function(x, filter_columns = TRUE, filter_to_HC = TRUE, renumber_sequences = TRUE, filter_after_counting = TRUE) {
  colname <- substitute(x)
  ## drop_na from v_call, j_call, v_identity, junction_aa
  x <- x %>% drop_na(v_call) %>% drop_na(j_call) %>% drop_na(junction_aa)
  ## this calculates SHM but depending on whether v_identity is from 0 to 1 or 0 to 100
  ## be careful with including if v_identity, if already processed will not have this!! - think we never invoke unless it is a unprocessed file BUT CHECK...
  ## idea - instead of if bcr_repertoire have if v_identity isn't blank...  allmisscols <- sapply(dt, function(x) all(is.na(x) | x == '' ))
  ## improvement on IG vs. TR.... if (substring(x$v_call[1],1,2) == "IG") { ... but would then need else if == "TR"...
  if (!is.na(mean(x$v_identity))) {
    x <- x %>% drop_na(v_identity)
    if (mean(x$v_identity) < 1) {
      x$shm <- (100 - (x$v_identity * 100))
    } else {
      x$shm <- (100 - x$v_identity)
    }
    x <- x %>% filter(vj_in_frame != "FALSE") %>%
      filter(vj_in_frame != "F")
    ## was splitting by '\', but we want to keep these
    # vgene0.pieces <- strsplit(x$v_call,"[*]")
    # x$vgene0 <- sapply(vgene0.pieces, "[", 1)
    # vgene.pieces1 <- strsplit(x$vgene0,"[/]")
    # x$vgene1 <- sapply(vgene.pieces1, "[", 1)
    # vgene.pieces <- strsplit(x$vgene1,",")
    # x$vgene <- sapply(vgene.pieces, "[", 1)
    vgene0.pieces <- strsplit(x$v_call,"[*]")
    x$vgene0 <- sapply(vgene0.pieces, "[", 1)
    vgene.pieces <- strsplit(x$vgene0,",")
    x$vgene <- sapply(vgene.pieces, "[", 1)
    ## vgene family is vgf, need to also delimit by * or / not just -
    vgf.pieces0 <- strsplit(x$vgene,"-")
    x$vgf0 <- sapply(vgf.pieces0, "[", 1)
    vgf.pieces1 <- strsplit(x$vgf0,"[*]")
    x$vgf1 <- sapply(vgf.pieces1, "[", 1)
    vgf.pieces2 <- strsplit(x$vgf1,"S")
    x$vgf2 <- sapply(vgf.pieces2, "[", 1)
    vgf.pieces <- strsplit(x$vgf2,"[/]")
    x$vgf <- sapply(vgf.pieces, "[", 1)
    ## reverting to original code for jgene...
    # jgene.pieces <- strsplit(x$j_call,"-")
    # x$jgene <- sapply(jgene.pieces, "[", 1)
    # next try
    # x$vgene <- getGene(x$v_call, first=TRUE, strip_d=TRUE)
    # x$vgf <- getFamily(x$v_call, first=TRUE, strip_d=TRUE)
    # x$jgene <- getGene(x$j_call, first=TRUE, strip_d=TRUE)
    # first try
    # x$vgf <- substring(x$vgene, 1,5)
    # x$vgene <- getGene(x$v_call, first=TRUE, strip_d=TRUE)
    # x$vgf <- substring(x$vgene, 1,5)
    x$jgene <- getGene(x$j_call, first=TRUE, strip_d=TRUE)
    x$jgene <- substring(x$jgene, 1,5)
  } else {
    x$shm <- 0
    ## was splitting by '\', but we want to keep these
    # vgene0.pieces <- strsplit(x$v_call,"[*]")
    # x$vgene0 <- sapply(vgene0.pieces, "[", 1)
    # vgene.pieces <- strsplit(x$vgene0,"[/]")
    # x$vgene <- sapply(vgene.pieces, "[", 1)
    vgene.pieces <- strsplit(x$v_call,"[*]")
    x$vgene <- sapply(vgene.pieces, "[", 1)
    ## vgene family is vgf, need to also delimit by * or / not just -
    vgf.pieces0 <- strsplit(x$vgene,"-")
    x$vgf0 <- sapply(vgf.pieces0, "[", 1)
    vgf.pieces1 <- strsplit(x$vgf0,"[*]")
    x$vgf1 <- sapply(vgf.pieces1, "[", 1)
    vgf.pieces <- strsplit(x$vgf1,"[/]")
    x$vgf <- sapply(vgf.pieces, "[", 1)
    ## reverting to original code for jgene...FOR TCR INSTEAD NEED TO REMOVE BY * AND -
    # x$jgene <- getGene(x$j_call, first=TRUE, strip_d=TRUE)
    # x$jgene <- substring(x$jgene, 1,5)
    jgene.pieces0 <- strsplit(x$j_call,"-")
    x$jgene0 <- sapply(jgene.pieces0, "[", 1)
    jgene.pieces <- strsplit(x$jgene0,"[*]")
    x$jgene <- sapply(jgene.pieces, "[", 1)
    ## next try
    # x$vgene <- getGene(x$v_call, first=TRUE, strip_d=TRUE)
    # x$vgf <- getFamily(x$v_call, first=TRUE, strip_d=TRUE)
    # x$jgene <- getGene(x$j_call, first=TRUE, strip_d=TRUE)
    ## first try
    # x$vgf <- substring(x$vgene, 1,5)
    # x$vgene <- x$v_gene
    # x$vgf <- x$v_subgroup
    # x$jgene <- x$j_subgroup
  }
  ## this removes columns with all NAs - moving to near bottom - tried moving up but instead back to near bottom and adding statements like this if needed: if (!is.na(mean(x$v_identity))) { 
  ## trying here again...latest - using new code
  x <- x %>% select(where(not_all_na))
  # x <- x[!map_lgl(x, ~ all(is.na(.)))]
  ## this makes new standard cdr3 column (sometimes already exists, but there should always be a junction_aa) by removing both ends of the junction_aa column
  x$cdr3_aa_imgt <- x$junction_aa
  str_sub(x$cdr3_aa_imgt, -1, -1) <- ""
  str_sub(x$cdr3_aa_imgt, 1, 1) <- ""
  ## this calculates the CDR3 length
  x$cdr3length_imgt <- nchar(x$cdr3_aa_imgt)
  ## removing non-productive, out of frame, stop codons, any X in CDR3 (was 'nnnn' in sequence)
  x <- x %>% filter(productive != "FALSE") %>%
    filter(productive != "F")
  x <- x[ grep("\\*", x$junction_aa, invert = TRUE) , ]
  x <- x[ grep("\\X", x$junction_aa, invert = TRUE) , ]
  ### removing all sequences with IMGT CDR3 less than 3
  x <- x %>% filter(cdr3length_imgt > 2.8)
  ## next lines create V gene family, J gene columns - making more universal and compatible with TCRs...see above
  # x$vgene <- unlist(strsplit(x$v_call, "[*]"))
  # x$vgf <- unlist(strsplit(x$v_call, "-"))
  # x$jgene <- unlist(strsplit(x$j_call, "-"))
  # x$vgene <- getGene(x$v_call, first=TRUE, strip_d=TRUE)
  # x$vgf <- substring(x$vgene, 1,5)
  # x$jgene <- getGene(x$j_call, first=TRUE, strip_d=TRUE)
  # x$jgene <- substring(x$jgene, 1,5)
  ## this creates new column vgf_jgene which is used in all shiny plots
  x <- x %>% unite(vgf_jgene, vgf, jgene, sep = "_", remove = FALSE, na.rm = TRUE)
  ## this removes any rows without CDR3, or with junctions that are not 3-mers
  x <- x %>% filter(!is.na(cdr3length_imgt)) %>%
    filter(is.wholenumber(cdr3length_imgt))
  # if there is a clone_id column this will make a count of reads_per_clone
  if ("clone_id" %in% names(x)) {
    x <- x %>% add_count(clone_id) %>%
      rename(reads_per_clone = n)
  }
  ## if no cregion column, make one  !(x %in% y) - moving substitution from IG to Ig outside of the if to always have it...
  ## adding first check to grab from c_call if it exists,   if (!is.na(mean(x$v_identity))) {  - better??  if(str == "HELLO" | str == "hello"){
  ## latest: first if no cregion, or if cregion is mostly NAs - removing second part...if (!("cregion" %in% names(x)) | (is.na(Mode(x$cregion)))) {
  if (!("cregion" %in% names(x))) {
    ## then use c_call if it exists & is not mostly NAs - if not, then just take first 3 characters of v_call
    ## if any NAs, need to replace!! (do last)
    if (("c_call" %in% names(x))) {
      x$cregion <- x$c_call
      # x$cregion <- x$c_call %>% replace_na(Mode(x$c_call))
      x$cregion <- gsub('IGH','IG',x$cregion)
      x$cregion <- gsub('TRAC','TRA',x$cregion)
      x$cregion <- gsub('TRBC','TRB',x$cregion)
      x$cregion <- gsub('TRDC','TRD',x$cregion)
      x$cregion <- gsub('TRGC','TRG',x$cregion)
      x$cregion <- gsub('IGMC','IGM',x$cregion)
      x$cregion <- gsub('IGDC','IGD',x$cregion)
      x$cregion <- gsub('IGGC','IGG',x$cregion)
      x$cregion <- gsub('IGAC','IGA',x$cregion)
      x$cregion <- gsub('IGKC','IGK',x$cregion)
      x$cregion <- gsub('IGLC','IGL',x$cregion)
      ## need to add a replace_na - instead using new function above
      # x$cregion <- x$cregion %>% replace_na(str_sub(x$v_call, end=3))
      # x$cregion <- gsub('IG','Ig',x$cregion)  ...substring(x$vgf[1],1,2) == "TR"
    } else {
      x$cregion <- str_sub(x$v_call, end=3)
    }
  }
  ## for new coalesce calculation below, make sure there is a locus field
  if (!("locus" %in% names(x))) {
    x$locus <- str_sub(x$v_call, end=3)
  }
  ## one option for replacing NAs - try coalesce per https://stackoverflow.com/questions/34071875/replace-a-value-na-with-the-value-from-another-column-in-r
  # dfABy %>% mutate(A = coalesce(A,B)) this is replacing all NAs in cregion with the value from locus
  x <- x %>% mutate(cregion = coalesce(cregion,locus))
  ## then changing names slightly next
  x$cregion <- gsub('IG','Ig',x$cregion)
  ## adding conversion to standardize light chain naming if necessary...
  x$cregion <- gsub("IgK","Kappa",x$cregion)
  x$cregion <- gsub("IgL","Lambda",x$cregion)
  ## another option for replacing NAs but loses TRA vs TRB etc.
  # if (substring(x$vgf[1],1,2) == "TR") {
  #   x$cregion <- x$cregion %>% replace_na("TR")
  # } else {
  #   x$cregion <- x$cregion %>% replace_na("IgH")
  # }
  ## after these defining steps of cregion, now replace any NAs after!
  ## making more important columns used in plotting, also a rounding step - later adding ncount_vgene as well
  x <- x %>%
    add_count(vgf_jgene,cdr3length_imgt) %>%
    rename(ncount = n) %>%
    ungroup() %>%
    add_count(vgene,cdr3length_imgt) %>%
    rename(ncount_vgene = n) %>%
    ungroup() %>%
    group_by(vgene,cdr3length_imgt) %>%
    mutate(shm_byvgene_mean = mean(shm, na.rm = TRUE)) %>%
    # NOTE ADDING MAX & MEAN SHM BY V-GENE
    mutate(shm_byvgene_max = max(shm, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(vgf_jgene,cdr3length_imgt) %>%
    mutate(shm_mean = mean(shm, na.rm = TRUE)) %>%
    # NOTE ADDIN MAX SHM AS WELL..
    mutate(shm_max = max(shm, na.rm = TRUE)) %>%
    mutate(shm_mean = na_if(shm_mean, "NaN")) %>%
    mutate(shm_max = na_if(shm_max, "-Inf")) %>%
    mutate(across(shm, round, 2)) %>%
    mutate(across(shm_max, round, 2)) %>%
    mutate(across(shm_mean, round, 2)) %>%
    mutate(shm_byvgene_mean = na_if(shm_byvgene_mean, "NaN")) %>%
    mutate(shm_byvgene_max = na_if(shm_byvgene_max, "-Inf")) %>%
    mutate(across(shm_byvgene_max, round, 2)) %>%
    mutate(across(shm_byvgene_mean, round, 2))
  ## this removes columns with all NAs - moving to near bottom - tried moving up but instead back to near bottom and adding statements like this if needed: if (!is.na(mean(x$v_identity))) { 
  ## trying further up again...
  # x <- x[!map_lgl(x, ~ all(is.na(.)))]
  ## this will filter the dataset if filter_columns option is set to true - note the any_of which allows columns to be missing - later adding various 10x column outputs
  # vars2 <- c("sequence_id", "binding", "neutralization", "cregion", "cdr3_aa_imgt","vgene", "vgf_jgene", "vgf","jgene", "cdr3length_imgt", "shm", "shm_max", "shm_mean", "ncount", "reads_per_clone")
  vars2 <- c("sequence_id", "nCount_RNA", "nFeature_RNA", "sample", "animal", "timepoint", "corrected_timepoint", "new_id", "binding", "neutralization", "cregion", "cdr3_aa_imgt","vgene", "vgf_jgene", "vgf","jgene", "cdr3length_imgt", "umi_count", "consensus_count", "duplicate_count", "clone_id", "shm", "shm_mean", "shm_max", "shm_byvgene_mean", "shm_byvgene_max", "ncount", "ncount_vgene", "reads_per_clone")
  if (filter_columns) {
    x <- x %>% select(any_of(vars2))
  }
  ## this will filter to heavy chains only......ADDING TRA & TRB
  if (filter_to_HC) {
    # x <- x %>% filter(cregion == "IgH" | cregion == "IgA" | cregion == "IgD" | cregion == "IgE" | cregion == "IgG" | cregion == "IgM")
    # x <- x %>% filter(cregion == "IgH" | cregion == "IgA" | cregion == "IgD" | cregion == "IgE" | cregion == "IgG" | cregion == "IgM" | cregion == "IgA1" | cregion == "IgA2" | cregion == "IgG1" | cregion == "IgG2" | cregion == "IgG3" | cregion == "IgG4" | cregion == "IgG2A" | cregion == "IgG2B" | cregion == "IgG2C" | cregion == "TRA" | cregion == "TRB")
    # filter(expr, cell_type %in% c("bj fibroblast", "hesc")) - changing to simpler negation......SHIFTING BACK IN CASE THIS IS CAUSING SHINY SERVER CODE TO FAIL...
    ## instead of exact matches to kappa or lambda, use grepl...https://stackoverflow.com/questions/22850026/filter-rows-which-contain-a-certain-string
    # x <- x %>% filter(!cregion %in% c("Kappa", "Lambda"))
    x <- x %>% dplyr::filter(!grepl('Kappa|Lambda', cregion))
    # x <- x %>% filter(cregion == "IgH" | cregion == "IgA" | cregion == "IgD" | cregion == "IgE" | cregion == "IgG" | cregion == "IgM" | cregion == "IgA1" | cregion == "IgA2" | cregion == "IgG1" | cregion == "IgG2" | cregion == "IgG3" | cregion == "IgG4" | cregion == "IgG2A" | cregion == "IgG2B" | cregion == "IgG2C" | cregion == "TRA" | cregion == "TRB")
  }
  ## this will remove all redundant sequences with same vgf/gene & cdr3 motif...note we count above so okay to collapse here!!
  if (filter_after_counting) {
    x <- x %>%
      group_by(cdr3_aa_imgt,vgf_jgene) %>%
      summarize_all(first) %>%
      rename(ncountfull = ncount) %>%
      ungroup() %>%
      add_count(vgf_jgene,cdr3length_imgt) %>%
      rename(ncount = n) %>%
      relocate(ncount, .before = shm)
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
  ## also adding as.character to the function - doesn't seem necessary here
  # x$sequence_id <- as.character(x$sequence_id)
  # x$cdr3_aa_imgt <- as.character(x$cdr3_aa_imgt)
  ## need to always check and remove underscores from all names
  x$sequence_id <- gsub("\\_","\\-",x$sequence_id)
  return(x)
}


##############################################################################################################################
##############################################################################################################################

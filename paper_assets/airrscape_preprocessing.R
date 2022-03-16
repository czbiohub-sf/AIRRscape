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

## login to iReceptor, search metadata for SubjectID = 7450 & Diagnosis = COVID-19
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
# fastq-dump --split-files --origfmt --gzip SRR8647465

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
# fastq-dump --split-files --origfmt --gzip SRR17417727
# fastq-dump --split-files --origfmt --gzip SRR17417728

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

OAS.clusters.all <- rbind(OAS.clusters.p148,OAS.clusters.p172,OAS.clusters.p194,OAS.clusters.p199,OAS.clusters.p203,OAS.clusters.p208,OAS.clusters.p232,OAS.clusters.p237,OAS.clusters.p238,OAS.clusters.p240,OAS.clusters.p249,OAS.clusters.p252,OAS.clusters.p255,OAS.clusters.p265,OAS.clusters.p275,OAS.clusters.p276,OAS.clusters.p287,OAS.clusters.p289,OAS.clusters.p299,OAS.clusters.p301,OAS.clusters.p307,OAS.clusters.p311,OAS.clusters.p320,OAS.clusters.p346,OAS.clusters.p376,OAS.clusters.p391,OAS.clusters.p422,OAS.clusters.p444,OAS.clusters.p455,OAS.clusters.p479,OAS.clusters.p481,OAS.clusters.p489,OAS.clusters.p500,OAS.clusters.p514,OAS.clusters.p515,OAS.clusters.p517,OAS.clusters.p520,OAS.clusters.p524,OAS.clusters.p529,OAS.clusters.p543,OAS.clusters.p551,OAS.clusters.p555,OAS.clusters.p558,OAS.clusters.p563,OAS.clusters.p569)

OAS.clusters.all$ANARCI_numbering <- NULL
OAS.clusters.all$ANARCI_status <- NULL
write.table(OAS.clusters.all, "OAS_sept21_germ-pass2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#den.bulk.OAS <- OAS.clusters.all
rm(list=ls(pattern="OAS."))
## this file "OAS_sept21_germ-pass2.tsv" is already available in the /paper_assets/intermediate_files folder
### further processing can be found in airrscape_processing.R script


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

### CODE FOR IMPORTING HIV-MT1214 DATA

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
OAS.clusters.mt1214all <- OAS.clusters.mt1214all %>% unite(sequence_id, dataset, obs, sep = "_", remove = TRUE, na.rm = TRUE)
hiv.bulk.mt1214t <- OAS.clusters.mt1214all
#write.table(OAS.clusters.mt1214all, "MT1214downloaded.tab", sep = "\t", row.names = FALSE, quote = FALSE)
rm(list=ls(pattern="OAS."))
## this file "MT1214downloaded.tab" is already available in the /paper_assets/intermediate_files folder
### further processing can be found in airrscape_processing.R script

##############################################################################################################################
##############################################################################################################################

### CODE FOR INTIAL IMPORTING SETLIFF 2018 HIV DATA FROM OAS AND COMBINING 100 RAW DATASETS INTO 12 DONOR TIMEPOINTS:
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
### further processing can be found in airrscape_processing.R script

##############################################################################################################################
##############################################################################################################################

### CODE FOR IMPORTING ANTIBODY SEQUENCES FROM CoV-AbDab
## first download full database as CSV (e.g. http://opig.stats.ox.ac.uk/webapps/covabdab/static/downloads/CoV-AbDab_xxyyzz.csv)
### version 1.0 of AIRRscape usees CoV-AbDab from CoV-AbDab_090721 database snapshot

## this list includes all CoV antibodies, and many non-human
## focus is protein sequence, and some entries are not complete (i.e. not full length sequence)
## list is first trimmed to antibodies with full length sequence, binding to SARS-CoV-2, and derived from humans or humanized mice

## for AIRRscape visualization of SHM, we searched Genbank using tblastn tosearch for 100% sequence match to nucleotide sequences 
## these were then checked for artificially high SHM due to codon optimization (e.g. Seydoux et al., 2020 (https://www.biorxiv.org/content/10.1101/2020.05.12.091298v1))


## these are only in SRA, so need full processing through Immcantation
## d13-enriched
## https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR17417727
## d13-stimulated
## https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR17417728
# fastq-dump --split-files --origfmt --gzip SRR17417727
# fastq-dump --split-files --origfmt --gzip SRR17417728

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

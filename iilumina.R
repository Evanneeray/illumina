rm(list=ls())

BiocManager::install('limma')
library(limma)
BiocManager::install('GEOquery')
library(GEOquery)
library(tidyverse)
library(readxl)

###### Download and processing of the GSE files ########
#list of all the GSEs we decided to include that use Illumina beadchip platform
GSE_beadchip <- read_excel('GEO summary table.xlsx') %>% 
  filter(str_detect(`Platform type`,'beadchip')) %>% 
  pull(GSE)

#downloading GSEs, un-taring them and chacking if all of them have the suppl files available
no_files <- c()
for (gse in GSE_beadchip){
  getGEOSuppFiles(gse)
  downloaded_GSEs <- list.files()
  if (gse %in% downloaded_GSEs)
    {untar(paste0(gse,'/',gse,'_RAW.tar'), exdir = paste0('./',gse))}
  else
    {no_files <- c(append(no_files,gse))}}
print(no_files)

#Un-zipping individual gz files to get text files
all_files <- c()
to_remove <- c()
for (gse in GSE_beadchip){
  all_files <- c(append(all_files,list.files(paste0('./',gse),full.names = TRUE)))
  to_remove <- c(append(to_remove,paste0('./',gse,'/', gse,'_RAW.tar')))}
files_to_unzip <- setdiff(all_files,to_remove)

for (file in files_to_unzip){
  gunzip(file, destname = (str_sub(file, 1, -4)),overwrite = FALSE,remove =FALSE)}
#####################


###### extracting raw data from individual GSMs #########

###### GSE101710 ########
GSE101710_rawData <- read.ilmn("GSE101710/GSE101710_Expression.raw.txt",
                               expr='SAMPLE',
                               sep='\t',
                               annotation='ID_REF')
nec(GSE101710_rawData) -> GSE101710_normalised #in-chip normalisation
log2(GSE101710_normalised$E) -> GSE101710_log2 #log2 transformation of expression data
data.frame(GSE101710_rawData$genes,GSE101710_log2) -> GSE101710_fixed #adding Illumina probe IDs
colnames(GSE101710_fixed)[1] <- 'Probe_Id' # changing Illumina probe ID column name for merging
getGEO('GSE101710') -> GSE101710_matrix # Loading matrix files with GSM sample names.
GSE101710_GSMs <- colnames(GSE101710_matrix$GSE101710_series_matrix.txt.gz)
colnames(GSE101710_fixed)[2:80] <- GSE101710_GSMs # add colnames (GSMs) from the matrix file

read.csv('GSE101710/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt',sep='\t',skip=8) -> R1
data.frame('Probe_Id'=R1$Probe_Id,
           'Unigene_ID'=R1$Unigene_ID,
           'Entrez_ID'=R1$Entrez_Gene_ID)->R1r
#####################

###### GSE65218 ########
ds<-read.csv('GSE65218/GSE65218_non-normalized.txt',sep='\t')
colnames(ds)
grep('TargetID|Detection',colnames(ds))->no
setdiff(c(1:303),no)->nos
colnames(ds)[nos]<-paste('SIGNAL',colnames(ds)[nos],sep=' ')
colnames(ds)[1] <- 'ID_REF'
write.csv(ds,'GSE65218/GSE65218_non-normalized_fixed.txt',quote = F,row.names = F)

GSE65218_rawData <- read.ilmn("GSE65218/GSE65218_non-normalized_fixed.txt",
                               expr='SIGNAL',
                               sep=',',
                               annotation='ID_REF')
nec(GSE65218_rawData) -> GSE65218_normalised 
log2(GSE65218_normalised$E) -> GSE65218_log2
data.frame(GSE65218_rawData$genes,GSE65218_log2) -> GSE65218_fixed
colnames(GSE65218_fixed)[1] <- 'Probe_Id'
getGEO('GSE65218') -> GSE65218_matrix
GSE65218_GSMs <- colnames(GSE65218_matrix$GSE65218_series_matrix.txt.gz)
colnames(GSE65218_fixed)[2:152] <- GSE65218_GSMs

#there is no platform file available and the probe IDs are very strange
###############################

###### GSE59714 ########
#only platform files provided, no raw expression data file available

read.csv('GSE59714/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt',sep='\t',skip=8) -> R1
data.frame('Probe_Id'=R1$Probe_Id,
           'Unigene_ID'=R1$Unigene_ID,
           'Entrez_ID'=R1$Entrez_Gene_ID)->R1r
###############################

###### GSE59654 ########
ds<-read.csv('GSE59654/GSE59654_PBMC.raw.corrected.txt',sep='\t')
colnames(ds)
grep('ID_REF|Detection',colnames(ds))->no
setdiff(c(1:313),no)->nos
colnames(ds)[nos]<-paste('SIGNAL',colnames(ds)[nos],sep='')
colnames(ds)[1] <- 'ID_REF'
write.csv(ds,'GSE59654/GSE59654_PBMC.raw.corrected_fixed.txt',quote = F,row.names = F)

GSE59654_rawData <- read.ilmn("GSE59654/GSE59654_PBMC.raw.corrected_fixed.txt",
                               expr='SIGNAL',
                               sep=',',
                               annotation='ID_REF')
nec(GSE59654_rawData) -> GSE59654_normalised
#not sure if the in-chip normalisation is required if the file was saved as "corrected"
log2(GSE59654_normalised$E) -> GSE59654_log2
data.frame(GSE59654_rawData$genes,GSE59654_log2) -> GSE59654_fixed 
colnames(GSE59654_fixed)[1] <- 'Probe_Id'
getGEO('GSE59654') -> GSE59654_matrix
GSE59654_GSMs <- colnames(GSE59654_matrix$GSE59654_series_matrix.txt.gz)
colnames(GSE59654_fixed)[2:157] <- GSE59654_GSMs 

read.csv('GSE59654/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt',sep='\t',skip=8) -> R1
data.frame('Probe_Id'=R1$Probe_Id,
           'Unigene_ID'=R1$Unigene_ID,
           'Entrez_ID'=R1$Entrez_Gene_ID)->R1r
###############################

###### GSE94497 ########
ds<-read.csv('GSE94497/GSE94497_non-normalized.txt',sep='\t')
ds<- select(ds, -'Array_Address_Id')
colnames(ds)
colnames(ds)<-paste('SIGNAL',colnames(ds),sep='')
colnames(ds)[1] <- 'ID_REF'
write.csv(ds,'GSE94497/GSE94497_non-normalized_fixed.txt',quote = F,row.names = F)

GSE94497_rawData <- read.ilmn("GSE94497/GSE94497_non-normalized_fixed.txt",
                              expr='SIGNAL',
                              sep=',',
                              annotation='ID_REF')
#can't be normalised because it's missing the p values
log2(GSE94497_rawData$E) -> GSE94497_log2
data.frame(GSE94497_rawData$genes,GSE94497_log2) -> GSE94497_fixed 
colnames(GSE94497_fixed)[1] <- 'Probe_Id' 
getGEO('GSE94497') -> GSE94497_matrix
GSE94497_GSMs <- colnames(GSE94497_matrix$GSE94497_series_matrix.txt.gz)
colnames(GSE94497_fixed)[2:55] <- GSE94497_GSMs

read.csv('GSE94497/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt',sep='\t',skip=8) -> R1
data.frame('Probe_Id'=R1$Probe_Id,
           'Unigene_ID'=R1$Unigene_ID,
           'Entrez_ID'=R1$Entrez_Gene_ID)->R1r
###############################

###### GSE94496 ########
ds<-read.csv('GSE94496/GSE94496_non-normalized.txt',sep='\t')
ds<- select(ds, -'Array_Address_Id')
colnames(ds)
colnames(ds)<-paste('SIGNAL',colnames(ds),sep='')
colnames(ds)[1] <- 'ID_REF'
write.csv(ds,'GSE94496/GSE94496_non-normalized_fixed.txt',quote = F,row.names = F)

GSE94496_rawData <- read.ilmn("GSE94496/GSE94496_non-normalized_fixed.txt",
                              expr='SIGNAL',
                              sep=',',
                              annotation='ID_REF')
#can't be normalised because it's missing the p values
log2(GSE94496_rawData$E) -> GSE94496_log2
data.frame(GSE94496_rawData$genes,GSE94496_log2) -> GSE94496_fixed 
colnames(GSE94496_fixed)[1] <- 'Probe_Id'
getGEO('GSE94496') -> GSE94496_matrix
GSE94496_GSMs <- colnames(GSE94496_matrix$GSE94496_series_matrix.txt.gz)
colnames(GSE94496_fixed)[2:209] <- GSE94496_GSMs
###############################


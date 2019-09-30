#### Diffex RNA iPS

rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

## ProteinCoding
inputCounts <- "~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/data/merged_pt6_hg38_mapping/mergedMapping_PC_pt6_hg38.s2.multi.Gencode27.Exon.sjdb.txt"
#### 1. read RNA data
rna.input <- read.delim(inputCounts,skip=1)
head(rna.input)
rna.abund <- rna.input[,7:length(rna.input[1,])]

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2) 
#install.packages("devtools")
library(devtools)
#install_github(repo = "perllb/deseqAbstraction",username = "perllb")
library(deseqAbstraction)

## FeatureCount raw file
path <- inputCounts
## Get the sampleFile name header
header <- read.delim(path,nrows = 1,skip=1,header=T)
## Define sample names, by removing pre- and suffix 
samples <- sub('.*\\_', '', names(header))

#### 1.3 Define metadata
##### Note that 'condition' and 'samples' are needed here. Conditions define the design of DESeq object, defines what can be tested for and how the data is normalized.

colDat <- read.delim(file = "data/merged_pt6_hg38_mapping/mergedMapping_pt6_hg38_colData.txt",sep="\t")
colDat

cbind(colDat,samples[-c(1:6)])

#######

#### ips only , remove h9#####
colDat.ips <- colDat[colDat$line!="HS_h9",]
colDat.ips$condition <- colDat.ips$species

rna.ips <- rna.input[,c(1:6,grep(paste(colDat.ips$samples,collapse = "|"),colnames(rna.input)))]
head(rna.ips)
header <- read.delim(file = path,header = F,nrows = 1)
path.ips <- "data/merged_pt6_hg38_mapping/merged_rna.raw.ips.txt"
write.table(header,file = path.ips,quote = F,row.names = F,col.names = F,sep="\t")
write.table(rna.ips,file = path.ips,quote = F,row.names = F,col.names = T,sep="\t",append = T)



dabs.ips <- deseqAbs$new("humanChimp",colData=colDat.ips,filename=path.ips,design = formula(~batch+condition))
head(dabs.ips$rawCounts)
dabs.ips$colData

dabs.ips$makeDiffex()

dabs.ips$test$Default

write.table(x = dabs.ips$test$Default,file = "data/DESeq_Human.vs.Chimp_allIPS.txt",quote = F,row.names = T,col.names = T,sep = "\t")

# get up
sign <- getSign(x = dabs.ips$test$Default,p = 0.05,l = 0)
sign.up <- sign$up
sign.down <- sign$down

write.table(sign.up,file = "data/DESeq_higherHuman_Human.vs.Chimp_allIPS.txt",quote = F,row.names = T,col.names = T,sep = "\t")
write.table(sign.down,file = "data/DESeq_higherChimp_Human.vs.Chimp_allIPS.txt",quote = F,row.names = T,col.names = T,sep = "\t")


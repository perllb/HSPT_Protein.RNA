### fix count table
# Chimp samples to pantro6
# human samples to hg38

# Unique

rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/KRABZNF-manuscript/Data/bulkRNA_allQuant/")

### 1. Read HG38 Multi RNAseq data ####
# Genome: hg38
# Mapping: unique

## Read human 
## FeatureCount raw file
path <- "~/Dropbox (MN)/Per/PhD/Projects/KRABZNF-manuscript/Data/bulkRNA_allQuant/hg38.s2.unique.Gencode27.Exon.sjdb.txt"

## Get the sampleFile name header
header <- read.delim(path,nrows = 1,skip=1,header=T)

## Define sample names, by removing pre- and suffix 
colNames <- gsub(pattern ="X.projects.fs1.medpvb.backup.projects.ChimpHuman.RNAseq.Aligned_hg38_STAR_unique.sjdb.genc.v27.hg38.unique.sjdb.genc.v27.",replacement = "",gsub(pattern = "Aligned.out.bam",replacement = "",x = names(header)))

#### 1.3 Define metadata
##### Note that 'condition' and 'samples' are needed here. Conditions define the design of DESeq object, defines what can be tested for and how the data is normalized.

## fix colData
line <- c(rep(c("PT_C6","HS_H6","HS_h9"),8),c("HS_H48","HS_H48","PT_CPT","PT_CPT"))
species <- c(rep(c("chimp","human","human"),8),c("human","human","chimp","chimp"))
time <- c(rep(c(rep("d13",3),rep("d14",3),rep("d15",3),rep("d16",3)),2),rep("d14",4))
cond <- paste(line,time,sep = "-")
batch <- c(rep(c("batch_1","batch_2"),each=12),rep("batch_3",4))
colDat <- data.frame(species=species,time=time,condition=species,line=line,samples=colNames[-c(1:6)],batch=batch)

human.input <- read.delim(path,skip=1)
head(human.input)
names(human.input) <- colNames
head(colDat)

##### READ pantro6 unique ###########

## FeatureCount raw file
path <- "~/Dropbox (MN)/Per/PhD/Projects/KRABZNF-manuscript/Data/bulkRNA_allQuant/pantro6.s2.unique.Gencode27.Exon.sjdb.txt"

## Get the sampleFile name header
header <- read.delim(path,nrows = 1,skip=1,header=T)
## Define sample names, by removing pre- and suffix 
colNames <- gsub(pattern ="X.projects.fs1.medpvb.backup.projects.ChimpHuman.RNAseq.Aligned_pantro6_STAR_unique.sjdb.pantro6.unique.sjdb.genc27.lifted.",replacement = "",gsub(pattern = "Aligned.out.bam",replacement = "",x = names(header)))

#### 1.3 Define metadata
##### Note that 'condition' and 'samples' are needed here. Conditions define the design of DESeq object, defines what can be tested for and how the data is normalized.

## fix colData
line <- c(rep(c("PT_C6","HS_H6","HS_h9"),8),c("HS_H48","HS_H48","PT_CPT","PT_CPT"))
species <- c(rep(c("chimp","human","human"),8),c("human","human","chimp","chimp"))
time <- c(rep(c(rep("d13",3),rep("d14",3),rep("d15",3),rep("d16",3)),2),rep("d14",4))
cond <- paste(line,time,sep = "-")
batch <- c(rep(c("batch_1","batch_2"),each=12),rep("batch_3",4))
colDat <- data.frame(species=species,time=time,condition=species,line=line,samples=colNames[-c(1:6)],batch=batch)

chimp.input <- read.delim(path,skip=1)
names(chimp.input) <- colNames
head(chimp.input)
head(colDat)

##### 3. get chimp samples from chimp input
head(chimp.input)

chimp.coldat <- colDat[colDat$species=="chimp",]
chimp.samples <- chimp.input[,names(chimp.input) %in% chimp.coldat$samples]
rownames(chimp.samples) <- chimp.input$Geneid
colnames(chimp.samples) <- paste("pt6_multi",chimp.coldat$line,chimp.coldat$samples,sep = "_")
head(chimp.samples)

##### 4. get human samples from human input
head(human.input)

human.coldat <- colDat[colDat$species=="human",]
human.samples <- human.input[,names(human.input) %in% human.coldat$samples]
rownames(human.samples) <- human.input$Geneid
colnames(human.samples) <- paste("hg38_multi",human.coldat$line,human.coldat$samples,sep = "_")
head(human.samples)

##### 5. Merge human and chimp
merged.hs.pt <- merge(x = human.samples,y = chimp.samples,by=0)
rownames(merged.hs.pt) <- merged.hs.pt$Row.names
merged.hs.pt <- merged.hs.pt[,-1]

dim(merged.hs.pt)
dim(human.samples)
dim(chimp.samples)

head(human.input)
head(merged.hs.pt)

merged.hs.pt.info <- merge(x = human.input[,1:6],y = merged.hs.pt,by.x=1,by.y=0)
rownames(merged.hs.pt.info) <- merged.hs.pt.info$Geneid
head(merged.hs.pt.info)

### 6. PC GENES ##########

# Get only protein coding genes 
#pcgenes <- read.csv(file = "~/Documents/bioinformatics/genomicData/hg38/gencode/gencode.protein_coding.v27.ID.txt",header=F)
pcgenes <- read.csv(file = "~/Documents/bioinformatics/genomicData/hg38/gencode/gencode.protein_coding.v27.ID.txt",header=F)
head(pcgenes)
pcgenes <- pcgenes[!duplicated(pcgenes$V1),]

merged.PC <- merge(pcgenes,merged.hs.pt.info,by.x=1,by.y=1)
head(merged.PC)
dim(merged.PC)
length(pcgenes)

merged.PC[grep("MT-",merged.PC$x),]

setdiff(rownames(human.samples),rownames(chimp.samples))

#### 7. Print table #####
## PC

header <- read.delim(file = path,header = F,nrows = 1)
path.merged.PC <- "merge/Merged_Unique_Hg38.PanTro6_ProteinCoding.txt"
head(merged.PC)
write.table(header,file = path.merged.PC,quote = F,row.names = F,col.names = F,sep="\t")
write.table(merged.PC,file = path.merged.PC,quote = F,row.names = F,col.names = T,sep="\t",append = T)

## All genes
header <- read.delim(file = path,header = F,nrows = 1)
path.merged <- "merge/Merged_Unique_Hg38.PanTro6_allGenes.txt"
write.table(header,file = path.merged,quote = F,row.names = F,col.names = F,sep="\t")
write.table(merged.hs.pt.info,file = path.merged,quote = F,row.names = F,col.names = T,sep="\t",append = T)


#### 8. print coldata ####
chimp.coldat <- colDat[colDat$species=="chimp",]
chimp.coldat$mapGenome <- "panTro6"
chimp.coldat$mapping <- "uniqueMap"
human.coldat <- colDat[colDat$species=="human",]
human.coldat$mapGenome <- "hg38"
human.coldat$mapping <- "uniqueMap"

write.table(x = rbind(human.coldat,chimp.coldat),file = "Merge/Merged_Unique_Hg38.PanTro6_ColData.txt",quote = F,row.names = F,col.names = T,sep = "\t")


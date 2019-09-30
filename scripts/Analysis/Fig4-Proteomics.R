# Fig 4 Proteomics ######

rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

# 1. Read proteome abundance ####
# Matched by geneName _ all _ normalized and imputed
proteome.input <- read.delim(file = "data/MatchedByGeneNames_HighOnly_normIMP_ModNames.txt",header=T)
head(proteome.input)
names(proteome.input)[1:3] <- c("H48","H48.1","H48.2")

prot.norm <- proteome.input[,1:12]
rownames(prot.norm) <- make.names(proteome.input$X.31,unique = T)

# Fix MT genes ID
prot.det <- as.character(rownames(prot.norm))
prot.det[grep("MT",prot.det)]
prot.det.mod <- gsub(pattern = "^MT\\.",replacement = "MT-",x = prot.det)
prot.det.mod <- gsub(pattern = "CORO7.PAM16",replacement = "CORO7-PAM16",x = prot.det.mod)

rownames(prot.norm) <- prot.det.mod
head(prot.norm)
dim(prot.norm)

### Testing according to limma (http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html)

# load limma and qvalue
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("qvalue")
library(limma)
library(qvalue)
library(wesanderson)
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

citation("ebayes")
# plotting norm data
boxplot(prot.norm)
head(prot.norm)

# Preparing data for limma
human.id <- colnames(prot.norm)[1:6]
chimp.id <- colnames(prot.norm)[7:12]

design <- model.matrix(~factor(c(2,2,2,2,2,2,1,1,1,1,1,1)))
colnames(design) <- c("Intercept","Diff")

fit <- lmFit(prot.norm[, c(human.id,chimp.id)], design)
fit$coefficients

res.eb <- eb.fit(prot.norm[, c(human.id,chimp.id)], design)
head(res.eb)
res.eb["PBX1",]
####

# colData prot
cnames.p <- colnames(prot.norm)
colDat.prot <- data.frame(
  colnames = cnames.p,
  species = ifelse(grepl("H",cnames.p),yes = "human",no = "chimp"),
  line = gsub("\\..*","",cnames.p)
)
colDat.prot

########


# 2. read RNA data ####
setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

## ProteinCoding
inputCounts <- "~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/data/merged_pt6_hg38_mapping/mergedMapping_PC_pt6_hg38.s2.multi.Gencode27.Exon.sjdb.txt"

rna.input <- read.delim(inputCounts,skip=1)
head(rna.input)
rna.abund <- rna.input[,7:length(rna.input[1,])]

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2) 
#install.packages("devtools")
library(devtools)
install_github(repo = "perllb/deseqAbstraction",username = "perllb")
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

# ips d14 only 
colDat.ips.d14 <- colDat[colDat$line!="HS_h9" & colDat$time=="d14",]
colDat.ips.d14

rna.ips.d14 <- rna.input[,c(1:6,grep(paste(colDat.ips.d14$samples,collapse = "|"),colnames(rna.input)))]
head(rna.ips.d14)
header <- read.delim(file = path,header = F,nrows = 1)
path.ips.d14 <- "data/merged_pt6_hg38_mapping/rna.raw.ips.d14.txt"
write.table(header,file = path.ips.d14,quote = F,row.names = F,col.names = F,sep="\t")
write.table(rna.ips.d14,file = path.ips.d14,quote = F,row.names = F,col.names = T,sep="\t",append = T)

dabs.ips.d14 <- deseqAbs$new("humanChimp",colData=colDat.ips.d14,filename=path.ips.d14,design = formula(~batch+condition))
head(dabs.ips.d14$rawCounts)
dabs.ips.d14$makeVST(blind=F)

ips.d14.vst.br <- limma::removeBatchEffect(assay(dabs.ips.d14$VST),batch = colDat.ips.d14$batch)
head(ips.d14.vst.br)

#######


### 3. Calc average on cell lines ####
colnames(prot.norm) <- colDat.prot$colnames
df <- prot.norm
prot.mean <- sapply( levels(colDat.prot$line), function(lvl) rowMeans( df[,colDat.prot$line == lvl] ) )
head(prot.mean)

# rna deseq
colDat.d14 <- colDat[colDat$time=="d14" & colDat$line!="h9",]
df <- ips.d14.vst.br
tail(df)
rna.mean.deseq <- sapply( levels(dabs.ips.d14$colData$line), function(lvl) rowMeans( df[,dabs.ips.d14$colData$line == lvl] ) )
head(rna.mean.deseq)
dim(rna.mean.deseq)

# Merge means 
colnames(rna.mean.deseq) <- paste("RNA",colnames(rna.mean.deseq),sep = "_")
colnames(prot.mean) <- paste("Prot",colnames(prot.mean),sep = "_")
head(rna.mean.deseq)
head(prot.mean)

merge.mean <- merge(rna.mean.deseq,prot.mean,by=0)
dim(merge.mean)
rownames(merge.mean) <- merge.mean$Row.names
merge.mean <- merge.mean[,-1]
head(merge.mean)

#####

# 4. Calc average species (d14 iPS lines only) #####
# protein
head(prot.norm)
colnames(prot.norm) <- colDat.prot$colnames
df <- prot.norm
prot.mean.species <- sapply( levels(colDat.prot$species), function(lvl) rowMeans( df[,colDat.prot$species == lvl] ) )

# rna deseq
colDat.d14 <- colDat[colDat$time=="d14" & colDat$line!="h9",]
df <- ips.d14.vst.br
tail(df)
rna.mean.species <- sapply( levels(dabs.ips.d14$colData$species), function(lvl) rowMeans( df[,dabs.ips.d14$colData$species == lvl] ) )
head(rna.mean.species)
dim(rna.mean.species)

### Merge means #
colnames(rna.mean.species) <- paste("RNA",colnames(rna.mean.species),sep = "_")
colnames(prot.mean.species) <- paste("Prot",colnames(prot.mean.species),sep = "_")
head(rna.mean.species)
head(prot.mean.species)

merge.mean.species <- merge(rna.mean.species,prot.mean.species,by=0)
rownames(merge.mean.species) <- merge.mean.species$Row.names
merge.mean.species <- merge.mean.species[,-1]
head(merge.mean.species)
dim(merge.mean.species)
########

plotDir <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig4-Proteomics/"


# Fig4B: Ma plot RNA - highlight prot det RNA
dabs.ips.d14$makeDiffex()

plotDat <- data.frame(baseMean = log2(dabs.ips.d14$test$Default$baseMean),log2FC=dabs.ips.d14$test$Default$log2FoldChange,protDet=factor(ifelse(rownames(dabs.ips.d14$test$Default) %in% rownames(prot.norm),yes = 1,no = 0)))
head(plotDat)

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig4-Proteomics/Fig4B_MaPlot_RNA_protDetMarked.pdf"
pdf(file = filename,width = 6,height = 6)
p <- ggplot(data = plotDat,mapping = aes(x = baseMean,y = log2FC,color=protDet)) +
  geom_point(alpha=.4,cex=.4) +
  geom_point(data = plotDat[plotDat$protDet==1,],alpha=.4) +
  theme_classic() +
  scale_color_manual(values = c("black",wes_palette("Darjeeling1")[4])) +
  ggtitle(label = "RNA changes - genes detected on protein level")
print(p)
dev.off()


#####

## Fig4C: Protein markers


# Define marker groups
Pluripotency <- data.frame(Genes=c("NANOG","POU5F1","KLF4","MYC","LIN28A"),Group="Pluripotency")
dForebrain <- data.frame(Genes=c("FOXG1","OTX2","PAX6","MEIS2"),Group="dorsal Forebrain")
vForebrain <- data.frame(Genes=c("PITX2","BARHL1","NKX2-1"),Group="ventral Forebrain")
Midbrain <-  data.frame(Genes=c("OTX1","EN1","LMX1A","CORIN","FOXA2","PAX8"),Group="Midbrain / Hindbrain")
Hindbrain <-  data.frame(Genes=c("HOXA2"),Group="Midbrain / Hindbrain")
Progenitors <- data.frame(Genes=c("SOX2","TUBB","DCX","TBR1","EOMES","CUX2","OLIG2","GFAP","S100G","EMX1"),Group="Fate specification")

fix.mat.marker = function(markers) {
  
  mat <- merge(x = markers,y = prot.norm,by.x=1,by.y=0,sort=F)
  rownames(mat) <- mat$x
  mat <- mat[,-1]
  
  notIn <- setdiff(markers,rownames(mat))
  for(gene in notIn){
    
    mat[gene,] <- rep(min(prot.norm),ncol(prot.norm))  
  }
  return(mat)
  
}

plurip.dat <- fix.mat.marker(as.character(Pluripotency$Genes))
prog.dat <- fix.mat.marker(as.character(Progenitors$Genes))
dforebrain.dat <- fix.mat.marker(as.character(dForebrain$Genes))
vforebrain.dat <- fix.mat.marker(as.character(vForebrain$Genes))
midbrain.dat <- fix.mat.marker(as.character(Midbrain$Genes))
hindbrain.dat   <- fix.mat.marker(as.character(Hindbrain$Genes))

mark.dat <- rbind(plurip.dat,dforebrain.dat,vforebrain.dat,midbrain.dat,hindbrain.dat,prog.dat)
allMarkers <- rbind(Pluripotency,dForebrain,vForebrain,Midbrain,Hindbrain,Progenitors)
rownames(allMarkers) <- allMarkers$Genes

colData <- colDat.prot

# species
a1col <- c( "darkolivegreen", 	"rosybrown")
# line
a2col <- rev(c("darkolivegreen3","darkolivegreen1","burlywood4","burlywood2"))
# time
a3col <- c("grey100","grey75","grey50","grey25")

mycolors1 <- a1col
names(mycolors1) <- c("human","chimp")

mycolors2 <- a2col
names(mycolors2) <- c("C6","CPT","H6","H48")

mycolors3 <- a3col
names(mycolors3) <- c("d13","d14","d15","d16")

#order  markdat by molecule -> species -> line
ordr <- order(colData$line)

#View(mark.dat)
library(RColorBrewer)
col_markersGr <- colorRampPalette(brewer.pal(length(unique(allMarkers$Group)),"Set2"))(length(unique(allMarkers$Group)))
names(col_markersGr) <- unique(allMarkers$Group)

myColors <- list(species = mycolors1,line = mycolors2,Group=col_markersGr)
myColors
annot <- colData[ordr,c(2,3)]
annot
colnames(mark.dat) <- make.names(colnames(mark.dat),unique = T)
rownames(annot) <- colnames(mark.dat)

annot.row <- data.frame(allMarkers$Group)
rownames(annot.row) <- allMarkers$Genes
colnames(annot.row) <- "Group"

library(pheatmap)

gaps <- which(duplicated(annot.row$Group) == F)-1

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig4-Proteomics/Fig4C_Markers_Protein_GroupedRows.pdf",width = 5,height = 8)
pheatmap(mark.dat,annotation_col = annot,annotation_row = annot.row,cluster_rows = F,
         annotation_colors = myColors,color = rev(colorRampPalette(brewer.pal(10,"RdBu"))(200)),
         cluster_cols = F,cutree_cols = 2,fontsize_col = 5,gaps_row = gaps,gaps_col = 6)
dev.off()



#### Fig 4D - RNA between cell lines, same species #####
# corr plot function
plotCor <- function(x,y,species,lineX,lineY,molecule){
  
  fit <- lm(y ~ x) 
  cor <- cor.test(x = x, y = y,method = "pearson",conf.level = 0.95)
  
  plot(x = x,y = y,pch=16,cex=.2,ylab=paste(lineY," abundance",sep=""),xlab = paste(lineX," abundance",sep=""))
  mtext(text = paste(molecule,": ",species,"\ncell lines (abundance)",sep=""),line = 1,cex = 1)
  abline(fit,lty=2,lwd=1.5,col="red")
  grid(nx = 10,ny = 10)
  points(x,y,cex=.3,pch=16)
  legend("topleft",bg = "white",
         legend=paste("Corr (pearson): ",format(cor$estimate,digits = 3)),
         cex = .7)
  
  return(cor$estimate)
}

# RNA
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig4-Proteomics/Fig4D_RNA_cellLine_sameSpecies_Means.pdf",height = 7,width = 4)
par(mfrow=c(2,1))
x <- merge.mean$RNA_HS_H6
y <- merge.mean$RNA_HS_H48
plotCor(x,y,"Human","human H6","human H48","RNA")

# chimp RNA
x <- merge.mean$RNA_PT_C6
y <- merge.mean$RNA_PT_CPT
plotCor(x,y,"Chimp","chimp C6","chimp CPT","RNA")
dev.off()


#### Fig 4E - Protein between cell lines, same species #####

# Protein
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig4-Proteomics/Fig4E_Protein_sameSpecies_cellLineMeans.pdf",height = 7,width = 4)
par(mfrow=c(2,1))
x <- merge.mean$Prot_H6
y <- merge.mean$Prot_H48
plotCor(x,y,"Human","human H6","human H48","Protein")

# chimp Protein
x <- merge.mean$Prot_C6
y <- merge.mean$Prot_CPT
plotCor(x,y,"Chimp","chimp C6","chimp CPT","Protein")
dev.off()



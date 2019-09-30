## Fig2 - bulkRNAseq characterization #####

# 1. Prepare data ####

## Load and prepare RNAseq data

rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

## ProteinCoding
inputCounts <- "~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/data/merged_pt6_hg38_mapping/merge/Merged_MultiMapping_Hg38.PanTro6_ProteinCoding.txt"

### read RNA data
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

#### Define metadata
##### Note that 'condition' and 'samples' are needed here. Conditions define the design of DESeq object, defines what can be tested for and how the data is normalized.

colDat <- read.delim(file = "data/merged_pt6_hg38_mapping/merge/Merged_MultiMapping_Hg38.PanTro6_ColData.txt",sep="\t")
colDat

cbind(colDat,samples[-c(1:6)])

#### ips H6 and CS only 
colDat.ips <- colDat[colDat$line!="HS_h9" & colDat$line!="HS_H48" & colDat$line!="PT_CPT",]
colDat.ips$condition <- paste(colDat.ips$time,colDat.ips$species,sep = "_")

rna.ips <- rna.input[,c(1:6,grep(paste(colDat.ips$samples,collapse = "|"),colnames(rna.input)))]
head(rna.ips)
header <- read.delim(file = path,header = F,nrows = 1)
path.ips <- "data/merged_pt6_hg38_mapping/merge/rna.raw.ips.PC.C6.H6.txt"
write.table(header,file = path.ips ,quote = F,row.names = F,col.names = F,sep="\t")
write.table(rna.ips,file = path.ips,quote = F,row.names = F,col.names = T,sep="\t",append = T)

dabs.ips <- deseqAbs$new("humanChimp",colData=colDat.ips,filename=path.ips,design = formula(~batch+condition))
head(dabs.ips$rawCounts)
dabs.ips$makeVST(blind=F)

dabs.ips.vst.br <- limma::removeBatchEffect(assay(dabs.ips$VST),batch = colDat.ips$batch)

# all samples 
dabs.all <- deseqAbs$new("humanChimp",colData=colDat,filename=path,design = formula(~batch+condition))
dabs.all$getAverage()
head(dabs.all$FPKMMean$Mean) 
head(dabs.all$FPKM)


## 
pcgenes <- read.csv(file = "~/Documents/bioinformatics/genomicData/hg38/gencode/gencode.protein_coding.v27.ID.txt",header=F)
head(pcgenes)
pcgenes <- pcgenes[!duplicated(pcgenes$V1),]

rna.input.PC <- merge(pcgenes,rna.input,by.x=1,by.y=1)
head(rna.input.PC)
dim(rna.input.PC)
length(pcgenes)

### Get data and coldata for CHIMP SANDRA only 
colDat.cSandra <- colDat[colDat$line=='PT_C6',]
colDat.cSandra$condition <- colDat.cSandra$time

rna.cSandra <- rna.input.PC[,c(1:6,grep(paste(colDat.cSandra$samples,collapse = "|"),colnames(rna.input.PC)))]
head(rna.cSandra)
header <- read.delim(file = path,header = F,nrows = 1)
path.cSandra.PC <- "data/rna.raw.cSandra.PC.txt"
write.table(header,file = path.cSandra.PC,quote = F,row.names = F,col.names = F,sep="\t")
write.table(rna.cSandra,file = path.cSandra.PC,quote = F,row.names = F,col.names = T,sep="\t",append = T)

dabs.cSandra <- deseqAbs$new("humanChimp",colData=colDat.cSandra,filename=path.cSandra.PC,design = formula(~batch+condition))
head(dabs.cSandra$rawCounts)

### Get data and coldata for human ips6 only 
colDat.hips6 <- colDat[colDat$line=='HS_H6',]
colDat.hips6$condition <- colDat.hips6$time

rna.hips6 <- rna.input.PC[,c(1:6,grep(paste(colDat.hips6$samples,collapse = "|"),colnames(rna.input.PC)))]
head(rna.hips6)
header <- read.delim(file = path,header = F,nrows = 1)
path.hips6.PC <- "data/rna.raw.hips6.PC.txt"
write.table(header,file = path.hips6.PC,quote = F,row.names = F,col.names = F,sep="\t")
write.table(rna.hips6,file = path.hips6.PC,quote = F,row.names = F,col.names = T,sep="\t",append = T)

dabs.hips6 <- deseqAbs$new("humanChimp",colData=colDat.hips6,filename=path.hips6.PC,design = formula(~batch+condition))
head(dabs.hips6$rawCounts)

### makeVST(), makeDiffex() 

dabs.hips6$makeVST(blind = F)
dabs.cSandra$makeVST(blind = F)

dabs.hips6$makeDiffex()
dabs.cSandra$makeDiffex()

dabs.hips6$test$Default
dabs.cSandra$test$Default

### Remove batch effect 
library(limma)

hips6.vst.brm <- limma::removeBatchEffect(x = assay(dabs.hips6$VST),batch = dabs.hips6$colData$batch)
csan.vst.brm <- limma::removeBatchEffect(x = assay(dabs.cSandra$VST),batch = dabs.cSandra$colData$batch)

### Fuse cSandra and h6 VST 
vst.brm.hips6.cSan <- merge(hips6.vst.brm,csan.vst.brm,by=0)
head(vst.brm.hips6.cSan)

a1 <- rep(dabs.cSandra$colData$condition,2)
a2 <- c(as.character(dabs.hips6$colData$line) ,as.character(dabs.cSandra$colData$line))

rownames(vst.brm.hips6.cSan) <- vst.brm.hips6.cSan$Row.names
vst.brm.hips6.cSan <- vst.brm.hips6.cSan[,-1]

colnames(vst.brm.hips6.cSan) <- paste(a2,colnames(vst.brm.hips6.cSan),sep = "_")

### Cut off low expressed genes
min(vst.brm.hips6.cSan)
max(vst.brm.hips6.cSan)

vst.brm.hips6.cSan 

#######


# 2. Fig1C: plot Markers ######
### Fuunction to make heatmap of data, markers
plotMarkers <- function(data,markers,colData,myColors=NULL,heatCol = "RdBu") {
  
  
  #order  markdat by molecule -> species -> line
  ordr <- order(colData$line,colData$time)
  
  mark.dat <- data[rownames(data) %in% markers,]
  notin <- setdiff(markers,rownames(mark.dat))
  idx <- nrow(mark.dat)+1
  
  mark.dat <- mark.dat[,ordr]
  
  #View(mark.dat)
  
  annot <- colData[ordr,rev(c(1,4,2))]
  annot
  myColors
  colnames(mark.dat) <- make.names(colnames(mark.dat),unique = T)
  rownames(annot) <- colnames(mark.dat)
  
  
  library(pheatmap)
  library(RColorBrewer)
  
  if(!is.null(myColors)) {
    return(pheatmap(mark.dat,annotation_col = annot,annotation_colors = myColors,color = rev(colorRampPalette(brewer.pal(10,heatCol))(200)),cluster_cols = F,cutree_cols = 2,fontsize_col = 5))
  }else{
    return(pheatmap(mark.dat,annotation_col = annot,color = rev(colorRampPalette(brewer.pal(10,heatCol))(200)),cluster_cols = F,cutree_cols = 2,fontsize_col = 5))
  }
}


##### Plot all iPS #
# Set annotation colors 
# species
a1col <- c( "darkolivegreen", 	"rosybrown")
# line
a2col <- rev(c("darkolivegreen3","burlywood4"))
# time
a3col <- c("grey100","grey75","grey50","grey25")

mycolors1 <- a1col
names(mycolors1) <- c("human","chimp")

mycolors2 <- a2col
names(mycolors2) <- c("PT_C6","HS_H6")

mycolors3 <- a3col
names(mycolors3) <- c("d13","d14","d15","d16")

# Define marker groups
Pluripotency <- data.frame(Genes=c("NANOG","POU5F1","KLF4","MYC","LIN28A"),Group="Pluripotency")
dForebrain <- data.frame(Genes=c("FOXG1","OTX2","PAX6","MEIS2"),Group="dorsal Forebrain")
vForebrain <- data.frame(Genes=c("PITX2","BARHL1","NKX2-1"),Group="ventral Forebrain")
Midbrain <-  data.frame(Genes=c("OTX1","EN1","LMX1A","CORIN","FOXA2","PAX8"),Group="Midbrain / Hindbrain")
Hindbrain <-  data.frame(Genes=c("HOXA2"),Group="Midbrain / Hindbrain")
Progenitors <- data.frame(Genes=c("SOX2","TUBB","DCX","VIM","HES1","TBR1","EOMES","CUX2","OLIG2","GFAP","S100G","EMX1"),Group="Fate specification")

fix.mat.marker = function(markers) {
  
  mat <- merge(x = markers,y = dabs.ips.vst.br,by.x=1,by.y=0,sort=F)
  rownames(mat) <- mat$x
  mat <- mat[,-1]
  
  notIn <- setdiff(markers,rownames(mat))
  for(gene in notIn){
    
    mat[gene,] <- rep(min(dabs.ips.vst.br),ncol(dabs.ips.vst.br))  
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

#order  markdat by molecule -> species -> line
ordr <- order(colDat.ips$line,colDat.ips$time)

mark.dat <- mark.dat[,ordr]
head(mark.dat)

#View(mark.dat)
library(RColorBrewer)
col_markersGr <- colorRampPalette(brewer.pal(length(unique(allMarkers$Group)),"Set2"))(length(unique(allMarkers$Group)))
names(col_markersGr) <- unique(allMarkers$Group)

myColors <- list(species = mycolors1,line = mycolors2,time=mycolors3,Group=col_markersGr)
myColors
annot <- colDat.ips[ordr,rev(c(1,4,2))]
annot
colnames(mark.dat) <- make.names(colnames(mark.dat),unique = T)
rownames(annot) <- colnames(mark.dat)

annot.row <- data.frame(allMarkers$Group)
rownames(annot.row) <- allMarkers$Genes
colnames(annot.row) <- "Group"

library(pheatmap)

gaps <- which(duplicated(annot.row$Group) == F)-1

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2A_Markers_Neuro_RNA.only_GroupedRows.pdf",width = 5,height = 8)
pheatmap(mark.dat,annotation_col = annot,annotation_row = annot.row,cluster_rows = F,
         annotation_colors = myColors,color = rev(colorRampPalette(brewer.pal(10,"RdBu"))(200)),
         cluster_cols = F,cutree_cols = 2,fontsize_col = 5,gaps_row = gaps,gaps_col = 8,show_colnames = F)
dev.off()

########

# 3. Fig2B: Dynamic Markers ######
library(ggplot2)
library(reshape2)

dabs.ips$makeFPKM()
cond.b <- paste(colDat.ips$line,colDat.ips$time,colDat.ips$batch,sep = "_")
colnames(dabs.ips$FPKM) <- cond.b
melted <- melt(data = dabs.ips$FPKM)
head(melted)

times <- read.table(text = as.character(melted$Var2), sep = "_", as.is = TRUE)$V3
melted$time <- times

lines <- read.table(text = as.character(melted$Var2), sep = "_", as.is = TRUE)$V2
melted$line <- lines

batches <- read.table(text = as.character(melted$Var2), sep = "_", as.is = TRUE)$V5
melted$batch <- paste("b",batches,sep = "")

head(melted)

#### Plot function
plotGeom <- function(gene,es=T) {
  
  curr <- melted[melted$Var1==gene,]
  ifelse(test = es,yes = curr.2 <- curr[curr$line!="CPT" & curr$line!="H48",],no = curr.2 <- curr[curr$line!="CPT" & curr$line!="H48" & curr$line!="h9",])
  curr.3 <- curr[curr$line=="CPT" | curr$line=="H48",]
  
  ggplot(data = curr.2,mapping = aes(x = time,y = value,col = line)) +
    geom_point(aes(shape=batch),size=2) +
    geom_smooth(aes(group=line),se=F,method = "loess") +
    ggtitle(label = gene) +
    ylab(label = "FPKM") +
    scale_colour_manual(values=c("burlywood3","darkolivegreen4","darkolivegreen2")) +
    theme_bw()
  
}

## Multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## plot without ES
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2B_Marker.Dynamics_LIN28A.NEUROD1.SYP.pdf",width = 5,height = 12)
p1 <- plotGeom("LIN28A",es = F)
p2 <- plotGeom("NEUROD1",es = F)
#p3 <- plotGeom("ASCL1",es = F)
p4 <- plotGeom("SYP",es = F)
multiplot(p1,p2,p4,cols=1)
dev.off()

######

# 4. Fig2C: d16 vs d13 heatmaps

plotHeatGenes <- function(data,genes,a1,a2,n1,n2,sdCut=sdCut,z=z,col=col,rev=F) {
  
  genes.exp <- getGenes(data = data,genes = genes)
  rownames(genes.exp) <- make.names(genes.exp[,1],unique = T)
  genes.exp <- genes.exp[,-1]
  
  # Set annotation colors (9 colors)
  a2col <- c( "darkolivegreen", 	"rosybrown")
  #  olive green, gold,ligth sea, forest green, rosy brown, moccasin  , red ,   dark blue, 
  a1col <- c("grey100","grey75","grey50","grey25")
  # antique white, saddle brown , ligth steel blue,  slate brue,  navy , corn flower blue, teal  ,  lime,      khaki
  ## Change to RdYlBu (if RedBlue == F)
  heatCol <- col
  
  ifelse(test = rev,yes = heatScaleCol <- rev(colorRampPalette(brewer.pal(9,heatCol))(200)),no = heatScaleCol <- colorRampPalette(brewer.pal(9,heatCol))(200))
  
  ## get standard deviation of each gene
  sd.exp <- apply(genes.exp,1,sd)
  
  # data to plot
  plotData <- genes.exp[sd.exp > sdCut,]
  
  ## Show rownames? only if less that 41 genes are plotted
  rowShow <- T; if(nrow(plotData)>100) { rowShow <- F }
  print(paste("There are ",nrow(plotData)," genes in you gene-set with sd > ",sdCut,".",sep=""))
  
  ## scale by row?
  scale <- 'none' ; if(z) { scale <- 'row' }
  
  df <- data.frame(Var1 = factor(a1),Var2 = factor(a2))
  rownames(df) <- colnames(data)
  colnames(df) <- c(n1,n2)
  names(plotData) <- rownames(df)
  
  
  
  ord <- order(df$line,df$day)
  df <- df[ord,]
  
  head(plotData[,ord])
  
  mycolors <- a1col[1:length(unique(a1))]
  names(mycolors) <- unique(a1)
  
  mycolors2 <- a2col[1:length(unique(a2))]
  names(mycolors2) <- unique(a2)
  
  mycolors <- list(a = mycolors,b = mycolors2)
  names(mycolors) <- c(n1,n2)
  
  pheatmap(plotData[,ord], annotation_col = df, annotation_colors = mycolors,border_color = NULL,cluster_rows = T, show_rownames = rowShow, cluster_cols = F,scale = scale,color = heatScaleCol)
  head(plotData[,ord])
  
  return(rownames(plotData))
}

pcut <- 0.01
fcl <- 0
sdCut <- .2

## UP human genes 
sign.hips6 <- getSignName(x = dabs.hips6$test$Default,p = pcut,l = fcl)
genes <- sign.hips6$up

# Human data
data <- hips6.vst.brm
head(data)

a1 <- dabs.hips6$colData$time
a2 <- dabs.hips6$colData$line

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2C_Heat_UpHuman_HumanData_d16.d13_padj.",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,height = 7,width = 5)
plotHeatGenes(col="BuPu",data = data,genes = genes,a1 = a1,a2=a2,n1="day",n2 = 'line',sdCut=sdCut,z=T)
dev.off()

# Chimp data
data <- csan.vst.brm
head(data)

a1 <- dabs.cSandra$colData$time
a2 <- dabs.cSandra$colData$line

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2C_Heat_UpHuman_ChimpData_d16.d13_padj.",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,height = 7,width = 5)
plotHeatGenes(col="BuPu",data = data,genes = genes,a1 = a1,a2=a2,n1="day",n2 = 'line',sdCut=sdCut,z=T)
dev.off()

## Down human genes
sign.hips6 <- getSignName(x = dabs.hips6$test$Default,p = pcut,l = fcl)
genes <- sign.hips6$down
sdCut <- sdCut

# Human data
data <- hips6.vst.brm
head(data)

a1 <- dabs.hips6$colData$time
a2 <- dabs.hips6$colData$line

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2C_Heat_DownHuman_HumanData_d16.d13_padj.0.01_FC.0.pdf",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,height = 7,width = 5)
plotHeatGenes(col="RdPu",data = data,genes = genes,a1 = a1,a2=a2,n1="day",n2 = 'line',sdCut=sdCut,z=T,rev=T)
dev.off()

# Chimp data
data <- csan.vst.brm
head(data)

a1 <- dabs.cSandra$colData$time
a2 <- dabs.cSandra$colData$line

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2C_Heat_DownHuman_ChimpData_d16.d13_padj.0.01_FC.0.pdf",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,height = 7,width = 5)
plotHeatGenes(col="RdPu",data = data,genes = genes,a1 = a1,a2=a2,n1="day",n2 = 'line',sdCut=sdCut,z=T,rev=T)
dev.off()


##########

## Plot correlated changes, d16-d13 human and chimp #####
rna.padj <- 1e-3
rna.fc <- 0
baseMean <- 200

combine <- merge(data.frame(dabs.hips6$test$Default),data.frame(dabs.cSandra$test$Default),by=0)
head(combine)
combine.exp <- combine[combine$baseMean.x>baseMean,]
col <- ifelse(test = combine.exp$padj.x<rna.padj & abs(combine.exp$log2FoldChange.x)>rna.fc,
              # if human is significant
              yes = ifelse(test = combine.exp$padj.y<rna.padj & abs(combine.exp$log2FoldChange.y)>rna.fc,
                           # prot sign and ALSO rna 
                           yes = ifelse(test = c(combine.exp$log2FoldChange.x>0)==c(combine.exp$log2FoldChange.y>0),
                                        # Both sign AND same direction: GREEN
                                        yes = "forestgreen",
                                        # Both sign NOT same direction: RED
                                        no = "red"),
                           # prot sign, but NOT RNA: BLUE
                           no = "blue"),
              # human is not sign
              no = ifelse(test = combine.exp$padj.y<rna.padj & abs(combine.exp$log2FoldChange.y)>rna.fc,
                          # Prot NOT sign, RNA sign : 
                          yes = "goldenrod",
                          # NO sign
                          no = "black"))
cex <- ifelse(col=="black",.2,.4)
table(col)

#  4. Plot FC - FC correlation + color sign. #####
filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2D_FCFCcorr_signCol_d16.d13_padj.",rna.padj,"_FC.",rna.fc,"_baseMean.",baseMean,".pdf",sep = "")
pdf(file = filename)
c1 <- combine.exp$log2FoldChange.x
c2 <- combine.exp$log2FoldChange.y
plot(c1,c2,ylim=c(-5,5),xlim=c(-5,5),pch=16,cex=cex,col=col,ylab="Human: log2fc(human/chimp)",xlab="Chimp: log2fc(human/chimp)")
cor <- cor.test( x = c1, y = c2,method = "pearson",conf.level = 0.95)
cor$estimate
fit<-lm(c1~c2)
mtext(text = paste("BaseMean: ",baseMean,", RNA padj: ",rna.padj," log2FC: RNA_",rna.fc,sep = ""),cex = .5)
mtext(text = "Correlation log2FC (d16 vs d13) Human vs Chimp",line=1)
grid(nx = 10,ny = 10)
abline(fit,lty=2,col="red",lwd=1)
legend("topleft",bty='n',legend=paste("r^2: ",format(summary(fit)$r.squared,digits = 3),"\nCorr (pearson): ",format(cor$estimate,digits = 3)))
legend("bottomright",legend = paste(c("non-signficant","Chimp-only","Both-sameChange","Human-only","Both-oppositeChange"),": ",table(col)),col = names(table(col)),pch=16,cex=.6)
points(combine.exp$log2fc.species,combine.exp$log2FoldChange,col=col,pch=16,cex=cex)
points(combine.exp$log2fc.species,combine.exp$log2FoldChange,col=col,pch=16,cex=ifelse(col=="forestgreen",yes = .4,no = 0))
points(combine.exp$log2fc.species,combine.exp$log2FoldChange,col=col,pch=16,cex=ifelse(col=="red",yes = .4,no = 0))
dev.off()


#  5. Plot FC - FC correlation without color sign. #####
combine <- merge(data.frame(dabs.hips6$test$Default),data.frame(dabs.cSandra$test$Default),by=0)
head(combine)
baseMean <- 200
combine.exp <- combine[combine$baseMean.x>baseMean,]

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2D_FCFCcorr_noCol.pdf_baseMean.",baseMean,".pdf",sep = "")
pdf(file = filename)
c1 <- combine.exp$log2FoldChange.x
c2 <- combine.exp$log2FoldChange.y
plot(c1,c2,pch=16,cex=.4,col="black",ylab="Human: log2fc(human/chimp)",xlab="Chimp: log2fc(human/chimp)")
cor <- cor.test( x = c1, y = c2,method = "pearson",conf.level = 0.95)
cor$estimate
fit<-lm(c1~c2)
mtext(text = "Correlation log2FC (d16 vs d13) Human vs Chimp",line=1)
grid(nx = 10,ny = 10)
abline(fit,lty=2,col="red",lwd=1)
legend("topleft",legend=paste("Corr (pearson): ",format(cor$estimate,digits = 3)))
dev.off()


### 6. Supp Fig: PCA ######

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

## ProteinCoding
inputCounts <- "~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/data/merged_pt6_hg38_mapping/mergedMapping_PC_pt6_hg38.s2.multi.Gencode27.Exon.sjdb.txt"
### read RNA data
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

#### Define metadata
##### Note that 'condition' and 'samples' are needed here. Conditions define the design of DESeq object, defines what can be tested for and how the data is normalized.

colDat <- read.delim(file = "data/merged_pt6_hg38_mapping/mergedMapping_pt6_hg38_colData.txt",sep="\t")
colDat

cbind(colDat,samples[-c(1:6)])

### All iPS 

colDat.ips <- colDat[colDat$line!="HS_h9",]
colDat.ips$condition <- paste(colDat.ips$time,colDat.ips$species,sep = "_")

rna.ips <- rna.input[,c(1:6,grep(paste(colDat.ips$samples,collapse = "|"),colnames(rna.input)))]
head(rna.ips)
header <- read.delim(file = path,header = F,nrows = 1)
write.table(header,file = "data/rna.raw.ips.txt",quote = F,row.names = F,col.names = F,sep="\t")
write.table(rna.ips,file = "data/rna.raw.ips.txt",quote = F,row.names = F,col.names = T,sep="\t",append = T)
path.ips <- "data/rna.raw.ips.txt"

dabs.ips <- deseqAbs$new("humanChimp",colData=colDat.ips,filename=path.ips,design = formula(~batch+condition))
head(dabs.ips$rawCounts)
dabs.ips$makeVST(blind=F)

dabs.ips.vst.br <- limma::removeBatchEffect(assay(dabs.ips$VST),batch = colDat.ips$batch)

### PCA batch corrected
pca <- prcomp(x = t(dabs.ips.vst.br))
x <- barplot(100*(pca$sdev/sum(pca$sdev)),ylab="% Variance")
axis(side = 1,at = x,labels = 1:length(x))

pch <- ifelse(test = dabs.ips$colData$species=="human" ,
              yes = ifelse(test = dabs.ips$colData$line=="HS_H6",
                           yes = 16,no = 17),
              no = ifelse(test = dabs.ips$colData$line=="PT_C6",
                          yes = 16,no = 17))
col <- ifelse(test = dabs.ips$colData$species=="chimp",
              yes = ifelse(test = dabs.ips$colData$time == "d13",yes = "lightskyblue1",
                           no = ifelse(test = dabs.ips$colData$time == "d14",yes = "lightskyblue",
                                       no = ifelse(test = dabs.ips$colData$time == "d15",yes = "slateblue1",
                                                   no = "slateblue4"))),
              no =  ifelse(test = dabs.ips$colData$time == "d13",yes = "peachpuff",
                           no = ifelse(test = dabs.ips$colData$time == "d14",yes = "sienna1",
                                       no = ifelse(test = dabs.ips$colData$time == "d15",yes = "red",
                                                   no = "red4"))))

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/FigSuppPCA_PC1.2_PC1.3_BatchCorr.pdf"
pdf(file = filename,width = 10,height = 6)
par(mfrow=c(1,2))
plot(pca$x,cex=2.2,col=col,pch=pch,main = "PCA (PC1 vs PC2) RNA-seq",
     xlab=c(paste("PC1 (",round(x = 100*pca$sdev[1]/sum(pca$sdev),digits = 2),"%)",sep = "")),
     ylab=c(paste("PC2 (",round(x = 100*pca$sdev[2]/sum(pca$sdev),digits = 2),"%)",sep = "")),lwd=2)
points(pca$x[dabs.ips$colData$batch=="batch_2",],col="black",pch=1,cex=2,lwd=1,lty=2)
legend(x = -20,y = 45,legend = c("Human - H6","Human - H48","Chimp - C6","Chimp - CPT"),col = rep(c("sienna1","lightskyblue"),each=2),pch=c(16,17,16,17),pt.cex = 2)

plot(pca$x[,c(1,3)],cex=2.2,col=col,pch=pch,main = "PCA (PC1 vs PC3) RNA-seq",
     xlab=c(paste("PC1 (",round(x = 100*pca$sdev[1]/sum(pca$sdev),digits = 2),"%)",sep = "")),
     ylab=c(paste("PC3 (",round(x = 100*pca$sdev[3]/sum(pca$sdev),digits = 2),"%)",sep = "")),lwd=2)
points(pca$x[dabs.ips$colData$batch=="batch_2",c(1,3)],col="black",pch=1,cex=2,lwd=1,lty=2)
dev.off()

### PCA NOT batch corrected
pca <- prcomp(x = t(assay(dabs.ips$VST)))
x <- barplot(100*(pca$sdev/sum(pca$sdev)),ylab="% Variance")
axis(side = 1,at = x,labels = 1:length(x))

pch <- ifelse(test = dabs.ips$colData$species=="human" ,
              yes = ifelse(test = dabs.ips$colData$line=="HS_H6",
                           yes = 16,no = 17),
              no = ifelse(test = dabs.ips$colData$line=="PT_C6",
                          yes = 16,no = 17))
col <- ifelse(test = dabs.ips$colData$species=="chimp",
              yes = ifelse(test = dabs.ips$colData$time == "d13",yes = "lightskyblue1",
                           no = ifelse(test = dabs.ips$colData$time == "d14",yes = "lightskyblue",
                                       no = ifelse(test = dabs.ips$colData$time == "d15",yes = "slateblue1",
                                                   no = "slateblue4"))),
              no =  ifelse(test = dabs.ips$colData$time == "d13",yes = "peachpuff",
                           no = ifelse(test = dabs.ips$colData$time == "d14",yes = "sienna1",
                                       no = ifelse(test = dabs.ips$colData$time == "d15",yes = "red",
                                                   no = "red4"))))

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/FigSuppPCA_PC1.2_PC1.3_notBatchCorr.pdf"
pdf(file = filename,width = 10,height = 6)
par(mfrow=c(1,2))
plot(pca$x,cex=2.2,col=col,pch=pch,main = "PCA (PC1 vs PC2) RNA-seq",
     xlab=c(paste("PC1 (",round(x = 100*pca$sdev[1]/sum(pca$sdev),digits = 2),"%)",sep = "")),
     ylab=c(paste("PC2 (",round(x = 100*pca$sdev[2]/sum(pca$sdev),digits = 2),"%)",sep = "")),lwd=2)
points(pca$x[dabs.ips$colData$batch=="batch_2",],col="black",pch=1,cex=2,lwd=1,lty=2)

plot(pca$x[,c(1,3)],cex=2.2,col=col,pch=pch,main = "PCA (PC1 vs PC3) RNA-seq",
     xlab=c(paste("PC1 (",round(x = 100*pca$sdev[1]/sum(pca$sdev),digits = 2),"%)",sep = "")),
     ylab=c(paste("PC3 (",round(x = 100*pca$sdev[3]/sum(pca$sdev),digits = 2),"%)",sep = "")),lwd=2)
points(pca$x[dabs.ips$colData$batch=="batch_2",c(1,3)],col="black",pch=1,cex=2,lwd=1,lty=2)
legend("bottomleft",legend = c("Human - H6","Human - H48","Chimp - C6","Chimp - CPT"),col = rep(c("sienna1","lightskyblue"),each=2),pch=c(16,17,16,17),pt.cex = 2)

dev.off()

#### BOxplots d16-d13 diff #### 
library(reshape2)
#dabs.hips6$getAverage()

colnames(hips6.vst.brm)
hips6.vst.brm.mean <- data.frame(d13=rowMeans(hips6.vst.brm[,grep("d13",colnames(hips6.vst.brm))]),
                                 d14=rowMeans(hips6.vst.brm[,grep("d14",colnames(hips6.vst.brm))]),
                                 d15=rowMeans(hips6.vst.brm[,grep("d15",colnames(hips6.vst.brm))]),
                                 d16=rowMeans(hips6.vst.brm[,grep("d16",colnames(hips6.vst.brm))]))
csan.vst.brm.mean <- data.frame(d13=rowMeans(csan.vst.brm[,grep("d13",colnames(csan.vst.brm))]),
                                 d14=rowMeans(csan.vst.brm[,grep("d14",colnames(csan.vst.brm))]),
                                 d15=rowMeans(csan.vst.brm[,grep("d15",colnames(csan.vst.brm))]),
                                 d16=rowMeans(csan.vst.brm[,grep("d16",colnames(csan.vst.brm))]))

## VST 
human.up <- getGenes(data = hips6.vst.brm.mean,genes = sign.hips6$up)
human.down <- getGenes(data = hips6.vst.brm.mean,genes = sign.hips6$down)

chimp.up <- getGenes(data = csan.vst.brm.mean,genes = sign.hips6$up)
chimp.down <- getGenes(data = csan.vst.brm.mean,genes = sign.hips6$down)

human.up.sc <- t(apply(X = as.matrix(human.up[,-1]),MARGIN = 1,FUN = scale))
rownames(human.up.sc) <- human.up$genes
colnames(human.up.sc) <- colnames(human.up[,-1])
human.down.sc <- t(apply(X = as.matrix(human.down[,-1]),MARGIN = 1,FUN = scale))
rownames(human.down.sc) <- human.down$genes
colnames(human.down.sc) <- colnames(human.up[,-1])
chimp.up.sc <- t(apply(X = as.matrix(chimp.up[,-1]),MARGIN = 1,FUN = scale))
rownames(chimp.up.sc) <- chimp.up$genes
colnames(chimp.up.sc) <- colnames(human.up[,-1])
chimp.down.sc <- t(apply(X = as.matrix(chimp.down[,-1]),MARGIN = 1,FUN = scale))
rownames(chimp.down.sc) <- chimp.down$genes
colnames(chimp.down.sc) <- colnames(human.up[,-1])

hs.up.melt <- melt(data = human.up.sc)
hs.down.melt <- melt(data = human.down.sc)
pt.up.melt <- melt(data = chimp.up.sc)
pt.down.melt <- melt(data = chimp.down.sc)

hs.up.melt$day <- read.table(text = as.character(hs.up.melt$Var2), sep = ".", as.is = TRUE)$V1
pt.up.melt$day <- read.table(text = as.character(pt.up.melt$Var2), sep = ".", as.is = TRUE)$V1
hs.down.melt$day <- read.table(text = as.character(hs.down.melt$Var2), sep = ".", as.is = TRUE)$V1
pt.down.melt$day <- read.table(text = as.character(pt.down.melt$Var2), sep = ".", as.is = TRUE)$V1

hs.up.melt$variable <- factor(x = hs.up.melt$Var1,levels = levels(hs.up.melt$Var1)[c(1,5,2,6,3,7,4,8)])
hs.down.melt$variable <- factor(x = hs.down.melt$Var1,levels = levels(hs.down.melt$Var1)[c(1,5,2,6,3,7,4,8)])
pt.up.melt$variable <- factor(x = pt.up.melt$Var1,levels = levels(pt.up.melt$Var1)[c(1,5,2,6,3,7,4,8)])
pt.down.melt$variable <- factor(x = pt.down.melt$Var1,levels = levels(pt.down.melt$Var1)[c(1,5,2,6,3,7,4,8)])

head(hs.up.melt)
head(pt.up.melt)

#install.packages("ggthemes") # Install 
library(ggthemes)

hs.up.melt$group <- paste(hs.up.melt$day,"HS",sep = "_")
pt.up.melt$group <- paste(pt.up.melt$day,"PT",sep = "_")

up.df <- rbind(hs.up.melt,pt.up.melt)
tail(up.df)
##### Boxplot 
filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2Ca_boxplot_zScore_humanUp_human_d16.d13_padj.",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,width = 6,height = 5)
p <- ggplot(data = up.df,mapping = aes(x = group,y = value,fill=group,group=group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90,hjust = 1)) +
  ylab(label = "Z-score") + 
  xlab(label = "day") +
  scale_fill_manual(values = rep(brewer.pal(n = 4,name = "BuPu"),each=2)) +
  ggtitle(label = "Z-score: genes up in human") 
print(p)
dev.off()

# human down
hs.down.melt$group <- paste(hs.down.melt$day,"HS",sep = "_")
pt.down.melt$group <- paste(pt.down.melt$day,"PT",sep = "_")

down.df <- rbind(hs.down.melt,pt.down.melt)
filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2Ca_boxplot_zScore_humanDOWN_human_d16.d13_padj.",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,width = 6,height = 5)
p <- ggplot(data = down.df,mapping = aes(x = group,y = value,fill=group,group=group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90,hjust = 1)) +
  ylab(label = "Z-score") + 
  xlab(label = "day") +
  scale_fill_manual(values = rep(brewer.pal(n = 4,name = "BuPu"),each=2)) +
  ggtitle(label = "Z-score: genes down in human") 
print(p)
dev.off()

# chimp up (human genes up)
filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2Ca_boxplot_humanUP_chimp_d16.d13_padj.",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,width = 6,height = 5)
p <- ggplot(data = pt.up.melt,mapping = aes(x = variable,y = value,fill=day)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90,hjust = 1)) +
  ylab(label = "VST") + 
  scale_fill_manual(values = c("grey100","grey75","grey50","grey25")) +
  ggtitle(label = "Abundance in Chimp: genes up in human") +
print(p)
dev.off()

# chimp down (human genes down)
filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig2-bulkRNA/Fig2Ca_boxplot_humanDOWN_chimp_d16.d13_padj.",pcut,"_FC.",fcl,".pdf",sep = "")
pdf(file = filename,width = 6,height = 5)
p <- ggplot(data = pt.down.melt,mapping = aes(x = variable,y = value,fill=day)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90,hjust = 1)) +
  ylab(label = "VST") + 
  scale_fill_manual(values = c("grey100","grey75","grey50","grey25")) +
  ggtitle(label = "Abundance in Chimp: genes down in human") +
  coord_cartesian(ylim = c(7,15))
print(p)
dev.off()

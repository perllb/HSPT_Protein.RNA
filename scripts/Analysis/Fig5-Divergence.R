## Figure 5 - Divergence human vs chimp #####

rm(list=ls())

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

# plotting norm data
boxplot(prot.norm)
head(prot.norm)

# Preparing data for limma
human.id <- colnames(prot.norm)[1:6]
chimp.id <- colnames(prot.norm)[7:12]

design <- model.matrix(~factor(c(2,2,2,2,2,2,1,1,1,1,1,1)))
colnames(design) <- c("Intercept","Diff")

res.eb <- eb.fit(prot.norm[, c(human.id,chimp.id)], design)
head(res.eb)

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

###################

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

########

plotDir <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/"

# Fig5A_ abundance human vs chimp RNA #####

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

par(mfrow=c(1,1))
# Fig5A: RNA: Human vs Chimp ####
filename <- paste(plotDir,"Fig5A_RNA_human.vs.chimp.pdf",sep = "")
pdf(filename,width = 6,height = 6)
x <- merge.mean.species$RNA_human
y <- merge.mean.species$RNA_chimp
plotCor(x,y,"","human (H6 and H48)","chimp (C6 and CPT)","RNA - det Protein")
dev.off()

########

# Fig5B_ maPlot RNA #####

prot.det <- as.character(rownames(prot.norm))
dabs.ips.d14$makeDiffex()


df.protdet <- dabs.ips.d14$test$Default[rownames(dabs.ips.d14$test$Default) %in% prot.det,]
head(df.protdet)
dim(df.protdet)

maLimit <- 0
p <- 1e-3
l <- 0
limit <- 0

df <- data.frame(baseMean = log2(df.protdet$baseMean),log2FC=df.protdet$log2FoldChange,
                      sign=factor(ifelse(test = !is.na(df.protdet$padj) & df.protdet$padj<p,yes = 1,no = 0)),
                      signUp=factor(ifelse(test = !is.na(df.protdet$padj) & df.protdet$padj<p & df.protdet$log2FoldChange>0 ,yes = 1,no = 0)),
                      signDown=factor(ifelse(test = !is.na(df.protdet$padj) & df.protdet$padj<p & df.protdet$log2FoldChange<0 ,yes = 1,no = 0)))
# 2 if sign up, 1 if sign Down
df$col <- factor(as.numeric(as.character(df$sign))+as.numeric(as.character(df$signUp)))
u <- nrow(df[df$signUp==1,])
d <- nrow(df[df$signDown==1,])
n <- nrow(df) - u - d

df$col <- ifelse(test = df$col==2,yes = paste("Up (",u,")",sep = ""),
                 no = ifelse(test = df$col==1,yes = paste("Down (",d,")",sep = ""),no = paste("Neg (",n,")",sep = "")))

col.RNA.up <- "darkorange2"
col.RNA.down <- "cyan3"
col.prot.up <- "darkmagenta"
col.prot.down <- "aquamarine3"

if(maLimit>0){
  df$log2FC <- ifelse(test = df$log2FC>maLimit,
                   yes = maLimit,no = ifelse(test = df$log2FC<(-maLimit),
                                             yes = (-maLimit),no = df$log2FC))
}else{
  df$log2FC <- df$log2FC
}



maxLim <- max(c(abs(df$log2FC[!is.na(df$log2FC)])))

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/Fig5B_RNA_maPlot_padj.",p,"_log2fc_",l,".pdf",sep = "")
pdf(file = filename,width = 6,height = 6)

pRNA <- ggplot(data = df,mapping = aes(x = baseMean,y = log2FC,color=col)) +
  geom_point(alpha=.4,cex=.4) +
  geom_point(data = df[df$signUp==1,],alpha=.4,col=col.RNA.up) +
  geom_point(data = df[df$signDown==1,],alpha=.4,col=col.RNA.down) +
  theme_classic() +
  ggtitle(label = "RNA changes - genes detected on protein level") +
  scale_color_manual(name="Changes",values = c(col.RNA.down,"black",col.RNA.up)) +
  guides(colour = guide_legend(override.aes = list(size=3,alpha=.9))) +
  geom_hline(yintercept = 0,lty=2,col="grey50")
  
print(pRNA)
dev.off()


######

# Fig5C Protien human vs chimp ####
filename <- paste(plotDir,"Fig5C_Protein_human.vs.chimp.pdf",sep = "")
pdf(filename,width = 6,height = 6)
x <- merge.mean.species$Prot_human
y <- merge.mean.species$Prot_chimp
plotCor(x,y,"","human (H6 and H48)","chimp (C6 and CPT)","Protein")
dev.off()

# Fig5D maPlot protein (limma) ####

limma.fc.padj <- res.eb[,c("logFC","q.mod")]
rownames(limma.fc.padj) <- rownames(res.eb)

prot.norm.shift <- prot.norm+abs(min(prot.norm))
baseMean.prot <- rowMeans(prot.norm.shift)
df.protdet <- merge(baseMean.prot,limma.fc.padj,by=0)

head(df.protdet)

colnames(df.protdet) <- c("ProteinName","baseMean","log2FoldChange","padj")
head(df.protdet)

maLimit <- 0
p <- 1e-3
l <- 0
limit <- 0

df <- data.frame(baseMean = df.protdet$baseMean,log2FC=df.protdet$log2FoldChange,
                 sign=factor(ifelse(test = !is.na(df.protdet$padj) & df.protdet$padj<p,yes = 1,no = 0)),
                 signUp=factor(ifelse(test = !is.na(df.protdet$padj) & df.protdet$padj<p & df.protdet$log2FoldChange>0 ,yes = 1,no = 0)),
                 signDown=factor(ifelse(test = !is.na(df.protdet$padj) & df.protdet$padj<p & df.protdet$log2FoldChange<0 ,yes = 1,no = 0)))

# 2 if sign up, 1 if sign Down
df$col <- factor(as.numeric(as.character(df$sign))+as.numeric(as.character(df$signUp)))
u <- nrow(df[df$signUp==1,])
d <- nrow(df[df$signDown==1,])
n <- nrow(df) - u - d

df$col <- ifelse(test = df$col==2,yes = paste("Up (",u,")",sep = ""),
                 no = ifelse(test = df$col==1,yes = paste("Down (",d,")",sep = ""),no = paste("Neg (",n,")",sep = "")))
head(df)

col.RNA.up <- "darkorange2"
col.RNA.down <- "cyan3"
col.prot.up <- "darkmagenta"
col.prot.down <- "aquamarine3"

if(maLimit>0){
  df$log2FC <- ifelse(test = df$log2FC>maLimit,
                      yes = maLimit,no = ifelse(test = df$log2FC<(-maLimit),
                                                yes = (-maLimit),no = df$log2FC))
}else{
  df$log2FC <- df$log2FC
}

maxLim <- max(c(abs(df$log2FC[!is.na(df$log2FC)])))

filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/Fig5D_Protein_maPlot_padj.",p,"_log2fc_",l,".pdf",sep = "")
pdf(file = filename,width = 6,height = 6)

pProt <- ggplot(data = df,mapping = aes(x = baseMean,y = log2FC,color=col)) +
  geom_point(alpha=.4,cex=.4) +
  geom_point(data = df[df$signUp==1,],alpha=.4,col=col.RNA.up) +
  geom_point(data = df[df$signDown==1,],alpha=.4,col=col.RNA.down) +
  theme_classic() +
  ggtitle(label = "Protein changes") +
  scale_color_manual(name="Changes",values = c(col.RNA.down,"black",col.RNA.up)) +
  guides(colour = guide_legend(override.aes = list(size=3,alpha=.9))) +
  geom_hline(yintercept = 0,lty=2,col="grey50")
print(pProt)
dev.off()

########################################

## Supplemental figures #####

# FigS2C - RNA within cell lines ########
head(ips.d14.vst.br)
colnames(ips.d14.vst.br) <- make.names(colDat.ips.d14$line,unique = T)
# human vs human
human.vs.human.corrs <- c()
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_allRepsRNA_human.vs.human.pdf",height = 6,width = 9)
par(mfrow=c(2,3))
for(i in 1:4){
  for(j in 1:4){
    if(i<j){
      print(i)
      print(j)
      x <- ips.d14.vst.br[,i]
      y <- ips.d14.vst.br[,j]
      corr <- plotCor(x,y,"human vs human",colnames(ips.d14.vst.br)[i],colnames(ips.d14.vst.br)[j],"RNA")
      human.vs.human.corrs <- c(human.vs.human.corrs,corr)
    }
  }
}
dev.off()

# human vs chimp
human.vs.chimp.corrs <- c()
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_allRepsRNA_human.vs.chimp.pdf",height = 12,width = 12)
par(mfrow=c(4,4))
for(i in 1:4){
  for(j in 5:8){
    print(i)
    print(j)
    x <- ips.d14.vst.br[,i]
    y <- ips.d14.vst.br[,j]
    corr <- plotCor(x,y,"human vs chimp",colnames(ips.d14.vst.br)[i],colnames(ips.d14.vst.br)[j],"RNA")
    human.vs.chimp.corrs <- c(human.vs.chimp.corrs,corr)
  }
}
dev.off()

chimp.vs.chimp.corrs <- c()
# chimp vs chimp
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_allRepsRNA_chimp.vs.chimp.pdf",height = 9,width = 15)
par(mfrow=c(2,3))
for(i in 5:8){
  for(j in 5:8){
    print(i)
    print(j)
    if(i<j){
      x <- ips.d14.vst.br[,i]
      y <- ips.d14.vst.br[,j]
      corr <- plotCor(x,y,"chimp.vs.chimp",colnames(ips.d14.vst.br)[i],colnames(ips.d14.vst.br)[j],"RNA")
      chimp.vs.chimp.corrs <- c(chimp.vs.chimp.corrs,corr)
    }
  }
}
dev.off()

# plot Pearson's corr
library(ggplot2)
library(wesanderson)
df.corr <- data.frame(corr=c(human.vs.human.corrs,human.vs.chimp.corrs,chimp.vs.chimp.corrs),group=c(rep("HS.vs.HS",length(human.vs.human.corrs)),rep("HS.vs.PT",length(human.vs.chimp.corrs)),rep("PT.vs.PT",length(chimp.vs.chimp.corrs))))

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_RNA_Correlation_within.cell.lines.pdf"
pdf(filename,width = 5,height = 6)
ggplot(data = df.corr,mapping = aes(x = group,y = corr,fill=group)) +
  geom_violin() +
  geom_boxplot(width=.1,fill="white") +
  scale_fill_manual(values = wes_palette("Royal1")) +
  coord_cartesian(ylim = c(0.7, 1)) +
  ggtitle(label = "Correlations within cell lines")
dev.off()

########


# Protein 
# human vs human
human.vs.human.corrs <- c()
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_allRepsProtein_human.vs.human.pdf",height = 9,width = 15)
par(mfrow=c(3,5))
for(i in 1:6){
  for(j in 1:6){
    if(i<j){
      print(i)
      print(j)
      x <- prot.norm[,i]
      y <- prot.norm[,j]
      corr <- plotCor(x,y,"human vs human",names(prot.norm)[i],names(prot.norm)[j],"Protein")
      human.vs.human.corrs <- c(human.vs.human.corrs,corr)
    }
  }
}
dev.off()

# human vs chimp
human.vs.chimp.corrs <- c()
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_allRepsProtein_human.vs.chimp.pdf",height = 15,width = 15)
par(mfrow=c(6,6))
for(i in 1:6){
  for(j in 7:12){
    print(i)
    print(j)
    x <- prot.norm[,i]
    y <- prot.norm[,j]
    corr <- plotCor(x,y,"human vs chimp",names(prot.norm)[i],names(prot.norm)[j],"Protein")
    human.vs.chimp.corrs <- c(human.vs.chimp.corrs,corr)
  }
}
dev.off()

chimp.vs.chimp.corrs <- c()
# chimp vs chimp
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_allRepsProtein_chimp.vs.chimp.pdf",height = 9,width = 15)
par(mfrow=c(3,6))
for(i in 7:12){
  for(j in 7:12){
    print(i)
    print(j)
    if(i<j){
      x <- prot.norm[,i]
      y <- prot.norm[,j]
      corr <- plotCor(x,y,"chimp vs chimp",names(prot.norm)[i],names(prot.norm)[j],"Protein")
      chimp.vs.chimp.corrs <- c(chimp.vs.chimp.corrs,corr)
    }
  }
}
dev.off()

# plot Pearson's corr
library(ggplot2)
library(wesanderson)
df.corr <- data.frame(corr=c(human.vs.human.corrs,human.vs.chimp.corrs,chimp.vs.chimp.corrs),group=c(rep("HS.vs.HS",length(human.vs.human.corrs)),rep("HS.vs.PT",length(human.vs.chimp.corrs)),rep("PT.vs.PT",length(chimp.vs.chimp.corrs))))

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_Protein_Correlation_within.cell.lines.pdf"
pdf(filename,width = 5,height = 6)
ggplot(data = df.corr,mapping = aes(x = group,y = corr,fill=group)) +
  geom_violin() +
  geom_boxplot(width=.1,fill="white") +
  scale_fill_manual(values = wes_palette("Royal1")) +
  coord_cartesian(ylim = c(0.7, 1)) +
  ggtitle(label = "Correlations within cell lines")
dev.off()

########

# Fig S2B - between cell lines - intraspecies - RNA #######
head(merge.mean)

# RNA
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_RNA_cellLineMeans.pdf",height = 9,width = 6)
par(mfrow=c(3,2))
x <- merge.mean$RNA_HS_H6
y <- merge.mean$RNA_HS_H48
plotCor(x,y,"Human vs Human","human H6","human H48","RNA")

x <- merge.mean$RNA_PT_C6
y <- merge.mean$RNA_PT_CPT
plotCor(x,y,"Chimp vs Chimp","chimp C6","chimp CPT","RNA")

x <- merge.mean$RNA_HS_H6
y <- merge.mean$RNA_PT_C6
plotCor(x,y,"Human vs Chimp","human H6","chimp C6","RNA")

x <- merge.mean$RNA_HS_H6
y <- merge.mean$RNA_PT_CPT
plotCor(x,y,"Human vs Chimp","human H6","chimp CPT","RNA")

x <- merge.mean$RNA_HS_H48
y <- merge.mean$RNA_PT_C6
plotCor(x,y,"Human vs Chimp","human H48","chimp C6","RNA")

x <- merge.mean$RNA_HS_H48
y <- merge.mean$RNA_PT_CPT
plotCor(x,y,"Human vs Chimp","human H48","chimp CPT","RNA")
dev.off()

# Protein
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_Protein_cellLineMeans.pdf",height = 9,width = 6)
par(mfrow=c(3,2))
x <- merge.mean$Prot_H6
y <- merge.mean$Prot_H48
plotCor(x,y,"Human","human H6","human H48","Protein")

# chimp Protein
x <- merge.mean$Prot_C6
y <- merge.mean$Prot_CPT
plotCor(x,y,"Chimp","chimp C6","chimp CPT","Protein")

# human vs chimp Protein
x <- merge.mean$Prot_H6
y <- merge.mean$Prot_C6
plotCor(x,y,"human vs chimp","human H6","chimp C6","Protein")

# human vs chimp Protein
x <- merge.mean$Prot_H6
y <- merge.mean$Prot_CPT
plotCor(x,y,"human vs chimp","human H6","chimp CPT","Protein")

# human vs chimp Protein
x <- merge.mean$Prot_H48
y <- merge.mean$Prot_C6
plotCor(x,y,"human vs chimp","human H48","chimp C6","Protein")

# human vs chimp Protein
x <- merge.mean$Prot_H48
y <- merge.mean$Prot_CPT
plotCor(x,y,"human vs chimp","human H48","chimp CPT","Protein")

dev.off()



## SupFig2: RNA vs Protien #####

# Human
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig5-Divergence_hs_pt/FigS_RNA.vs.Protein.pdf",height = 5,width = 10)
par(mfrow=c(1,2))
x <- merge.mean.species$Prot_human
y <- merge.mean.species$RNA_human
plotCor(x,y,"RNA vs Protein","Protein","RNA","Human")

# Chimp
x <- merge.mean.species$Prot_chimp
y <- merge.mean.species$RNA_chimp
plotCor(x,y,"RNA vs Protein","Protein","RNA","Chimp")
dev.off()

#######





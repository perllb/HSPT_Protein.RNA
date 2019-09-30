## Fig 6 - dependency of RNA ####

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

# plotting norm data
#boxplot(prot.norm)
head(prot.norm)

# Preparing data for limma
human.id <- colnames(prot.norm)[1:6]
chimp.id <- colnames(prot.norm)[7:12]

design <- model.matrix(~factor(c(2,2,2,2,2,2,1,1,1,1,1,1)))
colnames(design) <- c("Intercept","Diff")

res.eb <- eb.fit(prot.norm[, c(human.id,chimp.id)], design)
head(res.eb)
dim(res.eb)
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
inputCounts <- "~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/data/merged_pt6_hg38_mapping/merge/Merged_MultiMapping_Hg38.PanTro6_ProteinCoding_iPS.txt"

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

colDat <- read.delim(file = "data/merged_pt6_hg38_mapping/merge/Merged_MultiMapping_Hg38.PanTro6_ColData.txt",sep="\t")
colDat <- colDat[colDat$line != "HS_h9",]

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

#### Prepare data #####
# protein
limma.fc.padj <- res.eb[,c("logFC","p.mod","q.mod")]
rownames(limma.fc.padj) <- rownames(res.eb)

prot.norm.shift <- prot.norm+abs(min(prot.norm))
baseMean.prot <- rowMeans(prot.norm.shift)
df.prot.limma <- merge(baseMean.prot,limma.fc.padj,by=0)

head(limma.fc.padj)
head(df.prot.limma)

colnames(df.prot.limma) <- c("ProteinName","baseMean","log2FoldChange","padj","q.mod")
head(df.prot.limma)

# RNA ####

dabs.ips.d14$makeDiffex()
dim(prot.norm)
prot.det <- rownames(prot.norm)
setdiff(prot.det,rownames(dabs.ips.d14$test$Default))

###################


df <- dabs.ips.d14$test$Default[rownames(dabs.ips.d14$test$Default) %in% prot.det.mod,]
head(df)
dim(df)

### Make plot 
p.rna <- 0.001
p.protein <- p.rna
l <- 0
RNA.up <- df[!is.na(df$padj) & df$padj<p.rna & df$log2FoldChange>l,]
RNA.down <- df[!is.na(df$padj) & df$padj<p.rna & df$log2FoldChange<l,]
Prot.up <- df.prot.limma[!is.na(df.prot.limma$padj) & df.prot.limma$padj<p.protein & df.prot.limma$log2FoldChange>l,]
Prot.down <- df.prot.limma[!is.na(df.prot.limma$padj) & df.prot.limma$padj<p.protein & df.prot.limma$log2FoldChange<l,]

dim(RNA.up)
dim(RNA.down)
dim(Prot.up)
dim(Prot.down)

Prot_RNA.up <- merge(rownames(RNA.up),df.prot.limma,by=1)
Prot_RNA.down <- merge(rownames(RNA.down),df.prot.limma,by=1)
RNA_prot.up <- merge(Prot.up$ProteinName,df,by.x=1,by.y=0)
RNA_prot.down <- merge(Prot.down$ProteinName,df,by.x=1,by.y=0)

## Test ### 
hist(RNA_prot.up$log2FoldChange,breaks = 50)
hist(RNA_prot.down$log2FoldChange,breaks = 50)
hist(Prot_RNA.down$log2FoldChange,breaks = 50)
hist(Prot_RNA.up$log2FoldChange,breaks = 50)

t.test(RNA_prot.down$log2FoldChange,RNA_prot.up$log2FoldChange)
t.test(Prot_RNA.down$log2FoldChange,Prot_RNA.up$log2FoldChange)

df.fc <- data.frame(log2FC=c(Prot_RNA.up$log2FoldChange,Prot_RNA.down$log2FoldChange,RNA_prot.up$log2FoldChange,RNA_prot.down$log2FoldChange))
df.fc$genes <- c(rep("RNA-up",length(Prot_RNA.up$log2FoldChange)),rep("RNA-down",length(Prot_RNA.down$log2FoldChange)),rep("Protein-up",length(RNA_prot.up$log2FoldChange)),rep("Protein-down",length(RNA_prot.down$log2FoldChange)))
df.fc$moleculeFC <- ifelse(grepl("Protein",df.fc$genes),yes = "RNA",no = "Protein")
df.fc$moleculeGenes <- ifelse(grepl("Protein",df.fc$genes),yes = "Protein",no = "RNA")

par(mfrow=c(1,1))

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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

library(ggplot2)

# Vioplot
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig6-RNAchanges_corrProt/Fig6A_Vioplot.Boxplot_Sign.Changes_FCotherMolecule.pdf"
pdf(file = filename,width = 4,height = 8)
p1 <- ggplot(data = df.fc,mapping = aes(x = genes,y = log2FC)) +
  geom_abline(slope = 0,intercept = 0,lty=2) +
  geom_violin(aes(fill=moleculeFC)) +
  geom_boxplot(outlier.size = .4,width=.05,color="black") +
  theme_classic() +
  scale_fill_manual(values = wes_palette(name = "BottleRocket2")) +
  ggtitle(label = "Dependence of RNA/Protein changes",subtitle = "Log2FC of one molecule, grouped on sign.changes in the other") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# boxplot without outliers
p2 <- ggplot(data = df.fc,mapping = aes(x = genes,y = log2FC)) +
  geom_abline(slope = 0,intercept = 0,lty=2) +
  geom_boxplot(outlier.shape = NA,aes(fill=moleculeFC)) +
  theme_classic() +
  coord_cartesian(ylim = c(-4.4,4)) +
  scale_fill_manual(values = wes_palette(name = "BottleRocket2")) +
  ggtitle(label = "Dependence of RNA/Protein changes",subtitle = "Log2FC of one molecule, grouped on sign.changes in the other") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

?geom_boxplot

multiplot(p1,p2)
dev.off()



### Plot FC - FC corr  RNA up ####
head(RNA.up)
head(Prot_RNA.up)
merge.RNA.up <- merge(data.frame(RNA.up),Prot_RNA.up,by.x=0,by.y=1)
merge.RNA.up$Dir <- "up"
merge.RNA.down <- merge(data.frame(RNA.down),Prot_RNA.down,by.x=0,by.y=1)
merge.RNA.down$Dir <- "down"
merge.allRNA <- merge(data.frame(dabs.ips.d14$test$Default),df.prot.limma,by.x=0,by.y=1)
merge.allRNA$Dir <- "all"
head(merge.allRNA)
head(merge.RNA.up)
df <- rbind(merge.RNA.down,merge.RNA.up,merge.allRNA)
head(df)

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig6-RNAchanges_corrProt/Fig6B_maPlot_Sign.RNA_Protein_levels.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df,mapping = aes(x = baseMean.y,y = log2FoldChange.y,col=Dir)) +
  geom_point(alpha=.4,size=.6,col="black") +
  geom_point(data = df[df$Dir=="down",],alpha=.7,size=4,col="lightskyblue") +
  geom_point(data = df[df$Dir=="up",],alpha=.7,size=4,col="red") +
  geom_point(data = df[df$Dir=="down",],alpha=.2,size=4,col="lightskyblue") +
  theme_classic() +
  xlab(label = "baseMean Protein") +
  ylab(label = "log2FC Protein")  +
  ggtitle(label = "Protein changes",subtitle = "Genes significant in RNA marked")
dev.off()

## Plot FC - FC corr  Protein up ####
merge.RNA.up <- merge(RNA_prot.up,data.frame(Prot.up),by=1)
merge.RNA.up$Dir <- "up"
merge.RNA.down <- merge(RNA_prot.down,data.frame(Prot.down),by=1)
merge.RNA.down$Dir <- "down"
merge.allRNA <- merge(data.frame(dabs.ips.d14$test$Default),df.prot.limma,by.x=0,by.y=1)
merge.allRNA$Dir <- "all"
colnames(merge.allRNA)[1] <- "x"
head(merge.allRNA)
head(merge.RNA.up)

df <- rbind(merge.RNA.down,merge.RNA.up,merge.allRNA)
head(df)
df$baseMean.x <- log2(df$baseMean.x)

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig6-RNAchanges_corrProt/Fig6B_maPlot_Sign.Protein_RNA_levels.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df,mapping = aes(x = baseMean.x,y = log2FoldChange.x,col=Dir)) +
  geom_point(alpha=.4,size=.6,col="black") +
  geom_point(data = df[df$Dir=="down",],alpha=.4,size=2,col="lightskyblue") +
  geom_point(data = df[df$Dir=="up",],alpha=.4,size=2,col="red") +
  geom_point(data = df[df$Dir=="down",],alpha=.04,size=2,col="lightskyblue") +
  theme_classic() +
  xlab(label = "baseMean RNA") +
  ylab(label = "log2FC RNA") +
  ggtitle(label = "RNA changes",subtitle = "Genes significant in protein marked")
dev.off()



####

df.fc <- data.frame(log2FC=c(Prot_RNA.up$log2FoldChange,Prot_RNA.down$log2FoldChange,RNA_prot.up$log2FoldChange,RNA_prot.down$log2FoldChange))
df.fc$genes <- c(rep("RNA-up",length(Prot_RNA.up$log2FoldChange)),rep("RNA-down",length(Prot_RNA.down$log2FoldChange)),rep("Protein-up",length(RNA_prot.up$log2FoldChange)),rep("Protein-down",length(RNA_prot.down$log2FoldChange)))
df.fc$moleculeFC <- ifelse(grepl("Protein",df.fc$genes),yes = "RNA",no = "Protein")
df.fc$moleculeGenes <- ifelse(grepl("Protein",df.fc$genes),yes = "Protein",no = "RNA")


wilcox.test(x = RNA_prot.down$log2FoldChange,y = RNA_prot.up$log2FoldChange)
wilcox.test(x = Prot_RNA.down$log2FoldChange,y = Prot_RNA.up$log2FoldChange)
t.test(x = RNA_prot.down$log2FoldChange,y = RNA_prot.up$log2FoldChange)
t.test(x = Prot_RNA.down$log2FoldChange,y = Prot_RNA.up$log2FoldChange)

head(df.fc)
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig6-RNAchanges_corrProt/Fig6D_densityPlot_Sign.Protein_RNA.levels.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.fc[df.fc$moleculeFC=="RNA",],mapping = aes(x = log2FC,fill=genes)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept = 0,lty=3) +
  ggtitle(label = "FC density of RNA",subtitle = "Proteins sign.changed")
dev.off()

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig6-RNAchanges_corrProt/Fig6D_densityPlot_Sign.Protein_RNA.levels_coordXlim.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.fc[df.fc$moleculeFC=="RNA",],mapping = aes(x = log2FC,fill=genes)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept = 0,lty=3) +
  coord_cartesian(xlim = c(-2.5,2.5)) +
  ggtitle(label = "FC density of RNA",subtitle = "Proteins sign.changed")
dev.off()

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig6-RNAchanges_corrProt/Fig6D_densityPlot_Sign.RNA_Protein.levels.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.fc[df.fc$moleculeFC=="Protein",],mapping = aes(x = log2FC,fill=genes)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept = 0,lty=3) +
  ggtitle(label = "FC density of Protein",subtitle = "RNA sign.changed")
dev.off()

t.test(RNA_prot.down$log2FoldChange,RNA_prot.up$log2FoldChange)
t.test(Prot_RNA.up$log2FoldChange,Prot_RNA.down$log2FoldChange)

ks.test(RNA_prot.down$log2FoldChange,RNA_prot.up$log2FoldChange)
ks.test(Prot_RNA.up$log2FoldChange,Prot_RNA.down$log2FoldChange)

wilcox.test(RNA_prot.down$log2FoldChange,RNA_prot.up$log2FoldChange)
wilcox.test(Prot_RNA.up$log2FoldChange,Prot_RNA.down$log2FoldChange)


##### Permutation tests #####


df <- dabs.ips.d14$test$Default[rownames(dabs.ips.d14$test$Default) %in% prot.det.mod,]
head(df)
dim(df)


RNA.up <- df[!is.na(df$padj) & df$padj<p.rna & df$log2FoldChange>l,]
RNA.down <- df[!is.na(df$padj) & df$padj<p.rna & df$log2FoldChange<l,]
Prot.up <- df.prot.limma[!is.na(df.prot.limma$padj) & df.prot.limma$padj<p.protein & df.prot.limma$log2FoldChange>l,]
Prot.down <- df.prot.limma[!is.na(df.prot.limma$padj) & df.prot.limma$padj<p.protein & df.prot.limma$log2FoldChange<l,]

dim(RNA.up)
dim(RNA.down)
dim(Prot.up)
dim(Prot.down)

Prot_RNA.up <- merge(rownames(RNA.up),df.prot.limma,by=1)
Prot_RNA.down <- merge(rownames(RNA.down),df.prot.limma,by=1)
RNA_prot.up <- merge(Prot.up$ProteinName,df,by.x=1,by.y=0)
RNA_prot.down <- merge(Prot.down$ProteinName,df,by.x=1,by.y=0)

mean(Prot_RNA.down$log2FoldChange)
median(Prot_RNA.down$log2FoldChange)
mean(Prot_RNA.up$log2FoldChange)
median(Prot_RNA.up$log2FoldChange)

mean(RNA_prot.up$log2FoldChange)
median(RNA_prot.up$log2FoldChange)
mean(RNA_prot.down$log2FoldChange)
median(RNA_prot.down$log2FoldChange)


###### PERMUTATION TESTS #######


perm.tester = function(geneSetIdx,FC,nperm=1000) {
  
  # Function to run permutations
  # calculate difference between selected group (same n as #genes prot up) and all other genes
  one.test = function(x,y) {
    
    xstar<-sample(x)
    mean(y[xstar==1])-mean(y[xstar==0])
    
  }
  
  ## Get difference in true data
  alt.diff <- mean(FC[geneSetIdx==1])-mean(FC[geneSetIdx==0])
  
  ## 1000 permutations
  many <- replicate(n = nperm,expr = one.test(x = geneSetIdx,y = FC))
  hist(many,xlim = c(min(min(many),alt.diff),max(max(many),alt.diff)))
  abline(v=alt.diff, lwd=2, col="purple")
  
  ## get p-value:
  # fraction of higher values in permuted data versus to true data
  return((1+mean(abs(many) > abs(alt.diff)))/nperm)
  
}

# Test RNA FC prot changed genes
head(df)
df.RNAfc <- data.frame(ID=rownames(df),FC=df$log2FoldChange)
df.RNAfc$FC[is.na(df.RNAfc$FC)] <- 0

## Test prot UP genes
prot.up.idx <- ifelse(df.RNAfc$ID %in% RNA_prot.up$x,yes = 1,no = 0)
table(prot.up.idx)
perm.tester(geneSetIdx = prot.up.idx,FC = df.RNAfc$FC,nperm=10000)

## Test prot Down genes
prot.down.idx <- ifelse(df.RNAfc$ID %in% RNA_prot.down$x,yes = 1,no = 0)
table(prot.down.idx)
perm.tester(geneSetIdx = prot.down.idx,FC = df.RNAfc$FC,nperm = 10000)

# Test Protein FC, RNA changed genes
head(df.prot.limma)
df.PROTfc <- data.frame(ID=df.prot.limma$ProteinName,FC=df.prot.limma$log2FoldChange)
df.PROTfc$FC[is.na(df.PROTfc$FC)] <- 0

## Test RNA UP genes
RNA.up.idx <- ifelse(df.PROTfc$ID %in% Prot_RNA.up$x,yes = 1,no = 0)
table(RNA.up.idx)
perm.tester(geneSetIdx = RNA.up.idx,FC = df.PROTfc$FC,nperm=10000)

## Test prot Down genes
RNA.down.idx <- ifelse(df.PROTfc$ID %in% Prot_RNA.down$x,yes = 1,no = 0)
table(RNA.down.idx)
perm.tester(geneSetIdx = RNA.down.idx,FC = df.PROTfc$FC,nperm = 10000)

######



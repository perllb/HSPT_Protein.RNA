## Single cell analysis #####
# Fig1 scRNAseq: chimp #####

# 1. Load and prepare data ####

#rm(list=ls())
setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/20180521_10xInhouse_scRNAseq/seurat/")

#install.packages('Seurat')

library(Seurat)
library(dplyr)

data <- Read10X(data.dir = "~/Dropbox (MN)/Per/PhD/Projects/Chimp/20180521_10xInhouse_scRNAseq/matrices/ciPSSandra_D14_filtered_gene_bc_matrices/GRCh38/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = data))
dense.size
sparse.size <- object.size(x = data)
sparse.size
dense.size/sparse.size

# create object 
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
data <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, 
                           project = "10X_data")

dim(data@data)

# UMI and mito 
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = data@data), value = TRUE)
percent.mito <- Matrix::colSums(data@raw.data[mito.genes, ])/Matrix::colSums(data@raw.data)
data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")

GenePlot(object = data, gene1 = "nUMI", gene2 = "percent.mito",cex.use = .4)
GenePlot(object = data, gene1 = "nUMI", gene2 = "nGene",cex.use = .4)

upper <- 5000
lower <- 2000
data <- FilterCells(object = data, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(lower, -Inf), high.thresholds = c(upper, 0.05))


dim(data@data)

# Normalized Data 
data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Find Variable genes 
par(mfrow=c(1,1))
data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

###  SCale data , run PCA, visualize##
length(x = data@var.genes)

data <- ScaleData(object = data, vars.to.regress = c("nUMI", "percent.mito"))

# run PCA
data <- RunPCA(object = data, pc.genes = data@var.genes,do.print = F)

## run TSNE
data <- RunTSNE(object = data, dims.use = 1:10, do.fast=T ,dim.embed = 3)


# Get and mark cells with NeuroD1 / G1 expression
neurod1 <- data@data["NEUROD1",]#[1:5,1:5]
neurog1 <- data@data["NEUROG1",]
dim(data@data)
length(data@ident)

cells.nd1 <- names(neurod1[neurod1>0])
cells.ng1 <- names(neurog1[neurog1>0])
all <- rownames(data@dr$pca@cell.embeddings)

tsne.d123 <- data.frame(data@dr$tsne@cell.embeddings[,1:3])
tsne.d123$ng1nd1 <- ifelse(rownames(tsne.d123) %in% cells.ng1,yes = ifelse(test = rownames(tsne.d123) %in% cells.nd1,yes = 3,no = 2),no = ifelse(rownames(tsne.d123) %in% cells.nd1,yes = 1,no = 0))
tsne.d123$ng1nd1 <- factor(tsne.d123$ng1nd1)
head(tsne.d123)

###########

# Get and mark cells with ANKRD1 / TAGLN expression
ankrd1 <- data@data["ANKRD1",]#[1:5,1:5]
ctgf <- data@data["CTGF",]
dim(data@data)
length(data@ident)

cells.ankrd1 <- names(ankrd1[ankrd1>0])
cells.ctgf <- names(ctgf[ctgf>0])
all <- rownames(data@dr$pca@cell.embeddings)

tsne.d123 <- data.frame(data@dr$tsne@cell.embeddings[,1:3])
tsne.d123$nd1 <- neurod1
tsne.d123$ng1 <- neurog1
tsne.d123$ankrd1 <- ankrd1
tsne.d123$ctgf <- ctgf
tsne.d123$ng1nd1 <- factor(ifelse(rownames(tsne.d123) %in% cells.ng1,yes = ifelse(test = rownames(tsne.d123) %in% cells.nd1,yes = 3,no = 2),no = ifelse(rownames(tsne.d123) %in% cells.nd1,yes = 1,no = 0)))
tsne.d123$ankr.ctgf <- factor(ifelse(rownames(tsne.d123) %in% cells.ankrd1,yes = ifelse(test = rownames(tsne.d123) %in% cells.ctgf,yes = 3,no = 2),no = ifelse(rownames(tsne.d123) %in% cells.ctgf,yes = 1,no = 0)))
head(tsne.d123)
######

plotDir <- "../../Proteome_RNAseq.comparison/Data_2018-10-19/Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig3-scRNAseq/"


# 2. Fig3E: Pie NeuroD1 G1 ######

### Plot NEUROD1 and G1 cells

library(reshape2)
PCA.d123 <- data.frame(data@dr$pca@cell.embeddings[,1:5])
PCA.d123$ng1nd1 <- ifelse(rownames(PCA.d123) %in% cells.ng1,yes = ifelse(test = rownames(PCA.d123) %in% cells.nd1,yes = 3,no = 2),no = ifelse(rownames(PCA.d123) %in% cells.nd1,yes = 1,no = 0))
head(PCA.d123)
PCA.d123$ng1nd1 <- factor(PCA.d123$ng1nd1)

#### Plot number of NeuroD1 / G1 + cells
nd1.only <- setdiff(cells.nd1,cells.ng1)
ng1.only <- setdiff(cells.ng1,cells.nd1)
both.nd1.g1 <- intersect(cells.nd1,cells.ng1)
all <- rownames(data@dr$pca@cell.embeddings)

tmp <- setdiff(all,nd1.only)
tmp1 <- setdiff(tmp,ng1.only)
tmp2 <- setdiff(tmp1,both.nd1.g1)

df <- data.frame(not = length(tmp2),nd1=length(nd1.only),ng1=length(ng1.only),both=length(both.nd1.g1))

filename <- paste(plotDir,"Fig3_chimp_pie_neuroD1.G1.pdf",sep = "")

pdf(file = filename,width = 7,height = 7)
par(mfrow=c(1,1))
pie(c(not = length(tmp2),nd1=length(nd1.only),ng1=length(ng1.only),both=length(both.nd1.g1)),
    labels = c(paste("Negative: ",round(x = 100*length(tmp2)/length(all),digits = 2),"%",sep = ""),
               paste("NeuroD1: ",round(x = 100*length(nd1.only)/length(all),digits = 2),"%",sep = ""),
               paste("NeuroG1: ",round(x = 100*length(ng1.only)/length(all),digits = 2),"%",sep = ""),
               paste("Both: ",round(x = 100*length(both.nd1.g1)/length(all),digits = 2),"%",sep = "")),
    main = "Proportion of NeuroD1 / NeuroG1+ cells",
    col = c("grey70","firebrick","goldenrod1","orange"))
dev.off()

########

#### Plot number of ANKRD1 + TAGLN + cells
ankrd1.only <- setdiff(cells.ankrd1,cells.ctgf)
ctgf.only <- setdiff(cells.ctgf,cells.ankrd1)
both.ankr.ctgf <- intersect(cells.ctgf,cells.ankrd1)
all <- rownames(data@dr$pca@cell.embeddings)

tmp <- setdiff(all,ankrd1.only)
tmp1 <- setdiff(tmp,ctgf.only)
tmp2 <- setdiff(tmp1,both.ankr.ctgf)

df <- data.frame(not = length(tmp2),ctgf=length(ctgf),ankrd1=length(ankrd1.only),both=length(both.ankr.ctgf))

filename <- paste(plotDir,"Fig3_chimp_pie_ANKRD1.CTGF.pdf",sep = "")

pdf(file = filename,width = 7,height = 7)
par(mfrow=c(1,1))
pie(c(not = length(tmp2),tagln=length(ctgf.only),ankrd1=length(ankrd1.only),both=length(both.ankr.ctgf)),
    labels = c(paste("Negative: ",round(x = 100*length(tmp2)/length(all),digits = 2),"%",sep = ""),
               paste("CTGF: ",round(x = 100*length(ctgf.only)/length(all),digits = 2),"%",sep = ""),
               paste("ANKRD1: ",round(x = 100*length(ankrd1.only)/length(all),digits = 2),"%",sep = ""),
               paste("Both: ",round(x = 100*length(both.ankr.ctgf)/length(all),digits = 2),"%",sep = "")),
    main = "Proportion of CTGF / ANKRD1+ cells",
    col = c("grey70","darkolivegreen3","blue","orange"))
dev.off()

# Fig3 D: tSNE- NeuroD1 G1 #######

### tSNE 1 and 2. NeuroD1 G1
filename <- paste(plotDir,"Fig3D_chimp_tSNE_NeuroD1_G1_tsne1.tsne2.pdf",sep = "")
pdf(filename,width = 6,height = 6)
par(mfrow=c(1,1))
ggplot(data = tsne.d123,mapping = aes(x = tSNE_1,y = tSNE_2)) +
  geom_point(col="black",alpha=.1,size=.7) +
  geom_point(data = tsne.d123,mapping = aes(size=neurod1,col="firebrick"),alpha=.5) + 
  geom_point(data = tsne.d123,mapping = aes(size=neurog1,col="goldenrod"),alpha=.5) + 
  geom_point(data = tsne.d123,mapping = aes(size=neurod1),alpha=.1,col="firebrick") + 
  scale_size_continuous(range = c(0,5),limits = c(0.01,4.5),guide_legend(title="log2(reads)")) +
  scale_color_manual(values = c('firebrick','goldenrod'),labels=c("NeuroD1","NeuroG1"))
dev.off()

max(tsne.d123$ctgf)
max(tsne.d123$ankrd1)
max(tsne.d123$nd1)
max(tsne.d123$ng1)
########

# Fig3 Da: tSNE- ANKRD1 CTGF #######

filename <- paste(plotDir,"Fig3D_chimp_tSNE_ANKRD1.CTGF_tsne1.tsne2.pdf",sep = "")
pdf(filename,width = 6,height = 6)
par(mfrow=c(1,1))
ggplot(data = tsne.d123,mapping = aes(x = tSNE_1,y = tSNE_2)) +
  geom_point(col="black",alpha=.1,size=.7) +
  geom_point(data = tsne.d123,mapping = aes(size=ankrd1,col="blue"),alpha=.3) + 
  geom_point(data = tsne.d123,mapping = aes(size=ctgf,col="darkolivegreen3"),alpha=.3) + 
  geom_point(data = tsne.d123,mapping = aes(size=ankrd1),alpha=.1,col="blue") + 
  scale_size_continuous(range = c(0,5),limits = c(.1,4.5),guide_legend(title="log2(reads)")) +
  scale_color_manual(values = c('blue','darkolivegreen3'),labels=c("Ankrd1","Ctgf")) 
dev.off()


## Tsne blak
filename <- paste(plotDir,"Fig3D_chimp_tSNE_black_tsne1.tsne2.pdf",sep = "")
pdf(filename,width = 6,height = 6)
par(mfrow=c(1,1))
ggplot(data = tsne.d123,mapping = aes(x = tSNE_1,y = tSNE_2)) +
  geom_point(col="black",alpha=.8,size=.7)
dev.off()


# Fig3 C: tSNE with FOXG1, PAX6 and OTX2 ####
foxg1 <- data@data["FOXG1",]#[1:5,1:5]
pax6 <- data@data["PAX6",]
otx2 <- data@data["OTX2",]
dim(data@data)
length(data@ident)
df.all <- data.frame(fox=foxg1,pax=pax6,otx=otx2)
df.melt <- melt(df.all)
names(df.melt) <- c("variable","log2(reads)")
head(df.melt)

cells.foxg1 <- names(foxg1[foxg1>0])
cells.pax6 <- names(pax6[pax6>0])
cells.otx2 <- names(otx2[otx2>0])
all <- rownames(data@dr$pca@cell.embeddings)
tsne.d123$fox <- foxg1
tsne.d123$pax <- pax6
tsne.d123$otx <- otx2


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


head(tsne.d123)
filename <- paste(plotDir,"Fig3C_chimp_tSNE_foxg1.pax6.otx2_tsne1.tsne2.pdf",sep = "")
pdf(filename,width = 12,height = 10)

p1 <- ggplot(data = tsne.d123,mapping = aes(x = tSNE_1,y = tSNE_2,size=pax)) +
  geom_point(col="black",alpha=.4,size=.7) +
  geom_point(data = tsne.d123,mapping = aes(size=pax,col="firebrick3"),alpha=.2) + 
  scale_size_continuous(range = c(0,3),limits = c(.1,max(c(tsne.d123$pax,tsne.d123$fox,tsne.d123$otx))+.5),guide_legend(title="log2(reads)")) +
  scale_color_manual(values = c('firebrick3'),labels=c("PAX6")) +
  ggtitle(label = "PAX6")

p2 <-ggplot(data = tsne.d123,mapping = aes(x = tSNE_1,y = tSNE_2)) +
  geom_point(col="black",alpha=.4,size=.7) +
  geom_point(data = tsne.d123,mapping = aes(size=otx,col="darkgreen"),alpha=.2) + 
  scale_size_continuous(range = c(0,3),limits = c(.1,max(c(tsne.d123$pax,tsne.d123$fox,tsne.d123$otx))+.5),guide_legend(title="log2(reads)")) +
  scale_color_manual(values = c('darkgreen'),labels=c("OTX2")) +
  ggtitle(label = "OTX2")

p3 <- ggplot(data = tsne.d123,mapping = aes(x = tSNE_1,y = tSNE_2)) +
  geom_point(col="black",alpha=.4,size=.7) +
  geom_point(data = tsne.d123,mapping = aes(size=fox,col="darkblue"),alpha=.2) + 
  scale_size_continuous(range = c(0,3),limits = c(.1,max(c(tsne.d123$pax,tsne.d123$fox,tsne.d123$otx))+.5),guide_legend(title="log2(reads)")) +
  scale_color_manual(values = c('darkblue'),labels=c("FOXG1")) +
  ggtitle(label = "FOXG1")

multiplot(p1,p2,p3,cols=2)
dev.off()
#######

# Fig3 A: PCA with Cell cycle color ####
# Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "~/Dropbox (MN)/Per/PhD/Projects/Chimp/20180521_10xInhouse_scRNAseq/archive (1 lane only)/seurat/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# 
data <- CellCycleScoring(object = data, s.genes = s.genes, g2m.genes = g2m.genes, 
                         set.ident = TRUE)
#data <- RunPCA(object = data, pc.genes = data@var.genes, pcs.print = 1:4, genes.print = 10)

## PC structure
PCA.d123 <- data.frame(data@dr$pca@cell.embeddings[,1:5])
PCA.d123$cc <- factor(data@meta.data$Phase)
tsne.d123$cc <- factor(data@meta.data$Phase)

filename <- paste(plotDir,"Fig3A_chimp_PCA_CCcolors.pdf",sep = "")
pdf(file = filename)
ggplot(data = PCA.d123,mapping = aes(x = PC1,y = PC2)) +
  geom_point(mapping = aes(col=cc),alpha=.4)
dev.off()

########

# Fig3 B: PCA- regress CC ####
data <- ScaleData(object = data, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = T)

# PCA on var genes - now with regressed out cell cycle differences
data <- RunPCA(object = data, pc.genes = data@var.genes, pcs.print = 1:4)
pca.cc.d123 <- data.frame(data@dr$pca@cell.embeddings[,1:5])
pca.cc.d123$cc <- factor(data@meta.data$Phase)

pdf(file = paste(plotDir,"Fig3B_chimp_PCA_varGenes.regressedCC_colorCC.pdf",sep=""))
ggplot(data = pca.cc.d123,mapping = aes(x = PC1,y = PC2)) +
  geom_point(alpha=.4)
dev.off()


########
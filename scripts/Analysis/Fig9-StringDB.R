## Get number of protein-interactions of proteins #####

rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

## Get protein data 

# Matched by geneName _ all _ normalized and imputed
proteome.input <- read.delim(file = "data/MatchedByGeneNames_HighOnly_normIMP_ModNames.txt",header=T)
head(proteome.input)
names(proteome.input)[1:3] <- c("H48","H48.1","H48.2")

prot.norm <- proteome.input[,1:12]
rownames(prot.norm) <- make.names(proteome.input$X.31,unique = T)

# Fix MT genes ID
prot.det <- as.character(rownames(prot.norm))
prot.det[grep("MT",prot.det)]
prot.det.mod <- gsub(pattern = "^MT\\.",replacement = "MT-",rm(x = prot.det)
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
library(wesanderson)
library(qvalue)
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
inputCounts <- "~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/data/merged_pt6_hg38_mapping/merge/Merged_MultiMapping_Hg38.PanTro6_ProteinCoding.txt"

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
colDat

cbind(colDat,samples[-c(1:6)])

# ips d14 only 
colDat.ips.d14 <- colDat[colDat$line!="HS_h9" & colDat$time=="d14",]
colDat.ips.d14$condition <- colDat.ips.d14$species

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

dabs.ips.d14$makeDiffex()
dabs.ips.d14$test$Default

#######

## Get RNA prot.det ####
prot.det <- rownames(prot.norm)

df.protdet <- dabs.ips.d14$test$Default[rownames(dabs.ips.d14$test$Default) %in% prot.det,]
head(df.protdet)
dim(df.protdet)

## Merge RNA and protein diffex ####
combine <- merge(x = data.frame(df.protdet),y=res.eb,by=0)
head(combine)

plot(log2(combine$baseMean),combine$log2FoldChange,pch=16,cex=.3)

# combine protein mean
protMean <- rowMeans(prot.norm)
combine <- merge(x = combine,y = protMean,by.x=1,by.y=0)
head(combine)

##########

### Read PPI stringdb data ####
library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Read PPI StringDB - freq as computed in unix - all ppis, all scores ####
prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.txt",sep = ",",stringsAsFactors = F,header = F)
## Convert ENSP to gene name #####
conversion.table <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"),
                          filters = "ensembl_peptide_id" , values = list(prot.1$V2),
                          mart = ensembl)
head(conversion.table)


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

plotPPI <- function(prot.1,conversion.table=conversion.table,desc="PPIs") {
  
  ## Merge hgnc symbol and freq table
  sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
  head(sym.freq)
  
  ### Merge protein and RNA stats with PPIs #####
  head(combine)
  head(sym.freq)
  combine.ppi <- merge(combine,sym.freq,by.x=1,by.y=3)
  head(combine.ppi)
  
  ######
  
  #### plot prot FC-FC corr with color on PPI #####
  library(ggplot2)
  
  df <- data.frame(baseMean.prot = combine.ppi$y+abs(min(combine.ppi$y)),
                   baseMean.RNA = log2(combine.ppi$baseMean),
                   log2FC.prot=combine.ppi$logFC,
                   log2FC.RNA=combine.ppi$log2FoldChange,
                   PPI=combine.ppi$V1)
  df <- df[order(df$PPI),]
  head(df)
  
  #filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig7-ProtAbundanceBin/Fig7_Protein_maPlot_bins.pdf"
  #pdf(file = filename,width = 6,height = 5)
  p1 <- ggplot(data = df,mapping = aes(x = log2FC.prot,y = log2FC.RNA,color=PPI)) +
    geom_point(aes(size=PPI,alpha=PPI),pch=16) +
    theme_classic() +
    ggtitle(label = "FC-FC corr",subtitle = paste("number of ",desc,sep = "")) +
    geom_hline(yintercept = 0,lty=2,col="grey50") +
    scale_color_continuous(low = "lightskyblue",high = "darkred") +
    scale_alpha_continuous(range = c(.2,1) ) +
    scale_size_continuous(range = c(.3,3.5))
  
  ## plot protein ma with color of PPI
  p2 <- ggplot(data = df,mapping = aes(x = baseMean.prot,y = log2FC.prot,color=PPI)) +
    geom_point(pch=16,aes(size=PPI,alpha=PPI)) +
    theme_classic() +
    ggtitle(label = "Protein expression",subtitle = paste("number of ",desc,sep = "")) +
    geom_hline(yintercept = 0,lty=2,col="grey50") +
    scale_color_continuous(low = "lightskyblue",high = "darkred") +
    scale_alpha_continuous(range = c(.2,1) ) +
    scale_size_continuous(range = c(.3,3.5))
  
  ## plot RNA ma with color of PPI
  p3 <- ggplot(data = df,mapping = aes(x = baseMean.RNA,y = log2FC.RNA,color=PPI)) +
    geom_point(pch=16,aes(size=PPI,alpha=PPI)) +
    theme_classic() +
    ggtitle(label = "RNA expression",subtitle = paste("number of ",desc,sep = "")) +
    geom_hline(yintercept = 0,lty=2,col="grey50") +
    scale_color_continuous(low = "lightskyblue",high = "darkred") +
    scale_alpha_continuous(range = c(.2,1) ) +
    scale_size_continuous(range = c(.3,3.5))
  
  multiplot(p1,p2,p3,cols = 3)
  
}

## Read PPI StringDB - freq as computed in unix - all ppis, all scores ####
prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/AllPPI.pdf"
pdf(file = filename,width = 15,height = 5)
plotPPI(prot.1 = prot.1,conversion.table = conversion.table,desc = "PPIs (all)")
dev.off()

## Read PPI StringDB - freq as computed in unix - all ppis, score 400 ####
prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/AllPPI_score400.pdf"
pdf(file = filename,width = 15,height = 5)
plotPPI(prot.1 = prot.1,conversion.table = conversion.table,desc = "PPIs (score400)")
dev.off()

## Read PPI StringDB - freq as computed in unix - binding  ppis, score 400 ####
prot.1 <- read.delim(file = "data/StringDB//9606.binding.Acting.Links.v10.5.Freq.Prot1_score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Binding._score400.pdf"
pdf(file = filename,width = 15,height = 5)
plotPPI(prot.1 = prot.1,conversion.table = conversion.table,desc = "Active binding interactions (score400)")
dev.off()

#### Plot PPI in bin (priotien FC) #####
head(combine)

## Merge hgnc symbol and freq table
prot.1 <- read.delim(file = "data/StringDB/9606.binding.Acting.Links.v10.5.Freq.Prot1_score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)
sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
head(sym.freq)
names(sym.freq) <- c("ENSP","PPI_bAct","hgnc_sym")

### Merge protein and RNA stats with PPIs #####
combine.ppi <- merge(combine,sym.freq[-1],by.x=1,by.y=2)
head(combine.ppi)

prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)
sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
head(sym.freq)
names(sym.freq) <- c("ENSP","PPI_alllink","hgnc_sym")

### Merge protein and RNA stats with PPIs #####
df.ppi <- merge(combine.ppi,sym.freq[,-1],by.x=1,by.y=2)
head(df.ppi)

##############

## bin Protein FC #########
df.ppi$fcBin.prot <- ifelse(test = abs(df.ppi$logFC)<=.5,yes = "0-0.5",
                       no = ifelse(test = abs(df.ppi$logFC)>.5 & abs(df.ppi$logFC)<=1,yes = "0.5-1",
                                   no = ifelse(test = abs(df.ppi$logFC)>1 & abs(df.ppi$logFC)<=1.5,yes = "1-1.5",
                                               no = ifelse(test = abs(df.ppi$logFC)>1.5 & abs(df.ppi$logFC)<=2,yes = "1.5-2",
                                                           no = ifelse(test = abs(df.ppi$logFC)>2 & abs(df.ppi$logFC)<=3,yes = "2-3",
                                                                       no = ifelse(test = abs(df.ppi$logFC)>3 & abs(df.ppi$logFC)<=4,yes = "3-4",
                                                                                   no = "4<"))))))
df.ppi$fcBin.prot

df.ppi <- df.ppi[order(abs(df.ppi$logFC)),]
head(df.ppi)

### 9 bins 
df.ppi$ProtfcBinN <- factor(c(rep(1:9,each=c(floor(3874/9))),9,9,9,9))

library(ggplot2)
df.ppi$absLog2Prot <- abs(df.ppi$logFC)
library(tidyverse)
minMax <- df.ppi %>% 
  select(absLog2Prot,ProtfcBinN) %>%
  group_by(ProtfcBinN) %>%
  summarise(min = min(absLog2Prot),max = max(absLog2Prot))

minMax.ppi <- merge(df.ppi$ProtfcBinN,minMax,by.x=1,by.y=1)

tail(df.ppi[df.ppi$ProtfcBinN==9,])
df.ppi$ProtfcBinN.min <- minMax.ppi$min
df.ppi$ProtfcBinN.max <- minMax.ppi$max
df.ppi$ProtfcBinN.minMax <- factor(paste("FC: (",round(df.ppi$ProtfcBinN.min,digits = 2),",",round(df.ppi$ProtfcBinN.max,digits = 2),")",sep = ""))

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_protFCbin9_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = ProtfcBinN.minMax,y = PPI_bAct,fill=ProtfcBinN,group=ProtfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
#  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$ProtfcBinN.minMax)) +
  ggtitle(label = "PPI: all proteins",subtitle = " Binding Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,250)) 
dev.off()

head(df.ppi)
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-allLinks400_protFCbin9_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = ProtfcBinN.minMax,y = PPI_alllink,fill=ProtfcBinN,group=ProtfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$ProtfcBinN.minMax)) +
  ggtitle(label = "PPI: all proteins",subtitle = "All PPI >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,600)) 
dev.off()


### 5 bins 
df.ppi$ProtfcBinN <- factor(c(rep(1:5,each=c(floor(3874/5))),5,5,5,5))
library(ggplot2)
df.ppi$absLog2Prot <- abs(df.ppi$logFC)
library(tidyverse)
minMax <- df.ppi %>% 
  select(absLog2Prot,ProtfcBinN) %>%
  group_by(ProtfcBinN) %>%
  summarise(min = min(absLog2Prot),max = max(absLog2Prot))

minMax.ppi <- merge(df.ppi$ProtfcBinN,minMax,by.x=1,by.y=1)

df.ppi$ProtfcBinN.min <- minMax.ppi$min
df.ppi$ProtfcBinN.max <- minMax.ppi$max
df.ppi$ProtfcBinN.minMax <- factor(paste("FC: (",round(df.ppi$ProtfcBinN.min,digits = 2),",",round(df.ppi$ProtfcBinN.max,digits = 2),")",sep = ""))

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_protFCbin5_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = ProtfcBinN.minMax,y = PPI_bAct,fill=ProtfcBinN,group=ProtfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$ProtfcBinN.minMax)) +
  ggtitle(label = "PPI: all proteins",subtitle = " Binding Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,200)) 
dev.off()

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-allPPI.400_protFCbin5_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = ProtfcBinN.minMax,y = PPI_alllink,fill=ProtfcBinN,group=ProtfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$ProtfcBinN.minMax)) +
  ggtitle(label = "PPI: all proteins",subtitle = "All PPI >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,600)) 
dev.off()

ggplot(data = df.ppi,mapping = aes(y = absLog2Prot,x = PPI_bAct)) +
  # geom_boxplot(outlier.size = .4) +
  geom_point(alpha=.2) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all proteins",subtitle = " Binding Links >400") +
  xlab(label = "# PPIs _ binding ") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

head(df.ppi)
ggplot(data = df.ppi,mapping = aes(y = absLog2Prot,x = PPI_alllink)) +
  # geom_boxplot(outlier.size = .4) +
  geom_point(alpha=.2) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all proteins",subtitle = " Binding Links >400") +
  xlab(label = "# PPIs _ all links sc400") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_protFCbin_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin.prot,y = PPI_bAct,fill=fcBin.prot,group=fcBin.prot)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all proteins",subtitle = " Binding Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-allLinks400_protFCbin_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin.prot,y = PPI_alllink,fill=fcBin.prot,group=fcBin.prot)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all proteins",subtitle = "All Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

######

## Bin on RNA changes #######
## bin

## bin Protein FC #########
df.ppi$log2FoldChange[is.na(df.ppi$log2FoldChange)] <- 0
df.ppi$fcBin.RNA <- ifelse(test =  abs(df.ppi$log2FoldChange)<=.25,yes = "0-0.25",
                            no = ifelse(test = abs(df.ppi$log2FoldChange)>.25 & abs(df.ppi$log2FoldChange)<=.5,yes = "0.25-0.5",
                                        no = ifelse(test = abs(df.ppi$log2FoldChange)>.5 & abs(df.ppi$log2FoldChange)<=1,yes = "0.5-1",
                                                    no = ifelse(test = abs(df.ppi$log2FoldChange)>1 & abs(df.ppi$log2FoldChange)<=1.5,yes = "1-1.5",
                                                                no = ifelse(test = abs(df.ppi$log2FoldChange)>1.5 & abs(df.ppi$log2FoldChange)<=2,yes = "1.5-2",
                                                                            no = "2<")))))
df.ppi$fcBin.RNA

df.ppi <- df.ppi[order(abs(df.ppi$log2FoldChange)),]
head(df.ppi)

### 9 bins 
df.ppi$RNAfcBinN <- factor(c(rep(1:9,each=c(floor(3874/9))),9,9,9,9))

df.ppi$absLog2RNA <- abs(df.ppi$log2FoldChange)

minMax <- df.ppi %>% 
  select(absLog2RNA,RNAfcBinN) %>%
  group_by(RNAfcBinN) %>%
  summarise(min = min(absLog2RNA),max = max(absLog2RNA))

head(minMax)

minMax.ppi <- merge(df.ppi$RNAfcBinN,minMax,by.x=1,by.y=1)
tail(minMax.ppi)

tail(df.ppi[df.ppi$RNAfcBinN==9,])

df.ppi$RNAfcBinN.min <- minMax.ppi$min
df.ppi$RNAfcBinN.max <- minMax.ppi$max
tail(df.ppi)
df.ppi$RNAfcBinN.minMax <- factor(paste("FC: (",round(df.ppi$RNAfcBinN.min,digits = 3),",",round(df.ppi$RNAfcBinN.max,digits = 3),")",sep = ""))

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_RNAFCbin9_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = RNAfcBinN.minMax,y = PPI_bAct,fill=RNAfcBinN,group=RNAfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$RNAfcBinN.minMax)) +
  ggtitle(label = "PPI: all RNAs",subtitle = " Binding Links >400") +
  xlab(label = "RNA log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,250)) 
dev.off()

head(df.ppi)
pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-allLinks400_RNAFCbin9_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = RNAfcBinN.minMax,y = PPI_alllink,fill=RNAfcBinN,group=RNAfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$RNAfcBinN.minMax)) +
  ggtitle(label = "PPI: all RNAs",subtitle = "All PPI >400") +
  xlab(label = "RNA log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,600)) 
dev.off()


### 5 bins 
df.ppi$absLog2RNA <- abs(df.ppi$log2FoldChange)
df.ppi <- df.ppi[order(df.ppi$absLog2RNA),]
df.ppi$RNAfcBinN <- factor(c(rep(1:5,each=c(floor(3874/5))),5,5,5,5))

minMax <- df.ppi %>% 
  select(absLog2RNA,RNAfcBinN) %>%
  group_by(RNAfcBinN) %>%
  summarise(min = min(absLog2RNA),max = max(absLog2RNA))

head(minMax)

minMax.ppi <- merge(df.ppi$RNAfcBinN,minMax,by.x=1,by.y=1)
tail(minMax.ppi)

df.ppi$RNAfcBinN.min <- minMax.ppi$min
df.ppi$RNAfcBinN.max <- minMax.ppi$max
df.ppi$RNAfcBinN.minMax <- factor(paste("FC: (",round(df.ppi$RNAfcBinN.min,digits = 3),",",round(df.ppi$RNAfcBinN.max,digits = 3),")",sep = ""))
tail(df.ppi)
head(df.ppi)

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_RNAFCbin5_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = RNAfcBinN.minMax,y = PPI_bAct,fill=RNAfcBinN,group=RNAfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$RNAfcBinN.minMax)) +
  ggtitle(label = "PPI: all RNAs",subtitle = " Binding Links >400") +
  xlab(label = "RNA log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,200)) 
dev.off()

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-allPPI.400_RNAFCbin5_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = RNAfcBinN.minMax,y = PPI_alllink,fill=RNAfcBinN,group=RNAfcBinN)) +
  geom_boxplot(outlier.shape = NA) +
  #  geom_violin() +
  theme_bw() +
  scale_fill_brewer(palette = "Blues",labels=unique(df.ppi$RNAfcBinN.minMax)) +
  ggtitle(label = "PPI: all RNAs",subtitle = "All PPI >400") +
  xlab(label = "RNA log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  coord_cartesian(ylim = c(0,600)) 
dev.off()

ggplot(data = df.ppi,mapping = aes(y = absLog2RNA,x = PPI_bAct)) +
  # geom_boxplot(outlier.size = .4) +
  geom_point(alpha=.2) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all RNAs",subtitle = " Binding Links >400") +
  xlab(label = "# PPIs _ binding ") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(data = df.ppi,mapping = aes(y = absLog2RNA,x = PPI_alllink)) +
  # geom_boxplot(outlier.size = .4) +
  geom_point(alpha=.2) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all RNAs",subtitle = " Binding Links >400") +
  xlab(label = "# PPIs _ all links sc400") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_RNAFCbin_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin.RNA,y = PPI_bAct,fill=fcBin.RNA,group=fcBin.RNA)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all RNAs",subtitle = " Binding Links >400") +
  xlab(label = "RNA log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-allLinks400_RNAFCbin_allGenes.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin.RNA,y = PPI_alllink,fill=fcBin.RNA,group=fcBin.RNA)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: all RNAs",subtitle = "All Links >400") +
  xlab(label = "RNA log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

###########





#### RNA changes - protein dependence and PPIs ######

## maProtein levels - sign RNA 
limma.fc.padj <- res.eb[,c("logFC","p.mod","q.mod")]
rownames(limma.fc.padj) <- rownames(res.eb)

prot.norm.shift <- prot.norm+abs(min(prot.norm))
baseMean.prot <- rowMeans(prot.norm.shift)
df.prot.limma <- merge(baseMean.prot,limma.fc.padj,by=0)

head(limma.fc.padj)
head(df.prot.limma)

colnames(df.prot.limma) <- c("ProteinName","baseMean","log2FoldChange","padj","q.mod")
head(df.prot.limma)
df.rna <- dabs.ips.d14$test$Default[rownames(dabs.ips.d14$test$Default) %in% prot.det.mod,]
head(df.rna)
dim(df.rna)

### Make plot 
p.rna <- 0.001
p.protein <- p.rna
l <- 0

RNA.up <- df.rna[!is.na(df.rna$padj) & df.rna$padj<p.rna & df.rna$log2FoldChange>l,]
RNA.down <- df.rna[!is.na(df.rna$padj) & df.rna$padj<p.rna & df.rna$log2FoldChange<l,]
Prot.up <- df.prot.limma[!is.na(df.prot.limma$padj) & df.prot.limma$padj<p.protein & df.prot.limma$log2FoldChange>l,]
Prot.down <- df.prot.limma[!is.na(df.prot.limma$padj) & df.prot.limma$padj<p.protein & df.prot.limma$log2FoldChange<l,]

dim(RNA.up)
dim(RNA.down)
dim(Prot.up)
dim(Prot.down)

Prot_RNA.up <- merge(rownames(RNA.up),df.prot.limma,by=1)
Prot_RNA.down <- merge(rownames(RNA.down),df.prot.limma,by=1)
RNA_prot.up <- merge(Prot.up$ProteinName,df.rna,by.x=1,by.y=0)
RNA_prot.down <- merge(Prot.down$ProteinName,df.rna,by.x=1,by.y=0)

#####

## read PPI 

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
df.RNAchanged <- rbind(merge.RNA.down,merge.RNA.up)
head(df.RNAchanged)

## Binding  links #####
## Merge hgnc symbol and freq table
prot.1 <- read.delim(file = "data/StringDB/9606.binding.Acting.Links.v10.5.Freq.Prot1_score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)                 

sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
colnames(sym.freq) <- c("ENSP","PPIbAct","sym")
head(sym.freq)

df.ppi <- merge(df.RNAchanged,sym.freq[,-1],by.x=1,by.y=2,all.x=T)
head(df.ppi)
df.ppi$PPIbAct[is.na(df.ppi$PPIbAct)] <- 0
df.ppi <- df.ppi[order(df.ppi$PPIbAct),]

### all binding 400 links## Merge hgnc symbol and freq table
prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)

sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
colnames(sym.freq) <- c("ENSP","PPIall","sym")
head(sym.freq)
head(df)

df.ppi <- merge(df.ppi,sym.freq[,-1],by.x=1,by.y=2,all.x=T)
head(df.ppi)
df.ppi$PPIall[is.na(df.ppi$PPIall)] <- 0
df.ppi$PPIbAct[is.na(df.ppi$PPIbAct)] <- 0

df.ppi <- df.ppi[order(df.ppi$PPIbAct),]

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9C_maPlot_Sign.RNA_Protein_levels_PPIcol_Binding_Facet.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.ppi,mapping = aes(x = baseMean.y,y = log2FoldChange.y,col=PPIbAct,size=PPIbAct)) +
  geom_point(alpha=.8) +
  theme_classic() +
  xlab(label = "baseMean Protein") +
  ylab(label = "log2FC Protein")  +
  ggtitle(label = "Protein changes - PPI",subtitle = "Genes significant in RNA - all binding >400") +
  scale_color_continuous(low = "lightskyblue",high = "darkred") +
  scale_size_continuous(range = c(.5,6)) +
  scale_alpha_continuous(range = c(.6,1) ) +
  geom_hline(yintercept = 0,lty=3,col="grey60") +
  facet_wrap(facets = ~Dir)
dev.off()

df.ppi <- df.ppi[order(df.ppi$PPIall),]
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9C_maPlot_Sign.RNA_Protein_levels_PPIcol_AllLinks400_Facet.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.ppi,mapping = aes(x = baseMean.y,y = log2FoldChange.y,col=PPIall,size=PPIall)) +
  geom_point(alpha=.8) +
  theme_classic() +
  xlab(label = "baseMean Protein") +
  ylab(label = "log2FC Protein")  +
  ggtitle(label = "Protein changes - PPI",subtitle = "Genes significant in RNA - all links >400") +
  scale_color_continuous(low = "lightskyblue",high = "darkred") +
  scale_size_continuous(range = c(.5,6)) +
  scale_alpha_continuous(range = c(.6,1) ) +
  geom_hline(yintercept = 0,lty=3,col="grey60") +
  facet_wrap(facets = ~Dir)
dev.off()

## Bin for PPI - changes in Protein
df.ppi$fcBin <- ifelse(test = abs(df.ppi$log2FoldChange.y)<=.25,yes = "0-0.25",
                       no = ifelse(test = abs(df.ppi$log2FoldChange.y)>.25 & abs(df.ppi$log2FoldChange.y)<=.5,yes = "0.25-0.5",
                                   no = ifelse(test = abs(df.ppi$log2FoldChange.y)>.5 & abs(df.ppi$log2FoldChange.y)<=1,yes = "0.5-1",
                                               no = ifelse(test = abs(df.ppi$log2FoldChange.y)>1 & abs(df.ppi$log2FoldChange.y)<=1.5,yes = "1-1.5",
                                                           no = ifelse(test = abs(df.ppi$log2FoldChange.y)>1.5 & abs(df.ppi$log2FoldChange.y)<=2,yes = "1.5-2",
                                                                       no = ifelse(test = abs(df.ppi$log2FoldChange.y)>2 & abs(df.ppi$log2FoldChange.y)<=4,yes = "3-4",
                                                                                  no = "4<"))))))
df.ppi$fcBin

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_protFCbin_GenesChangedOnRNA.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin,y = PPIbAct,fill=fcBin,group=fcBin)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  ggtitle(label = "PPI: Genes changed on RNA",subtitle = " Binding Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-AllLinks400_protFCbin_GenesChangedOnRNA.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin,y = PPIall,fill=fcBin,group=fcBin)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  ggtitle(label = "PPI: Genes changed on RNA",subtitle = "All Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

######


#### Bin RNA changes - genes changed on Protein
merge.PROT.up <- merge(data.frame(Prot.up),RNA_prot.up,by.x=1,by.y=1)
merge.PROT.up$Dir <- "up"
merge.PROT.down <- merge(data.frame(Prot.down),RNA_prot.down,by.x=1,by.y=1)
merge.PROT.down$Dir <- "down"
merge.allPROT <- merge(df.prot.limma,data.frame(dabs.ips.d14$test$Default),by.x=1,by.y=0)
merge.allPROT$Dir <- "all"
head(merge.allPROT)
head(merge.PROT.up)
df.PROTchanged <- rbind(merge.PROT.down,merge.PROT.up)
head(df.PROTchanged)

## Binding  links #####
## Merge hgnc symbol and freq table
prot.1 <- read.delim(file = "data/StringDB/9606.binding.Acting.Links.v10.5.Freq.Prot1_score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)                 

sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
colnames(sym.freq) <- c("ENSP","PPIbAct","sym")
head(sym.freq)

df.ppi <- merge(df.PROTchanged,sym.freq[,-1],by.x=1,by.y=2,all.x=T)
head(df.ppi)
df.ppi$PPIbAct[is.na(df.ppi$PPIbAct)] <- 0
df.ppi <- df.ppi[order(df.ppi$PPIbAct),]

### all binding 400 links## Merge hgnc symbol and freq table
prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)

sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
colnames(sym.freq) <- c("ENSP","PPIall","sym")
head(sym.freq)
head(df)

df.ppi <- merge(df.ppi,sym.freq[,-1],by.x=1,by.y=2,all.x=T)
head(df.ppi)
df.ppi$PPIall[is.na(df.ppi$PPIall)] <- 0
df.ppi$PPIbAct[is.na(df.ppi$PPIbAct)] <- 0

df.ppi <- df.ppi[order(df.ppi$PPIbAct),]
head(df.ppi)

df.ppi$baseMean.y <- log2(df.ppi$baseMean.y)

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9C_maPlot_Sign.Protein_RNA_levels_PPIcol_Binding_Facet.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.ppi,mapping = aes(x = baseMean.y,y = log2FoldChange.y,col=PPIbAct,size=PPIbAct)) +
  geom_point(alpha=.8) +
  theme_classic() +
  xlab(label = "baseMean RNA") +
  ylab(label = "log2FC RNA")  +
  ggtitle(label = "RNA changes - PPI",subtitle = "Genes significant on Protein - all binding >400") +
  scale_color_continuous(low = "lightskyblue",high = "darkred") +
  scale_size_continuous(range = c(.5,6)) +
  scale_alpha_continuous(range = c(.6,1) ) +
  geom_hline(yintercept = 0,lty=3,col="grey60") +
  facet_wrap(facets = ~Dir)
dev.off()

df.ppi <- df.ppi[order(df.ppi$PPIall),]
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9C_maPlot_Sign.RNA_Protein_levels_PPIcol_AllLinks400_Facet.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.ppi,mapping = aes(x = baseMean.y,y = log2FoldChange.y,col=PPIall,size=PPIall)) +
  geom_point(alpha=.8) +
  theme_classic() +
  xlab(label = "baseMean RNA") +
  ylab(label = "log2FC RNA")  +
  ggtitle(label = "RNA changes - PPI",subtitle = "Genes significant on Protein - all links >400") +
  scale_color_continuous(low = "lightskyblue",high = "darkred") +
  scale_size_continuous(range = c(.5,6)) +
  scale_alpha_continuous(range = c(.6,1) ) +
  geom_hline(yintercept = 0,lty=3,col="grey60") +
  facet_wrap(facets = ~Dir)
dev.off()

## Bin for PPI - changes in RNA
head(df.ppi)
df.ppi$fcBin <- ifelse(test = abs(df.ppi$log2FoldChange.y)<=.25,yes = "0-0.25",
                       no = ifelse(test = abs(df.ppi$log2FoldChange.y)>.25 & abs(df.ppi$log2FoldChange.y)<=.5,yes = "0.25-0.5",
                                   no = ifelse(test = abs(df.ppi$log2FoldChange.y)>.5 & abs(df.ppi$log2FoldChange.y)<=1,yes = "0.5-1",
                                               no = ifelse(test = abs(df.ppi$log2FoldChange.y)>1 & abs(df.ppi$log2FoldChange.y)<=1.5,yes = "1-1.5",
                                                           no = ifelse(test = abs(df.ppi$log2FoldChange.y)>1.5 & abs(df.ppi$log2FoldChange.y)<=2,yes = "1.5-2",
                                                                       no = ifelse(test = abs(df.ppi$log2FoldChange.y)>2 & abs(df.ppi$log2FoldChange.y)<=4,yes = "3-4",
                                                                                   no = "4<"))))))
df.ppi$fcBin

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_RNAFCbin_GenesChangedOnProtein.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin,y = PPIbAct,fill=fcBin,group=fcBin)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  ggtitle(label = "PPI: Genes changed on Protein",subtitle = " Binding Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-AllLinks400_RNAFCbin_GenesChangedOnProtein.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin,y = PPIall,fill=fcBin,group=fcBin)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  ggtitle(label = "PPI: Genes changed on Protein",subtitle = "All Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()




















### ma plot PPI Protien sign. ####

merge.Prot.up <- merge(RNA_prot.up,data.frame(Prot.up),by=1)
merge.Prot.up$Dir <- "up"
merge.Prot.down <- merge(RNA_prot.down,data.frame(Prot.down),by=1)
merge.Prot.down$Dir <- "down"
merge.allRNA <- merge(data.frame(dabs.ips.d14$test$Default),df.prot.limma,by.x=0,by.y=1)
merge.allRNA$Dir <- "all"
colnames(merge.allRNA)[1] <- "x"
head(merge.allRNA)
head(merge.Prot.up)

df.protChange <- rbind(merge.Prot.down,merge.Prot.up)
head(df.protChange)
df.protChange$baseMean.x <- log2(df.protChange$baseMean.x)

## Merge hgnc symbol and freq table
prot.1 <- read.delim(file = "data/StringDB/9606.protein.links.v10.5.Freq.Prot1.score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)

sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
colnames(sym.freq) <- c("ENSP","PPIall","sym")
head(sym.freq)
head(df)

df.ppi <- merge(df.protChange,sym.freq[,-1],by.x=1,by.y=2,all.x=T)
head(df.ppi)
df.ppi$PPIall[is.na(df.ppi$PPIall)] <- 0

prot.1 <- read.delim(file = "data/StringDB/9606.binding..Links.v10.5.Freq.Prot1_score400.txt",sep = ",",stringsAsFactors = F,header = F)
head(prot.1)                 

sym.freq <- merge(prot.1,conversion.table,by.x=2,by.y=1)
colnames(sym.freq) <- c("ENSP","PPIbAct","sym")
head(sym.freq)

df.ppi <- merge(df.ppi,sym.freq[,-1],by.x=1,by.y=2,all.x=T)
head(df.ppi)
df.ppi$PPIbAct[is.na(df.ppi$PPIbAct)] <- 0
df.ppi <- df.ppi[order(df.ppi$PPIbAct),]

df.ppi$PPIbAct[is.na(df.ppi$PPIbAct)] <- 0

df.ppi <- df.ppi[order(df.ppi$PPIbAct),]

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9C_maPlot_Sign.Protein_RNA_levels_PPIcol_binding400_Facet.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.ppi,mapping = aes(x = baseMean.x,y = log2FoldChange.x,color=PPIbAct,size=PPIbAct,alpha=PPIbAct)) +
  geom_point() +
  theme_classic() +
  xlab(label = "baseMean RNA") +
  ylab(label = "log2FC RNA")  +
  ggtitle(label = "RNA changes - PPI",subtitle = "Genes significant in Protein - binding >400") +
  scale_color_continuous(low = "lightskyblue",high = "darkred") +
  scale_size_continuous(range = c(.2,6)) +
  scale_alpha_continuous(range = c(.4,1) ) +
  geom_hline(yintercept = 0,lty=3,col="grey60") +
  facet_wrap(facets = ~Dir)
dev.off()

df.ppi <- df.ppi[order(df.ppi$PPIall),]
filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9C_maPlot_Sign.Protein_RNA_levels_PPIcol_AllLinks400_Facet.pdf"
pdf(file = filename,width = 6,height = 6)
ggplot(data = df.ppi,mapping = aes(x = baseMean.x,y = log2FoldChange.x,color=PPIall,size=PPIall,alpha=PPIall)) +
  geom_point() +
  theme_classic() +
  xlab(label = "baseMean RNA") +
  ylab(label = "log2FC RNA")  +
  ggtitle(label = "RNA changes - PPI",subtitle = "Genes significant in Protein - bindingAll >400") +
  scale_color_continuous(low = "lightskyblue",high = "darkred") +
  scale_size_continuous(range = c(.2,6)) +
  scale_alpha_continuous(range = c(.4,1) ) +
  geom_hline(yintercept = 0,lty=3,col="grey60") +
  facet_wrap(facets = ~Dir)
dev.off()

## Bind for PPI
df.ppi$fcBin <- ifelse(test = abs(df.ppi$log2FoldChange.y)<=.5,yes = "0-0.5",
                       no = ifelse(test = abs(df.ppi$log2FoldChange.y)>.5 & abs(df.ppi$log2FoldChange.y)<=1,yes = "0.5-1",
                                   no = ifelse(test = abs(df.ppi$log2FoldChange.y)>1 & abs(df.ppi$log2FoldChange.y)<=1.5,yes = "1-1.5",
                                               no = ifelse(test = abs(df.ppi$log2FoldChange.y)>1.5 & abs(df.ppi$log2FoldChange.y)<=2,yes = "1.5-2",
                                                           no = ifelse(test = abs(df.ppi$log2FoldChange.y)>2 & abs(df.ppi$log2FoldChange.y)<=3,yes = "2-3",
                                                                       no = ifelse(test = abs(df.ppi$log2FoldChange.y)>3 & abs(df.ppi$log2FoldChange.y)<=4,yes = "3-4",
                                                                                   no = "4<"))))))
df.ppi$fcBin

pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-binding400_protFCbin_GenesChangedOnProtein.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin,y = PPIbAct,fill=fcBin,group=fcBin)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: Genes changed on Protein",subtitle = " Binding Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()


pdf(file = "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig9-StringDB/Fig9D_Boxplot_PPIint-AllLinks400_protFCbin_GenesChangedOnProtein.pdf",width = 4,height = 4)
ggplot(data = df.ppi,mapping = aes(x = fcBin,y = PPIall,fill=fcBin,group=fcBin)) +
  geom_boxplot(outlier.size = .4) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  ggtitle(label = "PPI: Genes changed on Protein",subtitle = "All Links >400") +
  xlab(label = "Protein log2 Fold-Change bin") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()


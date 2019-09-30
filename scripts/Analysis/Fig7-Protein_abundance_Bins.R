## High med low abundance protein #####

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

### Plot ma Protein

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

df.protdet <- df.protdet[order(-df.protdet$baseMean),]
df <- data.frame(baseMean = df.protdet$baseMean,log2FC=df.protdet$log2FoldChange)
df$bin <- c(rep("1-high",(nrow(df.protdet)/4+1)),rep("2-mid-high",(nrow(df.protdet)/4+1)),rep("3-mid-low",(nrow(df.protdet)/4+1)),rep("4-low",nrow(df.protdet)/4))

vline.val <- df$baseMean[which(!duplicated(df$bin))]
maxLim <- max(c(abs(df$log2FC[!is.na(df$log2FC)])))

filename <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig7-ProtAbundanceBin/Fig7_Protein_maPlot_bins.pdf"
pdf(file = filename,width = 6,height = 5)
pProt <- ggplot(data = df,mapping = aes(x = baseMean,y = log2FC,color=bin)) +
  geom_point(alpha=.8,cex=.4) +
  theme_classic() +
  ggtitle(label = "Protein bins",subtitle = paste("Four bins of ",ceiling((nrow(df.protdet)/4+1))," genes each",sep="")) +
  guides(colour = guide_legend(override.aes = list(size=3,alpha=.9))) +
  geom_hline(yintercept = 0,lty=2,col="grey50") +
  scale_color_manual(values = brewer.pal(n = 4, name = "RdBu")) +
  geom_vline(xintercept = vline.val[2],lty=3,col="grey50") +
  geom_vline(xintercept = vline.val[3],lty=3,col="grey50") +
  geom_vline(xintercept = vline.val[4],lty=3,col="grey50")
print(pProt)
dev.off()

#######

### Plot FC FC corr in bins (PRotien bins) #######

protFDR <- 0.001
fc.prot <- 0
fc.rna <- 0
rna.padj <- 0.001

RNAlimit <- NA
protLimit <- NA

head(combine)
plotFCcorr = function(protLimit=NA,RNAlimit=NA,maLimit=6) {
  
  if(!is.na(protLimit)){
    combine <- combine[order(-combine$y),]
    combine.exp <- combine[protLimit[1]:protLimit[2],]
  }
  
  if(!is.na(RNAlimit)){
    combine.exp <- combine.exp[combine.exp$baseMean>RNAlimit,]
  }
  
  if(is.na(RNAlimit) && is.na(protLimit)){
    combine.exp <- combine[combine$baseMean>baseMean,]
  }
  
  head(combine.exp)
  
  # 3. color data according to significance ######
  col <- ifelse(test = combine.exp$p.mod<protFDR & abs(combine.exp$logFC)>fc.prot,
                # if prot is significant
                yes = ifelse(test = combine.exp$padj<rna.padj & abs(combine.exp$log2FoldChange)>fc.rna,
                             # prot sign and ALSO rna 
                             yes = ifelse(test = c(combine.exp$logFC>0)==c(combine.exp$log2FoldChange>0),
                                          # Both sign AND same direction: GREEN
                                          yes = "forestgreen",
                                          # Both sign NOT same direction: RED
                                          no = "red"),
                             # prot sign, but NOT RNA: BLUE
                             no = "skyblue3"),
                # Prot is not sign
                no = ifelse(test = combine.exp$padj<rna.padj & abs(combine.exp$log2FoldChange)>fc.rna,
                            # Prot NOT sign, RNA sign : 
                            yes = "goldenrod2",
                            # NO sign
                            no = "black"))
  cex <- ifelse(col=="black",yes = .4,no = .8)
  table(col)
  
  
  #  4. Plot FC - FC correlation + color sign. #####
  #filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Final Figure/Figures/Figure3/Fig3A_ProtLim.",protLimit[1],"-",protLimit[2],"_RNAlim.",RNAlimit,"_FCcorr_RNA_maPlot_padj.",rna.padj,"_log2fc_",fc.prot,".pdf",sep = "")
  #print(filename)
  #pdf(file = filename)
  
  c1 <- combine.exp$logFC
  c2 <- combine.exp$log2FoldChange
  pch <- 16
  
  if(!is.na(maLimit)){
    c1[c1>=maLimit] <- maLimit
    c2[c2>=maLimit] <- maLimit
    c1[c1<=-maLimit] <- -maLimit
    c2[c2<=-maLimit] <- -maLimit
    
    pch <- ifelse(test = abs(c1)==maLimit | abs(c2)==maLimit,yes = 2,no = 16)
  }
  
  plot(c1,c2,pch=pch,cex=cex,col=col,ylab="RNA: log2fc(human/chimp)",xlab="PROTEIN: log2fc(human/chimp)",ylim=c(-maLimit,maLimit),xlim=c(-maLimit,maLimit))
  cor <- cor.test( x = c1, y = c2,method = "pearson",conf.level = 0.95)
  cor$estimate
  fit<-lm(c1~c2)
  mtext(text = paste("BaseMean: ",RNAlimit,", ProteinPadj: ",protFDR,", RNAseq padj: ",rna.padj," log2FC: RNA_",fc.rna,"Prot_",fc.prot,sep = ""),cex = .5)
  mtext(text = paste("proteinRange: ",protLimit[1],"-",protLimit[2],sep = ""),cex = 1,line = 1)
  grid(nx = 10,ny = 10)
  abline(fit,lty=2,col="red",lwd=1)
  legend("topleft",bty='n',legend=paste("Pearson's: ",format(cor$estimate,digits = 3),sep = ""))
  colLeg <- table(col)[c("black","skyblue3","forestgreen","goldenrod2","red")]
 # legend("topright",legend = paste(c("non-signficant","Protein-only","Both-sameChange","RNA-only","Both-oppositeChange"),": ",colLeg,sep = ""),col = names(colLeg),pch=16,cex=.9)
  points(combine.exp$logFC,combine.exp$log2FoldChange,col=col,pch=16,cex=cex)
  points(combine.exp$logFC,combine.exp$log2FoldChange,col=col,pch=16,cex=ifelse(col=="forestgreen",yes = .4,no = 0))
  points(combine.exp$logFC,combine.exp$log2FoldChange,col=col,pch=16,cex=ifelse(col=="red",yes = .4,no = 0))
  
  colLeg <- colLeg[!is.na(names(colLeg))]
  pie(x = colLeg,col = names(colLeg),labels = "")
  legend("topright",legend = paste(c("non-signficant","Protein-only","Both-sameChange","RNA-only","Both-oppositeChange"),": ",colLeg,sep = ""),col = names(colLeg),pch=16,cex=1)
  
}


### Plot based on Protein abundance bins

plotDir <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig7-ProtAbundanceBin/"
filename <- paste(plotDir,"/Fig7_FCcorr_div4bin_prot.pdf",sep = "")
RNAlimit <- NA

pdf(file = filename,width = 10,height = 20)
par(mfrow=c(4,2))
protLimit <- c(1,1239)
plotFCcorr(protLimit = protLimit,RNAlimit = RNAlimit,maLimit = 6)

protLimit <- c(1239+1,1240+1239)
plotFCcorr(protLimit = protLimit,RNAlimit = RNAlimit,maLimit = 6)

protLimit <- c(1240+1239+1,1240+1239+1+1239)
plotFCcorr(protLimit = protLimit,RNAlimit = RNAlimit,maLimit = 6)

protLimit <- c(1240+1239+1+1239,1240+1239+1+1239+1239)
plotFCcorr(protLimit = protLimit,RNAlimit = RNAlimit,maLimit = 6)
dev.off()

#########

### Bin RNA
head(combine)
plotFCcorr.RNAbin = function(protLimit=NA,RNAlimit=NA,maLimit=10) {
  
  if(!is.na(protLimit)){
    combine.exp <- combine[combine$y>protLimit,]
  }
  
  if(!is.na(RNAlimit)){
    combine <- combine[order(-combine$baseMean),]  
    combine.exp <- combine[RNAlimit[1]:RNAlimit[2],]
  }
  
  if(is.na(RNAlimit) && is.na(protLimit)){
    combine.exp <- combine[combine$baseMean>baseMean,]
  }
  
  head(combine.exp)
  
  # 3. color data according to significance ######
  col <- ifelse(test = combine.exp$p.mod<protFDR & abs(combine.exp$logFC)>fc.prot,
                # if prot is significant
                yes = ifelse(test = combine.exp$padj<rna.padj & abs(combine.exp$log2FoldChange)>fc.rna,
                             # prot sign and ALSO rna 
                             yes = ifelse(test = c(combine.exp$logFC>0)==c(combine.exp$log2FoldChange>0),
                                          # Both sign AND same direction: GREEN
                                          yes = "forestgreen",
                                          # Both sign NOT same direction: RED
                                          no = "red"),
                             # prot sign, but NOT RNA: BLUE
                             no = "skyblue3"),
                # Prot is not sign
                no = ifelse(test = combine.exp$padj<rna.padj & abs(combine.exp$log2FoldChange)>fc.rna,
                            # Prot NOT sign, RNA sign : 
                            yes = "goldenrod2",
                            # NO sign
                            no = "black"))
  cex <- ifelse(col=="black",yes = .4,no = .8)
  table(col)
  
  
  #  4. Plot FC - FC correlation + color sign. #####
  #filename <- paste("Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Final Figure/Figures/Figure3/Fig3A_ProtLim.",protLimit[1],"-",protLimit[2],"_RNAlim.",RNAlimit,"_FCcorr_RNA_maPlot_padj.",rna.padj,"_log2fc_",fc.prot,".pdf",sep = "")
  #print(filename)
  #pdf(file = filename)
  
  c1 <- combine.exp$logFC
  c2 <- combine.exp$log2FoldChange
  pch <- 16
  
  if(!is.na(maLimit)){
    c1[c1>=maLimit] <- maLimit
    c2[c2>=maLimit] <- maLimit
    c1[c1<=-maLimit] <- -maLimit
    c2[c2<=-maLimit] <- -maLimit
    
    pch <- ifelse(test = abs(c1)==maLimit | abs(c2)==maLimit,yes = 2,no = 16)
  }
  
  plot(c1,c2,pch=pch,cex=cex,col=col,ylab="RNA: log2fc(human/chimp)",xlab="PROTEIN: log2fc(human/chimp)",ylim=c(-maLimit,maLimit),xlim=c(-maLimit,maLimit))
  cor <- cor.test( x = c1, y = c2,method = "pearson",conf.level = 0.95)
  cor$estimate
  fit<-lm(c1~c2)
  mtext(text = paste("ProteinPadj: ",protFDR,", RNAseq padj: ",rna.padj," log2FC: RNA_",fc.rna,"Prot_",fc.prot,sep = ""),cex = .5)
  mtext(text = paste("RNA-range: ",RNAlimit[1],"-",RNAlimit[2],sep = ""),cex = 1,line = 1)
  grid(nx = 10,ny = 10)
  abline(fit,lty=2,col="red",lwd=1)
  legend("topleft",bty='n',legend=paste("Pearson's: ",format(cor$estimate,digits = 3),sep = ""))
  colLeg <- table(col)[c("black","skyblue3","forestgreen","goldenrod2","red")]
  # legend("topright",legend = paste(c("non-signficant","Protein-only","Both-sameChange","RNA-only","Both-oppositeChange"),": ",colLeg,sep = ""),col = names(colLeg),pch=16,cex=.9)
  points(combine.exp$logFC,combine.exp$log2FoldChange,col=col,pch=16,cex=cex)
  points(combine.exp$logFC,combine.exp$log2FoldChange,col=col,pch=16,cex=ifelse(col=="forestgreen",yes = .4,no = 0))
  points(combine.exp$logFC,combine.exp$log2FoldChange,col=col,pch=16,cex=ifelse(col=="red",yes = .4,no = 0))

  colLeg <- colLeg[!is.na(names(colLeg))]
  pie(x = colLeg,col = names(colLeg),labels = "")
  legend("topright",legend = paste(c("non-signficant","Protein-only","Both-sameChange","RNA-only","Both-oppositeChange"),": ",colLeg,sep = ""),col = names(colLeg),pch=16,cex=1)
  
}


### Plot based on Protein abundance bins

plotDir <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig7-ProtAbundanceBin/"
filename <- paste(plotDir,"/Fig7_FCcorr_div4bin_RNA.pdf",sep = "")
protLimit <- NA

pdf(file = filename,width = 10,height = 20)
par(mfrow=c(4,2))
RNAlimit <- c(1,1239)
plotFCcorr.RNAbin(RNAlimit = RNAlimit,protLimit = NA,maLimit = 6)

RNAlimit <- c(1239+1,1240+1239)
plotFCcorr.RNAbin(RNAlimit = RNAlimit,protLimit = NA,maLimit = 6)

RNAlimit <- c(1240+1239+1,1240+1239+1+1239)
plotFCcorr.RNAbin(RNAlimit = RNAlimit,protLimit = NA,maLimit = 6)

RNAlimit <- c(1240+1239+1+1239,1240+1239+1+1239+1239)
plotFCcorr.RNAbin(RNAlimit = RNAlimit,protLimit = NA,maLimit = 6)

dev.off()

#####


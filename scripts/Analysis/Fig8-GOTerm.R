### GO term analysis ######

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

## Read GO Terms #####

data.dir <- "Plots_manuscript_scripts_figs/Merged_pt6_hg38_counts/Figures Final/Fig8-GOterm/GO/"

## read the plain stats - # genes in each term
prot.BP <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/BP.txt",sep=""),header=F)
prot.MF <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/MF.txt",sep=""),header=F)
prot.CC <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/CC.txt",sep=""),header=F)
prot.ProtC <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/ProteinClass.txt",sep=""),header=F)
prot.Path <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/Pathway.txt",sep=""),header=F)

RNA.BP <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/BP.txt",sep=""),header=F)
RNA.MF <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/MF.txt",sep=""),header=F)
RNA.CC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/CC.txt",sep=""),header=F)
RNA.ProtC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/ProteinClass.txt",sep=""),header=F)
RNA.Path <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/Pathway.txt",sep=""),header=F)

both.BP <- read.delim(file = paste(data.dir,"/Changed.both/Changed/BP.txt",sep=""),header=F)
both.MF <- read.delim(file = paste(data.dir,"/Changed.both/Changed/MF.txt",sep=""),header=F)
both.CC <- read.delim(file = paste(data.dir,"/Changed.both/Changed/CC.txt",sep=""),header=F)
both.ProtC <- read.delim(file = paste(data.dir,"/Changed.both/Changed/ProteinClass.txt",sep=""),header=F)
both.Path <- read.delim(file = paste(data.dir,"/Changed.both/Changed/Pathway.txt",sep=""),header=F)


## Read the tests - up / down 
protUp.BP <- read.delim(file = paste(data.dir,"/Changed.protein.only/Protup/ProtUp_SBP.txt",sep=""),header=F,skip=12)
protUp.MF <- NULL
protUp.CC <- read.delim(file = paste(data.dir,"/Changed.protein.only/Protup/ProtUp_CC.txt",sep=""),header=F,skip=12)
protUp.ProtC <- read.delim(file = paste(data.dir,"/Changed.protein.only/Protup/ProtUp_ProtC.txt",sep=""),header=F,skip=12)
protUp.Path <- read.delim(file = paste(data.dir,"/Changed.protein.only/Protup/ProtUp_pathways.txt",sep=""),header=F,skip=12)

protDown.BP <- read.delim(file = paste(data.dir,"/Changed.protein.only/ProtDown/ProtDown_SBP.txt",sep=""),header=F,skip=12)
protDown.MF <- read.delim(file = paste(data.dir,"/Changed.protein.only/ProtDown/ProtDown_SMF.txt",sep=""),header=F,skip=12)
protDown.CC <- read.delim(file = paste(data.dir,"/Changed.protein.only/ProtDown/ProtDown_CC.txt",sep=""),header=F,skip=12)
protDown.ProtC <- read.delim(file = paste(data.dir,"/Changed.protein.only/ProtDown/ProtDown_ProtC.txt",sep=""),header=F,skip=12)
protDown.Path <- read.delim(file = paste(data.dir,"/Changed.protein.only/ProtDown/ProtDown_pathways.txt",sep=""),header=F,skip=12)

RNAup.BP <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAup/RNAup_SBP.txt",sep=""),header=F,skip=12)
RNAup.MF <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAup/RNAup_SMF.txt",sep=""),header=F,skip=12)
RNAup.CC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAup/RNAup_CC.txt",sep=""),header=F,skip=12)
RNAup.ProtC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAup/RNAup_ProtC.txt",sep=""),header=F,skip=12)
RNAup.Path <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAup/RNAup_pathways.txt",sep=""),header=F,skip=12)

RNAdown.BP <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAdown/RNAdown_SBP.txt",sep=""),header=F,skip=12)
RNAdown.MF <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAdown/RNAdown_SMF.txt",sep=""),header=F,skip=12)
RNAdown.CC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAdown/RNAdown_CC.txt",sep=""),header=F,skip=12)
RNAdown.ProtC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAdown/RNAdown_ProtC.txt",sep=""),header=F,skip=12)
RNAdown.Path <- NULL #read.delim(file = paste(data.dir,"/Changed.RNA.only/RNAdown/RNAdown_pathways.txt",sep=""),header=F,skip=12)

BothUp.BP <- read.delim(file = paste(data.dir,"/Changed.both/BothUp/BothUp_SBP.txt",sep=""),header=F,skip=12)
BothUp.MF <- read.delim(file = paste(data.dir,"/Changed.both/BothUp/BothUp_SMF.txt",sep=""),header=F,skip=12)
BothUp.CC <- NULL# read.delim(file = paste(data.dir,"/Changed.both/BothUp/BothUp_CC.txt",sep=""),header=F,skip=12)
BothUp.ProtC <- read.delim(file = paste(data.dir,"/Changed.both/BothUp/BothUp_ProtC.txt",sep=""),header=F,skip=12)
BothUp.Path <- read.delim(file = paste(data.dir,"/Changed.both/BothUp/BothUp_pathways.txt",sep=""),header=F,skip=12)

BothDown.BP <- read.delim(file = paste(data.dir,"/Changed.both/BothDown/BothDown_SBP.txt",sep=""),header=F,skip=12)
BothDown.MF <- read.delim(file = paste(data.dir,"/Changed.both/BothDown/BothDown_SMF.txt",sep=""),header=F,skip=12)
BothDown.CC <- read.delim(file = paste(data.dir,"/Changed.both/BothDown/BothDown_CC.txt",sep=""),header=F,skip=12)
BothDown.ProtC <- read.delim(file = paste(data.dir,"/Changed.both/BothDown/BothDown_ProtC.txt",sep=""),header=F,skip=12)
BothDown.Path <- read.delim(file = paste(data.dir,"/Changed.both/BothDown/BothDown_pathways.txt",sep=""),header=F,skip=12)


## Read the tests - up / down fused
prot.f.BP <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/TestEnrichment/Prot_BP.txt",sep=""),header=F,skip=12)
prot.f.MF <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/TestEnrichment/Prot_BP.txt",sep=""),header=F,skip=12)
prot.f.CC <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/TestEnrichment/Prot_CC.txt",sep=""),header=F,skip=12)
prot.f.ProtC <- read.delim(file = paste(data.dir,"/Changed.protein.only/Changed/TestEnrichment/Prot_ProtC.txt",sep=""),header=F,skip=12)
prot.f.Path <- NULL

RNA.f.BP <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/TestEnrichment/RNA_BP.txt",sep=""),header=F,skip=12)
RNA.f.MF <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/TestEnrichment/RNA_BP.txt",sep=""),header=F,skip=12)
RNA.f.CC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/TestEnrichment/RNA_CC.txt",sep=""),header=F,skip=12)
RNA.f.ProtC <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/TestEnrichment/RNA_ProtC.txt",sep=""),header=F,skip=12)
RNA.f.Path <- read.delim(file = paste(data.dir,"/Changed.RNA.only/Changed/TestEnrichment/RNA_Path.txt",sep=""),header=F,skip=12)


both.f.BP <- read.delim(file = paste(data.dir,"/Changed.both/Changed/TestEnrichment/Both_BP.txt",sep=""),header=F,skip=12)
both.f.MF <- read.delim(file = paste(data.dir,"/Changed.both/Changed/TestEnrichment/Both_BP.txt",sep=""),header=F,skip=12)
both.f.CC <- read.delim(file = paste(data.dir,"/Changed.both/Changed/TestEnrichment/Both_CC.txt",sep=""),header=F,skip=12)
both.f.ProtC <- read.delim(file = paste(data.dir,"/Changed.both/Changed/TestEnrichment/Both_ProtC.txt",sep=""),header=F,skip=12)
both.f.Path <- read.delim(file = paste(data.dir,"/Changed.both/Changed/TestEnrichment/Both_Path.txt",sep=""),header=F,skip=12)

########## 

## Plot GO Terms #####
library(ggplot2)
library(wesanderson)

### Plot enrichment tests up/down separate ######

##### C C #######
term <- "Cellular Compartment (Slim)"
p.up <- protUp.CC
p.down <- protDown.CC
r.up <- RNAup.CC
r.down <- RNAdown.CC
b.up <- BothUp.CC
b.down <- BothDown.CC

p.up$Group <- "Prot only"
p.down$Group <- "Prot only"
r.up$Group <- "RNA only"
r.down$Group <- "RNA only"
b.up$Group <- "Both"
b.down$Group <- "Both"

p.up$Dir <- "Up"
p.down$Dir <- "Down"
r.up$Dir <- "Up"
r.down$Dir <- "Down"
b.up$Dir <- "Up"
b.down$Dir <- "Down"

df.cc <- rbind(r.up,r.down,p.up,p.down,b.down)
df.cc$V1 <- sub(" \\(.*","",df.cc$V1)
df.cc$V6 <- log2(df.cc$V6)
df.cc$V1 <- factor(x = df.cc$V1,levels = unique(df.cc$V1))

p.CC <- ggplot(data = df.cc,mapping = aes(x=reorder(V1,Group),y=V6,fill=Dir)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("burlywood3","darkolivegreen4")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",position=position_dodge(width=1.5)) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.CC

filename <- paste(data.dir,"/Plots/Up.down.sep_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 4)
print(p.CC)
dev.off()

#######

##### B P  #######
term <- "Biological Process (Slim)"

p.up <- protUp.BP
p.down <- protDown.BP
r.up <- RNAup.BP
r.down <- RNAdown.BP
b.up <- BothUp.BP
b.down <- BothDown.BP

p.up$Group <- "Prot only"
p.down$Group <- "Prot only"
r.up$Group <- "RNA only"
r.down$Group <- "RNA only"
b.up$Group <- "Both"
b.down$Group <- "Both"

p.up$Dir <- "Up"
p.down$Dir <- "Down"
r.up$Dir <- "Up"
r.down$Dir <- "Down"
b.up$Dir <- "Up"
b.down$Dir <- "Down"

df.BP <- rbind(r.up,r.down,p.up,p.down,b.up,b.down)
df.BP$V1 <- sub(" \\(.*","",df.BP$V1)
df.BP$V6 <- log2(df.BP$V6)
df.BP$V1 <- factor(df.BP$V1,levels = unique(df.BP$V1))

p.BP <- ggplot(data = df.BP,mapping = aes(x=reorder(V1,Group),y=V6,fill=Dir)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("burlywood3","darkolivegreen4")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",position=position_dodge(width=1.5)) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.BP

filename <- paste(data.dir,"/Plots/Up.down.sep_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 7)
print(p.BP)
dev.off()

#####

##### M F  #######
term <- "Molecular Function (Slim)"

p.up <- protUp.MF
p.down <- protDown.MF
r.up <- RNAup.MF
r.down <- RNAdown.MF
b.up <- BothUp.MF
b.down <- BothDown.MF

p.up$Group <- "Prot only"
p.down$Group <- "Prot only"
r.up$Group <- "RNA only"
r.down$Group <- "RNA only"
b.up$Group <- "Both"
b.down$Group <- "Both"

p.up$Dir <- "Up"
p.down$Dir <- "Down"
r.up$Dir <- "Up"
r.down$Dir <- "Down"
b.up$Dir <- "Up"
b.down$Dir <- "Down"

df.MF <- rbind(r.up,r.down,p.down,b.up,b.down)
df.MF$V1 <- sub(" \\(.*","",df.MF$V1)
df.MF$V6 <- log2(df.MF$V6)
df.MF$V1 <- factor(df.MF$V1,levels = unique(df.MF$V1))

p.MF <- ggplot(data = df.MF,mapping = aes(x=reorder(V1,Group),y=V6,fill=Dir)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("burlywood3","darkolivegreen4")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",position=position_dodge(width=1.5)) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.MF

filename <- paste(data.dir,"/Plots/Up.down.sep_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 4)
print(p.MF)
dev.off()

######

##### path #######
term <- "Panther Pathway"

p.up <- protUp.Path
p.down <- protDown.Path
r.up <- RNAup.Path
r.down <- RNAdown.Path
b.up <- BothUp.Path
b.down <- BothDown.Path

p.up$Group <- "Prot only"
p.down$Group <- "Prot only"
r.up$Group <- "RNA only"
r.down$Group <- "RNA only"
b.up$Group <- "Both"
b.down$Group <- "Both"

p.up$Dir <- "Up"
p.down$Dir <- "Down"
r.up$Dir <- "Up"
r.down$Dir <- "Down"
b.up$Dir <- "Up"
b.down$Dir <- "Down"

df.Pathway <- rbind(r.up,p.down,b.up,b.down)
df.Pathway$V1 <- sub(" \\(.*","",df.Pathway$V1)
df.Pathway$V6 <- log2(df.Pathway$V6)
df.Pathway$V1 <- factor(df.Pathway$V1,levels = unique(df.Pathway$V1))

p.Pathway <- ggplot(data = df.Pathway,mapping = aes(x=reorder(V1,Group),y=V6,fill=Dir)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("burlywood3","darkolivegreen4")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",position=position_dodge(width=1.5)) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.Pathway

filename <- paste(data.dir,"/Plots/Up.down.sep_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 4)
print(p.Pathway)
dev.off()

######


##### ProtC #######
term <- "Protein Class"

p.up <- protUp.ProtC
p.down <- protDown.ProtC
r.up <- RNAup.ProtC
r.down <- RNAdown.ProtC
b.up <- BothUp.ProtC
b.down <- BothDown.ProtC

p.up$Group <- "Prot only"
p.down$Group <- "Prot only"
r.up$Group <- "RNA only"
r.down$Group <- "RNA only"
b.up$Group <- "Both"
b.down$Group <- "Both"

p.up$Dir <- "Up"
p.down$Dir <- "Down"
r.up$Dir <- "Up"
r.down$Dir <- "Down"
b.up$Dir <- "Up"
b.down$Dir <- "Down"

df.ProtC <- rbind(r.up,p.down,b.up,b.down)
df.ProtC$V1 <- sub(" \\(.*","",df.ProtC$V1)
df.ProtC$V6 <- log2(df.ProtC$V6)
df.ProtC$V1 <- factor(df.ProtC$V1,levels = unique(df.ProtC$V1))

p.ProtC <- ggplot(data = df.ProtC,mapping = aes(x=reorder(V1,Group),y=V6,fill=Dir)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("burlywood3","darkolivegreen4")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",position=position_dodge(width=1.5)) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.ProtC

filename <- paste(data.dir,"/Plots/Up.down.sep_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 4.5)
print(p.ProtC)
dev.off()

######
###################################

### Plot enrichment tests up/down combined ######

##### C C #######
term <- "Cellular Compartment (Slim)"

p <- prot.f.CC
r <- RNA.f.CC
b <- both.f.CC

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.cc <- rbind(r,p,b)
df.cc$V1 <- sub(" \\(.*","",df.cc$V1)
df.cc$V6 <- log2(df.cc$V6)
df.cc$V1 <- factor(x = df.cc$V1,levels = unique(df.cc$V1))

df.cc$fold <- ifelse(df.cc$V6>0,yes = "Pos",no = "Neg")
df.cc <- df.cc[as.character(df.cc$V1)!="Unclassified",]
p.CC <- ggplot(data = df.cc,mapping = aes(x=reorder(V1,Group),y=V6,fill=fold)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("Firebrick3","Darkblue")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",hjust=-.5) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.CC

filename <- paste(data.dir,"/Plots/UPdownFused_",term,".pdf",sep = "")
pdf(file = filename,width = 6,height = 4)
print(p.CC)
dev.off()

#######

##### B P  #######
term <- "Biological Process (Slim)"

p <- prot.f.BP
r <- RNA.f.BP
b <- both.f.BP

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.BP <- rbind(r,p,b)
df.BP$V1 <- sub(" \\(.*","",df.BP$V1)
df.BP$V6 <- log2(df.BP$V6)
df.BP$V1 <- factor(x = df.BP$V1,levels = unique(df.BP$V1))

df.BP$fold <- ifelse(df.BP$V6>0,yes = "Pos",no = "Neg")
df.BP <- df.BP[as.character(df.BP$V1)!="Unclassified",]
p.BP <- ggplot(data = df.BP,mapping = aes(x=reorder(V1,Group),y=V6,fill=fold)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("Firebrick3","Darkblue")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",hjust=-.5) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.BP

filename <- paste(data.dir,"/Plots/UPdownFused_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 6)
print(p.BP)
dev.off()
#####

##### M F  #######
term <- "Molecular Function (Slim)"

p <- prot.f.MF
r <- RNA.f.MF
b <- both.f.MF

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.MF <- rbind(r,p,b)
df.MF$V1 <- sub(" \\(.*","",df.MF$V1)
df.MF$V6 <- log2(df.MF$V6)
df.MF$V1 <- factor(x = df.MF$V1,levels = unique(df.MF$V1))

df.MF$fold <- ifelse(df.MF$V6>0,yes = "Pos",no = "Neg")
df.MF <- df.MF[as.character(df.MF$V1)!="Unclassified",]
p.MF <- ggplot(data = df.MF,mapping = aes(x=reorder(V1,Group),y=V6,fill=fold)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("Firebrick3","Darkblue")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",hjust=-.5) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.MF

filename <- paste(data.dir,"/Plots/UPdownFused_",term,".pdf",sep = "")
pdf(file = filename,width = 6,height = 6)
print(p.MF)
dev.off()

######

##### path #######
term <- "Panther Pathway"

p <- prot.f.Path
r <- RNA.f.Path
b <- both.f.Path

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.Path <- rbind(r,b)
df.Path$V1 <- sub(" \\(.*","",df.Path$V1)
df.Path$V6 <- log2(df.Path$V6)
df.Path$V1 <- factor(x = df.Path$V1,levels = unique(df.Path$V1))

df.Path$fold <- ifelse(df.Path$V6>0,yes = "Pos",no = "Neg")
df.Path <- df.Path[as.character(df.Path$V1)!="Unclassified",]

p.Path <- ggplot(data = df.Path,mapping = aes(x=reorder(V1,Group),y=V6,fill=fold)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("Darkblue")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",hjust=-.5) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.Path

filename <- paste(data.dir,"/Plots/UPdownFused_",term,".pdf",sep = "")
pdf(file = filename,width = 6,height = 3)
print(p.Path)
dev.off()

######


##### ProtC #######
term <- "Protein Class"

p <- prot.f.ProtC
r <- RNA.f.ProtC
b <- both.f.ProtC

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.ProtC <- rbind(r,p,b)
df.ProtC$V1 <- sub(" \\(.*","",df.ProtC$V1)
df.ProtC$V6 <- log2(df.ProtC$V6)
df.ProtC$V1 <- factor(x = df.ProtC$V1,levels = unique(df.ProtC$V1))

df.ProtC$fold <- ifelse(df.ProtC$V6>0,yes = "Pos",no = "Neg")
df.ProtC <- df.ProtC[as.character(df.ProtC$V1)!="Unclassified",]

p.ProtC <- ggplot(data = df.ProtC,mapping = aes(x=reorder(V1,Group),y=V6,fill=fold)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") + 
  coord_flip() +
  ggtitle(label = term) +
  ylab(label = "log2(Fold enrichment)") +
  xlab(label = "") +
  theme_classic() +
  scale_fill_manual(values = c("Firebrick3","Darkblue")) +
  geom_hline(yintercept = 0,lty=2) +
  facet_wrap(~Group) +
  geom_text(aes(label=V3),size=3,color="black",hjust=-.5) +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank())
p.ProtC

filename <- paste(data.dir,"/Plots/UPdownFused_",term,".pdf",sep = "")
pdf(file = filename,width = 6,height = 6)
print(p.ProtC)
dev.off()

######
##################
### Plot dry stats ######

## CC ####
term <- "Cellular Compartment (Slim)"
p <- prot.CC
r <- RNA.CC
b <- both.CC

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.cc <- rbind(r,p,b)
df.cc$V2 <- sub(" \\(.*","",df.cc$V2)
df.cc$V4 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.cc$V4)))
df.cc$V5 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.cc$V5)))
df.cc$V2 <- factor(x = df.cc$V2,levels = unique(df.cc$V2))

p.CC <- ggplot(data = df.cc,mapping = aes(x=Group,y=V5,fill=Group)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") +
  scale_fill_manual(values = c("forestgreen","skyblue3","goldenrod2")) +
  ggtitle(label = term) +
  facet_wrap(~V2,scales = "free_y") +
  xlab(label = "") +
  ylab(label = "% of genes") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p.CC

filename <- paste(data.dir,"/Plots/Changed_PercentageTerms_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 8)
print(p.CC)
dev.off()

#######

##### B P  #######
term <- "Biological Process (Slim)"

p <- prot.BP
r <- RNA.BP
b <- both.BP

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.BP <- rbind(r,p,b)
df.BP$V2 <- sub(" \\(.*","",df.BP$V2)
df.BP$V4 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.BP$V4)))
df.BP$V5 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.BP$V5)))
df.BP$V2 <- factor(x = df.BP$V2,levels = unique(df.BP$V2))

p.BP <- ggplot(data = df.BP,mapping = aes(x=Group,y=V5,fill=Group)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") +
  scale_fill_manual(values = c("forestgreen","skyblue3","goldenrod2")) +
  ggtitle(label = term) +
  facet_wrap(~V2,scales = "free_y") +
  xlab(label = "") +
  ylab(label = "% of genes") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p.BP

filename <- paste(data.dir,"/Plots/Changed_PercentageTerms_",term,".pdf",sep = "")
pdf(file = filename,width = 10,height = 8)
print(p.BP)
dev.off()


#####

##### M F  #######
term <- "Molecular Function (Slim)"

p <- prot.MF
r <- RNA.MF
b <- both.MF

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.MF <- rbind(r,p,b)
df.MF$V2 <- sub(" \\(.*","",df.MF$V2)
df.MF$V4 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.MF$V4)))
df.MF$V5 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.MF$V5)))
df.MF$V2 <- factor(x = df.MF$V2,levels = unique(df.MF$V2))

p.MF <- ggplot(data = df.MF,mapping = aes(x=Group,y=V5,fill=Group)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") +
  scale_fill_manual(values = c("forestgreen","skyblue3","goldenrod2")) +
  ggtitle(label = term) +
  facet_wrap(~V2,scales = "free_y") +
  xlab(label = "") +
  ylab(label = "% of genes") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p.MF

filename <- paste(data.dir,"/Plots/Changed_PercentageTerms_",term,".pdf",sep = "")
pdf(file = filename,width = 8,height = 8)
print(p.MF)
dev.off()

######

##### path #######
term <- "Panther Pathway"

p <- prot.Path
r <- RNA.Path
b <- both.Path

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.Path <- rbind(r,p,b)
df.Path$V2 <- sub(" \\(.*","",df.Path$V2)
df.Path$V4 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.Path$V4)))
df.Path$V5 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.Path$V5)))
df.Path <- df.Path[df.Path$V5>2,]
df.Path$V2 <- factor(x = df.Path$V2,levels = unique(df.Path$V2))

p.Path <- ggplot(data = df.Path,mapping = aes(x=Group,y=V5,fill=Group)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") +
  scale_fill_manual(values = c("forestgreen","skyblue3","goldenrod2")) +
  ggtitle(label = term) +
  facet_wrap(~V2,scales = "free_y") +
  xlab(label = "") +
  ylab(label = "% of genes") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p.Path

filename <- paste(data.dir,"/Plots/Changed_PercentageTerms_",term,".pdf",sep = "")
pdf(file = filename,width = 20,height = 12)
print(p.Path)
dev.off()

######

##### ProtC #######
term <- "Protein Class"

p <- prot.ProtC
r <- RNA.ProtC
b <- both.ProtC

p$Group <- "Prot only"
r$Group <- "RNA only"
b$Group <- "Both"

df.ProtC <- rbind(r,p,b)
df.ProtC$V2 <- sub(" \\(.*","",df.ProtC$V2)
df.ProtC$V4 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.ProtC$V4)))
df.ProtC$V5 <- as.numeric(gsub(pattern = "\\%",replacement="",x=as.character(df.ProtC$V5)))
df.ProtC$V2 <- factor(x = df.ProtC$V2,levels = unique(df.ProtC$V2))

p.ProtC <- ggplot(data = df.ProtC,mapping = aes(x=Group,y=V5,fill=Group)) +
  geom_bar(stat="identity",alpha=.9,position="dodge") +
  scale_fill_manual(values = c("forestgreen","skyblue3","goldenrod2")) +
  ggtitle(label = term) +
  facet_wrap(~V2,scales = "free_y") +
  xlab(label = "") +
  ylab(label = "% of genes") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray50",linetype = 3, size = 0.3),panel.grid.major.x = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))
p.ProtC

filename <- paste(data.dir,"/Plots/Changed_PercentageTerms_",term,".pdf",sep = "")
pdf(file = filename,width = 16,height = 10)
print(p.ProtC)
dev.off()

######

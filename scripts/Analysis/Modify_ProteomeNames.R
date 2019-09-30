## Change rownames protein ####
rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/Data_2018-10-19/")

proteome.input <- read.delim(file = "data/MatchedByGeneNames_HighOnly_normIMP.txt",skip=4,header=T)
head(proteome.input)

prot.det <- as.character(proteome.input$X.31)

# Change CIP2A to KIAA1524
prot.det.mod <- gsub(pattern = "CIP2A",replacement = "KIAA1524",x = prot.det.mod)
# Change VPS35L to C16orf62
prot.det.mod <- gsub(pattern = "VPS35L",replacement = "C16orf62",x = prot.det.mod)
# Change RTF2 to RTFDC1
prot.det.mod <- gsub(pattern = "RTF2",replacement = "RTFDC1",x = prot.det.mod)
# Change ATP5F1A to ATP5A1
prot.det.mod <- gsub(pattern = "ATP5F1A",replacement = "ATP5A1",x = prot.det.mod)
# Change CORO7-.AM16 to CORO7-PAM16
prot.det.mod <- gsub(pattern = "CORO7.PAM16",replacement = "CORO7-PAM16",x = prot.det.mod)
# Change RTRAF to C14orf166
prot.det.mod <- gsub(pattern = "RTRAF",replacement = "C14orf166",x = prot.det.mod)
# Change ODR4 to C1orf27
prot.det.mod <- gsub(pattern = "ODR4",replacement = "C1orf27",x = prot.det.mod)
# Change  ATP5F1B to ATP5B
prot.det.mod <- gsub(pattern = "ATP5F1B",replacement = "ATP5B",x = prot.det.mod)
# Change  PIP4P1 to TMEM55B
prot.det.mod <- gsub(pattern = "PIP4P1",replacement = "TMEM55B",x = prot.det.mod)
# Change  PIP4P2 to TMEM55A
prot.det.mod <- gsub(pattern = "PIP4P2",replacement = "TMEM55A",x = prot.det.mod)
# Change  MTREX to TMEM55A
prot.det.mod <- gsub(pattern = "MTREX",replacement = "SKIV2L2",x = prot.det.mod)
# Change  ECPAS to KIAA0368
prot.det.mod <- gsub(pattern = "ECPAS",replacement = "KIAA0368",x = prot.det.mod)
# Change  ATP5F1C to ATP5C1
prot.det.mod <- gsub(pattern = "ATP5F1C",replacement = "ATP5C1",x = prot.det.mod)
# Change  MRPL41 to MRPL27
prot.det.mod <- gsub(pattern = "MRPL41",replacement = "MRPL27",x = prot.det.mod)
# Change  LLGL1 to DLG4
prot.det.mod <- gsub(pattern = "LLGL1",replacement = "DLG4",x = prot.det.mod)
# Change  DGLUCY to C14orf159
prot.det.mod <- gsub(pattern = "DGLUCY",replacement = "C14orf159",x = prot.det.mod)

setdiff(prot.det.mod,prot.det)

proteome.input$X.31 <- prot.det.mod

write.table(x = proteome.input,file = "data/MatchedByGeneNames_HighOnly_normIMP_ModNames.txt",quote = F,append = F,sep = "\t",col.names = T,row.names = F)

library(colorspace)
library(MASS)
setwd("~/Documents/Projects/CCLE/R")
rm(list = ls())

CellLines <- read.table("./data/CCLE_sample_info_file_2012-10-18.txt", head=TRUE, sep = "\t", colClasses = c("character"))

RNAseq.file.name <- "./data/CCLE_RNAseq_081117.rpkm.gct"
RNAseq.dat <- read.table(RNAseq.file.name, head=TRUE, sep = "\t", skip = 2, na.strings = "null")
GenesIDs <- as.factor(sapply(RNAseq.dat$Name, function(x) unlist(strsplit(as.character(x), "[.]"))[1]))
RNAseq.dat$Name <- GenesIDs

# Filter cell lines to just include skin. 57 cell lines
SCL <- CellLines[CellLines$Histology %in% c("malignant_melanoma"),]
summary(as.factor(SCL$Hist.Subtype1))
summary(as.factor(SCL$Site.Primary))

MM.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% SCL$CCLE.name])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]
rm(RNAseq.dat)

##################################################################
#
# EMT signature correlates
#
##################################################################

#Skin "Differentiated" signature - 15 of 26 are shared with Breast Cancer Epithelial signature, which are:
# mitf: ALDH3B2 ERBB3 FXYD3 HPGD TMC6 
# not targeted by mitf: ARAP2 CKMT1A HOXC13 OVOL2 TUBBP5 VAV3
# down in WISP1 KO: ARAP2 VAV3
# up in WISP1 KO: CEACAM1
tmp <- matrix(data = c("BIK", -1.214385992, "GPR56", 4.184052198, "LLGL2", -0.599194087, 
                       "CTLA4", -2.047112291, "CGN", -0.553973887, "SP5", -4.908638476, 
                       "B3GAT1", -3.064132142, "DLL3", 1.310888085, "MITF", 3.893020066, 
                       "HPGD", -1.890261447, "MYH14", 0.89827663, "CCL3", -3.327471986, 
                       "ERBB3", 3.919007641, "MTUS1", 1.013823696, "TUBBP5", -3.399741752, 
                       "ARAP2", -2.054768322, "LEF1", 2.810160842, "FXYD3", 1.942820673, 
                       "ALDH3B2", -3.928680325, "CKMT1A", -1.93437492, "ESRP1", 0.371534458, 
                       "EDNRB", 3.299837111, "TMC6", 1.355537014, "FOXD3", 0.574941242, 
                       "CEACAM1", 1.030761688, "EN2", 0.443996142), ncol = 2, byrow = TRUE)

Esig <- data.frame(Genes = tmp[,1], NormFactor = as.numeric(tmp[,2]))

MM.RNAseq.E <- MM.RNAseq.dat[match(Esig$Genes, MM.RNAseq.dat$Description), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
log_EMM <- log2(apply(MM.RNAseq.E[,c(3:ncol(MM.RNAseq.E))], 2, mean) + 0.0001)
EMM <- apply(MM.RNAseq.E[,c(3:ncol(MM.RNAseq.E))]/2^Esig$NormFactor, 2, mean)/nrow(Esig)

#Check for correlation between epidermal signature genes and others
#cor(cbind(EBrCa, t(log2(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))] + 0.0001))))

EMM_Cor <- apply(log2(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))] + 0.0001), 1, function(x) cor(x,EMM))
Epi_Cor <- data.frame(GeneName = MM.RNAseq.dat$Description, Corr = EMM_Cor)

Epi_Cor <- Epi_Cor[abs(Epi_Cor$Corr) > 0.7 & !is.na(Epi_Cor$Corr),]

plot(density(EMM_Cor, adj = 0.2, na.rm = TRUE))

write.csv(Epi_Cor, file = "Skin-Epi-GeneCorrelation.csv")

##############################################################################
#Skin "De-differentiated" signature - 47 of these 84 are in BrCa mesenchymal signature, which are:
# POSTN COL3A1 GREM1 SERPINE1 FN1 CDH11 LOX COL6A2 FBN1 MMP2 SULF1   
# THY1COL6A3 VCAN COL6A1 ACTA2 PDGFRB COL1A1 TPM2 ADAM12 MYL9 FST     
# COL5A1 RCN3 FHL1 BGN MFAP5 SPOCK1 COL5A2 GJA1 IGFBP3 C1S THBS2  
# MAP1B ITGA5 CFH NID2 OLFML2B FGF1 ANKRD1 PCOLCE VEGFC PDGFC TFPI 
# TCF4 GLT8D2 WNT5A
tmp <- matrix(data = c("GREM1", -1.108527873, "COMP", -6.246580334, "COL1A1", 1.676171493, 
                       "CXCL12", -4.798731108, "MFAP5", -3.694334236, "POSTN", -1.026899727, 
                       "COL3A1", 0.037824991, "CDH11", -3.467002472, "THY1", -0.925747975, 
                       "WNT2", -8.204493412, "BGN", -0.426363339, "TWIST2", -4.001769723, 
                       "SERPINE1", 2.72622737, "COL5A1", 0.394575514, "DCN", -1.884229509, 
                       "SULF1", -0.576282599, "COL6A3", -0.592377762, "FOXC2", -3.773146543, 
                       "ASPN", -6.979012919, "PDGFRB", -1.590726383, "ITGBL1", -1.371222812, 
                       "MYL9", 2.069425824, "MXRA5", -4.885822053, "TNXB", -3.701116123, 
                       "LOX", 0.972368963, "GJA1", -0.941854563, "FBN1", 1.412787735, 
                       "SFRP4", -5.346991593, "COL6A2", 3.660198484, "PDGFRA", -2.318195624, 
                       "WISP1", -3.895709006, "TGFBI", 5.235744991, "NID2", -1.263634929, 
                       "ANKRD1", -0.554710348, "PAPPA", -2.050976984, "SERPINB2", 0.278276625, 
                       "CYP1B1", -1.275153916, "VEGFC", 0.106498104, "WNT5B", -1.107594273, 
                       "S100A4", 4.18650662, "VCAN", 1.501111971, "LOXL2", 4.250608876, 
                       "LRRC15", -1.637595346, "MMP2", 3.38648475, "PTGS1", -2.49454194, 
                       "INHBA", 0.400507034, "MALL", -2.577325698, "TPM2", 3.457653178, 
                       "EDNRA", -4.542401401, "ACTA2", 1.605676164, "CFB", -5.415954129, 
                       "TCF4", -1.615332903, "PRRX1", 0.66793088, "SPOCK1", 1.887874543, 
                       "LGR5", -5.435470896, "PLAU", -0.074808145, "FST", 2.365913639, 
                       "IL1R1", -0.498358101, "COL6A1", 5.066474543, "KRT14", -0.522035394, 
                       "ADAM12", 0.459366641, "COL5A2", 2.979012363, "NOTCH3", -0.378284442, 
                       "TFPI", 0.263926106, "WNT5A", 1.078965692, "RCN3", 3.112863026, 
                       "IFITM2", 3.009403733, "NT5E", 4.287225115, "FHL1", 1.477540026, 
                       "DES", -6.050789234, "RHOD", -2.420698019, "PDGFC", -0.079613497, 
                       "ZEB1", 0.465061124, "FGF1", -0.760882096, "GLT8D2", -1.374378, 
                       "IGFBP3", 5.044839966, "EGFR", -0.175279354, "FGF2", 0.741189508, 
                       "C1S", 1.482263317, "FAP", 1.929267289, "PCOLCE", 3.145042554, 
                       "PTRF", 5.236560361, "GADD45G", -5.405788006, "CLU", 0.69253925, 
                       "CFH", -0.810478538, "KRT7", 0.132794067, "KRT16", -3.456949727, 
                       "FN1", 8.745250548, "THBS2", 3.160974134, "EPS8L2", -0.744166138), ncol = 2, byrow = TRUE)

Msig <- data.frame(Genes = tmp[,1], NormFactor = as.numeric(tmp[,2]))

MM.RNAseq.M <- MM.RNAseq.dat[match(Msig$Genes, MM.RNAseq.dat$Description), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
log_MMM <- log2(apply(MM.RNAseq.M[,c(3:ncol(MM.RNAseq.M))], 2, mean) + 0.0001)
MMM <- apply(MM.RNAseq.M[,c(3:ncol(MM.RNAseq.M))]/2^Msig$NormFactor, 2, mean)/nrow(Msig)

#Check for correlation between mesenchymal signature genes and others
#cor(cbind(MBrCa, t(log2(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))] + 0.0001))))

MMM_Cor <- apply(log2(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))] + 0.0001), 1, function(x) cor(x,MMM))
Mesen_Cor <- data.frame(GeneName = MM.RNAseq.dat$Description, Corr = MMM_Cor)

Mesen_Cor <- Mesen_Cor[abs(Mesen_Cor$Corr) > 0.7 & !is.na(Mesen_Cor$Corr),]

plot(density(MMM_Cor, adj = 0.2, na.rm = TRUE))

write.csv(Mesen_Cor, file = "Skin-Mesen-GeneCorrelation.csv")

MM_State <- data.frame(CellLineName = colnames(MM.RNAseq.M)[c(3:ncol(MM.RNAseq.M))], Epithelial = EMM, Mesenchymal = MMM)
write.csv(MM_State, file = "Skin-EMT-State.csv")

MM_EMT <- read.csv(file = "Skin-EMT-State.csv")
##################################################################
#
#
#
##################################################################
#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
ColorDot <- c("black")
SlimCL <- lapply(as.character(MM_EMT$CellLineName), function(x) substr(x, 1, nchar(x)-5))

pdf("PCA-analysis-Skin-CellLinesProjection2.pdf", width = 7, height = 7)
plot(log2(MM_EMT$Mesenchymal + 0.0001), log2(MM_EMT$Epithelial + 0.0001), type = "p", pch = 19, col = ColorDot, xlab = "De-differentiatedl Signal", ylab = "Differentiated Signal", ylim = c(-10,-8), xlim = c(-2,3))

#legend(x="topright", legend = c("Basal", "Claudin Low", "HER2", "Luminal A", "Luminal B"), col = c("red", "yellow", "pink", "blue", "black"), pch = 19, cex = 1)
text(log2(MM_EMT$Mesenchymal + 0.0001), log2(MM_EMT$Epithelial + 0.0001)+0.5, labels = SlimCL, cex = 0.5)
dev.off() #



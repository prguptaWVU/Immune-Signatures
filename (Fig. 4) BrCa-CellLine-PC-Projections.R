library(colorspace)
library(MASS)
setwd("~/Documents/Projects/CCLE/R")
rm(list = ls())

CellLines <- read.table("./data/CCLE_sample_info_file_2012-10-18.txt", head=TRUE, sep = "\t", colClasses = c("character"))

RNAseq.file.name <- "./data/CCLE_RNAseq_081117.rpkm.gct"
RNAseq.dat <- read.table(RNAseq.file.name, head=TRUE, sep = "\t", skip = 2, na.strings = "null")

BrCaTypes.file.name <- "./data/BrCaCellLines-IntrinsicSubtypes.csv"
tmp <- read.csv(file = BrCaTypes.file.name)
BrCa_Types <- data.frame(CellName = tmp$X, Type = tmp$BrCaType)

#CNV.file.name <- "./data/CCLE_copynumber_byGene_2013-12-03.txt"
#CNV.dat <- read.table(CNV.file.name, head=TRUE, sep = "\t", na.strings = "null")

#InDel.file.name <- "./data/1650_HC_plus_RD_indels.maf.annotated"
#InDel.dat <- read.delim(InDel.file.name, head=TRUE, sep = "\t", na.strings = "null")

#Muts.file.name <- "./data/1650_HC_plus_RD_muts.maf.annotated"
#Muts.dat <- read.delim(Muts.file.name, head=TRUE, sep = "\t", na.strings = "null")

#Header.dat <- scan(CNV.file.name, sep = "\t", what = "character", nlines = 2)

# Filter cell lines to just include breast. 57 cell lines
SCL <- CellLines[CellLines$Site.Primary %in% c("breast"),]
summary(as.factor(SCL$Hist.Subtype1))
summary(as.factor(SCL$Histology))

#jnk <- colnames(RNAseq.dat)
#SCL$CCLE.name[!(SCL$CCLE.name %in% colnames(RNAseq.dat))]

#MM.CNV.dat <- cbind(CNV.dat[,c(1:5)], CNV.dat[,colnames(CNV.dat) %in% MMCL$CCLE.name])
BC.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% SCL$CCLE.name])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))])]

# # List from Tan et al. EMBO J 2014
# # List from Kaiser et al. Biotech Prog 2016
# List from references cited in Tan et al. EMBO J 2014 - MSigDB and others
# "C10orf116", "C12orf24", "C13orf15", "C16orf57", "C16orf61", "C19orf21", "C7orf68", "C9orf140", "C9orf30", "CXCR7", "GAGE8", 
# "GALNTL2", "HNT", "KLRK1", "LOC100507424", "LOC727820", "LOC729082", "LOC731429", "MALAT1", "PLOCE", "RPS27AP11", "RRS1", 
# "S1A4", "S1A8", "SERPINA3", "SIP1", "SRSF8", "TCF8", "TCFE2A", "TMEM30B", "ZFHX1B", 
EMTsig <- c("ABCC3", "ACSL1", "ACTA1", "ACTA2", "ADA", "ADAM12", "ADAP1", "ADH1C", "ADSSL1", "AEBP1", "AGAP1", "AGR2", "AIM1", "AKAP12", 
             "AKAP2", "AKT1", "AKT2", "AKT3", "ALDH3B2", "ALG6", "ALOX5", "ANK2", "ANK3", "ANKRD1", "ANLN", "ANXA4", "ANXA8L2", "ANXA9", 
             "AP1M2", "AP1S2", "APLP2", "AQP3", "ARAP2", "AREG", "ARHGAP8", "ARHGDIB", "ARHGEF18", "ARHGEF5", "ARPC1B", "ARPC5L", "ASF1B", 
             "ASNS", "ASPM", "ASPN", "ATF3", "ATP1A3", "ATP1B1", "ATP2C2", "AXIN2", "AXL", "B2M", "B3GAT1", "B3GNT3", "BAG2", "BARD1", 
             "BCL2L11", "BCL3", "BCL6", "BGN", "BIK", "BIRC3", "BLM", "BLNK", "BMP2", "BMP4", "BMP7", "BNIP3L", "BRCA1", "BRIX1", "BSPRY", 
             "BTG1", "BUB1", "C1orf106", "C1orf116", "C1orf63", 
             "C1orf74", "C1QTNF3", "C1S", "C4orf19", "C7orf10", "CA2", "CALD1", "CALM3", "CAMK2N1", 
             "CAPN1", "CAPZB", "CBLC", "CCL2", "CCL21", "CCL3", "CCNA2", "CCNB1", "CCNB2", "CCNG1", "CCR2", "CCT4", "CD164", "CD46", 
             "CD63", "CD68", "CD9", "CDC45", "CDC6", "CDC7", "CDCA5", "CDCA7", "CDH1", "CDH10", "CDH11", "CDH2", "CDH3", "CDK1", "CDK2", 
             "CDK2AP1", "CDKN1C", "CDS1", "CEACAM1", "CEACAM6", "CENPA", "CENPE", "CENPW", "CEP170", "CEP55", "CFB", "CFH", "CGN", "CHN1", 
             "CIT", "CITED2", "CKMT1A", "CKS1B", "CLDN1", "CLDN4", "CLDN7", "CLIC4", "CLU", "CNKSR1", "COL10A1", "COL11A1", "COL1A1", 
             "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COMP", "COMT", "COPZ2", "CORT", "CPT1A", "CREG1", "CRISPLD2", 
             "CRLF3", "CRY1", "CTLA4", "CTNNB1", "CTSB", "CTSD", "CTSH", "CTSK", "CTSZ", "CX3CR1", "CXCL1", "CXCL12", "CXCL3", "CXCL5", 
             "CXCR4", "CYB561", "CYBRD1", "CYP1B1", "CYP2A13", "CYP4B1", "CYP4F3", "DAB2", "DAPK1", "DCK", "DCN", "DDR1", "DDR2", 
             "DENND2D", "DENND5A", "DES", "DHFR", "DHH", "DIAPH3", "DLG5", "DLL1", "DLL3", "DLL4", "DSC2", "DSC3", "DSCC1", "DSG2", "DSP", 
             "DST", "DTL", "DTX4", "E2F1", "ECT2", "EDN2", "EDNRA", "EDNRB", "EEF2", "EFNA1", "EGF", "EGFR", "EGR1", "EHF", "EIF1AD", 
             "EIF4B", "EIF5A2", "ELF3", "ELMO3", "EMP1", "EMP3", "EN2", "ENC1", "ENG", "ENO1", "ENOPH1", "EOMES", "EPAS1", "EPCAM", 
             "EPHA1", "EPN3", "EPS8L1", "EPS8L2", "EPYC", "ERBB2", "ERBB3", "ERMP1", "ERRFI1", "ESPL1", "ESRP1", "ESRP2", "ETS1", "ETV1", 
             "EVPL", "EXO1", "EXPH5", "F11R", "F3", "FA2H", "FABP5", "FAP", "FBN1", "FBP1", "FBXO32", "FBXO41", "FBXO5", "FERMT2", "FGF1", 
             "FGF18", "FGF2", "FGF4", "FGFBP1", "FGFR1", "FGFR2", "FGG", "FHL1", "FHL2", "FJX1", "FKBP1A", "FMO1", "FN1", "FOS", "FOSB", 
             "FOSL1", "FOSL2", "FOXA1", "FOXC2", "FOXD3", "FST", "FSTL1", "FXYD3", "FZD1", "FZD7", "GADD45G", "GALK1", "GALNT3", 
             "GAS1", "GBP2", "GDA", "GEM", "GFPT2", "GGCT", "GJA1", "GLI1", "GLI2", "GLI3", "GLT8D2", "GLUD1", "GLYR1", "GMNN", 
             "GPR56", "GPX2", "GRB7", "GREM1", "GRHL2", "GSC", "H2AFX", "H2AFY", "HELLS", "HES1", "HGF", "HIF1A", "HIST1H2AC", "HIST1H4B", 
             "HIST1H4C", "HIST1H4J", "HIST2H2BE", "HLA-E", "HMGA2", "HMMR", "HNF1A", "HNMT", "HNRNPAB", "HOXC13", "HPGD", "HSPE1", 
             "HTRA1", "ICA1", "ID2", "IDH2", "IFI16", "IFI44", "IFIT3", "IFITM2", "IFITM3", "IFNGR1", "IGFBP3", "IHH", "IL11", "IL1R1", 
             "IL1RN", "IL20RA", "ILK", "INHBA", "INSIG1", "IRF6", "IRF7", "ITGA5", "ITGAV", "ITGB1", "ITGB4", "ITGB6", "ITGBL1", "JAG1", 
             "JAG2", "JUN", "JUP", "KCNK1", "KDELC1", "KIF11", "KIF23", "KLF5", "KLF6", "KLK6", "KLRC3", "KNTC1", "KRT13", 
             "KRT14", "KRT15", "KRT16", "KRT18", "KRT19", "KRT7", "KRT8", "KYNU", "LAD1", "LAMB2", "LCN2", "LEF1", "LEPRE1", "LETM1", 
             "LGALS1", "LGR5", "LHFP", "LLGL2", "LOX", "LOXL2", "LPAR6", "LRRC15", 
             "LRRC8C", "LSR", "LUM", "LY75", "MAD2L1", "MAFB", "MALL", "MANSC1", "MAOA", "MAP1B", "MAP2K5", "MAP3K1", "MAP3K5", 
             "MAP7", "MAPK1", "MAPK13", "MAPK14", "MAPK7", "MAPK8", "MBOAT7", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MET", 
             "MFAP1", "MFAP5", "MFHAS1", "MGST2", "MITF", "MKI67", "MKNK2", "MLF1IP", "MME", "MMP11", "MMP12", "MMP13", "MMP14", "MMP2", 
             "MMP25", "MMP3", "MMP7", "MMP9", "MPZL2", "MSH2", "MSN", "MST1R", "MSX1", "MSX2", "MT2A", "MTHFD2", "MTUS1", "MVP", "MXI1", 
             "MXRA5", "MXRA7", "MYB", "MYBL1", "MYC", "MYCBP", "MYH14", "MYL9", "MYO1D", "MYO5C", "MYO6", "NCK1", "NDC80", "NEAT1", 
             "NFE2L1", "NFE2L2", "NFIL3", "NFKB1", "NFKBIA", "NID2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NOX4", "NRCAM", "NT5E", 
             "NUAK1", "NUP155", "OAS1", "OCLN", "ODC1", "OLFML2B", "OR7E14P", "ORAOV1", "ORC6", "OVOL2", "PAK6", "PAPPA", "PARD3", 
             "PARD6A", "PC", "PCDH9", "PCNA", "PCOLCE", "PDGFA", "PDGFB", "PDGFC", "PDGFD", "PDGFRA", "PDGFRB", "PDZK1IP1", "PERP", 
             "PFN1", "PGK1", "PHF19", "PHLDA1", "PIK3C2B", "PIK3CA", "PITX2", "PKMYT1", "PKP3", "PLA2G7", "PLAU", "PLAUR", "PLCG2", 
             "PLEK2", "PLEKHO2", "PLK1", "PLK4", "PLS1", "PLXNB1", "PLXNB2", "PMP22", "PNRC1", "POF1B", "POLD3", "POLR3E", 
             "POLR3K", "POPDC3", "POSTN", "PPARD", "PPFIBP2", "PPIC", "PPL", "PRC1", "PRIM1", "PRKCA", "PRKCZ", "PROCR", "PRR4", "PRRX1", 
             "PRSS8", "PSAP", "PSME1", "PTEN", "PTGS1", "PTGS2", "PTK6", "PTPN22", "PTRF", "PYCARD", "RAB11FIP1", "RAB1A", "RAB20", 
             "RAB25", "RAB31", "RABGAP1L", "RAC3", "RAD51", "RAD51B", "RAN", "RARRES1", "RBL1", "RBM14", "RBM47", "RCN3", "RECK", 
             "RECQL4", "RELA", "RELT", "RFC2", "RFC3", "RFC4", "RFC5", "RGS2", "RHOA", "RHOB", "RHOD", "RHOU", "RNASET2", "RNF145", 
             "ROS1", "RPL37", "RPLP1", "RPS19BP1", "RPS27A", "RRAGA", "RRAS", "RRM2", "RTEL1", "RUNX1", "S100A14", 
             "S100A4", "S100A8", "S100P", "SAA2", "SACS", "SAT1", "SCD", "SCEL", "SCNN1A", "SCPEP1", "SDC1", "SDC2", 
             "SDCBP", "SDF2L1", "SEMA3B", "SEMA3C", "SEPP1", "SERINC5", "SERPINB2", "SERPINE1", "SERPINE2", "SERPINF1", 
             "SFN", "SFRP4", "SH2B3", "SH2D3A", "SH3KBP1", "SH3YL1", "SHANK2", "SHH", "SKA3", "SLC20A1", "SLC37A1", "SLC3A2", 
             "SLC9A3R1", "SLPI", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMARCA1", "SMTN", "SMURF1", "SNAI1", "SNAI2", "SNAI3",  
             "SNCAIP", "SND1", "SOAT1", "SOBP", "SOD2", "SORL1", "SOX1", "SOX9", "SP5", "SPA17", "SPARC", "SPINT1", "SPINT2", "SPOCK1", 
             "SPRR1B", "SPRR3", "SRM", "SRPX", "SRXN1", "SSH3", "ST14", "STAMBPL1", "STAP2", "STAT1", "STAT3", "STEAP1", "STIL", 
             "STK17A", "STMN1", "STOM", "STRA6", "SULF1", "SYNGR2", "TACC3", "TACSTD2", "TCF3", "TCF4", "TFDP1", 
             "TFPI", "TGFA", "TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBI", "TGFBR1", "TGFBR2", "TGFBR3", "THBS2", "THBS3", "THY1", 
             "TIMM10", "TIMP3", "TJP2", "TKT", "TMC6", "TMEFF1", "TMEM158", "TMEM59", "TMPO", "TMPRSS2", "TMPRSS4", "TNC", 
             "TNFAIP3", "TNFAIP6", "TNFRSF10A", "TNFRSF12A", "TNFRSF19", "TNFSF13", "TNK1", "TNXB", "TOB1", "TOP2A", "TOX3", "TP53INP1", 
             "TP63", "TPM2", "TRIM28", "TRIM29", "TRIM31", "TRPC1", "TSC22D3", "TSPAN1", "TSPAN13", "TSPAN15", "TTC39A", "TUBA1A", 
             "TUBA3C", "TUBA3D", "TUBA4A", "TUBA8", "TUBB", "TUBB3", "TUBB4B", "TUBB6", "TUBBP5", "TUBG1", "TWIST1", "TWIST2", "TYMS", 
             "UBA52", "UBE2H", "UBE2S", "UBE2T", "UCK2", "UGGT2", "UGT1A1", "UHRF1", "UPP1", "USP15", "VAMP8", "VAV3", "VCAN", "VDAC1P2", 
             "VEGFA", "VEGFC", "VGLL1", "VIM", "VLDLR", "VNN1", "WISP1", "WNT1", "WNT10A", "WNT10B", "WNT11", "WNT16", "WNT2", "WNT2B", 
             "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "XBP1", "YEATS4", 
             "YRDC", "YWHAH", "ZBTB16", "ZCCHC6", "ZEB1", "ZNF165", "ZNF367", "ZNF608", "ZWINT")


BC.RNAseq.EMT <- BC.RNAseq.dat[BC.RNAseq.dat$Description %in% EMTsig, ]

#Convert from RPKM to log2(RPKM + 0.0001)
SeqRes <- log2(BC.RNAseq.EMT[,c(3:ncol(BC.RNAseq.EMT))] + 0.0001)

# Principal component analysis
opar <- par(mfrow = c(2,2))
# Number of genes x Number of Cell lines
# t(SeqRes) gives PCA clustering of cell lines
fs1PCA <- prcomp(t(SeqRes), retx = TRUE, scale. = FALSE)

#Scree plot
EigenVals <- fs1PCA$sdev^2
TotVar <- sum(EigenVals)
PrPC <- c(1:10)
#PrPC[1] <- EigenVals[1]/
for ( i in 1:10 )
{
  PrPC[i] <- EigenVals[i]/TotVar
}

load.Matrix <- data.frame(fs1PCA$rotation)
# n_genes x n_PC          =   n_genes x n_pts      *  n_pts x n_PC
# (square n_genes = n_PC) = (here n_genes < n_pts) * (n_pts > n_PC)
PCA.data <- as.matrix(t(SeqRes)) %*% fs1PCA$rotation

BCType <- BrCa_Types$Type[match(colnames(BC.RNAseq.EMT)[c(3:ncol(BC.RNAseq.EMT))], BrCa_Types$CellName)]
DotColors <- c("red", "yellow", "pink", "blue", "black")
ColorDot <- DotColors[as.factor(BCType)]
SlimCL <- lapply(as.character(colnames(BC.RNAseq.EMT)[c(3:ncol(BC.RNAseq.EMT))]), function(x) substr(x, 1, nchar(x)-7))

pdf("PCA-analysis-BrCa-CellLines.pdf", width = 7, height = 7)
opar <- par(mfrow = c(1,1))
mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,40), ylab = "Variance (%)", xlab = "Principal Component")

plot(PCA.data[,c(1,2)], type = "n", xlim = c(-80, 100), ylim = c(-55, 55))
points(PCA.data[,c(1,2)], pch = 19, col = ColorDot, cex = 1)
legend(x="bottomright", legend = c("Basal", "Claudin Low", "HER2", "Luminal A", "Luminal B"), col = c("red", "yellow", "pink", "blue", "black"), pch = 19, cex = 1)
text(PCA.data[,1], PCA.data[,2]+0.1, labels = SlimCL, cex = 0.5)

dev.off() # End plot

# # Number of Cell lines x Number of genes
# # SeqRes gives PCA clustering of genes
# # t(SeqRes) gives PCA clustering of cell lines
# fs2PCA <- prcomp(SeqRes, retx = TRUE, scale. = FALSE)
# 
# # Eigenvals 819.3 172.5 32.4 20.1 12.1 11.0 10.18 8.73 8.05
# #Scree plot
# EigenVals <- fs2PCA$sdev^2
# TotVar <- sum(EigenVals)
# PrPC <- c(1:10)
# #PrPC[1] <- EigenVals[1]/
# for ( i in 1:10 )
# {
#   PrPC[i] <- EigenVals[i]/TotVar
# }
# 
# load.Matrix <- data.frame(fs2PCA$rotation)
# # n_genes x n_PC          =   n_genes x n_pts      *  n_pts x n_PC
# # (square n_genes = n_PC) = (here n_genes < n_pts) * (n_pts > n_PC)
# PCA2.data <- as.matrix(SeqRes) %*% fs2PCA$rotation
# 
# pdf("PCA-analysis-BrCa-Genes-Mar19.pdf", width = 7, height = 7)
# opar <- par(mfrow = c(2,2))
# mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,70), ylab = "Variance (%)", xlab = "Principal Component")
# lines(c(0,2,3,4,5,6,7,8,9,12), c(3.4, 3.2, 3.19, 3.04, 3.01, 2.89, 2.8, 2.74, 2.67, 2.63), col="red", lwd = 2, lty = 2)
# 
# mp <- barplot(EigenVals[c(1:10)], names.arg = c(1:10), ylim = c(0,900), ylab = "Variance", xlab = "Principal Component")
# lines(c(0,2,3,4,5,6,7,8,9,12), c(34.3, 34.1, 33.5, 32.2, 31.7, 31.4, 30.0, 29.9, 29.4, 28.1), col="red", lwd = 2, lty = 2)
# 
# # Keep gene for further analysis if it is expressed above the threshold in more than 
# # 5% of the samples
# GeneExpressed <- apply(SeqRes, 1, function(x) length(x[x > -0.304])/length(x) > 0.05 )
# rainColor <- colorRampPalette(c("red", "yellow", "blue"))(100)
# meanSeqRes <- rowMeans(SeqRes)
# PtCol <- rainColor[cut(meanSeqRes,breaks = 100)]
# #PtCol <- ifelse(meanSeqRes > -0.304, "black", "red")
# #PtCol <- ifelse(GeneExpressed, "black", "red")
# 
# opar <- par(mfrow = c(1,1))
# plot(-1*PCA2.data[,1], PCA2.data[,2],type = "n", xlim = c(-110, 80), ylim = c(-40, 40))
# #points(PCA2.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
# text(-1*PCA2.data[,1], PCA2.data[,2], labels = BC.RNAseq.EMT$Description, cex = 1, col = PtCol)
# points(-1*PCA2.data[BC.RNAseq.EMT$Description == "CDH1",1], PCA2.data[BC.RNAseq.EMT$Description == "CDH1",2], pch = 19, col = "blue", cex = 1)
# points(-1*PCA2.data[BC.RNAseq.EMT$Description == "VIM",1], PCA2.data[BC.RNAseq.EMT$Description == "VIM",2], pch = 19, col = "yellow", cex = 1)
# #lines(c(-100, 70), c(11.9, 11.9), col="red", lwd=2, lty=2)
# #lines(c(-100, 70), c(-11.9, -11.9), col="red", lwd=2, lty=2)
# #lines(c(30, 30), c(-100, 70), col="blue", lwd=2, lty=2)
# 
# KeepGene <- GeneExpressed
# # KeepGene <- meanSeqRes > -0.304
# # KeepGene <- PCA2.data[,1] < 30 
# plot(PCA2.data[KeepGene,c(2,3)], type = "n", xlim = c(-30, 30), ylim = c(-30, 30))
# #points(PCA2.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
# text(PCA2.data[KeepGene,c(2,3)], labels = BC.RNAseq.EMT$Description[KeepGene], cex = 0.5, col = "black" ) #PtCol[KeepGene]
# points(PCA2.data[BC.RNAseq.EMT$Description == "CDH1",2], PCA2.data[BC.RNAseq.EMT$Description == "CDH1",3], pch = 19, col = "blue", cex = 1)
# points(PCA2.data[BC.RNAseq.EMT$Description == "VIM",2], PCA2.data[BC.RNAseq.EMT$Description == "VIM",3], pch = 19, col = "yellow", cex = 1)
# lines(c(-110, 80), c(11.9, 11.9), col="red", lwd=2, lty=2)
# lines(c(-110, 80), c(-11.9, -11.9), col="red", lwd=2, lty=2)
# lines(c(-11.9, -11.9), c(-100, 70), col="blue", lwd=2, lty=2)
# lines(c(11.9, 11.9), c(-100, 70), col="blue", lwd=2, lty=2)
# 
# dev.off() # End plot 

##################################################################
#
# EMT signature correlates
#
##################################################################

#Epithelial signature
# tmp <- matrix(data = c("CYP4B1", -2.446785084, "ARHGAP8", -2.26795788, "ALDH3B2", 1.435970192, 
#                        "FA2H", -0.5285652, "BSPRY", 1.792107746, "F11R", -3.65098147, "GRHL2", 2.272778131, 
#                        "OVOL2", 0.060031635, "CKMT1A", 1.912022287, "F11R", -6.502068123, "BLNK", -1.318945039, 
#                        "FXYD3", 3.329906423, "FOXA1", 2.694653058, "RAB25", 3.799895055, "FBP1", 1.101376263, 
#                        "TMPRSS2", -0.278656984, "WNT3A", -4.758677792, "ESRP1", 3.205173557, "EPN3", 2.478648163, 
#                        "TUBBP5", -1.176766827, "PRSS8", 3.62237251, "EVPL", 2.310660775, "S100A14", 4.668678552, 
#                        "C4orf19", -3.344019889, "OR7E14P", -0.610740935, "BIK", 1.175798652, "AGR2", 3.328402466, 
#                        "GRB7", 2.928559204, "ATP2C2", 1.407978993, "EPHA1", 1.463343712, "LAD1", 2.832460693, 
#                        "ST14", 3.857103038, "MYO5C", 1.827949913, "ESRP2", 3.149686814, "PDGFB", 1.634909829, 
#                        "LLGL2", 3.503815163, "AP1M2", 4.059617591, "TTC39A", 1.093368035, "SPINT1", 4.28176267, 
#                        "ICA1", 1.255665985, "SCNN1A", 2.045415355, "ELF3", 2.930432283, "MYH14", 2.611742248, 
#                        "POF1B", -1.061477381, "IRF6", 2.512716819, "CLDN4", 4.439361084, "EPCAM", 5.035654841, 
#                        "SEPP1", 0.99912883, "EHF", 2.218922409, "CLDN7", 4.421022987, "HOXC13", 2.023664943, 
#                        "VAV3", 1.255238988, "DSC2", 1.588680098, "TSPAN15", 3.95722022, "SLC37A1", 1.632150536, 
#                        "CDH1", 3.978846105, "ERBB3", 3.496068379, "IL1RN", -2.209158624, "VAMP8", 5.014175672, 
#                        "CX3CR1", -5.547077396, "PTK6", 1.411144861, "EDN2", -0.110772193, "CGN", 2.421659511, 
#                        "TMC6", 2.429819726, "CEACAM6", 0.388508527, "SPINT2", 6.814914712, "MSX2", 0.707879411, 
#                        "EFNA1", 3.376318109, "MAP7", 2.769827993, "C1orf106", 0.724307937, "TSPAN1", 2.87527624, 
#                        "SHANK2", 0.378890092, "ANXA9", 2.246682225, "WNT4", -1.178900939, "B3GAT1", -5.761322743, 
#                        "BMP7", -0.636426934, "CDS1", 2.462592987, "PPL", 2.507467951, "MST1R", 1.579844441, 
#                        "SORL1", 1.279696831, "IL20RA", -2.822257097, "ANK3", 0.481008272, "KRT8", 8.088252672, 
#                        "GADD45G", 0.527432699, "DENND2D", 0.616714594, "EPS8L1", 2.770283683, "JUP", 5.202752406, 
#                        "CNKSR1", 1.661710701, "SH2D3A", 2.276145054, "ARAP2", -0.499625902, "CXCR4", -0.115041413, 
#                        "HPGD", -2.759089613, "LSR", 4.538511283, "EXPH5", -0.065494998, "CEACAM1", -0.99607925, 
#                        "GALNT3", 2.556610315, "RBM47", 2.886676807, "MYB", 0.970250216, "PKP3", 3.954563786, 
#                        "WNT7B", 1.672477704, "OCLN", 2.162243252, "GPX2", -0.062285999), ncol = 2, byrow = TRUE)

tmp <- matrix(data = c("CYP4B1", -2.930755502, "ARHGAP8", -2.737894731, "ALDH3B2", 0.925513024, "FA2H", -0.968204588,
                        "BSPRY", 1.369225794, "F11R", -4.070036907, "GRHL2", 1.892342621, "OVOL2", -0.30078643, 
                       "CKMT1A", 1.562032689, "F11R", -6.931774755, "BLNK", -1.673881837, "FXYD3", 2.967212661, 
                       "FOXA1", 2.288768997, "RAB25", 3.446292345, "FBP1", 0.724376998, "TMPRSS2", -0.624672513, 
                       "WNT3A", -5.04273598, "ESRP1", 2.890557471, "EPN3", 2.105650044, "TUBBP5", -1.502082985, 
                       "PRSS8", 3.315563824, "EVPL", 1.970005485, "S100A14", 4.349447482, "C4orf19", -3.712172017, 
                       "OR7E14P", -0.937264451, "BIK", 0.87253461, "AGR2", 2.965597498, "GRB7", 2.594647753, 
                       "ATP2C2", 1.08865336, "EPHA1", 1.159691669, "LAD1", 2.54637928, "ST14", 3.581367197, 
                       "MYO5C", 1.519849397, "ESRP2", 2.845976426, "PDGFB", 1.328508982, "LLGL2", 3.238111023, 
                       "AP1M2", 3.773188064, "TTC39A", 0.782003305, "SPINT1", 4.03245131, "ICA1", 0.960966321, 
                       "SCNN1A", 1.827126949, "ELF3", 2.701826979, "MYH14", 2.385552519, "POF1B", -1.282651427, 
                       "IRF6", 2.289520448, "CLDN4", 4.187521817, "EPCAM", 4.808169937, "SEPP1", 0.731504121, 
                       "EHF", 1.968944877, "CLDN7", 4.182235248, "HOXC13", 1.815842128, "VAV3", 1.007058023, 
                       "DSC2", 1.387721981, "TSPAN15", 3.72918076, "SLC37A1", 1.409167363, "CDH1", 3.792806331, 
                       "ERBB3", 3.224466404, "IL1RN", -2.388289369, "VAMP8", 4.78755525, "CX3CR1", -5.765761708, 
                       "PTK6", 1.169819059, "EDN2", -0.32606694, "CGN", 2.201078372, "TMC6", 2.232191228, 
                       "CEACAM6", 0.14387677, "SPINT2", 6.609323447, "MSX2", 0.437999837, "EFNA1", 3.177855185, 
                       "MAP7", 2.57947976, "C1orf106", 0.585614563, "TSPAN1", 2.653517551, "SHANK2", 0.144306741, 
                       "ANXA9", 2.008795445, "WNT4", -1.359998531, "B3GAT1", -5.898091275, "BMP7", -0.811401435, 
                       "CDS1", 2.279520695, "PPL", 2.330982964, "MST1R", 1.378381813, "SORL1", 1.096615392, 
                       "IL20RA", -2.950370686, "ANK3", 0.274722527, "KRT8", 7.912717037, "GADD45G", 0.354412218, 
                       "DENND2D", 0.400939541, "EPS8L1", 2.589570102, "JUP", 5.052623026, "CNKSR1", 1.506935172, 
                       "SH2D3A", 2.103163142, "ARAP2", -0.653338607, "CXCR4", -0.315141164, "HPGD", -2.807736224, 
                       "LSR", 4.419957081, "EXPH5", -0.235253428, "CEACAM1", -1.086964078, "GALNT3", 2.401095422, 
                       "RBM47", 2.718693097, "MYB", 0.762850161, "PKP3", 3.806420702, "WNT7B", 1.552310838, 
                       "OCLN", 2.000057913, "GPX2", -0.157037495), ncol = 2, byrow = TRUE)

Esig <- data.frame(Genes = tmp[,1], NormFactor = as.numeric(tmp[,2]))

BC.RNAseq.E <- BC.RNAseq.dat[match(Esig$Genes, BC.RNAseq.dat$Description), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
log_EBrCa <- log2(apply(BC.RNAseq.E[,c(3:ncol(BC.RNAseq.E))], 2, mean) + 0.0001)
EBrCa <- apply(BC.RNAseq.E[,c(3:ncol(BC.RNAseq.E))]/2^Esig$NormFactor, 2, mean)

#Mesenchymal signature
# tmp <- matrix(data = c("SACS", 0.290062684, "SH3KBP1", 1.466203798, "BAG2", 0.49636556, "CLIC4", 5.159359081, 
#               "C7orf10", -2.51888532, "CXCL3", -3.423618245, "S100A4", 3.778864441, "SDC2", 1.237691836, 
#               "PHLDA1", 2.195399561, "WNT5A", 0.111286371, "CTSB", 4.9741377, "CYBRD1", 1.525324739, 
#               "GLT8D2", -1.47153745, "TCF4", -2.335071698, "TRPC1", -0.827329743, "PITX2", -3.782760584, 
#               "PLAUR", 1.563004224, "MXRA7", 1.509509711, "LEPRE1", 2.79209761, "TFPI", 0.228503595, 
#               "ASPN", -6.336878371, "PMP22", 2.066637029, "GFPT2", -0.323637449, "PDGFC", 0.450379972, 
#               "ITGB1", 6.585941966, "IFITM3", 4.96794764, "GEM", -0.126991774, "COPZ2", -0.429158273, 
#               "VEGFC", 0.509088251, "CD68", 2.005694138, "MME", -1.977989143, "ANK2", -3.891723139, 
#               "SMARCA1", 1.403981074, "TGFB1", 1.995928221, "AKAP2", -4.329024029, "SFRP4", -5.868109735, 
#               "PCOLCE", 1.747514167, "ANKRD1", -1.004379089, "LRRC15", -4.410930186, "WNT2", -8.414143146, 
#               "FGF1", -3.573058761, "OLFML2B", -2.023325465, "HMGA2", -2.507777882, "PROCR", 1.002464008, 
#               "TNC", 1.4768736, "AEBP1", -1.794438431, "TWIST1", -1.265599572, "EDNRA", -5.015670378, 
#               "COMP", -5.895696801, "DAB2", 0.311142128, "TMEM158", -0.813402427, "MMP3", -4.774672962, 
#               "AKAP12", -0.516119842, "CCL2", -1.388376675, "GAS1", -2.551990724, "NID2", -1.282312728, 
#               "WISP1", -6.762209988, "CFH", -2.157779226, "LHFP", 0.23020315, "TUBB6", 4.06764218, 
#               "DDR2", -0.75369598, "TIMP3", 2.43646936, "GLI2", -4.47772844, "LGALS1", 7.82056313, 
#               "PAPPA", -3.189433741, "WNT5B", -1.566230649, "AKT3", -0.7568643, "ITGA5", 2.145434418, 
#               "TMEFF1", -7.602217201, "MAP1B", 0.486147967, "THBS2", -1.820681347, "TGFB1I1", -0.249204358, 
#               "SERPINB2", -3.771357136, "SRPX", -0.647854286, "C1S", 0.286522657, "FOSL1", 1.575046264, 
#               "IGFBP3", 3.712917933, "RECK", -1.442349154, "GJA1", -0.346756644, "COL5A2", 1.216170119, 
#               "SPOCK1", -0.183911903, "MFAP5", -3.363210987, "BGN", 0.316022331, "FHL1", -0.979230602, 
#               "HTRA1", 1.268220144, "RCN3", 0.754204477, "COL5A1", 0.29846138, "ZEB1", -2.478131852, 
#               "FST", -0.991784444, "MYL9", 3.176709968, "EMP3", 1.618830309, "MMP14", 1.941605043, 
#               "LUM", -2.379093261, "ADAM12", -2.253021028, "CALD1", 2.704644634, "TPM2", 3.062659226, 
#               "COL1A1", 3.467409961, "PDGFRB", -3.71542655, "ACTA2", 0.335708898, "FOXC2", -5.681890132, 
#               "AXL", 0.449483087, "COL6A1", 2.630076832, "VCAN", -1.614362995, "CDH2", -0.225036023, 
#               "DCN", -4.730347576, "COL6A3", -1.540084326, "THY1", -2.697321101, "SERPINE2", 1.735635925, 
#               "SULF1", -2.00066345, "PRRX1", -4.144931635, "PDGFRA", -5.205499978, "FSTL1", 2.419411186, 
#               "MMP2", -0.265964685, "FBN1", -0.322073696, "LOXL2", 1.876172565, "COL6A2", 1.499796853, 
#               "LOX", 0.197294329, "CDH11", -3.293231657, "FAP", -3.665707164, "FN1", 4.354138118, 
#               "SERPINE1", 2.091363265, "SPARC", 2.448393824, "TWIST2", -5.342141994, "GREM1", -2.307420799, 
#               "TNFAIP6", -6.738435166, "VIM", 3.384074296, "COL3A1", -0.391285122, "POSTN", -2.389306523), ncol = 2, byrow = TRUE)

tmp <- matrix(data = c("SACS", 0.472970206, "SH3KBP1", 1.666039892, "BAG2", 0.678533429, "CLIC4", 5.357695319, 
                       "C7orf10", -2.278524021, "CXCL3", -3.158784335, "S100A4", 3.979290912, "SDC2", 1.438793936, 
                       "PHLDA1", 2.426884449, "WNT5A", 0.328355237, "CTSB", 5.196116557, "CYBRD1", 1.706871501, 
                       "GLT8D2", -1.220461532, "TCF4", -2.147816288, "TRPC1", -0.603392817, "PITX2", -3.554157614, 
                       "PLAUR", 1.791729101, "MXRA7", 1.713568899, "LEPRE1", 3.009541309, "TFPI", 0.411657212, 
                       "ASPN", -6.17208793, "PMP22", 2.257755866, "GFPT2", -0.080852634, "PDGFC", 0.685190065, 
                       "ITGB1", 6.814265856, "IFITM3", 5.240077831, "GEM", 0.093927677, "COPZ2", -0.211689908, 
                       "VEGFC", 0.725773099, "CD68", 2.258521576, "MME", -1.681553355, "ANK2", -3.649720552, 
                       "SMARCA1", 1.635389435, "TGFB1", 2.269683222, "AKAP2", -4.041856426, "SFRP4", -5.640378342, 
                       "PCOLCE", 2.000280195, "ANKRD1", -0.728067557, "LRRC15", -4.242573529, "WNT2", -8.226283667, 
                       "FGF1", -3.356119985, "OLFML2B", -1.760393004, "HMGA2", -2.200766566, "PROCR", 1.2849028, 
                       "TNC", 1.752573207, "AEBP1", -1.535094887, "TWIST1", -0.981666108, "EDNRA", -4.725167387, 
                       "COMP", -5.579854026, "DAB2", 0.566387224, "TMEM158", -0.566031985, "MMP3", -4.475255139, 
                       "AKAP12", -0.200232147, "CCL2", -1.076341517, "GAS1", -2.244825991, "NID2", -0.991067623, 
                       "WISP1", -6.472539191, "CFH", -1.835546882, "LHFP", 0.556957735, "TUBB6", 4.385684555, 
                       "DDR2", -0.499571201, "TIMP3", 2.708720847, "GLI2", -4.162141483, "LGALS1", 8.131331323, 
                       "PAPPA", -2.864530248, "WNT5B", -1.25325296, "AKT3", -0.435085869, "ITGA5", 2.436332189, 
                       "TMEFF1", -7.29831855, "MAP1B", 0.793154502, "THBS2", -1.526213764, "TGFB1I1", 0.08152135, 
                       "SERPINB2", -3.387563228, "SRPX", -0.303884222, "C1S", 0.657286522, "FOSL1", 1.921391474, 
                       "IGFBP3", 4.107856925, "RECK", -1.055382469, "GJA1", -0.027502484, "COL5A2", 1.54700311, 
                       "SPOCK1", 0.176900432, "MFAP5", -3.044749983, "BGN", 0.689247627, "FHL1", -0.623010179, 
                       "HTRA1", 1.651777818, "RCN3", 1.102992153, "COL5A1", 0.691017249, "ZEB1", -2.173006491, 
                       "FST", -0.616215239, "MYL9", 3.500876038, "EMP3", 1.970144774, "MMP14", 2.333526183, 
                       "LUM", -2.019433847, "ADAM12", -1.893251343, "CALD1", 3.105059652, "TPM2", 3.424432511, 
                       "COL1A1", 3.774726107, "PDGFRB", -3.35270898, "ACTA2", 0.731992976, "FOXC2", -5.204224579, 
                       "AXL", 0.827110047, "COL6A1", 3.044217945, "VCAN", -1.246998063, "CDH2", 0.167713451, 
                       "DCN", -4.36665853, "COL6A3", -1.20107819, "THY1", -2.3790873, "SERPINE2", 2.164131944, 
                       "SULF1", -1.634809894, "PRRX1", -3.771012534, "PDGFRA", -4.811990706, "FSTL1", 2.856280855, 
                       "MMP2", 0.203812802, "FBN1", 0.049060769, "LOXL2", 2.302632429, "COL6A2", 1.945111732, 
                       "LOX", 0.585673513, "CDH11", -2.867068173, "FAP", -3.173643284, "FN1", 4.795853899, 
                       "SERPINE1", 2.518784945, "SPARC", 2.892494616, "TWIST2", -4.787716139, "GREM1", -1.869417265, 
                       "TNFAIP6", -6.164812282, "VIM", 3.834933599, "COL3A1", 0.039456497, "POSTN", -1.952092743), ncol = 2, byrow = TRUE)

Msig <- data.frame(Genes = tmp[,1], NormFactor = as.numeric(tmp[,2]))

BC.RNAseq.M <- BC.RNAseq.dat[match(Msig$Genes, BC.RNAseq.dat$Description), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
log_MBrCa <- log2(apply(BC.RNAseq.M[,c(3:ncol(BC.RNAseq.M))], 2, mean) + 0.0001)
MBrCa <- apply(BC.RNAseq.M[,c(3:ncol(BC.RNAseq.M))]/2^Msig$NormFactor, 2, mean)

BCType <- BrCa_Types$Type[match(colnames(BC.RNAseq.M)[c(3:ncol(BC.RNAseq.M))], BrCa_Types$CellName)]
BrCa_State <- data.frame(CellLineName = colnames(BC.RNAseq.M)[c(3:ncol(BC.RNAseq.M))], Epithelial = EBrCa, Mesenchymal = MBrCa, Intrinsic_Subtype = BCType)
write.csv(BrCa_State, file = "BrCa-EMT-State.csv")

BrCa_EMT <- read.csv(file = "BrCa-EMT-State.csv")
##################################################################
#
#
#
##################################################################
#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
DotColors <- c("red", "yellow", "pink", "blue", "black")
ColorDot <- DotColors[as.factor(BrCa_EMT$Intrinsic_Subtype)]
SlimCL <- lapply(as.character(BrCa_EMT$CellLineName), function(x) substr(x, 1, nchar(x)-7))

pdf("PCA-analysis-BrCa-CellLinesProjection-March2019.pdf", width = 25, height = 25)
plot(log2(BrCa_EMT$Mesenchymal + 0.0001), log2(BrCa_EMT$Epithelial + 0.0001), type = "p", pch = 19, col = ColorDot, 
     xlab = "Mesenchymal Signal", ylab = "Epithelial Signal", xlim = c(-2, 10), ylim = c(-6,6))

legend(x="bottomleft", legend = c("Basal", "Claudin Low", "HER2", "Luminal A", "Luminal B"), col = c("red", "yellow", "pink", "blue", "black"), pch = 19, cex = 1)
text(log2(BrCa_EMT$Mesenchymal + 0.0001), log2(BrCa_EMT$Epithelial + 0.0001)+0.1, labels = SlimCL, cex = 0.5)
dev.off() #

jnk <- as.character(BrCa_EMT$CellLineName)

#Multidimensional Scaling Approach
# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- dist(PCA.data) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results


library(colorspace)
library(MASS)
setwd("~/Documents/Projects/CCLE/R")
rm(list = ls())

#Reads are in TPM, so don't need to normalize to gene length as it's already done
RNAseq.file.name <- "./data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
RNAseq.dat <- read.table(RNAseq.file.name, head=TRUE, sep = "\t", skip = 0, na.strings = "null", stringsAsFactors = FALSE)

RNAseq.file.name2 <- "./data/GSE118389_tpm_rsem.txt"
RNAseq.dat2 <- read.table(RNAseq.file.name2, head=TRUE, sep = "\t", skip = 0, na.strings = "null", stringsAsFactors = FALSE)
RNAseq.dat2a <- RNAseq.dat2["WISP1",]
write.csv(RNAseq.dat2a, file = "scBrCa-WISP1-April19-2.csv")

SampleInfo.file.name <- "./data/GSE75688_final_sample_information.txt"
Info.dat <- read.table(SampleInfo.file.name, head=TRUE, sep = "\t", skip = 0, na.strings = "null", stringsAsFactors = FALSE)
KeepSample <- as.character(Info.dat$sample[Info.dat$type == "SC" & Info.dat$index == "Tumor"])
KeepSampleB <- as.character(Info.dat$sample[Info.dat$type == "Bulk" & Info.dat$index == "Tumor"])
#Header.dat <- scan(CNV.file.name, sep = "\t", what = "character", nlines = 2)

scRNAseq.dat <- RNAseq.dat[,c("gene_id", "gene_name", "gene_type", KeepSample)]
scRNAseq.dat2 <- RNAseq.dat[,c("gene_id", "gene_name", "gene_type", KeepSampleB)]
test1 <- RNAseq.dat[match("WISP1", RNAseq.dat$gene_name), ]
write.csv(test1, file = "scBrCa-WISP1-April19.csv")

KeepSampleC <- which(test1 > 0)

test2 <- RNAseq.dat[match("WISP1", RNAseq.dat$gene_name), c("gene_id", "gene_name", "gene_type", KeepSampleC)]

tmp <- read.csv("./data/hsa-GeneAnnotation.csv")
Genelength <- tmp[tmp$Transcript.type == "protein_coding",]
##################################################################
#
# EMT signature correlates
#
##################################################################

#Epithelial signature
# Changes SEPP1:SELENOP, C1orf106:INAVA
tmp <- matrix(data = c("CYP4B1", -5, "ARHGAP8", -5, "ALDH3B2", 0, 
                       "FA2H", 0, "BSPRY", 0, "F11R", -5, "GRHL2", 0, 
                       "OVOL2", 0, "CKMT1A", 0, "BLNK", 0, 
                       "FXYD3", 0, "FOXA1", 0, "RAB25", 0, "FBP1", 0, 
                       "TMPRSS2", 0, "WNT3A", -3, "ESRP1", 0, "EPN3", 0, 
                       "TUBBP5", 0, "PRSS8", 0, "EVPL", 0, "S100A14", 2, 
                       "C4orf19", -4, "OR7E14P", -1, "BIK", 0, "AGR2", 0, 
                       "GRB7", 0, "ATP2C2", 0, "EPHA1", 0, "LAD1", 0, 
                       "ST14", 0, "MYO5C", 0, "ESRP2", 1, "PDGFB", 0, 
                       "LLGL2", 3, "AP1M2", 2, "TTC39A", 0, "SPINT1", 2, 
                       "ICA1", 0, "SCNN1A", 0, "ELF3", 2, "MYH14", 0, 
                       "POF1B", -3, "IRF6", 1, "CLDN4", 2, "EPCAM", 2, 
                       "SEPP1", 0, "EHF", 0, "CLDN7", 0, "HOXC13", 0, 
                       "VAV3", 0, "DSC2", 0, "TSPAN15", 2, "SLC37A1", 0, 
                       "CDH1", 1, "ERBB3", 2, "IL1RN", 0, "VAMP8", 2, 
                       "CX3CR1", 0, "PTK6", 0, "EDN2", 0, "CGN", 0, 
                       "TMC6", 0, "CEACAM6", 0, "SPINT2", 5, "MSX2", 0, 
                       "EFNA1", 0, "MAP7", 0, "C1orf106", 0, "TSPAN1", 0, 
                       "SHANK2", 0, "ANXA9", 0, "WNT4", 0, "B3GAT1", -5, 
                       "BMP7", 0, "CDS1", 0, "PPL", 0, "MST1R", 0, 
                       "SORL1", 0, "IL20RA", 0, "ANK3", 0, "KRT8", 5, 
                       "GADD45G", 0, "DENND2D", 0, "EPS8L1", 0, "JUP", 3, 
                       "CNKSR1", 0, "SH2D3A", 0, "ARAP2", -1, "CXCR4", -1, 
                       "HPGD", -1, "LSR", 2, "EXPH5", 0, "CEACAM1", 0, 
                       "GALNT3", 1, "RBM47", 1, "MYB", 0, "PKP3", 1, 
                       "WNT7B", 1, "OCLN", 0, "GPX2", 0), ncol = 2, byrow = TRUE)
Esig <- data.frame(Genes = tmp[,1], KD = as.numeric(tmp[,2]))

BC.RNAseq.E <- scRNAseq.dat[match(Esig$Genes, scRNAseq.dat$gene_name), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
#EBrCa <- apply(BC.RNAseq.E[,c(4:ncol(BC.RNAseq.E))]/(Egl * 2^Esig$NormFactor), 2, function(x) mean(x, na.rm = TRUE))
EBrCa <- apply(BC.RNAseq.E[,c(4:ncol(BC.RNAseq.E))]/(BC.RNAseq.E[,c(4:ncol(BC.RNAseq.E))] + 2^Esig$KD), 2, function(x) mean(x, na.rm = TRUE))

# Mesenchymal signature
# change C7orf10:SUGCT, LEPRE1:P3H1, LHFP:LHFPL6
tmp <- matrix(data = c("SACS", 1, "SH3KBP1", 1, "BAG2", 0, "CLIC4", 5, 
                       "C7orf10", -1, "CXCL3", -2, "S100A4", 4, "SDC2", 1, 
                       "PHLDA1", 2, "WNT5A", 0, "CTSB", 5, "CYBRD1", 0, 
                       "GLT8D2", 0, "TCF4", -1, "TRPC1", 0, "PITX2", -1, 
                       "PLAUR", 1, "MXRA7", 1, "LEPRE1", 4, "TFPI", 0, 
                       "ASPN", -2, "PMP22", 2, "GFPT2", 0, "PDGFC", 0, 
                       "ITGB1", 7, "IFITM3", 6, "GEM", 2, "COPZ2", 0, 
                       "VEGFC", 1, "CD68", 3, "MME", 0, "ANK2", -3, 
                       "SMARCA1", 2, "TGFB1", 2, "AKAP2", -1, "SFRP4", -1, 
                       "PCOLCE", 4, "ANKRD1", 0, "LRRC15", -3, "WNT2", -5, 
                       "FGF1", 0, "OLFML2B", 0, "HMGA2", 0, "PROCR", 2, 
                       "TNC", 3, "AEBP1", 0, "TWIST1", 0, "EDNRA", -1, 
                       "COMP", -1, "DAB2", 2, "TMEM158", 0, "MMP3", -2, 
                       "AKAP12", 1, "CCL2", 0, "GAS1", 0, "NID2", 0, 
                       "WISP1", 0, "CFH", 0, "LHFP", 2, "TUBB6", 4, 
                       "DDR2", 3, "TIMP3", 4, "GLI2", -1, "LGALS1", 8, 
                       "PAPPA", -1, "WNT5B", 0, "AKT3", 0, "ITGA5", 3, 
                       "TMEFF1", -3, "MAP1B", 2, "THBS2", 0, "TGFB1I1", 1, 
                       "SERPINB2", 0, "SRPX", 2, "C1S", 3, "FOSL1", 3, 
                       "IGFBP3", 5, "RECK", 0, "GJA1", 2, "COL5A2", 3, 
                       "SPOCK1", 2, "MFAP5", 0, "BGN", 3, "FHL1", 0, 
                       "HTRA1", 3, "RCN3", 2, "COL5A1", 3, "ZEB1", 0, 
                       "FST", 1, "MYL9", 2, "EMP3", 3, "MMP14", 3, 
                       "LUM", 2, "ADAM12", 1, "CALD1", 2, "TPM2", 3, 
                       "COL1A1", 2, "PDGFRB", 0, "ACTA2", 2, "FOXC2", -1, 
                       "AXL", 1, "COL6A1", 3, "VCAN", 1, "CDH2", 2, 
                       "DCN", 0, "COL6A3", 1, "THY1", 1, "SERPINE2", 3, 
                       "SULF1", 2, "PRRX1", 0, "PDGFRA", 0, "FSTL1", 3, 
                       "MMP2", 2, "FBN1", 2, "LOXL2", 3, "COL6A2", 4, 
                       "LOX", 4, "CDH11", 1, "FAP", 1, "FN1", 5, 
                       "SERPINE1", 5, "SPARC", 5, "TWIST2", 0, "GREM1", 2, 
                       "TNFAIP6", 0, "VIM", 5, "COL3A1", 3, "POSTN", 2), ncol = 2, byrow = TRUE)
Msig <- data.frame(Genes = as.character(tmp[,1]), KD = as.numeric(tmp[,2]))

BC.RNAseq.M <- scRNAseq.dat[match(Msig$Genes, scRNAseq.dat$gene_name), ]
write.csv(BC.RNAseq.M, file = "scBrCa-Genes-April19.csv")
#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
#MBrCa <- apply(BC.RNAseq.M[,c(4:ncol(BC.RNAseq.M))]/(Mgl * 2^Msig$NormFactor), 2, function(x) mean(x, na.rm = TRUE))
MBrCa <- apply((BC.RNAseq.M[,c(4:ncol(BC.RNAseq.M))])/(BC.RNAseq.M[,c(4:ncol(BC.RNAseq.M))] + 2^Msig$KD), 2, function(x) mean(x, na.rm = TRUE))

BrCa_State <- data.frame(CellName = colnames(BC.RNAseq.M)[c(4:ncol(BC.RNAseq.M))], Epithelial = EBrCa, Mesenchymal = MBrCa)
write.csv(BrCa_State, file = "scBrCa-EMT-State-April4.csv")

BrCa_EMT <- read.csv(file = "scBrCa-EMT-State-April4.csv")
##################################################################
#
#
#
##################################################################
#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
DotColors <- c("blue", "blue", "black", "black", "pink", "pink", "pink", "red", "red", "red", "red", "red")
DotShape <- c(21, 22, 21, 22, 21, 22, 23, 21, 22, 23, 24, 25)
SlimCL <- sapply(as.character(BrCa_EMT$CellName), function(x) substr(x, 1, nchar(x)-3))
ColorDot <- DotColors[as.factor(SlimCL)]
ShapeDot <- DotShape[as.factor(SlimCL)]

xvals <- 10^seq(-6,0,len = 1000)
yvals <- 1 - xvals

pdf("PCA-analysis-BrCa-scProjection-April4.pdf", width = 5, height = 5)
plot(log2(BrCa_EMT$Mesenchymal + 0.0001), log2(BrCa_EMT$Epithelial + 0.0001), type = "p", pch = ShapeDot, col = ColorDot, lwd = 2,
     xlab = "Mesenchymal Signal", ylab = "Epithelial Signal", xlim = c(-6, 0), ylim = c(-6,0))
lines(log2(xvals+0.0001), log2(yvals+0.0001), lty = 2, col = "black")

legend(x="bottomleft", legend = names(summary(as.factor(SlimCL))), col = DotColors, pch = DotShape, cex = 0.7)
#text(log2(BrCa_EMT$Mesenchymal + 0.0001), log2(BrCa_EMT$Epithelial + 0.0001)+0.1, labels = SlimCL, cex = 0.5)
dev.off() #


# # List from Tan et al. EMBO J 2014
# #genes not present: TMEM30B, F11R, C19orf21, C10orf116, C12orf24, UGT1A1, MIR200A, MIR200C, MIR429, MIR141, MIR34A, CTLA4, EPYC, FGF4, HNF1A, MYCBP, WNT1, WNT8A, WNT9B, PLAY, PLOCE, NOXA, HNT, TMEFF1, AEBP12, 
# EMTsigA <- c("WISP1", "CDH1", "AGR2", "EPCAM", "KRT19", "RAB25", "TACSTD2", "S100P", "CEACAM6", "GALNT3", "FXYD3", "SPINT2", "SCNN1A", 
#               "ST14", "ESRP1", "S100A14", "CLDN7", "ERBB3", "RBM47", "SPINT1", "ELF3", "CLDN4", "PRSS8", "SH3YL1", "EHF", "LCN2", "JUP", 
#               "VAMP8", "KRT8", "C1orf106", "KRT7", "DSP", "CDS1", "ITGB4", "TMPRSS4", "LSR", "SORL1", "GRHL2", "PPL", "C1orf116", "TSPAN1", 
#               "MAP7", "SLPI", "TOX3", "ARHGAP8", "LAD1", "GPX2", "CTSH", "GPR56", "FA2H", "KLF5", "AREG", "KRT18", "SCEL", "CDH3", 
#               "MPZL2", "AIM1", "OVOL2", "LLGL2", "ESRP2", "MYO5C", "DDR1", "VGLL1", "IRF6", "SFN", "TSPAN13", "KCNK1", "MYO1D", 
#               "PKP3", "ITGB6", "LY75", "MAPK13", "TTC39A", "ELMO3", "CEACAM1", "DTX4", "ERBB2", "RAB11FIP1", "ATP2C2", "MST1R", "AP1M2", 
#               "TGFA", "MYO6", "PTK6", "OAS1", "FBP1", "AQP3", "CBLC", "EPHA1", "BSPRY", "SH2D3A", "EPS8L1", "GRB7", "C4orf19", "KLK6", 
#               "TJP2", "PLS1", "DENND2D", "EPS8L2", "IL20RA", "HES1", "IL1RN", "EXPH5", "ARHGDIB", "CAMK2N1", "HPGD", "SYNGR2", 
#               "PERP", "MANSC1", "DSC2", "POF1B", "SERINC5", "BIK", "ANXA9", "MALL", "EPN3", "STAP2", "FOXA1", "PYCARD", "ZNF165", 
#               "SLC37A1", "ANK3", "TSPAN15", "HNMT", "ABCC3", "SDC1", "CKMT1A", "TOB1", "B3GNT3", "TMC6", "CD9", "ADAP1", "ATP1B1", "SHANK2", 
#               "CYB561", "ERMP1", "RAB20", "MYH14", "CAPN1", "ALDH3B2", "TRIM31", "ARAP2", "SSH3", "ICA1", "ARHGEF5", "ALOX5", "RHOD", "TMPRSS2", 
#               "MTUS1", "CYP4F3", "PPFIBP2", "RABGAP1L", "PLXNB2", "MGST2", "OR7E14P", "EVPL", "CD46", "KRT15", "CNKSR1", "BLNK", "COMT", "ANXA4", 
#               "TNFSF13", "OCLN", "SLC9A3R1", "XBP1", "PTRF", "RECK", "GFPT2", "FSTL1", "MXRA7", "MYL9", "BAG2", "KDELC1", "TRPC1", 
#               "POPDC3", "CEP170", "LHFP", "COL5A2", "SOAT1", "SERPINE1", "TGFB1I1", "AKAP12", "SOBP", "ETV1", "LEPRE1", "CHN1", "SH2B3", "ANK2", 
#               "SRPX", "TUBB6", "TPM2", "DENND5A", "GJA1", "AP1S2", "MAP1B", "GLYR1", "MSN", "PMP22", "LGALS1", "CALD1", "TMEM158", 
#               "TUBA1A", "FERMT2", "FHL1", "SPARC", "LOXL2", "AXL", "SACS", "EMP3", "ZEB1", "VIM")
# 
# # List from Kaiser et al. Biotech Prog 2016
# EMTsigB <- c("ACTA2", "ADAM12", "ASPN", "AXIN2", "BGN", "C1QTNF3", "C7orf10", "CDH11", "CLDN1", "COL10A1", "COL11A1", "COL1A1", "COL3A1", 
#              "COL5A1", "COL6A2", "COL6A3", "COMP", "COPZ2", "CRISPLD2", "CTSK", "DCN", "EDNRA", "EN2", "FAP", "FBN1", "FGF18", 
#              "FN1", "FOSL1", "FST", "FZD7", "GLT8D2", "GREM1", "ID2", "INHBA", "ITGBL1", "JUN", "LEF1", "LGR5", "LOX", "LRRC15", 
#              "LUM", "MFAP5", "MITF", "MMP11", "MMP2", "MMP7", "MXRA5", "MYC", "NID2", "NRCAM", "NUAK1", "OLFML2B", "PDGFRB", 
#              "PITX2", "POSTN", "PPARD", "PRRX1", "RAB31", "RCN3", "RHOU", "SERPINF1", "SFRP4", "SNAI2", "SP5", "SPOCK1", "STRA6", 
#              "SULF1", "TCF4", "THBS2", "THY1", "TIMP3", "TNFAIP6", "TNFRSF19", "VEGFA", "WNT10A", "WNT10B", "WNT11", "WNT16", "WNT2", 
#              "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8B", "WNT9A") 
# 
# EMTsig <- c(EMTsigA, EMTsigB)

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

#just first time through
FirstTime <- TRUE
if(FirstTime) FtestInit <- rep(0, length(EMTsig)) else FtestDist <- rep(0, length(EMTsig))
colR <- rep("red", length(EMTsig))

pdf("BrCa-EMT-Genes-Big.pdf", width = 7, height = 7)
#summary(BC.RNAseq.dat$Description == "MIR34A") #"Mir200", "Mir34"
opar <- par(mfrow = c(3,3))
for (i in 1:length(EMTsig))
{
  xval <- as.numeric(BC.RNAseq.EMT[BC.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(BC.RNAseq.EMT))])
  yval <- as.numeric(BC.RNAseq.EMT[match(EMTsig[i], BC.RNAseq.EMT$Description), c(3:ncol(BC.RNAseq.EMT))])
  WH <- xval > 0.81
  WLGH <- length(yval[yval > 0.81 & !WH])
  WLGL <- length(yval[yval <= 0.81 & !WH])
  WHGH <- length(yval[yval > 0.81 & WH])
  WHGL <- length(yval[yval <= 0.81 & WH])
  
  Ftest <- fisher.test(matrix(c(WLGH, WHGH, WLGL, WHGL), nrow = 2))
  if(FirstTime) FtestInit[i] <- Ftest$p.value else FtestDist[i] <- Ftest$p.value
  colR[i] <- ifelse(Ftest$p.value < 0.05*length(FtestInit[FtestInit <= Ftest$p.value])/length(FtestInit), "black", "red") #Benjamini-Hochberg FDR correction
  plot(log2(xval + 0.001), log2(yval + 0.001), xlim = c(-10,6), ylim = c(-10,10), xlab = "WISP1 (log2 RPKM)", 
       main = list(cex = 0.8, col = colR[i], paste("EMT:", EMTsig[i], " p-value =", sprintf("%6.4f", Ftest$p.value))))
  #,"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.81, 0.81)), c(-10,14), col = "blue", lty = 2)
  lines(c(-10,6), log2(c(0.81, 0.81)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

FtestInit <- sort(FtestInit, decreasing = FALSE)
FirstTime <- FALSE
outGCA <- data.frame(Genename = EMTsig, BHPvalue = FtestDist, Sig = colR)
write.csv(outGCA, file = "BrCa-WISP1GeneRelation.csv")
dev.off() # End plot 


xval <- as.numeric(BC.RNAseq.EMT[BC.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(BC.RNAseq.EMT))])
yval <- as.numeric(BC.RNAseq.EMT[BC.RNAseq.EMT$Description == "CDH1", c(3:ncol(BC.RNAseq.EMT))]) + as.numeric(BC.RNAseq.EMT[BC.RNAseq.EMT$Description == "VIM", c(3:ncol(BC.RNAseq.EMT))])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), ylim = c(-13,15), main = paste("EMT:", EMTsig[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
lines(log2(c(0.81, 0.81)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.81, 0.81)), col = "blue", lty = 2)

#Convert from RPKM to log2(RPKM + 0.0001)
SeqRes <- log2(BC.RNAseq.EMT[,c(3:ncol(BC.RNAseq.EMT))] + 0.0001)

# # Principal component analysis
# opar <- par(mfrow = c(2,2))
# # Number of genes x Number of Cell lines
# # t(SeqRes) gives PCA clustering of cell lines
# fs1PCA <- prcomp(t(SeqRes), retx = TRUE, scale. = FALSE)
# 
# #Scree plot
# EigenVals <- fs1PCA$sdev^2
# TotVar <- sum(EigenVals)
# PrPC <- c(1:10)
# #PrPC[1] <- EigenVals[1]/
# for ( i in 1:10 )
# {
#   PrPC[i] <- EigenVals[i]/TotVar
# }
# 
# load.Matrix <- data.frame(fs1PCA$rotation)
# # n_genes x n_PC          =   n_genes x n_pts      *  n_pts x n_PC
# # (square n_genes = n_PC) = (here n_genes < n_pts) * (n_pts > n_PC)
# PCA.data <- as.matrix(t(SeqRes)) %*% fs1PCA$rotation
# 
# pdf("PCA-analysis-BrCa2.pdf", width = 7, height = 7)
# opar <- par(mfrow = c(2,3))
# mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,60), ylab = "Variance (%)", xlab = "Principal Component")
# 
# CDH1 <- as.numeric(BC.RNAseq.EMT[BC.RNAseq.EMT$Description == "CDH1",c(3:ncol(BC.RNAseq.EMT))])
# VIM <- as.numeric(BC.RNAseq.EMT[BC.RNAseq.EMT$Description == "VIM",c(3:ncol(BC.RNAseq.EMT))])
# EMTsignal = 0.5*(1 - (CDH1 - 0.3)/200 + (VIM - 0.6)/1000)
# rainColor <- colorRamp(c("red", "blue"), space = "rgb", interpolate = "linear")
# rainColor <- colorRampPalette(c("red", "yellow", "blue"))(100)
# PtCol <- rainColor[cut(EMTsignal,breaks = 100)]
# 
# plot(EMTsignal, rep(1, length(EMTsignal)), pch = 19, col = PtCol)
# 
# cbind(CDH1, VIM, EMTsignal)
# 
# opar <- par(mfrow = c(1,1))
# plot(PCA.data[,c(1,2)], type = "n", xlim = c(-60, 60), ylim = c(-35, 30))
# points(PCA.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
# legend("topright",legend=c("Epithelial (CDH1+ VIM-)", "Mesenchymal (CDH1- VIM+)"),col =rainColor[c(1,100)],pch=20)
# 
# dev.off() # End plot 

# Number of Cell lines x Number of genes
# SeqRes gives PCA clustering of genes
# t(SeqRes) gives PCA clustering of cell lines
fs2PCA <- prcomp(SeqRes, retx = TRUE, scale. = FALSE)

#Scree plot
EigenVals <- fs2PCA$sdev^2
TotVar <- sum(EigenVals)
PrPC <- c(1:10)
#PrPC[1] <- EigenVals[1]/
for ( i in 1:10 )
{
  PrPC[i] <- EigenVals[i]/TotVar
}

load.Matrix <- data.frame(fs2PCA$rotation)
# n_genes x n_PC          =   n_genes x n_pts      *  n_pts x n_PC
# (square n_genes = n_PC) = (here n_genes < n_pts) * (n_pts > n_PC)
PCA2.data <- as.matrix(SeqRes) %*% fs2PCA$rotation

pdf("PCA-analysis-BrCa-Genes5.pdf", width = 7, height = 7)
opar <- par(mfrow = c(2,2))
mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,70), ylab = "Variance (%)", xlab = "Principal Component")
lines(c(0,2,3,4,5,6,7,8,9,12), c(3.4, 3.2, 3.19, 3.04, 3.01, 2.89, 2.8, 2.74, 2.67, 2.63), col="red", lwd = 2, lty = 2)

# Keep gene for further analysis if it is expressed above the threshold in more than 
# 5% of the samples
GeneExpressed <- apply(SeqRes, 1, function(x) length(x[x > -0.304])/length(x) > 0.05 )
rainColor <- colorRampPalette(c("red", "yellow", "blue"))(100)
meanSeqRes <- rowMeans(SeqRes)
#PtCol <- rainColor[cut(meanSeqRes,breaks = 100)]
#PtCol <- ifelse(meanSeqRes > -0.304, "black", "red")
PtCol <- ifelse(GeneExpressed, "black", "red")

opar <- par(mfrow = c(1,1))
plot(-1*PCA2.data[,1], PCA2.data[,2],type = "n", xlim = c(-110, 80), ylim = c(-40, 40))
#points(PCA2.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
text(-1*PCA2.data[,1], PCA2.data[,2], labels = BC.RNAseq.EMT$Description, cex = 1, col = PtCol)
points(-1*PCA2.data[BC.RNAseq.EMT$Description == "CDH1",1], PCA2.data[BC.RNAseq.EMT$Description == "CDH1",2], pch = 19, col = "blue", cex = 1)
points(-1*PCA2.data[BC.RNAseq.EMT$Description == "VIM",1], PCA2.data[BC.RNAseq.EMT$Description == "VIM",2], pch = 19, col = "yellow", cex = 1)
lines(c(-100, 70), c(11.9, 11.9), col="red", lwd=2, lty=2)
lines(c(-100, 70), c(-11.9, -11.9), col="red", lwd=2, lty=2)
#lines(c(30, 30), c(-100, 70), col="blue", lwd=2, lty=2)

KeepGene <- GeneExpressed
# KeepGene <- meanSeqRes > -0.304
# KeepGene <- PCA2.data[,1] < 30 
plot(PCA2.data[KeepGene,c(3,4)], type = "n", xlim = c(-30, 30), ylim = c(-30, 30))
#points(PCA2.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
text(PCA2.data[KeepGene,c(3,4)], labels = BC.RNAseq.EMT$Description[KeepGene], cex = 0.5, col = "black" ) #PtCol[KeepGene]
points(PCA2.data[BC.RNAseq.EMT$Description == "CDH1",3], PCA2.data[BC.RNAseq.EMT$Description == "CDH1",4], pch = 19, col = "blue", cex = 1)
points(PCA2.data[BC.RNAseq.EMT$Description == "VIM",3], PCA2.data[BC.RNAseq.EMT$Description == "VIM",4], pch = 19, col = "yellow", cex = 1)
lines(c(-110, 80), c(11.9, 11.9), col="red", lwd=2, lty=2)
lines(c(-110, 80), c(-11.9, -11.9), col="red", lwd=2, lty=2)
lines(c(-11.9, -11.9), c(-100, 70), col="blue", lwd=2, lty=2)
lines(c(11.9, 11.9), c(-100, 70), col="blue", lwd=2, lty=2)

dev.off() # End plot 

outPCA <- data.frame(name = BC.RNAseq.EMT$Description, GeneMean = meanSeqRes, KeepGene = GeneExpressed, PC1 = PCA2.data[,1], PC2 = PCA2.data[,2], PC3 = PCA2.data[,3])
write.csv(outPCA, file = "GenePCAresult4.csv")
plot(density(PCA2.data[,3], adj = 0.25))
L1 <- density(c(NCx, NCy), adjust = 0.25, na.rm=TRUE, from = -30, to = 30)
lines(L1, col = "red")

Delta <- 0.5
x <- seq(-30,20, by=Delta)
# Here P is the same as dnorm (probability density function for normal
# distribution), but other functions could be tried here.
P <- function(x, mean, sd)
{
  variance <- sd^2
  exp(-(x-mean)^2/(2*variance)) / sqrt(2*pi*variance)
}
KeepGene2 <- PCA2.data[,3] < 5 & PCA2.data[,3] > -10
y <- P(x, median(PCA2.data[KeepGene2,3]), sd(PCA2.data[KeepGene2,3]))/1.2
lines(x, y, col = "red")
#####################################################
#Summary results
Ntot = 50
NCx = vector('double',783*Ntot)
NCy = vector('double',783*Ntot)

for ( i in 1:Ntot )
{
  #Negative control 
  RS <- sample(simplify2array(SeqRes), 783*57, replace=T)
  A2 <- matrix(c(RS), nrow = 783)
  Rdat <- SeqRes
  Rdat[c(1:783),c(1:57)] <- A2

  NC_PCA <- prcomp(Rdat, retx = TRUE, scale. = FALSE)
  PCA2.data <- as.matrix(Rdat) %*% NC_PCA$rotation
  NCx[783*(i-1)+c(1:783)] <- PCA2.data[,1]
  NCy[783*(i-1)+c(1:783)] <- PCA2.data[,2]
}

#Scree plot
EigenVals <- NC_PCA$sdev^2
TotVar <- sum(EigenVals)
PrPC <- c(1:10)
#PrPC[1] <- EigenVals[1]/
for ( i in 1:10 )
{
  PrPC[i] <- EigenVals[i]/TotVar
}
# Results for PrPC
#  [1] 0.03444504 0.03232179 0.03194559 0.03037966 0.03014023 0.02890219 0.02836817 0.02742888 0.02667278
# [10] 0.02630945

mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,60), ylab = "Variance (%)", xlab = "Principal Component")

NC2d <- kde2d(NCx, NCy, n = 25)
contour(NC2d, xlab = "X", ylab = "Y", levels = c(0.00001, 0.0001, 0.001, 0.002))
image(NC2d, zlim = c(0, 0.0001))

contour(NC2d, levels=c(3.7e-01), add=TRUE, col="red", labels="P1", lwd=3)
N.total <- sum(NC2d$z)
100*sum(NC2d$z[NC2d$z > 3.7e-01]) / N.total #95%

L1 <- density(c(NCx, NCy), adjust = 0.25, na.rm=TRUE, from = -30, to = 30)
N.total <- sum(L1$y)
100*sum(L1$y[L1$y > 0.167e-1]) / N.total #95%
L1$x[L1$y < 0.167e-1]
#NCx NCy 95% Limits: -11.8 11.9
#NCx NCy 90% Limits: -10.2 9.8

plot(L1)

##################################################################
#
# EMT signature correlates
#
##################################################################

#Old Epithelial signature
 # Esig <- c("CYP4B1", "ARHGAP8", "ALDH3B2", "FA2H", "BSPRY", "F11R", "GRHL2", "OVOL2", "CKMT1A", 
 #          "BLNK", "FXYD3", "FOXA1", "RAB25", "FBP1", "TMPRSS2", "ESRP1", "EPN3", "TUBBP5", 
 #          "PRSS8", "EVPL", "S100A14", "C4orf19", "OR7E14P", "BIK", "AGR2", "GRB7", "ATP2C2", 
 #          "EPHA1", "LAD1", "ST14", "MYO5C", "ESRP2", "PDGFB", "LLGL2", "AP1M2", "TTC39A", 
 #          "SPINT1", "ICA1", "MYH14", "POF1B", "IRF6", "CLDN4", "EPCAM", "SEPP1", "EHF", 
 #          "CLDN7", "HOXC13", "VAV3", "TSPAN15", "SLC37A1", "CDH1", "ERBB3", "VAMP8", "PTK6", 
 #          "EDN2", "CGN", "TMC6", "CEACAM6", "SPINT2", "MSX2", "EFNA1", "MAP7", "C1orf106", 
 #          "TSPAN1", "SHANK2", "ANXA9", "WNT4", "BMP7", "CDS1", "PPL", "MST1R", "SORL1", 
 #          "IL20RA", "ANK3", "KRT8", "GADD45G", "DENND2D", "EPS8L1", "JUP", "CNKSR1", "SH2D3A", 
 #          "ARAP2", "CXCR4", "HPGD", "LSR", "EXPH5", "GALNT3", "RBM47", "MYB", "PKP3", "WNT7B", 
 #          "OCLN", "GPX2")

tmp <- matrix(data = c("CYP4B1", -2.446785084, "ARHGAP8", -2.26795788, "ALDH3B2", 1.435970192, 
                       "FA2H", -0.5285652, "BSPRY", 1.792107746, "F11R", -3.65098147, "GRHL2", 2.272778131, 
                       "OVOL2", 0.060031635, "CKMT1A", 1.912022287, "F11R", -6.502068123, "BLNK", -1.318945039, 
                       "FXYD3", 3.329906423, "FOXA1", 2.694653058, "RAB25", 3.799895055, "FBP1", 1.101376263, 
                       "TMPRSS2", -0.278656984, "WNT3A", -4.758677792, "ESRP1", 3.205173557, "EPN3", 2.478648163, 
                       "TUBBP5", -1.176766827, "PRSS8", 3.62237251, "EVPL", 2.310660775, "S100A14", 4.668678552, 
                       "C4orf19", -3.344019889, "OR7E14P", -0.610740935, "BIK", 1.175798652, "AGR2", 3.328402466, 
                       "GRB7", 2.928559204, "ATP2C2", 1.407978993, "EPHA1", 1.463343712, "LAD1", 2.832460693, 
                       "ST14", 3.857103038, "MYO5C", 1.827949913, "ESRP2", 3.149686814, "PDGFB", 1.634909829, 
                       "LLGL2", 3.503815163, "AP1M2", 4.059617591, "TTC39A", 1.093368035, "SPINT1", 4.28176267, 
                       "ICA1", 1.255665985, "SCNN1A", 2.045415355, "ELF3", 2.930432283, "MYH14", 2.611742248, 
                       "POF1B", -1.061477381, "IRF6", 2.512716819, "CLDN4", 4.439361084, "EPCAM", 5.035654841, 
                       "SEPP1", 0.99912883, "EHF", 2.218922409, "CLDN7", 4.421022987, "HOXC13", 2.023664943, 
                       "VAV3", 1.255238988, "DSC2", 1.588680098, "TSPAN15", 3.95722022, "SLC37A1", 1.632150536, 
                       "CDH1", 3.978846105, "ERBB3", 3.496068379, "IL1RN", -2.209158624, "VAMP8", 5.014175672, 
                       "CX3CR1", -5.547077396, "PTK6", 1.411144861, "EDN2", -0.110772193, "CGN", 2.421659511, 
                       "TMC6", 2.429819726, "CEACAM6", 0.388508527, "SPINT2", 6.814914712, "MSX2", 0.707879411, 
                       "EFNA1", 3.376318109, "MAP7", 2.769827993, "C1orf106", 0.724307937, "TSPAN1", 2.87527624, 
                       "SHANK2", 0.378890092, "ANXA9", 2.246682225, "WNT4", -1.178900939, "B3GAT1", -5.761322743, 
                       "BMP7", -0.636426934, "CDS1", 2.462592987, "PPL", 2.507467951, "MST1R", 1.579844441, 
                       "SORL1", 1.279696831, "IL20RA", -2.822257097, "ANK3", 0.481008272, "KRT8", 8.088252672, 
                       "GADD45G", 0.527432699, "DENND2D", 0.616714594, "EPS8L1", 2.770283683, "JUP", 5.202752406, 
                       "CNKSR1", 1.661710701, "SH2D3A", 2.276145054, "ARAP2", -0.499625902, "CXCR4", -0.115041413, 
                       "HPGD", -2.759089613, "LSR", 4.538511283, "EXPH5", -0.065494998, "CEACAM1", -0.99607925, 
                       "GALNT3", 2.556610315, "RBM47", 2.886676807, "MYB", 0.970250216, "PKP3", 3.954563786, 
                       "WNT7B", 1.672477704, "OCLN", 2.162243252, "GPX2", -0.062285999), ncol = 2, byrow = TRUE)
Esig <- data.frame(Genes = tmp[,1], NormFactor = as.numeric(tmp[,2]))

BC.RNAseq.E <- BC.RNAseq.dat[match(Esig$Genes, BC.RNAseq.dat$Description), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
log_EBrCa <- log2(apply(BC.RNAseq.E[,c(3:ncol(BC.RNAseq.E))], 2, mean) + 0.0001)
EBrCa <- apply(BC.RNAseq.E[,c(3:ncol(BC.RNAseq.E))]/2^Esig$NormFactor, 2, mean)

#Check for correlation between epidermal signature genes and others
#cor(cbind(EBrCa, t(log2(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))] + 0.0001))))

EBrCa_Cor <- apply(log2(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))] + 0.0001), 1, function(x) cor(x,EBrCa))
Epi_Cor <- data.frame(GeneName = BC.RNAseq.dat$Description, Corr = EBrCa_Cor)

Epi_Cor <- Epi_Cor[abs(Epi_Cor$Corr) > 0.7 & !is.na(Epi_Cor$Corr),]

plot(density(EBrCa_Cor, adj = 0.2, na.rm = TRUE))

write.csv(Epi_Cor, file = "BrCa-Epi-GeneCorrelation.csv")

#Old Mesenchymal signature
# Msig <- c("POSTN", "COL3A1", "VIM", "GREM1", "SPARC", "SERPINE1", "FN1", "CDH11", "LOX", 
#           "COL6A2", "FBN1", "MMP2", "SULF1", "SERPINE2", "THY1", "COL6A3", "CDH2", "VCAN", 
#           "COL6A1", "AXL", "ACTA2", "PDGFRB", "COL1A1", "TPM2", "ADAM12", "LUM", "EMP3", 
#           "MYL9", "FST", "ZEB1", "COL5A1", "RCN3", "HTRA1", "FHL1", "BGN", "MFAP5", "SPOCK1", 
#           "COL5A2", "GJA1", "RECK", "IGFBP3", "C1S", "SRPX", "TGFB1I1", "THBS2", "MAP1B", 
#           "ITGA5", "AKT3", "PAPPA", "LGALS1", "TIMP3", "DDR2", "TUBB6", "LHFP", "CFH", 
#           "NID2", "GAS1", "CCL2", "AKAP12", "TMEM158", "DAB2", "TWIST1", "AEBP1", "TNC", 
#           "PROCR", "HMGA2", "OLFML2B", "FGF1", "ANKRD1", "PCOLCE", "TGFB1", "SMARCA1", 
#           "ANK2", "MME", "CD68", "VEGFC", "COPZ2", "GEM", "ITGB1", "PDGFC", "GFPT2", 
#           "PMP22", "TFPI", "LEPRE1", "MXRA7", "PLAUR", "PITX2", "TRPC1", "TCF4", "GLT8D2", 
#           "CYBRD1", "CTSB", "WNT5A", "PHLDA1", "SDC2", "S100A4", "CXCL3", "C7orf10", 
#           "CLIC4", "BAG2", "SH3KBP1", "SACS")

tmp <- matrix(data = c("SACS", 0.290062684, "SH3KBP1", 1.466203798, "BAG2", 0.49636556, "CLIC4", 5.159359081, 
              "C7orf10", -2.51888532, "CXCL3", -3.423618245, "S100A4", 3.778864441, "SDC2", 1.237691836, 
              "PHLDA1", 2.195399561, "WNT5A", 0.111286371, "CTSB", 4.9741377, "CYBRD1", 1.525324739, 
              "GLT8D2", -1.47153745, "TCF4", -2.335071698, "TRPC1", -0.827329743, "PITX2", -3.782760584, 
              "PLAUR", 1.563004224, "MXRA7", 1.509509711, "LEPRE1", 2.79209761, "TFPI", 0.228503595, 
              "ASPN", -6.336878371, "PMP22", 2.066637029, "GFPT2", -0.323637449, "PDGFC", 0.450379972, 
              "ITGB1", 6.585941966, "IFITM3", 4.96794764, "GEM", -0.126991774, "COPZ2", -0.429158273, 
              "VEGFC", 0.509088251, "CD68", 2.005694138, "MME", -1.977989143, "ANK2", -3.891723139, 
              "SMARCA1", 1.403981074, "TGFB1", 1.995928221, "AKAP2", -4.329024029, "SFRP4", -5.868109735, 
              "PCOLCE", 1.747514167, "ANKRD1", -1.004379089, "LRRC15", -4.410930186, "WNT2", -8.414143146, 
              "FGF1", -3.573058761, "OLFML2B", -2.023325465, "HMGA2", -2.507777882, "PROCR", 1.002464008, 
              "TNC", 1.4768736, "AEBP1", -1.794438431, "TWIST1", -1.265599572, "EDNRA", -5.015670378, 
              "COMP", -5.895696801, "DAB2", 0.311142128, "TMEM158", -0.813402427, "MMP3", -4.774672962, 
              "AKAP12", -0.516119842, "CCL2", -1.388376675, "GAS1", -2.551990724, "NID2", -1.282312728, 
              "WISP1", -6.762209988, "CFH", -2.157779226, "LHFP", 0.23020315, "TUBB6", 4.06764218, 
              "DDR2", -0.75369598, "TIMP3", 2.43646936, "GLI2", -4.47772844, "LGALS1", 7.82056313, 
              "PAPPA", -3.189433741, "WNT5B", -1.566230649, "AKT3", -0.7568643, "ITGA5", 2.145434418, 
              "TMEFF1", -7.602217201, "MAP1B", 0.486147967, "THBS2", -1.820681347, "TGFB1I1", -0.249204358, 
              "SERPINB2", -3.771357136, "SRPX", -0.647854286, "C1S", 0.286522657, "FOSL1", 1.575046264, 
              "IGFBP3", 3.712917933, "RECK", -1.442349154, "GJA1", -0.346756644, "COL5A2", 1.216170119, 
              "SPOCK1", -0.183911903, "MFAP5", -3.363210987, "BGN", 0.316022331, "FHL1", -0.979230602, 
              "HTRA1", 1.268220144, "RCN3", 0.754204477, "COL5A1", 0.29846138, "ZEB1", -2.478131852, 
              "FST", -0.991784444, "MYL9", 3.176709968, "EMP3", 1.618830309, "MMP14", 1.941605043, 
              "LUM", -2.379093261, "ADAM12", -2.253021028, "CALD1", 2.704644634, "TPM2", 3.062659226, 
              "COL1A1", 3.467409961, "PDGFRB", -3.71542655, "ACTA2", 0.335708898, "FOXC2", -5.681890132, 
              "AXL", 0.449483087, "COL6A1", 2.630076832, "VCAN", -1.614362995, "CDH2", -0.225036023, 
              "DCN", -4.730347576, "COL6A3", -1.540084326, "THY1", -2.697321101, "SERPINE2", 1.735635925, 
              "SULF1", -2.00066345, "PRRX1", -4.144931635, "PDGFRA", -5.205499978, "FSTL1", 2.419411186, 
              "MMP2", -0.265964685, "FBN1", -0.322073696, "LOXL2", 1.876172565, "COL6A2", 1.499796853, 
              "LOX", 0.197294329, "CDH11", -3.293231657, "FAP", -3.665707164, "FN1", 4.354138118, 
              "SERPINE1", 2.091363265, "SPARC", 2.448393824, "TWIST2", -5.342141994, "GREM1", -2.307420799, 
              "TNFAIP6", -6.738435166, "VIM", 3.384074296, "COL3A1", -0.391285122, "POSTN", -2.389306523), ncol = 2, byrow = TRUE)
Msig <- data.frame(Genes = tmp[,1], NormFactor = as.numeric(tmp[,2]))

BC.RNAseq.M <- BC.RNAseq.dat[match(Msig$Genes, BC.RNAseq.dat$Description), ]

#Average across genes for cell line and convert from RPKM to log2(RPKM + 0.0001)
log_MBrCa <- log2(apply(BC.RNAseq.M[,c(3:ncol(BC.RNAseq.M))], 2, mean) + 0.0001)
MBrCa <- apply(BC.RNAseq.M[,c(3:ncol(BC.RNAseq.M))]/2^Msig$NormFactor, 2, mean)

#Check for correlation between mesenchymal signature genes and others
#cor(cbind(MBrCa, t(log2(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))] + 0.0001))))

MBrCa_Cor <- apply(log2(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))] + 0.0001), 1, function(x) cor(x,MBrCa))
Mesen_Cor <- data.frame(GeneName = BC.RNAseq.dat$Description, Corr = MBrCa_Cor)

Mesen_Cor <- Mesen_Cor[abs(Mesen_Cor$Corr) > 0.7 & !is.na(Mesen_Cor$Corr),]

plot(density(MBrCa_Cor, adj = 0.2, na.rm = TRUE))

write.csv(Mesen_Cor, file = "BrCa-Mesen-GeneCorrelation.csv")

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

pdf("PCA-analysis-BrCa-CellLinesProjection3.pdf", width = 25, height = 25)
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

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS",	type="p")
PtType <- unlist(lapply(SitePrimary.k1, function(x) ifelse(x == "skin" || x == "breast", 19, 1)))
points(x, y, pch = PtType, col = colorTest, cex = 1)
#text(x, y, labels = row.names(PCA.data.k1), cex=.7)

plot(x, y, ylim = c(0, 150), xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS",	type="n")
text(c(1), 5*c(1:length(levels(as.factor(SitePrimary.k1)))), labels = levels(as.factor(SitePrimary.k1)), col = colorMap[c(1:length(levels(as.factor(SitePrimary.k1))))])
#Take a smaller subset
# n_genes x n_PC          =   n_genes x n_pts      *  n_pts x n_PC
# (n_genes > 2 PCs) = (here n_genes < n_pts) * (n_pts > 2 PCs) * 2x2 rotation matrix
theta = 45*pi/180
rotMat <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE)
rotLine1 <- line1 %*% rotMat

plot(PCA.data[,c(2,3)]%*% rotMat, type = "n", xlim = c(-15, 10), ylim = c(-10, 10))
text(PCA.data[,c(2,3)]%*% rotMat, labels=rownames(PCA.data), cex = 0.5)
lines(rotLine1[,1], rotLine1[,2])

PCA.data2 <- as.matrix(tot.seqRes) %*% fs1PCA$rotation
plot(PCA.data2[,c(2,3)] %*% rotMat, type = "n", xlim = c(-22, 15), ylim = c(-15, 15))
text(PCA.data2[,c(2,3)] %*% rotMat, labels=rownames(PCA.data2), cex = 0.5)
lines(rotLine1[,1], rotLine1[,2])

rPCA.data2 <- PCA.data2[,c(2,3)] %*% rotMat
plot(density(rPCA.data2[,1], adj = 0.1))

write.csv(rPCA.data2, file = "RotatedGenePCA.csv")


EMTsigDN <- c("WISP1", "CAV2", "CDH1", "CDH3", "DSP", "FGFBP1", "IL1RN", "KRT19", "MITF", "MST1R", "NUDT13", "OCLN", "RGS2", "SPP1", 
              "TFPI2", "TSPAN13") # "PPPDE2", 

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigDN, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigDN))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigDN[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT DN:", EMTsigDN[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

EMTsigTF <- c("WISP1", "CTNNB1", "ESR1", "FOXC2", "GSC", "MITF", "MSX1", "NOTCH1", "SMAD2", "SNAI1", "SNAI2", "SNAI3", "SOX10", 
              "STAT3", "TCF3", "TCF4", "TWIST1", "ZEB1", "ZEB2") # "SIP1", 

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigTF, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigTF))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigTF[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT TF:", EMTsigTF[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}
#source("filename")
opar <- par(mfrow = c(1,1))
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH2", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNNB1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNND1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNNA1", c(3:ncol(MM.RNAseq.dat))]) 
xval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1"  , c(3:ncol(MM.RNAseq.dat))])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), ylim = c(5,10), main = "CDH1 + CDH2 + CTNNB1 + CTNND1 + CTNNA1", col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)

text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))

pCDH2 <- function(x, max = 128, kd = 1) 
{
  y = max * x^2/(x + kd)^2 + 0.06
  y
}

pCDH1 <- function(x, max = 128, kd = 1) 
{
  y = max * (1 - x^2/(x + kd)^2) + 0.06
  y
}

opar <- par(mfrow = c(1,1))
rainbowColor <- rainbow(n = 256, start = 0, end = 1)
rainbowColor <- heat.colors(n = 256)
rainbowColor <- colorRampPalette(c("red", "yellow", "green", "blue"))(256)
#WISP1min <- min(log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))]))) 
#WISP1max <- max(log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))]))) 
#WISP1shade <- ceiling(255*(1- (log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])) - WISP1min) /(WISP1max - WISP1min))) 

pdf("CDH1vsCDH2-BrCa.pdf", width = 7, height = 7)
#Just plot for breast cancer cell lines
SCL <- CellLines[CellLines$Site.Primary %in% c("breast"),]
MM.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% SCL$CCLE.name])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]
Name <- CellLines$CCLE.name[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]

WISP1min <- min(log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))]))) 
WISP1max <- max(log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))]))) 
WISP1shade <- ceiling(255*(1- (log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])) - WISP1min) /(WISP1max - WISP1min))) 
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH1", c(3:ncol(MM.RNAseq.dat))])  
xval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH2"  , c(3:ncol(MM.RNAseq.dat))])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-7,9), ylim = c(-5,9), main = "CDH2 (x) vs CDH1 (y)", col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
points(log2(xval + 0.0001), log2(yval + 0.0001), lwd = 2, pch = 21, col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)], bg = rainbowColor[1+WISP1shade])
lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
#points(rep(10,length(xval)), 5+log2(0.002 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])), lwd = 2, pch = 21, col = c("black"), bg = rainbowColor[1+WISP1shade]) #rainbowColor[25.6*c(0:10)]

WISP1val <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])
points(log2(pCDH2(WISP1val, max = 64, kd = 0.1)), log2(pCDH1(WISP1val, kd = 0.01)), lwd = 1, pch = 24, col = "black", bg =rainbowColor[1+WISP1shade])

Cellxy <- cbind(log2(xval + 0.0001), log2(yval + 0.0001))
Modxy <- cbind(log2(pCDH2(WISP1val, max = 64, kd = 0.1)), log2(pCDH1(WISP1val, kd = 0.01)))
Distxy <- diag(distmat(Cellxy, Modxy))
xval.r <- xval[Distxy >= 5.5]
yval.r <- yval[Distxy >= 5.5]
Name.r <- Name[Distxy >= 5.5]
BrCa.Name.keep <- Name[Distxy < 5.5]
text(log2(xval.r + 0.0001), log2(yval.r + 0.0001)+0.2, Name.r, cex = 0.5, col = "black")
dev.off() # End plot 

#Use filtered breast cancer cell lines
MM.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% BrCa.Name.keep])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]
EMTsig <- c("WISP1", "SNAI1", "SNAI2", "ZEB2", "FN1", "VIM", "CDH2", "CDH1", "ZEB1")

pdf("EMT-Genes-BrCa-Filtered.pdf", width = 7, height = 7)
EMTsigUP <- c("WISP1", "AHNAK", "BMP1", "CALD1", "CAMK2N1", "CDH2", "CDH6", "COL1A2", "COL3A1", "COL5A2", "FN1", "FOXC2", "GNG11", "GSC", 
              "IGFBP4", "ITGA5", "ITGAV", "MMP2", "MMP3", "MMP9", "MSN", "SERPINE1", "SNAI1", "SNAI2", "SNAI3", "SOX10", 
              "SPARC", "STEAP1", "TCF4", "TIMP1", "TMEFF1", "TMEM132A", "TWIST1", "VCAN", "VIM", "VPS13A", "WNT5A", "WNT5B")

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigUP, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigUP))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigUP[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT UP:", EMTsigUP[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

EMTsigDN <- c("WISP1", "CAV2", "CDH1", "CDH3", "DSP", "FGFBP1", "IL1RN", "KRT19", "MITF", "MST1R", "NUDT13", "OCLN", "RGS2", "SPP1", 
              "TFPI2", "TSPAN13") # "PPPDE2", 

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigDN, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigDN))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigDN[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT DN:", EMTsigDN[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

EMTsigTF <- c("WISP1", "CTNNB1", "ESR1", "FOXC2", "GSC", "MITF", "MSX1", "NOTCH1", "SMAD2", "SNAI1", "SNAI2", "SNAI3", "SOX10", 
              "STAT3", "TCF3", "TCF4", "TWIST1", "ZEB1", "ZEB2") # "SIP1", 

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigTF, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigTF))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigTF[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT TF:", EMTsigTF[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}
#source("filename")
opar <- par(mfrow = c(1,1))
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH2", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNNB1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNND1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNNA1", c(3:ncol(MM.RNAseq.dat))]) 
xval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1"  , c(3:ncol(MM.RNAseq.dat))])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), ylim = c(5,10), main = "CDH1 + CDH2 + CTNNB1 + CTNND1 + CTNNA1", col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)

text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))

dev.off() # End plot 

#Just plot for skin cancer cell lines
pdf("CDH1vsCDH2-Mela.pdf", width = 7, height = 7)
#Just plot for skin cancer cell lines
SCL <- CellLines[CellLines$Site.Primary %in% c("skin"),]
MM.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% SCL$CCLE.name])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]
Name <- CellLines$CCLE.name[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]

WISP1min <- min(log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))]))) 
WISP1max <- max(log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))]))) 
WISP1shade <- ceiling(255*(1- (log2(0.0001 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])) - WISP1min) /(WISP1max - WISP1min))) 
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH1", c(3:ncol(MM.RNAseq.dat))])  
xval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH2"  , c(3:ncol(MM.RNAseq.dat))])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-7,9), ylim = c(-5,9), main = "CDH2 (x) vs CDH1 (y)", col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
points(log2(xval + 0.0001), log2(yval + 0.0001), lwd = 2, pch = 21, col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)], bg = rainbowColor[1+WISP1shade])
lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
#points(rep(10,length(xval)), 5+log2(0.002 + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])), lwd = 2, pch = 21, col = c("black"), bg = rainbowColor[1+WISP1shade]) #rainbowColor[25.6*c(0:10)]

WISP1val <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1", c(3:ncol(MM.RNAseq.dat))])
points(log2(pCDH2(WISP1val, max = 25, kd = 0.01)), log2(pCDH1(WISP1val, kd = 0.003)), lwd = 1, pch = 24, col = "black", bg =rainbowColor[1+WISP1shade])

Cellxy <- cbind(log2(xval + 0.0001), log2(yval + 0.0001))
Modxy <- cbind(log2(pCDH2(WISP1val, max = 64, kd = 0.1)), log2(pCDH1(WISP1val, kd = 0.01)))
Distxy <- diag(distmat(Cellxy, Modxy))
xval.r <- xval[Distxy >= 11]
yval.r <- yval[Distxy >= 11]
Name.r <- Name[Distxy >= 11]
Melan.Name.keep <- Name[Distxy < 11]
text(log2(xval.r + 0.0001), log2(yval.r + 0.0001)+0.2, Name.r, cex = 0.5, col = "black")
dev.off() # End plot 

plot(density(Distxy, adj = 0.1))

#Use filtered skin cancer cell lines
MM.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% Melan.Name.keep])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(MM.RNAseq.dat[,c(3:ncol(MM.RNAseq.dat))])]
EMTsig <- c("WISP1", "SNAI1", "SNAI2", "ZEB2", "FN1", "VIM", "CDH2", "CDH1", "ZEB1")

pdf("EMT-Genes-Melan-Filtered.pdf", width = 7, height = 7)
EMTsigUP <- c("WISP1", "AHNAK", "BMP1", "CALD1", "CAMK2N1", "CDH2", "CDH6", "COL1A2", "COL3A1", "COL5A2", "FN1", "FOXC2", "GNG11", "GSC", 
              "IGFBP4", "ITGA5", "ITGAV", "MMP2", "MMP3", "MMP9", "MSN", "SERPINE1", "SNAI1", "SNAI2", "SNAI3", "SOX10", 
              "SPARC", "STEAP1", "TCF4", "TIMP1", "TMEFF1", "TMEM132A", "TWIST1", "VCAN", "VIM", "VPS13A", "WNT5A", "WNT5B")

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigUP, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigUP))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigUP[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT UP:", EMTsigUP[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

EMTsigDN <- c("WISP1", "CAV2", "CDH1", "CDH3", "DSP", "FGFBP1", "IL1RN", "KRT19", "MITF", "MST1R", "NUDT13", "OCLN", "RGS2", "SPP1", 
              "TFPI2", "TSPAN13") # "PPPDE2", 

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigDN, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigDN))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigDN[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT DN:", EMTsigDN[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

EMTsigTF <- c("WISP1", "CTNNB1", "ESR1", "FOXC2", "GSC", "MITF", "MSX1", "NOTCH1", "SMAD2", "SNAI1", "SNAI2", "SNAI3", "SOX10", 
              "STAT3", "TCF3", "TCF4", "TWIST1", "ZEB1", "ZEB2") # "SIP1", 

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigTF, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigTF))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigTF[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT TF:", EMTsigTF[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}

EMTsigMela <- c("WISP1", "MITF", "TYR", "PMEL", "MLANA", "CDH1", "IRF4", "AXL", "NGFR", "SERPINE1", "FN1", "S100B", "CNTN6", "L1CAM", "FYN", "MAP2", "NCAM1")  

MM.RNAseq.EMT <- MM.RNAseq.dat[MM.RNAseq.dat$Description %in% EMTsigMela, ]

opar <- par(mfrow = c(3,3))
for (i in 2:length(EMTsigMela))
{
  xval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == "WISP1"  , c(3:ncol(MM.RNAseq.EMT))])
  yval <- as.numeric(MM.RNAseq.EMT[MM.RNAseq.EMT$Description == EMTsigMela[i], c(3:ncol(MM.RNAseq.EMT))])
  plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = paste("EMT Melanoma:", EMTsigMela[i]), sub = paste("Cor=", sprintf("%6.3f", cor.test(xval, yval)$estimate),"; P-val=", sprintf("%6.4f", cor.test(xval, yval)$p.value)), col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
  lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
  lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)
  #text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
  #text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
}
#source("filename")
opar <- par(mfrow = c(1,1))
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CDH2", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNNB1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNND1", c(3:ncol(MM.RNAseq.dat))]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "CTNNA1", c(3:ncol(MM.RNAseq.dat))]) 
xval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1"  , c(3:ncol(MM.RNAseq.dat))])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), ylim = c(5,10), main = "CDH1 + CDH2 + CTNNB1 + CTNND1 + CTNNA1", col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)

text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))

dev.off() # End plot 


opar <- par(mfrow = c(1,1))
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "TMED10", c(3:119)]) 
yval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "PSEN1", c(3:119)]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "APH1B", c(3:119)]) *  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "APH1B", c(3:119)]) /(as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "APH1A", c(3:119)]) + as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "APH1B", c(3:119)])) +
  #   as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "PDLIM1", c(3:119)]) +
  as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "NCSTN", c(3:119)]) 
xval <- as.numeric(MM.RNAseq.dat[MM.RNAseq.dat$Description == "WISP1"  , c(3:119)])
plot(log2(xval + 0.0001), log2(yval + 0.0001), xlim = c(-13,6), main = "CDH1 + CTNNB1 + CTNND1 + CTNNA1", col = c("red", "black", "black", "black", "black", "blue")[as.factor(Site)])
lines(log2(c(0.2, 0.2)), c(-5,14), col = "blue", lty = 2)
lines(c(-13,6), log2(c(0.2, 0.2)), col = "blue", lty = 2)

text(2.5, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
text(2.5, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))

text(7, log2(min(yval + 0.0001)) + 0.95*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
text(7, log2(min(yval + 0.0001)) + 0.9*(log2(max(yval)) - log2(min(yval + 0.0001))), paste("P-value = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
dev.off() # End plot 


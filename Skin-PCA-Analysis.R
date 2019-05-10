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

#MM.CNV.dat <- cbind(CNV.dat[,c(1:5)], CNV.dat[,colnames(CNV.dat) %in% MMCL$CCLE.name])
BC.RNAseq.dat <- cbind(RNAseq.dat[,c(1:2)], RNAseq.dat[,colnames(RNAseq.dat) %in% SCL$CCLE.name])
Site <- CellLines$Site.Primary[CellLines$CCLE.name %in% colnames(BC.RNAseq.dat[,c(3:ncol(BC.RNAseq.dat))])]
rm(RNAseq.dat)

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
rm(BC.RNAseq.dat)
gc()

#Convert from RPKM to log2(RPKM + 0.0001)
SeqRes <- log2(BC.RNAseq.EMT[,c(3:ncol(BC.RNAseq.EMT))] + 0.0001)

# # Principal component analysis
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

opar <- par(mfrow = c(2,2))
pdf("PCA-analysis-Skin-PCA-Scree-Nov18.pdf", width = 5, height = 5)
mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,80), cex.names = 1.5, cex.axis = 1.5, ylab = "Variance (%)", xlab = "Principal Component")
lines(c(0,2,3,4,5,6,7,8,9,12), c(3.4, 3.2, 3.19, 3.04, 3.01, 2.89, 2.8, 2.74, 2.67, 2.63), col="red", lwd = 2, lty = 2)
dev.off() # End plot 

# Keep gene for further analysis if it is expressed above the threshold in more than 
# 5% of the samples
GeneExpressed <- apply(SeqRes, 1, function(x) length(x[x > -0.304])/length(x) > 0.05 )
rainColor <- colorRampPalette(c("red", "yellow", "blue"))(100)
meanSeqRes <- rowMeans(SeqRes)
PtCol <- rainColor[cut(meanSeqRes,breaks = 100)]
#PtCol <- ifelse(meanSeqRes > -0.304, "black", "red")
#PtCol <- ifelse(GeneExpressed, "black", "red")

pdf("PCA-analysis-Skin-PC1vsPC2-Nov18.pdf", width = 7, height = 7)
opar <- par(mfrow = c(1,1))
plot(-1*PCA2.data[,1], PCA2.data[,2], type = "n")#, xlim = c(-100, 70), ylim = c(-35, 40))
#points(PCA2.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
text(-1*PCA2.data[,1], PCA2.data[,2], labels = BC.RNAseq.EMT$Description, col = PtCol)
points(-1*PCA2.data[BC.RNAseq.EMT$Description == "CDH1",1], PCA2.data[BC.RNAseq.EMT$Description == "CDH1",2], pch = 19, col = "blue", cex = 1)
points(-1*PCA2.data[BC.RNAseq.EMT$Description == "VIM",1], PCA2.data[BC.RNAseq.EMT$Description == "VIM",2], pch = 19, col = "red", cex = 1)
#lines(c(-100, 70), c(19, 19), col="red", lwd=2, lty=2)
#lines(c(-100, 70), c(-19, -19), col="red", lwd=2, lty=2)
lines(c(50, 50), c(-100, 70), col="blue", lwd=2, lty=2)
dev.off() # End plot 

pdf("PCA-analysis-Skin-PC2vsPC3-Nov18.pdf", width = 7, height = 7)
KeepGene <- GeneExpressed
# KeepGene <- PCA2.data[,1] > -30 
plot(PCA2.data[KeepGene,c(3,2)], type = "n", xlim = c(-22, 22), ylim = c(-34, 26))
#points(PCA2.data[,c(1,2)], pch = 19, col = PtCol, cex = 1)
text(PCA2.data[KeepGene,c(3,2)], labels = BC.RNAseq.EMT$Description[KeepGene], cex = 1.0, col = "black" ) #PtCol[KeepGene]
points(PCA2.data[BC.RNAseq.EMT$Description == "CDH1",3], PCA2.data[BC.RNAseq.EMT$Description == "CDH1",2], pch = 19, col = "blue", cex = 1)
points(PCA2.data[BC.RNAseq.EMT$Description == "VIM",3], PCA2.data[BC.RNAseq.EMT$Description == "VIM",2], pch = 19, col = "red", cex = 1)
lines(c(-100, 70), c(11.9, 11.9), col="red", lwd=2, lty=2)
lines(c(-100, 70), c(-11.9, -11.9), col="red", lwd=2, lty=2)
lines(c(-11.9, -11.9), c(-100, 70), col="blue", lwd=2, lty=2)
lines(c(11.9, 11.9), c(-100, 70), col="blue", lwd=2, lty=2)
dev.off() # End plot 

outPCA <- data.frame(name = BC.RNAseq.EMT$Description, GeneMean = meanSeqRes, KeepGene = GeneExpressed, PC1 = PCA2.data[,1], PC2 = PCA2.data[,2], PC3 = PCA2.data[,3])
write.csv(outPCA, file = "Skin-GenePCAresult-Nov18.csv")
plot(density(PCA2.data[,3], adj = 0.25))
L1 <- density(c(NCx, NCy), adjust = 0.25, na.rm=TRUE, from = -30, to = 30)
lines(L1, col = "red")

#######
# Plot relationship between transcript length and GC content
TLGC.file.name <- "./data/hsa-GeneLengthAndGC.csv"
TLGC.dat <- read.table(TLGC.file.name, head=TRUE, sep = ",", skip = 0, na.strings = "null")

TxLi <- as.numeric(TLGC.dat$transcript_length[match(as.character(BC.RNAseq.EMT$Name), as.character(TLGC.dat$ensembl_gene_id))])
GCi <- as.numeric(TLGC.dat$percentage_gene_gc_content[match(as.character(BC.RNAseq.EMT$Name), as.character(TLGC.dat$ensembl_gene_id))])

opar <- par(mfrow = c(3,1))
plot(cbind(PCA2.data[,1], TxLi), type = "n")#, xlim = c(-100, 70), ylim = c(-35, 40))
text(cbind(PCA2.data[,1], TxLi), labels = BC.RNAseq.EMT$Description, col = PtCol)

plot(cbind(PCA2.data[,1], GCi), type = "n")#, xlim = c(-100, 70), ylim = c(-35, 40))
text(cbind(PCA2.data[,1], GCi), labels = BC.RNAseq.EMT$Description, col = PtCol)

plot(cbind(PCA2.data[,1], meanSeqRes), type = "n")#, xlim = c(-100, 70), ylim = c(-35, 40))
text(cbind(PCA2.data[,1], meanSeqRes), labels = BC.RNAseq.EMT$Description, col = PtCol)

#######

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
#Control PCA loadings
#Summary results
Ntot = 50
NCx = vector('double',783*Ntot)
NCy = vector('double',783*Ntot)

for ( i in 1:Ntot )
{
  #Negative control 
  RS <- sample(simplify2array(SeqRes), 783*56, replace=T)
  A2 <- matrix(c(RS), nrow = 783)
  Rdat <- SeqRes
  Rdat[c(1:783),c(1:56)] <- A2

  NC_PCA <- prcomp(Rdat, retx = TRUE, scale. = FALSE)
  NC_PCA.data <- as.matrix(Rdat) %*% NC_PCA$rotation
  NCx[783*(i-1)+c(1:783)] <- NC_PCA.data[,1]
  NCy[783*(i-1)+c(1:783)] <- NC_PCA.data[,2]
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


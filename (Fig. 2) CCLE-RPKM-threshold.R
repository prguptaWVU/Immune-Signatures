library(colorspace)
#library(robustbase)
#library(MASS)
setwd("~/Documents/Projects/CCLE/R")
rm(list = ls())

CellLines <- read.table("./data/CCLE_sample_info_file_2012-10-18.txt", head=TRUE, sep = "\t", colClasses = c("character"))

# Read in RNAseq data from CCLE
RNAseq.file.name <- "./data/CCLE_RNAseq_081117.rpkm.gct"
RNAseq.dat <- read.table(RNAseq.file.name, head=TRUE, sep = "\t", skip = 2, na.strings = "null")
GenesIDs <- as.factor(sapply(RNAseq.dat$Name, function(x) unlist(strsplit(as.character(x), "[.]"))[1]))
RNAseq.dat$Name <- GenesIDs

# Read in RPPA data from MD Anderson
RPPA.file.name <- "./data/MCLP-v1.1-Level4.txt"
RPPA.dat <- read.table(RPPA.file.name, head=TRUE, sep = "\t", skip = 1, na.strings = "NA")
RPPA.Genes <- colnames(RPPA.dat)[-c(1,2)]

# Filter cell lines to just include skin. 57 cell lines and breast. 57 cell lines
#SCL <- CellLines[(CellLines$Histology %in% c("malignant_melanoma")) | (CellLines$Site.Primary %in% c("breast")),]
#SCL <- CellLines[(CellLines$Histology %in% c("malignant_melanoma")),]
#CommonCL <- SCL$CCLE.name[(SCL$CCLE.name %in% RPPA.dat$Sample_Name) & (SCL$CCLE.name %in% colnames(RNAseq.dat))]
#CommonCL <- SCL$CCLE.name[!(SCL$CCLE.name %in% colnames(RNAseq.dat))]

CommonCL <- CellLines$CCLE.name[(CellLines$CCLE.name %in% RPPA.dat$Sample_Name) & (CellLines$CCLE.name %in% colnames(RNAseq.dat))]

# Genes that have a good relationship between RPPA and RNAseq measures
CorrGenes <- c("LCK", "RAB25", "PTGS2", "BCL2L11", "CASP8", "FN1", "IGFBP2", "IRS1", 
               "GATA2", "MME", "KIT", "INPP4B", "SMAD1", "ITGA2", "PDGFRB", "CDKN1A", "MACC1", 
               "EGFR", "SERPINE1", "CAV2", "CAV1", "DUSP4", "NRG1", "MYC", "EPPK1", "CD274",
               "CDKN2A", "ANXA1", "SYK", "GATA3", "MGMT", "CD44", "ETS1", "ERBB3", "RB1", 
               "CDH1", "CLDN7", "TP53", "STAT5A", "PRKCA", "CDH2", "BCL2", "NOTCH3", "AXL",
               "JUN", "PARP1", "MSH6", "REL", "GLS", "CTNNA1", "SQSTM1", "TP63", "ASNS", 
               "MET", "EIF4EBP1", "PTK2", "NOTCH1", "LDHA", "GAB2", "KRT8", "PXN", "TYRO3", 
               "ERBB2", "FASN", "KEAP1", "CCNE1", "JAG1", "PYGB", "SRC", "TGM2", "PREX1", 
               "SYP", "AR", "G6PD", "PARK7", "MTOR", "CASP9", )


x_Threshold <- rep(0, length(CorrGenes))
var_Threshold <- rep(0, length(CorrGenes))
pdf("RNAseq-vs-RPPA-CorrGenes-rlm2.pdf", width = 7, height = 7)

dotcol <- sapply(CommonCL, function(x) ifelse(grepl("SKIN", x), "red", ifelse(grepl("BREAST", x), "blue", "black")))
opar <- par(mfrow = c(3,3))
for (i in 1:length(CorrGenes))
{
  xval <- as.numeric(RNAseq.dat[RNAseq.dat$Description == CorrGenes[i],CommonCL])
  yval <- as.numeric(RPPA.dat[match(CommonCL, RPPA.dat$Sample_Name), CorrGenes[i]])

  CalcStats <- length(yval[!is.na(yval)]) > 2
  if(CalcStats)
  {
    #  lmres <- lm(yval ~ xval)
  #  lmres <- lmrob(yval ~ xval)
    lmres <- glm(yval ~ xval, family = binomial())
  #  lmres <- rlm(yval ~ xval) # MASS package
    summary(lmres)
    plot(log10(xval+0.01), yval, col = dotcol, xlim = c(-2.5,4), ylim = c(-5,3), main = CorrGenes[i], xlab = "RNAseq", ylab = "RPPA")
    xsam <- seq(-2.5,3,0.1)
    lines(xsam, lmres$coefficients[1] + lmres$coefficients[2]*(10^xsam), lty = 2)
    x_Threshold[i] <- var(lmres$residuals)/lmres$coefficients[2]
    lines(log10(c(x_Threshold[i], x_Threshold[i])), c(-5,3), lty = 2, col = "red")
    text(2, -4.3, paste("Thresh = ", sprintf("%6.1f", x_Threshold[i])))
    text(-1, 2.5, paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
    text(-1, 1.8, paste("P-val = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
  }
  plot(lmres)
}

plot(density(log10(x_Threshold), adj = 0.5), xlab = "Threshold Value", ylab = "Density", main = list(cex = 0.8, paste("Median Threshold =", sprintf("%6.1f", median(x_Threshold)))))

dev.off() # End plot 

# Test gene
Tgene <- "STMN1"

xval <- as.numeric(RNAseq.dat[match(Tgene, RNAseq.dat$Description),CommonCL])
yval <- as.numeric(RPPA.dat[match(CommonCL, RPPA.dat$Sample_Name), Tgene])
  
lmres <- nls(yval ~ int + (10^Max)*xval/(xval + KD), start = list(int = -1.2, KD = 11, Max = 0.5), control = list(warnOnly = TRUE), trace = TRUE)
    
plot(log10(xval+0.01), yval, col = dotcol, xlim = c(-2.5,4), ylim = c(-5,3), main = CommonGenes[i], xlab = "RNAseq", ylab = "RPPA")
xsam <- seq(-2.5,4,0.1)
lines(xsam, coef(lmres)[1] + (10^coef(lmres)[3])*10^xsam/(10^xsam + coef(lmres)[2]), lty = 2)

x_Thresh <- var(yval - (coef(lmres)[1] + (10^coef(lmres)[3])*xval/(xval + coef(lmres)[2])), na.rm = TRUE)*coef(lmres)[2]/(10^coef(lmres)[3])
lines(log10(c(x_Thresh, x_Thresh)), c(-5,3), lty = 2, col = "red")
text(2, -4.3, paste("Thresh = ", sprintf("%6.1f", x_Thresh)))

text(-1, 2.5, paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
text(-1, 1.8, paste("P-val = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))

# End Test gene 

pdf("RNAseq-vs-RPPA-CorrGenes-nls-0p1b.pdf", width = 7, height = 7)
x_Threshold <- rep(NA, length(CorrGenes))
var_Threshold <- rep(NA, length(CorrGenes))
dotcol <- sapply(CommonCL, function(x) ifelse(grepl("SKIN", x), "red", ifelse(grepl("BREAST", x), "blue", "black")))
opar <- par(mfrow = c(3,3))
for (i in 1:length(CorrGenes))
{
  xval <- as.numeric(RNAseq.dat[match(CorrGenes[i], RNAseq.dat$Description),CommonCL])
  yval <- as.numeric(RPPA.dat[match(CommonCL, RPPA.dat$Sample_Name), CorrGenes[i]])
  
  #  lmres <- lm(yval ~ xval)
  lmres <- nls(yval ~ int + Max*xval/(xval + KD), start = list(int = 0, KD = 11, Max = 3), 
              control = list(warnOnly = TRUE))

  summary(lmres)
  plot(log10(xval+0.01), yval, col = dotcol, xlim = c(-2.5,4), ylim = c(-5,3), main = CorrGenes[i], xlab = "RNAseq", ylab = "RPPA")
  xsam <- seq(-2.5,4,0.1)
  lines(xsam, coef(lmres)[1] + coef(lmres)[3]*10^xsam/(10^xsam + coef(lmres)[2]), lty = 2)
  if(coef(lmres)[3] > 2){
#    x_Threshold[i] <- var(yval - (coef(lmres)[1] + coef(lmres)[3]*xval/(xval + coef(lmres)[2])), na.rm = TRUE)*coef(lmres)[2]/coef(lmres)[3]
    x_Threshold[i] <- 0.1*coef(lmres)[2]/(coef(lmres)[3]-0.1) #0.291 is the median residual variance
    var_Threshold[i] <- var(yval - (coef(lmres)[1] + coef(lmres)[3]*xval/(xval + coef(lmres)[2])), na.rm = TRUE)
    lines(log10(c(x_Threshold[i], x_Threshold[i])), c(-5,3), lty = 2, col = "red")
    text(2, -4.3, paste("Thresh = ", sprintf("%6.1f", x_Threshold[i])))
  }
  CalcStats <- length(yval[!is.na(yval)]) > 2
#  if(CalcStats)
#  {
#    text(-1, 2.5, paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
#    text(-1, 1.8, paste("P-val = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
#  }
  plot(log10(xval+0.01),yval - (coef(lmres)[1] + coef(lmres)[3]*xval/(xval + coef(lmres)[2])), xlim = c(-2.5,4), ylim = c(-3,3))
}

plot(density(log10(x_Threshold), adj = 0.5, na.rm = TRUE), xlab = "Threshold Value", ylab = "Density", main = list(cex = 0.8, paste("Median Threshold =", sprintf("%6.1f", median(x_Threshold, na.rm = TRUE)))))

dev.off() # End plot 

#################
#
# Do analysis for all common genes
#

CommonGenes <- as.character(RNAseq.dat$Description[as.character(RNAseq.dat$Description) %in% RPPA.Genes])
BadActors <- c("KRT5", "GAPDH", "BRCA1", "SCD", "PRDX1", "KIAA1324")

pdf("RNAseq-vs-RPPA-CommonGenes-nls-CorrGT0p36.pdf", width = 7, height = 7)
x_Threshold <- rep(NA, length(CommonGenes))
dotcol <- sapply(CommonCL, function(x) ifelse(grepl("SKIN", x), "red", ifelse(grepl("BREAST", x), "blue", "black")))
opar <- par(mfrow = c(3,3))
for (i in 1:length(CommonGenes))
{
  xval <- as.numeric(RNAseq.dat[match(CommonGenes[i], RNAseq.dat$Description),CommonCL])
  yval <- as.numeric(RPPA.dat[match(CommonCL, RPPA.dat$Sample_Name), CommonGenes[i]])
  
  CalcStats <- (length(yval[!is.na(yval)]) > 15 & length(xval[xval > 0 & !is.na(yval)]) > 15)
  if(CalcStats)
  {
    #  lmres <- lm(yval ~ xval)
    lmres <- nls(yval ~ int + (10^Max)*xval/(xval + KD), start = list(int = -1.2, KD = 11, Max = 0.5), control = nls.control(warnOnly = TRUE), trace = TRUE)
  
    #summary(lmres)
    plot(log10(xval+0.01), yval, col = dotcol, xlim = c(-2.5,4), ylim = c(-5,3), main = CommonGenes[i], xlab = "RNAseq", ylab = "RPPA")
    xsam <- seq(-2.5,4,0.1)
    lines(xsam, coef(lmres)[1] + (10^coef(lmres)[3])*10^xsam/(10^xsam + coef(lmres)[2]), lty = 2)
    if(cor.test(xval, yval)$estimate > 0.36 & !(CommonGenes[i] %in% BadActors)){
      x_Threshold[i] <- 0.025*coef(lmres)[2]/0.975 #10% above baseline
#      var_Threshold[i] <- var(yval - (coef(lmres)[1] + (10^coef(lmres)[3])*xval/(xval + coef(lmres)[2])), na.rm = TRUE)*coef(lmres)[2]/(10^coef(lmres)[3])
      lines(log10(c(x_Threshold[i], x_Threshold[i])), c(-5,3), lty = 2, col = "red")
      text(2, -4.3, paste("Thresh = ", sprintf("%6.1f", x_Threshold[i])))
    }
    text(-1, 2.5, paste("Corr = ", sprintf("%6.4f", cor.test(xval, yval)$estimate)))
    text(-1, 1.8, paste("P-val = ", sprintf("%6.4f", cor.test(xval, yval)$p.value)))
#    plot(log10(xval+0.01),yval - (coef(lmres)[1] + coef(lmres)[3]*xval/(xval + coef(lmres)[2])), xlim = c(-2.5,4), ylim = c(-3,3))
  }
}

plot(density(log10(x_Threshold), adj = 0.5, from = -3, to = 4, na.rm = TRUE), xlab = "Threshold Value", xlim = c(-2,2), ylab = "Density", lwd = 2, main = list(cex = 0.8, paste("Median Threshold =", sprintf("%6.1f", median(x_Threshold, na.rm = TRUE)))))

dev.off() # End plot 

##################################################################
#
# EMT signature correlates
#
##################################################################
# Filter cell lines to just include skin. 57 cell lines and breast. 57 cell lines
SCL <- CellLines[CellLines$Site.Primary %in% c("breast"),]
SCL <- CellLines[(CellLines$Histology %in% c("malignant_melanoma")) | (CellLines$Site.Primary %in% c("breast")),]
#SCL <- CellLines[(CellLines$Histology %in% c("malignant_melanoma")),]

summary(as.factor(SCL$Hist.Subtype1))
summary(as.factor(SCL$Site.Primary))

CommonCL <- SCL$CCLE.name[(SCL$CCLE.name %in% RPPA.dat$Sample_Name) & (SCL$CCLE.name %in% colnames(RNAseq.dat))]
dotcol <- sapply(CommonCL, function(x) ifelse(grepl("SKIN", x), "red", ifelse(grepl("BREAST", x), "blue", "black")))



RNAseq.dat[RNAseq.dat$Description == "CDH1", CommonCL]
RPPA.dat[as.factor(CommonCL), c("CDH1")]

RPPA.dat[, c("CDH1")]
CellLines$CCLE.name[CellLines$Site.Primary %in% c("kidney")]

# Filter cell lines to just include breast. 57 cell lines
SCL <- CellLines[CellLines$Site.Primary %in% c("breast"),]
summary(as.factor(SCL$Hist.Subtype1))
summary(as.factor(SCL$Histology))

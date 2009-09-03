# preprocessing Affy arrays on UNIX server
#
# TODO:
# - change pathways to files in a more generally valid way
# - make last part also directly applicable for gcrma as well as farms?
# - here, a rather exceptionl case is used, namely the 133B array that is not implemented in the pipeline
# a good idea?
#
# Author: Willem Talloen, ...
###############################################################################


#==================================================================
# reading in data
library(affy)
library(farms)
library(gcrma)

# setwd('/home/wtalloen/projects/2009/oncology/Hdm2/signatureWork/p53activity/wd')
setwd('/projects/bix/Public_Non_Incyte/Microarray/Projects/2009/Oncology/Hdm2/signatureWork/caseStudies/p53activity/wd')
pDat <- read.delim("../data/pDataMiller2005A_B.txt", header = TRUE)
pDat <- pDat[!is.na(pDat$ GEO.Sample.Accession..),]
rownames(pDat) <- paste(pDat$ GEO.Sample.Accession.., '.CEL', sep='')

selCel <- rownames(pDat)[which(pDat$Affy.platform == 'HG-U133B')]

setwd("/dvl/bix/current/chip/docs/external/GSE3494/cel_files")
listf <- list.files()

CELs2 <- intersect(listf,selCel)

pDat2 <- pDat[CELs2,]
data.raw <- ReadAffy(filenames = CELs2, phenoData = pDat2[CELs2,])
setwd('/projects/bix/Public_Non_Incyte/Microarray/Projects/2009/Oncology/Hdm2/signatureWork/caseStudies/p53activity/wd')
#==================================================================


#==================================================================
# preprocessing

load('expressionSetGcrma18days.Rda')
#-----------------------------------------
# using Entrez based Probe Set definitions
.libPaths("~/RlibraryUnix")
library(hgu133bhsentrezgprobe)
library(hgu133bhsentrezgcdf)
data.raw@cdfName <- "hgu133bhsentrezg"
#-----------------------------------------

#-----------------------------------------
# GCRMA
GcrmaEntrez <- gcrma(data.raw)
pData(GcrmaEntrez) <- pDat2
#-----------------------------------------

#-----------------------------------------
# FARMS and I/NI calls
FarmsEntrez <- q.farms(data.raw)
INIs <- INIcalls(FarmsEntrez)
summary(INIs)
iFarmsEntrez <- getI_Eset(INIs)
pData(iFarmsEntrez) <- pDat2
save(iFarmsEntrez, file = "iFarmsEntrez.rda")
#-----------------------------------------

#-----------------------------------------
# add featureData
library(a4)
library(hgu133plus2hsentrezgJnJ.db)

GcrmaEntrez <- GcrmaEntrez[-grep('AFFX-',featureNames(GcrmaEntrez)),]
GcrmaEntrez <- addGeneInfo_hgu133plus2hsentrezgJnJ(GcrmaEntrez)
#-----------------------------------------

#-----------------------------------------
# Remove the .CEL from the sampleNames 
sampleNames(GcrmaEntrez) <- sub('.CEL','',sampleNames(GcrmaEntrez))
colnames(pData(GcrmaEntrez)) <- sub('.CEL','',colnames(pData(GcrmaEntrez)))
#-----------------------------------------

save(GcrmaEntrez, file='eset.rda')
#==================================================================

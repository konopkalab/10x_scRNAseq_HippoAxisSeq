rm(list = ls())
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(reshape2)


#load data
P57.data <- Read10X(data.dir = "P57/outs/filtered_feature_bc_matrix")
P57 <- CreateSeuratObject(counts = P57.data, project ="P57")

#filter
P57 <- PercentageFeatureSet(object = P57, pattern = "^MT-", col.name = "percent.mt")

P57Filt<-subset(x = P57, subset = nCount_RNA< 10000 & nFeature_RNA >300& percent.mt<5)



#remove X, Y, M genes
mGenes <- scan("/work/psychiatry/s175390/resources/database/FILTER_MITO_X_Y/cellranger_hg19_genes_ChrM.txt", what = "", sep = "\n") 
xGenes <- scan("/work/psychiatry/s175390/resources/database/FILTER_MITO_X_Y/cellranger_hg19_genes_ChrX.txt", what = "", sep = "\n") 
yGenes <- scan("/work/psychiatry/s175390/resources/database/FILTER_MITO_X_Y/cellranger_hg19_genes_ChrY.txt", what = "", sep = "\n")


mxyGenes <- unique(sort(c(mGenes,xGenes,yGenes))) #1187 genes
keepGenes <- unique(sort(setdiff(row.names(P57Filt), mxyGenes)))
P57Filt_data<-GetAssayData(object = P57Filt)
tempP571 <- as.matrix(P57Filt_data)
keepCells <- colnames(tempP571)
P57Filt_rawdata <- GetAssayData(object = P57Filt, slot = "counts" )
tempP572 <- as.matrix(P57Filt_rawdata)
P57.newdata <- tempP572[keepGenes,keepCells]
P57New <- CreateSeuratObject(counts = P57.newdata, project = "P57_woMTXY")
metaAll <- as.data.frame(P57@meta.data)
P57pMito <- metaAll[keepCells, "percent.mt"]
names(P57pMito) <- row.names(metaAll[keepCells,])
P57New$percent.mt <- P57pMito
save(P57New, P57.newdata, file = "P57_FILT_woMTXY.RData")


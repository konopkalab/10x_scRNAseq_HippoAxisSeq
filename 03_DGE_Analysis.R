rm(list = ls())
library(tidyverse)
library(lme4)
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(MAST)
library(data.table)


#subset clusters, save normalized counts and metadata for DEG analysis

load("processing_qc/UMAP_2000_64_04.RData")
for(test.cluster in unique(hippo.combined@meta.data$seurat_clusters)){
	cluster<-subset(x=hippo.com, idents=test.cluster);
	cluster <- NormalizeData(object = cluster);
	df<-as.data.frame(cluster[["RNA"]]@data);
	metadata<-as.data.frame(cluster@meta.data);
	meta=metadata;
	meta<-meta[,c(2,3,4,5,6,7,8,9,10,11,21)];
	df <- df[,match(rownames(metadata),colnames(df))];
	data1<-as.matrix(df);
	save(data1,meta,file=paste("DGE/S_",test.cluster,".RData",sep=""));
}

#Make SingleCellAssay object to input in MAST
	
group1="anterior"
group2="posterior"

load("DGE/S_1.RData")
test.tsne=meta
test.data=data1[,rownames(test.tsne)]
rows=rowSums(as.matrix(test.data)> 0) > ncol(test.data)*0.10
test.data = test.data[rows,]
print(paste("Testing",nrow(test.data),"genes"))
sca=FromMatrix(as.matrix(test.data),as.data.frame(colnames(test.data)),as.data.frame(rownames(test.data)))
groups=as.character(test.tsne[,4])
axis<-factor(groups)
axis<-relevel(axis,group1)
cdr = colSums(assay(sca)>0)
age=test.tsne$age
sex=as.character(test.tsne$sex)
version=as.character(test.tsne$version)
batch=as.character(test.tsne$batch)
ind=as.character(test.tsne$donor)
dur=test.tsne$dur
mito_perc=test.tsne$percent.mt
colData(sca)$axis<-axis
colData(sca)$cngeneson <- scale(cdr)
colData(sca)$age<-scale(age)
colData(sca)$sex<-sex
colData(sca)$batch<-batch
colData(sca)$version<-version
colData(sca)$ind<-ind
colData(sca)$mito_perc<-scale(mito_perc)
colData(sca)$dur<-scale(dur)
save(sca, file="DGE/SCA_1.RData")

#LMM using MAST to detect DEGs across the axis in a given cluster(cluster 1 here)

form=as.formula("~axis + (1|ind) + cngeneson + age + sex + dur + mito_perc +version")
zlmCond = zlm(form, sca, method = "glmer", ebayes = F, silent=T)
summaryCond <- summary(zlmCond, doLRT='axisposterior')
summaryDt <- summaryCond$datatable
save(summaryDt,file="DGE/LMM_C1.R")
fcHurdle <- merge(summaryDt[contrast=='axisposterior' & component=='H',.(primerid,`Pr(>Chisq)`)],summaryDt[contrast=='axisposterior' & component=='logFC', .(primerid, coef,ci.hi, ci.lo)], by='primerid');
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')];
fcHurdleSig <- merge(fcHurdle[fdr<0.05& abs(coef)>0.3], as.data.table(mcols(sca)),by='primerid');
setorder(fcHurdleSig, fdr);
gene.IDs<-as.character(fcHurdleSig$primerid);
FCs<-as.character(fcHurdleSig$coef);
FDRs<-as.character(fcHurdleSig$fdr);
out=cbind(gene.IDs,FCs,FDRs);
write.table(out,file="DGE/LMM_C1.txt",sep="\t")


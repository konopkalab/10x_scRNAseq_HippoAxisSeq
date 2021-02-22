suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(emmeans)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(magrittr)
library(broom)
library(muscat)
library(Seurat)
library(clustree)
library(leiden)
})
source("input_data/Utils.R")

#plan("multiprocess", workers = 48)
#options(future.globals.maxSize = 5000 * 1024^2)

dir.create("Reanalysis_Scaled")

sce <- readRDS("input_data/HippoAxis.rds")

SeizFreq <- data.frame(orig.ident = c("A56_woXYMT","P56_woXYMT","A57_woXYMT","P57_woXYMT",
                                         "A67_woXYMT","P67_woXYMT","A76_woXYMT","P76_woXYMT",
                                         "A82_woXYMT","P82_woXYMT"),
                          SeizFreq = c(0.5,0.5,3,3,1.5,1.5,1.0,1.0,4.0,4.0))

temporary <- colData(sce) %>%
                as.data.frame() %>%
                rownames_to_column("TMP") %>%
                select(TMP,orig.ident,nCount_RNA,nFeature_RNA, percent.mt, group, age, sex,dur,batch) %>%
                dplyr::rename(Axis = group,nUMI = nCount_RNA,nGene=nFeature_RNA) %>%
                mutate(Axis = as.factor(Axis), sex = as.factor(sex),batch = as.factor(batch),orig.ident = as.factor(orig.ident)) %>%
                left_join(SeizFreq) %>%
                column_to_rownames("TMP") %>% 
                DataFrame()

colData(sce) <- temporary


# Create a Seurat Object
sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

processing_seurat_scaling <- function(object, vars_to_regress = NULL, npcs = 50, res = 0.8) {
  	object %>% 
  	NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%  
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>% 
    ScaleData(vars.to.regress = vars_to_regress, model.use = "linear") %>% 
    RunPCA(features=NULL, weight.by.var = TRUE, npcs = npcs, reduction.name = "pca") %>%
    JackStraw(dims = npcs) %>% 
    ScoreJackStraw(dims = 1:npcs) %>%
    FindNeighbors(reduction = "pca", dims = 1:npcs, nn.eps = 0.5) %>%  
 	  FindClusters(resolution = res, algorithm = 1,n.iter = 1000,save.SNN = TRUE) %>%
    RunTSNE(reduction = "pca", dims = 1:npcs) %>% 
    RunUMAP(reduction = "pca", dims = 1:npcs) 
}


Hippo_Scaled <- processing_seurat_scaling(sce.seurat, 
											vars_to_regress = c("nUMI","percent.mt","age","dur","batch","SeizFreq"),
											npcs = 30, 
											res = 0.8)

save(Hippo_Scaled, file = "Reanalysis_Scaled/Hippo_Data_Scaling.RData")

sessionInfo()

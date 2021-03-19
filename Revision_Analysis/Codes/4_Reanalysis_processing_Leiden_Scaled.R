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

# Load object
load("Reanalysis_Scaled/Hippo_Data_Scaling.RData")

# Calculate Leiden
membership <- leiden(Hippo_Scaled@graphs$RNA_snn, resolution_parameter = 0.8)
Idents(object = Hippo_Scaled) <- as.factor(membership)

# Select resolution and run UMAP
Hippo_Scaled <- RunUMAP(object = Hippo_Scaled, 
										reduction = "pca", 
										dims = 1:30)

pdf("Reanalysis_Scaled/Hippo_Scaled_UMAP_Leiden.pdf", width = 7, height = 6)
DimPlot(object = Hippo_Scaled, reduction = "umap", group.by = "ident", label = TRUE, pt.size = 0.5)
dev.off()

save(Hippo_Scaled,file = "Reanalysis_Scaled/Hippo_Scaled_Leiden.RData")
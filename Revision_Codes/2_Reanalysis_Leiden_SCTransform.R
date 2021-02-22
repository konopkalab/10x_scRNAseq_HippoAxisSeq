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
load("Reanalysis_SCTransform/Hippo_Data_Integrated.RData")

# Calculate Leiden
membership <- leiden(Hippo_Integrated@graphs$integrated_snn, resolution_parameter = 0.6)
Idents(object = Hippo_Integrated) <- as.factor(membership)

# Select resolution and run UMAP
Hippo_Integrated <- RunUMAP(object = Hippo_Integrated, 
										reduction = "pca", 
										dims = 1:30)

pdf("Reanalysis_SCTransform/Hippo_Integrated_UMAP_Leiden.pdf", width = 7, height = 6)
DimPlot(object = Hippo_Integrated, reduction = "umap", group.by = "ident", label = TRUE, pt.size = 0.5)
dev.off()

save(Hippo_Integrated,file = "Reanalysis_SCTransform/Hippo_Integrated_Leiden.RData")
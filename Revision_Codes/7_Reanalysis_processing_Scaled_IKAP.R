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
library(ggplot2)
library(ggpubr)
library(ggrastr)
})
source("input_data/Utils.R")
source("input_data/IKAP_Seurat3.R")

# Load the data
load("Reanalysis_Scaled/Hippo_Data_Scaling_30pcs.RData")

Hippo_IKAP <- Hippo_Scaled_30

Hippo_IKAP <- IKAP(Hippo_IKAP, out.dir = "Reanalysis_Scaled/IKAP")

save(Hippo_IKAP, file = "Reanalysis_Scaled/Hippo_Scaled_30_IKAP.RData")

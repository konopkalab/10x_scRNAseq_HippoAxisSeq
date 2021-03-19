suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(emmeans)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(patchwork)
library(miloR)
})

sce <- readRDS("input_data/HippoAxis_Filt.rds")
reducedDimNames(sce) <- c("pca.corrected","umap")
assays(sce) <- assays(sce)["counts"]


processing_milo <- function(object, k = 30, d = 30, prop = 0.1 ) {
  	object %>% 
	Milo() %>%
	buildGraph(k = k , d = d, reduced.dim = "pca.corrected") %>%
	makeNhoods(prop = prop, k = k, d=d, refined = TRUE, reduced_dims = "pca.corrected") %>%
	countCells(meta.data = data.frame(colData(.)), sample="orig.ident") %>%
	calcNhoodDistance(d=d, reduced.dim = "pca.corrected")
}


sce_milo <- processing_milo(sce,k=30,d=30,prop=0.1) 


sce_design <- data.frame(colData(sce_milo))[,c("orig.ident","Axis", "batch","age","sex","dur","version")] %>%
				dplyr::rename(ID = orig.ident) %>%
				mutate(batch = as.factor(batch), Axis = as.factor(Axis),sex = as.factor(sex),version = as.factor(version)) %>%
				distinct() %>%
				remove_rownames() %>%
				column_to_rownames("ID")

da_results <- testNhoods(sce_milo, design = ~ Axis, design.df = sce_design)

sce_milo <- buildNhoodGraph(sce_milo)

umap_pl <- plotReducedDim(sce_milo, dimred = "umap", colour_by="Axis", text_by = "ident", text_size = 3) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sce_milo, da_results, layout="umap",alpha=0.2)

pdf("Reviewer1/MILO_Results.pdf",width=12,height=8)
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
dev.off()

pdf("Reviewer1/MILO_Results_Graph.pdf",width=8,height=8)
print(nh_graph_pl)
dev.off()

da_results <- annotateNhoods(sce_milo, da_results, coldata_col = "ident") 

openxlsx::write.xlsx(da_results, file = "Reviewer1/Cell_Milo_Results.xlsx", colNames = TRUE, borders = "columns")

pdf("Reviewer1/MILO_Results_Beeswarm.pdf",width=6,height=8)
plotDAbeeswarm(da_results, group.by = "ident",alpha = 0.2)
dev.off()









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

dir.create("Reanalysis_SCTransform")

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

sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

sce.seurat <- SplitObject(sce.seurat, split.by = "batch")
#sce.seurat <- sce.seurat[c("Batch2","Batch3","Batch5","Batch6")]

# Transform and regress out covariates
for (i in 1:length(sce.seurat)) {
    sce.seurat[[i]] <- SCTransform(sce.seurat[[i]], 
				    						vars.to.regress = c("nUMI","percent.mt","age","dur","SeizFreq"), 
											verbose = FALSE)
    }

integ_features <- SelectIntegrationFeatures(object.list = sce.seurat, 
											nfeatures = 3000) 

sce.seurat <- PrepSCTIntegration(object.list = sce.seurat, 
											anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = sce.seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

Hippo_Integrated <- IntegrateData(
								anchorset = integ_anchors,
								new.assay.name = "integrated",
								normalization.method = "SCT",
								dims = 1:30,
								k.weight = 100,
								sd.weight = 1,
								do.cpp = TRUE,
								eps = 0.5,
								verbose = TRUE
								)

Hippo_Integrated <- RunPCA(object = Hippo_Integrated, 
								features=NULL, 
								weight.by.var = TRUE, 
								ndims.print = 1:5, 
								nfeatures.print = 30, 
								npcs = 30, 
								reduction.name = "pca")

Hippo_Integrated <- FindNeighbors(object = Hippo_Integrated, 
										reduction = "pca", 
										dims = 1:30, 
										nn.eps = 0.5)

Hippo_Integrated <- FindClusters(object = Hippo_Integrated, 
										resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4), 
										algorithm = 1,
										n.iter = 1000)

pdf("Reanalysis_SCTransform/Hippo_Integrated_Clustree.pdf", width = 12, height = 6)
clustree(Hippo_Integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Select resolution and run UMAP
#Idents(object = Hippo_Integrated) <- "integrated_snn_res.0.8"

Hippo_Integrated <- RunUMAP(object = Hippo_Integrated, 
										reduction = "pca", 
										dims = 1:30)

pdf("Reanalysis_SCTransform/Hippo_Integrated_UMAP.pdf", width = 7, height = 6)
DimPlot(object = Hippo_Integrated, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

pdf("Reanalysis_SCTransform/Hippo_Integrated_UMAP_Leuv_MultiRes.pdf", width = 20, height = 20)
cowplot::plot_grid(ncol = 3,
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.0.2")+ggtitle("louvain_0.2"),
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.0.4")+ggtitle("louvain_0.4"),
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.0.6")+ggtitle("louvain_0.6"),
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.0.8")+ggtitle("louvain_0.8"),
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.1")+ggtitle("louvain_1.0"),
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.1.2")+ggtitle("louvain_1.2"),
  DimPlot(Hippo_Integrated, reduction = "umap", group.by = "integrated_snn_res.1.4")+ggtitle("louvain_1.4")
)
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(Hippo_Integrated) <- "RNA"
Hippo_Integrated <- NormalizeData(object = Hippo_Integrated, 
						normalization.method = "LogNormalize", 
						scale.factor = 10000)

save(Hippo_Integrated, file = "Reanalysis_SCTransform/Hippo_Data_Integrated.RData")


# Chose the best resolution for the data
Hippo_Integrated2 <- SCTransform_ChooseRes(Hippo_Integrated,n.pcs=30,res.low=0.1,res.high=1.5,res.n=20)
save(Hippo_Integrated2, file = "Reanalysis_SCTransform/Hippo_Data_Integrated_ChooseRes.RData")

sessionInfo()

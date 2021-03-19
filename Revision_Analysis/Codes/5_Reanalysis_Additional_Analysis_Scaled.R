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
source("utils/Utils.R")

# Load the data
load("Reanalysis_Scaled/Hippo_Data_Scaling_30pcs.RData")

Hippo_Scaled <- Hippo_Scaled_30

# Remove res columns and clusters
Hippo_Scaled@meta.data$RNA_snn_res.0.8 <- NULL
Hippo_Scaled@meta.data$seurat_clusters <- NULL


# Multiple resolution
Hippo_Scaled <- FindClusters(object = Hippo_Scaled, 
								resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8), 
								algorithm = 1,
								n.iter = 1000)

pdf("Reanalysis_Scaled/Hippo_Scaled_Clustree.pdf", width = 12, height = 6)
clustree(Hippo_Scaled@meta.data, prefix = "RNA_snn_res.", node_colour = "sc3_stability")
dev.off()


pdf("Reanalysis_Scaled/Hippo_Scaled_JackStrawPlot.pdf", width = 6, height = 5)
JackStrawPlot(Hippo_Scaled, dims = 1:30, reduction = "pca")
dev.off()

pdf("Reanalysis_Scaled/Hippo_Scaled_ElbowPlot.pdf", width = 4, height = 4)
ElbowPlot(Hippo_Scaled, reduction = "pca", ndims = 30)
dev.off()


pdf("Reanalysis_Scaled/Hippo_Scaled_UMAP_Leuv_MultiRes.pdf", width = 20, height = 20)
cowplot::plot_grid(ncol = 3,
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.0.2")+ggtitle("louvain_0.2"),
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.0.4")+ggtitle("louvain_0.4"),
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.0.6")+ggtitle("louvain_0.6"),
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.0.8")+ggtitle("louvain_0.8"),
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.1")+ggtitle("louvain_1.0"),
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.1.2")+ggtitle("louvain_1.2"),
  DimPlot(Hippo_Scaled, reduction = "umap", group.by = "RNA_snn_res.1.4")+ggtitle("louvain_1.4")
)
dev.off()


# Select Resolution
Idents(object = Hippo_Scaled) <- "RNA_snn_res.0.8"

Hippo_Scaled@meta.data$seurat_clusters <- Hippo_Scaled@meta.data$RNA_snn_res.0.8

# Visualized multiple UMAPs
umap <- as.data.frame(Embeddings(Hippo_Scaled, reduction = "umap"))
meta <- as.data.frame(Hippo_Scaled@meta.data)

metadata <- cbind(umap,meta)%>% 
	group_by(seurat_clusters) %>% 
	mutate(N = n()) %>%
  	ungroup() %>% 
  	as.data.frame()

a <- ggplot(metadata, aes(x=umap_1, y=umap_2)) +
      ggrastr::geom_point_rast(aes(colour = age),size=0.2,alpha=0.3) +
      theme_classic() + 
      gradient_color("red")


b <- ggplot(metadata, aes(x=umap_1, y=umap_2)) +
      ggrastr::geom_point_rast(aes(colour = batch),size=0.2,alpha=0.3) +
      theme_classic() #+ 
      #gradient_color("blue")


c <- ggplot(metadata, aes(x=umap_1, y=umap_2)) +
      ggrastr::geom_point_rast(aes(colour = sex),size=0.2,alpha=0.3) +
      theme_classic()


d <- ggplot(metadata, aes(x=umap_1, y=umap_2)) +
      ggrastr::geom_point_rast(aes(colour = Axis),size=0.2,alpha=0.3)+
      theme_classic() + 
      scale_colour_manual(values = c("gray60", "purple"))

plot2by2 <- cowplot::plot_grid(a, b,c,d,labels=c("A", "B","C","D"), ncol = 2)
cowplot::save_plot("Reanalysis_Scaled/Hippo_Data_Scaling_UMAPs_Covariates.pdf", plot2by2, ncol = 2,base_height=6,base_width=4)


# Visualize
pdf("Reanalysis_Scaled/Hippo_Data_Scaling_UMAPs_Clusters.pdf", width = 5, height = 5)

label <- data.frame(seurat_clusters=levels(metadata$seurat_clusters),label=levels(metadata$seurat_clusters))

label_2 <- metadata %>% 
  group_by(seurat_clusters) %>% 
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2),N = n()) %>% 
  left_join(label) %>%
  as.data.frame()

ggplot(metadata, aes(x=umap_1, y=umap_2)) +
ggrastr::geom_point_rast(aes(colour = seurat_clusters),size=0.5) +
ggrepel::geom_text_repel(data = label_2, aes(label = label),
							color = "black",
							#fontface = 'bold',
							segment.colour = "grey60",
						    box.padding = unit(0.25, "lines"),
						    point.padding = unit(0.5, "lines"),
						    nudge_x = .15,
						    nudge_y = 1,
						    size = 2.5) + 
		theme_minimal() +
		theme(legend.position="none")
dev.off()

# UMAP Alternative by Axis
pdf("Reanalysis_Scaled/Hippo_Data_Scaling_UMAP_Axis.pdf", width = 6, height = 6)
ggplot(metadata, aes(x=umap_1, y=umap_2, color=Axis)) +
ggrastr::geom_point_rast(size=0.1) +
theme_minimal() +
theme(legend.position="none") + 
scale_colour_manual(values = c("gray60", "purple"))
dev.off()

# UMAP Alternative by Subject
pdf("Reanalysis_Scaled/Hippo_Data_Scaling_UMAPs_BySubject.pdf", width = 6, height = 6)
ggplot(metadata, aes(x=umap_1, y=umap_2, color=orig.ident)) +
ggrastr::geom_point_rast(size=0.5) +
facet_wrap(.~orig.ident,ncol=3) +
theme_classic() +
theme(legend.position="none")
dev.off()

# UMAP Alternative by Axis
pdf("Reanalysis_Scaled/Hippo_Data_Scaling_UMAPs_ByAxis.pdf", width = 6, height = 4)
ggplot(metadata, aes(x=umap_1, y=umap_2, color=orig.ident)) +
ggrastr::geom_point_rast(size=0.5) +
facet_wrap(.~Axis,ncol=3) +
theme_classic() +
theme(legend.position="none")
dev.off()

pdf("Reanalysis_Scaled/Hippo_Data_Scaling_UMAP_StackedBar.pdf", width = 8, height = 4)
ggplot(metadata, aes(x=seurat_clusters, fill=Axis)) + 
geom_bar(position = "fill")+
scale_y_continuous(labels = scales::percent_format())+
theme_classic() +
rotate_x_text(angle = 45) + 
geom_hline(yintercept=0.5, linetype="dashed", color = "black") + 
scale_fill_manual(values = c("gray60", "purple")) +
ylab("") + xlab("")
dev.off()

pdf("Reanalysis_Scaled/Hippo_Data_Scaling_UMAP_StackedBar_Ordered.pdf", width = 8, height = 4)
metadata %>% 
  mutate(x = forcats::fct_reorder(seurat_clusters, as.numeric(Axis), mean)) %>% 
ggplot(aes(x=x, fill=Axis)) + 
geom_bar(position = "fill")+
scale_y_continuous(labels = scales::percent_format())+
theme_classic() +
rotate_x_text(angle = 45) + 
geom_hline(yintercept=0.75, linetype="dashed", color = "black") + 
geom_hline(yintercept=0.5, linetype="dashed", color = "black") + 
geom_hline(yintercept=0.25, linetype="dashed", color = "black") + 
scale_fill_manual(values = c("gray60", "purple")) +
ylab("") + xlab("")
dev.off()



genes <- c("RBFOX3",
           "SLC17A7",
           "MAML2",
           "SV2B",
           "SATB2",
           "HS3ST4",
           "TYRO3",
           "PFKP",
           "FN1",
           "GAD1",
           "GAD2")

pdf("Reanalysis_Scaled/Markers_Violin.pdf",width=8,height=12)
StackedVlnPlot(Hippo_Scaled,features = genes) 
dev.off()













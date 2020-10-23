#integrate_anchor, norm, scale to nUMI, percent "age", "sex", "dur","percent.mt", "batch", "nCount_RNA"

rm(list=ls()) 
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

load("A76_sing.RData")
load("P76_sing.RData")
load("P82_sing.RData")
load("A82_sing.RData")
load("P56_sing.RData")
load("A56_sing.RData")
load("P67_sing.RData")
load("A67_sing.RData")
load("P57_sing.RData")
load("A57_sing.RData")


combined57 <- merge(x=P57New, y=A57New, add.cell.ids = c("P57", "A57"))
combined67 <- merge(x=P67New, y=A67New, add.cell.ids = c("P67", "A67"))
combined56 <- merge(x=P56New, y=A56New, add.cell.ids = c("P56", "A56"))
combined76 <- merge(x=P76New_S, y=A76New_S, add.cell.ids = c("P76", "A76"))
combined82 <- merge(x=P82New_S, y=A82New_S, add.cell.ids = c("P82", "A82"))

 
combined57 <- NormalizeData(object = combined57)
combined57 <- FindVariableFeatures(object = combined57, selection.method = "vst", nfeatures = 2000)

combined67  <- NormalizeData(object = combined67)
combined67  <- FindVariableFeatures(object = combined67, selection.method = "vst", nfeatures = 2000)

combined56 <- NormalizeData(object = combined56)
combined56 <- FindVariableFeatures(object = combined56, selection.method = "vst", nfeatures = 2000)

combined76 <- NormalizeData(object = combined76)
combined76 <- FindVariableFeatures(object = combined76, selection.method = "vst", nfeatures = 2000)

combined82 <- NormalizeData(object = combined82)
combined82 <- FindVariableFeatures(object = combined82, selection.method = "vst", nfeatures = 2000)

combined86 <- NormalizeData(object = combined86)
combined86 <- FindVariableFeatures(object = combined86, selection.method = "vst", nfeatures = 2000)


hippo.anchors <- FindIntegrationAnchors(object.list = list(combined86, combined57, combined67, combined56, combined76, combined82), dims = 1:30)
hippo.combined <- IntegrateData(anchorset = hippo.anchors, dims = 1:30)

DefaultAssay(hippo.combined) <- "integrated"


hippo.combined <- ScaleData(object = hippo.combined, verbose = FALSE, vars.to.regress = c("age", "sex", "dur","percent.mt", "batch", "nCount_RNA", "version"))

hippo.combined <- RunPCA(object = hippo.combined, npcs = 100, verbose = TRUE)

hippo.combined<- JackStraw(object = hippo.combined, num.replicate = 100, dims = 100)
hippo.combined<- ScoreJackStraw(object = hippo.combined, dims = 1:100)

plotJS <-JackStrawPlot(object = hippo.combined, dims = 1:100)
ggsave("integ/Jackstaw.pdf", plot = plotJS, width = 12, height =8, units = "in", dpi = 300)

save(hippo.combined, file = "integ/integ_JS.RData")












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
})

source("Utils.R")
dir.create("Reviewer2")

sce <- readRDS("input_data/HippoAxis_Filt.rds")
#reducedDimNames(sce) <- c("pca.corrected","umap")

################### 
# Add information #
###################

metadata <- colData(sce) %>%
        as.data.frame() %>%
        mutate(PC1 = reducedDim(sce, "PCA")[,1], PC2 = reducedDim(sce, "PCA")[,2]) %>%
        mutate(UMAP1 = reducedDim(sce, "UMAP")[,1], UMAP2 = reducedDim(sce, "UMAP")[,2]) %>%
        select(PC1,PC2, UMAP1, UMAP2,orig.ident, Axis,age,sex, dur,ident,batch,version,fre) %>%
        mutate(ident = as.character(ident)) %>%
        mutate(ident = if_else(ident == "29", "CApyr2",ident)) %>%
        dplyr::rename(Subjects = orig.ident, Age = age, Sex = sex, EpDur = dur, Cell=ident, Batch = batch,Version=version,Frequency=fre)



a <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
      ggrastr::geom_point_rast(aes(colour = Age),size=0.2,alpha=0.3) +
      theme_classic() + 
      gradient_color("red")


b <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
      ggrastr::geom_point_rast(aes(colour = EpDur),size=0.2,alpha=0.3) +
      theme_classic() + 
      gradient_color("blue")


c <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
      ggrastr::geom_point_rast(aes(colour = Sex),size=0.2,alpha=0.3) +
      theme_classic() + 
      scale_colour_manual(values = c("black", "red"))


d <- ggplot(metadata, aes(x=UMAP1, y=UMAP2)) +
      ggrastr::geom_point_rast(aes(colour = Axis),size=0.2,alpha=0.3)+
      theme_classic() + 
      scale_colour_manual(values = c("gray60", "purple"))

plot2by2 <- cowplot::plot_grid(a, b,c,d,labels=c("A", "B","C","D"), ncol = 2)
cowplot::save_plot("Reviewer2/UMAPs_Covariates.pdf", plot2by2, ncol = 2,base_height=6,base_width=4)


# 
pdf("Reviewer2/Heatmap_CellPerAge.pdf",width=6,height=6)
table(metadata$Cell,metadata$Age,metadata$Axis) %>%
as.data.frame() %>%
group_by(Var3) %>%
mutate(Scale = (Freq - min(Freq))/diff(range(Freq))) %>%
ggplot(aes(x = Var2,y = Var1)) + 
    geom_tile(aes(fill = Scale)) + 
    geom_text(aes(y=Var1,label=Freq)) +
        xlab("") + ylab("")+
    #coord_equal(ratio = 5) + 
    scale_fill_gradient(low="white", high="firebrick3")+
    #viridis::scale_fill_viridis(option="plasma") + 
    ggthemes::theme_tufte(base_family="Helvetica") +
    facet_wrap(~Var3) + theme(legend.position = "none") 
dev.off()

pdf("Reviewer2/Heatmap_CellPerSex.pdf",width=5,height=6)
table(metadata$Cell,metadata$Sex,metadata$Axis) %>%
as.data.frame() %>%
group_by(Var3) %>%
mutate(Scale = (Freq - min(Freq))/diff(range(Freq))) %>%
ggplot(aes(x = Var2,y = Var1)) + 
    geom_tile(aes(fill = Scale)) +
        geom_text(aes(y=Var1,label=Freq)) +
        xlab("") + ylab("")+
    #coord_equal(ratio = 5) + 
    scale_fill_gradient(low="white", high="blue")+
    #viridis::scale_fill_viridis(option="plasma") + 
    ggthemes::theme_tufte(base_family="Helvetica") +
    facet_wrap(~Var3)
dev.off()


# Variance explained by pseudobulk analysis
# Muscat aggregate data
sce_aggregate <- aggregateData(sce, 
    				assay = "counts", fun = "sum",
    				by = "orig.ident")

names(assays(sce_aggregate)) <- "counts"

dataExp <- as.data.frame(assays(sce_aggregate)$counts)
cpm <- apply(dataExp, 2, function(x) x/sum(as.numeric(x)) * 10^6)
logcpm <- log2(cpm+1)


# Create a metadata
# Group
tmp2 <- colData(sce) %>%
        as.data.frame()  %>% 
			group_by(Axis, orig.ident) %>% 
			tally() %>%
			dplyr::rename(NofCells = n)		

tmp3 <- colData(sce) %>%
        as.data.frame()  %>%
        select(orig.ident, batch,age,sex,dur,version,fre) %>%
        group_by(orig.ident, batch,age,sex,dur,version,fre) %>%
        summarise()

metadata_pseudobulk <- Reduce(dplyr::full_join, list(tmp2,tmp3)) %>%
		as.data.frame() %>%
		dplyr::rename(Subjects = orig.ident) %>%
	    mutate(batch = as.factor(batch), sex = as.factor(sex),version=as.factor(version)) %>%
	    column_to_rownames("Subjects")

# variance explained by covaraites
var <- VarExp(logcpm,metadata_pseudobulk,10,FALSE)
pdf("Reviewer2/Variance_Explained_Pseudobulk.pdf",width=5,height=4)
plotVarExp(var,"Variance Explained")
dev.off()

## Variance explained by cell
sce_aggregate2 <- aggregateData(sce, 
    				assay = "counts", fun = "sum",
    				by = c("Cell","orig.ident"))


names <- names(assays(sce_aggregate2))



dataExp <- list()
logcpm <- list()
for(i in 1:length(names)){
dataExp[[i]] <- as.data.frame(assays(sce_aggregate2)[[i]])
logcpm[[i]] <- log2(apply(dataExp[[i]], 2, function(x) x/sum(as.numeric(x)) * 10^6) + 1) %>% as.data.frame()
}

names(logcpm) <- names


# Group
tmp1 <- colData(sce) %>%
        as.data.frame() %>% 
			group_by(Axis, Cell,orig.ident) %>% 
			tally() %>%
			dplyr::rename(total = n)		 	

tmp2 <- colData(sce) %>%
        as.data.frame() %>% 
			group_by(Axis, orig.ident) %>% 
			tally() %>%
			dplyr::rename(NofCells = n)		

tmp3 <- colData(sce) %>%
        as.data.frame() %>%
        select(orig.ident, batch,age,sex,dur,version,fre) %>%
        group_by(orig.ident, batch,age,sex,dur,version,fre) %>%
        summarise()

metadata_pseudobulk_cell <- Reduce(dplyr::full_join, list(tmp1, tmp2,tmp3)) %>%
							mutate(Others = NofCells-total) %>% 
							as.data.frame() %>%
							dplyr::rename(Subjects = orig.ident) %>%							
					    	mutate(batch = as.factor(batch), sex = as.factor(sex),version=as.factor(version)) %>%
					    	select(Axis, Cell, Subjects,NofCells,batch, age,sex,dur,version,fre) %>%
					    	split(.,.$Cell)

metadata_pseudobulk_cell <- map(metadata_pseudobulk_cell, ~ (.x %>% `rownames<-`(.x$Subjects) %>% select(-Subjects,-Cell))) 

# remove olig6 because only 2 subjects
metadata_pseudobulk_cell["Olig5"] <- NULL
logcpm["Olig5"] <- NULL
metadata_pseudobulk_cell["Den.Gyr3"] <- NULL
logcpm["Den.Gyr3"] <- NULL

var <- list()
logcpm_tmp <- list()
for(i in 1:length(logcpm)){
logcpm_tmp[[i]] <- logcpm[[i]][,colnames(logcpm[[i]]) %in% rownames(metadata_pseudobulk_cell[[i]])]
var[[i]] <- VarExp(logcpm_tmp[[i]],metadata_pseudobulk_cell[[i]],10,FALSE)
}

names(var) <- names(logcpm)

df <- do.call(rbind,var) %>% 
		as.data.frame() %>%
		rownames_to_column("Cell") %>%
		pivot_longer(!Cell, values_to="WAPV",names_to="Covariates")

pdf("Reviewer2/Variance_Explained_Pseudobulk_ByCell.pdf",width=6,height=6)
 ggbarplot(df, "Covariates", "WAPV",
   fill = "Cell", color = "Cell") +
   facet_wrap(.~Cell) +
	ggthemes::theme_tufte(base_family="Helvetica") + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none")  
dev.off()


pdf("Reviewer2/Variance_Explained_Pseudobulk_ByCell_boxplot.pdf",width=5,height=4)
 ggboxplot(df, "Covariates", "WAPV",
    color = "Covariates",
    add = "jitter") +
    ggthemes::theme_tufte(base_family="Helvetica") + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none") 
dev.off()

# Fold Change Distribution
fc <- list()
logcpm_tmp <- list()
for(i in 1:length(logcpm)){
logcpm_tmp[[i]] <- logcpm[[i]][,colnames(logcpm[[i]]) %in% rownames(metadata_pseudobulk_cell[[i]])]
fc[[i]] <- log2(rowMeans(logcpm_tmp[[i]][metadata_pseudobulk_cell[[i]]$Axis == "anterior"])/rowMeans(logcpm_tmp[[i]][metadata_pseudobulk_cell[[i]]$Axis == "posterior"]))
}

names(fc) <- names(logcpm)

x <- do.call(cbind,fc) %>%
        na.omit() %>% 
        as.data.frame() %>%
        rownames_to_column("Gene") %>%
        pivot_longer(!Gene,names_to="Cell",values_to="FoldChange")

pdf("Reviewer2/Fold_Change.pdf",width=10,height=6)
gghistogram(x, x = "FoldChange", fill = "Cell",rug = TRUE) + 
            facet_wrap(.~Cell,scales="free") +
    ggthemes::theme_tufte(base_family="Helvetica") + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none") +
    geom_vline(xintercept = 0.3, linetype="dashed", color = "black", size=1) +
    geom_vline(xintercept = -0.3, linetype="dashed", color = "black", size=1) +
    xlab("fc(anterior/posterior)")
dev.off()



# Get variance per gene
vars <- getVarianceExplained(sce, variables=c("Axis", "age","sex","dur","batch","Cell","fre"))
openxlsx::write.xlsx(vars, file = "Reviewer2/Variance_Explained_Statistics_AllCells.xlsx", colNames = TRUE, rowNames = TRUE,borders = "columns")



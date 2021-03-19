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
	    select(-NofCells) %>%
	    column_to_rownames("Subjects")

# variance explained by covaraites
var <- VarExp(logcpm,metadata_pseudobulk,10,FALSE)
pdf("Reviewer2/Variance_Explained_Pseudobulk_NoCell.pdf",width=5,height=4)
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
					    	select(Axis, Cell, Subjects,batch, age,sex,dur,version,fre) %>%
					    	split(.,.$Cell)

metadata_pseudobulk_cell <- map(metadata_pseudobulk_cell, ~ (.x %>% `rownames<-`(.x$Subjects) %>% select(-Subjects,-Cell))) 

# remove olig6 because only 2 subjects
metadata_pseudobulk_cell["Olig5"] <- NULL
logcpm["Olig5"] <- NULL
metadata_pseudobulk_cell["Den.Gyr3"] <- NULL
logcpm["Den.Gyr3"] <- NULL
metadata_pseudobulk_cell["Micro3"] <- NULL
logcpm["Micro3"] <- NULL
metadata_pseudobulk_cell["OPC2"] <- NULL
logcpm["OPC2"] <- NULL

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

pdf("Reviewer2/Variance_Explained_Pseudobulk_ByCell_NoCell.pdf",width=6,height=6)
 ggbarplot(df, "Covariates", "WAPV",
   fill = "Cell", color = "Cell") +
   facet_wrap(.~Cell) +
	ggthemes::theme_tufte(base_family="Helvetica") + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none")  
dev.off()


pdf("Reviewer2/Variance_Explained_Pseudobulk_ByCell_NoCell_boxplot.pdf",width=5,height=4)
 ggboxplot(df, "Covariates", "WAPV",
    color = "Covariates",
    add = "jitter") +
    ggthemes::theme_tufte(base_family="Helvetica") + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none") 
dev.off()



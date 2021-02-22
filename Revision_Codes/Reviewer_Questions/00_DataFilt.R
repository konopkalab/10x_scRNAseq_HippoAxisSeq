suppressPackageStartupMessages({
library(tidyverse)
library(SingleCellExperiment)
library(BiocParallel)
library(scater)
library(Seurat)
library(data.table)
})

load("input_data/Fig1_final.RData")


sce_hippo <- as.SingleCellExperiment(hippo.com)

################### 
# Add information #
###################
temporary <- colData(sce_hippo) %>%
				as.data.frame() %>%
                rownames_to_column("TMP") %>%
                select(-DF.classifications_0.25_0.01_3770,-DF.classifications_0.25_0.005_3794,-DF.classifications_0.25_0.005_2981,
                        -DF.classifications_0.25_0.01_1582,-pANN_0.25_0.01_3770,-pANN_0.25_0.005_2981,
                        -pANN_0.25_0.005_3794,-pANN_0.25_0.01_1582,-integrated_snn_res.1,
                        -integrated_snn_res.0.4,-integrated_snn_res.0.6,-integrated_snn_res.0.8,
                        -integrated_snn_res.1.2,-integrated_snn_res.1.4,-integrated_snn_res.1.6) %>%
                mutate(Cell = as.factor(ident)) %>%
                mutate(Class = case_when(grepl("Pyr|Den", Cell) ~ "Glutamatergic", 
                        grepl("In", Cell) ~ "Gabaergic",
                        grepl("Olig|Astro|OPC|Endo|Micro", Cell) ~ "Glia")) %>%
                mutate(Definition = case_when(grepl("Pyr", Cell) ~ "PyramidalN",
                        grepl("Den", Cell) ~ "GranuleN",
                        grepl("Olig", Cell) ~ "Oligodendrocytes",
                		grepl("OPC", Cell) ~ "OPC",
                        grepl("In", Cell) ~ "Inhibitory",
                        grepl("Astro", Cell) ~ "Astrocytes",
                        grepl("Endo", Cell) ~ "Endothelial",
                        grepl("Micro", Cell) ~ "Microglia")) %>%
                dplyr::rename(Axis = group) %>%
                column_to_rownames("TMP") %>% 
                DataFrame()

colData(sce_hippo) <- temporary

#assays(sce_hippo) <- assays(sce_hippo)["logcounts"]

saveRDS(sce_hippo,"input_data/HippoAxis_Filt.rds")


####
# Neurons

load("input_data/All_neurons_renamed.RData")


sce_neurons <- as.SingleCellExperiment(Neurons)

################### 
# Add information #
###################
temporary <- colData(sce_neurons) %>%
                                as.data.frame() %>%
                rownames_to_column("TMP") %>%
                select(-DF.classifications_0.25_0.01_3770,-DF.classifications_0.25_0.005_3794,-DF.classifications_0.25_0.005_2981,
                        -DF.classifications_0.25_0.01_1582,-pANN_0.25_0.01_3770,-pANN_0.25_0.005_2981,
                        -pANN_0.25_0.005_3794,-pANN_0.25_0.01_1582,-integrated_snn_res.1,
                        -integrated_snn_res.0.4,-integrated_snn_res.0.6,-integrated_snn_res.0.5) %>%
                mutate(Cell = as.factor(ident)) %>%
                mutate(Class = case_when(grepl("CA|Gra", Cell) ~ "Glutamatergic", 
                        grepl("In", Cell) ~ "Gabaergic",
                        grepl("Sub", Cell) ~ "Subiculum")) %>%
                mutate(Definition = case_when(grepl("CA1", Cell) ~ "CA1",
                        grepl("CA2", Cell) ~ "CA2",
                        grepl("CA3", Cell) ~ "CA3",
                                grepl("Gra", Cell) ~ "GranuleN",
                        grepl("In", Cell) ~ "Inhibitory",
                        grepl("Sub", Cell) ~ "Subiculum")) %>%
                dplyr::rename(Axis = group) %>%
                column_to_rownames("TMP") %>% 
                DataFrame()

colData(sce_neurons) <- temporary

#assays(sce_hippo) <- assays(sce_hippo)["logcounts"]

saveRDS(sce_neurons,"input_data/HippoAxis_Neurons_Filt.rds")



library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(WGCNA)
library(tidyverse)

files=list.files(pattern="STATISTICS_10kb")
tmp=as.data.frame(lapply(files,read.table,sep="\t",header=T)[[1]])

df <- tmp %>%
      select(VARIABLE,Sample,BETA,P)

#tmp <- tmp[c(1,7,8)]
#df <- melt(tmp)
df$log <- -log10(df$P)
df$Sample <- factor(df$Sample,levels=rev(levels(df$Sample)))


df <- df %>%
      group_by(VARIABLE) %>%
      mutate(FDR = p.adjust(P,method="BH"),BONF = p.adjust(P,method="bonferroni")) %>%
      mutate(logFDR = -log10(FDR),logBONF = -log10(BONF)) %>%
      as.data.frame()
df$log[df$log < 1.3] <- NA
df$logFDR[df$logFDR < 1.3] <- NA
df$logBONF[df$logBONF < 1.3] <- NA

openxlsx::write.xlsx(df, 
                     file = "Magma_Statistics_All_Update.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="P-values")

df <- df %>% 
      separate(VARIABLE, c("Cluster","Class"), "_",remove=FALSE) %>%
      filter(grepl("Cluster",VARIABLE))

df$VARIABLE <- factor(df$VARIABLE,levels = c(
"Cluster0_Anterior", "Cluster0_Posterior",  "Cluster11_Anterior", "Cluster11_Posterior", "Cluster3_Anterior", "Cluster3_Posterior",  "Cluster13_Anterior", "Cluster13_Posterior", "Cluster16_Anterior", "Cluster16_Posterior", "Cluster1_Anterior", "Cluster1_Posterior", "Cluster10_Anterior", "Cluster10_Posterior", 
"Cluster15_Anterior", "Cluster15_Posterior", "Cluster14_Anterior", "Cluster14_Posterior", "Cluster6_Anterior", "Cluster6_Posterior", "Cluster2_Anterior", "Cluster2_Posterior", "Cluster7_Anterior", "Cluster7_Posterior", "Cluster8_Anterior", "Cluster8_Posterior", "Cluster9_Anterior", "Cluster9_Posterior"
))


df$Sample <- factor(df$Sample,levels = rev(c(
"EPILEPSY_2018", "ADHD_2017", "ASD_2017", "BIP_2018", "MDD_2018", "ALZ_2019", "SZ_2018", "tCognFun", 
"tEduATT", "tIntelligence")))

pdf("MAGMA_BUBBLE_10kb_Clusters.pdf",width=8,height=3.5)
ggscatter(df, 
  			x = "Sample",
  			y = "VARIABLE",
   			size="log",
   			color="Class",
        palette=c("grey60","purple"),
   			alpha = 0.8,
   			xlab = "",ylab = "",) +
   			theme_minimal() + 
   			#gradient_color(c("red", "darkred"))+
   			rotate_x_text(angle = 45)+
        coord_flip()+
        scale_size(range = c(2, 8))

dev.off()



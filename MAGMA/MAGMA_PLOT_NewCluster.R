library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(WGCNA)
library(tidyverse)


files=list.files(pattern="_STATISTICS_NewClusters")
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
                     file = "Magma_Statistics_Clusters_Update.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="P-values")

df$Beta<- ifelse(is.na(df$log), df$log, df$BETA)

df1 <- df %>%
        filter(!(grepl('z', Sample)))


df$VARIABLE <- factor(df$VARIABLE, levels = c("Cluster_4", "Cluster_0", "Cluster_12", "Cluster_11", "Cluster_3", "Cluster_13", "Cluster_16", "Cluster_1", 
"Cluster_10","Cluster_5", "Cluster_15", "Cluster_14", "Cluster_17", "Cluster_6", "Cluster_2", "Cluster_7", "Cluster_8","Cluster_9"
))

df$Sample <- factor(df$Sample,levels = rev(c(
"EPILEPSY_2018", "ADHD_2017", "ASD_2017", "BIP_2018", "MDD_2018", "ALZ_2019", "SZ_2018", "tCognFun", 
"tEduATT", "tIntelligence", "zBMI", "zCHD", "zDIAB", "zHGT", "zOSTEO")))

pdf("MAGMA_BUBBLE_10kb_NewCluster.pdf",width=6,height=4)
ggscatter(df, 
  			x = "Sample",
  			y = "VARIABLE",
   			size="Beta",
   			color="log",
   			alpha = 0.8,
   			xlab = "",
            ylab = "",) +
   			theme_minimal() + 
   			gradient_color(c("red", "darkred"))+
   			rotate_x_text(angle = 45)+
        coord_flip()#+
            #scale_size(range = c(2, 10))
            #gradient_color(c("lightblue","darkblue"))
dev.off()


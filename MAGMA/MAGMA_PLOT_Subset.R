library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(WGCNA)
library(tidyverse)
files=list.files(pattern="STATISTICS_10kb")
tmp=as.data.frame(lapply(files,read.table,sep="\t",header=T)[[1]])

tmp <- tmp[c(1,7,8)]
tmp <- tmp[-grep("z",tmp$Sample),]
df <- melt(tmp)
f$log <- -log10(df$P)
df$Sample <- factor(df$Sample,levels=rev(levels(df$Sample)))


df <- df %>%
      group_by(VARIABLE) %>%
      mutate(FDR = p.adjust(P,method="BH"),BONF = p.adjust(P,method="bonferroni")) %>%
      mutate(logFDR = -log10(FDR),logBONF = -log10(BONF)) %>%
      as.data.frame()
df$log[df$log < 1.3] <- NA
df$logFDR[df$logFDR < 1.3] <- NA
df$logBONF[df$logBONF < 1.3] <- NA

int <- c("CA1_Anterior","CA1_Posterior",
                             "DG_Anterior","DG_Posterior",
                             "IN_Anterior","IN_Posterior")

df <- df %>% 
      separate(VARIABLE, c("Cluster","Class"), "_",remove=FALSE) %>%
      filter(VARIABLE %in% int)

df$VARIABLE <- factor(df$VARIABLE, levels = c("DG_Anterior", "DG_Posterior", "CA1_Anterior", "CA1_Posterior", "IN_Anterior", "IN_Posterior"))

df$Sample <- factor(df$Sample,levels = rev(c(
"EPILEPSY_2018", "ADHD_2017", "ASD_2017", "BIP_2018", "MDD_2018", "ALZ_2019", "SZ_2018", "tCognFun", 
"tEduATT", "tIntelligence", "zBMI", "zCHD", "zDIAB", "zHGT", "zOSTEO")))




pdf("MAGMA_BUBBLE_10kb_Subset.pdf",width=4,height=3.5)
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



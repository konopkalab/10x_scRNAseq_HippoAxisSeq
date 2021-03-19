suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(emmeans)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(data.table)
library(robustlmm)
})

dir.create("Reviewer1")

sce <- readRDS("input_data/HippoAxis_Filt.rds")

################### 
# Add information #
###################
metadata <- colData(sce) %>%
				as.data.frame()

pdf("Reviewer1/Heatmap_CellClassPerDonor.pdf",width=6,height=6)
table(metadata$Definition,metadata$donor,metadata$Axis) %>%
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
    facet_wrap(~Var3) + theme(legend.position = "none") +
    rotate_x_text(angle = 45)
dev.off()

pdf("Reviewer1/Heatmap_CellClusterPerDonor.pdf",width=6,height=6)
table(metadata$Cell,metadata$donor,metadata$Axis) %>%
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
    facet_wrap(~Var3) + theme(legend.position = "none") +
    rotate_x_text(angle = 45)
dev.off()





# Group
tmp1 <- metadata %>% 
			group_by(Axis, Cell,orig.ident) %>% 
			tally() %>%
			dplyr::rename(total = n)		 	

tmp2 <- metadata %>% 
			group_by(Axis, orig.ident) %>% 
			tally() %>%
			dplyr::rename(total_by_cell = n)		

tmp3 <- metadata %>%
        select(orig.ident, batch,age,sex,dur,fre,version) %>%
        group_by(orig.ident, batch,age,sex,dur,fre,version) %>%
        summarise()

Comb <- Reduce(dplyr::full_join, list(tmp1, tmp2,tmp3)) %>%
		mutate(Others = total_by_cell-total, Ratio = total/total_by_cell) %>% 
		as.data.frame() %>%
    mutate(batch = as.factor(batch), sex = as.factor(sex),version = as.factor(version))

# Model the cell abundance and gets statistics between Posterior/Anterior
formula = log2(total) ~ Axis*Cell + Axis*version + Axis*fre + (1|orig.ident)
model1 <- robustlmm::rlmer(formula = formula, method="DASvar", data = Comb)

emm1 <- emmeans(model1, specs = revpairwise ~ Axis | Cell)
res_model <- emm1$contrasts %>%
        summary(infer = TRUE, type = 'response') %>%
        rbind() %>%
        as.data.frame() %>%
        arrange(desc(estimate)) %>%
        mutate(FDR = p.adjust(p.value,method="BH")) %>%
        mutate(Sign = if_else(p.value < 0.05, "Sign","NotSign"))

openxlsx::write.xlsx(res_model, file = "Reviewer1/cell_proportion_MixModStats.xlsx", colNames = TRUE, borders = "columns")

# Visualize the odds ratios vs -log10(FDR)
pdf("Reviewer1/AnteriorVsPosterior_CellAbundance.pdf",width=5,height=5)
ggplot(aes(x = estimate, y = -log10(p.value)), data = res_model) +
  geom_point() + 
  #geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.value))) +
  geom_text_repel(aes(label = Cell), data = res_model %>% filter(p.value < 0.05)) + 
  #scale_x_log10() + 
  theme_minimal() + 
  labs(x = 'Anterior / Posterior', title = 'Cell type proportion change') + 
  geom_vline(xintercept = 0, linetype="dashed", color = "blue", size=1.5) +
  xlim(-15,15)
dev.off()

# Visualize the abundance by boxplot
pdf("Reviewer1/AnteriorVsPosterior_CellAbundance_Boxplot.pdf",width=8,height=8)
Comb$Ratio <- Comb$total/Comb$Others
Comb$log <- log10(Comb$total)
ggboxplot(Comb, x = "Axis", y = "log",
          color = "Axis", palette = c("gray60", "purple"),
          add = "jitter") +
          theme_minimal() +  
          #scale_y_log10() +
          facet_wrap(.~Cell) + 
          ylab("log10(Count)") +
  		theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# By Class
# Group
tmp1 <- metadata %>% 
      group_by(Axis, Definition,orig.ident) %>% 
      tally() %>%
      dplyr::rename(total = n)      

tmp2 <- metadata %>% 
      group_by(Axis, orig.ident) %>% 
      tally() %>%
      dplyr::rename(total_by_cell = n)    

tmp3 <- metadata %>%
        select(orig.ident, batch,age,sex,dur,fre,version) %>%
        group_by(orig.ident, batch,age,sex,dur,fre,version) %>%
        summarise()

Comb <- Reduce(dplyr::full_join, list(tmp1, tmp2,tmp3)) %>%
    mutate(Others = total_by_cell-total, Ratio = total/total_by_cell) %>% 
    as.data.frame() %>%
    mutate(batch = as.factor(batch), sex = as.factor(sex),version = as.factor(version))

# Model the cell abundance and gets statistics between Posterior/Anterior
formula = log2(total) ~ Axis*Definition + Axis*version + Axis*fre + (1|orig.ident)
model1 <- robustlmm::rlmer(formula = formula, method="DASvar", data = Comb)

emm1 <- emmeans(model1, specs = revpairwise ~ Axis | Definition)
res_model <- emm1$contrasts %>%
        summary(infer = TRUE, type = 'response') %>%
        rbind() %>%
        as.data.frame() %>%
        arrange(desc(estimate)) %>%
        mutate(FDR = p.adjust(p.value,method="BH")) %>%
        mutate(Sign = if_else(p.value < 0.05, "Sign","NotSign"))

openxlsx::write.xlsx(res_model, file = "Reviewer1/cell_class_proportion_MixModStats.xlsx", colNames = TRUE, borders = "columns")

# Visualize the odds ratios vs -log10(FDR)
pdf("Reviewer1/AnteriorVsPosterior_CellClass_Abundance.pdf",width=5,height=5)
ggplot(aes(x = estimate, y = -log10(p.value)), data = res_model) +
  geom_point() + 
  #geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.value))) +
  geom_text_repel(aes(label = Definition), data = res_model %>% filter(p.value < 0.05)) + 
  #scale_x_log10() + 
  theme_minimal() + 
  labs(x = 'Anterior / Posterior', title = 'Cell type proportion change') + 
  geom_vline(xintercept = 0, linetype="dashed", color = "blue", size=1.5) +
  xlim(-15,15)
dev.off()

# Visualize the abundance by boxplot
pdf("Reviewer1/AnteriorVsPosterior_CellClass_Abundance_Boxplot.pdf",width=8,height=8)
Comb$Ratio <- Comb$total/Comb$Others
Comb$log <- log10(Comb$total)
ggboxplot(Comb, x = "Axis", y = "log",
          color = "Axis", palette = c("gray60", "purple"),
          add = "jitter") +
          theme_minimal() +  
          #scale_y_log10() +
          facet_wrap(.~Definition) + 
          ylab("log10(Count)") +
      theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# No specific cells
tmp1 <- metadata %>% 
      group_by(Axis, Cell,orig.ident) %>% 
      tally() %>%
      dplyr::rename(total = n)      

tmp2 <- metadata %>% 
      group_by(Axis, orig.ident) %>% 
      tally() %>%
      dplyr::rename(total_by_cell = n)    

tmp3 <- metadata %>%
        select(orig.ident, batch,age,sex,dur,fre,version) %>%
        group_by(orig.ident, batch,age,sex,dur,fre,version) %>%
        summarise()

Comb <- Reduce(dplyr::full_join, list(tmp1, tmp2,tmp3)) %>%
    mutate(Others = total_by_cell-total, Ratio = total/total_by_cell) %>% 
    as.data.frame() %>%
    filter(!Cell %in% c("Olig5","Den.Gyr3","OPC2")) %>%
    mutate(batch = as.factor(batch), sex = as.factor(sex),version = as.factor(version))

# Model the cell abundance and gets statistics between Posterior/Anterior
formula = log2(total) ~ Axis*Cell + Axis*version + Axis*fre + (1|orig.ident)
model1 <- robustlmm::rlmer(formula = formula, method="DASvar", data = Comb)

emm1 <- emmeans(model1, specs = revpairwise ~ Axis | Cell)
res_model <- emm1$contrasts %>%
        summary(infer = TRUE, type = 'response') %>%
        rbind() %>%
        as.data.frame() %>%
        arrange(desc(estimate)) %>%
        mutate(FDR = p.adjust(p.value,method="BH")) %>%
        mutate(Sign = if_else(p.value < 0.05, "Sign","NotSign"))

openxlsx::write.xlsx(res_model, file = "Reviewer1/cell_proportion_MixModStats_ExceptSomeCells.xlsx", colNames = TRUE, borders = "columns")

# Visualize the odds ratios vs -log10(FDR)
pdf("Reviewer1/AnteriorVsPosterior_CellAbundance_ExceptSomeCells.pdf",width=5,height=5)
ggplot(aes(x = estimate, y = -log10(p.value)), data = res_model) +
  geom_point() + 
  #geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.value))) +
  geom_text_repel(aes(label = Cell), data = res_model %>% filter(p.value < 0.05)) + 
  #scale_x_log10() + 
  theme_minimal() + 
  labs(x = 'Anterior / Posterior', title = 'Cell type proportion change') + 
  geom_vline(xintercept = 0, linetype="dashed", color = "blue", size=1.5) +
  xlim(-15,15)
dev.off()

# Visualize the abundance by boxplot
pdf("Reviewer1/AnteriorVsPosterior_CellAbundance_Boxplot_ExceptSomeCells.pdf",width=8,height=8)
Comb$Ratio <- Comb$total/Comb$Others
Comb$log <- log10(Comb$total)
ggboxplot(Comb, x = "Axis", y = "log",
          color = "Axis", palette = c("gray60", "purple"),
          add = "jitter") +
          theme_minimal() +  
          #scale_y_log10() +
          facet_wrap(.~Cell) + 
          ylab("log10(Count)") +
      theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

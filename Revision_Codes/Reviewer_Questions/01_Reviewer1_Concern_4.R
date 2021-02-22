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
})

dir.create("Reviewer1")

sce <- readRDS("input_data/HippoAxis_Filt.rds")

################### 
# Add information #
###################
metadata <- colData(sce) %>%
				as.data.frame()


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
        select(orig.ident, batch,age,sex,dur) %>%
        group_by(orig.ident, batch,age,sex,dur) %>%
        summarise()

Comb <- Reduce(dplyr::full_join, list(tmp1, tmp2,tmp3)) %>%
		mutate(Others = total_by_cell-total) %>% 
		as.data.frame() %>%
    mutate(batch = as.factor(batch), sex = as.factor(sex))

# Model the cell abundance and gets statistics between Posterior/Anterior
formula = cbind(total, Others) ~ Cell * Axis + batch*Axis + sex*Axis + (1|orig.ident)
model1 <- lme4::glmer(formula = formula, family = 'binomial', data = Comb)

emm1 <- emmeans(model1, specs = revpairwise ~ Axis | Cell)
res_model <- emm1$contrasts %>%
  			summary(infer = TRUE, type = 'response') %>%
  			rbind() %>%
  			as.data.frame() %>%
  			arrange(desc(odds.ratio)) %>%
  			mutate(FDR = p.adjust(p.value,method="bonferroni")) %>%
        mutate(Sign = if_else(FDR < 0.05, "Sign","NotSign"))

openxlsx::write.xlsx(res_model, file = "Reviewer1/cell_proportion_statistics.xlsx", colNames = TRUE, borders = "columns")

# Visualize the odds ratios vs -log10(FDR)
pdf("Reviewer1/AnteriorVsPosterior_CellAbundance.pdf",width=5,height=5)
ggplot(aes(x = odds.ratio, y = -log10(FDR)), data = res_model) +
  geom_point() + 
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(FDR))) +
  geom_text_repel(aes(label = Cell), data = res_model %>% filter(FDR < 0.01)) + 
  scale_x_log10() + 
  theme_minimal() + 
  labs(x = 'Anterior / Posterior', title = 'Cell type proportion change') + 
  geom_vline(xintercept = 1, linetype="dashed", color = "blue", size=1.5)
dev.off()

# Visualize the abundance by boxplot
pdf("Reviewer1/AnteriorVsPosterior_CellAbundance_Boxplot.pdf",width=8,height=8)
Comb$Ratio <- Comb$total/Comb$Others
Comb$log <- log10(Comb$Ratio)
ggboxplot(Comb, x = "Axis", y = "Ratio",
          color = "Axis", palette = "jco",
          add = "jitter") +
          theme_minimal() +  
          scale_y_log10() +
          facet_wrap(.~Cell) + 
          ylab("log10(Count)") +
  		theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# now with cell propeller (speckle)

prop <- propeller(clusters = metadata$Cell, sample = metadata$orig.ident, group = metadata$Axis,robust = FALSE)
openxlsx::write.xlsx(prop, file = "Reviewer1/cell_proportion_statistics_speckle.xlsx", colNames = TRUE, borders = "columns")





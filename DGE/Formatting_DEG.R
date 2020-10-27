library(tidyverse)

files = list.files(pattern = '*.txt')
names <- gsub( "LMM_|_filtered_0.3fc_0.25percent.txt", "", files )
GeneSets = lapply(files, read.table,header=T,sep="\t")
names(GeneSets) <- names


df <- GeneSets %>% 
		map2_df(names,~mutate(.x,Class=.y)) %>%
		mutate(Direction= case_when(FC > 0 ~ "Posterior", FC < 0  ~ "Anterior")) %>%
		unite("ID", Class:Direction, remove = FALSE)

write.table(df,"HippAxisSeq_Degs.txt",sep="\t",quote=F)
openxlsx::write.xlsx(df, file = "HippAxisSeq_Degs.xlsx", colNames = TRUE, borders = "columns")


df <- GeneSets %>% 
		map2_df(names,~mutate(.x,Class=.y)) %>%
		mutate(Direction= case_when(FC > 0 ~ "Posterior", FC < 0  ~ "Anterior")) %>%
		unite("ID", Class:Direction, remove = FALSE) %>%
		select(Gene,ID)

write.table(df,"HippAxisSeq_Degs_MagmaFormat.txt",sep="\t",quote=F,row.names=F,col.names=F)
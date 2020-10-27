files <- list.files(pattern = ".gsa.out", full.names = FALSE)

ReadIn_ReadDistrbution <- function(FileName){
   in.df <- read.table(FileName, skip = 3, stringsAsFactors = FALSE, fill=TRUE, header=T)
   in.df <- in.df[complete.cases(in.df),]
   in.df <- transform(in.df,
                      Sample = gsub(".gsa.out","\\1", FileName)
   )
 return(in.df)
}

l <- lapply(files, function(x) ReadIn_ReadDistrbution(x))
data <- as.data.frame(do.call(rbind, l))
write.table(data,"MAGMA_STATISTICS.txt",sep="\t",quote=F)

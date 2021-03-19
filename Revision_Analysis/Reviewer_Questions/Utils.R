##########################################################
## Functions correlation and brain waves         	##
## FastCor = multicore rapid correlation		##
## VarExp = variance explained by PCAs 			##
## plotVarExp = Plot variance explained  		##
##########################################################


#### Rapid correlation
# exp = gene x sample matrix
# sme = sample x oscillation (one by one as loop)
FastCor <- function(exp,sme,method="pearson",alternative="two.sided",cores=1,override=FALSE,...){
	suppressPackageStartupMessages(library(parallel))
	mat = as.matrix(exp)
	nr = nrow(mat)
	if(nr > 1000){
		if(override){
			show("Wait for it!")
		} else {
			stop("Set override to TRUE")
		}
	}
	if(is.null(rownames(mat))){
		nm = as.character(1:nr)
	} else { 
		nm = rownames(mat)
	}
	corrAndP = mclapply(1:nr,function(i){
			res = cor.test(mat[i,],as.numeric(sme),method=method,alternative=alternative,exact=FALSE)
			cbind(res$estimate,res$p.value)
		},mc.cores=cores,...)
	Rho = unlist(sapply(corrAndP,function(i){i[,1]}))
	Pval  = unlist(sapply(corrAndP,function(i){i[,2]}))
	results = cbind(Rho,Pval)
	rownames(results) <- rownames(exp)
	return(results)
}

#### Variance Explained
# counts = gene x sample matrix
# meta = sample x covariate matrix
# threshold = number of PCA to consider (e.g. 5)
# inter = interaction between covariates (e.g. Age:Sex)

VarExp <- function(counts, meta, threshold, inter){
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(optimx))
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)

  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}

  pred.list <- colnames(meta)
  meta <- droplevels(meta)

  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}

  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit,control = lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}

# Bar plot for data visualization
plotVarExp <- function(pvca.res, title){
  suppressPackageStartupMessages(library(ggplot2))
  plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
  p <- ggplot2::ggplot(plot.dat, aes(x=eff, y=prop))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::geom_bar(stat="identity", fill="steelblue", colour="steelblue")
  p <- p + ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4)
  p <- p + ggplot2::scale_x_discrete(limits=names(pvca.res))
  p <- p + ggplot2::scale_y_continuous(limits = c(0,1))
  p <- p + ggplot2::labs(x= "Effects", y= "WAPV")
  p <- p + ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}

# Fisher Exact multi-factor
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

### For Single Cell

##-------------------------------------------------------
## FUNCTION TO SEURATIFY DATA
##-------------------------------------------------------
seuratify <- function(countData, sampleName)
 {
 seuratObject <- CreateSeuratObject(counts = countData, project = sampleName)
 # save(seuratObject, countData, file = paste(sampleName, "_DATA.RData", sep = ""))
 return(seuratObject)
 }

##-------------------------------------------------------
## FUNCTION TO CREATE BIN HISTOGRAM
##-------------------------------------------------------
plotbinhist <- function(seuratObject)
 {
 brks <- unique(c(seq(0, 1000, by = 100), seq(1000, 25000, by = 1000)))
 # brks <- c(0, 200, 400, 600, 800, 1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 100000)
 umiBins <- cut(seuratObject@meta.data$nCount_RNA, brks, include.lowest = T, right = FALSE)

 binData <- cbind(summary(umiBins))
 colnames(binData) <- c("nUMI")
 binData <- as.data.frame(binData)
 binData$bins <- row.names(binData)
 binData2 <- melt(binData)
 binData2$bins <- factor(binData2$bins, levels = row.names(binData))

 plotBins <- ggplot(binData2, aes(bins, value, fill = variable)) +
             geom_bar(stat = "identity", position = "dodge") +
             # facet_grid(rows = vars(variable)) +
             facet_grid(. ~ variable) +
             # coord_flip() +
             theme_bw() +
             theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
             geom_text_repel(aes(label = value), vjust = -1, angle = 90) +
             ylim(0, 50000) +
             labs(title = seuratObject@project.name, x = "Bins", y = "Number of Cells") +
             NULL
 # plotBins
 ggsave(paste(seuratObject@project.name, "_UMI_HIST.pdf", sep = ""), plot = plotBins, width = 16, height = 8, units = "in", dpi = 300)  
 }

##-------------------------------------------------------
## MITO GENES AND QC PLOTS
##-------------------------------------------------------
plotqcmito <- function(seuratObject)
 {
 ## Identify mitochondrial genes and plot percent mito genes along with nGenes and nUMI
 mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObject@assays$RNA@data), value = TRUE)
 percent.mito <- Matrix::colSums(seuratObject@assays$RNA@data[mito.genes, ])/Matrix::colSums(seuratObject@assays$RNA@data)
 seuratObject <- AddMetaData(object = seuratObject, metadata = percent.mito, col.name = "pMito")
 
 ## QC Plots before removing mito genes
 plot_nU <- VlnPlot(object = seuratObject, features = "nCount_RNA", pt.size = 0)
 plot_nG <- VlnPlot(object = seuratObject, features = "nFeature_RNA", pt.size = 0)
 plot_pM <- VlnPlot(object = seuratObject, features = "pMito", pt.size = 0)
 plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
 ggsave(paste(seuratObject@project.name, "_QC_1.pdf", sep = ""), plot = plotQC1, width = 12, height = 4, units = "in", dpi = 300)

 plot_nUpM <- FeatureScatter(object = seuratObject, feature1 = "nCount_RNA", feature2 = "pMito", cex.use = 1)
 plot_nUnG <- FeatureScatter(object = seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cex.use = 1)
 plotQC2 <- grid.arrange(plot_nUpM, plot_nUnG, ncol = 2)
 ggsave(paste(seuratObject@project.name, "_QC_2.pdf", sep = ""), plot = plotQC2, width = 10, height = 4, units = "in", dpi = 300)

 return(seuratObject)
 }

##-------------------------------------------------------
## FUNCTION TO CREATE HISTOGRAM
##-------------------------------------------------------
plotnumihist <- function(all.data, meta.data, prefix)
 {
 ## Generate histograms for UMI data
 dataHist <- as.data.frame(colSums(all.data))
 colnames(dataHist) <- "UMI_Count"

 hist.data <- merge(dataHist, meta.data, by = "row.names")
 row.names(hist.data) <- hist.data$Row.names
 hist.data$Row.names <- NULL

 ## Plot histograms UMI Distribution (auto-scale)
 plotHist <- ggplot(hist.data, aes(x = log10(UMI_Count), fill = LIBRARY)) + geom_histogram(bins = 100) + facet_grid(. ~ LIBRARY) + labs(title = prefix, x = "Log10(Number of UMIs)", y = "Number of Cells") + theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 90))
 ggsave(paste(prefix, "_HIST.pdf", sep = ""), plot = plotHist, width = 12, height = 4, units = "in", dpi = 150)
 }


## Stacked violin
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) +
        rotate_x_text(angle = 45)
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y)) 

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


SCTransform_ChooseRes <- function(
    input.srobj, n.pcs, sample.name =  format(Sys.time(), "%a-%b-%d-%X-%Y-%Z"),
    res.low = .01, res.high=3, res.n = 40, bias = "over", figdir=F) {

    ######## step 1: save the input seurat object as a new temporary object, 
    ########         dont want to overwrite or change the original one with all of the parameter scans
    
    srobj.tmp <- input.srobj 
    # in case there have been other things calculated in the metadata, just cut down to simplify/avoid errors
    srobj.tmp@meta.data = srobj.tmp@meta.data[,c(2:3)] # should just be the nUMI and nGene 
    

    ######## step 2: calculate the FindClusters over a large range of resolutions
    print("Performing parameter scan over multiple resolutions...")

    set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=3)

    DefaultAssay(srobj.tmp) <- "integrated"

    srobj.tmp = FindClusters(srobj.tmp, 
                               resolution = set.res,
                               algorithm = 1,
                               n.iter = 1000)

    #for(i in 2:length(set.res)){
    #    srobj.tmp = FindClusters(
    #        srobj.tmp,
    #        resolution = set.res[i], 
    #        algorithm = 1,
    #        n.iter = 1000)
    #    print(paste("          ", round(100*i/length(set.res)), "% done with parameter scan", sep=""))
    #}


    ######## step 3: output plot of how the resolution changes the number of clusters you get
    n.clusters = vector(mode="numeric", length=length(set.res))
    names(n.clusters) = set.res
    for(i in 1:length(n.clusters)){
        n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[,paste("integrated_snn_res.", names(n.clusters)[i], sep="")])))
    }
    
    ######## step 4: calculate the silhouette width for each resolution
    print("Computing a silhouette width for each cell, for each resolution...")
    require(cluster)

    dist.temp = cor(t(srobj.tmp@reductions$pca@cell.embeddings[,n.pcs]), method="pearson")
    dist.temp[is.na(dist.temp)] <- 0
    random.cells.choose = sample(1:nrow(dist.temp), round(nrow(dist.temp)/10, digits=0))
    dist.temp.downsample = dist.temp[random.cells.choose, random.cells.choose]
    sil.all.matrix = matrix(data=NA, nrow=nrow(dist.temp.downsample), ncol=0)
    
    for(i in 1:length(set.res)){
        clusters.temp = as.numeric(as.vector(
            srobj.tmp@meta.data[random.cells.choose,paste("integrated_snn_res.", set.res[i], sep="")]))
        if(length(table(clusters.temp))>1){
            sil.out = silhouette(clusters.temp, as.dist(1-as.matrix(dist.temp.downsample)))
            sil.all.matrix = cbind(sil.all.matrix, sil.out[,3])
        }
        if(length(table(clusters.temp))==1){
            sil.all.matrix = cbind(sil.all.matrix, rep(0, length(clusters.temp)))
        }
        print(paste("          ", round(100*i/length(set.res)), "% done with silhouette calculation", sep=""))
    
    }


    ######## step 5: calculate summary metric to compare the silhouette distributions,
    ########         average has worked well so far... could get fancier

    print("Identifying a best resolution to maximize silhouette width")
    sil.average = colMeans(sil.all.matrix)
    names(sil.average) = set.res


    ######## step 6: automate choosing resolution that maximizes the silhouette 
    hist.out = hist(sil.average, length(sil.average)/1.2,  plot=FALSE)
    
    #  take the ones that fall into the top bin, 
    #  and the max OR MIN of those  ******* can change this to under vs over cluster
    if(bias=="over"){
        resolution.choice = as.numeric(max(
            names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
    }
    if(bias=="under"){
        resolution.choice = as.numeric(min(
            names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
    }
    
    # get the silhouette of the best resolution: 
    silhouette.best = as.numeric(sil.average[paste(resolution.choice)])

    print(paste("Best Resolution Choice: ", resolution.choice, ", with average silhouette score of: ",
            round(silhouette.best, digits=3), ", giving ", as.numeric(n.clusters[paste(resolution.choice)]),
            " clusters", sep=""))


    ######### step 7: output plot and data 
    
    if (figdir) {
        
        print(paste0("Ouptutting summary statistics and returning seurat object... ",
              "This will create a pdf in your output directory,",
              " and will return your input seurat object ammended with the best choice",
              " for clusters (found as Best.Clusters in the meta.data matrix, and set to your new ident)..."))

        pdf(paste(figdir, "/", sample.name, ".pdf", sep=""),
        width=10, height=4, useDingbats=FALSE)
        par(mfrow=c(1,3))
        # Resolution vs # of Clusters
        plot(set.res, n.clusters, col="black", pch=19,
             type="p", xlab="Resolution", ylab="# Clusters",
             main="Resolution vs. # Clusters")
        # Resolution vs Average Silhouette
        plot(set.res, sil.average, col="black", pch=19,
             type="p", xlab="Resolution", ylab="Average Silhouette",
             main="Resolution vs. Average Silhouette")
        abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
        abline(v=resolution.choice, col="dodgerblue2", lty=2)

        # N Clusters vs Average Silhouette
        plot(n.clusters, sil.average, col="black", pch=19,
             type="p", xlab="# Clusters", ylab="Average Silhouette",
             main="# Clusters vs. Average Silhouette")
        abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
        abline(v=as.numeric(n.clusters[paste(resolution.choice)]), col="dodgerblue2", lty=2)
        dev.off()
    }
    
    ######## step 8: return the original seurat object, with the metadata containing a 
    ########         concatenated vector with the clusters defined by the best choice here,
    ########         as well as the ident set to this new vector
    
    Best.Clusters = srobj.tmp@meta.data[,paste("integrated_snn_res.", resolution.choice, sep="")]
    
    input.srobj$Best.Clusters = Best.Clusters
    Idents(object = input.srobj) = input.srobj$Best.Clusters
    input.srobj@misc$resolution.choice <- resolution.choice
    return(input.srobj)
}




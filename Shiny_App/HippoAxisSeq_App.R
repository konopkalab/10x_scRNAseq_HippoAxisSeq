library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("tidyverse")
library("BiocParallel")
library("scater")


# Download the data from here: "https://www.dropbox.com/s/6j0wroa41hnr5m9/HippoAxis.rds?dl=1"
# The data is ~700mb. 
sce <- readRDS("HippoAxis.rds")
sce <- scater::addPerCellQC(sce)
sce <- scater::addPerFeatureQC(sce,BPPARAM = SerialParam())

################### 
# Add information #
###################
temporary <- colData(sce) %>%
				as.data.frame() %>%
                rownames_to_column("TMP") %>%
                #select(TMP,orig.ident,nCount_RNA,nFeature_RNA, percent.mt, group, age, sex,dur,batch, donor,ident) %>%
                mutate(ident = as.character(ident)) %>%
                mutate(ident = if_else(ident == "29", "CApyr2",ident)) %>%
                mutate(Cell = as.factor(ident)) %>%
                mutate(Class = case_when(grepl("CA|Den", Cell) ~ "Glutamatergic", 
                        grepl("In", Cell) ~ "Gabaergic",
                        grepl("Astro|Olig|OPC|Endo|Micro", Cell) ~ "NonNeuronal")) %>%
                mutate(Definition = case_when(grepl("CA", Cell) ~ "CA.Neurons",
                		grepl("Den", Cell) ~ "GranuleN",
                        grepl("In", Cell) ~ "Inhibitory",
                        grepl("Astro", Cell) ~ "Astrocytes", 
                        grepl("Olig", Cell) ~ "Oligodendrocytes",
                        grepl("OPC", Cell) ~ "OPC",
                        grepl("Micro", Cell) ~ "Microglia",
                        grepl("Endo", Cell) ~ "Endothelial")) %>%
                rename(Axis = group) %>%
                column_to_rownames("TMP") %>% 
                DataFrame()

colData(sce) <- temporary

# Add info about % of expressed per cell
splitit <- function(x) split(seq(along=x),x)
cellType.idx <- splitit(colData(sce)[["Cell"]])
    rowdat.sce <- rowData(sce)
    for(i in names(cellType.idx)){
        message(Sys.time(), " computing propNucleiExprs for ", i)
        rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
            assay(sce, "counts")[, cellType.idx[[i]]],
            1,
            function(x){
                mean(x != 0)
            }
        )
    }

rowData(sce) <- rowdat.sce

# Initialize the iSEE app.

initial <- list()

###############sce#################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "UMAP", XAxis = 1L, YAxis = 2L,
    ColorByColumnData = "Definition", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Axis",
    SizeByColumnData = "sum", FacetByRow = "---", FacetByColumn = "---",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "CADPS2", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "---",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionEffect = "Transparent",
    SelectionColor = "#FF0000", SelectionAlpha = 0.1, ZoomData = numeric(0),
    BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = c("Color",
        "Shape"), ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
    CustomRowsText = "LINGO2
CADPS2
GRIA4
SHISA9
WIPF3
EPHA7", 
	ClusterRows = FALSE,
    ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE, VisualChoices = "Annotations",
    ColumnData = c("Definition", "Class","Axis"),
    RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
    UpperBound = NA_real_, AssayCenterRows = FALSE, AssayScaleRows = FALSE,
    DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
    LegendPosition = "Bottom", LegendDirection = "Horizontal",
    VisualBoxOpen = FALSE, SelectionEffect = "Color", SelectionColor = "#FF0000",
    PanelId = 1L, PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "Definition", XAxisFeatureName = "CADPS2",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "CADPS2", YAxisFeatureSource = "RowDataTable",
    YAxisFeatureDynamicSource = TRUE, ColorByColumnData = "Axis",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "Axis", SizeByColumnData = "sum", FacetByRow = "---",
    FacetByColumn = "Axis", ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "CADPS2", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "P57_AAAGTAGGTCCAGTAT",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionEffect = "Transparent",
    SelectionColor = "#FF0000", SelectionAlpha = 0.1, ZoomData = numeric(0),
    BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = "Color",
    ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())


################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Search = "", SearchColumns = c("",
    "", "", "", "", "", "", "", "", "", "", "", "", "", ""), PanelId = c(RowDataTable = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Search = "", SearchColumns = c("",
    "", "", "", "", "", "", "", "", "", "", "", "", "", ""), PanelId = c(RowDataTable = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

shiny::runApp(iSEE(sce, appTitle = "Ahyan et al. 2020, HippoAxis scRNA-seq", initial = initial))
Hippocampus Axis data analysis
==========================

This repository contains analysis code for the Hippocampus Axis scRNA-seq project carried out by researchers at the [Konopka Lab, UTSW](http://konopkalab.org/).

## Cite this

If you use anything in this repository please cite the following publication:

Paper: https://www.cell.com/neuron/fulltext/S0896-6273(21)00329-9

## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`processing_qc`](processing_qc/) | Output data from initial processing and quality check. | 01_Data_processing_QC.R |
| [`DGE`](DGE/) | Output data from DGE analysis. | 02_DGE_Analysis.R |
| [`MAGMA`](MAGMA/) | Output data from MAGMA gene set analysis. | 03_MAGMA.sh |
| [`Shiny_App`](Shiny_App/) | Input data to visualize the data. | HippoAxisSeq_App.R |

## Explore the data

We have provided an interactive web app that allow you to explore the data at single nucleus level. 

* https://human-hippo-axis.cells.ucsc.edu/

* Shiny app (Please download the data and build it using the HippoAxisSeq_App.R code).

![](HippoAxisSeq.gif)]


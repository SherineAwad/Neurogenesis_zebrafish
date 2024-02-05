library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)


args <- commandArgs(trailingOnly = TRUE)

mysample = args[1]
myRDS <- paste(mysample, "_subsetClusters.rds", sep="")

mysample
myRDS

myObject <- readRDS(myRDS)


DefaultAssay(myObject) <- "RNA"

myObject <- SCTransform(myObject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'UMAP_') 

myObject <- FindNeighbors(myObject, dims = 1:10)
myObject <- FindClusters(myObject, resolution = 1.2)
myObject <- RunUMAP(object = myObject, dims = 1:10) 

figure_name <- "" 
figure_name <- paste(mysample, "SubsetClusters.pdf", sep="") 
pdf(file =figure_name, width =12) 
DimPlot(myObject, reduction = "umap")
dev.off()

myRDS <- paste(mysample, "_cleaned.rds", sep="")
saveRDS(myObject, file = myRDS)






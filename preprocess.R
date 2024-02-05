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
myRDS <- paste(mysample, ".rds", sep="")

mysample
myRDS

myObject <-readRDS(myRDS) 

#idents <- Idents(myObject) <- factor(x = Idents(myObject), levels = sort(levels(myObject)))
#Idents(myObject) <- "sample"
#idents = Idents(myObject) 
#pdf(file="labelclusters.pdf", width =20)
#plot <- DimPlot(object = myObject)
#LabelClusters(plot = plot, id = 'ident')
#dev.off()


#Idents(myObject) <- "sample"
#cell_values <- c("LD_36hr","LD_54hr","LD_72hr","LD_96hr","LD_7D","LD_14D")
#mySubset <- subset(myObject, idents = cell_values, invert = FALSE)
#pdf(file="subsetclusters.pdf", width =20)
#plot <- DimPlot(object = myObject)
#LabelClusters(plot = plot, id = 'ident')
#dev.off()


DefaultAssay(myObject) <- "RNA"

myObject <- FindNeighbors(myObject, dims = 1:10)
myObject <- FindClusters(myObject, resolution = 1.6) 

figure_name <- ""
figure_name <- paste(mysample, "zebra_clusters.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


myRDS <- paste(mysample, "_Clusters.rds", sep="")
saveRDS(myObject, file = myRDS)



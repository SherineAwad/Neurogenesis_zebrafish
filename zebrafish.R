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

myObject$sample 
#myObject$new_celltype 
#Name of clusters
#myObject$tag1
#works well
#idents <- Idents(object = myObject)

print("new_celltype")
myObject$new_celltype
idents <- Idents(myObject) <- factor(x = Idents(myObject), levels = sort(levels(myObject)))
head(idents)  
#works well 
#mySubset <- subset(myObject, idents = c(64,278)) 

#AC is one of those new_celltype 
#mySubset <- subset(myObject, subset = new_celltype == 'AC') 
#mySubset <- subset(myObject, subset = sample == "LD_36hr")

Idents(myObject) <- "sample"
cell_values <- c("LD_36hr","LD_54hr","LD_72hr","LD_96hr","LD_7D","LD_14D")
mySubset <- subset(myObject, idents = cell_values, invert = FALSE)
#head(mySubset) 

pdf(file="labelclusters.pdf", width =20)
plot <- DimPlot(object = myObject)
LabelClusters(plot = plot, id = 'ident')
dev.off()

pdf(file="subsetlabelclusters.pdf", width =20)
plot <- DimPlot(object = mySubset)
LabelClusters(plot = plot, id = 'ident')
dev.off()



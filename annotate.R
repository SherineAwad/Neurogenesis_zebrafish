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
myRDS <- paste(mysample, "_Clusters.rds", sep="")

mysample
myRDS

myObject <-readRDS(myRDS) 

myObject <- RenameIdents(
   object = myObject,
   "0" = 'BC (M)',  
   "1" = 'BC (M)', 
   "2" = 'Rod(M)',
   "3" = 'Rod(M)',
   "4" = 'Rod(M)',
   "5" = 'Rod(M)',
   "6" = 'HC', 
   "7" = 'Cone (M)', 
   "8" = 'AC(M)',
   "9" = 'BC (M)',
   "10" = 'Cone (M)',
   "11" = 'AC pre', 
   "12" = 'AC(M)',
   "13" = 'BC (M)',
   "14" = 'BC (M)',
   "15" = 'RGC(M)', 
   "16" = 'AC(M)',
   "17" = 'BC (M)',
   "18" = 'REST MG', 
   "19" = 'BC pre',
   "20" = 'Cone pre', 
   "21" = 'AC(M)', 
   "22" = 'BC (M)',
   "23" = 'Rod pre',
   "24" = 'AC(M)',
   "25" = 'Cone (M)',
   "26" = 'AC(M)', 
   "27" = 'Act MG',  
   "28" = 'Cone (M)',
   "29" = 'Rod(M)',
   "30" = 'MGPCs', 
   "31" = 'BC (M)', 
   "32" = 'MGPCs',
   "33" = 'BC (M)',
   "34" = 'BC (M)',
   "35" = 'AC(M)',
   "36" = 'AC(M)',
   "37" = 'HC',  
   "38" = 'AC(M)',
   "39" = 'RGC pre',  
   "40" = 'Microglia/Macrophages',
   "41" = 'AC(M)',
   "42" = 'Microglia/Macrophages',
   "43" = 'AC(M)')


figure_name <- ""
figure_name <- paste(mysample, "zebra_clusters_renamed.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap", group.by = "orig.ident",  repel = TRUE) + ggtitle("UMAP")
DimPlot(myObject, reduction = "umap", label=TRUE, repel = TRUE) + ggtitle("UMAP")
dev.off()


cell_values <- c("REST MG", "MGPCs", "AC pre", "BC pre", "Rod pre", "Cone pre", "RGC pre", "Act MG") 
mySubset <- subset(myObject, idents = cell_values, invert = FALSE)
head(mySubset)


figure_name <- ""
figure_name <- paste(mysample, "SubsetUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(mySubset, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
dev.off()




myRDS <- paste(mysample, "_subsetClusters.rds", sep="")
saveRDS(mySubset, file = myRDS)





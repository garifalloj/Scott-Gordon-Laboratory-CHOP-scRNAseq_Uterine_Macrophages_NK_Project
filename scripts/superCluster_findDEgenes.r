#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
#options(Seurat.object.assay.version = 'v3')

setwd(".")

merged_EndoFT_combinedSCT_res9 <- readRDS("merged_Endo_FT.combined.sct.rds")

#setIdent to the annotations if not active ident

###Find DEGs across Endo vs FT superclusters
###########set super cluster names#######
#find super cluster markers

super_cluster_names <- c(
    "Macrophage",
"ILC",
"T cell",
"ILC",
"T cell",
"Macrophage",
"Macrophage",
"ILC",
"ILC",
"DC",
"ILC",
"DC",
"ILC",
"Monocyte",
"Monocyte",
"B cell",
"Macrophage",
"Monocyte",
"T cell",
"ILC",
"ILC",
"T cell",
"ILC",
"ILC",
"Gran-Fibroblast",
"T cell",
"DC",
"Gran-Fibroblast",
"DC",
"B cell",
)

names(super_cluster_names) <- levels(merged_EndoFT_combinedSCT_res9)
merged_EndoFT_combinedSCT_res9 <- RenameIdents(merged_EndoFT_combinedSCT_res9, super_cluster_names)

#add idents back to object here
merged_EndoFT_combinedSCT_res9$supercluster <- Idents(merged_EndoFT_combinedSCT_res9)

#find DEGs for major super cluster markers
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
merged_EndoFT_combinedSCT_res9_super_cluster_markers <- FindAllMarkers(merged_EndoFT_combinedSCT_res9, only.pos = TRUE)



###add meta for superCluster_condition
merged_EndoFT_combinedSCT_res9$superCluster_condition <- paste(merged_EndoFT_combinedSCT_res9$supercluster, merged_EndoFT_combinedSCT_res9$condition, sep = "_")
Idents(merged_EndoFT_combinedSCT_res9) <- "superCluster_condition"

###next run comparisons
#should be 9
macrophage_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Macrophage_Endometrium", ident.2 = "Macrophage_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(macrophage_superClust_Endo_vs_FT_markers, file = "macrophage_superClust_Endo_vs_FT_markers.csv")

#DC
DC_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "DC_Endometrium", ident.2 = "DC_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DC_superClust_Endo_vs_FT_markers, file = "DC_superClust_Endo_vs_FT_markers.csv")

#Monocyte
Monocyte_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Monocyte_Endometrium", ident.2 = "Monocyte_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Monocyte_superClust_Endo_vs_FT_markers, file = "Monocyte_superClust_Endo_vs_FT_markers.csv")

#Tcell
Tcell_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Tcell_Endometrium", ident.2 = "Tcell_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tcell_superClust_Endo_vs_FT_markers, file = "Tcell_superClust_Endo_vs_FT_markers.csv")


#ILC
ILC_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "ILC_Endometrium", ident.2 = "ILC_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ILC_superClust_Endo_vs_FT_markers, file = "ILC_superClust_Endo_vs_FT_markers.csv")

#Bcell
Bcell_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Bcell_Endometrium", ident.2 = "Bcell_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Bcell_superClust_Endo_vs_FT_markers, file = "Bcell_superClust_Endo_vs_FT_markers.csv")



#Gran-Fibroblast
Gran-Fibroblast_superClust_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Gran-Fibroblast_Endometrium", ident.2 = "Gran-Fibroblast_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Gran-Fibroblast_superClust_Endo_vs_FT_markers, file = "Gran-Fibroblast_superClust_Endo_vs_FT_markers.csv")



###need to subset all supercluster marker gene lists for pvalue < 0.05 and logfc > .5 and logfc < -1
perSuperCluster_Endo_vs_FT_markers <- list.files(pattern = '*_superClust_Endo_vs_FT_markers.csv')
for (i in perSuperCluster_Endo_vs_FT_markers){
    markers <- read.csv(i, row.names = 1)
    #pvalue_adjusted
    markers_subset <- subset(markers, p_val_adj < 0.05)
    #log2FC > 0.5 OR log2FC < -1
    markers_subset <- subset(markers_subset, avg_log2FC > .5 | avg_log2FC < -1)
    assign(i, markers_subset)
}

#then write all to new csv files
write.csv(macrophage_superClust_Endo_vs_FT_markers.csv, file ="subset_pval_log2fc_macrophage_superClust_Endo_vs_FT_markers.csv")
write.csv(DC_superClust_Endo_vs_FT_markers.csv, file ="subset_pval_log2fc_DC_superClust_Endo_vs_FT_markers.csv")
write.csv(Monocyte_superClust_Endo_vs_FT_markers.csv, file ="subset_pval_log2fc_Monocyte_superClust_Endo_vs_FT_markers.csv")
write.csv(Tcell_superClust_Endo_vs_FT_markers.csv, file ="subset_pval_log2fc_Tcell_superClust_Endo_vs_FT_markers.csv")
write.csv(ILC_superClust_Endo_vs_FT_markers.csv, file ="subset_pval_log2fc_ILC_superClust_Endo_vs_FT_markers.csv")
write.csv(Bcell_superClust_Endo_vs_FT_markers.csv, file ="subset_pval_log2fc_Bcell_superClust_Endo_vs_FT_markers.csv" )
write.csv(Gran-Fibroblast_superClust_Endo_vs_FT_markers, file ="subset_pval_log2fc_Gran-Fibroblast_superClust_Endo_vs_FT_markers.csv")





#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
#options(Seurat.object.assay.version = 'v3')

setwd(".")

merged_EndoFT_combinedSCT_res9 <- readRDS("merged_Endo_FT.combined.sct.rds")



#####generate per seurat cluster conserved marker genes
cluster0_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 0, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster0_conserved_markers_thresholds_res9, file = "cluster0_conserved_markers_thresholds_res9.csv")

cluster1_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 1, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster1_conserved_markers_thresholds_res9, file = "cluster1_conserved_markers_thresholds_res9.csv")


cluster2_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 2, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster2_conserved_markers_thresholds_res9, file = "cluster2_conserved_markers_thresholds_res9.csv")


cluster3_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 3, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster3_conserved_markers_thresholds_res9, file = "cluster3_conserved_markers_thresholds_res9.csv")



cluster4_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 4, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster4_conserved_markers_thresholds_res9, file = "cluster4_conserved_markers_thresholds_res9.csv")


cluster5_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 5, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster5_conserved_markers_thresholds_res9, file = "cluster5_conserved_markers_thresholds_res9.csv")


cluster6_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 6, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster6_conserved_markers_thresholds_res9, file = "cluster6_conserved_markers_thresholds_res9.csv")


cluster7_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 7, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster7_conserved_markers_thresholds_res9, file = "cluster7_conserved_markers_thresholds_res9.csv")


cluster8_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 8, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster8_conserved_markers_thresholds_res9, file = "cluster8_conserved_markers_thresholds_res9.csv")


cluster9_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 9, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster9_conserved_markers_thresholds_res9, file = "cluster9_conserved_markers_thresholds_res9.csv")


cluster10_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 10, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster10_conserved_markers_thresholds_res9, file = "cluster10_conserved_markers_thresholds_res9.csv")


cluster11_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 11, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster11_conserved_markers_thresholds_res9, file = "cluster11_conserved_markers_thresholds_res9.csv")


cluster12_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 12, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster12_conserved_markers_thresholds_res9, file = "cluster12_conserved_markers_thresholds_res9.csv")


cluster13_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 13, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster13_conserved_markers_thresholds_res9, file = "cluster13_conserved_markers_thresholds_res9.csv")



cluster14_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 14, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster14_conserved_markers_thresholds_res9, file = "cluster14_conserved_markers_thresholds_res9.csv")


cluster15_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 15, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster15_conserved_markers_thresholds_res9, file = "cluster15_conserved_markers_thresholds_res9.csv")


cluster16_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 16, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster16_conserved_markers_thresholds_res9, file = "cluster16_conserved_markers_thresholds_res9.csv")


cluster17_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 17, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster17_conserved_markers_thresholds_res9, file = "cluster17_conserved_markers_thresholds_res9.csv")


cluster18_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 18, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster18_conserved_markers_thresholds_res9, file = "cluster18_conserved_markers_thresholds_res9.csv")

cluster19_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 19, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster19_conserved_markers_thresholds_res9, file = "cluster19_conserved_markers_thresholds_res9.csv")


cluster20_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 20, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster20_conserved_markers_thresholds_res9, file = "cluster20_conserved_markers_thresholds_res9.csv")

cluster21_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 21, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster21_conserved_markers_thresholds_res9, file = "cluster21_conserved_markers_thresholds_res9.csv")


cluster22_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 22, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster22_conserved_markers_thresholds_res9, file = "cluster22_conserved_markers_thresholds_res9.csv")

cluster23_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 23, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster23_conserved_markers_thresholds_res9, file = "cluster23_conserved_markers_thresholds_res9.csv")

cluster24_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 24, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster24_conserved_markers_thresholds_res9, file = "cluster24_conserved_markers_thresholds_res9.csv")

cluster25_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 25, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster25_conserved_markers_thresholds_res9, file = "cluster25_conserved_markers_thresholds_res9.csv")

cluster26_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 26, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster26_conserved_markers_thresholds_res9, file = "cluster26_conserved_markers_thresholds_res9.csv")


cluster27_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 27, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster27_conserved_markers_thresholds_res9, file = "cluster27_conserved_markers_thresholds_res9.csv")

cluster28_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 28, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster28_conserved_markers_thresholds_res9, file = "cluster28_conserved_markers_thresholds_res9.csv")

cluster29_conserved_markers_thresholds_res9 <- FindConservedMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = 29, grouping.var = "condition",
    verbose = TRUE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster29_conserved_markers_thresholds_res9, file = "cluster29_conserved_markers_thresholds_res9.csv")




###add cell annotations after investigation of conserved marker genes from each cluster
#cluster0 = M2_Macs
#cluster1 = NK2a
#cluster2 = CD8_Tcells
#cluster3 = NK3_ILC1
#cluster4 = DN_Tcells
#cluster5 = M1a_Macs
#cluster6 = M3_Macs
#cluster7 = NK1b
#cluster8 = NK1a
#cluster9 = DC2
#cluster10 = prolif_NKb
#cluster11 = DC1
#cluster12 = NK2b
#cluster13 = prolif_Mono
#cluster14 = cMono
#cluster15 = Bcells
#cluster16 = M1b_Macs
#cluster17 = ncMono
#cluster18 = stressed_T_cell
#cluster19 = prolif_NKa
#cluster20 = ILC3
#cluster21 = prolif_T_cell_a
#cluster22 = CD16+NK
#cluster23 = prolif_NKc
#cluster24 = fibroblast
#cluster25 = prolif_T_cell_b
#cluster26 = pDC
#cluster27 = GranMast
#cluster28 = Tol_DC
#cluster29 = plasma_Bcells

#replace names
new.cluster.ids <- c("M2_Macs",
"NK2a",
"CD8_Tcells",
"NK3_ILC1",
"DN_Tcells",
"M1a_Macs",
"M3_Macs",
"NK1b",
"NK1a",
"DC2",
"prolif_NKb",
"DC1",
"NK2b",
"prolif_Mono",
"cMono",
"Bcells",
"M1b_Macs",
"ncMono",
"stressed_T_cell",
"prolif_NKa",
"ILC3",
"prolif_T_cell_a",
"CD16+NK",
"prolif_NKc",
"fibroblast",
"prolif_T_cell_b",
"pDC",
"GranMast",
"Tol_DC",
"plasma_Bcells")
names(new.cluster.ids) <- levels(merged_EndoFT_combinedSCT_res9)
merged_EndoFT_combinedSCT_res9 <- RenameIdents(merged_EndoFT_combinedSCT_res9, new.cluster.ids)

#add idents back to object here
merged_EndoFT_combinedSCT_res9$annotation <- Idents(merged_EndoFT_combinedSCT_res9)

DimPlot(merged_EndoFT_combinedSCT_res9, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend()


#check cell cycle proliferating clusters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_EndoFT_combinedSCT_res9 <- CellCycleScoring(merged_EndoFT_combinedSCT_res9, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
DimPlot(merged_EndoFT_combinedSCT_res9, group.by = 'Phase')

cell_annotations_breakdown <- as.data.frame(table(Idents(merged_EndoFT_combinedSCT_res9)))
write.csv(cell_annotations_breakdown, file = "merged_EndoFT_annotated_cell_type_table.csv", row.names = FALSE)



#######Find DEGs across Endo vs FT within each cell type 
merged_EndoFT_combinedSCT_res9$annotation.condition <- paste(merged_EndoFT_combinedSCT_res9$annotation, merged_EndoFT_combinedSCT_res9$condition,
    sep = "_")
Idents(merged_EndoFT_combinedSCT_res9) <- "annotation_condition"

M2_Macs_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "M2_Macs_Endometrium", ident.2 = "M2_Macs_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(M2_Macs_Endo_vs_FT_markers, file = "M2_Macs_Endo_vs_FT_markers.csv")

NK2a_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "NK2a_Endometrium", ident.2 = "NK2a_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NK2a_Endo_vs_FT_markers, file = "NK2a_Endo_vs_FT_markers.csv")

CD8_Tcells_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "CD8_Tcells_Endometrium", ident.2 = "CD8_Tcells_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CD8_Tcells_Endo_vs_FT_markers, file = "CD8_Tcells_Endo_vs_FT_markers.csv")

NK3_ILC1_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "NK3_ILC1_Endometrium", ident.2 = "NK3_ILC1_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NK3_ILC1_Endo_vs_FT_markers, file = "NK3_ILC1_Endo_vs_FT_markers.csv")

DN_Tcells_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "DN_Tcells_Endometrium", ident.2 = "DN_Tcells_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DN_Tcells_Endo_vs_FT_markers, file = "DN_Tcells_Endo_vs_FT_markers.csv")

M1a_Macs_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "M1a_Macs_Endometrium", ident.2 = "M1a_Macs_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(M1a_Macs_Endo_vs_FT_markers, file = "M1a_Macs_Endo_vs_FT_markers.csv")

M3_Macs_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "M3_Macs_Endometrium", ident.2 = "M3_Macs_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(M3_Macs_Endo_vs_FT_markers, file = "M3_Macs_Endo_vs_FT_markers.csv")

NK1b_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "NK1b_Endometrium", ident.2 = "NK1b_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NK1b_Endo_vs_FT_markers, file = "NK1b_Endo_vs_FT_markers.csv")

NK1a_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "NK1a_Endometrium", ident.2 = "NK1a_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NK1a_Endo_vs_FT_markers, file = "NK1a_Endo_vs_FT_markers.csv")

DC2_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "DC2_Endometrium", ident.2 = "DC2_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DC2_Endo_vs_FT_markers, file = "DC2_Endo_vs_FT_markers.csv")

prolif_NKb_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "prolif_NKb_Endometrium", ident.2 = "prolif_NKb_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(prolif_NKb_Endo_vs_FT_markers, file = "prolif_NKb_Endo_vs_FT_markers.csv")

DC1_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "DC1_Endometrium", ident.2 = "DC1_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DC1_Endo_vs_FT_markers, file = "DC1_Endo_vs_FT_markers.csv")

NK2b_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "NK2b_Endometrium", ident.2 = "NK2b_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NK2b_Endo_vs_FT_markers, file = "NK2b_Endo_vs_FT_markers.csv")

prolif_Mono_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "prolif_Mono_Endometrium", ident.2 = "prolif_Mono_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(prolif_Mono_Endo_vs_FT_markers, file = "prolif_Mono_Endo_vs_FT_markers.csv")

cMono_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "cMono_Endometrium", ident.2 = "cMono_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cMono_Endo_vs_FT_markers, file = "cMono_Endo_vs_FT_markers.csv")

Bcells_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Bcells_Endometrium", ident.2 = "Bcells_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Bcells_Endo_vs_FT_markers, file = "Bcells_Endo_vs_FT_markers.csv")



M1b_Macs_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "M1b_Macs_Endometrium", ident.2 = "M1b_Macs_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(M1b_Macs_Endo_vs_FT_markers, file = "M1b_Macs_Endo_vs_FT_markers.csv")

ncMono_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "ncMono_Endometrium", ident.2 = "ncMono_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ncMono_Endo_vs_FT_markers, file = "ncMono_Endo_vs_FT_markers.csv")

stressed_T_cell_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "stressed_T_cell_Endometrium", ident.2 = "stressed_T_cell_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(stressed_T_cell_Endo_vs_FT_markers, file = "stressed_T_cell_Endo_vs_FT_markers.csv")

prolif_NKa_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "prolif_NKa_Endometrium", ident.2 = "prolif_NKa_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(prolif_NKa_Endo_vs_FT_markers, file = "prolif_NKa_Endo_vs_FT_markers.csv")

ILC3_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "ILC3_Endometrium", ident.2 = "ILC3_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ILC3_Endo_vs_FT_markers, file = "ILC3_Endo_vs_FT_markers.csv")

prolif_T_cell_a_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "prolif_T_cell_a_Endometrium", ident.2 = "prolif_T_cell_a_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(prolif_T_cell_a_Endo_vs_FT_markers, file = "prolif_T_cell_a_Endo_vs_FT_markers.csv")

CD16posNK_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "CD16+NK_Endometrium", ident.2 = "CD16+NK_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CD16posNK_Endo_vs_FT_markers, file = "NK_Endo_vs_FT_markers.csv")

prolif_NKc_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "prolif_NKc_Endometrium", ident.2 = "prolif_NKc_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(prolif_NKc_Endo_vs_FT_markers, file = "prolif_NKc_Endo_vs_FT_markers.csv")

fibroblast_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "fibroblast_Endometrium", ident.2 = "fibroblast_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(fibroblast_Endo_vs_FT_markers, file = "fibroblast_Endo_vs_FT_markers.csv")

prolif_T_cell_b_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "prolif_T_cell_b_Endometrium", ident.2 = "prolif_T_cell_b_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(prolif_T_cell_b_Endo_vs_FT_markers, file = "prolif_T_cell_b_Endo_vs_FT_markers.csv")

pDC_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "pDC_Endometrium", ident.2 = "pDC_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pDC_Endo_vs_FT_markers, file = "pDC_Endo_vs_FT_markers.csv")

GranMast_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "GranMast_Endometrium", ident.2 = "GranMast_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(GranMast_Endo_vs_FT_markers, file = "GranMast_Endo_vs_FT_markers.csv")

Tol_DC_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "Tol_DC_Endometrium", ident.2 = "Tol_DC_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tol_DC_Endo_vs_FT_markers, file = "Tol_DC_Endo_vs_FT_markers.csv")

plasma_Bcells_Endo_vs_FT_markers  <- FindMarkers(merged_EndoFT_combinedSCT_res9, assay = "SCT", ident.1 = "plasma_Bcells_Endometrium", ident.2 = "plasma_Bcells_First_Trimester", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(plasma_Bcells_Endo_vs_FT_markers, file = "plasma_Bcells_Endo_vs_FT_markers.csv")

####perform filtering of marker files
perCelltype_Endo_vs_FT_markers <- list.files(pattern = '*Endo_vs_FT_markers.csv')

#subset for 
#pvalue_adju < 0.05
#log2FC > 0.5 (up in Endo)
#log2FC < -1 (up in FT)


perCelltype_Endo_vs_FT_markers <- list.files(pattern = '*Endo_vs_FT_markers.csv')
for (i in perCelltype_Endo_vs_FT_markers){
    markers <- read.csv(i, row.names = 1)
    #pvalue_adjusted
    markers_subset <- subset(markers, p_val_adj < 0.05)
    #log2FC > 0.5 OR log2FC < -1
    markers_subset <- subset(markers_subset, avg_log2FC > .5 | avg_log2FC < -1)
    assign(i, markers_subset)
}


####write each one to a csv file
write.csv(M2_Macs_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_M2_Macs_Endo_vs_FT_markers.csv")
write.csv(NK2a_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_NK2a_Endo_vs_FT_markers.csv")
write.csv(CD8_Tcells_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_CD8_Tcells_Endo_vs_FT_markers.csv")
write.csv(NK3_ILC1_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_NK3_ILC1_Endo_vs_FT_markers.csv")
write.csv(DN_Tcells_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_DN_Tcells_Endo_vs_FT_markers.csv")
write.csv(M1a_Macs_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_M1a_Macs_Endo_vs_FT_markers.csv")
write.csv(M3_Macs_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_M3_Macs_Endo_vs_FT_markers.csv")
write.csv(NK1b_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_NK1b_Endo_vs_FT_markers.csv")
write.csv(NK1a_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_NK1a_Endo_vs_FT_markers.csv")
write.csv(DC2_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_DC2_Endo_vs_FT_markers.csv")
write.csv(prolif_NKb_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_prolif_NKb_Endo_vs_FT_markers.csv")
write.csv(DC1_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_DC1_Endo_vs_FT_markers.csv")
write.csv(NK2b_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_NK2b_Endo_vs_FT_markers.csv")
write.csv(prolif_Mono_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_prolif_Mono_Endo_vs_FT_markers.csv")
write.csv(cMono_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_cMono_Endo_vs_FT_markers.csv")
write.csv(Bcells_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_Bcells_Endo_vs_FT_markers.csv")
write.csv(M1b_Macs_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_M1b_Macs_Endo_vs_FT_markers.csv")
write.csv(ncMono_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_ncMono_Endo_vs_FT_markers.csv")
write.csv(stressed_T_cell_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_stressed_T_cell_Endo_vs_FT_markers.csv")
write.csv(prolif_NKa_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_prolif_NKa_Endo_vs_FT_markers.csv")
write.csv(ILC3_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_ILC3_Endo_vs_FT_markers.csv")
write.csv(prolif_T_cell_a_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_prolif_T_cell_a_Endo_vs_FT_markers.csv")
write.csv(CD16+NK_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_CD16+NK_Endo_vs_FT_markers.csv")
write.csv(prolif_NKc_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_prolif_NKc_Endo_vs_FT_markers.csv")
write.csv(fibroblast_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_fibroblast_Endo_vs_FT_markers.csv")
write.csv(prolif_T_cell_b_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_prolif_T_cell_b_Endo_vs_FT_markers.csv")
write.csv(pDC_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_pDC_Endo_vs_FT_markers.csv")
write.csv(GranMast_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_GranMast_Endo_vs_FT_markers.csv")
write.csv(Tol_DC_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_Tol_DC_Endo_vs_FT_markers.csv")
write.csv(plasma_Bcells_Endo_vs_FT_markers.csv, file = "subset_pval_log2fc_plasma_Bcells_Endo_vs_FT_markers.csv")

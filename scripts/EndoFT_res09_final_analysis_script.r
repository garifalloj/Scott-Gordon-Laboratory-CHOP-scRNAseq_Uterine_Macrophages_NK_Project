#!/usr/bin/env Rscript


library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
#options(Seurat.object.assay.version = 'v3')
#using v3 seurat assays by default


setwd(".")

for (file in c("Endo3", "FT4", "Endo2", "Endo4", "FT2", "FT3", "FT1", "Endo1")){
    seurat_data <- Read10X(data.dir = paste0(file, "/outs/filtered_feature_bc_matrix"))
    seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                    min.cells = 3,
                                    project = file)
    assign(file, seurat_obj)
}



#check cell numbers, add condition, add mito content. then filter
FT4[["condition"]] <- "First_Trimester"
Endo4[["condition"]] <- "Endometrium"
FT3[["condition"]] <- "First_Trimester"
Endo1[["condition"]] <- "Endometrium"
FT1[["condition"]] <- "First_Trimester"
FT2[["condition"]] <- "First_Trimester"
Endo3[["condition"]] <- "Endometrium"
Endo2[["condition"]] <- "Endometrium"

#merge all objects, violin plot QC
merged_Endo_FT <- merge(Endo1, y = c(Endo2, Endo4, Endo3, FT3, FT2, FT4, FT1))

##add mito content and Run QC
merged_Endo_FT[["percent.mt"]] <- PercentageFeatureSet(merged_Endo_FT, pattern = "^MT-")
pdf("VlnPlotQC_merged_Endo_FT_cellranger.pdf", height = 6, width = 10)
VlnPlot(merged_Endo_FT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Create .RData object to load at any time
#saveRDS(merged_Endo_FT, file="merged_Endo_FT_cellranger_rawSeur.rds")


#perform basic QC filtering
merged_Endo_FT <- subset(merged_Endo_FT, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
#81,234 cells down to 62,229 post filtering


####split object and find anchors here.
#split object into Endo and FT
merged_Endo_FT.list <- SplitObject(merged_Endo_FT, split.by = "condition")

Endo_split <- merged_Endo_FT.list[["Endometrium"]]
FT_split <- merged_Endo_FT.list[["First_Trimester"]]


#########run this on BOTH Endo and FT splits individually
# normalize and run dimensionality reduction on control dataset
#CANNOT use orig.ident in this context
#code based on Seurat integration vignette stim vs control
Endo_split <- SCTransform(Endo_split, vst.flavor = "v2", verbose = TRUE, vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30, verbose = TRUE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = TRUE) %>%
    FindClusters(resolution = 0.7, verbose = TRUE)

#check dimplot by orig.ident AND by percent.mt
DimPlot(Endo_split, reduction = 'umap')
DimPlot(Endo_split, reduction = 'umap', group.by = "orig.ident")

FeaturePlot(Endo_split, features = 'percent.mt')
FeaturePlot(Endo_split, features = 'nCount_SCT')

#####repeat process for FT_split
FT_split <- SCTransform(FT_split, vst.flavor = "v2", verbose = TRUE, vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30, verbose = TRUE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = TRUE) %>%
    FindClusters(resolution = 0.7, verbose = TRUE)

#check dimplot by orig.ident AND by percent.mt
DimPlot(FT_split, reduction = 'umap')
DimPlot(FT_split, reduction = 'umap', group.by = "orig.ident")

FeaturePlot(FT_split, features = 'percent.mt')
FeaturePlot(FT_split, features = 'nCount_SCT')



####then integrate across both separate lists using residuals
merged_Endo_FT.list <- list(Endo_split = Endo_split, FT_split = FT_split)
features <- SelectIntegrationFeatures(object.list = merged_Endo_FT.list, nfeatures = 3000)
merged_Endo_FT.list <- PrepSCTIntegration(object.list = merged_Endo_FT.list, anchor.features = features)
#To integrate the two datasets, we use the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
##CAN IGNORE WARNINGS regarding different cells and/or features from existing assay SCT, https://github.com/satijalab/seurat/issues/7145


##this step can take a LONG time and lot of memory!!!
merged_Endo_FT.anchors <- FindIntegrationAnchors(object.list = merged_Endo_FT.list, normalization.method = "SCT",
    anchor.features = features)
merged_Endo_FT.combined.sct <- IntegrateData(anchorset = merged_Endo_FT.anchors, normalization.method = "SCT")
#Perform an integrated analysis
#Now we can run a single integrated analysis on all cells:

merged_Endo_FT.combined.sct <- RunPCA(merged_Endo_FT.combined.sct, verbose = FALSE)
merged_Endo_FT.combined.sct <- RunUMAP(merged_Endo_FT.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
merged_Endo_FT.combined.sct <- FindNeighbors(merged_Endo_FT.combined.sct, reduction = "pca", dims = 1:30)
merged_Endo_FT.combined.sct <- FindClusters(merged_Endo_FT.combined.sct, resolution = 0.9)
#check final resolution, should be 0.9 after testing


#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(merged_Endo_FT.combined.sct, reduction = "umap", split.by = "condition")

#generate integrated dimplots

DimPlot(merged_Endo_FT.combined.sct, reduction = 'umap')
DimPlot(merged_Endo_FT.combined.sct, reduction = 'umap', split.by = "condition")
DimPlot(merged_Endo_FT.combined.sct, reduction = 'umap', group.by = "condition")

FeaturePlot(merged_Endo_FT.combined.sct, features = "percent.mt")

#need to run this before finding marker gene comparisons
merged_Endo_FT.combined.sct <- PrepSCTFindMarkers(merged_Endo_FT.combined.sct)



saveRDS(merged_Endo_FT.combined.sct, file = "merged_Endo_FT.combined.sct.rds")


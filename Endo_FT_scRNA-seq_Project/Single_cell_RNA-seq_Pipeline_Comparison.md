Single-cell RNA-seq Pipeline Comparison Analysis

Pipelines analyzed:

#for each mapping and cell calling algorithm, the 10x Genomics human GRCh38 default reference was utilized. In cases of pseudo-alignment (kallisto-BUS, alevin), the same reference was utilized to create an appropriate reference index.

1_cellranger_standard_seurat_anchor_integration_byCondition
2_starSolo_standard_seurat_anchor_integration_byCondition
3_kallistoBus_standard_seurat_anchor_integration_byCondition
4_alevin_standard_seurat_anchor_integration_byCondition

5_cellranger_harmony_integration_byCondition
6_starSolo_harmony_integration_byCondition
7_kallistoBus_harmony_integration_byCondition
8_alevin_harmony_integration_byCondition

#pipeline outputs were used to create cell x gene barcode matrices and used as input into seurat v5. Matrices were filtered to only keep cells expressing between 200-4000 genes and less than 20% mitochondrial content. 

#for the kallisto-BUS pipeline, the dropUtils algorithm was used to determine which barcodes belonged to cells (FDR <= 0.01 for isCell). Transcripts from the pseudo-alignment were converted to their appropriate gene symbols and used in the seurat object downstream.


After each seurat object was processed using either the seurat-anchor integration or the harmony integration method by condition (Endo vs. FT), singleR R package was used to determine the likely cell annotation for each cell based on the published Vento-tormo dataset. Cell annotation labels were added to each seurat object as meta data columns.

Pending: all cell annotated seurat objects were merged in R to assess differences in cell proportions and their significance. Note: gene expression and other assessments on the combined dataset were not possible as the same dataset was replicated and would not be meaningful.

Propellor R package was then used to assess differences in cell proportions and their significance across the varying pipelines.

Questions:
sample FT149 seems to have strong independent clustering. Examine the #UMIs and possible regress out differences in UMI.
alternative: seurat SCT normalization before integration methods would regress out UMIs naturally.

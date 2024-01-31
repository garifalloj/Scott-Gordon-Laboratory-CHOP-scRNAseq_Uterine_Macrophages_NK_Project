# Uterine macrophages and NK cells exhibit population-and gene-level changes after implantation but maintain pro-invasive properties

Abstract here

### Bioinformatics Strategy

Raw data from sequencing was processed using the cell ranger pipeline (10x Genomics) to demultiplex data into per library cell barcode matrices. Cell barcode matrices were used as input into Seurat, filtered for quality thresholds and aggregated. Aggregated data was normalized using the sctransform normalization method. To account for differences between the Endometrium and first trimester samples, an integrated method seurat anchor integration, was utilized to ensure robustness of clustering into appropriate cell type-clusters and account for possible batch effects.  

### Code Availability

This repository will contain all associated scripts to replicate the analyses provided in the manuscript. All analyses were performed in the R statistical computing environment. Libraries utilized:

```
library(Seurat)
library(dplyr)
library(patchwork)
library(UCell)
library(singleR)
library(singleCellExperiment)
library(monocle)
library(ggplot2)
library(Rcolorbrewer)

```

R version 4.2.3 was utilized for this analysis. Additional sessionInfo() is stored as sessionInfo.txt in this repository.





### Individual R scripts are in the scripts subfolder

Scripts organized by process.




### Processed Outputs and Supplemental figures

Per annotated clusters marker genes files and additional figures will be described here as needed.




### Public Data

Data is available through the NCBI GEO portal, accession number GEO1241532532. GEO uploads in the fastq files and the processed count matrices for each library. 






### Additional Processed Data

Additional processed data for the fresh and frozen analysis can be located here [URL].



 
 
 ### Notes and Information
 

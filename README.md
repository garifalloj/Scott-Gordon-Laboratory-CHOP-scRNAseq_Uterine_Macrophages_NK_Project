# Uterine macrophages and NK cells exhibit population-and gene-level changes after implantation but maintain pro-invasive properties

Prior to pregnancy, hormonal changes lead to cellular adaptations in the endometrium allowing for embryo implantation. Critical for successful pregnancy establishment, innate immune cells constitute a significant proportion of uterine cells prior to arrival of the embryo and throughout the first trimester in humans and animal models. Abnormal uterine immune cell function during implantation is believed to play a role in multiple adverse pregnancy outcomes. Current work in humans has focused on uterine immune cells present after pregnancy establishment, and limited in vitro models exist to explore unique functions of these cells. With single-cell RNA-sequencing (scRNAseq), we comprehensively compared the human uterine immune landscape of the endometrium during the window of implantation and the decidua during the first trimester of pregnancy. We uncovered global and cell-type-specific gene signatures for each timepoint. Immune cells in the endometrium prior to implantation expressed genes associated with immune metabolism, division, and activation. In contrast, we observed widespread interferon signaling during the first trimester of pregnancy. We also provide evidence of specific inflammatory pathways enriched in pre- and post-implantation macrophages and natural killer (NK) cells in the uterine lining. Using our novel implantation-on-a-chip (IOC) to model human implantation ex vivo, we demonstrate for the first time that uterine macrophages strongly promote invasion of extravillous trophoblasts (EVTs), a process essential for pregnancy establishment. Pre- and post-implantation uterine macrophages promoted EVT invasion to a similar degree as pre- and post-implantation NK cells on the IOC. This work provides a foundation for further investigation of the individual roles of uterine immune cell subtypes present prior to embryo implantation and during early pregnancy, which will be critical for our understanding of pregnancy complications associated with abnormal trophoblast invasion and placentation.

### Experimental Design

Sort-purified CD45+ cells from mid-secretory (implantation window) endometrium were compared to Sort-purified CD45+ cells from first-trimester elective termination decidua.

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

R version 4.2.3 was utilized for this analysis.





### Individual R scripts are in the scripts subfolder

Scripts organized by process.




### Processed Outputs and Supplemental figures

Per annotated clusters marker genes files and additional figures will be described here as needed.




### Public Data

Data is available through the NCBI GEO portal, accession number GSE255282 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255282). Data is uploaded as the fastq files and the processed count matrices for each library. 






 
 
 ### Additional Notes and Information
 

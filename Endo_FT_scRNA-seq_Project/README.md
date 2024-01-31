### Uterine macrophages and NK cells exhibit population-and gene-level changes after implantation but maintain pro-invasive properties

Abstract here

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

```

# Individual R scripts are in the scripts subfolder






 
 8 total samples. The experimental design was to look at sorted CD45+ immune cells from: 

endometrial biopsies from the “window of implantation” (around the time of ovulation, pre pregnancy): Endo_63, Endo_69, Endo_71, Endo_72 

the same tissue, now called decidua, from first trimester elective terminations of pregnancy: FT_130, FT_131, FT_135, FT_149 

 ### Biobox online platform used for intitial data analysis pipeline and clustering.
 

 

I did spend a ton of time assigning cluster identity based on the Vento-Tormo dataset from first-trimester pregnancy and other resources and came up with the attached excel spreadsheet of cluster identities. The tab names go Supercluster_subcluster. Hopefully using the same parameters, we can get the same clusters. Regardless, I’ll have to look through and refine the cluster assignments a bit. 

 

I’m not too sure with this number of cells/samples how much higher we can/should push the cluster resolution. Do you have ideas about when we get too many clusters with too little power? I worry we get to a point where the subclusters get ridiculous and don’t help us biologically. 

 

The primary objectives are to generate figures showing the following. Let me know if you like particular graph styles more than others, seems like there are many good options from this paper and others (https://www.nature.com/articles/s41467-021-21892-z): 

A single master map showing superclusters of the basic cell types considering all 8 samples: 

T cells 

B cells 

NK cells and other innate lymphoid cells 

Granulocytes 

Monocytes 

Macrophages 

Dendritic cells 

 

Kind of like this, though some of the yellows and oranges are brutal. 

 

 

A single master map showing all subclusters 

Same as above, just all subclusters shown. 

I will have to confirm cell type annotations 

 

Kind of like this one from the Vento-Tormo paper: 

 

 

A monocle analysis showing possible developmental relationships among the subclusters. Might help me name some of the less clear monocyte and macrophage subsets. 

 

A figure like this one showing the top few genes identifying each subcluster. It would get crazy showing all the subclusters together, so I’d like to make one panel each on 

Just T cell subclusters 

Just NK cell and ILC subclusters  

Just monocyte subclusters 

Just dendritic cell subclusters 

Just macrophage subclusters 

 

 

 

The integrated analysis comparing the two UMAPs for pre-implantation and first-trimester like my home cooked one above. This is going to be the meat of the paper, basically how things change quantitatively and qualitatively between pre- and post-implantation of the embryo. From the software, I have a list of DEGs between the Endo and FT clusters. Hopefully these will stay more or less the same, but we’ll have to see. Once I have those, we can run some thru gene ontology or other grouping algorithm. 

We need one table comparing the number of events in each cluster between the endo samples and the FT samples, plus a statistical test. Likely a table with number of cells in one column, fold change of FT/endo in the next column, p value in the next column, and each row designating the subcluster identity. Ideally, this would be organized by fold change, from most enriched in FT to most enriched in Endo.  

I really like this layout in C, showing representation of each cluster by category, in our case it would be Endo and FT. I think the T cells, B cells, ILCs, and NK cells (“lymphoid cells”) should be one bar plot and the rest of the clusters (“myeloid cells”) should be a second. (https://www.researchgate.net/figure/Integration-analysis-of-single-cell-RNA-sequencing-datasets-from-SLE-patients-and_fig1_360084554):   

A separate table, similar layout as above, showing number of genes significantly differing between each cluster in Endo vs FT. (x many genes significantly up in Endo, y genes significantly up in FT) 

 

One heatmap each showing the top significant genes enriched in Endo samples: 

Just T cell and B cell subclusters 

Just NK cell and ILC subclusters  

Just monocyte subclusters 

Just dendritic cell subclusters 

Just macrophage subclusters 

 

One heatmap each showing the top significant genes enriched in FT samples: 

Just T cell and B cell subclusters 

Just NK cell and ILC subclusters  

Just monocyte subclusters 

Just dendritic cell subclusters 

Just macrophage subclusters 

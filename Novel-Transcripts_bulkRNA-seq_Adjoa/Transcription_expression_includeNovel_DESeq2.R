#!bin/sh

#1st step use prepDE.py3 to extract read count information from stringtie -e outputs (in each folder)
#prepDE.py derives hypothetical read counts for each transcript from the coverage values estimated by StringTie for each transcript, by using this simple formula: #reads_per_transcript = coverage * transcript_len / read_len

#generate input file in directory, 
ll *.stringtie_output_GRCm38_wMerged.gtf > inputList_forDESeq2Prep.txt

#input = text file with sampleID and path to its gtf
# -l average read length = 150bp
conda activate
python3 prepDE.py3 -i inputList_forDESeq2Prep.txt -l 150 -g gene_count_matrix_fromMerged.csv -t transcript_count_matrix_fromMerged.csv
#successful, took 3 mins does not need extra memory


#These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR (using the DESeqDataSetFromMatrix and DGEList functions, respectively).
#DESeqDataSetFromMatrix function


#example workflow
setwd("/Users/garifalloj/Desktop/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts/ballgown_DTE_analysis_06222023")

library("DESeq2")
#Load gene(/transcript) count matrix and labels
countData <- as.matrix(read.csv("transcript_count_matrix_fromMerged.csv", row.names="transcript_id"))
colData <- read.csv("design_stringtie_cond_treat.csv", header = TRUE, sep = ",", row.names=1)
#colData <- read.csv(PHENO_DATA, sep="\t", row.names=1)
#Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
#[1] TRUE
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
[1] TRUE
#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData, design = ~ condition)
#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
res <- results(dds)
#Sort by adjusted p-value and display
(resOrdered <- res[order(res$padj), ])

write.csv(resOrdered, file = "WTNK_vs_CD122mut_DEseq2_transcripts.csv")

#possibly filter to remove low transcript outliers
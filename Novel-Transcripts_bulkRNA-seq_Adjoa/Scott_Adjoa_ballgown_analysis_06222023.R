library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)


setwd("/Users/garifalloj/Desktop/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts/ballgown_DTE_analysis_06222023")

#R 4.2.1 on local Rstudio

#add sample groups phenotype data
pheno_data_cond_treat <- read.csv("design_stringtie_cond_treat.csv", header = TRUE, sep = ",")
#condition = WT, IFNARKO, CD122mut, or WT_NK
#treatment = no_stim, stim_10pIFNa


##must create manually pData 
#check sample order with:
#sampleNames(bg)


bg = ballgown(dataDir=".", samplePattern='HTS034', meas='all', pData=pheno_data_cond_treat)
bg
#163612 transcripts, 22 samples
#try to read in at same time with pData

##now filter:
##possibly repeat from here with filtering and additional parameters
bg = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
#this removes all genes with a variance across samples less than 1



transcript_fpkm = texpr(bg, 'FPKM')
transcript_cov = texpr(bg, 'cov')
whole_tx_table = texpr(bg, 'all')
exon_mcov = eexpr(bg, 'mcov')
junction_rcount = iexpr(bg)
whole_intron_table = iexpr(bg, 'all')
gene_expression = gexpr(bg)



exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

results_transcripts = stattest(bg, feature='transcript', meas='FPKM', covariate='condition')
#head(stat_results)


##needs to be named "bg" for stattest to work correctly??

results_genes = stattest(bg, feature="gene", meas="FPKM", covariate='condition')


#Add gene names and gene IDs to the results_transcripts data frame:
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg), geneIDs=ballgown::geneIDs(bg), results_transcripts)


#Sort the results from the smallest P value to the largest:

results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)


#Write the results to a csv file that can be shared and distributed:

write.csv(results_transcripts, "transcript_results.csv", row.names=FALSE)

write.csv(results_genes, "gene_results.csv", row.names=FALSE)


#Identify transcripts and genes with a q value <0.05:

subset(results_transcripts,results_transcripts$qval<0.05)

subset(results_genes,results_genes$qval<0.05)


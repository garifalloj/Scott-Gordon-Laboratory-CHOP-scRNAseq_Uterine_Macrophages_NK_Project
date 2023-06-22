Project Description:

Examine the role of macrophages and NK cells after stim and no stim for DE genes and DE transcripts. Examine the role of novel transcripts within the CD122mutants

bulkRNAseq libraries performed with TruSeq total RNAseq with ribo depletion. 150bp x 150bp read length

GRCm38/mm10 reference for previous qPCR experiments.

Bioinformatics Updates:
libraries aligned with STAR 2.7.10 (ensembl GRCm38 with 150bp read length index) and salmon gencode for pseudo alignment/quantification.

gene counting in Rsubread with featureCounts algorithm.

Downstream DE in DEseq2. 


Novel transcript assembly:

Overall pipeline:
https://www.nature.com/articles/nprot.2016.095#Tab1

Transcripts assembled by stringtie (JHU)

Need to identify transcripts from the stringtie based data with isoformAnalyzerR package: issue with annotation fasta/references, can't annotate ORFs
Could possibly identify differentially expressed transcripts across the CD122 mutants in DEseq2 and/or visualize those transcripts in Ballgown R package.

Maybe SCISSOR for transcript shape changes? Would be applicable for particular genes. Maybe Il2rb? or any others that are highly DE in CD122mutants?
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7804101/
Can it be run on multiple genes and samples??


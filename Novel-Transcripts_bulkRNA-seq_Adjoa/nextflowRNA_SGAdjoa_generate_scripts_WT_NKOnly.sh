#!/bin/sh
#must be full pwd to main fastqs
fastq_directory=/mnt/isilon/neonatology_bioinformatics/projects/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts

#cd into project directory that has per sample folders fastqs
cd $fastq_directory

for fastq_file in HTS034*_WT_NK_*R1_001.fastq.gz; do
    (
        sn=$(echo $fastq_file | cut -f 1-6 -d'_');
        echo -e "sample,fastq_1,fastq_2,strandedness
$sn,$fastq_directory/$fastq_file,$fastq_directory/$fastq_file,auto" > /mnt/isilon/neonatology_bioinformatics/projects/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts/nf_RNAseq_salmon/"$sn"_nfcoreRNA_samplesheet.csv;
    )
    done

for fastq_file in HTS034*_WT_NK_*R1_001.fastq.gz; do
    (
        sn=$(echo $fastq_file | cut -f 1-6 -d'_');
        echo -e "#!/bin/sh
#SBATCH -t 7:00:00
#SBATCH --mem=40G
#SBATCH -c 3
module load Java/15.0.1
module load singularity/3.7.0

export NXF_SINGULARITY_CACHEDIR="/mnt/isilon/neonatology_bioinformatics/programs/singularity3.7_images_temp/singularity_image"

nextflow run nf-core/rnaseq --input /mnt/isilon/neonatology_bioinformatics/projects/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts/nf_RNAseq_salmon/"$sn"_nfcoreRNA_samplesheet.csv --outdir /mnt/isilon/neonatology_bioinformatics/projects/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts/nf_RNAseq_salmon/$sn --pseudo_aligner salmon --salmon_index /mnt/isilon/neonatology_bioinformatics/reference_genomes/bulk_RNAseq/salmon_GRCm38_vM25_decoyAware_index --fasta /mnt/isilon/neonatology_bioinformatics/reference_genomes/bulk_RNAseq/GRCm38.primary_assembly.genome.fa.gz --gtf /mnt/isilon/neonatology_bioinformatics/reference_genomes/bulk_RNAseq/gencode.vM25.primary_assembly.annotation.gtf --gencode  -profile singularity" > /mnt/isilon/neonatology_bioinformatics/projects/Scott_Gordon_Adjoa_bulkRNAseq_novel_transcripts/nf_RNAseq_salmon/"$sn"_nfcoreRNA_salmon_gencode.sh;
    )
    done

##after running, need to change second fastq input file to "R2"
#also need to delete WT_NK*.sh, WT_NK*.csv and repeat for cut -f 1-6 for sn command
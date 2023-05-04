##examine BAM files in IGV web browswer
#subset for only chr15: 78363456-78396569
chr15: 78300000-78400000

HTS034_02_WT_Mouse2_NSAligned.sortedByCoord.out.bam
HTS034_02_WT_Mouse2_NSAligned.toTranscriptome.out.bam

module load sam-bcf-tools/1.6
samtools index HTS034_02_WT_Mouse2_NSAligned.sortedByCoord.out.bam

##need to index all BAM and transcriptome BAM files before viewing in IGV.


#try just chr15

samtools view -b HTS034_02_WT_Mouse2_NSAligned.sortedByCoord.out.bam > chr15_WT_mouse2_NS_subset.bam
chr15:78,300,000-78400000 > chr15_Il2rb_region.bam
chr15:78,477,256-78,513,621
###goal, index all BAM and transcriptome files before viewing IGV

#interactive node 
srun -t 5:00:00 --mem=30G --pty bash

module load sam-bcf-tools/1.6




#example:
#samtools view input.bam "Chr10:18000-45500" > output.bam


#index all BAM files NOT transcriptome files, then subset for region of interest. Merge appropriate replicates, then visualize in IGV:
for i in *.Aligned.sortedByCoord.out.bam; do samtools index $i; done

for i in *.Aligned.sortedByCoord.out.bam; do samtools view -b $i "15:78477256-78513621" > $i.chr15_Il2rb_region.bam; done

#try manually:
samtools view HTS034_13_CD122mut_Mouse1_NS.Aligned.sortedByCoord.out.bam "15:78477256-78513621" > test_13_CD122_NS_reads.txt

samtools view -c HTS034_13_CD122mut_Mouse1_NS.Aligned.sortedByCoord.out.bam "15:78477256-78513621"

need to index all

for i in *Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam; do samtools index $i; done


#merge all subsetted BAMs for each group, 3 per group. then pull screen shots.

samtools merge WT_no_stim_merged_Il2rb_locus.bam HTS034_*_WT_Mouse*_NS.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index WT_no_stim_merged_Il2rb_locus.bam
samtools view -c WT_no_stim_merged_Il2rb_locus.bam "15:78477256-78513621"
#then index this BAM, view in IGV, save screenshot


##repeat for WT_stim: HTS034_04_WT_Mouse1_10pIFNa.Aligned.sortedByCoord
samtools merge WT_stim_10pIFNa_merged_Il2rb_locus.bam HTS034_*_WT_Mouse*_10pIFNa.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index WT_stim_10pIFNa_merged_Il2rb_locus.bam
samtools view -c WT_stim_10pIFNa_merged_Il2rb_locus.bam "15:78477256-78513621"


##repeat for IFNARKO NS
samtools merge IFNARKO_no_stim_merged_Il2rb_locus.bam HTS034_*_IFNARKO_Mouse*_NS.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index IFNARKO_no_stim_merged_Il2rb_locus.bam
samtools view -c IFNARKO_no_stim_merged_Il2rb_locus.bam "15:78477256-78513621"

#repeat for IFNARKO stim HTS034_10_IFNARKO_Mouse1_10pIFNa.Aligned.sortedByCoord.out.bam 
samtools merge IFNARKO_stim_10pIFNa_merged_Il2rb_locus.bam HTS034_*_IFNARKO_Mouse*_10pIFNa.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index IFNARKO_stim_10pIFNa_merged_Il2rb_locus.bam
samtools view -c IFNARKO_stim_10pIFNa_merged_Il2rb_locus.bam "15:78477256-78513621"
#spans exon1-2 and exon10, depending on which transcript


##repeat for CD122mut stim HTS034_16_CD122mut_Mouse1_10pIFNa.Aligned.sortedByCoord.out.bam
samtools merge CD122mut_stim_10pIFNa_merged_Il2rb_locus.bam HTS034_*_CD122mut_Mouse*_10pIFNa.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index CD122mut_stim_10pIFNa_merged_Il2rb_locus.bam
samtools view -c CD122mut_stim_10pIFNa_merged_Il2rb_locus.bam "15:78477256-78513621"


##repeat for CD122mut no stim HTS034_14_CD122mut_Mouse2_NS.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools merge CD122mut_no_stim_merged_Il2rb_locus.bam HTS034_*_CD122mut_Mouse*_NS.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index CD122mut_no_stim_merged_Il2rb_locus.bam
samtools view -c CD122mut_no_stim_merged_Il2rb_locus.bam "15:78477256-78513621"


###for CD122mutants, no reads in exons 1-2, spread out 

#produces output for every bp coordinate, not needed
samtools depth -r 15:78477256-78513621 CD122mut_no_stim_merged_Il2rb_locus.bam


####################repeat for NK cell subset
for i in *.Aligned.sortedByCoord.out.bam; do samtools view -b $i "15:78477256-78513621" > $i.chr15_Il2rb_region.bam; done
#need to index all
for i in *Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam; do samtools index $i; done

###repeat merge replicates for stim/no stim
####WT_NK_NS: HTS034_19_WT_NK_Mouse1_NS.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools merge WT_NK_no_stim_merged_Il2rb_locus.bam HTS034_*_WT_NK_Mouse*_NS.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index WT_NK_no_stim_merged_Il2rb_locus.bam
samtools view -c WT_NK_no_stim_merged_Il2rb_locus.bam "15:78477256-78513621"
#70,122 reads

###WT_NK_stim: HTS034_21_WT_NK_Mouse1_10pIFNa.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools merge WT_NK_stim_10pIFNa_merged_Il2rb_locus.bam HTS034_*_WT_NK_Mouse*_10pIFNa.Aligned.sortedByCoord.out.bam.chr15_Il2rb_region.bam
samtools index WT_NK_stim_10pIFNa_merged_Il2rb_locus.bam
samtools view -c WT_NK_stim_10pIFNa_merged_Il2rb_locus.bam "15:78477256-78513621"
#71,065
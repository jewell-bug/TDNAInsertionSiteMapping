#!/bin/bash 

# Set GTF path
GTF_PATH="../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/braker.gtf"

#Sort the mutants BED file
bedtools sort -i mutants.bed > mutants_sorted.bed

#Sort the GTF file
bedtools sort -i $GTF_PATH > braker_sorted.gtf

# direct intersections 
bedtools intersect -a mutants_sorted.bed -b braker_sorted.gtf -wa -wb > direct_intersections.txt

# closest features search 
bedtools closest -a mutants_sorted.bed -b braker_sorted.gtf -d > closest_features.txt



PFAM_FILE="../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/gene_pfam_all.tsv"

# Search for g5275, g6713, and g6273 in the Pfam table
grep -E "g5275|g6713|g6273" $PFAM_FILE

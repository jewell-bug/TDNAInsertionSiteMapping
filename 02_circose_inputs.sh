#!/bin/bash 
#Load the software modules
module load samtools/1.19
module load bedtools/2.31.1
# Extract the chromosome/contig sizes from  BAM file
samtools view -H bwa_align_wt_tdna/B7_combined.bam | grep "@SQ" | sed 's/@SQ\tSN://g' | sed 's/\tLN:/\t/g' > genome_sizes.txt

# Create 1,000 base-pair  windows across  genome
bedtools makewindows -g genome_sizes.txt -w 1000 > windows.bed

# Calculate the average read depth for each 1kb window
bedtools coverage -a windows.bed -b bwa_align_wt_tdna/B7_combined.bam > circos_coverage_input.txt


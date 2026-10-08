#!/bin/bash 



module load blast/2.13.0
 Make a BLAST database of WT genome assembly
#makeblastdb -in ../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/CF40_WT_contigs.fasta -dbtype nucl -out wt_db

# Query  mutant insertion coordinates against the WT database

module load samtools/1.19
samtools faidx spades/B7_combined_R1_001-fastp.fastq.gz/B7_contigs.fasta  NODE_7_length_898147_cov_22.871728:430000-445000 > mutant_query.fasta

samtools faidx spades/C2_combined_R1_001-fastp.fastq.gz/C2_contigs.fasta NODE_19_length_525580_cov_23.332626:70000-82000 > mutant_C2_query.fasta
# Extract 15kb chunk (bp 1 to 15,000) around the insertion on E3 contig

samtools faidx spades/E3_combined_R1_001-fastp.fastq.gz/E3_contigs.fasta \
  NODE_2_length_1038953_cov_20.436220:1-15000 > mutant_E3_query.fasta

samtools faidx spades/F8_combined_R1_001-fastp.fastq.gz/F8_contigs.fasta \
  NODE_2_length_1176253_cov_22.155932:1165000-1175000 > mutant_F8_query.fasta

samtools faidx spades/G9_combined_R1_001-fastp.fastq.gz/G9_contigs.fasta \
  NODE_2_length_1041895_cov_24.673380:1030000-1040000 > mutant_G9_query.fasta

blastn -query mutant_G9_query.fasta -db wt_db -outfmt 6 > G9_mutant_map_wt.txt

blastn -query mutant_E3_query.fasta -db wt_db -outfmt 6 > E3_mutant_map_wt.txt

blastn -query mutant_C2_query.fasta -db wt_db -outfmt 6 > C2_mutant_map_wt.txt

blastn -query mutant_F8_query.fasta -db wt_db -outfmt 6 > F8_mutant_map_wt.txt

#B7


awk '$1=="NODE_3_length_997933_cov_20.759324" && $3=="transcript" && $4>=535000 && $5<=545000' ../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/braker.gtf

#C2

awk '$1=="NODE_37_length_267140_cov_20.161194" && $3=="transcript" && $4>=185000 && $5<=200000' ../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/braker.gtf

#F8
awk '$1=="NODE_2_length_1035884_cov_20.579413" && $3=="transcript" && $4<=1035884 && $5>=1030994' ../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/braker.gtf


#E3
awk '$1=="NODE_2_length_1035884_cov_20.579413" && $3=="transcript" && $4<=1033932 && $5>=1027248' ../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/braker.gtf


#G9 
awk '$1=="NODE_2_length_1035884_cov_20.579413" && $4<=1035884 && $5>=1030352' ../../CF40_seq_27_08_25/spades/CF40_165_028_R1_001-fastp.fastq.gz/braker/braker.gtf

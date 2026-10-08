#!/bin/bash 


# Load  BLAST 
module load blast/2.13.0

# database out of  SPAdes assembly file
makeblastdb -in spades/E3_combined_R1_001-fastp.fastq.gz/contigs.fasta -dbtype nucl -out spades_db_E3

# Align T-DNA sequence against the assembly database
blastn -query tdna_genes.fa -db spades_db_E3 -outfmt 6 -out gene_tdna_integration_results_E3.txt

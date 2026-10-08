#!/bin/bash 
module load bwa-mem2/2.1
module load bwa-mem2/2.1
module load samblaster/0.1.24
module load samtools/1.19

# set directories
INDEXDIR=bwa_index_wt_withtdna
mkdir -p $INDEXDIR

GENOME=CF40_WT_contigs.fasta


# Create index for genome
    # requires significant memory for indexing
bwa-mem2 index \
   -p $INDEXDIR/AP \
   $GENOME


#set directories
SAMPDIR=fastp_data

OUTDIR=bwa_align_wt_tdna
mkdir -p $OUTDIR

INDEX=bwa_index_wt_withtdna/AP

# sample ID list
SAMPLELIST=(E3_combined C2_combined)

# extract one sample ID
SAMPLE=${SAMPLELIST[$SLURM_ARRAY_TASK_ID]} 

# create read group string
RG=$(echo \@RG\\tID:$SAMPLE\\tSM:$SAMPLE)

# execute the alignment pipe:
bwa-mem2 mem -t 7 -R ${RG} ${INDEX} ${SAMPDIR}/${SAMPLE}_R1_001-fastp.fastq.gz $SAMPDIR/${SAMPLE}_R2_001-fastp.fastq.gz | \
	samblaster | \
	samtools view -S -h -u - | \
	samtools sort -T ${OUTDIR}/${SAMPLE}.temp -O BAM >$OUTDIR/${SAMPLE}.bam 

# index alignment file
samtools index ${OUTDIR}/${SAMPLE}.bam


##coverage 
samtools coverage bwa_align_wt_tdna/B7_combined.bam

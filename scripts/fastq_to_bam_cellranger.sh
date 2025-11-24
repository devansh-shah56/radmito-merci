#!/bin/bash

# ============================
# Single-cell BAM Pipeline
# FASTQ.gz → BAM using Cell Ranger
# Author: RADmito
# ============================


# Update sample IDs here
SAMPLES=("SRR15339529" "SRR15339530","SRR15339531","SRR15339532")

# Base directories
PROJECT_DIR="/NAS/ravikantv/dh607_project"
FASTQ_DIR="$PROJECT_DIR/fastq_files"
REF_DIR="$PROJECT_DIR/ref/refdata-gex-GRCh38-2024-A"

# Resources (your server: 128 cores, 1 TB RAM)
CORES=32
MEMORY=256

# -------- RUN PIPELINE --------
echo "Starting Cell Ranger BAM Pipeline..."
echo "Using $CORES cores and $MEMORY GB RAM"

for SAMPLE in "${SAMPLES[@]}"
do
    echo "----------------------------------------"
    echo "Processing Sample: $SAMPLE"
    echo "----------------------------------------"
############################################
# before running cellranger rename .fastq.gz file like SRR15339529_S1_L001_I1_001.fastq.gz
# before running cellranger rename .fastq.gz file like SRR15339529_S1_L001_R1_001.fastq.gz
# before running cellranger rename .fastq.gz file like SRR15339529_S1_L001_R2_001.fastq.gz

# before running cellranger rename .fastq.gz file like SRR15339530_S1_L001_I1_001.fastq.gz
# before running cellranger rename .fastq.gz file like SRR15339530_S1_L001_R1_001.fastq.gz
# before running cellranger rename .fastq.gz file like SRR15339530_S1_L001_R2_001.fastq.gz

# so on for other two samples SRR15339531, SRR15339532
###########################################
    cellranger count \
        --id=${SAMPLE}_bamrun \
        --transcriptome=$REF_DIR \
        --fastqs=$FASTQ_DIR/$SAMPLE \
        --sample=$SAMPLE \
        --nosecondary \
        --localcores=$CORES \
        --localmem=$MEMORY \
        --create-bam=true

    echo "Finished Sample: $SAMPLE"
    echo "BAM file saved at:"
    echo "$PROJECT_DIR/${SAMPLE}_bamrun/outs/possorted_genome_bam.bam"
    echo "----------------------------------------"
done

echo "All samples completed."

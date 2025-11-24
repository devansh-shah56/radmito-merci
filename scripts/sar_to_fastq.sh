

##############################################################################
# SRA to FASTQ Conversion Pipeline
# Generates 3 FASTQ files per SRA dataset (paired reads + unpaired)
##############################################################################

# Set up error handling
set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure



# Define your SRA accession numbers
SRA_IDS=(
    "SRR15339529"
    "SRR15339530"
    "SRR15339531"
    "SRR15339532"
)

# Define directories
BASE_DIR="${HOME}/sra_analysis"
SRA_DIR="${BASE_DIR}/sra_files"
FASTQ_DIR="${BASE_DIR}/fastq_files"

# Number of threads for parallel-fastq-dump (adjust based on your CPU)
THREADS=4


echo "=========================================="
echo "SRA to FASTQ Conversion Pipeline"
echo "=========================================="
echo "Start time: $(date)"
echo ""

# Activate conda environment
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate sra_to_fastq

# Create directories if they don't exist
mkdir -p "${SRA_DIR}"
mkdir -p "${FASTQ_DIR}"

# Verify tools are available
echo "Checking required tools..."
command -v prefetch >/dev/null 2>&1 || { echo "ERROR: prefetch not found"; exit 1; }
command -v fastq-dump >/dev/null 2>&1 || { echo "ERROR: fastq-dump not found"; exit 1; }
echo "✓ All tools available"
echo ""


# MAIN PIPELINE


for SRA_ID in "${SRA_IDS[@]}"; do
    echo "=========================================="
    echo "Processing: ${SRA_ID}"
    echo "=========================================="
    
    # Step 1: Download SRA file
    echo "[1/3] Downloading SRA file..."
    cd "${SRA_DIR}"
    
    if [ -f "${SRA_ID}" ]; then
        echo "✓ SRA file already exists: ${SRA_ID}"
    else
        # Try direct download first (faster, bypasses SSL issues)
        echo "Attempting direct download from AWS S3..."
        if wget -q --show-progress "https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRA_ID}/${SRA_ID}"; then
            echo "✓ Direct download successful"
        else
            echo "Direct download failed, trying prefetch..."
            prefetch "${SRA_ID}" --max-size 50G
            mv "${HOME}/ncbi/public/sra/${SRA_ID}.sra" "${SRA_DIR}/${SRA_ID}" 2>/dev/null || true
        fi
    fi
    
    # Verify SRA file exists
    if [ ! -f "${SRA_DIR}/${SRA_ID}" ]; then
        echo "ERROR: SRA file not found for ${SRA_ID}"
        continue
    fi
    
    FILE_SIZE=$(ls -lh "${SRA_DIR}/${SRA_ID}" | awk '{print $5}')
    echo "✓ SRA file size: ${FILE_SIZE}"
    echo ""
    
    # Step 2: Convert SRA to FASTQ with --split-3
    echo "[2/3] Converting SRA to FASTQ..."
    cd "${FASTQ_DIR}"
    
    # Check if FASTQ files already exist
    if [ -f "${SRA_ID}_1.fastq.gz" ] && [ -f "${SRA_ID}_2.fastq.gz" ]; then
        echo "✓ FASTQ files already exist for ${SRA_ID}"
    else
        echo "Running fastq-dump with --split-3 option..."
        fastq-dump --split-3 \
                   --gzip \
                   --skip-technical \
                   --readids \
                   --read-filter pass \
                   --dumpbase \
                   --clip \
                   "${SRA_DIR}/${SRA_ID}"
        
        echo "✓ Conversion complete"
    fi
    echo ""
    
    # Step 3: Verify output files
    echo "[3/3] Verifying output files..."
    
    EXPECTED_FILES=(
        "${SRA_ID}_1.fastq.gz"   # Forward reads (Read 1)
        "${SRA_ID}_2.fastq.gz"   # Reverse reads (Read 2)
        "${SRA_ID}.fastq.gz"     # Unpaired reads (if any)
    )
    
    echo "Checking for output files:"
    for FILE in "${EXPECTED_FILES[@]}"; do
        if [ -f "${FASTQ_DIR}/${FILE}" ]; then
            SIZE=$(ls -lh "${FASTQ_DIR}/${FILE}" | awk '{print $5}')
            # Count reads (approximate)
            READS=$(zcat "${FASTQ_DIR}/${FILE}" | wc -l)
            READ_COUNT=$((READS / 4))
            echo "  ✓ ${FILE} (${SIZE}, ~${READ_COUNT} reads)"
        else
            echo "  ✗ ${FILE} (not found - may indicate no unpaired reads)"
        fi
    done
    
   
    echo ""
    echo "Sample of first read from ${SRA_ID}_1.fastq.gz:"
    zcat "${FASTQ_DIR}/${SRA_ID}_1.fastq.gz" | head -4
    echo ""
    
    echo "✓ ${SRA_ID} processing complete"
    echo ""
done



echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Output directory: ${FASTQ_DIR}"
echo ""
echo "Summary of generated files:"
ls -lh "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null || echo "No FASTQ files found"
echo ""
echo "Next steps:"
echo "1. Check file sizes to ensure complete downloads"
echo "2. Run FastQC for quality control"
echo "3. Proceed with alignment using BWA"
echo ""
echo "=========================================="
#!/bin/bash
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0:20:00
#SBATCH -J sam_covered_bases
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/samtools/samtools-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/samtools/samtools-%j.err

set -e  # Exit on error

#####################################
# Parse arguments
#####################################
if [ $# -lt 3 ] || [ $# -gt 4 ]; then
    echo "ERROR: Incorrect number of arguments"
    echo "Usage: $0 <sampleID> <bamfile> <min_depth> [max_depth]"
    echo ""
    echo "Arguments:"
    echo "  sampleID    - Sample identifier (e.g., P20702_102)"
    echo "  bamfile     - Path to BAM file"
    echo "  min_depth   - Minimum depth threshold (e.g., 100)"
    echo "  max_depth   - Maximum depth threshold (optional, e.g., 200)"
    echo ""
    echo "Examples:"
    echo "  $0 P20702_102 sample.bam 100 200      # Depth range 100-200"
    echo "  $0 P20702_102 sample.bam 300          # Depth ≥300 (no upper limit)"
    exit 1
fi

SAMPLE_ID=$1
BAM_FILE=$2
MIN_DEPTH=$3
MAX_DEPTH=${4:-}  # Optional 4th argument, empty string if not provided

#####################################
# Paths
#####################################
CODEDIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/Friberg_2025_TranscriptError_eriks_copy/TE
REFINDEX=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/genome.fa.fai
GENES=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/genes_dm6.bed
INTERGENIC=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/intergenic_dm6.bed

# Determine output directory based on depth range
if [ -z "${MAX_DEPTH}" ]; then
    DEPTH_LABEL="${MIN_DEPTH}"
    DEPTH_DESC="≥${MIN_DEPTH}x"
else
    DEPTH_LABEL="${MIN_DEPTH}-${MAX_DEPTH}"
    DEPTH_DESC="${MIN_DEPTH}-${MAX_DEPTH}x"
fi

DEPTH_DIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/Depth_dir/${DEPTH_LABEL}_Q73

#####################################
# Load modules
#####################################
module load PDC
module load miniconda3/25.3.1-1-cpeGNU-24.11
module load samtools

#####################################
# Print job info
#####################################
echo "========================================"
echo "DEPTH CALCULATION - SLURM Job ${SLURM_JOB_ID}"
echo "========================================"
echo "Started at: $(date)"
echo ""
echo "=== Configuration ==="
echo "Sample ID: ${SAMPLE_ID}"
echo "BAM file: ${BAM_FILE}"
echo "Minimum depth: ${MIN_DEPTH}"
if [ -z "${MAX_DEPTH}" ]; then
    echo "Maximum depth: none (no upper limit)"
else
    echo "Maximum depth: ${MAX_DEPTH}"
fi
echo "Depth range: ${DEPTH_DESC}"
echo "Output directory: ${DEPTH_DIR}"
echo "====================="
echo ""

#####################################
# Validate inputs
#####################################
echo "=== Validating Inputs ==="

if [ ! -f "${BAM_FILE}" ]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi
echo "✓ BAM file exists"

if [ ! -f "${GENES}" ]; then
    echo "ERROR: Genes BED file not found: ${GENES}"
    exit 1
fi
echo "✓ Genes BED file exists"

if [ ! -f "${INTERGENIC}" ]; then
    echo "ERROR: Intergenic BED file not found: ${INTERGENIC}"
    exit 1
fi
echo "✓ Intergenic BED file exists"

if [ ! -d "${DEPTH_DIR}" ]; then
    echo "ERROR: Output directory does not exist: ${DEPTH_DIR}"
    echo "       Run covered_bases_cohort.sh first to create it!"
    exit 1
fi
echo "✓ Output directory exists"

echo ""

#####################################
# Calculate depth with samtools
#####################################
echo "=== Running samtools depth ==="
echo "Quality filters: -Q 1 -q 73"
echo "NOTE: samtools outputs all depths; filtering will be done in Python"
echo ""

GENES_DEPTH="${DEPTH_DIR}/${SAMPLE_ID}.genes.depth"
INTERGENIC_DEPTH="${DEPTH_DIR}/${SAMPLE_ID}.intergenic.depth"

echo "Processing genes..."
samtools depth -Q 1 -q 73 -b ${GENES} ${BAM_FILE} > ${GENES_DEPTH}
if [ $? -eq 0 ]; then
    GENES_SIZE=$(wc -l < ${GENES_DEPTH})
    echo "✓ Genes depth calculated (${GENES_SIZE} positions)"
else
    echo "ERROR: samtools depth failed for genes"
    exit 1
fi

echo ""
echo "Processing intergenic regions..."
samtools depth -Q 1 -q 73 -b ${INTERGENIC} ${BAM_FILE} > ${INTERGENIC_DEPTH}
if [ $? -eq 0 ]; then
    INTERGENIC_SIZE=$(wc -l < ${INTERGENIC_DEPTH})
    echo "✓ Intergenic depth calculated (${INTERGENIC_SIZE} positions)"
else
    echo "ERROR: samtools depth failed for intergenic regions"
    exit 1
fi

echo ""

#####################################
# Calculate covered bases
#####################################
echo "=== Calculating Covered Bases ==="
echo "Script: depth2covered_comp.py"
echo "Minimum depth: ${MIN_DEPTH}"
if [ -z "${MAX_DEPTH}" ]; then
    echo "Maximum depth: none (no upper limit)"
else
    echo "Maximum depth: ${MAX_DEPTH}"
fi
echo ""

COVERED_BASES_FILE="${DEPTH_DIR}/${SAMPLE_ID}.covered_bases.txt"

# Build Python command with optional max_depth
PYTHON_CMD="python ${CODEDIR}/depth2covered_comp.py \
    --sampleid ${SAMPLE_ID} \
    --DG ${GENES_DEPTH} \
    --DIG ${INTERGENIC_DEPTH} \
    --O ${COVERED_BASES_FILE} \
    --depth_limit ${MIN_DEPTH}"

# Add max_depth only if provided
if [ -n "${MAX_DEPTH}" ]; then
    PYTHON_CMD="${PYTHON_CMD} --max_depth ${MAX_DEPTH}"
fi

# Execute Python script
eval ${PYTHON_CMD}

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Covered bases calculated"
    echo ""
    echo "=== Output ==="
    cat ${COVERED_BASES_FILE}
else
    echo "ERROR: depth2covered_comp.py failed"
    exit 1
fi

#####################################
# Job completion
#####################################
echo ""
echo "========================================"
echo "Job completed successfully!"
echo "Finished at: $(date)"
echo ""
echo "Output files:"
echo "  - ${COVERED_BASES_FILE}"
echo "  - ${GENES_DEPTH}"
echo "  - ${INTERGENIC_DEPTH}"
echo "========================================"
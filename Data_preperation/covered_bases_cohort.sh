#!/bin/bash
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 3:30:00
#SBATCH -J covered_bases_wrapper

# Output and error files
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/covered_bases_cohort-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/covered_bases_cohort-%j.err

set -e

#####################################
# CONFIGURATION - CHANGE THESE VALUES
#####################################
MIN_DEPTH=500         # Minimum depth threshold (e.g., 100, 200, 300)
#MAX_DEPTH=500         # Maximum depth threshold (e.g., 200, 300, 400) - set to empty string for no max
                      # Example intervals: 100-200, 200-300, 300-400, 400-500

#####################################
# Paths
#####################################
CODEDIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/Friberg_2025_TranscriptError_eriks_copy/TE
BAMDIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/bam

# Create descriptive folder name based on depth range
if [ -z "${MAX_DEPTH}" ] || [ "${MAX_DEPTH}" == "none" ]; then
    DEPTH_LABEL="${MIN_DEPTH}"
    DEPTH_DESC="≥${MIN_DEPTH}x"
else
    DEPTH_LABEL="${MIN_DEPTH}-${MAX_DEPTH}"
    DEPTH_DESC="${MIN_DEPTH}-${MAX_DEPTH}x"
fi

DEPTH_DIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/Depth_dir/${DEPTH_LABEL}_Q73

#####################################
# Create output directory
#####################################
echo "========================================"
echo "DEPTH CALCULATION - Wrapper Script"
echo "========================================"
echo "Started at: $(date)"
echo ""
echo "=== Configuration ==="
echo "Minimum depth: ${MIN_DEPTH}"
echo "Maximum depth: ${MAX_DEPTH:-none (no upper limit)}"
echo "Depth range: ${DEPTH_DESC}"
echo "BAM directory: ${BAMDIR}"
echo "Output directory: ${DEPTH_DIR}"
echo "====================="
echo ""

# Create depth directory if it doesn't exist
if [ ! -d "${DEPTH_DIR}" ]; then
    echo "Creating output directory: ${DEPTH_DIR}"
    mkdir -p ${DEPTH_DIR}
else
    echo "Output directory already exists: ${DEPTH_DIR}"
    echo "WARNING: Existing files may be overwritten!"
fi

echo ""

#####################################
# Count BAM files
#####################################
BAM_COUNT=$(ls ${BAMDIR}/P20702_*.Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)

if [ ${BAM_COUNT} -eq 0 ]; then
    echo "ERROR: No BAM files found in ${BAMDIR}"
    echo "Expected pattern: P20702_*.Aligned.sortedByCoord.out.bam"
    exit 1
fi

echo "Found ${BAM_COUNT} BAM files to process"
echo ""

#####################################
# Submit jobs for each BAM file
#####################################
echo "=== Submitting Jobs ==="
echo ""

SUBMITTED_COUNT=0
SKIPPED_COUNT=0

for bamfile in ${BAMDIR}/P20702_*.Aligned.sortedByCoord.out.bam; do 
    sampleID=$(basename $bamfile | awk -F "." '{ print $1 }')  # P20702_102
    
    # Check if output already exists
    OUTPUT_FILE="${DEPTH_DIR}/${sampleID}.covered_bases.txt"
    if [ -f "${OUTPUT_FILE}" ]; then
        echo "[${SUBMITTED_COUNT}] SKIPPING ${sampleID} (output already exists)"
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        continue
    fi
    
    echo "[${SUBMITTED_COUNT}] Submitting ${sampleID}"
    echo "    BAM: ${bamfile}"
    
    # Submit job with both min and max depth parameters
    if [ -z "${MAX_DEPTH}" ] || [ "${MAX_DEPTH}" == "none" ]; then
        JOB_ID=$(sbatch --parsable ${CODEDIR}/covered_bases.sh ${sampleID} ${bamfile} ${MIN_DEPTH})
    else
        JOB_ID=$(sbatch --parsable ${CODEDIR}/covered_bases.sh ${sampleID} ${bamfile} ${MIN_DEPTH} ${MAX_DEPTH})
    fi
    
    echo "    Job ID: ${JOB_ID}"
    echo ""
    
    SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
    sleep 30
done

echo "========================================"
echo "Summary:"
echo "  - Total BAM files: ${BAM_COUNT}"
echo "  - Jobs submitted: ${SUBMITTED_COUNT}"
echo "  - Files skipped: ${SKIPPED_COUNT}"
echo "========================================"
echo ""
echo "Monitor jobs with: squeue -u $USER"
echo "Check output directory: ${DEPTH_DIR}"
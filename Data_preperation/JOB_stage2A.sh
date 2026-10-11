#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH --mem=24G
#SBATCH -J 2A_te_stage2a
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage2a/s2a-%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage2a/s2a-%A_%a.err
#
# Stage 2A - per-sample depth (array job, one task per sample)
# ============================================================
# Task N processes line N of the sample list: runs samtools depth on that
# sample's STAR BAM (MAPQ >= 1, base quality >= 73) and writes the depth
# denominators.
# Submitted by submit_pipeline.sh once Stages 0A and 0C have finished.
#
# Outputs (in OUT_DIR): <sample>.aggregates.json, <sample>.depth.npz,
#                       <sample>.callable.bed
#
set -euo pipefail
SAMPLES_LIST="${1:?usage: sbatch --array=1-N JOB_stage2A.sh samples.txt}"

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
CANON_DIR="${BASE_DIR}/refgenome/canonical"
FAI="${CANON_DIR}/genome.canonical.fa.fai"
ANALYSIS_BED="${CANON_DIR}/analysis.bed"
# CDS of the representative transcripts (Stage 0C)
CDS_BED="${CANON_DIR}/ns_sites.cds.bed"

BAM_DIR="${BASE_DIR}/transcription_errors/bam"
OUT_DIR="${BASE_DIR}/transcription_errors/stage2a/per_sample"
SCRIPT="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE/stage2a_depth_sample.py"

MIN_DEPTH=100
MAPQ=1
BASEQ=73
MITO=M

mkdir -p "${OUT_DIR}"

module load samtools
ml PDC
conda activate transcription_errors

SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_LIST}")
BAM="${BAM_DIR}/${SAMPLE_ID}.Aligned.sortedByCoord.out.bam"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE_ID}  bam=${BAM}"
[ -f "${BAM}" ] || { echo "ERROR: missing ${BAM}"; exit 1; }

python "${SCRIPT}" \
    --sample-id "${SAMPLE_ID}" \
    --bam "${BAM}" \
    --fai "${FAI}" \
    --analysis-bed "${ANALYSIS_BED}" \
    --cds-bed "${CDS_BED}" \
    --out-prefix "${OUT_DIR}/${SAMPLE_ID}" \
    --min-depth "${MIN_DEPTH}" --mapq "${MAPQ}" --baseq "${BASEQ}" \
    --mito-contig "${MITO}"

echo "Done: ${SAMPLE_ID}"
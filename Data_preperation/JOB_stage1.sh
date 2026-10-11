#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0:45:00
#SBATCH --mem=8G
#SBATCH -J 1_te_stage1
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage1/stage1-%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage1/stage1-%A_%a.err
#
# Stage 1 - per-sample SNV extraction (array job, one task per sample)
# ====================================================================
# Task N processes line N of the sample list: filters that sample's freebayes
# VCF to primitive SNVs passing the mean alt-quality threshold.
# Submitted by submit_pipeline.sh with --array=1-<number of samples>.
#
# Outputs (in OUT_DIR): <sample>.snv.parquet, <sample>.stage1.json
#
set -euo pipefail

SAMPLES_LIST="${1:?usage: sbatch --array=1-N JOB_stage1.sh samples.txt}"

#####################################
# Paths
#####################################
VCF_DIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/concatenated_freebayes_clean
OUT_DIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/stage1
# Canonical reference from Stage 0A, used to check the VCF REF bases.
FASTA=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/refgenome/canonical/genome.canonical.fa
SCRIPT=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/Friberg_2025_TranscriptError_eriks_copy/TE/stage1_extract_snvs.py

CHROMS="2L,2R,3L,3R,4,X,Y,M"
MITO=M
QA=74          # minimum mean alt base quality (QA/AO)

mkdir -p "${OUT_DIR}"

#####################################
# Environment
#####################################
ml PDC
conda activate transcription_errors

#####################################
# Pick this task's sample
#####################################
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_LIST}")
VCF="${VCF_DIR}/${SAMPLE_ID}.concatenated.vcf"

echo "Task ${SLURM_ARRAY_TASK_ID}: sample=${SAMPLE_ID}"
echo "VCF=${VCF}"
[ -f "${VCF}" ] || { echo "ERROR: missing ${VCF}"; exit 1; }

python "${SCRIPT}" \
    --sample-id "${SAMPLE_ID}" \
    --vcf "${VCF}" \
    --out-prefix "${OUT_DIR}/${SAMPLE_ID}" \
    --chromosomes "${CHROMS}" \
    --mito-contig "${MITO}" \
    --qa-threshold "${QA}" \
    --fasta "${FASTA}"

echo "Done: ${SAMPLE_ID}"
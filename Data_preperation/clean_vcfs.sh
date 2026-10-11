#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0:20:00
#SBATCH --mem=2G
#SBATCH -J te_vcf_clean
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage1/clean-%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage1/clean-%A_%a.err
set -euo pipefail

SAMPLES_LIST="${1:?usage: sbatch --array=1-N clean_vcfs.sh samples.txt}"

RAW_DIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/concatenated_freebayes
CLEAN_DIR=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/concatenated_freebayes_clean

mkdir -p "$CLEAN_DIR"

SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
RAW="${RAW_DIR}/${SAMPLE_ID}.concatenated.vcf"
CLEAN="${CLEAN_DIR}/${SAMPLE_ID}.concatenated.vcf"

[ -f "$RAW" ] || { echo "ERROR: missing $RAW"; exit 1; }

awk 'BEGIN{d=0} /^#/{ if(!d) print; next } { d=1; print }' "$RAW" > "$CLEAN"

n=$(grep -c '^#CHROM' "$CLEAN")
[ "$n" -eq 1 ] || { echo "ERROR: $n header blocks left in $CLEAN"; exit 1; }
echo "cleaned ${SAMPLE_ID}: $(grep -vc '^#' "$CLEAN") records -> $CLEAN"
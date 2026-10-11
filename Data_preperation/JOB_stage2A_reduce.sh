#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH --mem=64G
#SBATCH -J 2areduce_te_stage2a_reduce
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage2a/reduce-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage2a/reduce-%j.err
#
# Stage 2A reduce - cohort per-base depth
# =======================================
# Runs once after all Stage 2A tasks have finished. Sums the per-sample depth
# into per-base cohort arrays (summed depth and number of samples with depth
# >= MIN_DEPTH), used by Stage 3.
#
# Outputs (in OUT_DIR): cohort_sum_depth.<contig>.npy,
#                       cohort_n_callable.<contig>.npy
#
set -euo pipefail
SAMPLES_LIST="${1:?usage: sbatch JOB_stage2A_reduce.sh samples.txt}"

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
CANON_DIR="${BASE_DIR}/refgenome/canonical"
FAI="${CANON_DIR}/genome.canonical.fa.fai"
PER_SAMPLE="${BASE_DIR}/transcription_errors/stage2a/per_sample"
OUT_DIR="${BASE_DIR}/transcription_errors/stage2a/cohort"
SCRIPT="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE/stage2a_reduce.py"

ml PDC
conda activate transcription_errors

python "${SCRIPT}" \
    --depth-dir "${PER_SAMPLE}" \
    --fai "${FAI}" \
    --samples-list "${SAMPLES_LIST}" \
    --out-dir "${OUT_DIR}" \
    --min-depth 100            # must equal MIN_DEPTH in JOB_stage2A.sh

echo "Reduce done -> ${OUT_DIR}"
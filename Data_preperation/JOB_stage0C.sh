#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --mem=32G
#SBATCH -J te_stage0c
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage0/0C_ns-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage0/0C_ns-%j.err
#
# Stage 0C - Nei-Gojobori N/S site table
# ======================================
# Computes synonymous and nonsynonymous site counts for every CDS position of
# the representative (longest-CDS) transcript of each gene, from the canonical
# reference built in Stage 0A. See stage0c_build_ns_sites.py for the method.
#
# Outputs (in refgenome/canonical/):
#   ns_sites.parquet           per-position n_site/s_site  (Stages 2B, 3)
#   ns_sites.channels.parquet  per-change weights          (Stage 3)
#   ns_sites.cds.bed           the CDS region              (Stages 2A, 3)
#
set -euo pipefail

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
CANON_DIR="${BASE_DIR}/refgenome/canonical"
SCRIPT_DIR="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE"

ml PDC
conda activate transcription_errors

python "${SCRIPT_DIR}/stage0c_build_ns_sites.py" \
    --gtf "${CANON_DIR}/annotation.canonical.gtf" \
    --fasta "${CANON_DIR}/genome.canonical.fa" \
    --out "${CANON_DIR}/ns_sites.parquet" \
    --mito-contig M
echo "Stage 0C done -> ${CANON_DIR}/ns_sites.parquet"
#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 2
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -J 2B_te_stage2b
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage2b/s2b-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage2b/s2b-%j.err
#
# Stage 2B - annotation
# =====================
#   1) collect the distinct SNVs of all samples into one VCF
#   2) annotate it with snpEff, using the custom database from Stage 0B
#   3) add consequence (representative transcripts only), transcribed strand
#      and spectrum channel to every sample's Stage 1 table
#
# Outputs: <OUT_DIR>/per_sample/<sample>.annotated.parquet
#
set -euo pipefail
SAMPLES_LIST="${1:?usage: sbatch JOB_stage2B.sh samples.txt}"

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
CANON_DIR="${BASE_DIR}/refgenome/canonical"
FAI="${CANON_DIR}/genome.canonical.fa.fai"
GTF="${CANON_DIR}/annotation.canonical.gtf"

STAGE1_DIR="${BASE_DIR}/transcription_errors/stage1"
NS_TABLE="${BASE_DIR}/refgenome/canonical/ns_sites.parquet"     # from Stage 0C
OUT_DIR="${BASE_DIR}/transcription_errors/stage2b"
WORK="${OUT_DIR}/work"
SCRIPT_DIR="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE"

# snpEff custom database (Stage 0B)
SNPEFF_CONFIG="${BASE_DIR}/refgenome/snpeff_custom.config"
SNPEFF_DATADIR="${BASE_DIR}/refgenome/snpeff_data"
GENOME="BDGP6.46.canonical"

mkdir -p "${WORK}" "${OUT_DIR}/per_sample"

ml PDC
conda activate transcription_errors

UNION_VCF="${WORK}/union_snvs.vcf"
ANN_VCF="${WORK}/union_snvs.ann.vcf"

echo "=== (1) build union VCF ==="
python "${SCRIPT_DIR}/stage2b_build_union_vcf.py" \
    --stage1-dir "${STAGE1_DIR}" --fai "${FAI}" --out-vcf "${UNION_VCF}" \
    --samples-list "${SAMPLES_LIST}"

echo "=== (2) snpEff annotate ==="
snpEff -config "${SNPEFF_CONFIG}" -dataDir "${SNPEFF_DATADIR}" \
    -noStats "${GENOME}" "${UNION_VCF}" > "${ANN_VCF}"

echo "=== (3) finalize (parse ANN + spectrum + enrich samples) ==="
python "${SCRIPT_DIR}/stage2b_finalize.py" \
    --annotated-vcf "${ANN_VCF}" \
    --stage1-dir "${STAGE1_DIR}" \
    --gtf "${GTF}" \
    --ns-table "${NS_TABLE}" \
    --fai "${FAI}" \
    --samples-list "${SAMPLES_LIST}" \
    --out-dir "${OUT_DIR}/per_sample"

echo "Stage 2B done -> ${OUT_DIR}/per_sample"
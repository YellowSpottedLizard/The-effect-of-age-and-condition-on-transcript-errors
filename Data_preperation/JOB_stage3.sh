#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH --mem=64G
#SBATCH -J 3_te_stage3
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage3/s3-%A_%a.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage3/s3-%A_%a.err
#
# Stage 3 - filtering and summary (array job, one task per SETTINGS entry)
# ========================================================================
# Each entry of SETTINGS below is one filter setting:
#
#     "label|filter|params_json"
#
#   label   output prefix: <OUT_DIR>/<label>.summary.tsv, .overview.tsv,
#           .overview_long.tsv, .params.json
#   filter  none         no germline disqualification (unfiltered baseline)
#           across_rate  disqualify a site-allele whose cohort-wide alt-read
#                        fraction (summed AO / summed depth) exceeds
#                        params_json "across_rate"
#
# submit_pipeline.sh counts the SETTINGS entries to size the array. Without an
# array index, all settings run one after another.
#
set -euo pipefail
SAMPLES_LIST="${1:?usage: sbatch [--array=1-N] JOB_stage3.sh samples.txt}"

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
CANON_DIR="${BASE_DIR}/refgenome/canonical"
TE="${BASE_DIR}/transcription_errors"
SCRIPT_DIR="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE"

STAGE2B_DIR="${TE}/stage2b/per_sample"
COHORT_DIR="${TE}/stage2a/cohort"
PER_SAMPLE_DIR="${TE}/stage2a/per_sample"
NS_TABLE="${CANON_DIR}/ns_sites.parquet"
NS_CHANNEL="${CANON_DIR}/ns_sites.channels.parquet"
FASTA="${CANON_DIR}/genome.canonical.fa"
FAI="${FASTA}.fai"
GENES_BED="${CANON_DIR}/genes.bed"          # gene spans: denominator for spec_*
CDS_BED="${CANON_DIR}/ns_sites.cds.bed"     # CDS: denominator for spec_*_cds
OUT_DIR="${TE}/stage3"
MIN_DEPTH=100        # must equal MIN_DEPTH in JOB_stage2A.sh
MIN_CALLABLE=1       # samples with depth >= MIN_DEPTH needed to analyse a position
MITO=M               # reported in separate '_mito' columns

mkdir -p "${OUT_DIR}"

ml PDC
set +u
conda activate transcription_errors
set -u

cd "${SCRIPT_DIR}"

# ============================ SETTINGS ====================================
SETTINGS=(
  "unfiltered|none|{}"
  "across_5pct|across_rate|{\"across_rate\":0.05}"
  "across_1pct|across_rate|{\"across_rate\":0.01}"
  "across_0p1pct|across_rate|{\"across_rate\":0.001}"
  "across_0p01pct|across_rate|{\"across_rate\":0.0001}"
)
# ==========================================================================

run_one() {
    local entry="$1"
    local label="${entry%%|*}"; local rest="${entry#*|}"
    local fname="${rest%%|*}"; local params="${rest#*|}"
    echo "=== ${label}  (${fname})  ${params} ==="
    python stage3_summary.py \
        --stage2b-dir "${STAGE2B_DIR}" --cohort-dir "${COHORT_DIR}" \
        --per-sample-dir "${PER_SAMPLE_DIR}" \
        --ns-table "${NS_TABLE}" --ns-channel "${NS_CHANNEL}" \
        --fasta "${FASTA}" --fai "${FAI}" --genes-bed "${GENES_BED}" --cds-bed "${CDS_BED}" \
        --samples-list "${SAMPLES_LIST}" \
        --filter "${fname}" --params "${params}" \
        --out-prefix "${OUT_DIR}/${label}" \
        --min-depth "${MIN_DEPTH}" --min-callable "${MIN_CALLABLE}" \
        --mito-contig "${MITO}"
}

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    idx=$((SLURM_ARRAY_TASK_ID - 1))
    if [ "${idx}" -ge "${#SETTINGS[@]}" ]; then
        echo "task ${SLURM_ARRAY_TASK_ID} > ${#SETTINGS[@]} settings; nothing to do"; exit 0
    fi
    run_one "${SETTINGS[$idx]}"
else
    echo "(no array index - running all ${#SETTINGS[@]} settings serially)"
    for entry in "${SETTINGS[@]}"; do run_one "${entry}"; done
fi

echo "Stage 3 task done -> ${OUT_DIR}"

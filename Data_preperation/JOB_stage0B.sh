#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 2
#SBATCH -t 0:40:00
#SBATCH --mem=16G
#SBATCH -J te_stage0_snpeff
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage0/snpeff-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage0/snpeff-%j.err
#
# Stage 0B - build a custom snpEff database from the canonical reference
# =======================================================================
#  
#
# Requires snpeff (bioconda) + a Java runtime (pulled in by the conda package).
#   conda install -c bioconda -c conda-forge snpeff
#
set -euo pipefail

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
CANON_DIR="${BASE_DIR}/refgenome/canonical"
SCRIPT_DIR="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE"
CANON_FASTA="${CANON_DIR}/genome.canonical.fa"
CANON_GTF="${CANON_DIR}/annotation.canonical.gtf"

GENOME="BDGP6.46.canonical"
DATA_DIR="${BASE_DIR}/refgenome/snpeff_data"     # snpEff -dataDir
CONFIG="${BASE_DIR}/refgenome/snpeff_custom.config"

ml PDC
conda activate transcription_errors

mkdir -p "${DATA_DIR}/${GENOME}"

# snpEff expects genes.gtf + sequences.fa under data/<genome>/
cp "${CANON_GTF}"   "${DATA_DIR}/${GENOME}/genes.gtf"
cp "${CANON_FASTA}" "${DATA_DIR}/${GENOME}/sequences.fa"

# ---- minimal custom config ----
# Mito (contig 'M') translated with the invertebrate mitochondrial code (table 5).
BUNDLED_CONFIG=$(ls ${CONDA_PREFIX}/share/snpeff*/snpEff.config 2>/dev/null | head -1)
[ -z "${BUNDLED_CONFIG}" ] && { echo "ERROR: no bundled snpEff.config under CONDA_PREFIX"; exit 1; }
cp "${BUNDLED_CONFIG}" "${CONFIG}"
cat >> "${CONFIG}" <<EOF

${GENOME}.genome : Drosophila_melanogaster_BDGP6.46_canonical
${GENOME}.M.codonTable : Invertebrate_Mitochondrial
EOF

echo "=== Building snpEff database: ${GENOME} ==="
BUILD_LOG="${DATA_DIR}/${GENOME}/build.log"
# -noCheckProtein / -noCheckCds: we do not supply Ensembl cds.fa / protein.fa.
# The two gates below replace what that validation would have caught.
snpEff build -gtf22 -v -noCheckProtein -noCheckCds \
    -config "${CONFIG}" -dataDir "${DATA_DIR}" "${GENOME}" 2>&1 | tee "${BUILD_LOG}"

#############################################################################
#  are the transcripts marked protein-coding?
#
# If the GTF has child rows before their transcript rows, snpEff reads every
# CDS but marks almost no transcript as coding. The build still exits 0, so
# the protein-coding fraction is read from the build summary instead.
#############################################################################
echo ""
echo "=== Gate 1: protein-coding content of the database ==="
N_TX=$(grep -m1 '^# Transcripts' "${BUILD_LOG}" | tr -dc '0-9')
N_PC=$(grep -m1 '^# Protein coding transcripts' "${BUILD_LOG}" | tr -dc '0-9')
N_PCG=$(grep -m1 '^# Protein coding genes' "${BUILD_LOG}" | tr -dc '0-9')
echo "  transcripts               : ${N_TX:-?}"
echo "  protein coding transcripts: ${N_PC:-?}"
echo "  protein coding genes      : ${N_PCG:-?}"
if [ -z "${N_TX}" ] || [ -z "${N_PC}" ] || [ "${N_TX}" -eq 0 ]; then
    echo "ERROR: could not read the build summary from ${BUILD_LOG}"; exit 1
fi
# Drosophila: ~30k of ~35k transcripts are protein coding. Anything under 20%
# means the marking step failed, not that the annotation is unusual.
if [ "$(( 100 * N_PC / N_TX ))" -lt 20 ]; then
    echo ""
    echo "ERROR: only ${N_PC} of ${N_TX} transcripts are protein-coding."
    echo "       snpEff read the CDS records but did not mark the transcripts."
    echo "       Check the GTF (Stage 0A should have normalised it) with:"
    echo "         python3 ${SCRIPT_DIR}/snpeff_gtf.py inspect --gtf ${CANON_GTF}"
    exit 1
fi
echo "  OK ($(( 100 * N_PC / N_TX ))% of transcripts are protein-coding)"

#############################################################################
#  does the database give coding annotations at real CDS positions?
#
# Builds SNVs from the middle of real CDS segments (REF read from the FASTA),
# annotates them, and requires at least 50% of nuclear probes to receive a
# coding consequence (synonymous, missense, stop gained/lost, start lost).
#############################################################################
echo ""
echo "=== Gate 2: coding-call probe on real CDS positions ==="
PROBE_VCF="${DATA_DIR}/${GENOME}/coding_probe.vcf"
PROBE_ANN="${DATA_DIR}/${GENOME}/coding_probe.ann.vcf"
python "${SCRIPT_DIR}/snpeff_gtf.py" make-testvcf \
    --gtf "${CANON_GTF}" --fasta "${CANON_FASTA}" --out "${PROBE_VCF}"
snpEff -config "${CONFIG}" -dataDir "${DATA_DIR}" -noStats "${GENOME}" \
    "${PROBE_VCF}" > "${PROBE_ANN}"
python "${SCRIPT_DIR}/snpeff_gtf.py" score --vcf "${PROBE_ANN}"


echo ""
echo "Database built: ${GENOME} in ${DATA_DIR}"

#!/bin/bash -l
#SBATCH -A naiss2025-22-1482
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 0:30:00
#SBATCH --mem=16G
#SBATCH -J te_stage0_ref
#SBATCH -o /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage0/ref-%j.out
#SBATCH -e /cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors/slurmout/stage0/ref-%j.err
#
# Stage 0A - canonical reference
# ==============================
# Builds the reference that every later stage uses, from the FASTA freebayes
# was run on (genome.fa) and its matching annotation (genes.gtf, the GTF the
# STAR index was built from).
#
# Steps:
#   1. Keep only 2L, 2R, 3L, 3R, 4, X, Y and the mitochondrion; rename headers
#      to bare contig names, with the mitochondrion renamed to 'M'.
#   2. Restrict and rename the GTF to the same contigs.
#   3. Normalise the GTF so snpEff parses coding transcripts correctly
#      (snpeff_gtf.py).
#   4. Verify that the GTF coordinates describe the FASTA by translating the
#      longest CDS of each gene (verify_ref_gtf.py). Stops on a mismatch.
#   5. Write genes.bed (gene spans, used in Stage 3) and analysis.bed (whole
#      canonical genome, used in Stage 2A).
#
# Outputs (in refgenome/canonical/):
#   genome.canonical.fa(.fai), annotation.canonical.gtf, reference_check.tsv,
#   genes.bed, analysis.bed
#
set -euo pipefail

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
REF_DIR="${BASE_DIR}/refgenome"
OUT_DIR="${REF_DIR}/canonical"
SCRIPT_DIR="${BASE_DIR}/Friberg_2025_TranscriptError_eriks_copy/TE"

# ---- INPUTS -----------------------------------------------------------------
# The FASTA used for variant calling.
IN_FASTA="${REF_DIR}/genome.fa"

# Annotation matching genome.fa. Can be overridden from the environment:
#   IN_GTF=/path/to/other.gtf sbatch JOB_stage0A.sh
IN_GTF="${IN_GTF:-${REF_DIR}/genes.gtf}"

# Nuclear contigs we keep, in canonical (output) naming.
KEEP_NUCLEAR=(2L 2R 3L 3R 4 X Y)
# The mito contig is auto-detected in both FASTA and GTF and renamed to this:
MITO_DST="M"

# Fail the job if fewer than this fraction of representative CDS are free of
# internal stop codons when read off the FASTA (catches a wrong assembly).
MIN_CDS_OK_FRAC=0.90
# Separate, lower floor for the start-codon and stop-codon tests (catches a
# small coordinate offset). Lower because real annotations contain some
# deliberately incomplete CDS (Ensembl cds_start_NF / cds_end_NF).
MIN_START_STOP_FRAC=0.80

ml PDC
module load bioinfo-tools 2>/dev/null || true
set +u
conda activate transcription_errors
set -u
module load samtools

mkdir -p "${OUT_DIR}"
OUT_FASTA="${OUT_DIR}/genome.canonical.fa"
OUT_GTF="${OUT_DIR}/annotation.canonical.gtf"

[ -f "${IN_FASTA}" ] || { echo "ERROR: missing ${IN_FASTA}"; exit 1; }
[ -f "${IN_GTF}" ]   || { echo "ERROR: missing ${IN_GTF}"; exit 1; }
[ -f "${IN_FASTA}.fai" ] || samtools faidx "${IN_FASTA}"

echo "=== INPUTS ==="
echo "  FASTA : ${IN_FASTA}"
echo "  GTF   : ${IN_GTF}"
echo ""

#############################################################################
# Resolve contig names in the FASTA and the GTF.
#
# We do NOT assume a naming scheme. For each canonical contig we look for the
# bare name first ('2L'), then a 'chr'-prefixed variant ('chr2L'). The mito is
# found by pattern (M / MT / chrM / *mitochond*). This makes the job work
# unchanged whether genome.fa uses Ensembl-style or UCSC-style names, and it
# fails loudly instead of silently dropping a chromosome.
#############################################################################
FA_NAMES=$(cut -f1 "${IN_FASTA}.fai")
GTF_NAMES=$(awk -F'\t' '!/^#/ && NF>=9 {print $1}' "${IN_GTF}" | sort -u)

find_contig() {   # $1 = canonical name, $2 = newline-separated name list
    local want="$1" list="$2"
    if printf '%s\n' "${list}" | grep -qx -- "${want}";    then printf '%s' "${want}";    return 0; fi
    if printf '%s\n' "${list}" | grep -qx -- "chr${want}"; then printf '%s' "chr${want}"; return 0; fi
    return 1
}

find_mito() {     # $1 = newline-separated name list
    printf '%s\n' "$1" | awk '
        tolower($0)=="m" || tolower($0)=="mt" || tolower($0)=="chrm" ||
        tolower($0)=="chrmt" || tolower($0) ~ /mitochond/ { print; exit }'
}

echo "=== Resolving contig names ==="
FA_SRC=(); GTF_SRC=()
for c in "${KEEP_NUCLEAR[@]}"; do
    s=$(find_contig "${c}" "${FA_NAMES}")  || { echo "ERROR: no contig for '${c}' in ${IN_FASTA}"; exit 1; }
    g=$(find_contig "${c}" "${GTF_NAMES}") || { echo "ERROR: no seqname for '${c}' in ${IN_GTF}"; exit 1; }
    FA_SRC+=("${s}"); GTF_SRC+=("${g}")
    printf '  %-4s FASTA=%-8s GTF=%s\n' "${c}" "${s}" "${g}"
done

FA_MITO=$(find_mito "${FA_NAMES}")
GTF_MITO=$(find_mito "${GTF_NAMES}")
[ -n "${FA_MITO}" ]  || { echo "ERROR: no mitochondrial contig found in ${IN_FASTA}"; exit 1; }
[ -n "${GTF_MITO}" ] || { echo "ERROR: no mitochondrial seqname found in ${IN_GTF}"; exit 1; }
printf '  %-4s FASTA=%-8s GTF=%s   -> renamed to %s\n' "mito" "${FA_MITO}" "${GTF_MITO}" "${MITO_DST}"

#############################################################################
# FASTA: extract the canonical contigs, normalising every header to the bare
# canonical name. samtools faidx keeps the full original header, so we rewrite
# the '>' line of each record explicitly.
#############################################################################
echo ""
echo "=== Extracting canonical contigs from FASTA ==="
: > "${OUT_FASTA}"
for i in "${!KEEP_NUCLEAR[@]}"; do
    dst="${KEEP_NUCLEAR[$i]}"; src="${FA_SRC[$i]}"
    echo "  ${src} -> ${dst}"
    samtools faidx "${IN_FASTA}" "${src}" | sed "1s/^>.*/>${dst}/" >> "${OUT_FASTA}"
done
echo "  ${FA_MITO} -> ${MITO_DST}"
samtools faidx "${IN_FASTA}" "${FA_MITO}" | sed "1s/^>.*/>${MITO_DST}/" >> "${OUT_FASTA}"

samtools faidx "${OUT_FASTA}"
echo "Canonical FASTA contigs (name, length):"
cut -f1,2 "${OUT_FASTA}.fai" | sed 's/^/  /'

#############################################################################
# GTF: keep only the canonical seqnames, rewriting column 1 to canonical names.
#############################################################################
echo ""
echo "=== Filtering + renaming GTF ==="
KEEP_MAP=""
for i in "${!KEEP_NUCLEAR[@]}"; do
    KEEP_MAP="${KEEP_MAP}${GTF_SRC[$i]}=${KEEP_NUCLEAR[$i]};"
done
KEEP_MAP="${KEEP_MAP}${GTF_MITO}=${MITO_DST};"

awk -v FS='\t' -v OFS='\t' -v mapstr="${KEEP_MAP}" '
    BEGIN {
        n = split(mapstr, pairs, ";")
        for (i = 1; i <= n; i++) {
            if (pairs[i] == "") continue
            split(pairs[i], kv, "=")
            m[kv[1]] = kv[2]
        }
    }
    /^#/ { print; next }
    NF >= 9 && ($1 in m) { $1 = m[$1]; print }
' "${IN_GTF}" > "${OUT_GTF}"

echo "GTF seqnames now:"
awk -F'\t' '!/^#/ && NF>=9 {print $1}' "${OUT_GTF}" | sort -u | sed 's/^/  /'
echo "GTF feature types now:"
awk -F'\t' '!/^#/ && NF>=9 {print $3}' "${OUT_GTF}" | sort | uniq -c | sort -k1,1nr | sed 's/^/  /'

N_CDS=$(awk -F'\t' '!/^#/ && NF>=9 && $3=="CDS"' "${OUT_GTF}" | wc -l)
echo "  CDS records: ${N_CDS}"
if [ "${N_CDS}" -eq 0 ]; then
    echo ""
    echo "ERROR: the canonical GTF has no CDS records."
    echo "       Stage 0c (N/S sites), the synonymous depth denominator and"
    echo "       snpEff's coding annotation all need CDS features with frames."
    echo "       Re-run with the Ensembl annotation instead:"
    echo "         IN_GTF=${REF_DIR}/Drosophila_melanogaster.BDGP6.46.112.gtf sbatch JOB_stage0A.sh"
    exit 1
fi

#############################################################################
# Normalise the GTF for snpEff.
#
# In genes.gtf, exon rows precede their transcript row and exon_ids are reused
# across sibling transcripts. snpEff reads the file in one pass, so it creates
# a stub transcript from the first exon it meets and rejects the real
# transcript row as a duplicate. The stub never gets transcript_biotype
# "protein_coding", so nearly all transcripts are treated as non-coding and
# coding SNVs are annotated as non_coding_transcript_exon_variant instead of
# synonymous/missense. No error is raised; the synonymous counts are simply
# near zero.
#
# snpeff_gtf.py fix puts parent rows first and makes exon_ids unique. It does
# not change any coordinates. 'inspect --assert-clean' stops the job if the
# problem is still present.
#############################################################################
echo ""
echo "=== Normalising GTF for snpEff (parent rows first, unique exon_ids) ==="
python "${SCRIPT_DIR}/snpeff_gtf.py" fix \
    --gtf "${OUT_GTF}" --out "${OUT_GTF}.norm"
mv "${OUT_GTF}.norm" "${OUT_GTF}"

# Refuse to continue on a GTF that would build a non-coding database.
python "${SCRIPT_DIR}/snpeff_gtf.py" inspect --gtf "${OUT_GTF}" --assert-clean

#############################################################################
# Verify that the GTF coordinates describe the FASTA sequence. For each gene's
# longest CDS: starts with a start codon, ends with a stop, no internal stops.
# Catches an annotation from a different assembly release, or an off-by-one
# coordinate convention error.
#############################################################################
echo ""
echo "=== Verifying GTF coordinates against the FASTA ==="
python "${SCRIPT_DIR}/verify_ref_gtf.py" \
    --fasta "${OUT_FASTA}" \
    --gtf "${OUT_GTF}" \
    --mito-contig "${MITO_DST}" \
    --min-ok-frac "${MIN_CDS_OK_FRAC}" \
    --min-start-stop-frac "${MIN_START_STOP_FRAC}" \
    --report "${OUT_DIR}/reference_check.tsv"

#############################################################################
# BEDs: gene spans (Stage 3 denominator for spec_*) and the whole canonical
# genome (Stage 2A analysis region).
#############################################################################
echo ""
echo "=== Building genes / analysis BEDs ==="
module load bedtools 2>/dev/null || true

GENES_BED="${OUT_DIR}/genes.bed"
ANALYSIS_BED="${OUT_DIR}/analysis.bed"

# Gene spans from 'gene' features, falling back to transcript, then exon, for
# GTFs without explicit gene rows.
GENE_FEAT="gene"
for feat in gene transcript exon; do
    if awk -v FS='\t' -v f="${feat}" '!/^#/ && NF>=9 && $3==f {found=1; exit} END{exit !found}' "${OUT_GTF}"; then
        GENE_FEAT="${feat}"; break
    fi
done
echo "  genic intervals derived from '${GENE_FEAT}' features"
awk -v FS='\t' -v OFS='\t' -v f="${GENE_FEAT}" '!/^#/ && NF>=9 && $3==f { print $1, $4-1, $5 }' "${OUT_GTF}" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "${GENES_BED}"

awk -v OFS='\t' '{ print $1, 0, $2 }' "${OUT_FASTA}.fai" | sort -k1,1 -k2,2n > "${ANALYSIS_BED}"

bases() { awk '{s += $3 - $2} END {printf "%d", s+0}' "$1"; }
echo "  genes      : $(wc -l < "${GENES_BED}") intervals, $(bases "${GENES_BED}") bp"
echo "  mito genic : $(awk -v m="${MITO_DST}" '$1==m {s += $3 - $2} END {printf "%d", s+0}' "${GENES_BED}") bp"

echo ""
echo "Done."
echo "  FASTA      : ${OUT_FASTA}"
echo "  GTF        : ${OUT_GTF}"
echo "  genes BED  : ${GENES_BED}"
echo "  analysis   : ${ANALYSIS_BED}"
echo "  ref check  : ${OUT_DIR}/reference_check.tsv"

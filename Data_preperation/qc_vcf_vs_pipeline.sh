#!/bin/bash -l
# =============================================================================
# qc_vcf_vs_pipeline.sh - reconcile the raw VCF against Stage 1 and Stage 3
# =============================================================================
# Run ON THE CLUSTER, where the VCFs live.
#
#   ./qc_vcf_vs_pipeline.sh samples.txt [N]
#
# N = how many samples to check (default 8). Every sample takes ~30 s, so use
# a handful routinely and the full 286 once before submission.
#
# THE CENTRAL IDENTITY
# --------------------
# For each compartment, the alleles the VCF offers and the alleles Stage 3
# accounts for must balance exactly:
#
#   awk_pass  -  rows_merged_by_trim  ==  total_alleles + excluded_alleles
#
#   awk_pass            TYPE=snp alleles with QA/AO >= 73, INCLUDING the
#                       REF-padded ones (Stage 1 left-aligns rather than
#                       dropping them, so they must be counted here too)
#   rows_merged_by_trim trimming can move two alleles onto the same
#                       (chrom,pos,ref,alt); those rows are collapsed and their
#                       AO summed, so the ALLELE count falls while the READ
#                       count is preserved
#   excluded_alleles    rows dropped for sitting outside the analysed region or
#                       outside this sample's own callable intervals
#
# If this balances, nothing was silently lost between the VCF and the summary.
# If it does not, the difference tells you which stage to look at.
#
# The read identity is stricter, because trimming preserves reads exactly:
#
#   awk_reads  ==  total_mutations + excluded_mutations
#
# =============================================================================
set -euo pipefail

SAMPLES="${1:?usage: ./qc_vcf_vs_pipeline.sh samples.txt [n_samples]}"
NCHECK="${2:-8}"

BASE=/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/transcription_errors
VCF_DIR="${BASE}/concatenated_freebayes_clean"
STAGE1_DIR="${BASE}/stage1"
STAGE3_DIR="${BASE}/stage3_Q74"
FILTER="f1_none_M"          # MUST be the unfiltered one: f8 removes alleles
QA_MIN=74

echo "VCF dir : ${VCF_DIR}"
echo "Filter  : ${FILTER}  (must be f1_none - an f8 file will not balance)"
echo ""

# ---- awk: count passing alleles and reads per compartment -------------------
count_vcf () {
  awk -F'\t' -v QMIN="${QA_MIN}" '$1!~/^#/ {
    n=split($8,info,";"); ao=""; qa=""; ty="";
    for(i=1;i<=n;i++){
      if(info[i]~/^AO=/)        ao=substr(info[i],4);
      else if(info[i]~/^QA=/)   qa=substr(info[i],4);
      else if(info[i]~/^TYPE=/) ty=substr(info[i],6);
    }
    na=split(ao,A,","); split(qa,Q,","); split(ty,T,","); split($5,ALT,",");
    c = ($1=="M") ? "mito" : "nuc";
    for(j=1;j<=na;j++){
      if(T[j]!="snp" || A[j]<=0) continue;
      if(Q[j]/A[j] < QMIN) continue;
      alt[c]++; rd[c]+=A[j];
      if(!(length($4)==1 && length(ALT[j])==1)) pad[c]++;
    }
  } END{
    printf "%d %d %d %d %d %d\n",
      alt["nuc"]+0, rd["nuc"]+0, pad["nuc"]+0,
      alt["mito"]+0, rd["mito"]+0, pad["mito"]+0
  }' "$1"
}

printf "%-14s %-6s %10s %10s %10s %10s  %s\n" \
  SAMPLE COMP VCF_ALLELE S3_ALLELE VCF_READS S3_READS STATUS

fail=0
head -n "${NCHECK}" "${SAMPLES}" | while read -r S; do
  [ -z "${S}" ] && continue
  VCF="${VCF_DIR}/${S}.concatenated.vcf"
  [ -f "${VCF}" ] || { echo "MISSING VCF ${S}"; continue; }
  set -- $(count_vcf "${VCF}")
  AN=$1; RN=$2; PN=$3; AM=$4; RM=$5; PM=$6

  MERGED=$(python3 -c "
import json,sys
d=json.load(open('${STAGE1_DIR}/${S}.stage1.json'))
print(d['counts'].get('rows_merged_by_trim',0))" 2>/dev/null || echo 0)

  S3="${STAGE3_DIR}/${FILTER}.summary.tsv" python3 - "$S" "$AN" "$RN" "$PN" "$AM" "$RM" "$PM" "$MERGED" <<'PY2'
import sys, csv, os
S, an, rn, pn, am, rm, pm, merged = sys.argv[1], *map(int, sys.argv[2:])
path = os.environ.get("S3", "")
row = None
with open(path) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["sample_id"] == S:
            row = r; break
if row is None:
    print(f"{S:<14} -- not in summary.tsv --"); sys.exit(0)

g = lambda k: int(float(row[k]))
for comp, sfx, va, vr, vp in (("nuc", "", an, rn, pn), ("mito", "_mito", am, rm, pm)):
    s_all = g(f"total_alleles{sfx}") + g(f"excluded_alleles{sfx}")
    s_rd  = g(f"total_mutations{sfx}") + g(f"excluded_mutations{sfx}")
    # merged rows are not attributed to a compartment in the log; allow the
    # allele identity to absorb up to `merged` on either side
    ok_a = abs((va - s_all)) <= merged
    ok_r = (vr == s_rd)
    status = "OK" if (ok_a and ok_r) else ("ALLELE" if not ok_a else "") + ("READS" if not ok_r else "")
    print(f"{S:<14} {comp:<6} {va:10d} {s_all:10d} {vr:10d} {s_rd:10d}  {status}"
          f"   (padded={vp}, merged<={merged})")
PY2
done
echo ""
echo "Interpretation:"
echo "  ALLELE mismatch  -> Stage 1 kept or dropped alleles the VCF scan disagrees about."
echo "                      Check alts_nonsnv / alts_trimmed / alts_below_qa in the"
echo "                      stage1.json; a large alts_nonsnv means trim_alleles is not"
echo "                      recovering the REF-padded records."
echo "  READS mismatch   -> AO is not being summed faithfully. Trimming preserves reads"
echo "                      exactly, so any difference here is a real bug."
echo "  Both OK          -> nothing was lost between the VCF and the summary table."
#!/usr/bin/env python3
"""
Stage 1 - per-sample VCF -> filtered SNV table
==============================================
For one sample's freebayes VCF:

  * Canonicalises contig names (contig_names.py) and keeps whitelisted
    contigs only. Per-contig counts of kept and skipped records go to the
    JSON log.
  * Processes each ALT allele separately. Padded alleles are trimmed to their
    minimal form (trim_alleles) and kept if they are SNVs (one base each and
    freebayes TYPE 'snp').
  * Keeps an allele if its mean alt base quality, QA / AO, is at least
    --qa-threshold (default 74). QA is freebayes' summed alt base quality, so
    QA / AO is the mean per read.
  * Records AO, RO, DP, QA, mean alt quality and VAF (AO/DP and AO/(AO+RO)).
  * Checks the (trimmed) REF base against the reference FASTA and records
    the result in ref_fasta_match. Mismatches are counted, not corrected or
    removed; the VCF REF is used downstream.
  * Merges rows that end up with the same (pos, ref, alt) after trimming.

Outputs:
  <out_prefix>.snv.parquet    one row per passing (pos, ref, alt)
  <out_prefix>.stage1.json    run log with counts

Consequence, gene, transcribed strand and spectrum channel are added in
Stage 2B, after snpEff annotation.

Requires: pysam, pandas, pyarrow.
"""
import argparse
import json
import os
import sys
from collections import Counter

import numpy as np
import pysam
import pandas as pd

import contig_names as CN

DEFAULT_CHROMS = "2L,2R,3L,3R,4,X,Y,M"


def trim_alleles(ref, alt):
    """
    Trim one padded freebayes allele to its minimal form.
    -> (pos_offset, ref, alt)

    freebayes pads REF to span the longest allele at a locus, so a single-base
    substitution can arrive as e.g. AAG>GAG when another allele at the same
    position is longer. Testing len(ref)==len(alt)==1 without trimming would
    discard such substitutions, and more often at high-depth loci (where more
    alleles are called), such as the mitochondrion.

    Strips the shared suffix (POS unchanged), then the shared prefix (POS
    moves right by the number of bases removed). Never trims below one base,
    so a multi-nucleotide substitution still fails the SNV test.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    off = 0
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        off += 1
    return off, ref, alt


def as_tuple(val, n):
    """Normalise a pysam INFO value to a length-n tuple aligned with ALTs."""
    if val is None:
        return (None,) * n
    if isinstance(val, (tuple, list)):
        return tuple(val)
    return (val,) * n  # scalar broadcast (shouldn't happen for Number=A fields)


def main():
    ap = argparse.ArgumentParser(
        description="Stage 1: filter a single sample's freebayes VCF to primitive SNVs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--sample-id", required=True,
                    help="Sample identifier, e.g. P20702_102")
    ap.add_argument("--vcf", required=True,
                    help="Path to the sample's *.concatenated.vcf (top folder)")
    ap.add_argument("--out-prefix", required=True,
                    help="Output path prefix (dir must exist), e.g. /.../stage1/P20702_102")
    ap.add_argument("--chromosomes", default=DEFAULT_CHROMS,
                    help="Comma-separated whitelist of canonical contig names to keep")
    ap.add_argument("--mito-contig", default=CN.CANONICAL_MITO,
                    help="Canonical name for the mitochondrial contig")
    ap.add_argument("--qa-threshold", type=float, default=74.0,
                    help="Minimum mean per-read alt quality (QA/AO) to keep an alt")
    ap.add_argument("--fasta", default=None,
                    help="Canonical reference FASTA (.fa with .fai) to flag REF mismatches")
    args = ap.parse_args()

    whitelist = {c.strip() for c in args.chromosomes.split(",") if c.strip()}

    fasta = pysam.FastaFile(args.fasta) if args.fasta else None
    fasta_contigs = set(fasta.references) if fasta else set()

    counts = Counter()
    per_contig = Counter()      # canonical contig -> passing alts
    skipped_contigs = Counter()  # raw contig name -> records dropped
    rows = []
    positions = set()  # (chrom, start0), for the unique-position count

    vcf = pysam.VariantFile(args.vcf)

    for rec in vcf:
        counts["records_total"] += 1
        # canonicalise before the whitelist test
        chrom = CN.canonicalise(rec.chrom, args.mito_contig)
        if chrom not in whitelist:
            counts["records_contig_skipped"] += 1
            skipped_contigs[rec.chrom] += 1
            continue
        per_contig[f"{chrom}_records"] += 1

        alts = rec.alts or ()
        n = len(alts)
        if n > 1:
            counts["records_multiallelic"] += 1

        ao = as_tuple(rec.info.get("AO"), n)
        qa = as_tuple(rec.info.get("QA"), n)
        vtype = as_tuple(rec.info.get("TYPE"), n)
        ro = rec.info.get("RO", 0) or 0
        dp = rec.info.get("DP", 0) or 0
        ref = rec.ref

        fcontig = chrom if chrom in fasta_contigs else None
        if fasta is not None and fcontig is None:
            counts["ref_fasta_contig_absent"] += 1

        for i, alt in enumerate(alts):
            counts["alts_total"] += 1

            # trim padding before the SNV test (see trim_alleles)
            off, tref, talt = trim_alleles(ref, alt)

            is_snv = (len(tref) == 1 and len(talt) == 1)
            if vtype[i] is not None:
                is_snv = is_snv and (vtype[i] == "snp")
            if not is_snv:
                counts["alts_nonsnv"] += 1
                continue
            counts["alts_snv"] += 1
            if off or len(ref) > 1:
                counts["alts_trimmed"] += 1
                per_contig[f"{chrom}_alts_trimmed"] += 1

            ao_i = ao[i] or 0
            qa_i = qa[i] or 0
            if ao_i <= 0:
                counts["alts_ao_zero"] += 1
                continue

            mean_alt_q = qa_i / ao_i
            if mean_alt_q < args.qa_threshold:
                counts["alts_below_qa"] += 1
                continue

            # coordinates AFTER trimming
            start0 = rec.start + off
            pos1 = rec.pos + off

            # REF check per ALT, since the trim offset is per allele; also
            # validates the trimming coordinates.
            ref_match = None
            if fasta is not None and fcontig is not None:
                fbase = fasta.fetch(fcontig, start0, start0 + 1).upper()
                ref_match = (fbase == tref.upper())
                if not ref_match:
                    counts["ref_fasta_mismatch"] += 1

            if ref_match is not None:
                counts["alts_ref_checked"] += 1
            counts["alts_pass"] += 1
            per_contig[f"{chrom}_alts_pass"] += 1
            positions.add((chrom, start0))

            rows.append({
                "sample_id": args.sample_id,
                "chrom": chrom,
                "pos": pos1,                    # 1-based POS, post-trim
                "ref": tref,
                "alt": talt,
                "type": vtype[i] if vtype[i] is not None else "snp",
                "ao": int(ao_i),
                "ro": int(ro),
                "dp": int(dp),
                "qa": float(qa_i),
                "mean_alt_qual": round(mean_alt_q, 3),
                "vaf_dp": round(ao_i / dp, 6) if dp > 0 else None,
                "vaf_ar": round(ao_i / (ao_i + ro), 6) if (ao_i + ro) > 0 else None,
                "ref_fasta_match": ref_match,
            })

    vcf.close()
    if fasta is not None:
        fasta.close()

    # ---- write outputs ----
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)) or ".", exist_ok=True)

    df = pd.DataFrame(rows)

    # Trimming can move an allele onto a position that already carries the same
    # (ref, alt) from another record (e.g. 'M 1541 GA>GT' becomes '1542 A>T').
    # Overlapping records may count the same reads, so the merged row takes
    # the record with the highest AO (and its QA); RO and DP take the maximum.
    # This avoids counting any read twice.
    counts["rows_before_merge"] = len(df)
    if len(df):
        key = ["sample_id", "chrom", "pos", "ref", "alt"]
        n_before = len(df)
        df = (df.sort_values(key + ["ao"], ascending=[True] * len(key) + [False])
                .groupby(key, as_index=False, sort=False)
                .agg({"type": "first", "ao": "first", "qa": "first",
                      "ro": "max", "dp": "max", "ref_fasta_match": "first"}))
        df["mean_alt_qual"] = (df["qa"] / df["ao"]).round(3)
        df["vaf_dp"] = np.where(df["dp"] > 0, (df["ao"] / df["dp"]).round(6), np.nan)
        df["vaf_ar"] = np.where((df["ao"] + df["ro"]) > 0,
                                (df["ao"] / (df["ao"] + df["ro"])).round(6), np.nan)
        counts["rows_merged_by_trim"] = n_before - len(df)

    parquet_path = f"{args.out_prefix}.snv.parquet"
    df.to_parquet(parquet_path, index=False)

    counts["positions_unique"] = len(positions)
    counts["rows_written"] = len(df)
    mito_pass = per_contig.get(f"{args.mito_contig}_alts_pass", 0)
    log = {
        "sample_id": args.sample_id,
        "vcf": args.vcf,
        "qa_threshold": args.qa_threshold,
        "chromosomes": sorted(whitelist),
        "mito_contig": args.mito_contig,
        "fasta_checked": fasta is not None,
        "counts": dict(counts),
        "per_contig": dict(per_contig),
        "skipped_contigs": dict(skipped_contigs.most_common(50)),
        "alts_pass_mito": mito_pass,
        "alts_pass_nuclear": counts.get("alts_pass", 0) - mito_pass,
    }
    if counts.get("alts_ref_checked"):
        log["ref_mismatch_pct"] = round(
            100.0 * counts.get("ref_fasta_mismatch", 0) / counts["alts_ref_checked"], 4)

    json_path = f"{args.out_prefix}.stage1.json"
    with open(json_path, "w") as fh:
        json.dump(log, fh, indent=2)

    # summary to stderr (SLURM .err)
    print(f"[STAGE1] sample={args.sample_id}", file=sys.stderr)
    for k in ("records_total", "records_contig_skipped", "records_multiallelic",
              "alts_total", "alts_snv", "alts_nonsnv", "alts_trimmed",
              "alts_below_qa", "alts_pass", "rows_merged_by_trim",
              "positions_unique", "ref_fasta_mismatch"):
        if k in counts:
            print(f"[STAGE1]   {k} = {counts[k]}", file=sys.stderr)
    print(f"[STAGE1]   alts_pass_nuclear = {log['alts_pass_nuclear']}", file=sys.stderr)
    print(f"[STAGE1]   alts_pass_mito ({args.mito_contig}) = {mito_pass}", file=sys.stderr)
    if mito_pass == 0:
        print(f"[STAGE1]   WARNING: no passing variants on the mitochondrion. "
              f"If that is unexpected, check the contig names actually present "
              f"in the VCF against the whitelist {sorted(whitelist)}.",
              file=sys.stderr)
    if skipped_contigs:
        top = ", ".join(f"{k}={v}" for k, v in skipped_contigs.most_common(5))
        print(f"[STAGE1]   skipped contigs (top 5): {top}", file=sys.stderr)
    print(f"[STAGE1] wrote {parquet_path} ({len(df)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
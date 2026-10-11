#!/usr/bin/env python3
"""
Stage 2B (part 2) - finalize annotation
=======================================
Takes the snpEff-annotated union VCF and, for every sample's Stage 1 table:

  1. Adds the consequence of each SNV on the representative transcript of its
     gene (the longest-CDS transcripts used for N and S in Stage 0C). If
     several representative transcripts carry a coding effect, the most severe
     is taken (stop_gained > stop_lost > start_lost > missense > synonymous).
     SNVs with no coding effect on a representative transcript are labelled
     'other', even if snpEff reports a coding effect on another isoform, so
     that synonymous and missense counts cover the same positions as N and S.
  2. Adds the transcribed strand, from the GTF gene (or, if absent,
     transcript) spans: '+', '-', or none where no feature or features on both
     strands cover the base.
  3. Adds the spectrum channel: ref>alt in transcript orientation (complemented
     for '-' strand), T written as U. No channel where the strand is none.

Output: <out-dir>/<sample>.annotated.parquet (Stage 1 columns plus consequence,
gene_id, transcript_id, effect, impact, gene_strand, spectrum_channel).

Stops if almost no annotation matches a representative transcript (transcript
ID mismatch), if no strand could be derived from the GTF, or if too few
nuclear SNVs get a coding consequence (snpEff treating transcripts as
non-coding).

Requires: pysam, pandas, pyarrow, numpy.
"""
import argparse
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
import pysam

import contig_names as CN

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
CODING = {"synonymous_variant": "synonymous",
          "missense_variant": "missense",
          "stop_gained": "stop_gained",
          "stop_lost": "stop_lost",
          "start_lost": "start_lost"}
# preference order when one variant hits several coding effects on one transcript
PRIORITY = ["stop_gained", "stop_lost", "start_lost", "missense_variant", "synonymous_variant"]

_attr_re = re.compile(r'(\w+) "([^"]*)"')

_ver_re = re.compile(r'\.\d+$')


def strip_version(tid):
    """'FBtr0306592.1' -> 'FBtr0306592'.

    snpEff appends transcript_version to its Feature_ID; the GTF (and so the
    Stage 0C table) uses the bare ID.
    """
    return _ver_re.sub("", tid) if tid else tid


def to_rna(b):
    return "U" if b == "T" else b


def spectrum_channel(ref, alt, strand):
    """Strand-aware RNA channel, or None if strand is unusable."""
    if strand == "+":
        f, t = ref, alt
    elif strand == "-":
        f, t = COMP.get(ref), COMP.get(alt)
    else:
        return None
    if f is None or t is None:
        return None
    return f"{to_rna(f)}>{to_rna(t)}"


def pick_gene_feature(gtf, preference=("gene", "transcript", "exon")):
    """
    First feature type in `preference` present in the GTF. genes.gtf has no
    'gene' rows; transcript spans give the same transcribed strand, since a
    gene's span is the union of its transcripts' spans.
    """
    seen = set()
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 3)
            if len(f) >= 3:
                seen.add(f[2])
    for feat in preference:
        if feat in seen:
            return feat, seen
    raise ValueError(f"{gtf}: none of {preference} present; feature types are "
                     f"{sorted(seen)}")


def build_strand_code(gtf, sizes, feature=None):
    """{contig: int8 array} 0=none 1=+ 2=- 3=ambiguous, from gene-level features."""
    if feature is None:
        feature, seen = pick_gene_feature(gtf)
        if feature != "gene":
            print(f"[STAGE2B] GTF has no 'gene' feature; deriving transcribed "
                  f"strand from '{feature}' instead (present types: "
                  f"{sorted(seen)})", file=sys.stderr)
    code = {c: np.zeros(n, dtype=np.int8) for c, n in sizes.items()}
    n_used = 0
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            c, s, e, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            if c not in code:
                continue
            val = 1 if strand == "+" else 2
            seg = code[c][s:e]
            opp = ((seg == 1) & (val == 2)) | ((seg == 2) & (val == 1)) | (seg == 3)
            seg[seg == 0] = val
            seg[opp] = 3
            code[c][s:e] = seg
            n_used += 1

    # an empty strand array would leave every spectrum channel empty
    stranded = sum(int((v != 0).sum()) for v in code.values())
    print(f"[STAGE2B] strand from '{feature}': {n_used} features, "
          f"{stranded} bases with a strand", file=sys.stderr)
    if stranded == 0:
        raise SystemExit(
            f"[STAGE2B] FATAL: no strand information derived from {gtf} using "
            f"feature '{feature}'. Every mutation would get an empty spectrum "
            f"channel. Check that the GTF's seqnames match the canonical .fai.")
    return code


def strand_of(code, chrom, pos1):
    if chrom not in code:
        return None
    v = code[chrom][pos1 - 1]
    return {1: "+", 2: "-"}.get(int(v), None)  # 0/3 -> None


def parse_ann(vcf_path, rep_transcripts):
    """
    (chrom, pos, ref, alt) -> dict(consequence, gene_id, transcript_id, effect,
    impact).

    The consequence is taken only from representative transcripts. SNVs with
    no coding effect on one are labelled 'other'; their gene/effect fields then
    come from snpEff's first annotation for the allele, for reference only.
    """
    lookup = {}
    counts = Counter()
    vf = pysam.VariantFile(vcf_path)
    for rec in vf:
        ann = rec.info.get("ANN")
        if ann is None:
            continue
        if isinstance(ann, str):
            ann = (ann,)
        for alt in (rec.alts or ()):
            chosen = None      # (priority_rank, fields): best representative-tx coding hit
            first = None       # snpEff's first annotation for this allele
            for entry in ann:
                parts = entry.split("|")
                if len(parts) < 7 or parts[0] != alt:
                    continue
                effect, impact, gene_id, feat_id = parts[1], parts[2], parts[4], parts[6]
                if first is None:
                    first = (effect, impact, gene_id, feat_id)
                if strip_version(feat_id) in rep_transcripts:
                    effs = effect.split("&")
                    for pr, term in enumerate(PRIORITY):
                        if term in effs:
                            if chosen is None or pr < chosen[0]:
                                chosen = (pr, (effect, impact, gene_id, feat_id))
                            break
            if chosen is not None:
                effect, impact, gene_id, feat_id = chosen[1]
                consequence = CODING[PRIORITY[chosen[0]]]
                counts["rep_coding"] += 1
            elif first is not None:
                effect, impact, gene_id, feat_id = first
                consequence = "other"
                if any(t in effect.split("&") for t in PRIORITY):
                    counts["other_coding_on_non_representative"] += 1
            else:
                continue
            lookup[(rec.chrom, rec.pos, rec.ref, alt)] = {
                "consequence": consequence, "gene_id": gene_id,
                "transcript_id": feat_id, "effect": effect, "impact": impact,
            }
            counts[consequence] += 1
    vf.close()
    return lookup, counts


def load_fai(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--annotated-vcf", required=True)
    ap.add_argument("--stage1-dir", required=True)
    ap.add_argument("--gtf", required=True, help="canonical GTF")
    ap.add_argument("--ns-table", required=True, help="Stage 0C ns_sites.parquet (representative transcripts)")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-coding-frac", type=float, default=0.10,
                    help="Stop if fewer than this fraction of nuclear SNVs get a "
                         "coding consequence (snpEff treating transcripts as non-coding)")
    ap.add_argument("--samples-list", required=True,
                    help="Samples to process (one ID per line)")
    args = ap.parse_args()

    sizes = load_fai(args.fai)
    os.makedirs(args.out_dir, exist_ok=True)

    rep_tx = {strip_version(t) for t in
              pd.read_parquet(args.ns_table, columns=["transcript_id"])["transcript_id"].unique()}
    print(f"[STAGE2B] representative transcripts: {len(rep_tx)}", file=sys.stderr)

    lookup, ann_counts = parse_ann(args.annotated_vcf, rep_tx)
    print(f"[STAGE2B] annotated SNVs: {len(lookup)} | {dict(ann_counts)}", file=sys.stderr)

    # If transcript IDs did not match, nothing would be classified as coding.
    n_rep_hit = ann_counts["rep_coding"]
    frac_rep = n_rep_hit / len(lookup) if lookup else 0.0
    print(f"[STAGE2B] coding on a representative transcript: "
          f"{n_rep_hit} ({100*frac_rep:.1f}%); coding only on another isoform "
          f"(labelled 'other'): {ann_counts['other_coding_on_non_representative']}",
          file=sys.stderr)
    if lookup and frac_rep < 0.01:
        raise SystemExit(
            "[STAGE2B] FATAL: essentially no snpEff annotation matched a "
            "representative transcript from the N/S table. The two files "
            "disagree about transcript-id format (version suffix?), so "
            "synonymous/missense would not correspond to N and S.")

    strand_code = build_strand_code(args.gtf, sizes)

    samples = [l.strip() for l in open(args.samples_list) if l.strip()]
    files = [os.path.join(args.stage1_dir, f"{s}.snv.parquet") for s in samples]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        sys.exit(f"[STAGE2B] FATAL: {len(missing)} Stage 1 table(s) missing, "
                 f"e.g. {missing[:5]}")
    gcount = Counter()
    ccount = Counter()
    for f in files:
        df = pd.read_parquet(f)
        cons, genes, txs, effs, imps, strands, chans = [], [], [], [], [], [], []
        for chrom, pos, ref, alt in zip(df["chrom"], df["pos"], df["ref"], df["alt"]):
            a = lookup.get((chrom, pos, ref, alt))
            if a:
                cons.append(a["consequence"]); genes.append(a["gene_id"])
                txs.append(a["transcript_id"]); effs.append(a["effect"]); imps.append(a["impact"])
            else:
                cons.append("other"); genes.append(None); txs.append(None)
                effs.append(None); imps.append(None)
            st = strand_of(strand_code, chrom, pos)
            strands.append(st)
            ch = spectrum_channel(ref, alt, st)
            chans.append(ch)
            gcount["total"] += 1
            gcount["with_channel" if ch else "no_channel"] += 1
        df["consequence"] = cons
        df["gene_id"] = genes
        df["transcript_id"] = txs
        df["effect"] = effs
        df["impact"] = imps
        df["gene_strand"] = strands
        df["spectrum_channel"] = chans
        sample = os.path.basename(f).replace(".snv.parquet", "")
        df.to_parquet(os.path.join(args.out_dir, f"{sample}.annotated.parquet"), index=False)

        # coding-call counts by compartment, for the guard below
        comp = df["chrom"].map(CN.compartment)
        cod = df["consequence"].isin(tuple(CODING.values()))
        for c in ("nuclear", "mito"):
            m = (comp == c)
            ccount[f"{c}_total"] += int(m.sum())
            ccount[f"{c}_coding"] += int((m & cod).sum())

    tot = gcount["total"] or 1
    print(f"[STAGE2B] global SNVs={gcount['total']} "
          f"with_channel={gcount['with_channel']} "
          f"no_channel={gcount['no_channel']} "
          f"({100*gcount['no_channel']/tot:.1f}% no strand)", file=sys.stderr)

    # ---- coding-call guard ----
    # If snpEff treats the nuclear transcripts as non-coding, nearly every SNV
    # becomes 'other' without any error; stop instead.
    for c in ("nuclear", "mito"):
        n, k = ccount[f"{c}_total"], ccount[f"{c}_coding"]
        frac = k / n if n else 0.0
        print(f"[STAGE2B] {c:8s} SNVs={n:>10,} coding={k:>9,} ({100*frac:5.2f}%)",
              file=sys.stderr)
    n_nuc = ccount["nuclear_total"]
    f_nuc = ccount["nuclear_coding"] / n_nuc if n_nuc else 0.0
    if n_nuc > 10000 and f_nuc < args.min_coding_frac:
        raise SystemExit(
            f"[STAGE2B] FATAL: only {100*f_nuc:.2f}% of nuclear SNVs got a coding "
            f"consequence (floor {100*args.min_coding_frac:.0f}%). snpEff is almost "
            f"certainly annotating the nuclear transcripts as non-coding - check "
            f"the Stage 0B build log and run\n"
            f"    python3 snpeff_gtf.py inspect --gtf <canonical gtf>\n"
            f"Override with --min-coding-frac 0 if this is genuinely expected.")


if __name__ == "__main__":
    main()
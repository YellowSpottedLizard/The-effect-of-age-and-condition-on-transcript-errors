#!/usr/bin/env python3
"""
Stage 2B (part 1) - union VCF of distinct SNVs
==============================================
A consequence depends only on (chrom, pos, ref, alt), not on the sample, so
the distinct SNVs of all samples' Stage 1 tables are collected and annotated
with snpEff once. Writes a minimal VCF sorted in .fai contig order.

Requires: pandas, pyarrow.
"""
import argparse
import os
import sys

import pandas as pd


def load_fai_order(fai):
    order = []
    with open(fai) as fh:
        for line in fh:
            order.append(line.split("\t")[0])
    return order


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--stage1-dir", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--samples-list", required=True,
                    help="Samples to include (one ID per line). Only their Stage 1 "
                         "tables are read, so other files in --stage1-dir are ignored.")
    args = ap.parse_args()

    samples = [l.strip() for l in open(args.samples_list) if l.strip()]
    files = [os.path.join(args.stage1_dir, f"{s}.snv.parquet") for s in samples]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        sys.exit(f"[STAGE2B] FATAL: {len(missing)} sample(s) in the samples "
                 f"list have no Stage 1 table, e.g. {missing[:5]}")

    frames = []
    for f in files:
        df = pd.read_parquet(f, columns=["chrom", "pos", "ref", "alt"])
        frames.append(df)
    allv = pd.concat(frames, ignore_index=True).drop_duplicates()

    order = load_fai_order(args.fai)
    rank = {c: i for i, c in enumerate(order)}
    allv["__r"] = allv["chrom"].map(rank)
    n_off = int(allv["__r"].isna().sum())
    if n_off:
        print(f"[STAGE2B] WARNING: {n_off} SNVs on contigs absent from the .fai dropped",
              file=sys.stderr)
    allv = allv.dropna(subset=["__r"]).sort_values(["__r", "pos"]).reset_index(drop=True)

    with open(args.out_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        with open(args.fai) as fai:
            for line in fai:
                c, n = line.split("\t")[:2]
                fh.write(f"##contig=<ID={c},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in zip(allv["chrom"], allv["pos"], allv["ref"], allv["alt"]):
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"[STAGE2B] {len(allv)} distinct SNVs -> {args.out_vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()
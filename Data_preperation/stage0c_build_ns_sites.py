#!/usr/bin/env python3
"""
Stage 0C - Nei-Gojobori N/S site table
======================================
For every coding position, computes the number of nonsynonymous (n_site) and
synonymous (s_site) sites. Stage 3 sums these over the CDS positions that are
callable in each sample to obtain the dN/dS denominators N and S.

Method (Nei & Gojobori 1986):
  At each codon position, consider the three possible single-nucleotide
  changes and translate the mutated codon.
    s_site = (number of synonymous changes) / 3
    n_site = (number of missense changes)   / 3
  Changes that create a stop codon are excluded from both, so n_site and
  s_site describe the same classes of events that Stage 3 counts as missense
  and synonymous. At positions where a change can create a stop codon,
  n_site + s_site < 1. Reference stop codons are skipped.

  * One representative transcript per gene (the longest CDS; ties broken by
    transcript ID), so each position is counted once.
  * The mitochondrial contig (default 'M') uses NCBI translation table 5
    (invertebrate mitochondrial); all other contigs use table 1.
  * Codons containing N are skipped.
  * If two genes' representative CDS cover the same position (e.g. antisense
    overlap), the first one encountered is kept and the collision is counted.

Outputs (named from --out, e.g. ns_sites.parquet):
  ns_sites.parquet           chrom, pos (1-based), n_site, s_site, gene_id,
                             transcript_id. Used by Stages 2B and 3.
  ns_sites.channels.parquet  one row per classified change: chrom, pos,
                             channel (transcript-strand ref>alt in U notation),
                             is_syn, weight (1/3). Used by Stage 3 for
                             spectrum-weighted N and S.
  ns_sites.cds.bed           the positions of ns_sites.parquet as merged
                             intervals; defines the CDS region in Stages 2A
                             and 3.

Requires: pysam, pandas, pyarrow, biopython.
"""
import argparse
import re
import sys
from collections import defaultdict, Counter
from functools import lru_cache

import numpy as np
import pysam
import pandas as pd
from Bio.Data import CodonTable

BASES = ("A", "C", "G", "T")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

_TBL = {
    1: CodonTable.unambiguous_dna_by_id[1],
    5: CodonTable.unambiguous_dna_by_id[5],
}

_attr_re = re.compile(r'(\w+) "([^"]*)"')


def revcomp(s):
    return "".join(COMP.get(b, "N") for b in reversed(s))


@lru_cache(maxsize=None)
def translate_codon(codon, table_id):
    """Return amino acid char, '*' for stop, or None if untranslatable."""
    codon = codon.upper()
    if len(codon) != 3 or "N" in codon:
        return None
    t = _TBL[table_id]
    if codon in t.stop_codons:
        return "*"
    return t.forward_table.get(codon)  # None if somehow ambiguous


def _rna(b):
    return "U" if b == "T" else b


@lru_cache(maxsize=None)
def site_fractions(codon, pos_in_codon, table_id):
    """
    For one position of a sense codon, return (n_site, s_site, changes), or
    None for a reference stop codon.

    `changes` holds (channel, is_syn, weight) for each synonymous or missense
    change, with channel the coding-strand ref>alt in U notation and weight
    1/3, so that the weights of the synonymous changes sum to s_site and those
    of the missense changes to n_site. Changes to a stop codon are omitted.
    """
    aa_ref = translate_codon(codon, table_id)
    if aa_ref is None or aa_ref == "*":
        return None
    ref_base = codon[pos_in_codon]
    changes = []
    n_syn = n_mis = 0
    for b in BASES:
        if b == ref_base:
            continue
        aa_mut = translate_codon(codon[:pos_in_codon] + b + codon[pos_in_codon + 1:], table_id)
        if aa_mut == "*":
            continue
        is_syn = aa_mut == aa_ref
        n_syn += is_syn
        n_mis += not is_syn
        changes.append((f"{_rna(ref_base)}>{_rna(b)}", is_syn, 1.0 / 3))
    return (n_mis / 3, n_syn / 3, tuple(changes))


def parse_attrs(field):
    return dict(_attr_re.findall(field))


def main():
    ap = argparse.ArgumentParser(
        description="Build the Nei-Gojobori per-CDS-position N/S site table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--gtf", required=True, help="Canonical GTF (mito named M)")
    ap.add_argument("--fasta", required=True, help="Canonical FASTA (.fa + .fai)")
    ap.add_argument("--out", required=True, help="Output parquet path")
    ap.add_argument("--mito-contig", default="M",
                    help="Contig that uses translation table 5")
    args = ap.parse_args()

    fasta = pysam.FastaFile(args.fasta)
    fa_contigs = set(fasta.references)

    # ---- collect CDS segments per transcript ----
    # tx -> dict(gene, chrom, strand, segs=[(start1,end1,frame)])
    tx = {}
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            chrom, start, end, strand, frame, attrs = f[0], int(f[3]), int(f[4]), f[6], f[7], f[8]
            a = parse_attrs(attrs)
            tid = a.get("transcript_id")
            gid = a.get("gene_id")
            if tid is None:
                continue
            fr = int(frame) if frame in ("0", "1", "2") else 0
            rec = tx.setdefault(tid, {"gene": gid, "chrom": chrom, "strand": strand, "segs": []})
            rec["segs"].append((start, end, fr))

    # ---- one representative transcript per gene: the longest CDS ----
    by_gene = defaultdict(list)
    for tid, rec in tx.items():
        rec["cds_len"] = sum(e - s + 1 for s, e, _ in rec["segs"])
        by_gene[rec["gene"]].append(tid)
    chosen = [sorted(tids, key=lambda t: (-tx[t]["cds_len"], t))[0]
              for tids in by_gene.values()]

    # ---- walk codons, emit per-position n/s ----
    counts = Counter()
    seen = {}  # (chrom,pos1) -> (n,s,gene,tx)  (first wins)
    channel_rows = []  # (chrom,pos1,channel,is_syn,weight)

    for tid in chosen:
        rec = tx[tid]
        chrom, strand = rec["chrom"], rec["strand"]
        if chrom not in fa_contigs:
            counts["tx_contig_absent"] += 1
            continue
        table_id = 5 if chrom == args.mito_contig else 1

        # order segments in transcript (5'->3') direction
        segs = sorted(rec["segs"], key=lambda x: x[0], reverse=(strand == "-"))

        # build coding sequence + parallel genomic coords (1-based)
        coding = []
        gcoord = []
        for (s1, e1, _fr) in segs:
            seq = fasta.fetch(chrom, s1 - 1, e1).upper()  # genomic 5'->3'
            if strand == "+":
                coding.append(seq)
                gcoord.extend(range(s1, e1 + 1))
            else:
                coding.append(revcomp(seq))
                gcoord.extend(range(e1, s1 - 1, -1))
        coding = "".join(coding)

        # trim leading phase of the FIRST segment in transcript order
        lead = segs[0][2]
        if lead:
            coding = coding[lead:]
            gcoord = gcoord[lead:]

        counts["tx_used"] += 1

        # iterate complete codons
        ncodons = len(coding) // 3
        for ci in range(ncodons):
            codon = coding[ci * 3: ci * 3 + 3]
            gc = gcoord[ci * 3: ci * 3 + 3]
            if "N" in codon:
                counts["codons_with_N"] += 1
                continue
            for j in range(3):
                frac = site_fractions(codon, j, table_id)
                if frac is None:
                    counts["ref_stop_sites"] += 1
                    continue
                n_site, s_site, changes = frac
                key = (chrom, gc[j])
                if key in seen:
                    counts["position_collisions"] += 1
                    continue
                seen[key] = (n_site, s_site, rec["gene"], tid)
                for channel, is_syn, w in changes:
                    channel_rows.append((chrom, gc[j], channel, is_syn, w))
                counts["sites_emitted"] += 1

    fasta.close()

    # ---- write ----
    rows = [
        {"chrom": c, "pos": p, "n_site": v[0], "s_site": v[1],
         "gene_id": v[2], "transcript_id": v[3]}
        for (c, p), v in seen.items()
    ]
    df = pd.DataFrame(rows).sort_values(["chrom", "pos"]).reset_index(drop=True)
    df.to_parquet(args.out, index=False)

    # channel-resolved table for spectrum-weighted dN/dS
    ch = pd.DataFrame(channel_rows, columns=["chrom", "pos", "channel", "is_syn", "weight"])
    ch_out = args.out.replace(".parquet", "") + ".channels.parquet"
    ch.to_parquet(ch_out, index=False)

    # ---- CDS BED from the same positions, so the CDS region used in Stages
    # 2A and 3 covers exactly the positions that have N/S sites ----
    bed_out = args.out.replace(".parquet", "") + ".cds.bed"
    n_iv = 0
    with open(bed_out, "w") as fh:
        for chrom, sub in df.groupby("chrom", sort=True):
            p0 = np.sort(sub["pos"].to_numpy()) - 1          # 0-based
            brk = np.where(np.diff(p0) > 1)[0]
            starts = np.concatenate(([p0[0]], p0[brk + 1]))
            ends = np.concatenate((p0[brk], [p0[-1]])) + 1
            for s, e in zip(starts, ends):
                fh.write(f"{chrom}\t{int(s)}\t{int(e)}\n")
                n_iv += 1
    print(f"[STAGE0C] wrote {bed_out} ({n_iv} intervals, "
          f"{len(df)} CDS positions)", file=sys.stderr)

    tot_n = df["n_site"].sum()
    tot_s = df["s_site"].sum()

    # nuclear and mito reported separately (different genetic codes)
    is_mito = df["chrom"] == args.mito_contig
    for label, sub in (("nuclear", df[~is_mito]), ("mito", df[is_mito])):
        if len(sub) == 0:
            print(f"[STAGE0C] {label:8s}: no CDS positions", file=sys.stderr)
            continue
        print(f"[STAGE0C] {label:8s}: {len(sub)} CDS positions, "
              f"N={sub['n_site'].sum():.1f} S={sub['s_site'].sum():.1f} "
              f"mean s_site={sub['s_site'].mean():.4f}", file=sys.stderr)

    print(f"[STAGE0C] transcripts used      : {counts['tx_used']}", file=sys.stderr)
    print(f"[STAGE0C] CDS positions emitted : {counts['sites_emitted']}", file=sys.stderr)
    print(f"[STAGE0C] position collisions   : {counts['position_collisions']}", file=sys.stderr)
    print(f"[STAGE0C] codons with N         : {counts['codons_with_N']}", file=sys.stderr)
    print(f"[STAGE0C] reference stop sites  : {counts['ref_stop_sites']}", file=sys.stderr)
    print(f"[STAGE0C] total N sites         : {tot_n:.1f}", file=sys.stderr)
    print(f"[STAGE0C] total S sites         : {tot_s:.1f}", file=sys.stderr)
    print(f"[STAGE0C] genome-wide N/S ratio : {tot_n/tot_s:.3f}", file=sys.stderr)
    print(f"[STAGE0C] wrote {args.out} ({len(df)} rows)", file=sys.stderr)
    print(f"[STAGE0C] wrote {ch_out} ({len(ch)} channel rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
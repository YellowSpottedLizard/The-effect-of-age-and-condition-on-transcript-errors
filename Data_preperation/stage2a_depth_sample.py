#!/usr/bin/env python3
"""
Stage 2A (per sample) - depth aggregates and callable regions
=============================================================
Runs samtools depth on one sample's BAM, counting bases with mapping quality
>= 1 and base quality >= 73, as the denominator for the error rates. (Mutation
alleles are kept in Stage 1 if their mean alt base quality is >= 74.)

Outputs, for one sample:
  <prefix>.aggregates.json       depth sums and positions with depth >=
                                 MIN_DEPTH, split total / CDS / non-CDS
  <prefix>.depth.npz             per-contig (pos0, depth) over the analysis
                                 region; used by Stage 2A reduce and Stage 3
  <prefix>.callable.bed          runs with depth >= MIN_DEPTH; used by Stage 3

Nuclear and mitochondrial figures are reported separately. The mitochondrion
is transcribed by its own polymerase and has 10-100x the nuclear depth, so
pooling would let a few thousand mito bases dominate a genome-wide rate.
Nuclear columns are unsuffixed, mito columns end in '_mito' and pooled ones in
'_all'.

BAM contig names are mapped to canonical names (contig_names.py); canonical
BEDs are translated to the BAM's names before being passed to samtools.

Requires: numpy, samtools on PATH.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

import contig_names as CN


def load_fai(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def load_bed_mask(bed_path, sizes, mito=CN.CANONICAL_MITO):
    """Return {contig: bool array} marking positions inside the BED (0-based)."""
    mask = {c: np.zeros(n, dtype=bool) for c, n in sizes.items()}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            c, s, e = CN.canonicalise(f[0], mito), int(f[1]), int(f[2])
            if c in mask:
                mask[c][s:e] = True
    return mask


# ---------------------------------------------------------------------------
# contig-name plumbing between the canonical reference and the BAM
# ---------------------------------------------------------------------------
def bam_contigs(bam):
    """Contig names as spelled in the BAM header."""
    out = subprocess.run(["samtools", "view", "-H", bam],
                         stdout=subprocess.PIPE, text=True, check=True).stdout
    names = []
    for line in out.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    names.append(field[3:])
    return names


def translate_bed(src_bed, canon_to_bam, dst_bed):
    """
    Rewrite a canonical-named BED into the BAM's contig naming.

    Returns (n_intervals_written, set_of_canonical_contigs_with_no_BAM_match),
    so the caller can complain about dropped contigs rather than lose them
    quietly.
    """
    written = 0
    missing = set()
    with open(src_bed) as fh, open(dst_bed, "w") as out:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            c = CN.canonicalise(f[0])
            bam_name = canon_to_bam.get(c)
            if bam_name is None:
                missing.add(c)
                continue
            f[0] = bam_name
            out.write("\t".join(f) + "\n")
            written += 1
    return written, missing


def aggregate_depth(lines, cds_mask, sizes, min_depth, mito_contig=CN.CANONICAL_MITO):
    """
    Core aggregation over `samtools depth` lines (chrom, pos1, depth).

    Returns (aggregates dict, {contig: (pos0_array, depth_array)}).
    cds_mask: {contig: bool array} of CDS positions; sizes: {contig: length}.

    Aggregates are accumulated per compartment ('nuclear' / 'mito') and the
    pooled '_all' figures are derived at the end, so the three views are
    consistent with one another by construction.
    """
    fields = ("total_depth_sum", "cds_depth_sum", "noncds_depth_sum",
              "positions_covered_total", "positions_covered_cds",
              "positions_covered_noncds")
    acc = {comp: {f: 0 for f in fields} for comp in ("nuclear", "mito")}

    pos_lists = {c: [] for c in sizes}
    dep_lists = {c: [] for c in sizes}

    for line in lines:
        if not line:
            continue
        chrom_raw, pos_s, dep_s = line.rstrip("\n").split("\t")
        depth = int(dep_s)
        if depth <= 0:
            continue
        chrom = CN.canonicalise(chrom_raw, mito_contig)
        if chrom not in sizes:
            continue
        pos0 = int(pos_s) - 1
        if pos0 < 0 or pos0 >= sizes[chrom]:
            continue
        is_cds = bool(cds_mask[chrom][pos0]) if chrom in cds_mask else False
        a = acc["mito" if chrom == mito_contig else "nuclear"]

        a["total_depth_sum"] += depth
        if is_cds:
            a["cds_depth_sum"] += depth
        else:
            a["noncds_depth_sum"] += depth

        if depth >= min_depth:
            a["positions_covered_total"] += 1
            if is_cds:
                a["positions_covered_cds"] += 1
            else:
                a["positions_covered_noncds"] += 1

        pos_lists[chrom].append(pos0)
        dep_lists[chrom].append(depth)

    # nuclear keeps the original column names; mito gets '_mito'; pooled '_all'
    agg = {}
    for f in fields:
        agg[f] = acc["nuclear"][f]
        agg[f + "_mito"] = acc["mito"][f]
        agg[f + "_all"] = acc["nuclear"][f] + acc["mito"][f]

    arrays = {}
    for c in sizes:
        if pos_lists[c]:
            arrays[c] = (np.asarray(pos_lists[c], dtype=np.uint32),
                         np.asarray(dep_lists[c], dtype=np.uint32))
    return agg, arrays


def runs_from_positions(pos0, depth, min_depth):
    """Yield (start0, end) BED runs where depth >= min_depth (pos0 sorted)."""
    keep = pos0[depth >= min_depth]
    if keep.size == 0:
        return
    keep = np.sort(keep)
    brk = np.where(np.diff(keep) > 1)[0]
    starts = np.concatenate(([keep[0]], keep[brk + 1]))
    ends = np.concatenate((keep[brk], [keep[-1]])) + 1
    for s, e in zip(starts, ends):
        yield int(s), int(e)


def samtools_depth(bam, bed, mapq, baseq):
    # -d 0: no per-position depth cap. Some samtools versions cap at 8000,
    # which would truncate the mitochondrion (~22,000x).
    cmd = ["samtools", "depth", "-d", "0", "-Q", str(mapq), "-q", str(baseq),
           "-b", bed, bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in proc.stdout:
        yield line
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"samtools depth failed: {' '.join(cmd)}")


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fai", required=True, help="canonical genome .fa.fai")
    ap.add_argument("--analysis-bed", required=True, help="whole callable genome (one interval per contig)")
    ap.add_argument("--cds-bed", required=True,
                    help="CDS positions from Stage 0C (ns_sites.cds.bed)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--mapq", type=int, default=1)
    ap.add_argument("--baseq", type=int, default=73)
    ap.add_argument("--mito-contig", default=CN.CANONICAL_MITO,
                    help="Canonical mitochondrial contig; its statistics are "
                         "reported separately from the nuclear genome")
    args = ap.parse_args()

    sizes = load_fai(args.fai)
    cds_mask = load_bed_mask(args.cds_bed, sizes, args.mito_contig)

    # ---- reconcile canonical contig names with the BAM's own naming ----
    bam_names = bam_contigs(args.bam)
    canon_to_bam = CN.invert_rename_map(bam_names, set(sizes), args.mito_contig)
    absent = sorted(set(sizes) - set(canon_to_bam))
    if absent:
        print(f"[STAGE2A] WARNING: canonical contigs with no match in the BAM header: "
              f"{absent}", file=sys.stderr)
    if args.mito_contig in sizes and args.mito_contig not in canon_to_bam:
        print("[STAGE2A] WARNING: the mitochondrion is in the reference but not in "
              "the BAM header; all mito statistics will be zero for this sample.",
              file=sys.stderr)

    tmpdir = tempfile.mkdtemp(prefix=f"stage2a_{args.sample_id}_")
    an_bed = os.path.join(tmpdir, "analysis.bamnames.bed")
    n_an, miss_an = translate_bed(args.analysis_bed, canon_to_bam, an_bed)
    if miss_an:
        print(f"[STAGE2A] WARNING: intervals dropped for contigs absent from the BAM: "
              f"{sorted(miss_an)}", file=sys.stderr)
    if n_an == 0:
        raise SystemExit("[STAGE2A] FATAL: no analysis intervals survived contig-name "
                         "translation - the BAM and the canonical reference share no "
                         "contigs. Check the BAM header against the canonical .fai.")

    # main pass over the analysis region
    agg, arrays = aggregate_depth(
        samtools_depth(args.bam, an_bed, args.mapq, args.baseq),
        cds_mask, sizes, args.min_depth, args.mito_contig,
    )
    agg["sample_id"] = args.sample_id
    agg["min_depth"] = args.min_depth
    agg["mito_contig"] = args.mito_contig

    # sparse per-contig depth (canonical contig names)
    npz = {}
    for c, (p, d) in arrays.items():
        npz[f"{c}__pos"] = p
        npz[f"{c}__dep"] = d
    np.savez_compressed(f"{args.out_prefix}.depth.npz", **npz)

    # callable BED (canonical contig names)
    with open(f"{args.out_prefix}.callable.bed", "w") as fh:
        for c in sizes:
            if c in arrays:
                p, d = arrays[c]
                for s, e in runs_from_positions(p, d, args.min_depth):
                    fh.write(f"{c}\t{s}\t{e}\n")

    with open(f"{args.out_prefix}.aggregates.json", "w") as fh:
        json.dump(agg, fh, indent=2)

    shutil.rmtree(tmpdir, ignore_errors=True)

    print(f"[STAGE2A] {args.sample_id}: covered>={args.min_depth} "
          f"nuclear total={agg['positions_covered_total']} "
          f"cds={agg['positions_covered_cds']} "
          f"noncds={agg['positions_covered_noncds']} | "
          f"mito total={agg['positions_covered_total_mito']} "
          f"cds={agg['positions_covered_cds_mito']}", file=sys.stderr)
    print(f"[STAGE2A] {args.sample_id}: depth sum "
          f"nuclear={agg['total_depth_sum']} mito={agg['total_depth_sum_mito']}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Stage 2A reduce - cohort per-base depth summary
===============================================
Combines the per-sample <sample>.depth.npz files from stage2a_depth_sample.py
into two per-base arrays for each contig:

  cohort_sum_depth.<contig>.npy   depth summed over all samples
  cohort_n_callable.<contig>.npy  number of samples with depth >= MIN_DEPTH

Stage 3 uses them for the across-sample occurrence filter (summed AO / summed
depth at a position) and to define the analysed region (n_callable).

Stops if any sample in the list has no depth.npz, since a missing sample would
silently lower the summed depth at every position.

Requires: numpy.
"""
import argparse
import os
import sys

import numpy as np


def load_fai(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def read_samples(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--depth-dir", required=True, help="dir with <sample>.depth.npz etc.")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--samples-list", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-depth", type=int, default=100)
    args = ap.parse_args()

    sizes = load_fai(args.fai)
    samples = read_samples(args.samples_list)
    os.makedirs(args.out_dir, exist_ok=True)

    missing = [s for s in samples
               if not os.path.exists(os.path.join(args.depth_dir, f"{s}.depth.npz"))]
    if missing:
        sys.exit(f"[REDUCE] FATAL: {len(missing)} sample(s) have no depth.npz, "
                 f"e.g. {missing[:5]}")

    # ---- cohort n_callable + sum_depth per contig ----
    for contig, n in sizes.items():
        sum_depth = np.zeros(n, dtype=np.uint64)
        n_callable = np.zeros(n, dtype=np.uint16)
        present = 0
        for s in samples:
            npz_path = os.path.join(args.depth_dir, f"{s}.depth.npz")
            with np.load(npz_path) as z:
                pk, dk = f"{contig}__pos", f"{contig}__dep"
                if pk not in z.files:
                    continue
                pos = z[pk].astype(np.int64)
                dep = z[dk].astype(np.uint64)
            sum_depth[pos] += dep
            n_callable[pos[dep >= args.min_depth]] += 1
            present += 1
        np.save(os.path.join(args.out_dir, f"cohort_sum_depth.{contig}.npy"), sum_depth)
        np.save(os.path.join(args.out_dir, f"cohort_n_callable.{contig}.npy"), n_callable)
        print(f"[REDUCE] {contig}: {present} samples, "
              f"max n_callable={int(n_callable.max()) if n else 0}", file=sys.stderr)


if __name__ == "__main__":
    main()
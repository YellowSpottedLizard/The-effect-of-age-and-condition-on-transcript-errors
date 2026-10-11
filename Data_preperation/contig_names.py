#!/usr/bin/env python3
"""
contig_names.py - one contig vocabulary for the whole pipeline.

The reference FASTA, GTF, freebayes VCFs and STAR BAMs may name contigs
differently. Every stage maps names through canonicalise() onto:

    2L 2R 3L 3R 4 X Y M

  'chr2L' -> '2L'                              (UCSC prefix stripped)
  'M', 'MT', 'chrM', 'chrMT', 'mito', ...       -> 'M'
  anything containing 'mitochond'               -> 'M'
                                                   (e.g. Ensembl
                                                   'mitochondrion_genome')

Other names (e.g. unplaced scaffolds) are returned unchanged and are dropped
by each stage's contig whitelist.
"""

CANONICAL_NUCLEAR = ("2L", "2R", "3L", "3R", "4", "X", "Y")
CANONICAL_MITO = "M"
CANONICAL_ALL = CANONICAL_NUCLEAR + (CANONICAL_MITO,)

_MITO_EXACT = {"m", "mt", "chrm", "chrmt", "mito", "chrmito"}


def canonicalise(name, mito=CANONICAL_MITO):
    """Map any spelling of a contig name onto the pipeline's canonical name."""
    if name is None:
        return None
    n = str(name).strip()
    low = n.lower()
    if low in _MITO_EXACT or "mitochond" in low:
        return mito
    if low.startswith("chr") and len(n) > 3:
        n = n[3:]
        if n.lower() in _MITO_EXACT:
            return mito
    return n


def invert_rename_map(names, wanted, mito=CANONICAL_MITO):
    """
    {canonical_name: original_name} restricted to `wanted` canonical names.

    Used when we have to hand a canonical BED to a tool that speaks the BAM's
    own dialect (samtools depth), so the interval file is translated into the
    BAM's names on the way in and the output is translated back on the way out.
    Raises on an ambiguous mapping (two source contigs collapsing onto the same
    canonical name) rather than picking one arbitrarily.
    """
    out = {}
    for n in names:
        c = canonicalise(n, mito)
        if c not in wanted:
            continue
        if c in out and out[c] != n:
            raise ValueError(
                f"ambiguous contig naming: both {out[c]!r} and {n!r} canonicalise to {c!r}")
        out[c] = n
    return out


def is_mito(name, mito=CANONICAL_MITO):
    return canonicalise(name, mito) == mito


def compartment(name, mito=CANONICAL_MITO):
    """'mito' or 'nuclear'."""
    return "mito" if is_mito(name, mito) else "nuclear"

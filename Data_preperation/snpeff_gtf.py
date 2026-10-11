#!/usr/bin/env python3
"""
snpeff_gtf.py - prepare a GTF for snpEff and check the resulting database.

snpEff's GTF22 parser reads the file in a single pass and resolves each row's
parent transcript by ID. If a child row (exon, CDS) appears before its
transcript row, snpEff creates a stub transcript and later rejects the real
transcript row as a duplicate. The stub never receives transcript_biotype
"protein_coding", so the transcript is treated as non-coding and coding SNVs
are annotated as non_coding_transcript_exon_variant. Exon IDs shared between
sibling transcripts cause the same problem.

Subcommands
-----------
  fix           rewrite the GTF with parent rows first and unique exon_ids
                (Stage 0A)
  inspect       report row ordering and exon_id collisions; with
                --assert-clean, exit 1 if either is present (Stage 0A)
  make-testvcf  write a VCF of SNVs at real CDS positions, nuclear and mito
                (Stage 0B)
  score         report the fraction of those SNVs that snpEff annotated with
                a coding consequence; exit 1 below --min-coding-frac (Stage 0B)

Requires: pysam for make-testvcf; the other subcommands use the standard
library only.
"""
import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

MITO = "M"
CODING_TERMS = {"synonymous_variant", "missense_variant",
                "stop_gained", "stop_lost", "start_lost"}

# snpEff's GTF22 parser resolves parents by ID as it streams the file, so a
# child feature appearing before its parent makes it invent the parent. Rank
# parents ahead of children; keep exon before CDS, which is the conventional
# order and the one snpEff is tested against.
FEATURE_RANK = {"gene": 0, "transcript": 1, "mRNA": 1, "exon": 2, "CDS": 3,
                "start_codon": 4, "stop_codon": 5, "UTR": 6,
                "five_prime_utr": 6, "three_prime_utr": 6, "Selenocysteine": 7}

_attr_re = re.compile(r'(\w+)\s+"([^"]*)"')
# Extract single attributes rather than parsing the full attribute field of
# every line, to keep memory low on large GTFs.
_gid_re = re.compile(r'gene_id "([^"]*)"')
_tid_re = re.compile(r'transcript_id "([^"]*)"')
_eid_re = re.compile(r'exon_id "([^"]*)"')


def _grab(rx, s):
    m = rx.search(s)
    return m.group(1) if m else None


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def attrs(s):
    return dict(_attr_re.findall(s))


def comp(c):
    return "mito" if str(c) == MITO else "nuclear"


def hr(t):
    print("\n" + "=" * 76)
    print(t)
    print("=" * 76)


def pct(a, b):
    return f"{100.0*a/b:6.2f}%" if b else "   n/a"


def read_gtf(path):
    """Yield (fields, attrs) for every data line."""
    with opener(path) as fh:
        for ln, line in enumerate(fh, 1):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            yield ln, f, attrs(f[8])


# --------------------------------------------------------------------------- #
def cmd_inspect(args):
    hr(f"GTF structure: {args.gtf}")

    feat = defaultdict(Counter)
    tx_first = {}                # tid -> first feature type seen
    tx_has_tx = set()            # tid with an explicit transcript/mRNA row
    tx_comp = {}
    exon_id_tx = defaultdict(set)
    n_lines = 0

    with opener(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9:
                continue
            n_lines += 1
            c = comp(f[0])
            ftype = f[2]
            feat[c][ftype] += 1
            attr = f[8]
            tid = _grab(_tid_re, attr)
            if not tid:
                continue
            tx_comp[tid] = c
            if tid not in tx_first:
                tx_first[tid] = ftype
            if ftype in ("transcript", "mRNA"):
                tx_has_tx.add(tid)
            elif ftype == "exon":
                eid = _grab(_eid_re, attr)
                if eid:
                    exon_id_tx[eid].add(tid)

    print(f"  {n_lines:,} data lines, {len(tx_comp):,} transcripts\n")
    print("  feature counts by compartment:")
    print(f"    {'feature':<22s} {'nuclear':>12s} {'mito':>10s}")
    for k in sorted({k for c in feat for k in feat[c]}):
        print(f"    {k:<22s} {feat['nuclear'][k]:>12,} {feat['mito'][k]:>10,}")

    hr("1. Transcripts with a child row before their transcript row")
    bad, tot = Counter(), Counter()
    for tid, ftype in tx_first.items():
        c = tx_comp[tid]
        tot[c] += 1
        if tid in tx_has_tx and ftype not in ("transcript", "mRNA", "gene"):
            bad[c] += 1
    for c in ("nuclear", "mito"):
        print(f"  {c:8s}: {bad[c]:>7,} / {tot[c]:>7,}  {pct(bad[c], tot[c])}")

    hr("2. exon_id values shared by more than one transcript")
    dup = {e: t for e, t in exon_id_tx.items() if len(t) > 1}
    print(f"  {len(exon_id_tx):,} distinct exon_id values; {len(dup):,} shared")
    for e, t in list(dup.items())[:3]:
        print(f"    {e}  ->  {sorted(t)[:4]}")

    if args.assert_clean:
        problems = []
        if bad["nuclear"]:
            problems.append(f"{bad['nuclear']:,} nuclear transcripts have a child "
                            f"row before their 'transcript' row")
        if dup:
            problems.append(f"{len(dup):,} exon_id values shared across transcripts")
        if problems:
            print("\n" + "=" * 76)
            for p in problems:
                print(f"  FAIL: {p}")
            print("  snpEff would not mark these transcripts as protein-coding.")
            print("  Run `snpeff_gtf.py fix` first.")
            sys.exit(1)
        print("\n  ASSERT CLEAN: ordering and exon_ids are fine for snpEff.")


# --------------------------------------------------------------------------- #
def cmd_fix(args):
    """Rewrite the GTF grouped by transcript, parents first. Streams in one pass."""
    hr(f"Rewriting {args.gtf}\n       -> {args.out}")

    header = []
    groups = defaultdict(list)      # (contig, gene_id, tid) -> [(rank, start, end, line)]
    gorder = []
    seen = set()
    n_lines = 0

    with opener(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            n_lines += 1
            key = (f[0], _grab(_gid_re, f[8]) or "", _grab(_tid_re, f[8]) or "")
            if key not in seen:
                seen.add(key)
                gorder.append(key)
            groups[key].append((FEATURE_RANK.get(f[2], 9), int(f[3]), int(f[4]), f))

    print(f"  read {n_lines:,} data lines in {len(groups):,} transcript groups")

    n_reordered = n_exon_fixed = 0
    out = []
    for key in sorted(gorder, key=lambda k: (k[0], min(r[1] for r in groups[k]))):
        grp = groups[key]
        before = [r[3][2] for r in grp]
        grp.sort(key=lambda r: (r[0], r[1], r[2]))
        if before != [r[3][2] for r in grp]:
            n_reordered += 1
        if args.fix_exon_ids:
            n_ex = 0
            for _, _, _, f in grp:
                if f[2] == "exon" and "exon_id" in f[8]:
                    n_ex += 1
                    new = f"{key[2] or 'tx'}-E{n_ex}"
                    if _grab(_eid_re, f[8]) != new:
                        f[8] = _eid_re.sub(f'exon_id "{new}"', f[8])
                        n_exon_fixed += 1
        out.extend("\t".join(r[3]) + "\n" for r in grp)

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.writelines(header)
        fh.writelines(out)

    print(f"  transcript groups reordered : {n_reordered:,}")
    print(f"  exon_id values rewritten    : {n_exon_fixed:,}")
    print(f"  lines written               : {len(out):,}")


# --------------------------------------------------------------------------- #
def cmd_make_testvcf(args):
    """Emit a VCF of real mid-CDS positions, so the probe cannot be intergenic."""
    hr("Building a coding-position probe VCF")
    try:
        import pysam
    except ImportError:
        sys.exit("make-testvcf needs pysam (conda activate transcription_errors)")

    fa = pysam.FastaFile(args.fasta)
    cds = defaultdict(list)      # (contig, tx) -> [(start,end,strand)]
    for _, f, a in read_gtf(args.gtf):
        if f[2] != "CDS":
            continue
        tid = a.get("transcript_id", "?")
        cds[(f[0], tid)].append((int(f[3]), int(f[4])))

    print(f"  {len(cds):,} transcripts with CDS in the GTF")
    per = defaultdict(list)
    for (c, tid), segs in cds.items():
        segs.sort()
        # take a base comfortably inside the longest segment
        s, e = max(segs, key=lambda x: x[1] - x[0])
        if e - s < 12:
            continue
        per[comp(c)].append((c, (s + e) // 2, tid))

    recs = []
    for c in ("nuclear", "mito"):
        want = args.n_mito if c == "mito" else args.n_nuclear
        pool = per[c]
        step = max(1, len(pool) // max(want, 1))
        for c_, pos, tid in pool[::step][:want]:
            try:
                ref = fa.fetch(c_, pos - 1, pos).upper()
            except Exception:
                continue
            if ref not in "ACGT":
                continue
            alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
            recs.append((c_, pos, ref, alt, tid))
        print(f"  {c:8s}: {len([r for r in recs if comp(r[0])==c]):,} probe positions")

    order = {}
    with open(args.fasta + ".fai") as fh:
        for i, line in enumerate(fh):
            order[line.split("\t")[0]] = i
    recs.sort(key=lambda r: (order.get(r[0], 99), r[1]))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        with open(args.fasta + ".fai") as fai:
            for line in fai:
                c_, n = line.split("\t")[:2]
                fh.write(f"##contig=<ID={c_},length={n}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c_, pos, ref, alt, tid in recs:
            fh.write(f"{c_}\t{pos}\t.\t{ref}\t{alt}\t100\t.\t.\n")
    print(f"\n  wrote {len(recs):,} records -> {args.out}")


# --------------------------------------------------------------------------- #
def cmd_score(args):
    hr(f"Scoring {args.vcf}")
    n = Counter()
    coding = Counter()
    hgvs_c = Counter()
    terms = defaultdict(Counter)
    with opener(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or "ANN=" not in f[7]:
                continue
            c = comp(f[0])
            n[c] += 1
            ann = f[7].split("ANN=", 1)[1].split(";", 1)[0]
            hit = False
            for entry in ann.split(","):
                p = entry.split("|")
                if len(p) < 10:
                    continue
                for t in p[1].split("&"):
                    terms[c][t] += 1
                    if t in CODING_TERMS:
                        hit = True
                if p[9].startswith("c."):
                    hgvs_c[c] += 1
            if hit:
                coding[c] += 1

    for c in ("nuclear", "mito"):
        print(f"\n  --- {c} ---  {n[c]:,} annotated records")
        print(f"  with a coding effect : {coding[c]:>8,}  {pct(coding[c], n[c])}")
        print(f"  with c. HGVS         : {hgvs_c[c]:>8,}")
        for t, v in terms[c].most_common(5):
            mark = "  <-- coding" if t in CODING_TERMS else ""
            print(f"      {t:<45s} {v:>8,}{mark}")

    hr("VERDICT")
    ok = coding["nuclear"] > args.min_coding_frac * max(n["nuclear"], 1)
    if ok:
        print("  PASS - nuclear CDS positions receive coding annotations.")
    else:
        print(f"  FAIL - fewer than {100*args.min_coding_frac:.0f}% of nuclear CDS "
              f"probe positions received a coding annotation.")
        sys.exit(1)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("inspect")
    p.add_argument("--gtf", required=True)
    p.add_argument("--assert-clean", action="store_true",
                   help="Exit 1 if the GTF would give non-coding transcripts in snpEff.")
    p.set_defaults(fn=cmd_inspect)

    p = sub.add_parser("fix")
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--no-fix-exon-ids", dest="fix_exon_ids", action="store_false",
                   help="Only reorder rows; leave exon_id values unchanged.")
    p.set_defaults(fn=cmd_fix)

    p = sub.add_parser("make-testvcf")
    p.add_argument("--gtf", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n-nuclear", type=int, default=2000)
    p.add_argument("--n-mito", type=int, default=200)
    p.set_defaults(fn=cmd_make_testvcf)

    p = sub.add_parser("score")
    p.add_argument("--vcf", required=True)
    p.add_argument("--min-coding-frac", type=float, default=0.5,
                   help="Exit 1 if fewer than this fraction of nuclear probe "
                        "positions get a coding call.")
    p.set_defaults(fn=cmd_score)

    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
verify_ref_gtf.py - check that a GTF's CDS coordinates describe a FASTA.

Contig lengths alone cannot detect an annotation from a different assembly
release (Release 5 and 6 of D. melanogaster share contig names but differ in
coordinates). Instead, for each gene we take its longest CDS, splice the
segments in transcript order, trim the leading phase, and test:

  * is the length a multiple of 3?
  * does it begin with a start codon?
  * does it end with a stop codon (inside the CDS, or in the three bases
    directly after it, depending on the GTF's convention)?
  * is it free of internal stop codons?

Internal stops detect a large coordinate offset: in the correct frame they are
essentially absent, whereas a wrong frame hits a stop about every 20 codons.
A shift of one or two bases can occasionally pass that test, but then the
start- and stop-codon tests fail, so each test has its own threshold.
Start/stop get a lower threshold because real annotations contain some
deliberately incomplete CDS (Ensembl cds_start_NF / cds_end_NF).

The mitochondrial contig is translated with NCBI table 5 (invertebrate
mitochondrial) and reported as its own row. Several Drosophila mitochondrial
genes end in an incomplete stop codon (T or TA, completed to TAA by
polyadenylation), so a lower stop-codon fraction on the mito row is expected.
The thresholds are applied to genome-wide totals, which the nuclear genes
dominate.

Exits with status 2 if any test falls below its threshold.

Requires: pysam, biopython.
"""
import argparse
import re
import sys
from collections import defaultdict

import pysam
from Bio.Data import CodonTable

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
_attr_re = re.compile(r'(\w+) "([^"]*)"')

_TBL = {
    1: CodonTable.unambiguous_dna_by_id[1],
    5: CodonTable.unambiguous_dna_by_id[5],
}


def revcomp(s):
    return "".join(COMP.get(b, "N") for b in reversed(s))


def parse_attrs(field):
    return dict(_attr_re.findall(field))


def next_codon(fasta, chrom, strand, gcoord, contig_len):
    """
    The three bases immediately after the CDS, in transcript orientation.

    Under the GTF2.2 convention the terminal stop codon lives here rather than
    inside the CDS. Reading it is also a sharper test of the 3' boundary than
    the in-CDS check: it only reads as a stop if the final CDS coordinate is
    exactly right.
    """
    if not gcoord:
        return None
    last = gcoord[-1]                      # 1-based coord of the final CDS base
    if strand == "+":
        if contig_len and last + 3 > contig_len:
            return None
        return fasta.fetch(chrom, last, last + 3).upper()      # 0-based [last, last+3)
    if last - 3 < 1:
        return None
    seq = fasta.fetch(chrom, last - 4, last - 1).upper()       # 1-based last-3..last-1
    return revcomp(seq)


def classify(coding, table_id, next_codon=None):
    """
    Boolean checks for one spliced CDS.

    `next_codon` is the three bases immediately DOWNSTREAM of the CDS in
    transcript orientation, or None if they run off the contig.

    Two stop-codon conventions exist and both are legitimate:

      * GTF2.2 / Ensembl / FlyBase: the CDS feature EXCLUDES the terminal stop
        codon, which is annotated separately as a 'stop_codon' feature. The last
        codon of the spliced CDS is the final sense codon, and the stop sits in
        the three bases just after it.
      * Some GFF3-derived annotations INCLUDE the stop codon in the CDS.

    We test both and let the caller take whichever the GTF actually uses.
    Demanding the included-stop convention would fail a perfectly good Ensembl
    or FlyBase GTF at exactly 0.00%, which tells you nothing about whether the
    coordinates are right.
    """
    t = _TBL[table_id]
    out = {"len_mod3": (len(coding) % 3 == 0)}
    codons = [coding[i:i + 3] for i in range(0, len(coding) - len(coding) % 3, 3)]
    if not codons:
        return {"len_mod3": out["len_mod3"], "start_ok": False,
                "stop_in_cds": False, "stop_after_cds": False,
                "no_internal_stop": False, "n_codons": 0}
    out["start_ok"] = codons[0] in t.start_codons
    out["stop_in_cds"] = codons[-1] in t.stop_codons
    out["stop_after_cds"] = bool(next_codon) and (next_codon in t.stop_codons)
    # If the CDS ends in a stop, that codon is terminal and not "internal".
    # If it does not (stop-excluded convention), every codon must be checked -
    # skipping the last one there would let a real internal stop through.
    internal = codons[:-1] if out["stop_in_cds"] else codons
    out["no_internal_stop"] = not any(c in t.stop_codons for c in internal)
    out["n_codons"] = len(codons)
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Verify that a GTF's CDS coordinates match a FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--mito-contig", default="M")
    ap.add_argument("--min-ok-frac", type=float, default=0.90,
                    help="Fail if the fraction of representative CDS with no "
                         "internal stop falls below this")
    ap.add_argument("--min-start-stop-frac", type=float, default=0.80,
                    help="Fail if the fraction of representative CDS beginning "
                         "with a start codon, or ending with a stop codon, falls "
                         "below this. Lower than --min-ok-frac because real "
                         "annotations contain some incomplete (cds_start_NF / "
                         "cds_end_NF) transcripts.")
    ap.add_argument("--report", default=None, help="Optional per-contig TSV")
    ap.add_argument("--max-examples", type=int, default=5,
                    help="How many failing transcripts to print")
    args = ap.parse_args()

    fasta = pysam.FastaFile(args.fasta)
    fa_contigs = set(fasta.references)
    contig_len = dict(zip(fasta.references, fasta.lengths))

    # ---- collect CDS segments per transcript ----
    tx = {}
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            chrom, start, end, strand, frame, attrs = (
                f[0], int(f[3]), int(f[4]), f[6], f[7], f[8])
            a = parse_attrs(attrs)
            tid = a.get("transcript_id")
            gid = a.get("gene_id", tid)
            if tid is None:
                continue
            fr = int(frame) if frame in ("0", "1", "2") else 0
            rec = tx.setdefault(tid, {"gene": gid, "chrom": chrom,
                                      "strand": strand, "segs": []})
            rec["segs"].append((start, end, fr))

    if not tx:
        sys.exit("[VERIFY] FATAL: no CDS records with transcript_id in the GTF")

    # one representative transcript per gene: the longest CDS
    by_gene = defaultdict(list)
    for tid, rec in tx.items():
        rec["cds_len"] = sum(e - s + 1 for s, e, _ in rec["segs"])
        by_gene[rec["gene"]].append(tid)
    chosen = [sorted(t, key=lambda x: (-tx[x]["cds_len"], x))[0]
              for t in by_gene.values()]

    stats = defaultdict(lambda: defaultdict(int))
    failures = []

    for tid in chosen:
        rec = tx[tid]
        chrom, strand = rec["chrom"], rec["strand"]
        if chrom not in fa_contigs:
            stats[chrom]["contig_absent"] += 1
            continue
        table_id = 5 if chrom == args.mito_contig else 1

        segs = sorted(rec["segs"], key=lambda x: x[0], reverse=(strand == "-"))
        coding = []
        gcoord = []
        for (s1, e1, _fr) in segs:
            seq = fasta.fetch(chrom, s1 - 1, e1).upper()
            if strand == "+":
                coding.append(seq)
                gcoord.extend(range(s1, e1 + 1))
            else:
                coding.append(revcomp(seq))
                gcoord.extend(range(e1, s1 - 1, -1))
        coding = "".join(coding)
        lead = segs[0][2]
        if lead:
            coding = coding[lead:]
            gcoord = gcoord[lead:]
        if "N" in coding:
            stats[chrom]["has_N"] += 1

        nxt = next_codon(fasta, chrom, strand, gcoord, contig_len.get(chrom, 0))
        c = classify(coding, table_id, nxt)
        stats[chrom]["n"] += 1
        for k in ("len_mod3", "start_ok", "stop_in_cds", "stop_after_cds",
                  "no_internal_stop"):
            if c[k]:
                stats[chrom][k] += 1
        if not c["no_internal_stop"] and len(failures) < args.max_examples:
            failures.append((tid, rec["gene"], chrom, strand, c["n_codons"]))

    fasta.close()

    # ---- report ----
    header = ["contig", "n_representative_cds", "frac_len_mod3", "frac_start_ok",
              "frac_stop_in_cds", "frac_stop_after_cds", "frac_no_internal_stop",
              "n_contig_absent", "n_with_N"]
    rows = []
    tot_n = 0
    tot = {"no_internal_stop": 0, "start_ok": 0, "stop_in_cds": 0,
           "stop_after_cds": 0, "len_mod3": 0}
    for chrom in sorted(stats):
        s = stats[chrom]
        n = s["n"]
        tot_n += n
        for k in tot:
            tot[k] += s[k]
        frac = (lambda k: (s[k] / n) if n else float("nan"))
        rows.append([chrom, n, round(frac("len_mod3"), 4), round(frac("start_ok"), 4),
                     round(frac("stop_in_cds"), 4), round(frac("stop_after_cds"), 4),
                     round(frac("no_internal_stop"), 4),
                     s["contig_absent"], s["has_N"]])

    widths = [max(len(str(r[i])) for r in ([header] + rows)) for i in range(len(header))]
    def fmt(r):
        return "  ".join(str(v).ljust(w) for v, w in zip(r, widths))
    print("[VERIFY] " + fmt(header), file=sys.stderr)
    for r in rows:
        print("[VERIFY] " + fmt(r), file=sys.stderr)

    if args.report:
        with open(args.report, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(str(v) for v in r) + "\n")

    fracs = {k: (v / tot_n if tot_n else 0.0) for k, v in tot.items()}

    # Which stop-codon convention does this GTF use? Take whichever scores
    # higher; they are mutually exclusive in practice.
    if fracs["stop_after_cds"] >= fracs["stop_in_cds"]:
        stop_frac = fracs["stop_after_cds"]
        stop_n = tot["stop_after_cds"]
        convention = "stop codon EXCLUDED from CDS (GTF2.2 / Ensembl / FlyBase)"
    else:
        stop_frac = fracs["stop_in_cds"]
        stop_n = tot["stop_in_cds"]
        convention = "stop codon INCLUDED in CDS"

    print(f"[VERIFY] representative CDS checked : {tot_n}", file=sys.stderr)
    print(f"[VERIFY] stop-codon convention      : {convention}", file=sys.stderr)
    for k, label in (("no_internal_stop", "free of internal stops   "),
                     ("start_ok", "begin with a start codon "),
                     ("len_mod3", "length a multiple of 3    ")):
        print(f"[VERIFY] {label}: {tot[k]} ({100*fracs[k]:.2f}%)", file=sys.stderr)
    print(f"[VERIFY] terminated by a stop codon: {stop_n} ({100*stop_frac:.2f}%)",
          file=sys.stderr)

    if failures:
        print("[VERIFY] examples with internal stops:", file=sys.stderr)
        for tid, gid, chrom, strand, ncod in failures:
            print(f"[VERIFY]   {tid} ({gid}) {chrom}{strand} {ncod} codons",
                  file=sys.stderr)

    # each test gets its own floor; ANY of them failing condemns the pairing
    problems = []
    if fracs["no_internal_stop"] < args.min_ok_frac:
        problems.append(f"internal-stop-free {100*fracs['no_internal_stop']:.2f}% "
                        f"< {100*args.min_ok_frac:.0f}%")
    if fracs["start_ok"] < args.min_start_stop_frac:
        problems.append(f"start-codon {100*fracs['start_ok']:.2f}% "
                        f"< {100*args.min_start_stop_frac:.0f}%")
    if stop_frac < args.min_start_stop_frac:
        problems.append(f"stop-codon {100*stop_frac:.2f}% "
                        f"< {100*args.min_start_stop_frac:.0f}% "
                        f"(best of both conventions: in-CDS "
                        f"{100*fracs['stop_in_cds']:.2f}%, after-CDS "
                        f"{100*fracs['stop_after_cds']:.2f}%)")

    if problems:
        print("", file=sys.stderr)
        print("[VERIFY] FATAL: the GTF and the FASTA do not agree.", file=sys.stderr)
        for p in problems:
            print(f"[VERIFY]   {p}", file=sys.stderr)
        print("[VERIFY] A correct pairing gives ~99% on all tests. A failure on", file=sys.stderr)
        print("[VERIFY] start/stop alone usually means a small coordinate offset", file=sys.stderr)
        print("[VERIFY] (e.g. a 0-based/1-based convention error); a failure on", file=sys.stderr)
        print("[VERIFY] internal stops means the GTF and FASTA are from different", file=sys.stderr)
        print("[VERIFY] assembly releases.", file=sys.stderr)
        sys.exit(2)

    print("[VERIFY] OK - annotation and sequence agree.", file=sys.stderr)


if __name__ == "__main__":
    main()

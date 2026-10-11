#!/usr/bin/env python3
"""
Stage 3 (summary) - per-sample summary table and cohort overview
=================================================================
Runs one filter setting end to end: builds the cohort site table and applies
the filter (stage3_filters), computes per-sample N/S sites (stage3_dnds),
counts the retained mutations and sums each sample's depth over the analysed
region, and writes the four output files below.

Regions
-------
Every count is made over the analysed region (positions callable in at least
--min-callable samples, minus positions where every observed allele was
disqualified) intersected with the sample's own callable intervals.

  all          the whole analysed region
  genespan     positions in --genes-bed
  cds          positions in --cds-bed
  synonymous   CDS positions weighted by s_site, the Nei-Gojobori synonymous
               fraction (0, 1/3, 2/3 or 1); counts in this region are sums of
               fractional site weights, not integers

Nuclear and mitochondrial contigs are reported separately; --mito-contig names
the latter.

Inputs
------
--stage2b-dir     <sample>.annotated.parquet, columns sample_id, chrom, pos,
                  ref, alt, type, ao, ro, dp, qa, mean_alt_qual, vaf_dp,
                  vaf_ar, ref_fasta_match, consequence, gene_id,
                  transcript_id, effect, impact, gene_strand, spectrum_channel
--cohort-dir      cohort_sum_depth.<contig>.npy, cohort_n_callable.<contig>.npy
--per-sample-dir  <sample>.depth.npz (arrays <contig>__pos, <contig>__dep),
                  <sample>.callable.bed
--ns-table        ns_sites.parquet: chrom, pos, n_site, s_site, gene_id,
                  transcript_id. Checked against the channel table below.
--ns-channel      ns_sites.channels.parquet: chrom, pos, channel, is_syn,
                  weight. The source of N, S and the synonymous depth weights.
--fasta, --fai    canonical reference and its index
--genes-bed       gene spans (genes.bed)
--cds-bed         CDS positions (ns_sites.cds.bed)
--samples-list    one sample ID per line
--filter, --params, --min-depth, --min-callable, --mito-contig, --out-prefix

Outputs
-------
<out-prefix>.summary.tsv - one row per sample. sample_id, then the columns
below for the nuclear genome and again with a '_mito' suffix:

  depth_sum_total, depth_sum_genespan, depth_sum_cds, depth_sum_noncds
                            summed read depth over the region
  depth_sum_synonymous      s_site-weighted depth over CDS positions
  positions_covered, positions_covered_genespan, positions_covered_cds,
  positions_covered_noncds  positions at depth >= --min-depth
  total_mutations, genespan_mutations, cds_mutations, synonymous_mutations,
  missense_mutations        alt reads (summed AO) of the retained SNVs
  total_alleles, genespan_alleles, cds_alleles
                            distinct (pos, ref, alt) of the same SNVs
  reads_per_allele          total_mutations / total_alleles
  excluded_mutations, excluded_alleles
                            retained SNVs outside the analysed region or
                            outside this sample's callable intervals
  rate_total_per_total_depth, rate_genespan_per_genespan_depth,
  rate_cds_per_cds_depth, rate_synonymous_per_syn_depth
                            each numerator over its matching denominator
  N_sites, S_sites          Nei-Gojobori potential sites, over the changes
                            the filter left observable
  pN, pS, dNdS              Dn/N, Ds/S and their ratio
  N_weighted, S_weighted, dNdS_weighted
                            the same, weighted by the sample's CDS spectrum
  spectrum_skipped, spectrum_skipped_pct, spectrum_skipped_cds
                            alt reads with no spectrum channel (no transcribed
                            strand at that position)
  spec_<channel>            alt reads per channel over gene spans; pairs with
                            depth_sum_genespan. The 12 channels are
                            A>C A>G A>U C>A C>G C>U G>A G>C G>U U>A U>C U>G,
                            written in transcript orientation.
                            sum(spec_*) == total_mutations - spectrum_skipped
  spec_<channel>_cds        the same channels restricted to CDS positions;
                            pairs with depth_sum_cds and with N_c/S_c.
                            sum(spec_*_cds) == cds_mutations
                                               - spectrum_skipped_cds

<out-prefix>.overview.tsv - one row for the cohort: filter, params,
mito_contig, n_samples, min_callable, then the statistics below for each
compartment x region, suffixed '' / '_genespan' / '_cds' / '_synonymous' for
the nuclear genome and '_mito' / '_mito_genespan' / '_mito_cds' /
'_mito_synonymous' for the mitochondrion:

  total_positions, total_read_depth_cohort, total_read_depth_per_sample_mean,
  avg_depth_per_position_cohort, avg_depth_per_position_per_sample,
  total_mutations, total_alleles, reads_per_allele, avg_mutations_per_sample,
  avg_mutations_per_MB, mutations_per_read_depth,
  positions_with_A, positions_with_C, positions_with_G, positions_with_U
                            reference bases over the region, on the reference
                            strand (not transcript orientation); T is reported
                            as U

<out-prefix>.overview_long.tsv - the same statistics with one row per
compartment x region: filter, compartment, region, then the unsuffixed
statistic names.

<out-prefix>.params.json - filter, params, mito_contig, min_depth, n_samples,
min_callable.

Requires: numpy, pandas, pyarrow, pysam. Imports stage3_filters, stage3_dnds,
contig_names.
"""
import argparse
import json
import os
from collections import Counter

import numpy as np
import pandas as pd
import pysam

import stage3_filters as F
import stage3_dnds as D
import contig_names as CN

COMPARTMENTS = ("nuclear", "mito")
REGIONS = ("all", "genespan", "cds", "synonymous")

# column-name suffix for each (compartment, region) cell of the overview
SUFFIX = {
    ("nuclear", "all"): "",
    ("nuclear", "genespan"): "_genespan",
    ("nuclear", "cds"): "_cds",
    ("nuclear", "synonymous"): "_synonymous",
    ("mito", "all"): "_mito",
    ("mito", "genespan"): "_mito_genespan",
    ("mito", "cds"): "_mito_cds",
    ("mito", "synonymous"): "_mito_synonymous",
}


def load_fai(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            sizes[f[0]] = int(f[1])
    return sizes


def load_bed_intervals(path):
    """{contig: [(start0,end), ...]} from a BED file."""
    d = {}
    if not os.path.exists(path):
        return d
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            d.setdefault(CN.canonicalise(f[0]), []).append((int(f[1]), int(f[2])))
    return d


def in_intervals(ivs, pos0):
    """Boolean mask: which 0-based positions fall inside sorted [start,end)."""
    if not ivs or pos0.size == 0:
        return np.zeros(pos0.size, dtype=bool)
    iv = np.asarray(sorted(ivs), dtype=np.int64)
    starts, ends = iv[:, 0], iv[:, 1]
    idx = np.searchsorted(starts, pos0, side="right") - 1
    ok = idx >= 0
    out = np.zeros(pos0.size, dtype=bool)
    out[ok] = pos0[ok] < ends[idx[ok]]
    return out


def load_bed_mask(path, sizes):
    mask = {c: np.zeros(n, dtype=bool) for c, n in sizes.items()}
    for c, ivs in load_bed_intervals(path).items():
        if c in mask:
            for a, b in ivs:
                mask[c][a:b] = True
    return mask


def build_syn_weight(ns_channel, sizes, disqualified=None):
    """
    {contig: float64 array} of synonymous site weight per 0-based position,
    0 outside CDS.

    Summed from the per-change rows of Stage 0C, excluding the changes the
    filter disqualified. With nothing disqualified this reproduces the s_site
    column of ns_sites.parquet.

    float64 (~1.2 GB for this genome) because float32 would make a whole-number
    site total read as 24.00000017.
    """
    w = {c: np.zeros(n, dtype=np.float64) for c, n in sizes.items()}
    dq_by_c = ({c: g for c, g in disqualified.groupby("chrom")}
               if disqualified is not None and len(disqualified) else {})
    n_set = 0
    for c, g in ns_channel.groupby("chrom"):
        c = CN.canonicalise(c)
        if c not in w:
            continue
        p0 = g["pos"].to_numpy() - 1
        cidx = np.array([D.CH_IDX.get(ch, -1) for ch in g["channel"]], dtype=np.int64)
        sel = ((p0 >= 0) & (p0 < w[c].size) & (cidx >= 0)
               & g["is_syn"].to_numpy().astype(bool)
               & D.live_mask(p0, cidx, dq_by_c.get(c)))
        np.add.at(w[c], p0[sel], g["weight"].to_numpy().astype(np.float64)[sel])
        n_set += int(np.count_nonzero(w[c]))
    return w, n_set


def build_analysed_mask(cohort_dir, sizes, fully_disq, min_nc):
    """Genome-wide analysed region: n_callable >= min_nc, minus fully-disqualified
    positions. min_nc is the same threshold as the filters' callability_pass."""
    mask = {}
    for contig, n in sizes.items():
        nc_path = os.path.join(cohort_dir, f"cohort_n_callable.{contig}.npy")
        if not os.path.exists(nc_path):
            raise FileNotFoundError(
                f"{nc_path} is missing; run Stage 2A reduce for every contig in "
                f"the .fai before Stage 3")
        nc = np.load(nc_path)
        if nc.size < n:
            nc = np.concatenate([nc, np.zeros(n - nc.size, dtype=nc.dtype)])
        mask[contig] = nc[:n] >= min_nc
    for (c, p1) in fully_disq:
        c = CN.canonicalise(c)
        if c in mask and 0 <= p1 - 1 < mask[c].size:
            mask[c][p1 - 1] = False
    return mask


# ---------------------------------------------------------------------------
# cohort-level stats, per compartment x region
# ---------------------------------------------------------------------------
def cohort_region_stats(fasta_path, cohort_dir, analysed, genespan_mask, cds_mask,
                        syn_weight, sizes, mito_contig):
    """
    -> {(compartment, region): {"positions": float, "depth": float,
                                "base": Counter of reference base -> weight}}

    'positions' and the base counts are integer counts for the 'all' and
    'cds' regions and sums of fractional s_site weights for 'synonymous'.
    'depth' is the cohort-summed read depth over the region, likewise weighted
    for 'synonymous'.
    """
    fa = pysam.FastaFile(fasta_path)
    fa_refs = set(fa.references)
    out = {(c, r): {"positions": 0.0, "depth": 0.0, "base": Counter()}
           for c in COMPARTMENTS for r in REGIONS}

    for contig, n in sizes.items():
        m = analysed.get(contig)
        if m is None or not m.any():
            continue
        comp = "mito" if contig == mito_contig else "nuclear"
        gsm = genespan_mask.get(contig)
        gsm = np.zeros(n, dtype=bool) if gsm is None else gsm
        gm = cds_mask.get(contig)
        gm = np.zeros(n, dtype=bool) if gm is None else gm
        sw = syn_weight.get(contig)
        sw = np.zeros(n, dtype=np.float64) if sw is None else sw

        sd_path = os.path.join(cohort_dir, f"cohort_sum_depth.{contig}.npy")
        if not os.path.exists(sd_path):
            raise FileNotFoundError(
                f"{sd_path} is missing; run Stage 2A reduce for every contig in "
                f"the .fai before Stage 3")
        sd = np.load(sd_path).astype(np.float64)
        if sd.size < n:
            sd = np.concatenate([sd, np.zeros(n - sd.size)])
        sd = sd[:n]

        sel_all = m
        sel_gsp = m & gsm
        sel_cds = m & gm
        # synonymous is a weight vector rather than a boolean selector: every
        # analysed CDS base contributes s_site of a synonymous site
        w_syn = np.where(m, sw, 0.0)

        out[(comp, "all")]["positions"] += float(sel_all.sum())
        out[(comp, "all")]["depth"] += float(sd[sel_all].sum())
        out[(comp, "genespan")]["positions"] += float(sel_gsp.sum())
        out[(comp, "genespan")]["depth"] += float(sd[sel_gsp].sum())
        out[(comp, "cds")]["positions"] += float(sel_cds.sum())
        out[(comp, "cds")]["depth"] += float(sd[sel_cds].sum())
        out[(comp, "synonymous")]["positions"] += float(w_syn.sum())
        out[(comp, "synonymous")]["depth"] += float((sd * w_syn).sum())

        if contig in fa_refs:
            seq = np.frombuffer(fa.fetch(contig).upper().encode(), dtype="S1")
            L = min(n, seq.size)
            for b in ("A", "C", "G", "T"):
                bm = seq[:L] == b.encode()
                out[(comp, "all")]["base"][b] += float((bm & sel_all[:L]).sum())
                out[(comp, "genespan")]["base"][b] += float((bm & sel_gsp[:L]).sum())
                out[(comp, "cds")]["base"][b] += float((bm & sel_cds[:L]).sum())
                out[(comp, "synonymous")]["base"][b] += float(w_syn[:L][bm].sum())

    fa.close()
    return out


# ---------------------------------------------------------------------------
# per-sample
# ---------------------------------------------------------------------------
def _blank_depth():
    return {"depth_total": 0, "depth_genespan": 0, "depth_cds": 0, "depth_noncds": 0,
            "depth_synonymous": 0.0,
            "cov_total": 0, "cov_genespan": 0, "cov_cds": 0, "cov_noncds": 0}


def per_sample_depth(npz_path, analysed_mask, genespan_mask, cds_mask, syn_weight,
                     min_depth, mito_contig):
    """
    Sum sample depth and count covered positions within the analysed region,
    split by compartment, and accumulate the s_site-weighted synonymous depth.
    """
    out = {c: _blank_depth() for c in COMPARTMENTS}
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"{npz_path} is missing; rerun Stage 2A for this sample")
    with np.load(npz_path) as z:
        contigs = sorted({k.split("__")[0] for k in z.files})
        for c in contigs:
            if f"{c}__pos" not in z.files:
                continue
            pos = z[f"{c}__pos"]
            dep = z[f"{c}__dep"].astype(np.int64)
            am = analysed_mask.get(c)
            if am is None:
                continue
            o = out["mito" if c == mito_contig else "nuclear"]
            gsm = genespan_mask.get(c)
            gm = cds_mask.get(c)
            sw = syn_weight.get(c)

            in_an = am[pos]
            pos, dep = pos[in_an], dep[in_an]
            if pos.size == 0:
                continue
            gsp = gsm[pos] if gsm is not None else np.zeros(pos.size, bool)
            cds = gm[pos] if gm is not None else np.zeros(pos.size, bool)
            o["depth_total"] += int(dep.sum())
            o["depth_genespan"] += int(dep[gsp].sum())
            o["depth_cds"] += int(dep[cds].sum())
            o["depth_noncds"] += int(dep[~cds].sum())
            if sw is not None:
                # each read at this base counts as s_site(base) of a
                # potentially-synonymous observation
                o["depth_synonymous"] += float(
                    np.dot(dep.astype(np.float64), sw[pos]))
            cov = dep >= min_depth
            o["cov_total"] += int(cov.sum())
            o["cov_genespan"] += int((cov & gsp).sum())
            o["cov_cds"] += int((cov & cds).sum())
            o["cov_noncds"] += int((cov & ~cds).sum())
    return out


def _mask_lookup(mask, chrom, pos, mito_contig):
    """Boolean array: which rows fall inside a {contig: bool array} mask."""
    out = np.zeros(chrom.size, dtype=bool)
    for c in pd.unique(chrom):
        gm = mask.get(CN.canonicalise(c, mito_contig))
        if gm is None:
            continue
        sel = np.flatnonzero(chrom == c)
        p0 = pos[sel] - 1
        ok = (p0 >= 0) & (p0 < gm.size)
        out[sel[ok]] = gm[p0[ok]]
    return out


def flag_rows(df, genespan_mask, cds_mask, mito_contig):
    """(is_mito, is_genespan, is_cds) boolean arrays for a mutation table."""
    if len(df) == 0:
        raise ValueError(
            "mutation table is empty: no SNV was retained by this filter, "
            "which should not happen on real data - check the Stage 2B "
            "annotated tables and the filter parameters")
    chrom = df["chrom"].to_numpy()
    pos = df["pos"].to_numpy()
    # canonicalise once per distinct contig rather than once per row
    canon = {c: CN.canonicalise(c, mito_contig) for c in pd.unique(chrom)}
    is_mito = (pd.Series(chrom).map(canon) == mito_contig).to_numpy()
    gsp = _mask_lookup(genespan_mask, chrom, pos, mito_contig)
    cds = _mask_lookup(cds_mask, chrom, pos, mito_contig)
    return is_mito, gsp, cds


def _blank_counts():
    d = {"total_mutations": 0, "total_alleles": 0,
         "genespan_mutations": 0, "genespan_alleles": 0,
         "cds_mutations": 0, "cds_alleles": 0,
         "synonymous_mutations": 0, "missense_mutations": 0,
         "spectrum_skipped": 0, "spectrum_skipped_cds": 0,
         "excluded_mutations": 0, "excluded_alleles": 0}
    # two spectra over the same mutations, differing only in region:
    #   spec_<ch>       pairs with depth_genespan
    #   spec_<ch>_cds   pairs with depth_cds
    d.update({f"spec_{ch}": 0 for ch in D.CHANNELS})
    d.update({f"spec_{ch}_cds": 0 for ch in D.CHANNELS})
    d["_chan_array"] = np.zeros(12)          # gene-span spectrum
    d["_chan_array_cds"] = np.zeros(12)      # CDS spectrum (used for dN/dS)
    return d


def analysable_rows(df, analysed_mask, callable_ivs, mito_contig):
    """
    Boolean mask: rows that sit both inside the filter's analysed region and
    inside this sample's own callable (depth >= min_depth) intervals, so that
    the mutation counts cover exactly the positions per_sample_depth() sums.
    freebayes called at raw DP >= 100 while Stage 2A measures depth at BQ >= 73,
    so without this a position could be in the numerator while contributing
    nothing to the denominator.
    """
    keep = np.zeros(len(df), dtype=bool)
    if len(df) == 0:
        return keep
    chrom = df["chrom"].to_numpy()
    pos = df["pos"].to_numpy()
    for c in pd.unique(chrom):
        cc = CN.canonicalise(c, mito_contig)
        am = analysed_mask.get(cc)
        if am is None:
            continue
        sel = np.flatnonzero(chrom == c)
        p0 = pos[sel] - 1
        ok = (p0 >= 0) & (p0 < am.size)
        m = np.zeros(sel.size, dtype=bool)
        m[ok] = am[p0[ok]]
        m &= in_intervals(callable_ivs.get(cc, []), p0)
        keep[sel] = m
    return keep


def per_sample_counts(sample_ret, genespan_mask, cds_mask, mito_contig,
                      analysed_mask, callable_ivs):
    """
    AO-weighted mutation counts and 12-channel spectrum for one sample's
    retained SNVs, split into nuclear and mitochondrial compartments.

    Counts are sums of AO (alt reads), matching the read-depth denominators.
    Distinct site-alleles are counted alongside as total_alleles / cds_alleles.
    """
    out = {c: _blank_counts() for c in COMPARTMENTS}
    if len(sample_ret) == 0:
        return out

    inreg = analysable_rows(sample_ret, analysed_mask, callable_ivs, mito_contig)
    is_mito_all, is_gsp_all, is_cds_all = flag_rows(
        sample_ret, genespan_mask, cds_mask, mito_contig)
    ao_all = sample_ret["ao"].to_numpy()

    for comp in COMPARTMENTS:
        comp_sel = is_mito_all if comp == "mito" else ~is_mito_all
        d = out[comp]

        # what the region restriction removed, reported rather than hidden
        dropped = comp_sel & ~inreg
        d["excluded_alleles"] = int(dropped.sum())
        d["excluded_mutations"] = int(ao_all[dropped].sum())

        sel = comp_sel & inreg
        if not sel.any():
            continue
        sub = sample_ret.loc[sel]
        ao = ao_all[sel]
        gsp = is_gsp_all[sel]
        cds = is_cds_all[sel]
        is_syn = (sub["consequence"] == "synonymous").to_numpy()
        is_mis = (sub["consequence"] == "missense").to_numpy()
        spec = sub["spectrum_channel"]
        no_chan = spec.isna().to_numpy()

        d["total_mutations"] = int(ao.sum())
        d["total_alleles"] = int(len(sub))
        d["genespan_mutations"] = int(ao[gsp].sum())
        d["genespan_alleles"] = int(gsp.sum())
        d["cds_mutations"] = int(ao[cds].sum())
        d["cds_alleles"] = int(cds.sum())
        d["synonymous_mutations"] = int(ao[is_syn].sum())
        d["missense_mutations"] = int(ao[is_mis].sum())
        d["spectrum_skipped"] = int(ao[no_chan].sum())
        d["spectrum_skipped_cds"] = int(ao[no_chan & cds].sum())

        # gene-span spectrum: every mutation assigned a transcribed strand;
        # sums to total_mutations - spectrum_skipped
        spec_np = spec.to_numpy()
        chan = {ch: int(ao[spec_np == ch].sum()) for ch in D.CHANNELS}
        # CDS spectrum: the same channels restricted to coding positions;
        # sums to cds_mutations - spectrum_skipped_cds
        chan_cds = {ch: int(ao[(spec_np == ch) & cds].sum()) for ch in D.CHANNELS}

        for ch in D.CHANNELS:
            d[f"spec_{ch}"] = chan[ch]
            d[f"spec_{ch}_cds"] = chan_cds[ch]
        d["_chan_array"] = np.array([chan[c] for c in D.CHANNELS])
        d["_chan_array_cds"] = np.array([chan_cds[c] for c in D.CHANNELS])
    return out


def _sample_block(counts, dpth, sites, suffix):
    """Assemble one compartment's worth of per-sample columns."""
    # Dn and Ds are AO-weighted, matching total_mutations. N and S are
    # potential-SITE counts, so pN and pS are read counts over site counts;
    # the ratio pN/pS is unaffected as long as reads-per-allele does not differ
    # systematically between missense and synonymous alleles.
    # The weighting spectrum is the CDS one: N_c and S_c are defined on CDS
    # positions, so a gene-span spectrum (intronic and UTR mutations included)
    # would mix regions inside one ratio.
    dd = D.dnds_from_components(counts["missense_mutations"],
                                counts["synonymous_mutations"],
                                sites["N"], sites["S"],
                                sites["N_c"], sites["S_c"],
                                counts["_chan_array_cds"], sites["opp_c"])
    tot = counts["total_mutations"] or 1
    syn_depth = dpth["depth_synonymous"]
    block = {
        f"depth_sum_total{suffix}": dpth["depth_total"],
        f"depth_sum_genespan{suffix}": dpth["depth_genespan"],
        f"depth_sum_cds{suffix}": dpth["depth_cds"],
        f"depth_sum_noncds{suffix}": dpth["depth_noncds"],
        # the denominator for a synonymous error rate
        f"depth_sum_synonymous{suffix}": round(syn_depth, 3),
        f"positions_covered{suffix}": dpth["cov_total"],
        f"positions_covered_genespan{suffix}": dpth["cov_genespan"],
        f"positions_covered_cds{suffix}": dpth["cov_cds"],
        f"positions_covered_noncds{suffix}": dpth["cov_noncds"],
        # AO-weighted counts (reads)
        f"total_mutations{suffix}": counts["total_mutations"],
        f"genespan_mutations{suffix}": counts["genespan_mutations"],
        f"cds_mutations{suffix}": counts["cds_mutations"],
        f"synonymous_mutations{suffix}": counts["synonymous_mutations"],
        f"missense_mutations{suffix}": counts["missense_mutations"],
        # distinct site-alleles, kept as a duplication / heteroplasmy diagnostic
        f"total_alleles{suffix}": counts["total_alleles"],
        f"genespan_alleles{suffix}": counts["genespan_alleles"],
        f"cds_alleles{suffix}": counts["cds_alleles"],
        f"reads_per_allele{suffix}":
            round(counts["total_mutations"] / counts["total_alleles"], 3)
            if counts["total_alleles"] else None,
        # rows dropped for sitting outside the analysed region or outside this
        # sample's callable intervals
        f"excluded_mutations{suffix}": counts["excluded_mutations"],
        f"excluded_alleles{suffix}": counts["excluded_alleles"],
        # convenience rates, so nobody has to remember which denominator pairs
        # with which numerator
        f"rate_total_per_total_depth{suffix}":
            (counts["total_mutations"] / dpth["depth_total"]) if dpth["depth_total"] else None,
        f"rate_genespan_per_genespan_depth{suffix}":
            (counts["genespan_mutations"] / dpth["depth_genespan"]) if dpth["depth_genespan"] else None,
        f"rate_cds_per_cds_depth{suffix}":
            (counts["cds_mutations"] / dpth["depth_cds"]) if dpth["depth_cds"] else None,
        f"rate_synonymous_per_syn_depth{suffix}":
            (counts["synonymous_mutations"] / syn_depth) if syn_depth > 0 else None,
        f"N_sites{suffix}": round(sites["N"], 3),
        f"S_sites{suffix}": round(sites["S"], 3),
        f"pN{suffix}": dd["pN"], f"pS{suffix}": dd["pS"], f"dNdS{suffix}": dd["dNdS"],
        f"N_weighted{suffix}": round(dd["N_w"], 3),
        f"S_weighted{suffix}": round(dd["S_w"], 3),
        f"dNdS_weighted{suffix}": dd["dNdS_weighted"],
        f"spectrum_skipped{suffix}": counts["spectrum_skipped"],
        f"spectrum_skipped_pct{suffix}": round(100 * counts["spectrum_skipped"] / tot, 2),
        f"spectrum_skipped_cds{suffix}": counts["spectrum_skipped_cds"],
    }
    block.update({f"spec_{ch}{suffix}": counts[f"spec_{ch}"] for ch in D.CHANNELS})
    block.update({f"spec_{ch}_cds{suffix}": counts[f"spec_{ch}_cds"] for ch in D.CHANNELS})
    return block


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--stage2b-dir", required=True)
    ap.add_argument("--cohort-dir", required=True)
    ap.add_argument("--per-sample-dir", required=True, help="stage2a per_sample dir")
    ap.add_argument("--ns-table", required=True)
    ap.add_argument("--ns-channel", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--genes-bed", required=True,
                    help="gene spans from Stage 0A (genes.bed); the denominator "
                         "region for the unsuffixed spec_* spectrum")
    ap.add_argument("--cds-bed", required=True,
                    help="CDS positions from Stage 0C (ns_sites.cds.bed); the "
                         "region the cds_* and synonymous columns share")
    ap.add_argument("--samples-list", required=True)
    ap.add_argument("--filter", required=True, choices=list(F.FILTERS))
    ap.add_argument("--params", required=True, help="JSON dict of thresholds")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--min-callable", type=int, default=1,
                    help="Minimum number of samples with depth >= --min-depth at a "
                         "position for it to be analysed")
    ap.add_argument("--mito-contig", default=CN.CANONICAL_MITO,
                    help="Contig reported separately from the nuclear genome; "
                         "its columns carry a '_mito' suffix")
    args = ap.parse_args()

    params = json.loads(args.params)
    sizes = load_fai(args.fai)
    samples = [l.strip() for l in open(args.samples_list) if l.strip()]
    genespan_mask = load_bed_mask(args.genes_bed, sizes)
    cds_mask = load_bed_mask(args.cds_bed, sizes)
    mito = args.mito_contig
    if mito not in sizes:
        print(f"[STAGE3] WARNING: mito contig {mito!r} is not in {args.fai}; "
              f"all '_mito' columns will be zero.")

    # --- filter ---
    mut = F.load_mutations(args.stage2b_dir, samples)
    n_samples = len(samples)

    site = F.build_site_table(mut, args.cohort_dir)
    labeled = F.apply_filter(site, args.filter, params, args.min_callable)
    fully_disq = F.fully_disqualified_positions(labeled)
    dq_changes = F.disqualified_changes(labeled)
    ret = F.retained_mutations(mut, labeled)

    analysed = build_analysed_mask(args.cohort_dir, sizes, fully_disq,
                                   args.min_callable)

    # --- per-change N/S opportunity, less the disqualified changes ---
    contig_groups = {c: ("mito" if c == mito else "nuclear") for c in sizes}
    ns_chan = pd.read_parquet(args.ns_channel)

    # the channel table must agree with the n_site/s_site columns
    nss = pd.read_parquet(args.ns_table, columns=["n_site", "s_site"])
    tot_chan = float(ns_chan["weight"].sum())
    tot_pos = float(nss["n_site"].sum() + nss["s_site"].sum())
    if not np.isclose(tot_chan, tot_pos, rtol=1e-9):
        raise SystemExit(
            f"[STAGE3] FATAL: {args.ns_channel} sums to {tot_chan:.4f} of "
            f"N+S opportunity but {args.ns_table} sums to {tot_pos:.4f}; "
            f"rerun Stage 0C so both files come from the same build")
    del nss

    syn_weight, n_cds_sites = build_syn_weight(ns_chan, sizes, dq_changes)
    print(f"[STAGE3] synonymous weights laid out for {n_cds_sites} CDS positions")
    print(f"[STAGE3] changes disqualified by the filter: {len(dq_changes)}")
    callable_paths = {s: os.path.join(args.per_sample_dir, f"{s}.callable.bed")
                      for s in samples}
    missing = [p for p in callable_paths.values() if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} sample(s) have no callable.bed, e.g. {missing[:5]}; "
            f"rerun Stage 2A for them")
    callable_by_sample = {s: load_bed_intervals(p) for s, p in callable_paths.items()}
    site_acc = D.compute_sites_per_sample(ns_chan, sizes, analysed,
                                          callable_by_sample, contig_groups,
                                          dq_changes)
    blank_sites = {"N": 0.0, "S": 0.0, "N_c": np.zeros(12),
                   "S_c": np.zeros(12), "opp_c": np.zeros(12)}

    # --- per-sample rows ---
    ret_by_sample = {s: g for s, g in ret.groupby("sample_id")}
    rows = []
    for s in samples:
        sr = ret_by_sample.get(s, ret.iloc[0:0])
        c = per_sample_counts(sr, genespan_mask, cds_mask, mito, analysed,
                              callable_by_sample.get(s, {}))
        dpth = per_sample_depth(os.path.join(args.per_sample_dir, f"{s}.depth.npz"),
                                analysed, genespan_mask, cds_mask, syn_weight,
                                args.min_depth, mito)
        row = {"sample_id": s}
        for comp in COMPARTMENTS:
            suffix = "" if comp == "nuclear" else "_mito"
            row.update(_sample_block(c[comp], dpth[comp],
                                     site_acc[s].get(comp, blank_sites), suffix))
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.out_prefix}.summary.tsv", sep="\t", index=False)

    # --- overview ---
    cstats = cohort_region_stats(args.fasta, args.cohort_dir, analysed,
                                 genespan_mask, cds_mask, syn_weight, sizes, mito)

    ret = ret.reset_index(drop=True)
    is_mito_ret, is_gsp_ret, is_cds_ret = flag_rows(ret, genespan_mask, cds_mask, mito)
    is_syn_ret = (ret["consequence"] == "synonymous").to_numpy()
    ao_ret = ret["ao"].to_numpy()

    # Same restriction as the per-sample rows: a mutation counts only where its
    # own sample contributes depth to the denominator. Without it the overview
    # numerator would cover a wider region than cohort_region_stats() sums, and
    # would not match the sum of the per-sample columns.
    in_region_ret = np.zeros(len(ret), dtype=bool)
    for smp, g in ret.groupby("sample_id", sort=False):
        in_region_ret[g.index.to_numpy()] = analysable_rows(
            g, analysed, callable_by_sample.get(smp, {}), mito)

    def _region_sel(comp, region):
        sel = (is_mito_ret if comp == "mito" else ~is_mito_ret) & in_region_ret
        if region == "genespan":
            sel = sel & is_gsp_ret
        elif region == "cds":
            sel = sel & is_cds_ret
        elif region == "synonymous":
            # the numerator that pairs with the s_site-weighted denominator
            sel = sel & is_syn_ret
        return sel

    def mutation_count(comp, region):
        """AO-weighted, matching the per-sample columns."""
        return int(ao_ret[_region_sel(comp, region)].sum())

    def allele_count(comp, region):
        return int(_region_sel(comp, region).sum())

    overview = {"filter": args.filter, "params": json.dumps(params),
                "mito_contig": mito, "n_samples": n_samples,
                "min_callable": args.min_callable}
    long_rows = []
    for comp in COMPARTMENTS:
        for region in REGIONS:
            st = cstats[(comp, region)]
            sfx = SUFFIX[(comp, region)]
            pos = st["positions"]
            dep = st["depth"]
            nmut = mutation_count(comp, region)
            nall = allele_count(comp, region)
            mb = pos / 1e6 if pos else 0.0
            stats = {
                "total_positions": pos,
                "total_read_depth_cohort": dep,
                "total_read_depth_per_sample_mean": dep / n_samples if n_samples else None,
                "avg_depth_per_position_cohort": (dep / pos) if pos else None,
                "avg_depth_per_position_per_sample": (dep / pos / n_samples) if (pos and n_samples) else None,
                "total_mutations": nmut,
                "total_alleles": nall,
                "reads_per_allele": (nmut / nall) if nall else None,
                "avg_mutations_per_sample": nmut / n_samples if n_samples else None,
                "avg_mutations_per_MB": (nmut / mb) if mb else None,
                "mutations_per_read_depth": (nmut / dep) if dep else None,
                "positions_with_C": st["base"].get("C", 0.0),
                "positions_with_G": st["base"].get("G", 0.0),
                "positions_with_A": st["base"].get("A", 0.0),
                "positions_with_U": st["base"].get("T", 0.0),
            }
            for k, v in stats.items():
                overview[f"{k}{sfx}"] = v
            long_rows.append({"filter": args.filter,
                              "compartment": comp,
                              "region": region, **stats})

    pd.DataFrame([overview]).to_csv(f"{args.out_prefix}.overview.tsv", sep="\t", index=False)
    pd.DataFrame(long_rows).to_csv(f"{args.out_prefix}.overview_long.tsv",
                                   sep="\t", index=False)
    with open(f"{args.out_prefix}.params.json", "w") as fh:
        json.dump({"filter": args.filter, "params": params,
                   "mito_contig": mito, "min_depth": args.min_depth,
                   "n_samples": n_samples,
                   "min_callable": args.min_callable}, fh, indent=2)

    print(f"[STAGE3:{args.filter}] nuclear positions={overview['total_positions']:.0f} "
          f"mutations={overview['total_mutations']} | "
          f"mito positions={overview['total_positions_mito']:.0f} "
          f"mutations={overview['total_mutations_mito']}")
    print(f"[STAGE3:{args.filter}] nuclear CDS positions="
          f"{overview['total_positions_cds']:.0f} "
          f"CDS mutations={overview['total_mutations_cds']}")
    print(f"[STAGE3:{args.filter}] reads/allele nuclear="
          f"{overview['reads_per_allele']} mito={overview['reads_per_allele_mito']} "
          f"(a large mito value means heteroplasmy or germline mito variants "
          f"are still in the set)")
    print(f"[STAGE3:{args.filter}] nuclear synonymous sites="
          f"{overview['total_positions_synonymous']:.1f} "
          f"syn-weighted depth={overview['total_read_depth_cohort_synonymous']:.1f} "
          f"syn mutations={overview['total_mutations_synonymous']}")
    print(f"[STAGE3:{args.filter}] -> {args.out_prefix}.summary.tsv")


if __name__ == "__main__":
    main()
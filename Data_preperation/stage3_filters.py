#!/usr/bin/env python3
"""
Stage 3 - cohort site table and germline filter
===============================================
Builds one site table from the Stage 2B per-sample annotated SNVs and the
Stage 2A cohort depth arrays, then applies a filter.

The filter labels every (chrom, pos, ref, alt) site-allele with:
  callability_pass       at least `min_callable` samples have depth >=
                         MIN_DEPTH at the position
  germline_disqualified  the allele looks germline/heteroplasmic, not an error
  retained               callability_pass & ~germline_disqualified

Filters
-------
  none         nothing is disqualified. Retains germline mitochondrial
               variants and heteroplasmy, which dominate the mito numerator.
  across_rate  disqualify a site-allele if its cohort-wide alt-read fraction,
               sum_ao / sum_depth, exceeds `across_rate`. sum_ao is the alt
               reads of this (ref, alt) summed over all samples; sum_depth is
               the depth at the position summed over all samples.

A position leaves the analysed region only if all its observed alleles are
disqualified (fully_disqualified_positions). Where only some are, the position
stays and disqualified_changes() lists the individual changes.

Requires: numpy, pandas, pyarrow.
"""
import os

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------
def load_mutations(stage2b_dir, samples):
    """Concatenate the listed samples' annotated SNVs into one long DataFrame."""
    files = [os.path.join(stage2b_dir, f"{s}.annotated.parquet") for s in samples]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        raise FileNotFoundError(f"{len(missing)} annotated table(s) missing, "
                                f"e.g. {missing[:5]}")
    return pd.concat((pd.read_parquet(f) for f in files), ignore_index=True)


def load_cohort_pos_features(cohort_dir, sites):
    """Attach n_callable and sum_depth to the unique (chrom,pos) of `sites`."""
    pos = sites[["chrom", "pos"]].drop_duplicates()
    n_call = np.zeros(len(pos), dtype=np.int64)
    s_dep = np.zeros(len(pos), dtype=np.int64)
    p0 = (pos["pos"].to_numpy() - 1)
    for contig in pos["chrom"].unique():
        nc_path = os.path.join(cohort_dir, f"cohort_n_callable.{contig}.npy")
        sd_path = os.path.join(cohort_dir, f"cohort_sum_depth.{contig}.npy")
        for f in (nc_path, sd_path):
            if not os.path.exists(f):
                raise FileNotFoundError(
                    f"{f} is missing; run Stage 2A reduce for every contig "
                    f"before Stage 3")
        nc = np.load(nc_path)
        sd = np.load(sd_path)
        m = (pos["chrom"].to_numpy() == contig)
        idx = p0[m]
        n_call[m] = nc[idx]
        s_dep[m] = sd[idx]
    pos = pos.assign(n_callable=n_call, sum_depth=s_dep)
    return pos


def build_site_table(mut_df, cohort_dir):
    """One row per (chrom,pos,ref,alt) with cohort aggregates."""
    g = mut_df.groupby(["chrom", "pos", "ref", "alt"], sort=False)
    site = g.agg(
        n_samples=("sample_id", "size"),
        sum_ao=("ao", "sum"),
        # constant within a (chrom, pos, ref, alt) group
        spectrum_channel=("spectrum_channel", "first"),
    ).reset_index()
    posf = load_cohort_pos_features(cohort_dir, site)
    site = site.merge(posf, on=["chrom", "pos"], how="left")
    site["sum_depth"] = site["sum_depth"].fillna(0).astype(np.int64)
    site["n_callable"] = site["n_callable"].fillna(0).astype(np.int64)
    return site


# ----------------------------------------------------------------------------
# filters: each returns germline_disqualified, a bool Series aligned to `site`
# ----------------------------------------------------------------------------
def filter_none(site, params):
    """No disqualification (unfiltered baseline)."""
    return pd.Series(False, index=site.index)


def filter_across_rate(site, params):
    """
    Disqualify a site-allele if sum_ao / sum_depth > params["across_rate"].

    Allele-specific: at a position where A>G exceeds the threshold, A>C is
    judged on its own pooled rate. Alleles at positions with zero summed depth
    are disqualified.
    """
    rate = site["sum_ao"] / site["sum_depth"].where(site["sum_depth"] > 0, np.nan)
    return (rate > params["across_rate"]) | rate.isna()


FILTERS = {
    "none": filter_none,
    "across_rate": filter_across_rate,
}


def apply_filter(site, name, params, min_callable):
    dq = FILTERS[name](site, params)
    out = site.copy()
    out["callability_pass"] = (site["n_callable"] >= min_callable).to_numpy()
    out["germline_disqualified"] = dq.to_numpy()
    out["retained"] = out["callability_pass"] & ~out["germline_disqualified"]
    return out


def fully_disqualified_positions(site_labeled):
    """Positions with >=1 observed alt, all disqualified, none retained."""
    g = site_labeled.groupby(["chrom", "pos"])
    summ = g.agg(n_alts=("alt", "size"),
                 n_retained=("retained", "sum"),
                 n_dq=("germline_disqualified", "sum")).reset_index()
    bad = summ[(summ["n_alts"] >= 1) & (summ["n_retained"] == 0) &
               (summ["n_dq"] == summ["n_alts"])]
    return set(zip(bad["chrom"], bad["pos"]))


def disqualified_changes(site_labeled):
    """
    DataFrame[chrom, pos, channel] of the changes the filter disqualified.

    `channel` is the transcript-orientation ref>alt from Stage 2B, the key the
    Stage 0C channel table uses. Disqualified alleles with no channel are
    dropped; they lie outside every transcript.
    """
    dq = site_labeled.loc[site_labeled["germline_disqualified"]
                          & site_labeled["spectrum_channel"].notna(),
                          ["chrom", "pos", "spectrum_channel"]]
    return (dq.rename(columns={"spectrum_channel": "channel"})
              .drop_duplicates()
              .reset_index(drop=True))


def retained_mutations(mut_df, site_labeled):
    """Subset the long mutation table to retained site-alleles."""
    keep = site_labeled.loc[site_labeled["retained"], ["chrom", "pos", "ref", "alt"]]
    return mut_df.merge(keep, on=["chrom", "pos", "ref", "alt"], how="inner")
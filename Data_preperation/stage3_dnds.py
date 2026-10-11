#!/usr/bin/env python3
"""
Stage 3 (dN/dS) - per-sample synonymous / nonsynonymous sites
=============================================================
Nei-Gojobori potential-site counts for each sample, summed over the CDS
positions that are both in the filter's analysed region and callable in that
sample:

  N, S        potential missense / synonymous sites, summed from the
              per-change rows of Stage 0C, excluding changes the filter
              disqualified. Changes that create a stop codon have no row and
              are in neither, so N pairs with Dn counted as missense only.
  N_c, S_c    the same, split by the 12 mutation channels
  opp_c       per channel, the number of positions carrying that channel's
              source base, less the positions where the filter disqualified
              that change

Nuclear and mitochondrial contigs are accumulated separately: they use
different genetic codes and have very different base composition, so their
potential-site counts are not interchangeable.

dnds_from_components() turns these into:
  pN = Dn/N, pS = Ds/S, dN/dS = pN/pS
  N_w, S_w, dNdS_weighted   the same sites weighted by the sample's own
                            spectrum, as a rate per opportunity (count_c /
                            opp_c), so that a uniform spectrum reproduces the
                            unweighted ratio

Dn and Ds are read counts (summed AO) while N and S are site counts; the ratio
pN/pS is unaffected as long as reads per allele do not differ systematically
between missense and synonymous alleles.

Method: per contig, scatter the per-change weights onto genome-length arrays
masked to the analysed region, take a cumulative sum, then sum over each
sample's callable intervals. One contig is held in memory at a time.

Requires: numpy.
"""
import numpy as np

CHANNELS = ["A>C", "A>G", "A>U", "C>A", "C>G", "C>U",
            "G>A", "G>C", "G>U", "U>A", "U>C", "U>G"]
CH_IDX = {c: i for i, c in enumerate(CHANNELS)}

# source base of each channel, in transcript orientation
BASES = ["A", "C", "G", "U"]
BASE_IDX = {b: i for i, b in enumerate(BASES)}
CH_SRC = np.array([BASE_IDX[c[0]] for c in CHANNELS])


def _change_keys(pos0, cidx):
    """One integer per (position, channel), for set membership tests."""
    return pos0.astype(np.int64) * 16 + cidx


def live_mask(pos0, cidx, dq_contig):
    """
    Boolean mask over change rows: True where the change was NOT disqualified.

    dq_contig: DataFrame[pos, channel] for this contig, or None.
    """
    if dq_contig is None or len(dq_contig) == 0:
        return np.ones(pos0.size, dtype=bool)
    dq_cidx = np.array([CH_IDX.get(c, -1) for c in dq_contig["channel"]],
                       dtype=np.int64)
    dq_keys = _change_keys(dq_contig["pos"].to_numpy() - 1, dq_cidx)
    return ~np.isin(_change_keys(pos0, cidx), dq_keys)


def _cumsum0(arr):
    """cumulative sum with a leading 0 so interval [a,b) = cum[b]-cum[a]."""
    out = np.zeros(arr.size + 1, dtype=np.float64)
    np.cumsum(arr, out=out[1:])
    return out


def compute_sites_per_sample(ns_channel, sizes, analysed_mask,
                             callable_by_sample, contig_groups=None,
                             disqualified=None):
    """
    ns_channel: DataFrame[chrom,pos,channel,is_syn,weight] from Stage 0C
                (pos 1-based)
    sizes: {contig: length}
    analysed_mask: {contig: bool array} (the filter's analysed region)
    callable_by_sample: {sample: {contig: [(start0,end), ...]}}
    contig_groups: {contig: group_label}; contigs absent from the mapping fall
        into 'nuclear'. Sites are accumulated per group, so the mitochondrion
        gets its own N/S.
    disqualified: DataFrame[chrom,pos,channel] from
        stage3_filters.disqualified_changes, or None for no exclusion.

    -> {sample: {group: {"N","S","N_c"(12,),"S_c"(12,),"opp_c"(12,)}}}
    """
    samples = list(callable_by_sample)
    if contig_groups is None:
        contig_groups = {c: "nuclear" for c in sizes}
    groups = sorted(set(contig_groups.get(c, "nuclear") for c in sizes))

    def _blank():
        return {"N": 0.0, "S": 0.0, "N_c": np.zeros(12), "S_c": np.zeros(12),
                "_opp_base": np.zeros(4), "_dq_c": np.zeros(12)}

    acc = {s: {g: _blank() for g in groups} for s in samples}

    ch_by_c = {c: g for c, g in ns_channel.groupby("chrom")}
    dq_by_c = ({c: g for c, g in disqualified.groupby("chrom")}
               if disqualified is not None and len(disqualified) else {})

    for contig, n in sizes.items():
        if contig not in ch_by_c:
            continue
        grp = contig_groups.get(contig, "nuclear")
        mask = analysed_mask.get(contig)
        if mask is None:
            mask = np.ones(n, dtype=bool)

        gc = ch_by_c[contig]
        p0 = gc["pos"].to_numpy() - 1
        cidx = np.array([CH_IDX.get(c, -1) for c in gc["channel"]], dtype=np.int64)
        issyn = gc["is_syn"].to_numpy().astype(bool)
        w = gc["weight"].to_numpy().astype(np.float64)

        in_region = mask[p0] & (cidx >= 0)
        keep = in_region & live_mask(p0, cidx, dq_by_c.get(contig))

        def _scatter(sel):
            na = np.zeros(n, dtype=np.float64)
            sa = np.zeros(n, dtype=np.float64)
            np.add.at(sa, p0[sel & issyn], w[sel & issyn])
            np.add.at(na, p0[sel & ~issyn], w[sel & ~issyn])
            return _cumsum0(na), _cumsum0(sa)

        cum_n, cum_s = _scatter(keep)

        cum_nc = [None] * 12
        cum_sc = [None] * 12
        for k in range(12):
            cum_nc[k], cum_sc[k] = _scatter(keep & (cidx == k))

        def _count(sel):
            flag = np.zeros(n, dtype=np.float64)
            flag[p0[sel]] = 1.0
            return _cumsum0(flag)

        # positions per source base, and per channel the positions where the
        # change was disqualified
        cum_opp = [_count(in_region & (CH_SRC[cidx] == bi)) for bi in range(4)]
        dq_sel = in_region & ~keep
        cum_dq = ([_count(dq_sel & (cidx == k)) for k in range(12)]
                  if dq_sel.any() else None)

        for s in samples:
            a_g = acc[s][grp]
            for (a, b) in callable_by_sample[s].get(contig, []):
                a_g["N"] += cum_n[b] - cum_n[a]
                a_g["S"] += cum_s[b] - cum_s[a]
                for k in range(12):
                    a_g["N_c"][k] += cum_nc[k][b] - cum_nc[k][a]
                    a_g["S_c"][k] += cum_sc[k][b] - cum_sc[k][a]
                    if cum_dq is not None:
                        a_g["_dq_c"][k] += cum_dq[k][b] - cum_dq[k][a]
                for bi in range(4):
                    a_g["_opp_base"][bi] += cum_opp[bi][b] - cum_opp[bi][a]

    for per_group in acc.values():
        for a_g in per_group.values():
            a_g["opp_c"] = a_g.pop("_opp_base")[CH_SRC] - a_g.pop("_dq_c")
    return acc


def dnds_from_components(Dn, Ds, N, S, N_c, S_c, spectrum_counts, opportunities):
    """
    spectrum_counts: array(12,) of the sample's observed channel counts.
    opportunities:   array(12,) of positions carrying each channel's source base.

    The spectrum weights are per-opportunity rates, count_c / opp_c, not raw
    channel counts: N_c already counts only positions whose base can undergo
    that change, so weighting by raw counts would apply base composition twice
    and a uniform mutation spectrum would not reproduce the unweighted ratio.
    Channels with no opportunities get zero weight.

    Returns dict with unweighted and spectrum-weighted dN/dS (None if undefined).
    """
    def ratio(a, b):
        return (a / b) if b and b > 0 else None

    pN = ratio(Dn, N)
    pS = ratio(Ds, S)
    dnds = ratio(pN, pS) if (pN is not None and pS) else None

    opp = np.asarray(opportunities, dtype=float)
    rate = np.divide(np.asarray(spectrum_counts, dtype=float), opp,
                     out=np.zeros(12), where=opp > 0)
    tot = rate.sum()
    if tot > 0:
        w = rate / tot                            # normalisation cancels in ratio
        N_w = float(np.dot(w, N_c))
        S_w = float(np.dot(w, S_c))
    else:
        N_w = S_w = 0.0
    pN_w = ratio(Dn, N_w)
    pS_w = ratio(Ds, S_w)
    dnds_w = ratio(pN_w, pS_w) if (pN_w is not None and pS_w) else None

    return {"pN": pN, "pS": pS, "dNdS": dnds,
            "N_w": N_w, "S_w": S_w, "dNdS_weighted": dnds_w}
This is the code repository for the detection of transcription errors.

An overview of the transcription error detection workflow is available in the figure workflow_overview.pdf.

The pipeline consist of the following main steps: 
-------------------------------------------------
1) Assembly of forward and reverse read with Pear:
pear -f P20702_139_R1.fastq -r P20702_139_R2.fastq -o P20702_139.assembled.fastq -y 2000M -c 80

2) Trimming the overhangs: 
cutadapt -q 60,60 -m 1 -o P20702-139.assembled.trimmed60.fastq P20702-139.assembled.fastq

3) FastQC of merged and trimemd reads

4) Mapping the merged reads to the reference genome: 
STAR --genomeDir starIndex \
--readFilesIn P20702-139.assembled.trimmed60.fastq \
--runThreadN 1 \
--outFileNamePrefix bam_tempdir/P20702-139 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard

5) Detecting variants:
freebayes -f ${ref} \
-r ${region} \
--haplotype-length 0 \
--min-alternate-count 1 \
--min-alternate-fraction 0 \
--pooled-continuous \
--min-coverage 100 \
--use-best-n-alleles 4 \
-i \
${base_dir}/bam/${sample}.Aligned.sortedByCoord.out.bam > ${base_dir}/outputs/${sample}_${region}.freebayes.vcf

6) Define the number of positions with "detectable variants" (sequence space):
samtools depth -Q 1 -q 73 -b $genes $bamfile > ${depth_dir}'/'${sampleID}.genes.depth
Followed by processing in python.



# Transcription-error calling pipeline

Variant calling and rate estimation for RNA transcription errors in
*Drosophila melanogaster*, from freebayes VCFs to per-sample error rates,
mutation spectra and dN/dS.

`submit_pipeline.sh` submits every stage to SLURM with the right dependencies:

```
./submit_pipeline.sh samples.txt              # run everything
./submit_pipeline.sh samples.txt --from 2B    # resume from one stage
```

`samples.txt` has one sample ID per line. Stages are `0A 0B 0C 1 2A 2AR 2B 3`.

```
0A ref --+--> 0B snpEff DB ---------------\
         +--> 0C N/S table ---------------+--> 2B annotate --\
         +--> 1 SNVs (array) -------------/                   +--> 3 (array)
         +--> 2A depth (array, also after 0C) --> 2A reduce --/
```

---


### Inputs the pipeline expects

| Path | What |
|---|---|
| `refgenome/genome.fa` (+ `.fai`) | the FASTA freebayes was run on |
| `refgenome/genes.gtf` | annotation matching that FASTA (the one the STAR index was built from) |
| `transcription_errors/concatenated_freebayes_clean/<sample>.concatenated.vcf` | per-sample freebayes calls, one `#CHROM` header per file |
| `transcription_errors/bam/<sample>.Aligned.sortedByCoord.out.bam` (+ index) | the STAR alignments those calls came from |

### Thresholds

| Setting | Value | Where |
|---|---|---|
| Mean alt base quality (QA/AO) to keep an allele | ≥ 74 | `QA` in `JOB_stage1.sh` |
| Base quality counted as depth | ≥ 73 | `BASEQ` in `JOB_stage2A.sh` |
| Mapping quality counted as depth | ≥ 1 | `MAPQ` in `JOB_stage2A.sh` |
| Depth for a position to be callable in a sample | ≥ 100 | `MIN_DEPTH` in `JOB_stage2A.sh`, `JOB_stage2A_reduce.sh` and `JOB_stage3.sh` — **keep all three equal** |
| Samples callable at a position for it to be analysed | ≥ 1 | `MIN_CALLABLE` in `JOB_stage3.sh` |
| Filter settings | 5 | `SETTINGS` in `JOB_stage3.sh`; `submit_pipeline.sh` counts them itself |

The allele threshold is a *mean* over the reads supporting an allele, while
depth is counted *per base*. freebayes was run without `--min-base-quality`,
so individual low-quality reads can contribute to an allele whose mean passes.

---

## Stage 0A — canonical reference

`JOB_stage0A.sh`, with `snpeff_gtf.py` and `verify_ref_gtf.py`

Restricts the FASTA and GTF to 2L, 2R, 3L, 3R, 4, X, Y and the mitochondrion
(renamed `M`), normalises the GTF so snpEff parses coding transcripts, checks
that the GTF coordinates describe the FASTA, and writes the region BEDs.

**Reads:** `refgenome/genome.fa` (+ `.fai`), `refgenome/genes.gtf`

**Writes,** all into `refgenome/canonical/`:

| File | Used by |
|---|---|
| `genome.canonical.fa` (+ `.fai`) | 0B, 0C, 1, 2A, 2A reduce, 2B, 3 |
| `annotation.canonical.gtf` | 0B, 0C, 2B |
| `genes.bed` — merged gene spans | 3 (the `genespan` region) |
| `analysis.bed` — the whole canonical genome | 2A (the region depth is measured over) |
| `reference_check.tsv` — per-contig CDS translation report | nothing; for inspection |

Stops if a chromosome is missing, if the GTF has no CDS features, if the GTF
still has rows snpEff would misparse, or if the CDS translations disagree with
the sequence (below 90% free of internal stops, or below 80% with a start or
stop codon).

## Stage 0B — snpEff database

`JOB_stage0B.sh`, with `snpeff_gtf.py` (`make-testvcf`, `score`)

Builds a snpEff database from the canonical reference, so annotation uses the
same coordinates, contig names and transcripts as everything else, with the
invertebrate mitochondrial codon table on `M`.

**Reads:** `refgenome/canonical/genome.canonical.fa`, `annotation.canonical.gtf`

**Writes:**

| File | Used by |
|---|---|
| `refgenome/snpeff_data/BDGP6.46.canonical/` — the database | 2B |
| `refgenome/snpeff_custom.config` | 2B |
| `build.log`, the probe VCF and its annotated copy | nothing; for inspection |

Two gates stop the job if the database is wrong: at least 20% of transcripts
must be protein-coding in the build log, and at least 50% of SNVs placed at
real CDS positions must come back with a coding consequence.

## Stage 0C — N/S site table

`JOB_stage0C.sh`, with `stage0c_build_ns_sites.py`

Nei-Gojobori potential-site counts for every CDS position of the longest-CDS
transcript of each gene. Changes that create a stop codon are excluded from
both N and S, so N pairs with a missense-only Dn.

**Reads:** `refgenome/canonical/genome.canonical.fa`, `annotation.canonical.gtf`

**Writes,** into `refgenome/canonical/`:

| File | Used by |
|---|---|
| `ns_sites.parquet` — `chrom, pos, n_site, s_site, gene_id, transcript_id` | 2B (which transcripts are representative), 3 (N and S, and the synonymous depth weights) |
| `ns_sites.channels.parquet` — `chrom, pos, channel, is_syn, weight` | 3 (per-channel N_c / S_c for the weighted dN/dS) |
| `ns_sites.cds.bed` — those positions as merged intervals | 2A and 3 (the `cds` region) |

## Stage 1 — per-sample SNV extraction

`JOB_stage1.sh` (array, one task per sample), with `stage1_extract_snvs.py`
and `contig_names.py`

Filters one sample's freebayes VCF to primitive SNVs: canonicalises contig
names, trims padded alleles, keeps single-base changes whose mean alt base
quality is at least 74, and checks each REF base against the reference.

**Reads:** `concatenated_freebayes_clean/<sample>.concatenated.vcf`,
`refgenome/canonical/genome.canonical.fa`

**Writes,** into `transcription_errors/stage1/`:

| File | Used by |
|---|---|
| `<sample>.snv.parquet` | 2B |
| `<sample>.stage1.json` — per-contig, trimming and merge counts | nothing; for inspection |

## Stage 2A — per-sample depth

`JOB_stage2A.sh` (array, one task per sample), with `stage2a_depth_sample.py`
and `contig_names.py`

Runs `samtools depth` over the analysis region and produces the denominators.

**Reads:** `transcription_errors/bam/<sample>.Aligned.sortedByCoord.out.bam`,
`refgenome/canonical/analysis.bed`, `ns_sites.cds.bed`, `genome.canonical.fa.fai`

**Writes,** into `transcription_errors/stage2a/per_sample/`:

| File | Used by |
|---|---|
| `<sample>.depth.npz` — per-base depth | 2A reduce, 3 |
| `<sample>.callable.bed` — runs at depth ≥ 100 | 3 |
| `<sample>.aggregates.json` — depth sums, total / CDS / non-CDS | nothing; for inspection |

## Stage 2A reduce — cohort depth

`JOB_stage2A_reduce.sh`, with `stage2a_reduce.py`

Runs once after the whole 2A array. Sums the per-sample depth into per-base
cohort arrays. Stops if any listed sample has no `depth.npz`.

**Reads:** every `<sample>.depth.npz`, `genome.canonical.fa.fai`

**Writes,** into `transcription_errors/stage2a/cohort/`:

| File | Used by |
|---|---|
| `cohort_sum_depth.<contig>.npy` — depth summed over all samples | 3 (the filter, and the overview denominators) |
| `cohort_n_callable.<contig>.npy` — samples callable at each base | 3 (the filter, and the analysed region) |

## Stage 2B — annotation

`JOB_stage2B.sh`, with `stage2b_build_union_vcf.py`, snpEff,
`stage2b_finalize.py` and `contig_names.py`

Collects the distinct SNVs of all samples into one VCF, annotates it once, and
adds the consequence, transcribed strand and spectrum channel back onto each
sample's table. A consequence is taken only from a representative transcript;
anything else is labelled `other`, so synonymous and missense counts cover the
same positions as N and S.

**Reads:** every `stage1/<sample>.snv.parquet`, the snpEff database and config
from 0B, `ns_sites.parquet` from 0C, `genome.canonical.fa.fai`

**Writes,** into `transcription_errors/stage2b/`:

| File | Used by |
|---|---|
| `per_sample/<sample>.annotated.parquet` | 3 |
| `work/union.vcf` and its annotated copy | nothing; for inspection |

Stops if almost nothing matches a representative transcript, if no strand
could be derived, or if too few nuclear SNVs get a coding consequence.

## Stage 3 — filtering and summary

`JOB_stage3.sh` (array, one task per filter setting), with `stage3_filters.py`,
`stage3_dnds.py`, `stage3_summary.py` and `contig_names.py`

Applies one filter setting, then counts mutations and depth over the analysed
region and writes the results. The filter disqualifies a site-allele whose
cohort-wide alt-read fraction (summed AO over summed depth) exceeds a
threshold; a position leaves the analysed region only if all of its alleles are
disqualified.

| Label | Filter | Threshold |
|---|---|---|
| `unfiltered` | none | — |
| `across_5pct` | across_rate | 5% |
| `across_1pct` | across_rate | 1% |
| `across_0p1pct` | across_rate | 0.1% |
| `across_0p01pct` | across_rate | 0.01% |

**Reads:** every `stage2b/per_sample/<sample>.annotated.parquet`, the cohort
arrays from 2A reduce, each sample's `depth.npz` and `callable.bed` from 2A,
`ns_sites.parquet` and `ns_sites.channels.parquet` from 0C,
`genome.canonical.fa` (+ `.fai`), `genes.bed`, `ns_sites.cds.bed`

**Writes,** into `transcription_errors/stage3/`, four files per label:

| File | Contents |
|---|---|
| `<label>.summary.tsv` | one row per sample: depth sums, positions covered, mutation counts, rates, N, S, pN, pS, dN/dS and the 12 spectrum channels, each with a `_mito` twin |
| `<label>.overview.tsv` | one row for the cohort, per compartment × region |
| `<label>.overview_long.tsv` | the same, one row per compartment × region |
| `<label>.params.json` | the exact settings used |

These are the inputs to the downstream R analysis. Every column is defined in
the docstring at the top of `stage3_summary.py`.

---

## Counting conventions

- **Mutations are reads, not positions.** Every `*_mutations` column is a sum
  of AO, the alt-supporting read count, so it pairs with a depth denominator.
  Distinct site-alleles are reported separately as `*_alleles`, and their ratio
  as `reads_per_allele`.
- **Numerator and denominator cover the same bases.** A mutation counts only
  where its own sample contributes depth: inside the analysed region and inside
  that sample's callable intervals. What this removes is reported as
  `excluded_mutations` / `excluded_alleles`.
- **Nuclear and mitochondrial results are separate throughout.** Unsuffixed
  columns are nuclear; `_mito` columns are the mitochondrion, which uses its
  own genetic code.
- **Synonymous rates use a weighted denominator.** `depth_sum_synonymous` is
  the depth at each CDS position weighted by `s_site`, the fraction of that
  base's possible changes that are synonymous.
- **Two spectra.** `spec_*` covers gene spans and pairs with
  `depth_sum_genespan`; `spec_*_cds` is restricted to coding positions and
  pairs with `depth_sum_cds`. Channels are written in transcript orientation
  with T as U.
- **`positions_with_A/C/G/U`** in the overview count reference bases on the
  **reference strand**, not in transcript orientation.

## Files

| File | Role |
|---|---|
| `submit_pipeline.sh` | submits every stage with its dependencies |
| `JOB_stage*.sh` | one SLURM job per stage |
| `snpeff_gtf.py` | GTF normalisation for snpEff and database checks (0A, 0B) |
| `verify_ref_gtf.py` | checks GTF coordinates against the FASTA (0A) |
| `stage0c_build_ns_sites.py` | N/S site table (0C) |
| `stage1_extract_snvs.py` | VCF to filtered SNV table (1) |
| `stage2a_depth_sample.py` | per-sample depth (2A) |
| `stage2a_reduce.py` | cohort depth (2A reduce) |
| `stage2b_build_union_vcf.py`, `stage2b_finalize.py` | annotation (2B) |
| `stage3_filters.py` | site table and germline filter (3) |
| `stage3_dnds.py` | per-sample N/S sites and dN/dS (3) |
| `stage3_summary.py` | counting, summary and overview (3) |
| `contig_names.py` | shared contig vocabulary (1, 2A, 2B, 3) |

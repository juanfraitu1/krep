# Handoff: krep vs RepeatMasker on T2T-CHM13

## The two things that were wrong

**1. K-mer counting was scoped per FASTA record, not per genome.**
`mask_fasta_graph` / `mask_sequence` rebuilt a CBF from whichever sequence was
passed in, so "how often does this k-mer occur" was answered against a 10 Mb
slice. A family with a few hundred copies genome-wide has an expected count
below 1 in that window and was structurally invisible. The 10 Mb slice was also
unrepresentative: 32.1% soft-masked vs 41.0% for chr1 and 40.3% genome-wide
(it is the p-arm subtelomere and misses the centromeric satellite arrays).

**2. The ground truth was not RepeatMasker.** The lowercase in NCBI's
`*_genomic.fna` covers 40.27% of the genome; the real RepeatMasker annotation
covers 54.19%. On chr1 only 83% of the lowercase is explained by RepeatMasker,
and the lowercase captures just 62% of RepeatMasker's calls. The NCBI mask is a
WindowMasker/TRF-style de novo mask, so the old F1 of 0.7737 measured agreement
with *another de novo masker*.

## What was added

- **`krep index`** — one streaming pass over the genome producing exact
  genome-wide k-mer counts. Hash sub-sampling (FracMinHash, `--sample`) keeps
  1 in N k-mers; because the decision depends only on the k-mer's own hash, a
  tracked k-mer is counted at every occurrence, so counts stay *exact* (unlike a
  Bloom filter, which inflates via collisions). Buffers spill to disk as sorted
  runs and are merged in one pass, so peak RAM is `--buffer`, not genome size.
- **`krep mask --index`** — streams the target and looks up genome-wide counts.
  Mask chr1 alone while using whole-genome context.
- **Streaming `compare-mask`** with per-record breakdown (the old one loaded
  both FASTAs, ~6 GB at genome scale).
- **`fasta::FastaStream` / `FastaWriter`** — record-at-a-time I/O.

## Bugs fixed along the way

- **BED chrom field contained the whole FASTA description line**
  (`NC_060925.1 Homo sapiens isolate CHM13 chromosome 1, ...`), which is invalid
  BED and cannot be joined against `rm.out`. Now uses the sequence ID.
- **`evaluate` ignored the chromosome column** — it swept all coordinates on one
  axis, so genome-wide BEDs would conflate chr1 and chr2 positions. Now grouped
  per chromosome.
- **Histogram bucketing used `63 - leading_zeros()` on a `u32`** (should be 31),
  putting every count in one bogus bucket.
- **`canonical()` recomputed the reverse complement with a k-iteration loop per
  k-mer** (5.4e10 ops genome-wide). `KmerIter` now rolls both strands
  incrementally; a test asserts equivalence with the old implementation.

## Current results

Genome index: `k=18 --sample 16 --min-count 2` over 3,117,275,501 bp →
**32 s**, 19.6M entries, 235 MB. K-mer count validates exactly
(3,117,275,093 = bases − 17 × 24 records).

Full chr1 masked in **1.8 s**. Against RepeatMasker (`chr1_rm.bed`):

| t | gap | P | R | F1 |
|---|-----|---|---|-----|
| 2 | 150 | 0.610 | 0.911 | 0.731 |
| 3 | 250 | 0.678 | 0.850 | 0.754 |
| **4** | **250** | **0.767** | **0.765** | **0.766** |
| 10 | 150 | 0.939 | 0.568 | 0.708 |

Against the NCBI lowercase on the old 10 Mb slice, genome-wide counting lifts
recall from 0.689 to 0.760 at equal F1 — without `--mismatch1` and despite 16x
sub-sampling.

## The ceiling

Per-family recall splits by family age. Nearly complete: Satellite (0.999),
SVA (0.993), Alu (0.978), ERVK (0.969), L1 (0.826). Largely missed: MIR (0.378),
L2 (0.294), CR1 (0.202), most DNA transposons (0.12-0.43).

**54% of missed bases on chr1 are ancient diverged families.** A 30%-diverged
MIR copy shares essentially no exact 18-mers with any other copy, so no
abundance method recovers it. Closing that gap needs consensus building plus
alignment, not tuning.

## Suggested next steps

1. **Re-tune on a representative target.** Use full chr1 (or several
   chromosomes) rather than the 10 Mb subtelomere slice, and score against
   `chr1_rm.bed`, not the fna lowercase.
2. **Spaced seeds** instead of contiguous k-mers — the single highest-value
   change for the ancient families. A weight-14 seed over an 18-base span
   tolerates mismatches in don't-care positions, which is exactly what MIR/L2
   need. Sub-sampling is compatible (hash the extracted seed).
3. **Multiple k in one index** (`k=15` alongside `k=18`) and take the union;
   smaller k is more divergence-tolerant but needs a higher threshold.
4. **`--mismatch1` is incompatible with hash sub-sampling** — a Hamming
   neighbour of a sampled k-mer is almost never itself sampled. Either run
   `--sample 1` for mismatch expansion, or use spaced seeds instead.
5. **Boundary refinement.** Seeds are ~16 bp apart, so region edges are fuzzy to
   that scale; trimming edges against a denser local scan would lift precision.
6. Genome-wide masking run (all 24 records) to confirm scaling end to end.

## Data locations

- Genome index: `C:\krep_work\chm13.k18.s16.kidx`
- RepeatMasker BED: `C:\krep_work\chr1_rm.bed` (merged),
  `C:\krep_work\chr1_rm_family.bed` (4-column, for per-family recall)
- Raw annotation: `~/krep_data/rm.out.gz` (from NCBI FTP, 198 MB)
- Scratch/temp: `C:\krep_work\tmp` (kept off OneDrive deliberately)

## Build

WSL cannot host the build (no Rust there, and only 3.9 GB RAM); the Windows
toolchain works from WSL:

```bash
export PATH="/mnt/c/mingw64/bin:$PATH"
/mnt/c/Users/jfris/.cargo/bin/cargo.exe build --release --target-dir target3
```

Pass Windows-style paths (`C:/...`) to the resulting `.exe`.

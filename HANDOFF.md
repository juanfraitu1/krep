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

## Round 2: spaced seeds, hysteresis, density gate (measured)

Goal was recall on the ancient families. Added, all behind flags with
defaults that preserve prior behaviour:

- `krep index --seed PATTERN` — symmetric spaced seeds (`kmer::SpacedSeed`,
  `SeedIter`); pattern stored in the index header (format v2; v1 indices are
  rejected with a rebuild message).
- `krep index --passes N` — hash-partitioned multi-pass counting so `--sample 1`
  fits in bounded temp/RAM. Test asserts byte-identical output for 1 vs 4 passes.
- `krep mask --index A --index B --index-threshold tA,tB` — union of several
  indices with per-index thresholds.
- `--index-threshold-low` (hysteresis) and `--min-hits` / `--min-density`
  (component density gate).

Results on full chr1 vs RepeatMasker (baseline k18@4 gap 250: F1 0.766):

- Hysteresis alone: **no gain** (best 0.762). Lowering the extension threshold
  raises recall but precision collapses because weak seeds *bridge* unique DNA
  between nearby repeats. Bridging, not nucleation, is the FP source.
- Spaced w16/s22 at `--sample 1`: at low thresholds it masks the whole
  chromosome (composition background: 47M distinct seeds with count >= 8). At
  t=64 with gap 100-150 it matches k18's F1 with higher P / lower R. Weight
  cannot drop below 16 at 3.1 Gb, so the sensitivity gain is only ~2x per copy
  and the required threshold increase eats it for mid-count families.
- **Union k18@4 ∪ spaced@64, gap 150, density 0.03: P 0.749 / R 0.802 /
  F1 0.7745** — best found. MIR 0.378→0.468, hAT-Charlie 0.507→0.593,
  L1 0.826→0.888, Simple_repeat 0.832→0.938. L2 flat (0.294→0.309).

Everything converges at F1 ≈ 0.775. That is the de novo copy-vs-copy ceiling
here; the remaining ~20% of RepeatMasker's bases are copies that share no seed
with any other copy at a usable weight.

## What would actually move the ceiling

1. **Consensus + alignment (RepeatScout-style).** Seed from high-count entries
   in the index, gather occurrences with flanks, extend a consensus by majority
   vote, then align consensi back (seed-chain-extend, banded) and mask hits.
   Helps L1M / MaLR / Charlie most; MIR/L2 need copies that align to each
   other, which is exactly what they lack.
2. **Library mode (Dfam consensi).** Only route to MIR/L2 parity. Stops being
   de novo. Index the ~1,500 human consensi with spaced seeds, scan the genome,
   extend hits.
3. **Density-gated bridging** (not implemented): allow a gap to be bridged only
   if hits inside it meet a local density, decoupling "extend through mutated
   stretch" from "jump across unique DNA". The density gate here is per
   component, not per gap.
4. Genome-wide run across all 24 records to confirm the streaming path end to
   end at scale.

## Data locations

- Genome indices: `C:\krep_work\chm13.k18.s16.kidx` (contiguous 18-mer,
  1/16 sampled, 235 MB) and `C:\krep_work\chm13.sp16.s1.kidx` (spaced
  w16/s22, dense, 560 MB)
- RepeatMasker BED: `C:\krep_work\chr1_rm.bed` (merged),
  `C:\krep_work\chr1_rm_family.bed` (4-column, for per-family recall)
- Raw annotation: `~/krep_data/rm.out.gz` (from NCBI FTP, 198 MB)
- Scratch/temp: `C:\krep_work\tmp` (kept off OneDrive deliberately)

## Build and test

Windows Application Control intermittently blocks freshly built executables
(os error 4551) — it hit the test harness and Cargo build scripts, while the
release `krep.exe` kept running. So:

- **Tests: WSL.** Rust (rustup, minimal) and `build-essential` are installed
  in WSL. Build on ext4, not `/mnt/c`:
  `~/.cargo/bin/cargo test --release --target-dir ~/krep_target`
- **Runs: Windows exe**, because it sees the host's 8 GB rather than WSL's
  3.9 GB. Build with the Windows toolchain from WSL and copy off OneDrive:

```bash
export PATH="/mnt/c/mingw64/bin:$PATH"
/mnt/c/Users/jfris/.cargo/bin/cargo.exe build --release --target-dir target3
cp target3/release/krep.exe /mnt/c/krep_work/krep.exe
```

Pass Windows-style paths (`C:/...`) to the `.exe`.

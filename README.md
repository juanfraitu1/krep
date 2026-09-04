# krep

A lightweight, de novo repeat masker that runs comfortably on a laptop. It uses a **counting Bloom filter** to find over-represented k-mers (WindowMasker-style), with optional **graph-based connected-component masking**, **Hamming-distance-1 neighborhood expansion**, and **de Bruijn graph assembly** to merge repeat signals. No RepeatMasker, no RepBase, no Dfam, no BLAST — just one small Rust binary.

## What it does

1. **`mock`** — generates a synthetic genome with known repeat families (including segmental duplications) and a ground-truth BED file.
2. **`mask`** — counts canonical k-mers in a counting Bloom filter, then masks regions dense in over-represented k-mers (window density) or merges linked high-count k-mers into coherent repeat blocks (graph / assembly modes). Outputs BED or GTF and optional soft/hard-masked FASTA.
3. **`demask`** — converts a soft-masked FASTA to an all-uppercase FASTA.
4. **`compare-mask`** — compares two soft-masked FASTAs base-by-base (e.g. RepeatMasker vs krep) and reports precision/recall/F1.
5. **`evaluate`** — compares the predicted BED to the ground-truth BED and reports precision, recall, F1, and per-family recall.

## Build

Rust (GNU toolchain) + MinGW-w64 linker are required on Windows. The project was built and tested with rustup + WinLibs MinGW-w64.

A convenient WinLibs package can be downloaded from the WinLibs GitHub releases:
<https://github.com/brechtsanders/winlibs_mingw/releases/tag/14.2.0posix-12.0.0-ucrt-r3>

```bash
cargo build --release
```

The resulting binary is `target/release/krep` (or `krep.exe` on Windows). If a local security policy blocks the binary in the default `target/` directory, build to an alternate directory:

```bash
cargo build --release --target-dir target2
```

## Quick start

```bash
# Generate a 10 Mb mock genome with ~5% divergence repeat copies
krep mock --size 10m --seed 42 --out genome.fa --bed repeats.bed

# Mask repetitive regions (fast default mode)
krep mask --genome genome.fa --out masked.bed

# High-accuracy mode on the mock genome
krep mask --genome genome.fa --graph --graph-gap 100 --k 18 --threshold 5 --mismatch1 --out masked.bed

# Evaluate
krep evaluate --truth repeats.bed --pred masked.bed

# De-mask a soft-masked FASTA (e.g. from RepeatMasker / NCBI)
krep demask --genome GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --out chm13_unmasked.fa

# Mask a real genome and output GTF as well as soft-masked FASTA
krep mask --genome chm13_unmasked.fa --graph --graph-gap 100 --k 18 --threshold 5 --mismatch1 --out-format gtf --out chm13_masked.gtf --soft chm13_krep_soft.fa

# Compare krep's soft mask to another soft-masked FASTA
krep compare-mask --reference GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --predicted chm13_krep_soft.fa
```

### Genome-wide counting (recommended for real genomes)

The commands above count k-mers *within whatever FASTA you pass in*. That is the
wrong scope for a real genome: a repeat family with a few hundred copies
genome-wide has an expected count below 1 inside a 10 Mb window, so it cannot be
detected there at all. Build an index over the whole genome once, then mask any
subset against it:

```bash
# One pass over the genome (~32 s for T2T-CHM13, ~600 MB peak RAM)
krep index --genome GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --k 18 --sample 16 --min-count 2 --out chm13.k18.s16.kidx

# Mask one chromosome using genome-wide counts (~1.8 s for chr1)
krep mask --genome chm13_chr1_unmasked.fa --index chm13.k18.s16.kidx \
  --index-threshold 4 --graph-gap 250 --min-len 30 --out chr1.bed

# Score against a RepeatMasker .out converted to BED
krep evaluate --truth chr1_rm.bed --pred chr1.bed
```

`--index-threshold` is a **genome-wide** occurrence count. It is on a completely
different scale from the per-slice `--threshold` and must be retuned, not
carried over: the same family that occurs 5 times in a 10 Mb slice occurs
~1500 times across 3.1 Gb.

Expected output with the high-accuracy mode on the 10 Mb mock genome:

```text
Truth bases:      416180
Predicted bases:  416054
Intersection:     414780
Precision:        0.9969
Recall:           0.9966
F1:               0.9968
Per-family recall:
  ALU: 0.9864
  LINE1: 0.9991
  LTR: 0.9941
  MICROSAT: 0.9834
  SAT: 0.9985
  SEG_DUP: 0.9995
  SINE: 0.9664
```

Runtime for this configuration is about **9 seconds** on a laptop.

## Algorithm

### Counting Bloom filter

- **u8 counters**, `k = 4` hash functions.
- Hashing uses double hashing from a single 64-bit key:  
  `h_i(x) = (h1(x) + i * h2(x)) mod m`, where `h1` and `h2` are splitmix64 outputs.
- Number of slots `m = 8 × number_of_kmers`, giving ~0.5 expected increments per counter.

### Masking modes

#### Window density (default)

1. **Insert** every canonical k-mer from the genome into the CBF.
2. **Query** every k-mer; k-mers with estimated count ≥ `threshold` are flagged as repetitive.
3. **Window scoring**: slide a sequence window of length `window` and require at least `density` fraction of the window's k-mers to be repetitive. Mark the whole window as masked.
4. Merge marked bases into BED regions and drop regions shorter than `min_len`.

#### Graph mode (`--graph`)

High-count k-mer positions are treated as nodes in an interval graph. Two positions are linked if they are within `graph_gap` bases of each other. Each connected component is expanded by `k - 1` bases and emitted as a masked region if it spans at least `min_len` bases.

Graph mode is especially good at bridging small low-count gaps inside repeat copies and merging fragmented k-mer hits into coherent repeat blocks, usually giving higher precision than window density.

#### Hamming-distance-1 expansion (`--mismatch1`)

A k-mer is treated as repetitive if it **or any of its Hamming-distance-1 canonical neighbors** is high-count. This recovers diverged repeat copies where mutation has changed one base. To keep the pass fast, high-count seed k-mers are collected into a hash set once; the per-position query then becomes a fast set-membership test instead of many CBF lookups. On the 10 Mb mock genome this adds only ~7 seconds compared with the no-mismatch graph mode.

#### Assembly mode (`--assembly`)

High-count k-mers are linked into connected components via a de Bruijn graph of consecutive genomic k-mers. Components whose total CBF abundance is below `assembly-abundance` are discarded, which can filter isolated noisy k-mers. This mode is provided as an alternative graph-theoretic approach; on the current mock genome the simpler graph-gap mode gives the best F1, but `--assembly` may help on noisier real data.

### Optimizations

- Coverage is tracked with a compact `BitVec` (1 bit per base).
- Window density is maintained with an O(n) sliding counter.
- K-mer counting and querying are parallelized with **rayon** per 1 Mb chunks.
- Mismatch-1 queries use a precomputed `HashSet<u64>` of high-count seed k-mers to avoid thousands of CBF lookups per base.

## CLI reference

### `mock`

| Flag | Default | Description |
|------|---------|-------------|
| `--size` | `10m` | Genome size (`k`/`m`/`g` suffixes supported) |
| `--gc` | `0.5` | Background GC fraction |
| `--seed` | `42` | Random seed |
| `--divergence` | `0.05` | Per-base divergence applied to each repeat copy |
| `--segdup-divergence` | `0.02` | Per-base divergence for segmental-duplication copies |
| `--out` | `genome.fa` | Output FASTA |
| `--bed` | `repeats.bed` | Ground-truth BED |

Repeat families in the mock genome: ALU, LINE1, SINE, LTR, MICROSAT, SAT, **SEG_DUP**.

### `mask`

| Flag | Default | Description |
|------|---------|-------------|
| `--genome` | required | Input FASTA |
| `--k` | `21` | k-mer length |
| `--k-list` | — | Comma-separated k values; union of masks from each k |
| `--threshold` | `3` | Minimum CBF count for a k-mer to be considered repetitive |
| `--window` | `40` | Sequence window length for density scoring |
| `--density` | `0.25` | Minimum fraction of repetitive k-mers in a window |
| `--min-len` | `30` | Shortest reported masked region |
| `--out` | `masked.bed` | Output BED |
| `--soft` | optional | Output soft-masked FASTA (lowercase for masked bases) |
| `--hard` | optional | Output hard-masked FASTA (`N` for masked bases) |
| `--graph` | false | Use graph-based connected-component masking |
| `--graph-gap` | `100` | Max gap (in bases) for linking k-mers in graph mode |
| `--mismatch1` | false | Allow one mismatch against high-count neighbors |
| `--assembly` | false | Use de Bruijn graph assembly of high-count k-mers |
| `--assembly-abundance` | `0` | Minimum component CBF abundance (0 = automatic) |
| `--cbf-factor` | `8` | CBF size multiplier (`slots = factor × kmers`). Lower = less RAM, more collisions |
| `--index` | — | Use a genome-wide index (from `krep index`) instead of counting within the input FASTA |
| `--index-threshold` | `10` | Genome-wide occurrence count required, used with `--index`. Different scale from `--threshold` |

### `index`

Build a genome-wide k-mer count index. Run once per genome, then mask any subset
against it.

| flag | default | meaning |
|------|---------|---------|
| `--genome` | – | input FASTA (the whole genome) |
| `--k` | 18 | k-mer length; needs k >= 17 at 3.1 Gb so random k-mers do not recur by chance |
| `--sample` | 16 | keep 1 in N k-mers by hash (power of two); counts stay exact |
| `--min-count` | 2 | drop k-mers rarer than this genome-wide |
| `--buffer` | 48000000 | k-mers held in RAM before spilling a run (8 bytes each) |
| `--tmp-dir` | `krep_tmp` | scratch space for sorted runs |
| `--out` | `genome.kidx` | output index |
| `--verbose` | off | per-record progress |

The command prints an occurrence-count histogram, which is the practical way to
pick `--index-threshold`.

### `demask`

| Flag | Default | Description |
|------|---------|-------------|
| `--genome` | required | Soft-masked input FASTA |
| `--out` | `demasked.fa` | Uppercase output FASTA |

### `compare-mask`

| Flag | Description |
|------|-------------|
| `--reference` | Original / truth soft-masked FASTA |
| `--predicted` | Predicted soft-masked FASTA |

### `evaluate`

| Flag | Description |
|------|-------------|
| `--truth` | Ground-truth BED (4-column BED enables per-family recall) |
| `--pred` | Predicted BED |

## Tuning tips

- **More repeats detected** (higher recall): lower `--threshold`, smaller `--k`, use `--k-list`, increase `--graph-gap`, or add `--mismatch1`.
- **Fewer false positives** (higher precision): raise `--threshold`, raise `--density`, raise `--min-len`, or use `--graph` with a moderate gap.
- The defaults are tuned for the mock genome with ~5% divergence. Real genomes may need different parameters.
- Segmental duplications are detected well with the defaults because they share long, nearly identical k-mer blocks.

## Performance

On a 10 Mb mock genome:

| Mode | Time | F1 |
|------|------|-----|
| Window density (default) | ~1.4 s | ~0.948 |
| Graph, `k=21` | ~1.3 s | ~0.98 |
| Graph + `--mismatch1`, `k=18`, `threshold=5` | ~9 s | **0.9968** |

- Memory: ~80 MB for the CBF plus a few auxiliary arrays (well under 200 MB total).
- Binary size: ~2.7 MB.

### Real-genome test (T2T-CHM13)

**A note on ground truth.** The lowercase in NCBI's `*_genomic.fna` is *not*
RepeatMasker output. Measured on this assembly:

| mask | genome coverage |
|------|-----------------|
| NCBI `genomic.fna` lowercase | 40.27% |
| RepeatMasker `_rm.out` | 54.19% |

On chr1 only 83% of the lowercase is explained by RepeatMasker, and the
lowercase captures just 62% of RepeatMasker's calls. The NCBI soft-mask is a
WindowMasker/TRF-style de novo mask. Scoring against it measures agreement with
a *different de novo masker*, not with RepeatMasker. Download
`GCF_..._rm.out.gz` and convert it to BED for a homology-based ground truth.

**Full chr1 (248 Mb), genome-wide index, scored against RepeatMasker:**

| `--index-threshold` | `--graph-gap` | precision | recall | F1 |
|---|---|---|---|---|
| 2 | 150 | 0.610 | 0.911 | 0.731 |
| 3 | 250 | 0.678 | 0.850 | 0.754 |
| **4** | **250** | **0.767** | **0.765** | **0.766** |
| 4 | 150 | 0.832 | 0.682 | 0.750 |
| 10 | 150 | 0.939 | 0.568 | 0.708 |

Masking full chr1 takes **1.8 s** once the index is built.

**Why recall plateaus.** Per-family recall at `t=4, gap=250` splits cleanly by
family age:

| recovered almost fully | recall | | largely missed | recall |
|---|---|---|---|---|
| Satellite/acro | 1.000 | | SINE/tRNA | 0.087 |
| Satellite/centr | 0.999 | | DNA/Crypton | 0.122 |
| Retroposon/SVA | 0.993 | | LINE/CR1 | 0.202 |
| SINE/Alu | 0.978 | | LINE/L2 | 0.294 |
| LTR/ERVK | 0.969 | | SINE/MIR | 0.378 |
| Simple_repeat | 0.832 | | LTR/ERVL | 0.424 |
| LINE/L1 | 0.826 | | DNA/hAT-Charlie | 0.507 |

**54% of all missed bases** on chr1 are ancient, highly diverged families
(L2, MIR, CR1, Helitron, DNA transposons). RepeatMasker finds these by aligning
to a curated Dfam consensus; a 30%-diverged MIR copy shares essentially no exact
18-mers with any other copy, so *no* k-mer-abundance method can recover it. This
is a structural ceiling, not a tuning problem — closing it requires consensus
building plus alignment, not a different threshold.

### Memory notes for large genomes

**With `krep index` (recommended).** Memory is bounded by `--buffer` (8 bytes
per buffered k-mer; the 48M default is ~384 MB) plus the largest single record,
regardless of genome size. Sorted runs spill to `--tmp-dir`; put that on a fast
local disk (under WSL, the ext4 root is far faster than `/mnt/c`). For
T2T-CHM13 at `k=18 --sample 16`: 32 s, ~2 GB of temporary files, a 235 MB index
holding 19.6M entries.

Sub-sampling is what makes this fit. Only k-mers with `hash(kmer) % sample == 0`
are tracked, which shrinks the table by `sample` without biasing counts: the
decision depends solely on the k-mer's own hash, so a tracked k-mer is counted at
every occurrence and its stored count is exact. `--graph-gap` must comfortably
exceed the sampling stride so sparse seeds still link (gap 250 with sample 16).

**Without an index (legacy per-record counting).** The CBF is sized
`genome_size x cbf_factor` bytes and is rebuilt per record, so a single human
chromosome needs several GB and a whole genome does not fit at all. Use
`krep index` for anything above a few tens of Mb.

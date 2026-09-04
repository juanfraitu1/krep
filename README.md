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

# Compare krep's soft mask to the original RepeatMasker soft mask
krep compare-mask --reference GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --predicted chm13_krep_soft.fa
```

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

### Real-genome test (T2T-CHM13 chr1, first 10 Mb)

Using the high-accuracy mode on an unmasked slice of the real T2T-CHM13 genome and comparing the krep soft mask back to the original RepeatMasker soft mask:

```text
Bases compared:        10000000
Reference masked:      3211510 (32.12%)
Predicted masked:      2508290 (25.08%)
Both masked (TP):      2212685
Both unmasked (TN):    6492885
Reference only (FN):   998825
Predicted only (FP):   295605
Precision:             0.8821
Recall:                0.6890
F1:                    0.7737
```

Runtime: ~8 seconds for the 10 Mb slice.

krep recovers ~69% of the bases that RepeatMasker soft-masked in this region, with ~88% precision. The de novo approach naturally misses some ancient/diverged repeats and library-specific elements (e.g. SINEs with little sequence identity), while also predicting some novel low-copy repeats RepeatMasker did not annotate. It is a fast, library-free complement to RepeatMasker rather than a byte-for-byte replacement.

### Memory notes for large genomes

The CBF size is roughly `genome_size × cbf_factor` bytes. A full human chromosome (e.g. chr1, ~250 Mb) with the default `--cbf-factor 8` needs ~2 GB of RAM for the CBF alone, plus the sequence and auxiliary structures. On an 8 GB laptop this can fail due to fragmented memory; reduce `--cbf-factor` (try `4` or `2`) or process smaller chromosomes / slices. For whole-genome runs, 16 GB or more is recommended.

Scales roughly linearly with genome size.

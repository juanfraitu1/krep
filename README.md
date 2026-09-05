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

# Optional second index with a spaced seed (weight 16 over 22 bases), densely
# sampled: ~12 min in 16 passes, 560 MB. Don't-care positions tolerate the
# substitutions that diverged copies carry.
krep index --genome GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --seed 1110110110110110110111 --sample 1 --min-count 8 --passes 16 \
  --out chm13.sp16.s1.kidx

# Union both indices, each with its own threshold, plus a density gate.
# Best configuration found on chr1 vs RepeatMasker: P 0.749 / R 0.802 / F1 0.775
krep mask --genome chm13_chr1_unmasked.fa \
  --index chm13.k18.s16.kidx --index chm13.sp16.s1.kidx --index-threshold 4,64 \
  --graph-gap 150 --min-hits 2 --min-density 0.03 --out chr1.bed
```

### Consensus building and library alignment (for diverged families)

K-mer abundance can only find a copy that still carries high-count k-mers. A
copy 30% diverged from its family consensus keeps 0.7^18 ≈ 0.16% of them, so
ancient families (MIR, L2, old L1, DNA transposons) are invisible at any
threshold. RepeatMasker sidesteps this by aligning each copy to a *consensus*,
which is far closer to every copy than copies are to each other. krep can now
do the same without a curated library:

```bash
# Build family consensi de novo from the index's abundant k-mers
# (RepeatScout-style greedy extension; writes a plain sequence dump once)
krep consensus --genome GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --index chm13.k18.s16.kidx --seq-dump chm13.kseq \
  --min-seed-count 50 --max-families 3000 --out chm13_consensi.fa

# Align the consensi back and mask hits; --index and --library union
krep mask --genome chm13_chr1_unmasked.fa --library chm13_consensi.fa \
  --index chm13.k18.s16.kidx --index-threshold 4 --graph-gap 150 --out chr1.bed
```

Library hits carry the consensus id and alignment score in BED column 4.

On a 10 Mb mock genome with copies 15% diverged, library masking alone reaches
P 0.998 / R 0.949 / F1 0.973 (ALU 0.995, LINE1 0.999, LTR 0.997, SINE 0.991),
where k-mer masking alone gets F1 0.48. At 25% divergence it recovers 93% of
ALU copies where k-mer masking recovers 9%.

On T2T-CHM13, `krep consensus` over the whole genome (min seed count 30) built
**1,124 consensi** (2.45 Mb, median 950 bp, longest 6,018 bp — a full L1) in
~30 min. Masking chr1 with that library alone, weight-9 seed and score floor 30,
scores **P 0.959 / R 0.742 / F1 0.836** against RepeatMasker — past the 0.775
ceiling of every k-mer configuration, at higher precision. Per family:

| family | k-mer best | library v2 | | family | k-mer best | library v2 |
|---|---|---|---|---|---|---|
| LINE/L1 | 0.888 | 0.819 | | DNA/hAT-Charlie | 0.593 | 0.491 |
| LTR/ERVK | 0.969 | 0.954 | | SINE/MIR | 0.468 | 0.334 |
| LTR/ERV1 | 0.802 | 0.769 | | LINE/L2 | 0.309 | 0.138 |
| LTR/ERVL | 0.426 | 0.503 | | LINE/CR1 | 0.232 | 0.009 |

The library wins on precision and on ERVL; the k-mer path still wins on the
oldest families because their consensi are the hardest to build (few seeds
survive at 35% divergence) and their copies are short fragments that clear
neither the seed chain nor the score floor. The two approaches union.

**Whole genome.** The final configuration — k18 index at threshold 8, gap 100,
unioned with the v2 library at weight-9 seed and score floor 30 — masks all
24 chromosomes of T2T-CHM13 in **35 minutes** on a 6-core laptop and writes a
soft-masked FASTA:

```bash
krep mask --genome GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --index chm13.k18.s16.kidx --index-threshold 8 --graph-gap 100 \
  --library chm13_consensi_v2.fa --lib-seed 11101001100111 --lib-min-score 30 \
  --min-len 30 --out chm13_krep.bed --soft chm13_krep_soft.fa
```

Scored base-by-base against RepeatMasker over the whole genome
(1,689,365,764 annotated bases):

| | value |
|---|---|
| masked | 1,382,100,035 bp (44.3%) |
| precision | **0.937** |
| recall | **0.767** |
| F1 | **0.844** |

| family | recall | | family | recall |
|---|---|---|---|---|
| SINE/Alu | 0.993 | | DNA/hAT-Charlie | 0.509 |
| Retroposon/SVA | 0.992 | | LTR/ERVL | 0.491 |
| Satellite/centr | 0.998 | | SINE/MIR | 0.330 |
| LTR/ERVK | 0.962 | | LINE/L2 | 0.146 |
| Simple_repeat | 0.884 | | LINE/CR1 | 0.022 |
| LINE/L1 | 0.834 | | RC/Helitron | 0.022 |
| LTR/ERV1 | 0.793 | | | |
| LTR/ERVL-MaLR | 0.750 | | | |

The remaining gap to RepeatMasker is concentrated in the oldest families (L2,
MIR, CR1, Helitron): their consensi are hard to build de novo and their copies
are short fragments at ~65% identity. Everything younger than ~150 My is
recovered at 75–99%.

### Hybrid mode: a curated library where de novo cannot reach

`--library` accepts any consensus FASTA, so the families krep cannot build can
come from Dfam (CC0). Feeding RepeatMasker's own 1,403 human consensi through
krep's aligner is also the cleanest diagnostic available: it separates
consensus quality from aligner sensitivity. On chr1:

| library | P | R | F1 | ERVL | Tip100 | Helitron | MIR | L2 |
|---|---|---|---|---|---|---|---|---|
| krep de novo (v2) | 0.965 | 0.738 | 0.837 | 0.503 | 0.167 | 0.020 | 0.334 | 0.138 |
| Dfam human | 0.994 | 0.806 | **0.890** | 0.805 | 0.618 | 0.526 | 0.376 | 0.273 |
| k18@32 + Dfam + tandem | 0.975 | 0.822 | 0.892 | | | | | |
| Dfam, single-hit gate −4,2 | 0.992 | 0.823 | 0.899 | 0.85 | 0.68 | 0.63 | 0.42 | 0.32 |
| **k18@32 + Dfam single-hit −4,2 + tandem + dust** | 0.971 | 0.838 | **0.899** | | | | 0.42 | 0.32 |

So for the mid-age families the consensus was the limit and Dfam fixes it;
for MIR and L2 it was not — RepeatMasker's own consensi barely move them
through krep's two-hit seed chain. `--lib-single-hit` triggers extension on
every gated seed hit instead: Dfam-only rises to **P 0.990 / R 0.834 /
F1 0.905**, with MIR 0.38→0.44, L2 0.27→0.35, CR1 0.21→0.30 and Helitron
0.53→0.67. Its cost is set by `--lib-gate`, the ungapped flank test that decides
whether a seed hit earns a banded alignment: the sensitivity comes from
short fragments whose flanks run off into random sequence, so a looser gate
finds more of them and runs more alignments that fail. At `-4,2` it costs
about 3× chained mode for most of the gain. Adding krep's de novo consensi on top of Dfam buys +1.7
recall for −2.7 precision (they carry segmental duplications and gene
families); the k-mer index at threshold 32 adds satellites and simple
repeats Dfam lacks at almost no precision cost.

```bash
# Dfam human consensi via the API (JSON-wrapped FASTA; unwrap the "body")
curl -s "https://dfam.org/api/families?clade=9606&clade_relatives=ancestors&format=fasta&limit=5000"
krep mask --genome chr1.fa --index chm13.k18.s16.kidx --index-threshold 32 --graph-gap 100 \
  --library dfam_human.fa --lib-seed 11101001100111 --lib-min-score 30 --tandem --out chr1.bed
```

**Whole genome, hybrid.** k18 index at threshold 32 + Dfam human consensi
in single-hit mode (gate −4,2) + tandem + DUST, all 24 chromosomes in
**41 minutes**, soft-masked FASTA written:

| | library-free (k18@8 + krep v2) | hybrid (k18@32 + Dfam single-hit + tandem + dust) |
|---|---|---|
| masked | 44.3% | 47.6% |
| precision | 0.937 | **0.967** |
| recall | 0.767 | **0.849** |
| F1 | 0.844 | **0.904** |

| family | library-free | hybrid |
|---|---|---|
| SINE/Alu | 0.993 | 0.995 |
| Retroposon/SVA | 0.992 | 0.991 |
| Satellite/centr | 0.998 | 0.997 |
| LTR/ERVK | 0.962 | 0.996 |
| LINE/L1 | 0.834 | 0.917 |
| LTR/ERV1 | 0.793 | 0.944 |
| LTR/ERVL-MaLR | 0.750 | 0.856 |
| Simple_repeat | 0.884 | 0.848 |
| DNA/TcMar-Tigger | 0.738 | 0.870 |
| LTR/ERVL | 0.491 | 0.819 |
| DNA/hAT-Charlie | 0.509 | 0.736 |
| DNA/hAT-Tip100 | 0.199 | 0.660 |
| RC/Helitron | 0.022 | 0.603 |
| SINE/MIR | 0.330 | 0.394 |
| LINE/L2 | 0.146 | 0.310 |
| LINE/CR1 | 0.022 | 0.253 |

### Learned hit filter

RepeatMasker's accept/reject decision is family- and divergence-dependent;
a single global score floor is not. `krep mask --lib-dump hits.tsv` writes
every accepted library hit with its features (score, forward/backward
scores, consensus coverage, flank gate scores, GC), and two scripts in
`scripts/` turn a dump labelled against a reference annotation into a model
for `--lib-model`:

```bash
# permissive dump on a training chromosome, then a test chromosome
krep mask --genome chr2.fa --library dfam_human.fa --lib-single-hit --lib-gate=-4,2 \
  --lib-min-score 15 --lib-dump chr2_hits.tsv --out /dev/null
python3 scripts/train_thresholds.py chr2_hits.tsv chr2_rm_family.bed model_thr.tsv   # per-consensus floors
python3 scripts/train_logistic.py chr2_hits.tsv chr2_rm_family.bed chr1_hits.tsv chr1_rm_family.bed model_thr.tsv model_logit.tsv
krep mask --genome chr1.fa --library dfam_human.fa --lib-single-hit --lib-gate=-4,2 \
  --lib-min-score 15 --lib-model model_logit.tsv --out chr1.bed
```

Trained on chr2, tested on chr1 (Dfam single-hit, gate −4,2):

| filter | P | R | F1 | MIR | L2 | CR1 |
|---|---|---|---|---|---|---|
| global floor 30 | 0.992 | 0.823 | 0.899 | 0.42 | 0.32 | 0.26 |
| per-consensus floors | 0.977 | 0.846 | 0.907 | 0.51 | 0.42 | 0.36 |
| logistic + per-consensus offset | 0.976 | 0.850 | **0.909** | 0.51 | 0.42 | 0.36 |

The learned floors say why: 1,008 of 1,402 consensi want a floor in the 20s,
48 want ≥50 and 30 are best dropped outright. The logistic model's weights
put most of the remaining signal in hit length and per-base score. The
training scripts are plain Python (no NumPy) and take a few minutes on ~1.6M
candidates.

**Whole genome with the learned filter** (hybrid + logistic model trained on
chr2): **P 0.9530 / R 0.8707 / F1 0.9100**, 2354.84 s — the final
configuration; see `HANDOFF.md` for the exact command and outputs.

Whether to use this is a policy choice: with Dfam, krep is no longer library-free.

### NCBI-style masking (WindowMasker + TRF + DUST)

NCBI's `genomic.fna` lowercase is not RepeatMasker; it is WindowMasker (k-mer
frequency), TRF (tandem repeats) and DUST (low complexity). krep's k-mer index
is the WindowMasker analogue; `--tandem` and `--dust` supply the other two:

```bash
krep mask --genome chr1.fa --index chm13.k18.s16.kidx --index-threshold 6 --graph-gap 100 \
  --tandem --dust --dust-threshold 3 --min-len 20 --out chr1_ncbi_style.bed --soft chr1_ncbi_style.fa
krep compare-mask --reference GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --predicted chr1_ncbi_style.fa
```

Against the NCBI lowercase on chr1 this scores P 0.853 / R 0.732 / F1 0.788
(the RepeatMasker-oriented configuration scores 0.757 against it). Tandem
detection is k-mer recurrence at a fixed period verified by periodic
self-identity; DUST is triplet-frequency skew over a 64-bp window. RepeatMasker's
Simple_repeat / Low_complexity classes are only half recovered by these
(0.50 / 0.30) because many are low-complexity with a *wobbling* period.

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
| `--index` | — | Use a genome-wide index (from `krep index`) instead of counting within the input FASTA. Repeatable; hits from all indices are unioned |
| `--index-threshold` | `10` | Genome-wide count at which a seed may *nucleate* a region. One value, or one per `--index` (comma-separated) — a spaced seed and a contiguous k-mer sit on different count scales |
| `--index-threshold-low` | = high | Lower count at which a seed may *extend or bridge* a region a stronger seed nucleated (hysteresis). Weak seeds never start a region |
| `--min-hits` | `1` | Minimum seed hits in a linked component |
| `--min-density` | `0` | Minimum hits per base over a linked component. With dense sampling, real copies sit well above 0.1; background clusters do not |
| `--library` | — | Consensus FASTA (from `krep consensus`, or any repeat library). Consensi are seeded with a weight-11 spaced seed on both strands; two hits on one diagonal trigger a banded X-drop alignment, and hits at or above `--lib-min-score` are masked. Unions with `--index` |
| `--lib-min-score` | `50` | Minimum alignment score (match +1, mismatch −1, gap −3) |
| `--lib-band` | `16` | Alignment band half-width |
| `--lib-xdrop` | `40` | Stop extending once the score falls this far below its best |
| `--lib-seed` | PatternHunter w11 | Spaced seed for library hits (weight 8–13, need not be symmetric). Weight 9 (`11101001100111`) is the sensitivity/speed knee on human |
| `--lib-single-hit` | off | Trigger extension on every gated seed hit instead of two hits on one diagonal; more sensitive to short diverged fragments |
| `--lib-dump` | — | Write every accepted library hit with its features to a TSV, for training `--lib-model` |
| `--lib-model` | — | Learned filter: per-consensus minimum scores (`name<TAB>score`) or a logistic model (`#logistic` header) from `scripts/` |
| `--lib-gate` | `4,6` | Single-hit gate "SUM,SIDE" on the ungapped flank score (32 bp each side of the seed, seed excluded; random flanks average −16 per side). This is the sensitivity/speed dial: on chr1 with Dfam, `4,6` 69 s / F1 0.890, `0,4` 101 s / 0.895, `-4,2` 190 s / 0.899, `-8,0` 432 s / 0.902, `-12,-2` 921 s / 0.904. Negative values need the `--lib-gate=-4,2` form |
| `--tandem` | off | Tandem repeats: k-mer (`--tandem-k`, 5) recurrence at a fixed period ≤ `--tandem-max-period` (500), kept if dense (`--tandem-density`, 0.25) and periodic identity ≥ `--tandem-identity` (0.7), length ≥ `--tandem-min-len` (20) |
| `--dust` | off | Low complexity: DUST triplet-skew score over `--dust-window` (64) bases above `--dust-threshold` (5; random ~0.5, (CA)n ~15, poly-A ~31) |

### `index`

Build a genome-wide k-mer count index. Run once per genome, then mask any subset
against it.

| flag | default | meaning |
|------|---------|---------|
| `--genome` | – | input FASTA (the whole genome) |
| `--k` | 18 | contiguous k-mer length; needs k >= 17 at 3.1 Gb so random k-mers do not recur by chance |
| `--seed` | – | spaced-seed pattern of 1/0 (care/don't-care), e.g. `1110110110110110110111`; must be symmetric; overrides `--k`. Weight must stay >= 16 at 3.1 Gb or composition background swamps the counts |
| `--passes` | 1 | stream the genome N times, each counting one hash partition; bounds temp disk and RAM for dense sampling |
| `--sample` | 16 | keep 1 in N k-mers by hash (power of two); counts stay exact |
| `--min-count` | 2 | drop k-mers rarer than this genome-wide |
| `--buffer` | 48000000 | k-mers held in RAM before spilling a run (8 bytes each) |
| `--tmp-dir` | `krep_tmp` | scratch space for sorted runs |
| `--out` | `genome.kidx` | output index |
| `--verbose` | off | per-record progress |

The command prints an occurrence-count histogram, which is the practical way to
pick `--index-threshold`.

### `consensus`

Build de novo repeat consensi (RepeatScout-style) from a contiguous k-mer
index. Seeds are taken most-abundant first; for each seed not already inside a
built consensus, up to `--max-occ` genomic occurrences are fetched, oriented by
strand, and the consensus is grown one base at a time by a banded
fit-alignment DP over all occurrences. Lanes are split into voters (choose the
base) and held-out judges (decide whether the total improved), which keeps the
extension from fitting noise; growth stops after `--lookahead` steps without
improvement. Consensi are filtered for length, support, tandem periodicity and
redundancy.

| flag | default | meaning |
|------|---------|---------|
| `--genome` | – | input FASTA (the whole genome) |
| `--index` | – | contiguous k-mer index (no `--seed`) for picking seeds |
| `--seq-dump` | `genome.kseq` | plain sequence dump for random access; created once, reused |
| `--min-seed-count` | 20 | only k-mers this frequent become seeds |
| `--max-families` | 2000 | stop after this many consensi |
| `--max-occ` | 100 | occurrences sampled per seed |
| `--flank` | 3000 | bases fetched on each side of a seed; caps consensus length |
| `--band` | 16 | alignment band half-width (tolerated net indel drift) |
| `--lookahead` | 100 | stop after this many steps without score improvement |
| `--min-len` | 50 | discard shorter consensi |
| `--min-support` | 10 | discard consensi with fewer copies fitting at ~62% identity |
| `--out` | `consensi.fa` | output FASTA |

Scoring is match +1 / mismatch −1 / gap −3. Two details decide whether this
works at all. The gap cost: at −2 a diagonal switch costs the same as a
mismatch and the DP drifts *upward* through random sequence (the linear phase
of alignment statistics), so extension never stops. And the voter/judge split:
choosing each base by plurality over the same lanes that score the step lets
the consensus fit noise — with few lanes the total creeps up through random
tails — so a third of the lanes are held out purely to judge improvement.

Seeds are checked against the Hamming-2 neighbourhood of every k-mer in
already-built consensi before extension; without that, subfamily seeds (one
or two mismatches from a consensus already built) rebuild the same family
thousands of times and the run never reaches low-count seeds.

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

**Closing the gap: what was tried, with numbers.** Three additions were
measured on full chr1 against RepeatMasker, all starting from the k18 baseline
(`t=4 gap=250`, F1 0.766):

| change | best config | P | R | F1 | verdict |
|---|---|---|---|---|---|
| Hysteresis (`--index-threshold-low`) alone | high 8 / low 3 / gap 250 | 0.702 | 0.834 | 0.762 | no gain: the FP source is *bridging*, not nucleation; weak seeds bridge as freely as strong ones |
| Spaced seed (w16/s22, `--sample 1`) alone | t 64 / gap 150 / density 0.05 | 0.832 | 0.724 | 0.774 | same F1 as k18, different profile; weight had to stay 16 and threshold rise 16x to beat composition background |
| **k18@4 ∪ spaced@64 + density gate** | gap 150 / `--min-density 0.03` | **0.749** | **0.802** | **0.775** | best found; +3.7 recall points |

Per-family recall, k18 baseline → union:

| family | before | after | | family | before | after |
|---|---|---|---|---|---|---|
| Simple_repeat | 0.832 | 0.938 | | SINE/MIR | 0.378 | 0.468 |
| SINE/Alu | 0.978 | 0.999 | | DNA/hAT-Charlie | 0.507 | 0.593 |
| LINE/L1 | 0.826 | 0.888 | | LINE/CR1 | 0.202 | 0.232 |
| DNA/TcMar-Tigger | 0.741 | 0.804 | | LINE/L2 | 0.294 | 0.309 |

Every configuration converges on F1 ≈ 0.775. That is the ceiling of
copy-vs-copy k-mer abundance on this genome, and it is structural:
RepeatMasker aligns each copy to a curated Dfam *consensus*, which is far
closer to every copy than any two copies are to each other. A 30%-diverged
MIR copy shares essentially no exact seed with any other copy at the weight a
3.1 Gb genome demands. L2 barely moved (0.294 → 0.309) no matter what. Closing
the remaining gap needs consensus building plus alignment — or a library — not
another threshold.

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

**Dense sampling.** `--sample 1` with a spaced seed produces 3.1e9 seed
instances; `--passes 16` keeps each pass at ~2.3 GB of temp and one buffer of
RAM. Retained entries are held in RAM across passes (12 bytes each), so pair
dense sampling with a higher `--min-count` (8 gave 47M entries / 560 MB).

**Without an index (legacy per-record counting).** The CBF is sized
`genome_size x cbf_factor` bytes and is rebuilt per record, so a single human
chromosome needs several GB and a whole genome does not fit at all. Use
`krep index` for anything above a few tens of Mb.

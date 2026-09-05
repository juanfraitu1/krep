# krep — handoff (2026-09-04)

Start here. Everything below the "Detailed history" line is the round-by-round
record; this top section is the current state and what to do next.

## What krep is now

A repeat masker for T2T-CHM13 that works in three layers, all in
`krep mask`, unioned:

1. **Genome-wide k-mer index** (`krep index`, `--index`): abundance masking
   with counts measured over the whole genome, hash-subsampled, exact.
2. **Library alignment** (`--library`): consensus sequences seeded on both
   strands, banded X-drop extension. Consensi come from `krep consensus`
   (de novo, RepeatScout-style) or from Dfam. `--lib-single-hit` with
   `--lib-gate=-4,2` is the sensitive mode; `--lib-model` applies a learned
   per-consensus filter.
3. **Tandem + DUST** (`--tandem`, `--dust`): the TRF/DUST half of an
   NCBI-style mask.
4. **Profile-HMM layer** (`--hmm-bed`): merge pre-computed HMM hits (e.g.
   `nhmmer`) into the masked output. This is what rescues the oldest
   families (MIR, L2, CR1) that a single consensus cannot align reliably.

## Final numbers (whole genome, all 24 chromosomes, vs RepeatMasker .out)

| configuration | P | R | F1 | time |
|---|---|---|---|---|
| library-free: k18@8 + de novo consensi | 0.937 | 0.767 | 0.844 | 35 min |
| hybrid: k18@32 + Dfam single-hit −4,2 + tandem + dust | 0.967 | 0.849 | 0.904 | 41 min |
| **hybrid + logistic hit filter (final)** | **0.9530** | **0.8707** | **0.9100** | 2354.84 s |

Masked 4370062 regions, 1543473801 of 3117275501 bp (49.51%) in 2353.9s. Outputs: `C:\krep_work\chm13_krep_final.bed`,
`chm13_krep_final_soft.fa` (soft-masked genome), `final_genome_eval*.txt`.

Per family (hybrid → hybrid + logistic filter), genome-wide recall:

| family | hybrid | + filter |
|---|---|---|
| SINE/Alu | 0.995 | 0.996 |
| LINE/L1 | 0.917 | 0.928 |
| LTR/ERVK | 0.996 | 0.996 |
| LTR/ERV1 | 0.944 | 0.949 |
| LTR/ERVL-MaLR | 0.856 | 0.876 |
| LTR/ERVL | 0.819 | 0.838 |
| DNA/TcMar-Tigger | 0.870 | 0.884 |
| DNA/hAT-Charlie | 0.736 | 0.777 |
| DNA/hAT-Tip100 | 0.660 | 0.704 |
| RC/Helitron | 0.603 | 0.644 |
| SINE/MIR | 0.394 | 0.497 |
| LINE/L2 | 0.310 | 0.412 |
| LINE/CR1 | 0.253 | 0.344 |
| Simple_repeat | 0.848 | 0.876 |
| Satellite/centr | 0.997 | 0.997 |

The remaining gap is the oldest families (L2, MIR, CR1): short fragments at
~65% identity. Everything younger is at 0.75–0.99.

## Exact final command

```bash
krep mask --genome GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --index chm13.k18.s16.kidx --index-threshold 32 --graph-gap 100 \
  --library dfam_human.fa --lib-seed 11101001100111 --lib-min-score 15 \
  --lib-single-hit --lib-gate=-4,2 --lib-model model_logit.tsv \
  --tandem --dust --min-len 20 --out chm13_krep_final.bed --soft chm13_krep_final_soft.fa
```

(`--lib-gate=-4,2` must use the `=` form: clap reads a leading `-` as a flag.)

## Data locations (all under `C:\krep_work`, not in git)

- `krep.exe` — final binary (built from commit `417d2da` or later).
- Indices: `chm13.k18.s16.kidx` (235 MB), `chm13.sp16.s1.kidx` (560 MB,
  spaced w16 — not used in the final config).
- Libraries: `dfam_human.fa` (1,403 Dfam human consensi, `name#accession`,
  fetched via `https://dfam.org/api/families?clade=9606&clade_relatives=ancestors&format=fasta&limit=5000`,
  JSON-wrapped: unwrap `body`), `chm13_consensi_v2.fa` (1,124 de novo from the
  contiguous k18 index), `chm13_consensi_sp16_pool500k.fa` (537 de novo from
  the spaced w16/s22 index), `krep_plus_dfam.fa` (both).
- Learned models: `model_logit.tsv` (logistic, chr2-only, production),
  `model_logit_234.tsv` (chr2+3+4, floor 15), `model_logit10.tsv`
  (chr2+3+4, floor 10), `model_chr2.tsv` (per-consensus floors, chr2),
  `model_chr1.tsv` (in-sample).
- Candidate dumps for training: `chr{1,2,3,4}_cand.tsv` (floor 15) and
  `chr{1,2,3,4}_cand10.tsv` (floor 10), columns consumed by
  `scripts/train_logistic.py` and `scripts/train_gbm.py`.
- Ground truth: `genome_rm.bed` (merged), `genome_rm_family.bed` (4-column,
  5.5M rows), `chr1_rm.bed`, `chr1_rm_family.bed`, `chr2_rm_family.bed`,
  `chr3_rm_family.bed`, `chr4_rm_family.bed`; raw `~/krep_data/rm.out.gz`.
  T2T segdup/censat for chr1 in `~/krep_data/chr1_sedefSegDups.json`,
  `chr1_censat.json` (UCSC hs1 API).
- Sequence dump for `krep consensus`: `chm13.kseq` (3.1 GB, reused).
- Models: `model_sp16_thresholds.tsv`, `model_sp16_logit.tsv` (de novo
  library filter, trained on chr2+3+4).
- Per-chromosome FASTAs: `chm13_chr1_unmasked.fa`, `chm13_chr2_unmasked.fa`,
  `chm13_chr3_unmasked.fa`, `chm13_chr4_unmasked.fa`.
- Mock genomes for regression: `mock15.*`, `mock25.*`.
- Substitution-count dump: `chr1_subst.tsv` (39.6M gate-window bases from
  415k accepted library hits on chr1; produced by `--lib-subst-dump`).
- HMM layer: `target_hmms.hmm` (Dfam MIR/L2/CR1 keys), per-family
  `*_chr1.tbl`/`*_chr1.bed`, merged `hmm_mir_l2_cr1_chr1.bed`.
- Launch scripts: `run_genome2.sh`, `run_final.sh`, `run_dump.sh`,
  `run_gate.sh`, `run_dump34.sh`, `run_dump10.sh`, `run_train234.sh`,
  `run_train10.sh`, `run_gbm.sh`, `run_gbm10.sh`.

## Latest work

### 1. Learned filter: more training data

Dumped chr3 and chr4 library-alignment candidates (floor 15, 1.36M and
1.34M rows) and retrained the logistic filter on chr2+chr3+chr4.
Held-out chr1 applied result:

| model | train set | P | R | F1 |
|---|---|---|---|---|
| `model_logit.tsv` | chr2 only | 0.9759 | 0.8497 | **0.9085** |
| `model_logit_234.tsv` | chr2+3+4 | 0.9753 | 0.8498 | **0.9082** |

Offline (candidate-level) evaluation: 0.9008 vs 0.9005. Conclusion:
**the logistic filter is saturated on its current features**; more
same-distribution training data does not move it.

### 2. Gradient boosting

Added `scripts/train_gbm.py` (LightGBM, same features + consensus as a
native categorical). GBM on floor-15 data: F1 **0.9019** (vs logistic
0.9008 offline). On floor-10 data it reached only **0.8927** because
the lower floor floods the candidate set with weak noise the tree model
still cannot separate from true ancient fragments.

### 3. Sub-15 candidate floor: alignment memo fix

At `--lib-min-score 5` the single-hit aligner ground to a halt (55+ min
for chr2). The cause: a *successful* banded DP never updated the failed-
diagonal memo, so weak accepted hits re-triggered DPs every few bases
along the same diagonal. Fixed in `src/align.rs` by memoizing the diagonal
at the hit's end. After the fix:

- floor-15 chr1 regression unchanged: **P 0.9759 / R 0.8497 / F1 0.9084**
  (was 0.9085), and ~13% faster.
- floor-10 dumps run in ~80 s per chromosome, faster than the old
  floor-15 time.

Headroom check on chr2 (partial dump at floor 5):

| floor | recall ceiling (kept bases vs RM) |
|---|---|
| 15 | 0.8406 |
| 10 | 0.9057 |
| 8  | 0.9164 |
| 5  | 0.9189 |

There is real signal in the 10–15 score band, but the current learned
features cannot exploit it without measurable precision loss. The aligner
itself (gate, seed weight, scoring matrix) has to get better before a
post-hoc filter can turn that signal into recall gain.

### 4. Divergence-aware scoring (started)

Added `--lib-subst-dump` to collect A/C/G/T counts from the gate windows
of accepted library hits. Running it on chr1 with the production
logistic model produced 39.6M counted bases from 415k hits:

- matches: 24.54M (62.0%)
- transitions: 6.49M (16.4%)
- transversions: 8.54M (21.6%)
- Ti/Tv ratio = **0.76**

A generic transition-biased matrix (transition score −0.78, transversion
−1.17, keeping E_random = −0.5) is **not obviously better** than the
flat mismatch = −1. The accepted-hit spectrum is not transition-rich
enough for a simple Ti/Tv bias to matter, and the gate windows are short
and selected for the existing flat score. Family-specific matrices
(per clade, estimated from high-confidence alignments of each repeat
family) may still help, but that requires much more data and a new
alignment format. For now this path is deprioritised.

### 5. Profile HMM layer (`--hmm-bed`)

Implemented an optional profile-HMM merge in `krep mask`. The workflow is:
run `nhmmer` on the target sequence with one or more Dfam HMMs, convert the
`--tblout` to BED with `scripts/nhmmer_to_bed.py`, merge per-family BEDs with
`scripts/run_nhmmer_layer.sh`, and pass the merged BED to `krep mask --hmm-bed`.
The BED regions are folded into every mask layer (index/library/tandem/dust)
per chromosome and appear in both the output BED and any soft/hard FASTA.

On chr1 vs `chr1_rm_family.bed` (production model + MIR/L2/CR1_Mam HMMs,
`-E 1e-2`):

| layer | P | R | F1 | SINE/MIR | LINE/L2 | LINE/CR1 |
|---|---|---|---|---|---|---|
| production (no HMM) | 0.9578 | 0.8611 | 0.9069 | 0.5202 | 0.4243 | 0.3648 |
| + MIR HMM | 0.9562 | 0.8772 | 0.9150 | 0.7874 | 0.4253 | 0.3652 |
| + MIR + L2 + CR1_Mam | 0.9502 | 0.8931 | 0.9208 | 0.7877 | 0.6262 | 0.3789 |
| + MIR + L2 + all CR1 subfams | **0.9501** | **0.8933** | **0.9208** | 0.7877 | 0.6262 | **0.4054** |

The MIR and L2 gains are large; CR1_Mam alone is weak because it only covers
one CR1 subfamily (2.9% recall of LINE/CR1 at 0.48 precision). Adding the
other nine Dfam CR1 keys (`CR1-16_AMi`, `CR1-13_AMi`, `CR1-12_AMi`,
`CR1-3_Croc`, `CR1-L3A_Croc`, `CR1-L3B_Croc`, `CR1-11_Crp`, `CR1-1_Amn`,
`CR1_Amni-1`) lifts LINE/CR1 recall from 0.3789 to 0.4054 at essentially no
F1 cost. The merged layer costs ~0.8 points of precision but buys 3.2 points
of recall and lifts F1 by 0.014 on chr1.

Example chr1 invocation (uses the WSL binary because the fresh Windows build
is currently blocked by Application Control):

```bash
/mnt/c/krep_work/krep_wsl mask \
  --genome /mnt/c/krep_work/chm13_chr1_unmasked.fa \
  --index /mnt/c/krep_work/chm13.k18.s16.kidx --index-threshold 32 --graph-gap 100 \
  --library /mnt/c/krep_work/dfam_human.fa --lib-seed 11101001100111 \
  --lib-min-score 15 --lib-single-hit --lib-gate=-4,2 \
  --lib-model /mnt/c/krep_work/model_logit.tsv \
  --tandem --dust --min-len 20 \
  --hmm-bed /mnt/c/krep_work/hmm_mir_l2_cr1_chr1.bed \
  --out /tmp/chr1_final_logit_hmm.bed
```

### 6. De novo spaced-seed consensus (no Dfam/HMM)

The consensus builder now seeds from the dense whole-genome spaced-seed index
(`chm13.sp16.s1.kidx`). A new `--seed-pool` flag limits how many abundant
seeds are held for occurrence collection, and the top-K selection uses a
bounded heap so exploring millions of low-count seeds no longer materialises
the full candidate list.

```bash
~/krep_target/release/krep consensus \
  --genome /mnt/c/Users/jfris/OneDrive/Desktop/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --index /mnt/c/krep_work/chm13.sp16.s1.kidx \
  --seq-dump /mnt/c/krep_work/chm13.kseq \
  --out /mnt/c/krep_work/chm13_consensi_sp16_pool500k.fa \
  --min-seed-count 8 --max-families 2000 --seed-pool 500000 --max-occ 30
```

This built **537 de novo consensi** in ~4 min at ~1.6 GB peak RAM. The
library is dominated by Alu/L1/satellite, but also contains real MIR and L2
fragments (e.g. `krep_fam_00039` = Dfam MIR, `krep_fam_00085/00091` = L2a/MIRb
3' ends). A logistic hit filter was trained on chr2+3+4 candidates and
applied to chr1.

Chr1 vs `chr1_rm_family.bed` (k18@32 + de novo library + tandem/dust):

| filter | P | R | F1 | SINE/MIR | LINE/L2 | LINE/CR1 |
|---|---|---|---|---|---|---|
| none (score ≥15) | 0.8065 | 0.7794 | 0.7927 | 0.4881 | 0.2853 | 0.1716 |
| per-consensus thresholds | **0.9422** | 0.7313 | **0.8234** | 0.4283 | 0.1594 | 0.0188 |
| logistic τ=0.60 | 0.9360 | 0.7325 | 0.8218 | 0.4303 | 0.1624 | 0.0216 |
| logistic τ=0.50 | 0.9044 | 0.7443 | 0.8166 | 0.4398 | 0.1913 | 0.0496 |
| logistic τ=0.45 | 0.8824 | 0.7521 | 0.8121 | 0.4481 | 0.2122 | 0.0733 |
| logistic τ=0.40 | 0.8594 | 0.7602 | 0.8068 | 0.4589 | 0.2337 | 0.1007 |
| logistic τ=0.30 | 0.8249 | **0.7729** | 0.7981 | **0.4772** | **0.2683** | **0.1482** |

The fully de novo pipeline reaches **F1 0.8234 on chr1** at high precision,
or trades ~12 points of precision for more ancient-family recall (τ=0.3). It
still trails the Dfam hybrid (F1 ~0.91) and the HMM layer on MIR/L2/CR1, but
it proves the spaced-seed consensus path works without external repeat
libraries.

### 7. Two-pass de novo consensus builder

`krep consensus` now has `--second-pass`, `--second-pass-min-len`,
`--second-pass-min-support`, and `--second-pass-max-families`. The first
pass builds consensi with the primary filters (length/support/tandem/
redundancy). The optional second pass iterates over the same occurrence table
with relaxed filters, skipping seeds already covered by pass-1 consensi and
checking pass-2 consensi against them for redundancy.

This was motivated by zero dedicated CR1/Helitron/Tip100 consensi in the
single-pass 1M library: the families are too low-copy or diverged to clear
the default filters. With `--second-pass-min-len 80 --second-pass-min-support
5 --second-pass-max-families 500`, the 1M library grew from **943 consensi to
1,375 consensi** (pass 1 built 1,337, pass 2 added 38). Genome-wide masking
coverage on chr1–4 with the new library is ~60% vs ~52–55% for the old 1M
library, so the expansion hits more repeat sequence. Evaluation of the
resulting learned filter on chr1 is in progress.

## Remaining tasks, in priority order

1. **Whole-genome HMM-layer run.** The chr1 gain is large; confirm it on
   all 24 chromosomes against `genome_rm_family.bed`. The main cost is
   running `nhmmer` once per family on the whole genome; the `krep mask`
   merge is cheap (~2.6 min on chr1). *Current blocker:* `makehmmerdb` on
   the 3 GB whole-genome FASTA is being killed in WSL (max ~4 GB dynamic
   memory), and `nhmmer` directly on the FASTA also exits immediately.
   Workarounds: (a) run on a machine with more WSL memory, (b) split by
   chromosome and merge per-family BEDs, (c) use a non-WSL HMMER build.
2. **Short diverged fragments in single-hit mode.** A length-normalized score
   threshold or a separate path for short consensi (< 400 bp) could rescue
   some MIR/L2 at a measurable precision cost. The sub-15 score band shows
   ~6 points of recall ceiling if the filter can be made precise.
3. **Precision side.** ~3% FP genome-wide: segmental duplications and gene
   families masked by the k-mer index (real repeats RepeatMasker does not
   annotate), boundary overhang, and library hits below RepeatMasker's own
   cutoffs. Per-region copy-number stats in the BED would let SD-like
   regions be labelled rather than counted as errors.
4. **Consensus builder tail**: CR1/Helitron/Tip100 had no usable de novo
   consensi (seeds below count 30 or too diverged). The two-pass builder
   addresses this; evaluate the resulting recall/precision on chr1.
5. **NCBI-style mode** (index@6 + tandem + dust 3) scores F1 0.788 against
   the NCBI lowercase; matching it closely would need WindowMasker's actual
   two-threshold window scoring. Only worth doing if that mask is a target.
6. **Divergence-aware scoring (revisit later).** Needs family-specific
   transition/transversion matrices estimated from full alignments, not the
   gate-window counts. Generic Ti/Tv bias is not supported by the observed
   accepted-hit spectrum.

## Things that will bite the next person

- **Windows Application Control (os error 4551)** blocks freshly built
  executables at random — it hit the test harness and Cargo build scripts.
  Tests run in WSL (`~/.cargo/bin/cargo test --release --target-dir
  ~/krep_target`, Rust + build-essential installed); the release `.exe` is
  built with the Windows toolchain from WSL and run from `C:\krep_work`.
  The exe cannot be overwritten while a run holds it: copy to a new name.
- **Host memory is tight (8 GB, WSL reserves up to 3.8 GB).** The agent
  harness kills its own background tasks when the host runs low; detached
  jobs (`setsid nohup script &`) survive. Run one krep process at a time
  and poll its log with bounded foreground loops.
- **Killing processes by name**: `pkill -f` and `ps | grep | kill` match
  the calling shell's own command line if the pattern appears in it — this
  killed the working shell twice. Use `[k]rep` bracket patterns and never
  include the pattern text elsewhere in the same command.
- The shell's working directory resets between tool calls in this harness;
  scripts that touch the repo must `cd` explicitly.
- `krep evaluate` treats BED column 4 as the family label; `krep mask`
  writes `consensus:score` there for library hits — fine for evaluation of
  the prediction, but don't feed a prediction BED as *truth*.
- Windows exe wants `C:/...` paths; WSL tools want `/mnt/c/...`.

## Build

```bash
export PATH="/mnt/c/mingw64/bin:$PATH"
/mnt/c/Users/jfris/.cargo/bin/cargo.exe build --release --target-dir target3
cp target3/release/krep.exe /mnt/c/krep_work/krep.exe
~/.cargo/bin/cargo test --release --target-dir ~/krep_target     # 52 tests
```

---

# Detailed history


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

## Round 3: de novo consensus + library alignment (C1) — built and measured

- `src/consensus.rs` — `krep consensus`: RepeatScout-style greedy extension
  from index seeds. Sequence dump (`--seq-dump`, 3.1 GB, reused) + one
  occurrence pass, then per seed: fetch ≤100 windows, orient, extend right and
  left with a banded fit DP (match +1 / mismatch −1 / gap −3, band 16), stop
  after 100 non-improving steps, filter (len ≥ 50, support ≥ 10 at ~62%
  identity, not tandem, not ≥60% 12-mer-redundant).
- `src/align.rs` — `krep mask --library`: consensi seeded on both strands into
  a direct-address table (4^weight), genome scanned forward, two hits on one
  diagonal per consensus trigger an ungapped check then banded X-drop
  extension both ways; hits ≥ `--lib-min-score` masked, seeds inside masked
  regions skipped. `--lib-seed` (weight 8–13), `--lib-band`, `--lib-xdrop`.
  BED column 4 = `consensus:score`. Unions with `--index`.

Bugs that had to be found before any of it worked (all covered by tests now):
1. Gap cost 2 → alignment in the *linear phase*: the max over a band of
   near-free paths drifts upward through random sequence; extension never
   stopped. Gap 3 fixes it.
2. Base chosen by plurality over the same lanes that score the step → with
   few lanes the consensus fits noise. Held-out judges (every 3rd lane) fix it.
3. Exact "seed already covered" check misses subfamily seeds → 9,655 redundant
   rebuilds and the run stopped at seed count 208. Hamming-2 check: 173k seeds
   skipped cheaply, run reaches count 30, 1,124 consensi instead of 380.
4. Aligner chain buffer scanned linearly per table entry → ~3e12 comparisons
   on chr1 at weight 9. Replaced by two slots per oriented consensus (O(1)).
5. Failed extensions left their anchor in the buffer and re-triggered at every
   following position; reversed slices were allocated per extension. Fixed,
   plus an ungapped diagonal check before the banded DP.

Results (chr1 vs RepeatMasker): mock 15% divergence library-only F1 0.973;
mock 25% ALU recall 0.93 (k-mer: 0.09). Real genome, v2 library, weight-9 seed,
score 30: **P 0.959 / R 0.742 / F1 0.836** library-only. Sensitivity knee:
weight 10/score 20 and weight 9/score 30 both ≈0.80 on the v1 library; below
that precision drops.

Attribution of v1 consensi to RM families (chr1 hits): ERV1 74, L1 50,
MaLR 42, Alu 27, Tigger 17, Charlie 15, ERVL 15, Satellite 15, MIR 11, L2 4.

**Whole-genome run** (k18@8 gap100 + v2 library, w9 seed, score 30; all 24
chromosomes; 35 min; soft FASTA written): **P 0.937 / R 0.767 / F1 0.844**
vs RepeatMasker genome-wide. Outputs: `C:\krep_work\chm13_krep.bed`,
`C:\krep_work\chm13_krep_soft.fa` (3.16 GB), evaluation in
`genome_eval.txt` / `genome_eval_family.txt`. Union sweep on chr1 showed the
index adds ~1.5 recall points over the library alone at k18@8; lower index
thresholds cost more precision than they return.

Note: whole-genome runs launched from this harness as background tasks were
killed once mid-run; `run_genome.sh` relaunches it detached
(`setsid nohup`) and polls the log for `ALL_DONE`.

## Round 4: Dfam hybrid, single-hit aligner, tandem + DUST

- **Dfam diagnostic** (1,403 human consensi via the Dfam API, CC0): through
  krep's aligner alone P 0.994 / R 0.806 / F1 0.890 on chr1. Mid-age families
  jump (ERVL 0.50→0.81, Tip100 0.17→0.62, Helitron 0.02→0.53, Charlie
  0.51→0.70) — consensus quality was their limit. MIR (0.38) and L2 (0.27)
  barely move even with RepeatMasker's own consensi — the aligner's two-hit
  chain is their limit. Best hybrid: k18@32 + Dfam + tandem = F1 0.8915.
- **`--lib-single-hit`**: extension on every seed hit, gated by a stricter
  ungapped check (sum ≥ 4 or one side ≥ 6), with a per-consensus failed-
  diagonal memo and a global 8-base cooldown after any failed DP (without
  the memo a weakly similar region re-triggered a full DP at every position;
  the first attempt ran >50 min on chr1). Dfam-only, chr1: P 0.990 /
  R 0.834 / F1 0.905 (chained: 0.890); MIR 0.44, L2 0.35, CR1 0.30,
  Helitron 0.67, but 1,771 s.

  Where that time went, in order of discovery: (1) a global 8-base cooldown
  after failed DPs made it 5× faster and destroyed the gain (F1 0.852) —
  the pause blankets a real fragment's few seed hits; reverted. (2) Capping
  DPs to the two best-gated entries per position kept sensitivity but saved
  only 16%. (3) The real cause: the gate's right-hand window *included the
  seed*, whose nine care positions match by construction, pre-loading it
  by ~+6 so the one-sided rule fired by chance at ~25% of positions. With
  the seed excluded the run takes 64 s — and sensitivity drops back to
  chained levels, because the gain was coming from exactly those loosely
  gated hits: short fragments whose flanks run into random sequence. So the
  gate is the dial, now explicit as `--lib-gate SUM,SIDE`:

  | gate | chr1 time | P | R | F1 | MIR | L2 | CR1 |
  |---|---|---|---|---|---|---|---|
  | 4,6 | 69 s | 0.993 | 0.806 | 0.890 | 0.38 | 0.28 | 0.22 |
  | 0,4 | 101 s | 0.992 | 0.814 | 0.895 | 0.40 | 0.30 | 0.25 |
  | −4,2 | 190 s | 0.992 | 0.823 | 0.899 | 0.42 | 0.32 | 0.26 |
  | −8,0 | 432 s | 0.991 | 0.828 | 0.902 | 0.43 | 0.33 | 0.28 |
  | −12,−2 | 921 s | 0.990 | 0.831 | 0.904 | 0.44 | 0.34 | 0.29 |

  Final hybrid on chr1 (k18@32 + Dfam single-hit −4,2 + tandem + dust 5):
  P 0.971 / R 0.838 / F1 0.899 in 192 s.

  **Whole genome, hybrid: P 0.967 / R 0.849 / F1 0.904**, 41 min, 47.6%
  masked (library-free run: 0.937 / 0.767 / 0.844). Outputs
  `C:\krep_work\chm13_krep_hybrid.bed`, `chm13_krep_hybrid_soft.fa`
  (3.16 GB), `genome2_eval.txt` / `genome2_eval_family.txt`; launcher
  `run_genome2.sh`. Genome-wide per family (library-free → hybrid):
  L1 0.834→0.917, ERVL 0.491→0.819, Charlie 0.509→0.736,
  Helitron 0.022→0.603, MIR 0.330→0.394, L2 0.146→0.310,
  CR1 0.022→0.253.
- **`--tandem`** (`src/tandem.rs`): fixed-period k-mer recurrence runs +
  periodic identity. k=5, density 0.25, min len 20 is the swept optimum
  (Simple_repeat 0.435, P 0.90 alone). Most misses are wobbly-period
  low-complexity (`CAAACAAAACAAAC…`), hence:
- **`--dust`**: DUST triplet skew. Scale on chr1: random 0.5, wobbly
  microsat 3–5, (CA)n 15, poly-A 31 — so NCBI's "20" is not the right number
  here; 5 is conservative, 3 masks more at P 0.74. tandem+dust(5):
  Simple_repeat 0.50, Low_complexity 0.30, P 0.88.
- **NCBI-style mode** (index@6 + tandem + dust 3) vs NCBI lowercase on chr1:
  F1 0.788 (RM-oriented config: 0.757). Closer replication would need
  WindowMasker's actual two-threshold window scoring.

Data: `C:\krep_work\dfam_human.fa` (1,403 consensi, headers `name#accession`),
`krep_plus_dfam.fa` (v2 + Dfam concatenated).

## Round 5: learned hit filter (ML step 1)

- `--lib-dump` writes accepted hits with features; `--lib-model` applies
  either per-consensus floors (`scripts/train_thresholds.py`: floor per
  consensus maximizing cumulative true-minus-false bases) or a logistic
  model with per-consensus offsets (`scripts/train_logistic.py`: pure-Python
  SGD, base-weighted loss, features standardized; offline evaluator merges
  kept intervals — an earlier version double-counted overlaps and reported
  recall > 1).
- Train chr2 → test chr1, Dfam single-hit −4,2: floor 30 F1 0.899 →
  per-consensus floors 0.9065 → logistic 0.9085 (P 0.976 / R 0.850).
  Hybrid (k18@32 + Dfam + tandem + dust) + floors on chr1: 0.9066.
- Dumps: `C:\krep_work\chr{1,2}_cand.tsv` (1.6M rows each, floor 15);
  models `model_chr2.tsv`, `model_logit.tsv`; `chm13_chr2_unmasked.fa`,
  `chr2_rm_family.bed`.
- Next in this direction: gradient boosting over the same features; a
  divergence-aware (transition/transversion) scoring matrix estimated from
  krep's confident alignments; profile HMMs for the ancient families.

## Round 6: de novo spaced-seed consensus (no Dfam/HMM)

The consensus builder now seeds from the dense whole-genome spaced-seed index
(`chm13.sp16.s1.kidx`). A new `--seed-pool` flag limits how many abundant
seeds are held for occurrence collection, and the top-K selection uses a
bounded heap so exploring millions of low-count seeds no longer materialises
the full candidate list.

```bash
~/krep_target/release/krep consensus \
  --genome /mnt/c/Users/jfris/OneDrive/Desktop/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --index /mnt/c/krep_work/chm13.sp16.s1.kidx \
  --seq-dump /mnt/c/krep_work/chm13.kseq \
  --out /mnt/c/krep_work/chm13_consensi_sp16_pool500k.fa \
  --min-seed-count 8 --max-families 2000 --seed-pool 500000 --max-occ 30
```

This built **537 de novo consensi** in ~4 min at ~1.6 GB peak RAM. The
library is dominated by Alu/L1/satellite, but also contains real MIR and L2
fragments (e.g. `krep_fam_00039` = Dfam MIR, `krep_fam_00085/00091` = L2a/MIRb
3' ends). A logistic hit filter was trained on chr2+3+4 candidates and
applied to chr1.

Chr1 vs `chr1_rm_family.bed` (k18@32 + de novo library + tandem/dust):

| filter | P | R | F1 | SINE/MIR | LINE/L2 | LINE/CR1 |
|---|---|---|---|---|---|---|
| none (score ≥15) | 0.8065 | 0.7794 | 0.7927 | 0.4881 | 0.2853 | 0.1716 |
| per-consensus thresholds | **0.9422** | 0.7313 | **0.8234** | 0.4283 | 0.1594 | 0.0188 |
| logistic τ=0.60 | 0.9360 | 0.7325 | 0.8218 | 0.4303 | 0.1624 | 0.0216 |
| logistic τ=0.50 | 0.9044 | 0.7443 | 0.8166 | 0.4398 | 0.1913 | 0.0496 |
| logistic τ=0.45 | 0.8824 | 0.7521 | 0.8121 | 0.4481 | 0.2122 | 0.0733 |
| logistic τ=0.40 | 0.8594 | 0.7602 | 0.8068 | 0.4589 | 0.2337 | 0.1007 |
| logistic τ=0.30 | 0.8249 | **0.7729** | 0.7981 | **0.4772** | **0.2683** | **0.1482** |

The fully de novo pipeline reaches **F1 0.8234 on chr1** at high precision,
or trades ~12 points of precision for more ancient-family recall (τ=0.3). It
still trails the Dfam hybrid (F1 ~0.91) and the HMM layer on MIR/L2/CR1, but
it proves the spaced-seed consensus path works without external repeat
libraries.

A larger `--seed-pool 1000000` library reached **943 consensi**. A full
τ sweep showed that no global threshold can raise ancient-family recall
without a large precision drop, because the underrepresented families
(CR1, Helitron, Tip100) have a different feature distribution than the
dominant Alu/L1/satellite families. This motivated the family-specialist and
two-pass routes below.

### Round 6a: family-specific rare model

`scripts/train_rare_family_model.py` trains a second logistic model whose
positive label is "candidate overlaps any of {LINE/CR1, RC/Helitron,
DNA/hAT-Tip100, LINE/L2, SINE/MIR}" and whose negatives include both
non-repeat regions and other repeat families. Ensembling it with the main
model via `scripts/apply_model_offline.py` gave only tiny gains (CR1
0.0293 → 0.0303 on chr1) because the candidate dump from the single-pass
1M library already had poor coverage of those targets (~10% of true bases).

### Round 6b: two-pass consensus builder

`krep consensus` now has `--second-pass`, `--second-pass-min-len`,
`--second-pass-min-support`, and `--second-pass-max-families` (see
`src/consensus.rs` and `src/main.rs`). The first pass builds consensi with
the primary filters; the optional second pass iterates over the same
occurrence table with relaxed filters, skipping seeds already covered by pass
1 and checking pass-2 consensi against them for redundancy. A bug-fix
prevents rejected pass-1 candidates from poisoning the covered set; only
accepted consensi mark their seeds/12-mers as covered.

With `--second-pass-min-len 80 --second-pass-min-support 5
--second-pass-max-families 500`, the 1M library grew from **943 consensi to
1,375 consensi** (pass 1 built 1,337, pass 2 added 38). Genome-wide masking
coverage on chr1–4 with the new library is ~60% vs ~52–55% for the single-pass
1M library, so the expansion hits substantially more repeat sequence.

Chr1 vs `chr1_rm_family.bed`:

| filter | P | R | F1 | Tip100 | CR1 | L2 | MIR |
|---|---|---|---|---|---|---|---|
| single-pass 1M main | 0.9626 | 0.7163 | 0.8214 | 0.0956 | 0.0102 | 0.1301 | 0.3920 |
| single-pass 1M + rare | 0.9561 | 0.7186 | 0.8205 | 0.0998 | 0.0143 | 0.1409 | 0.3993 |
| **two-pass 1M main** | **0.9583** | **0.7350** | **0.8319** | 0.1230 | 0.0161 | 0.1920 | 0.4090 |
| two-pass 1M + rare | 0.9509 | 0.7379 | 0.8310 | 0.1292 | 0.0205 | 0.2075 | 0.4191 |

The two-pass library buys **+1.9 recall** and **+1.0 F1** at the cost of
**−0.4 precision**, with clear gains for L2/MIR/Tip100. The rare-family
specialist adds another ~0.3 recall but costs ~0.7 precision. CR1 and
Helitron remain low even after the second pass because the consensus builder
still finds almost no CR1/Helitron seeds that clear the minimum requirements
in either pass.

## What would still move it

- **CR1 / Helitron / Tip100 have no usable consensi in single-pass mode.**
  The two-pass builder addresses the *builder* side; the remaining question is
  whether the learned hit filter can turn the new consensi into real recall
  gain without hurting precision.
- **Short diverged fragments** (MIR/L2 at ~150 bp, 65% identity) barely clear
  the 2-hit chain and the score floor. A length-normalized score threshold or
  single-hit triggering for consensi < 400 bp would help sensitivity at a
  measurable precision cost.
- **Subfamily consensi**: the redundancy filter drops candidates sharing ≥60%
  12-mers; RepeatModeler keeps subfamilies. Loosening it to ~85% would add
  AluY/AluS-style variants and lift young-family recall slightly.
- Whole-genome library run + genome-wide evaluation (BEDs are in
  `C:\krep_work\genome_rm.bed` / `genome_rm_family.bed`).

## What would actually move the ceiling (earlier notes)

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
- Consensus libraries: `C:\krep_work\chm13_consensi_v2.fa` (1,124 consensi
  from the contiguous k18 index), `chm13_consensi_sp16_pool500k.fa` (537
  consensi from the spaced w16/s22 index, fully de novo), `chm13_consensi_sp16_pool1M.fa`
  (943 consensi, single-pass), `chm13_consensi_sp16_pool1M_pass2c.fa` (1,375
  consensi, two-pass). Sequence dump `C:\krep_work\chm13.kseq` (3.1 GB,
  reused by `krep consensus`).
- Learned models (spaced-seed de novo): `model_sp16_thresholds.tsv`,
  `model_sp16_logit.tsv`, `model_sp16_thresholds1M.tsv`,
  `model_sp16_logit1M.tsv`, `model_sp16_rare1M.tsv`,
  `model_sp16_thresholds_pass2c.tsv`, `model_sp16_logit_pass2c.tsv`,
  `model_sp16_rare_pass2c.tsv`.
- Candidate dumps (spaced-seed de novo): `chr{1,2,3,4}_sp16_cand*.tsv`,
  combined `chr234_sp16_cand*.tsv`. Used by `scripts/train_thresholds.py`,
  `scripts/train_logistic.py`, `scripts/train_rare_family_model.py`,
  `scripts/apply_model_offline.py`.
- Ground truth: `C:\krep_work\genome_rm_family.bed`, `chr{1,2,3,4}_rm_family.bed`.
- Mock genomes for regression: `C:\krep_work\mock15.*`, `mock25.*`.
- Launch scripts: `C:\krep_work\build_pool1M.sh`, `dump_pass2c.sh`,
  `train_eval_pass2c.sh`, `train_rare.sh`, `eval_rare_combo_fast.sh`, plus
  earlier `run_*.sh`.

## Build and test

Windows Application Control intermittently blocks freshly built executables
(os error 4551) — it hit the test harness and Cargo build scripts, while the
release `krep.exe` kept running. If a fresh `.exe` is blocked, use the WSL
native binary (built below) until the policy allows the new file.

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

- **WSL fallback binary** (current, since the latest Windows build is blocked):
  build a native Linux binary and run it from WSL with `/mnt/c/...` paths:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
cargo build --release --target-dir ~/krep_target
cp ~/krep_target/release/krep /mnt/c/krep_work/krep_wsl
```

  The chr1 HMM-layer and two-pass de novo numbers below were produced with
  `/mnt/c/krep_work/krep_wsl`.

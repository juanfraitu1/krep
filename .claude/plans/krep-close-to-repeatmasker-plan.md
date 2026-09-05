# Plan: make krep much closer to RepeatMasker on CHM13 chr1

## Goal
Push the base-level agreement between krep's de novo soft mask and the T2T-CHM13 RepeatMasker soft mask as high as possible on the first 10 Mb of chr1, while making segmental-duplication simulation an optional `krep mock` mode.

## Current baseline (chr1 10 Mb slice, k=18, threshold=5, graph-gap=100, mismatch1)

```text
Precision: 0.8821
Recall:    0.6890
F1:        0.7737
```

The main issue is recall: krep misses ~31% of RepeatMasker-masked bases.

## Steps

### 1. Make SEG_DUP optional in mock mode
- Add a `--segdup` boolean flag to `Commands::Mock`.
- Only include `SEG_DUP` copies when the flag is set.
- Default behaviour: no SEG_DUP (cleaner mock for advisor-style comparisons).
- Update README mock table and examples.

### 2. Add `krep tune-mask` subcommand
- New subcommand that runs a grid search over masking parameters on a real genome slice and reports the best F1 against a reference soft-masked FASTA.
- Parameters to sweep:
  - `--k-list` candidates: 15, 17, 18, 19, 21, 25, 31, and small unions like `17,21`, `15,21,31`.
  - `--threshold`: 2, 3, 4, 5.
  - `--graph-gap`: 100, 250, 500, 1000.
  - `--mismatch1`: on/off.
  - `--min-len`: 30, 50, 100.
- It will internally call the mask logic and `compare_mask` logic without writing intermediate files, or write them to temp and clean up.
- Output: ranked parameter combos with precision/recall/F1 and runtime.

### 3. Run the tuning on the chr1 10 Mb slice
- Use the original soft-masked slice as reference.
- Identify the best single k, best gap, best threshold, and whether mismatch1 helps.
- Try k-list unions that may catch both young (large k) and old/diverged (small k) repeats.

### 4. Adopt and document the winning real-genome command
- Update README quick-start with the tuned command.
- Add a "Real genome (CHM13)" section showing precision/recall/F1.

### 5. Memory note for full chr1
- Full chr1 (~248 Mb) cannot allocate the default CBF on an 8 GB laptop.
- Keep `--cbf-factor` tunable; document that 16 GB is recommended for whole-genome runs and that the tool was tuned on a 10 Mb slice.
- Optionally run a 50 Mb slice if memory allows after tuning.

## Expected outcome
- SEG_DUP becomes optional.
- A `tune-mask` command exists for reusable parameter search.
- Real-genome F1 improves meaningfully from 0.7737, ideally into the high 0.8s or low 0.9s.
- README reflects the new best real-genome parameters.

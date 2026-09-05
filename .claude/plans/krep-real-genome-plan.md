# Plan: run krep on the real CHM13 genome and add de-mask / GTF support

## Goal
1. Add a `demask` subcommand that converts a soft-masked FASTA back to an all-uppercase FASTA.
2. Add `--out-format bed|gtf` to the `mask` subcommand so masked regions can be emitted as BED or GTF.
3. Run krep on the real `../GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` (3 GB soft-masked T2T-CHM13) and compare the re-soft-masked output with the original to see how well the de novo signal reproduces RepeatMasker soft masking.

## Why these features are useful
- `demask`: many genomes are distributed soft-masked, but some downstream tools need uppercase sequence. A dedicated command makes this discoverable.
- GTF: annotation pipelines often expect GTF/GFF rather than BED.
- Real-genome test: validates that krep scales beyond the synthetic 10 Mb benchmark and gives a concrete overlap metric against a widely used reference masking.

## Implementation steps

### 1. Add `demask` subcommand
- Add `Commands::Demask { genome: PathBuf, out: PathBuf }` to the CLI enum.
- Add `run_demask` that calls `fasta::read` (already uppercases) and `fasta::write`.
- Wire it into `main()`.
- Update README.

### 2. Add `--out-format` to `mask`
- Add `#[derive(ValueEnum, Clone)] enum OutputFormat { Bed, Gtf }`.
- Add `--out-format` argument to `Commands::Mask` (default `bed`).
- In `run_mask`, choose between:
  - `mask::regions_to_bed(&regions)`
  - `mask::regions_to_gtf(&regions)` which writes 1-based inclusive coordinates with a simple attribute line (`repeat_id "{i}"`).
- Update README.

### 3. Test on CHM13
- Run `krep demask` on the 3 GB input to produce `chm13_unmasked.fa`.
- Run `krep mask` on the unmasked FASTA using the current best mock-genome parameters (`--graph --graph-gap 100 --k 18 --threshold 5 --mismatch1`) and output both BED and soft-masked FASTA.
- Compare the krep soft-masked output to the original soft-masked file base-by-base, reporting:
  - total bases compared
  - agreement on masked (both lowercase)
  - agreement on unmasked (both uppercase ACGT)
  - disagreement rate
  - note: the files will not be identical because krep is de novo, but we expect high agreement in repeat-dense regions.

## Risks and mitigations
- **Memory**: the CHM13 FASTA is ~3 GB and the largest chromosome CBF will be ~2 GB. Peak RAM could reach ~5–6 GB. This fits a 16 GB laptop but may be tight on 8 GB. If it fails, we will fall back to masking a single large chromosome (e.g., chr1) for the proof-of-concept.
- **Runtime**: linear extrapolation from the 10 Mb benchmark suggests ~30–60 minutes for the full genome. We will run it in the background and report progress.
- **Soft-masked input**: `fasta::read` already uppercases sequence, so k-mer counting will work on the original file. The explicit `demask` step is still useful as a standalone command and as a clean intermediate file for the comparison.

## Expected outcome
- A new `demask` subcommand and `--out-format gtf` option in the CLI.
- A BED/GTF of krep-predicted repeats for CHM13.
- A quantitative base-level comparison between krep's soft masking and the original CHM13 soft masking.

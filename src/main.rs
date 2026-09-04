mod bitvec;
mod cbf;
mod evaluate;
mod fasta;
mod index;
mod kmer;
mod mask;
mod mock;
mod rng;

use clap::{Parser, Subcommand, ValueEnum};
use std::fs;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "krep")]
#[command(about = "Lightweight de novo repeat masking with a counting Bloom filter")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate a mock genome with known repetitive elements.
    Mock {
        /// Genome size (supports suffixes: k, m, g). Default 10m.
        #[arg(short, long, default_value = "10m")]
        size: String,

        /// Background GC fraction (0.0–1.0).
        #[arg(short, long, default_value_t = 0.5)]
        gc: f64,

        /// Random seed.
        #[arg(short, long, default_value_t = 42)]
        seed: u64,

        /// Per-base divergence for repeat copies (0.0–1.0).
        #[arg(short, long, default_value_t = 0.05)]
        divergence: f64,

        /// Per-base divergence for segmental-duplication copies (0.0–1.0).
        #[arg(long, default_value_t = 0.02)]
        segdup_divergence: f64,

        /// Include a segmental-duplication family in the mock genome.
        #[arg(long, default_value_t = false)]
        segdup: bool,

        /// Output FASTA file.
        #[arg(short, long, default_value = "genome.fa")]
        out: PathBuf,

        /// Output ground-truth BED file.
        #[arg(short, long, default_value = "repeats.bed")]
        bed: PathBuf,
    },

    /// Mask repetitive regions de novo using a counting Bloom filter.
    Mask {
        /// Input FASTA file.
        #[arg(short, long)]
        genome: PathBuf,

        /// k-mer length.
        #[arg(short, long, default_value_t = 21, group = "kmer")]
        k: usize,

        /// Comma-separated list of k values. The union of masked regions from each
        /// k is reported. Useful for catching both old/diverged repeats (small k)
        /// and recent precise repeats (large k).
        #[arg(long, value_delimiter = ',', group = "kmer")]
        k_list: Option<Vec<usize>>,

        /// Count threshold: a k-mer is considered repetitive when its CBF count is ≥ this.
        #[arg(short, long, default_value_t = 3)]
        threshold: u8,

        /// Minimum masked region length.
        #[arg(short = 'l', long, default_value_t = 30)]
        min_len: usize,

        /// Sequence window length used for density scoring.
        #[arg(short = 'w', long, default_value_t = 40)]
        window: usize,

        /// Minimum fraction of high-count k-mers required in a window.
        #[arg(short = 'd', long, default_value_t = 0.25)]
        density: f64,

        /// Output file for masked regions.
        #[arg(short, long, default_value = "masked.bed")]
        out: PathBuf,

        /// Output format: BED or GTF.
        #[arg(long, value_enum, default_value_t = OutputFormat::Bed)]
        out_format: OutputFormat,

        /// Output soft-masked FASTA (lowercase for masked bases).
        #[arg(long, value_name = "FILE")]
        soft: Option<PathBuf>,

        /// Output hard-masked FASTA (N for masked bases).
        #[arg(long, value_name = "FILE")]
        hard: Option<PathBuf>,

        /// Use graph-based connected-component masking instead of sliding-window
        /// density. High-count k-mer positions are linked if they are within
        /// `graph-gap` bases; connected components spanning ≥ min_len are masked.
        #[arg(long, default_value_t = false)]
        graph: bool,

        /// Maximum gap (in bases) between high-count k-mers for them to be linked
        /// in graph mode.
        #[arg(long, default_value_t = 100)]
        graph_gap: usize,

        /// Allow a k-mer to be considered repetitive if any Hamming-distance-1
        /// neighbor is high-count. This helps recover diverged repeat copies at
        /// the cost of more CBF queries.
        #[arg(long, default_value_t = false)]
        mismatch1: bool,

        /// Use a de Bruijn graph of consecutive high-count k-mers to build repeat
        /// contigs. Only components whose total CBF abundance passes the
        /// `assembly-abundance` filter are masked. Useful for very noisy or
        /// ancient repeat copies.
        #[arg(long, default_value_t = false)]
        assembly: bool,

        /// Minimum total CBF abundance for a de Bruijn graph component to be
        /// kept. 0 means automatic (threshold * 5, at least 20).
        #[arg(long, value_name = "COUNT", default_value_t = 0)]
        assembly_abundance: u64,

        /// Use a prebuilt genome-wide index (see `krep index`) instead of
        /// counting k-mers within the input FASTA. `k` comes from the index;
        /// `--threshold` is then a genome-wide occurrence count, which is a
        /// completely different scale from the per-slice counts, so it must be
        /// retuned rather than carried over.
        #[arg(long, value_name = "FILE")]
        index: Option<PathBuf>,

        /// Genome-wide occurrence threshold used with `--index`.
        #[arg(long, default_value_t = 10)]
        index_threshold: u32,

        /// Counting Bloom filter size multiplier: number of CBF slots is
        /// `factor × number_of_kmers`. Lower values reduce memory but increase
        /// hash collisions. Default 8; try 4 on memory-constrained laptops.
        #[arg(long, value_name = "FACTOR", default_value_t = 8)]
        cbf_factor: usize,
    },

    /// Build a genome-wide k-mer count index.
    ///
    /// Run this once over the full genome, then mask any subset against it.
    /// Counting genome-wide is the whole point: a repeat family with a few
    /// hundred copies genome-wide has an expected count below 1 inside a 10 Mb
    /// window, so per-slice counting cannot see it at all.
    Index {
        /// Input FASTA (the whole genome).
        #[arg(short, long)]
        genome: PathBuf,

        /// k-mer length. At genome scale k must be large enough that a random
        /// k-mer is not expected to recur by chance (k >= 17 for 3.1 Gb).
        #[arg(short, long, default_value_t = 18)]
        k: usize,

        /// Keep 1 in SAMPLE k-mers, chosen by hash. Must be a power of two.
        /// Sampling shrinks the table without biasing counts: a sampled k-mer
        /// is counted at every occurrence, so stored counts stay exact.
        #[arg(long, default_value_t = 16)]
        sample: u64,

        /// Drop k-mers occurring fewer than this many times genome-wide.
        /// Singletons dominate the table; discarding them is what makes the
        /// index loadable on a small machine.
        #[arg(long, default_value_t = 2)]
        min_count: u32,

        /// K-mers buffered in RAM before spilling a sorted run to disk.
        /// Each buffered k-mer costs 8 bytes, so 48M ~ 384 MB.
        #[arg(long, default_value_t = 48_000_000)]
        buffer: usize,

        /// Directory for temporary sorted runs. Put this on a fast local disk;
        /// under WSL, the ext4 root is far faster than /mnt/c.
        #[arg(long, default_value = "krep_tmp")]
        tmp_dir: PathBuf,

        /// Output index file.
        #[arg(short, long, default_value = "genome.kidx")]
        out: PathBuf,

        /// Print per-record progress while counting.
        #[arg(long, default_value_t = false)]
        verbose: bool,
    },

    /// Convert a soft-masked FASTA to an all-uppercase FASTA.
    Demask {
        /// Input soft-masked FASTA.
        #[arg(short, long)]
        genome: PathBuf,

        /// Output uppercase FASTA.
        #[arg(short, long, default_value = "demasked.fa")]
        out: PathBuf,
    },

    /// Compare two FASTA soft masks base-by-base (e.g. RepeatMasker vs krep).
    CompareMask {
        /// Reference / original soft-masked FASTA.
        #[arg(short, long)]
        reference: PathBuf,

        /// Predicted soft-masked FASTA.
        #[arg(short, long)]
        predicted: PathBuf,
    },

    /// Grid-search masking parameters against a reference soft mask and print
    /// the best combinations.
    TuneMask {
        /// Input uppercase FASTA to mask.
        #[arg(short, long)]
        genome: PathBuf,

        /// Reference soft-masked FASTA to compare against.
        #[arg(short, long)]
        reference: PathBuf,

        /// Comma-separated k values to try (e.g. 17,21). Single k values are tried
        /// individually and as combined k-lists.
        #[arg(long, value_delimiter = ',', default_value = "17,19,21,25")]
        ks: Vec<usize>,

        /// Comma-separated thresholds to try.
        #[arg(long, value_delimiter = ',', default_value = "2,3,4,5")]
        thresholds: Vec<u8>,

        /// Comma-separated graph-gap values to try. Use 0 to disable graph mode
        /// (window-density mode is used instead).
        #[arg(long, value_delimiter = ',', default_value = "100,250,500")]
        graph_gaps: Vec<usize>,

        /// Try with and without --mismatch1.
        #[arg(long, default_value_t = true)]
        try_mismatch: bool,

        /// Minimum masked region length.
        #[arg(long, default_value_t = 30)]
        min_len: usize,

        /// Window length for window-density mode.
        #[arg(long, default_value_t = 40)]
        window: usize,

        /// Density threshold for window-density mode.
        #[arg(long, default_value_t = 0.25)]
        density: f64,

        /// CBF size multiplier.
        #[arg(long, default_value_t = 8)]
        cbf_factor: usize,

        /// Print at most this many top results.
        #[arg(long, default_value_t = 10)]
        top: usize,
    },

    /// Evaluate predicted BED against a ground-truth BED.
    Evaluate {
        /// Ground-truth BED (4-column supported for per-family recall).
        #[arg(short, long)]
        truth: PathBuf,

        /// Predicted BED.
        #[arg(short, long)]
        pred: PathBuf,
    },
}

fn parse_size(s: &str) -> Result<usize, String> {
    let s = s.trim().to_ascii_lowercase();
    let (num_part, mult) = if s.ends_with('k') {
        (&s[..s.len() - 1], 1_000usize)
    } else if s.ends_with('m') {
        (&s[..s.len() - 1], 1_000_000usize)
    } else if s.ends_with('g') {
        (&s[..s.len() - 1], 1_000_000_000usize)
    } else {
        (&s[..], 1usize)
    };
    let n: f64 = num_part
        .parse()
        .map_err(|_| format!("invalid size: {}", s))?;
    let size = (n * mult as f64) as usize;
    if size == 0 {
        return Err("genome size must be > 0".into());
    }
    Ok(size)
}

fn run_mock(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    let Commands::Mock {
        size,
        gc,
        seed,
        divergence,
        segdup_divergence,
        segdup,
        out,
        bed,
    } = args
    else {
        unreachable!()
    };

    if !(0.0..=1.0).contains(gc) {
        return Err("GC must be between 0.0 and 1.0".into());
    }
    if !(0.0..=1.0).contains(divergence) {
        return Err("divergence must be between 0.0 and 1.0".into());
    }
    if !(0.0..=1.0).contains(segdup_divergence) {
        return Err("segdup_divergence must be between 0.0 and 1.0".into());
    }

    let size = parse_size(size)?;
    let mut rng = rng::Rng::new(*seed);
    let families = if *segdup {
        mock::families_with_segdup()
    } else {
        mock::default_families()
    };
    let (genome_seq, repeats) = mock::build_mock_genome_with_segdup(
        size,
        *gc,
        &families,
        *divergence,
        *segdup_divergence,
        &mut rng,
    );

    fasta::write(&[fasta::Record::new("chr1", genome_seq)], out)?;
    fs::write(bed, mock::repeats_to_bed("chr1", &repeats))?;

    let repeat_bp: usize = repeats.iter().map(|r| r.len()).sum();
    println!(
        "Generated chr1 ({} bp) with {} repeat copies covering {} bp ({:.2}%).",
        size,
        repeats.len(),
        repeat_bp,
        100.0 * repeat_bp as f64 / size as f64
    );
    Ok(())
}

fn run_mask(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    let Commands::Mask {
        genome,
        k,
        k_list,
        threshold,
        min_len,
        window,
        density,
        out,
        out_format,
        soft,
        hard,
        graph,
        graph_gap,
        mismatch1,
        assembly,
        assembly_abundance,
        cbf_factor,
        index: index_path,
        index_threshold,
    } = args
    else {
        unreachable!()
    };

    // Index-driven path: stream the target, look up genome-wide counts.
    if let Some(idx_path) = index_path {
        return run_mask_indexed(
            genome,
            idx_path,
            *index_threshold,
            *graph_gap,
            *min_len,
            out,
            out_format,
            soft.as_deref(),
        );
    }

    let records = fasta::read_uppercase(genome)?;
    let seq_records: Vec<(String, Vec<u8>)> =
        records.into_iter().map(|r| (r.header, r.seq)).collect();

    let ks: Vec<usize> = k_list
        .clone()
        .unwrap_or_else(|| vec![*k])
        .iter()
        .copied()
        .filter(|x| *x > 0)
        .collect();
    if ks.is_empty() {
        return Err("at least one positive k value is required".into());
    }

    let mismatches = if *mismatch1 { 1 } else { 0 };
    let cbf_factor = (*cbf_factor).max(1);
    let regions = if *assembly {
        mask::mask_fasta_assembly(
            &seq_records,
            &ks,
            *threshold,
            *min_len,
            mismatches,
            *assembly_abundance,
            cbf_factor,
        )
    } else if *graph {
        mask::mask_fasta_graph(
            &seq_records, &ks, *threshold, *graph_gap, *min_len, mismatches, cbf_factor)
    } else {
        mask::mask_fasta_union(
            &seq_records, &ks, *threshold, *window, *density, *min_len, mismatches, cbf_factor)
    };
    let output_text = match out_format {
        OutputFormat::Bed => mask::regions_to_bed(&regions),
        OutputFormat::Gtf => mask::regions_to_gtf(&regions),
    };
    fs::write(out, output_text)?;

    // Optional masked FASTA outputs.
    if let Some(path) = soft {
        let soft_records = apply_mask(&seq_records, &regions, MaskMode::Soft);
        fasta::write(&soft_records, path)?;
    }
    if let Some(path) = hard {
        let hard_records = apply_mask(&seq_records, &regions, MaskMode::Hard);
        fasta::write(&hard_records, path)?;
    }

    let masked_bp: usize = regions.iter().map(|r| r.len()).sum();
    println!("Masked {} regions ({} bp).", regions.len(), masked_bp);
    Ok(())
}

enum MaskMode {
    Soft,
    Hard,
}

fn apply_mask(
    records: &[(String, Vec<u8>)],
    regions: &[mask::MaskedRegion],
    mode: MaskMode,
) -> Vec<fasta::Record> {
    let mut out_records = Vec::with_capacity(records.len());
    for (header, seq) in records {
        let mut new_seq = seq.clone();
        for r in regions.iter().filter(|r| r.chrom == fasta::seq_id(header)) {
            for i in r.start..r.end.min(new_seq.len()) {
                new_seq[i] = match mode {
                    MaskMode::Soft => new_seq[i].to_ascii_lowercase(),
                    MaskMode::Hard => b'N',
                };
            }
        }
        out_records.push(fasta::Record::new(header.clone(), new_seq));
    }
    out_records
}

fn run_demask(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    let Commands::Demask { genome, out } = args else {
        unreachable!()
    };
    let mut records = fasta::read(genome)?;
    for rec in &mut records {
        for b in &mut rec.seq {
            *b = b.to_ascii_uppercase();
        }
    }
    fasta::write(&records, out)?;
    println!("Demasked {} sequence records to {}.", records.len(), out.display());
    Ok(())
}

#[derive(ValueEnum, Clone, Debug)]
enum OutputFormat {
    Bed,
    Gtf,
}

fn run_compare_mask(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    let Commands::CompareMask { reference, predicted } = args else {
        unreachable!()
    };

    // Streamed one record at a time. Loading both FASTAs in full costs ~6 GB at
    // genome scale, which does not fit on a typical laptop.
    let mut ref_stream = fasta::FastaStream::open(reference)?;
    let mut pred_stream = fasta::FastaStream::open(predicted)?;

    let mut total = 0u64;
    let mut both_masked = 0u64;
    let mut both_unmasked = 0u64;
    let mut ref_masked_only = 0u64;
    let mut pred_masked_only = 0u64;
    let mut compared_records = 0usize;

    let mut pending_ref = ref_stream.next_record(false)?;

    println!(
        "{:<18} {:>13} {:>9} {:>9} {:>8} {:>8} {:>8}",
        "record", "bases", "ref%", "pred%", "prec", "recall", "F1"
    );

    while let Some(pred_rec) = pred_stream.next_record(false)? {
        // Allow the predicted file to cover a subset of the reference records
        // (e.g. chr1 only, scored against the whole-genome FASTA).
        loop {
            match &pending_ref {
                Some(r) if fasta::seq_id(&r.header) == fasta::seq_id(&pred_rec.header) => break,
                Some(_) => pending_ref = ref_stream.next_record(false)?,
                None => {
                    return Err(format!(
                        "record {} not found in reference",
                        fasta::seq_id(&pred_rec.header)
                    )
                    .into())
                }
            }
        }
        let ref_rec = pending_ref.take().expect("matched above");
        pending_ref = ref_stream.next_record(false)?;
        compared_records += 1;

        let len = ref_rec.seq.len().min(pred_rec.seq.len());
        let (mut tp, mut fp, mut fn_, mut tn) = (0u64, 0u64, 0u64, 0u64);
        for i in 0..len {
            // Treat non-ACGT or lowercase as masked in the reference.
            let ref_masked = !matches!(ref_rec.seq[i], b'A' | b'C' | b'G' | b'T');
            let pred_masked = pred_rec.seq[i].is_ascii_lowercase();
            match (ref_masked, pred_masked) {
                (true, true) => tp += 1,
                (false, true) => fp += 1,
                (true, false) => fn_ += 1,
                (false, false) => tn += 1,
            }
        }

        let (p, r, f) = prf(tp, fp, fn_);
        println!(
            "{:<18} {:>13} {:>8.2}% {:>8.2}% {:>8.4} {:>8.4} {:>8.4}",
            fasta::seq_id(&ref_rec.header),
            len,
            100.0 * (tp + fn_) as f64 / len.max(1) as f64,
            100.0 * (tp + fp) as f64 / len.max(1) as f64,
            p,
            r,
            f
        );

        total += len as u64;
        both_masked += tp;
        pred_masked_only += fp;
        ref_masked_only += fn_;
        both_unmasked += tn;
    }

    if compared_records == 0 {
        return Err("no records compared".into());
    }

    let (precision, recall, f1) = prf(both_masked, pred_masked_only, ref_masked_only);
    println!("\n--- totals over {} record(s) ---", compared_records);
    println!("Bases compared:        {}", total);
    println!(
        "Reference masked:      {} ({:.2}%)",
        both_masked + ref_masked_only,
        100.0 * (both_masked + ref_masked_only) as f64 / total as f64
    );
    println!(
        "Predicted masked:      {} ({:.2}%)",
        both_masked + pred_masked_only,
        100.0 * (both_masked + pred_masked_only) as f64 / total as f64
    );
    println!("Both masked (TP):      {}", both_masked);
    println!("Both unmasked (TN):    {}", both_unmasked);
    println!("Reference only (FN):   {}", ref_masked_only);
    println!("Predicted only (FP):   {}", pred_masked_only);
    println!("Precision:             {:.4}", precision);
    println!("Recall:                {:.4}", recall);
    println!("F1:                    {:.4}", f1);
    Ok(())
}

fn prf(tp: u64, fp: u64, fn_: u64) -> (f64, f64, f64) {
    let precision = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let recall = if tp + fn_ > 0 { tp as f64 / (tp + fn_) as f64 } else { 0.0 };
    let f1 = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
    } else {
        0.0
    };
    (precision, recall, f1)
}

fn run_tune_mask(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    use std::time::Instant;

    let Commands::TuneMask {
        genome,
        reference,
        ks,
        thresholds,
        graph_gaps,
        try_mismatch,
        min_len,
        window,
        density,
        cbf_factor,
        top,
    } = args
    else {
        unreachable!()
    };

    let ref_records = fasta::read(reference)?;
    let genome_records = fasta::read_uppercase(genome)?;
    if ref_records.len() != genome_records.len() {
        return Err("reference and genome FASTAs must have the same number of records".into());
    }

    let seq_records: Vec<(String, Vec<u8>)> = genome_records
        .into_iter()
        .map(|r| (r.header, r.seq))
        .collect();

    // Precompute the per-record reference mask once.
    let ref_masks: Vec<Vec<bool>> = ref_records
        .iter()
        .zip(seq_records.iter())
        .map(|(ref_rec, (_, seq))| {
            let len = ref_rec.seq.len().min(seq.len());
            (0..len)
                .map(|i| !matches!(ref_rec.seq[i], b'A' | b'C' | b'G' | b'T'))
                .collect()
        })
        .collect();

    let total_ref_masked: u64 = ref_masks.iter().map(|m| m.iter().filter(|b| **b).count() as u64).sum();

    // Build candidate k configurations: each single k, plus a few unions.
    let mut k_configs: Vec<Vec<usize>> = ks.iter().map(|k| vec![*k]).collect();
    if ks.len() >= 2 {
        k_configs.push(ks.iter().copied().collect());
    }
    if ks.len() >= 3 {
        for i in 0..ks.len() {
            for j in (i + 1)..ks.len() {
                k_configs.push(vec![ks[i], ks[j]]);
            }
        }
    }

    let mismatch_options: Vec<bool> = if *try_mismatch { vec![false, true] } else { vec![false] };

    #[derive(Debug)]
    struct Trial {
        ks: Vec<usize>,
        threshold: u8,
        graph_gap: usize,
        mismatch: bool,
        precision: f64,
        recall: f64,
        f1: f64,
        ms: u128,
    }

    let mut trials: Vec<Trial> = Vec::new();

    for k_cfg in &k_configs {
        for threshold in thresholds {
            for &graph_gap in graph_gaps {
                for mismatch in &mismatch_options {
                    let start = Instant::now();
                    let regions = if graph_gap == 0 {
                        mask::mask_fasta_union(
                            &seq_records,
                            k_cfg,
                            *threshold,
                            *window,
                            *density,
                            *min_len,
                            if *mismatch { 1 } else { 0 },
                            *cbf_factor,
                        )
                    } else {
                        mask::mask_fasta_graph(
                            &seq_records,
                            k_cfg,
                            *threshold,
                            graph_gap,
                            *min_len,
                            if *mismatch { 1 } else { 0 },
                            *cbf_factor,
                        )
                    };
                    let ms = start.elapsed().as_millis();

                    // Build predicted mask.
                    let mut pred_masks: Vec<Vec<bool>> =
                        seq_records.iter().map(|(_, seq)| vec![false; seq.len()]).collect();
                    for r in &regions {
                        if let Some(idx) = seq_records.iter().position(|(h, _)| h == &r.chrom) {
                            let seq_len = pred_masks[idx].len();
                            for p in r.start..r.end.min(seq_len) {
                                pred_masks[idx][p] = true;
                            }
                        }
                    }

                    let mut tp = 0u64;
                    let mut fp = 0u64;
                    let mut fn_ = 0u64;
                    for (rm, pm) in ref_masks.iter().zip(pred_masks.iter()) {
                        for (r, p) in rm.iter().zip(pm.iter()) {
                            match (*r, *p) {
                                (true, true) => tp += 1,
                                (false, true) => fp += 1,
                                (true, false) => fn_ += 1,
                                _ => {}
                            }
                        }
                    }

                    let precision = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
                    let recall = if total_ref_masked > 0 {
                        tp as f64 / total_ref_masked as f64
                    } else {
                        0.0
                    };
                    let f1 = if precision + recall > 0.0 {
                        2.0 * precision * recall / (precision + recall)
                    } else {
                        0.0
                    };

                    trials.push(Trial {
                        ks: k_cfg.clone(),
                        threshold: *threshold,
                        graph_gap,
                        mismatch: *mismatch,
                        precision,
                        recall,
                        f1,
                        ms,
                    });
                }
            }
        }
    }

    trials.sort_by(|a, b| b.f1.partial_cmp(&a.f1).unwrap());

    println!("Ranked parameter combinations by F1:");
    println!(
        "{:>3} {:>10} {:>3} {:>5} {:>5} {:>8} {:>8} {:>8} {:>8}",
        "#", "ks", "thr", "gap", "mis", "prec", "rec", "f1", "ms"
    );
    for (i, t) in trials.iter().take(*top).enumerate() {
        println!(
            "{:>3} {:>10} {:>3} {:>5} {:>5} {:>8.4} {:>8.4} {:>8.4} {:>8}",
            i + 1,
            t.ks.iter()
                .map(|k| k.to_string())
                .collect::<Vec<_>>()
                .join(","),
            t.threshold,
            t.graph_gap,
            if t.mismatch { "yes" } else { "no" },
            t.precision,
            t.recall,
            t.f1,
            t.ms
        );
    }

    Ok(())
}

fn run_evaluate(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    let Commands::Evaluate { truth, pred } = args else {
        unreachable!()
    };

    let truth_ivs = evaluate::read_bed(truth)?;
    let pred_ivs = evaluate::read_bed(pred)?;
    let metrics = evaluate::evaluate(&truth_ivs, &pred_ivs);
    print!("{}", evaluate::format_metrics(&metrics));
    Ok(())
}

fn run_index(args: &Commands) -> Result<(), Box<dyn std::error::Error>> {
    use std::time::Instant;

    let Commands::Index {
        genome,
        k,
        sample,
        min_count,
        buffer,
        tmp_dir,
        out,
        verbose,
    } = args
    else {
        unreachable!()
    };

    let t0 = Instant::now();
    println!(
        "Indexing {} with k={}, sample=1/{}, min_count={} ...",
        genome.display(),
        k,
        sample,
        min_count
    );
    let stats = index::build(
        genome, *k, *sample, *min_count, *buffer, tmp_dir, out, *verbose,
    )?;

    println!("\nRecords:            {}", stats.records);
    println!("Bases:              {}", stats.bases);
    println!("K-mers seen:        {}", stats.total_kmers);
    println!(
        "K-mers sampled:     {} (1/{:.1})",
        stats.sampled_kmers,
        stats.total_kmers as f64 / stats.sampled_kmers.max(1) as f64
    );
    println!("Sorted runs:        {}", stats.runs);
    println!(
        "Index entries:      {} ({:.0} MB resident)",
        stats.entries,
        stats.entries as f64 * 12.0 / 1e6
    );
    print!("{}", index::format_histogram(&stats));
    println!("\nWrote {} in {:.1}s", out.display(), t0.elapsed().as_secs_f64());
    Ok(())
}

fn run_mask_indexed(
    genome: &PathBuf,
    idx_path: &PathBuf,
    threshold: u32,
    max_gap: usize,
    min_len: usize,
    out: &PathBuf,
    out_format: &OutputFormat,
    soft: Option<&std::path::Path>,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::time::Instant;

    let t0 = Instant::now();
    let idx = index::KmerIndex::load(idx_path)?;
    println!(
        "Loaded index: k={}, sample=1/{}, {} entries ({:.0} MB), min_count={}",
        idx.k,
        idx.sample,
        idx.len(),
        idx.memory_bytes() as f64 / 1e6,
        idx.min_count
    );
    if threshold < idx.min_count {
        eprintln!(
            "warning: --index-threshold {} is below the index min_count {}; \
             k-mers rarer than {} were discarded at build time and read as count 0",
            threshold, idx.min_count, idx.min_count
        );
    }
    if max_gap < idx.sample as usize {
        eprintln!(
            "warning: --graph-gap {} is smaller than the sampling stride {}; \
             seeds will not link into regions",
            max_gap, idx.sample
        );
    }

    let mut stream = fasta::FastaStream::open(genome)?;
    let mut soft_writer = match soft {
        Some(p) => Some(fasta::FastaWriter::create(p)?),
        None => None,
    };

    let mut all_regions = Vec::new();
    let mut total_bases = 0u64;
    let mut total_masked = 0u64;

    while let Some(mut rec) = stream.next_record(true)? {
        let regions =
            mask::mask_sequence_indexed(&rec.header, &rec.seq, &idx, threshold, max_gap, min_len);
        let masked: usize = regions.iter().map(|r| r.end - r.start).sum();
        println!(
            "  {:<18} {:>12} bp  {:>8} regions  {:>6.2}% masked",
            fasta::seq_id(&rec.header),
            rec.seq.len(),
            regions.len(),
            100.0 * masked as f64 / rec.seq.len().max(1) as f64
        );
        total_bases += rec.seq.len() as u64;
        total_masked += masked as u64;

        if let Some(w) = soft_writer.as_mut() {
            mask::apply_soft_mask(&mut rec.seq, &regions);
            w.write_record(&rec.header, &rec.seq)?;
        }
        all_regions.extend(regions);
    }

    if let Some(w) = soft_writer {
        w.finish()?;
    }

    let text = match out_format {
        OutputFormat::Bed => mask::regions_to_bed(&all_regions),
        OutputFormat::Gtf => mask::regions_to_gtf(&all_regions),
    };
    fs::write(out, text)?;

    println!(
        "\nMasked {} regions, {} of {} bp ({:.2}%) in {:.1}s",
        all_regions.len(),
        total_masked,
        total_bases,
        100.0 * total_masked as f64 / total_bases.max(1) as f64,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() {
    let cli = Cli::parse();
    let result = match &cli.command {
        cmd @ Commands::Mock { .. } => run_mock(cmd),
        cmd @ Commands::Index { .. } => run_index(cmd),
        cmd @ Commands::Mask { .. } => run_mask(cmd),
        cmd @ Commands::Demask { .. } => run_demask(cmd),
        cmd @ Commands::CompareMask { .. } => run_compare_mask(cmd),
        cmd @ Commands::TuneMask { .. } => run_tune_mask(cmd),
        cmd @ Commands::Evaluate { .. } => run_evaluate(cmd),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

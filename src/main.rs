mod bitvec;
mod cbf;
mod evaluate;
mod fasta;
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

        /// Counting Bloom filter size multiplier: number of CBF slots is
        /// `factor × number_of_kmers`. Lower values reduce memory but increase
        /// hash collisions. Default 8; try 4 on memory-constrained laptops.
        #[arg(long, value_name = "FACTOR", default_value_t = 8)]
        cbf_factor: usize,
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
    } = args
    else {
        unreachable!()
    };

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
        for r in regions.iter().filter(|r| r.chrom == *header) {
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

    let ref_records = fasta::read(reference)?;
    let pred_records = fasta::read(predicted)?;

    if ref_records.len() != pred_records.len() {
        return Err("reference and predicted FASTAs must have the same number of records".into());
    }

    let mut total = 0u64;
    let mut both_masked = 0u64;
    let mut both_unmasked = 0u64;
    let mut ref_masked_only = 0u64;
    let mut pred_masked_only = 0u64;

    for (ref_rec, pred_rec) in ref_records.iter().zip(pred_records.iter()) {
        if ref_rec.header != pred_rec.header {
            return Err(format!(
                "header mismatch: {} vs {}",
                ref_rec.header, pred_rec.header
            )
            .into());
        }
        let len = ref_rec.seq.len().min(pred_rec.seq.len());
        for i in 0..len {
            let rb = ref_rec.seq[i];
            let pb = pred_rec.seq[i];
            // Treat non-ACGT or lowercase as masked in the reference.
            let ref_masked = !matches!(rb, b'A' | b'C' | b'G' | b'T');
            // krep only produces lowercase ACGT for soft-masked bases.
            let pred_masked = pb.is_ascii_lowercase();
            total += 1;
            match (ref_masked, pred_masked) {
                (true, true) => both_masked += 1,
                (false, false) => both_unmasked += 1,
                (true, false) => ref_masked_only += 1,
                (false, true) => pred_masked_only += 1,
            }
        }
    }

    let tp = both_masked;
    let fp = pred_masked_only;
    let fn_ = ref_masked_only;
    let precision = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let recall = if tp + fn_ > 0 { tp as f64 / (tp + fn_) as f64 } else { 0.0 };
    let f1 = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
    } else {
        0.0
    };

    println!("Bases compared:        {}", total);
    println!("Reference masked:      {} ({:.2}%)", both_masked + ref_masked_only, 100.0 * (both_masked + ref_masked_only) as f64 / total as f64);
    println!("Predicted masked:      {} ({:.2}%)", both_masked + pred_masked_only, 100.0 * (both_masked + pred_masked_only) as f64 / total as f64);
    println!("Both masked (TP):      {}", both_masked);
    println!("Both unmasked (TN):    {}", both_unmasked);
    println!("Reference only (FN):   {}", ref_masked_only);
    println!("Predicted only (FP):   {}", pred_masked_only);
    println!("Precision:             {:.4}", precision);
    println!("Recall:                {:.4}", recall);
    println!("F1:                    {:.4}", f1);
    Ok(())
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

fn main() {
    let cli = Cli::parse();
    let result = match &cli.command {
        cmd @ Commands::Mock { .. } => run_mock(cmd),
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

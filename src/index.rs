//! Genome-wide k-mer count index.
//!
//! The masker originally counted k-mers per FASTA record, so "how often does
//! this k-mer occur" was answered against whichever slice was passed in. A
//! repeat family with a few hundred copies genome-wide has an *expected* count
//! below 1 in a 10 Mb window, so it was invisible at that scope. This module
//! answers the question against the whole genome instead.
//!
//! Doing that exactly would need ~3.1e9 counters, which does not fit in the
//! few GB available. Two ideas make it fit:
//!
//! 1. **Hash sub-sampling (FracMinHash).** A k-mer is tracked only when
//!    `hash(kmer) % sample == 0`. Because the decision depends solely on the
//!    k-mer's own hash, a tracked k-mer is tracked at *every* one of its
//!    occurrences, so the stored count is its exact genomic occurrence count —
//!    unlike a Bloom filter, which inflates counts through collisions. With
//!    `sample = 16` the table shrinks 16x while seeds still land every ~16 bp
//!    inside a repeat, which is dense enough for gap-linking to reassemble.
//! 2. **External merge sort.** Sampled k-mers are buffered, sorted, run-length
//!    encoded and spilled to disk; the runs are then merged in a single
//!    streaming pass. Peak memory is the buffer, not the genome.

use crate::fasta::FastaStream;
use crate::kmer::KmerIter;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 8] = b"KREPIDX1";
const VERSION: u32 = 1;
const HEADER_LEN: usize = 56;

/// Salt for the sampling hash, so sampling is independent of the k-mer's own
/// bit pattern (raw 2-bit encodings are highly structured).
const SAMPLE_SALT: u64 = 0x9E37_79B9_7F4A_7C15;

/// Bits of k-mer prefix used for the lookup table that narrows binary search.
const PREFIX_BITS: u32 = 22;

/// Bases per parallel block during counting.
const BLOCK: usize = 32 << 20;

#[inline]
pub fn splitmix64(state: u64) -> u64 {
    let mut z = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// True when this canonical k-mer falls in the retained hash sample.
#[inline]
pub fn is_sampled(kmer: u64, sample_mask: u64) -> bool {
    splitmix64(kmer ^ SAMPLE_SALT) & sample_mask == 0
}

#[derive(Debug)]
pub struct BuildStats {
    pub records: usize,
    pub bases: u64,
    pub total_kmers: u64,
    pub sampled_kmers: u64,
    pub entries: u64,
    pub runs: usize,
    /// `histogram[b]` counts distinct sampled k-mers whose count is in
    /// `[2^b, 2^(b+1))`. Useful for picking a threshold at genome scope.
    pub histogram: [u64; 33],
}

// `[u64; 33]` is past the array size that gets a blanket `Default` impl.
impl Default for BuildStats {
    fn default() -> Self {
        Self {
            records: 0,
            bases: 0,
            total_kmers: 0,
            sampled_kmers: 0,
            entries: 0,
            runs: 0,
            histogram: [0u64; 33],
        }
    }
}

fn run_path(tmp_dir: &Path, idx: usize) -> PathBuf {
    tmp_dir.join(format!("krep_run_{:05}.bin", idx))
}

/// Sort, run-length encode and spill one buffer of sampled k-mers.
fn flush_run(buf: &mut Vec<u64>, tmp_dir: &Path, idx: usize) -> io::Result<()> {
    buf.par_sort_unstable();
    let mut w = BufWriter::with_capacity(1 << 20, File::create(run_path(tmp_dir, idx))?);
    let mut i = 0usize;
    while i < buf.len() {
        let kmer = buf[i];
        let mut j = i + 1;
        while j < buf.len() && buf[j] == kmer {
            j += 1;
        }
        w.write_all(&kmer.to_le_bytes())?;
        w.write_all(&((j - i) as u32).to_le_bytes())?;
        i = j;
    }
    w.flush()?;
    buf.clear();
    Ok(())
}

/// One sorted (kmer, count) run on disk, read as a cursor for k-way merging.
struct RunReader {
    reader: BufReader<File>,
    current: Option<(u64, u32)>,
}

impl RunReader {
    fn open(path: &Path) -> io::Result<Self> {
        let mut r = Self {
            reader: BufReader::with_capacity(1 << 20, File::open(path)?),
            current: None,
        };
        r.advance()?;
        Ok(r)
    }

    fn advance(&mut self) -> io::Result<()> {
        let mut buf = [0u8; 12];
        match self.reader.read_exact(&mut buf) {
            Ok(()) => {
                let kmer = u64::from_le_bytes(buf[0..8].try_into().unwrap());
                let count = u32::from_le_bytes(buf[8..12].try_into().unwrap());
                self.current = Some((kmer, count));
            }
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => self.current = None,
            Err(e) => return Err(e),
        }
        Ok(())
    }
}

/// Build a genome-wide k-mer count index.
///
/// * `sample` — keep 1 in `sample` k-mers by hash; must be a power of two.
/// * `min_count` — drop k-mers occurring fewer than this many times. Since the
///   index exists to find *repeats*, singletons are the bulk of the table and
///   dropping them is what keeps the result loadable.
pub fn build(
    genome: &Path,
    k: usize,
    sample: u64,
    min_count: u32,
    buffer_kmers: usize,
    tmp_dir: &Path,
    out: &Path,
    verbose: bool,
) -> io::Result<BuildStats> {
    if !(11..=32).contains(&k) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("k must be between 11 and 32 (got {})", k),
        ));
    }
    if !sample.is_power_of_two() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--sample must be a power of two (got {})", sample),
        ));
    }
    fs::create_dir_all(tmp_dir)?;

    let sample_mask = sample - 1;
    let mut stats = BuildStats::default();
    let mut buf: Vec<u64> = Vec::with_capacity(buffer_kmers);

    // ---- Pass 1: stream the genome, extract sampled k-mers, spill runs. ----
    let mut stream = FastaStream::open(genome)?;
    while let Some(rec) = stream.next_record(true)? {
        stats.records += 1;
        stats.bases += rec.seq.len() as u64;
        let seq = &rec.seq;
        let n_starts = seq.len().saturating_sub(k - 1);

        let mut s = 0usize;
        while s < n_starts {
            let e = (s + BLOCK).min(n_starts);

            // Disjoint k-mer *start* ranges: block [s, e) reads the slice
            // seq[s .. e+k-1] so every k-mer is emitted exactly once.
            let mut bounds = Vec::new();
            let sub = (e - s).div_ceil(rayon::current_num_threads().max(1));
            let sub = sub.max(1 << 20);
            let mut b = s;
            while b < e {
                bounds.push((b, (b + sub).min(e)));
                b += sub;
            }

            let parts: Vec<(u64, Vec<u64>)> = bounds
                .into_par_iter()
                .map(|(bs, be)| {
                    let hi = (be + k - 1).min(seq.len());
                    let mut local = Vec::new();
                    let mut total = 0u64;
                    for (_pos, kmer) in KmerIter::new(&seq[bs..hi], k) {
                        total += 1;
                        if is_sampled(kmer, sample_mask) {
                            local.push(kmer);
                        }
                    }
                    (total, local)
                })
                .collect();

            for (total, local) in parts {
                stats.total_kmers += total;
                stats.sampled_kmers += local.len() as u64;
                buf.extend_from_slice(&local);
            }

            if buf.len() >= buffer_kmers {
                flush_run(&mut buf, tmp_dir, stats.runs)?;
                stats.runs += 1;
                if verbose {
                    eprintln!(
                        "  spilled run {} ({} records seen, {} Mbp)",
                        stats.runs,
                        stats.records,
                        stats.bases / 1_000_000
                    );
                }
            }
            s = e;
        }
        if verbose {
            eprintln!("  counted {} ({} bp)", rec.header, rec.seq.len());
        }
    }
    if !buf.is_empty() {
        flush_run(&mut buf, tmp_dir, stats.runs)?;
        stats.runs += 1;
    }
    drop(buf);

    // ---- Pass 2: k-way merge the runs, summing counts for equal k-mers. ----
    let mut readers: Vec<RunReader> = (0..stats.runs)
        .map(|i| RunReader::open(&run_path(tmp_dir, i)))
        .collect::<io::Result<_>>()?;

    let tmp_k = tmp_dir.join("krep_merge_kmers.bin");
    let tmp_c = tmp_dir.join("krep_merge_counts.bin");
    {
        let mut wk = BufWriter::with_capacity(1 << 20, File::create(&tmp_k)?);
        let mut wc = BufWriter::with_capacity(1 << 20, File::create(&tmp_c)?);

        let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
        for (i, r) in readers.iter().enumerate() {
            if let Some((kmer, _)) = r.current {
                heap.push(Reverse((kmer, i)));
            }
        }

        while let Some(Reverse((kmer, _))) = heap.peek().copied() {
            let mut total: u64 = 0;
            // Drain every run positioned on this k-mer.
            while let Some(&Reverse((top, idx))) = heap.peek() {
                if top != kmer {
                    break;
                }
                heap.pop();
                total += readers[idx].current.map(|(_, c)| c as u64).unwrap_or(0);
                readers[idx].advance()?;
                if let Some((next, _)) = readers[idx].current {
                    heap.push(Reverse((next, idx)));
                }
            }

            let capped = total.min(u32::MAX as u64) as u32;
            // `capped` is u32, so leading_zeros() is 0..=31 for a non-zero value.
            let bucket = (31 - capped.max(1).leading_zeros()) as usize;
            stats.histogram[bucket.min(32)] += 1;

            if capped >= min_count {
                wk.write_all(&kmer.to_le_bytes())?;
                wc.write_all(&capped.to_le_bytes())?;
                stats.entries += 1;
            }
        }
        wk.flush()?;
        wc.flush()?;
    }
    drop(readers);

    // ---- Assemble the final file: header, then k-mers, then counts. ----
    {
        let mut w = BufWriter::with_capacity(1 << 20, File::create(out)?);
        w.write_all(MAGIC)?;
        w.write_all(&VERSION.to_le_bytes())?;
        w.write_all(&(k as u32).to_le_bytes())?;
        w.write_all(&sample.to_le_bytes())?;
        w.write_all(&min_count.to_le_bytes())?;
        w.write_all(&0u32.to_le_bytes())?; // padding
        w.write_all(&stats.entries.to_le_bytes())?;
        w.write_all(&stats.total_kmers.to_le_bytes())?;
        w.write_all(&stats.sampled_kmers.to_le_bytes())?;
        for p in [&tmp_k, &tmp_c] {
            let mut r = BufReader::with_capacity(1 << 20, File::open(p)?);
            io::copy(&mut r, &mut w)?;
        }
        w.flush()?;
    }

    for i in 0..stats.runs {
        let _ = fs::remove_file(run_path(tmp_dir, i));
    }
    let _ = fs::remove_file(&tmp_k);
    let _ = fs::remove_file(&tmp_c);

    Ok(stats)
}

/// A loaded genome-wide k-mer count index.
pub struct KmerIndex {
    pub k: usize,
    pub sample: u64,
    pub sample_mask: u64,
    pub min_count: u32,
    pub total_kmers: u64,
    pub sampled_kmers: u64,
    kmers: Vec<u64>,
    counts: Vec<u32>,
    /// `prefix[p]` is the first index whose top `PREFIX_BITS` are >= `p`, so a
    /// lookup binary-searches a small bucket instead of the whole table.
    prefix: Vec<u32>,
    prefix_shift: u32,
}

fn read_exact_vec<R: Read>(r: &mut R, n: usize, width: usize) -> io::Result<Vec<u8>> {
    let mut buf = vec![0u8; n * width];
    r.read_exact(&mut buf)?;
    Ok(buf)
}

impl KmerIndex {
    pub fn load(path: &Path) -> io::Result<Self> {
        let mut r = BufReader::with_capacity(1 << 22, File::open(path)?);
        let mut head = [0u8; HEADER_LEN];
        r.read_exact(&mut head)?;
        if &head[0..8] != MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "not a krep index file",
            ));
        }
        let version = u32::from_le_bytes(head[8..12].try_into().unwrap());
        if version != VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported index version {}", version),
            ));
        }
        let k = u32::from_le_bytes(head[12..16].try_into().unwrap()) as usize;
        let sample = u64::from_le_bytes(head[16..24].try_into().unwrap());
        let min_count = u32::from_le_bytes(head[24..28].try_into().unwrap());
        let n = u64::from_le_bytes(head[32..40].try_into().unwrap()) as usize;
        let total_kmers = u64::from_le_bytes(head[40..48].try_into().unwrap());
        let sampled_kmers = u64::from_le_bytes(head[48..56].try_into().unwrap());

        let raw_k = read_exact_vec(&mut r, n, 8)?;
        let kmers: Vec<u64> = raw_k
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        drop(raw_k);
        let raw_c = read_exact_vec(&mut r, n, 4)?;
        let counts: Vec<u32> = raw_c
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        drop(raw_c);

        let prefix_shift = (2 * k as u32).saturating_sub(PREFIX_BITS);
        let buckets = 1usize << PREFIX_BITS;
        let mut prefix = vec![0u32; buckets + 1];
        // Count per bucket, then prefix-sum into start offsets.
        for &kmer in &kmers {
            let p = (kmer >> prefix_shift) as usize;
            prefix[p.min(buckets - 1) + 1] += 1;
        }
        for i in 1..=buckets {
            prefix[i] += prefix[i - 1];
        }

        Ok(Self {
            k,
            sample,
            sample_mask: sample - 1,
            min_count,
            total_kmers,
            sampled_kmers,
            kmers,
            counts,
            prefix,
            prefix_shift,
        })
    }

    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    /// Approximate resident size of the loaded index, in bytes.
    pub fn memory_bytes(&self) -> usize {
        self.kmers.len() * 8 + self.counts.len() * 4 + self.prefix.len() * 4
    }

    /// Genomic occurrence count of a canonical k-mer.
    ///
    /// Returns 0 for k-mers outside the hash sample and for sampled k-mers that
    /// fell below `min_count` at build time — callers must therefore test
    /// `is_sampled` before treating a 0 as "not repetitive".
    #[inline]
    pub fn count(&self, kmer: u64) -> u32 {
        let bucket = (kmer >> self.prefix_shift) as usize;
        let last = self.prefix.len() - 2;
        let b = bucket.min(last);
        let lo = self.prefix[b] as usize;
        let hi = self.prefix[b + 1] as usize;
        match self.kmers[lo..hi].binary_search(&kmer) {
            Ok(i) => self.counts[lo + i],
            Err(_) => 0,
        }
    }

    #[inline]
    pub fn sampled(&self, kmer: u64) -> bool {
        is_sampled(kmer, self.sample_mask)
    }
}

/// Render the build histogram as a threshold-selection aid.
pub fn format_histogram(stats: &BuildStats) -> String {
    use std::fmt::Write as _;
    let mut s = String::new();
    let _ = writeln!(s, "\nDistinct sampled k-mers by occurrence count:");
    let mut cumulative_above = 0u64;
    let totals: Vec<u64> = stats.histogram.to_vec();
    for b in (0..33).rev() {
        if totals[b] == 0 {
            continue;
        }
        cumulative_above += totals[b];
        let lo = 1u64 << b;
        let hi = (1u64 << (b + 1)) - 1;
        let _ = writeln!(
            s,
            "  count {:>10}-{:<10} {:>12}   (>= {}: {})",
            lo, hi, totals[b], lo, cumulative_above
        );
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sampling_is_deterministic_and_roughly_uniform() {
        let mask = 16 - 1;
        let kept = (0u64..200_000).filter(|&x| is_sampled(x, mask)).count();
        // Expect ~1/16 of 200k = 12500; allow generous slack.
        assert!(
            (10_000..15_000).contains(&kept),
            "sampled {} of 200000, expected ~12500",
            kept
        );
        // Same k-mer must always give the same answer.
        for x in 0u64..1000 {
            assert_eq!(is_sampled(x, mask), is_sampled(x, mask));
        }
    }

    #[test]
    fn build_and_query_roundtrip() {
        let dir = std::env::temp_dir().join("krep_index_test");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        // A 400 bp unit repeated 12 times, embedded in unique flanks.
        let unit: Vec<u8> = (0..400)
            .map(|i| b"ACGT"[(splitmix64(i as u64) % 4) as usize])
            .collect();
        let mut seq = Vec::new();
        for i in 0..3000u64 {
            seq.push(b"ACGT"[(splitmix64(i ^ 0xfeed) % 4) as usize]);
        }
        for _ in 0..12 {
            seq.extend_from_slice(&unit);
        }
        for i in 0..3000u64 {
            seq.push(b"ACGT"[(splitmix64(i ^ 0xbeef) % 4) as usize]);
        }

        let fa = dir.join("t.fa");
        {
            let mut f = File::create(&fa).unwrap();
            writeln!(f, ">chrT").unwrap();
            for c in seq.chunks(80) {
                f.write_all(c).unwrap();
                f.write_all(b"\n").unwrap();
            }
        }

        let idx_path = dir.join("t.kidx");
        // sample=1 keeps every k-mer, so counts must be exact.
        let stats = build(&fa, 15, 1, 2, 1 << 16, &dir, &idx_path, false).unwrap();
        assert_eq!(stats.records, 1);

        let idx = KmerIndex::load(&idx_path).unwrap();
        assert_eq!(idx.k, 15);
        assert_eq!(idx.sample, 1);

        // Interior k-mers of the repeat unit occur exactly 12 times.
        let mid = 3000 + 400 * 6 + 100;
        let (_, kmer) = KmerIter::new(&seq[mid..mid + 15], 15).next().unwrap();
        assert_eq!(idx.count(kmer), 12, "repeat k-mer should occur 12 times");

        // A k-mer from the unique flank is below min_count and absent.
        let (_, uniq) = KmerIter::new(&seq[100..115], 15).next().unwrap();
        assert_eq!(idx.count(uniq), 0, "unique k-mer should not be stored");

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn histogram_buckets_are_log2() {
        let bucket = |c: u32| (31 - c.max(1).leading_zeros()) as usize;
        assert_eq!(bucket(1), 0);
        assert_eq!(bucket(2), 1);
        assert_eq!(bucket(3), 1);
        assert_eq!(bucket(4), 2);
        assert_eq!(bucket(12), 3);
        assert_eq!(bucket(255), 7);
        assert_eq!(bucket(256), 8);
        assert!(bucket(u32::MAX) <= 32);
    }

    #[test]
    fn rejects_non_power_of_two_sample() {
        let dir = std::env::temp_dir();
        let r = build(
            &dir.join("missing.fa"),
            18,
            7,
            2,
            1 << 16,
            &dir,
            &dir.join("x.kidx"),
            false,
        );
        assert!(r.is_err());
    }
}

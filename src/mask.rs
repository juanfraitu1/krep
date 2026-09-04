//! De novo repeat masking using a counting Bloom filter.
//!
//! The core masking pass is parallelized with rayon: the genome is split into
//! overlapping chunks, each chunk builds a local counting Bloom filter and
//! local coverage bitset, and the results are merged.

use crate::bitvec::BitVec;
use crate::cbf::CountingBloomFilter;
use crate::kmer::{self, KmerIter, SeedIter};
use crate::rng::Rng;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// A predicted repeat region.
#[derive(Debug, Clone)]
pub struct MaskedRegion {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    /// Family label for library hits ("consensus:score"); empty otherwise.
    pub name: String,
}

#[allow(dead_code)]
impl MaskedRegion {
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}

/// Disjoint-set (union-find) with path compression and union by size.
struct UnionFind {
    parent: Vec<usize>,
    size: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.size[ra] < self.size[rb] {
            self.parent[ra] = rb;
            self.size[rb] += self.size[ra];
        } else {
            self.parent[rb] = ra;
            self.size[ra] += self.size[rb];
        }
    }
}

/// Suggested chunk size for parallel processing. Chosen to keep per-chunk
/// overhead small while giving enough work per thread.
const CHUNK_SIZE: usize = 1_000_000;

/// Split a sequence into chunks suitable for parallel k-mer processing.
/// Each chunk overlaps the next by `overlap` bases so that k-mers crossing
/// chunk boundaries are not lost.
fn chunk_bounds(len: usize, overlap: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return Vec::new();
    }
    if len <= CHUNK_SIZE + overlap {
        return vec![(0, len)];
    }
    let mut bounds = Vec::new();
    let mut start = 0;
    while start < len {
        let end = (start + CHUNK_SIZE + overlap).min(len);
        bounds.push((start, end));
        if end == len {
            break;
        }
        start = end - overlap;
    }
    bounds
}

/// Count k-mers in parallel across chunks and merge into a single CBF.
fn parallel_count_kmers(
    seq: &[u8],
    k: usize,
    hash_seed: u64,
    cbf_factor: usize,
) -> CountingBloomFilter {
    let n_kmers = KmerIter::new(seq, k).count();
    let overlap = k - 1;
    let bounds = chunk_bounds(seq.len(), overlap);

    let local_cbfs: Vec<CountingBloomFilter> = bounds
        .into_par_iter()
        .map(|(start, end)| {
            let mut local =
                CountingBloomFilter::with_factor(n_kmers, 4, hash_seed, cbf_factor);
            for (_pos, kmer) in KmerIter::new(&seq[start..end], k) {
                local.insert(kmer);
            }
            local
        })
        .collect();

    // Merge local CBFs by saturating addition of counters.
    let mut cbf = CountingBloomFilter::with_factor(n_kmers, 4, hash_seed, cbf_factor);
    for local in local_cbfs {
        cbf.merge(&local);
    }
    cbf
}

/// Build a set of canonical k-mer keys whose CBF count is ≥ threshold.
fn build_high_kmer_set(
    seq: &[u8],
    k: usize,
    cbf: &CountingBloomFilter,
    threshold: u8,
) -> HashSet<u64> {
    let mut set = HashSet::new();
    for (_pos, kmer) in KmerIter::new(seq, k) {
        if cbf.count(kmer) >= threshold {
            set.insert(kmer);
        }
    }
    set
}

/// Query k-mers in parallel and build a global bitset of high-count k-mer starts.
///
/// If `mismatches > 0`, a k-mer position is marked if the k-mer or any of its
/// Hamming-distance-`mismatches` canonical neighbors has count ≥ threshold.
/// The optional `high_kmers` set must be supplied when `mismatches > 0` and
/// contains all k-mer keys with count ≥ threshold; it replaces slow per-neighbor
/// CBF queries with fast hash-set membership tests.
fn parallel_query_kmers(
    seq: &[u8],
    k: usize,
    cbf: &CountingBloomFilter,
    threshold: u8,
    mismatches: usize,
    high_kmers: Option<&HashSet<u64>>,
) -> BitVec {
    let overlap = k - 1;
    let bounds = chunk_bounds(seq.len(), overlap);
    let global_len = seq.len();

    let local_highs: Vec<(usize, BitVec)> = bounds
        .into_par_iter()
        .map(|(start, end)| {
            let local_len = end - start;
            let mut local_high = BitVec::new(local_len);
            for (pos, kmer) in KmerIter::new(&seq[start..end], k) {
                let hit = if mismatches == 0 {
                    cbf.count(kmer) >= threshold
                } else {
                    let set = high_kmers.expect("high_kmers required when mismatches > 0");
                    if set.contains(&kmer) {
                        true
                    } else {
                        kmer::neighbors(kmer, k, mismatches, true)
                            .iter()
                            .any(|n| set.contains(n))
                    }
                };
                if hit {
                    local_high.set(pos);
                }
            }
            (start, local_high)
        })
        .collect();

    let mut high = BitVec::new(global_len);
    for (start, local_high) in local_highs {
        for i in 0..local_high.len() {
            if local_high.get(i) {
                high.set(start + i);
            }
        }
    }
    high
}

/// Count k-mers, then query high positions, building a mismatch-expanded seed
/// set automatically when needed.
fn query_high_positions(
    seq: &[u8],
    k: usize,
    cbf: &CountingBloomFilter,
    threshold: u8,
    mismatches: usize,
) -> BitVec {
    let high_kmers = if mismatches > 0 {
        Some(build_high_kmer_set(seq, k, cbf, threshold))
    } else {
        None
    };
    parallel_query_kmers(seq, k, cbf, threshold, mismatches, high_kmers.as_ref())
}

/// Query k-mer starts against an explicit seed set.
///
/// A position is marked if its canonical k-mer (or any Hamming-distance-
/// `mismatches` neighbor) is contained in `seed_set`.
fn parallel_query_with_set(
    seq: &[u8],
    k: usize,
    seed_set: &HashSet<u64>,
    mismatches: usize,
) -> BitVec {
    let overlap = k - 1;
    let bounds = chunk_bounds(seq.len(), overlap);
    let global_len = seq.len();

    let local_highs: Vec<(usize, BitVec)> = bounds
        .into_par_iter()
        .map(|(start, end)| {
            let local_len = end - start;
            let mut local_high = BitVec::new(local_len);
            for (pos, kmer) in KmerIter::new(&seq[start..end], k) {
                let hit = if mismatches == 0 {
                    seed_set.contains(&kmer)
                } else {
                    if seed_set.contains(&kmer) {
                        true
                    } else {
                        kmer::neighbors(kmer, k, mismatches, true)
                            .iter()
                            .any(|n| seed_set.contains(n))
                    }
                };
                if hit {
                    local_high.set(pos);
                }
            }
            (start, local_high)
        })
        .collect();

    let mut high = BitVec::new(global_len);
    for (start, local_high) in local_highs {
        for i in 0..local_high.len() {
            if local_high.get(i) {
                high.set(start + i);
            }
        }
    }
    high
}

/// Run the de novo masker on one chromosome sequence.
///
/// * `chrom` — chromosome / contig name.
/// * `seq` — uppercase ACGT sequence (Ns are ignored).
/// * `k` — k-mer length.
/// * `threshold` — a k-mer is considered repetitive when its CBF count is ≥ this.
/// * `window` — length of the sequence window used to decide if a region is repetitive.
/// * `density` — minimum fraction of high-count k-mers required inside a window.
/// * `min_len` — drop predicted regions shorter than this.
pub fn mask_sequence(
    chrom: &str,
    seq: &[u8],
    k: usize,
    threshold: u8,
    window: usize,
    density: f64,
    min_len: usize,
    mismatches: usize,
    cbf_factor: usize,
) -> Vec<MaskedRegion> {
    if seq.len() < k || window == 0 || window > seq.len() {
        return Vec::new();
    }

    let hash_seed = Rng::new(k as u64).next_u64();

    // Pass 1: count k-mers in parallel.
    let cbf = parallel_count_kmers(seq, k, hash_seed, cbf_factor);

    // Pass 2: query k-mers in parallel and mark high-count starts.
    let high = query_high_positions(seq, k, &cbf, threshold, mismatches);

    // Sliding-window density (single-threaded; very fast).
    let mut covered = BitVec::new(seq.len());
    let mut high_count = 0usize;
    let n_window_kmers = window.saturating_sub(k - 1);
    let min_high = ((n_window_kmers as f64) * density).ceil() as usize;

    let max_start = seq.len() - window;
    for i in 0..=max_start {
        let entering = i;
        let leaving = i + window - k + 1;

        if i == 0 {
            for p in entering..leaving {
                if high.get(p) {
                    high_count += 1;
                }
            }
        } else {
            if high.get(entering - 1) {
                high_count -= 1;
            }
            if leaving > 0 && high.get(leaving - 1) {
                high_count += 1;
            }
        }

        if high_count >= min_high {
            for p in i..i + window {
                covered.set(p);
            }
        }
    }

    // Convert covered runs into regions.
    let mut regions = Vec::new();
    for (start, end) in covered.set_runs() {
        if end - start >= min_len {
            regions.push(MaskedRegion {
                chrom: crate::fasta::seq_id(chrom).to_string(),
                start,
                end,
                name: String::new(),
            });
        }
    }

    regions
}

/// Convenience: mask a multi-record FASTA, returning a flat list of regions.
#[allow(dead_code)]
pub fn mask_fasta(
    records: &[(String, Vec<u8>)],
    k: usize,
    threshold: u8,
    window: usize,
    density: f64,
    min_len: usize,
    mismatches: usize,
    cbf_factor: usize,
) -> Vec<MaskedRegion> {
    mask_fasta_union(
        records,
        &[k],
        threshold,
        window,
        density,
        min_len,
        mismatches,
        cbf_factor,
    )
}

/// Mask with multiple k-mer sizes and return the union of all masked regions.
/// Each k is processed in parallel across chunks; different k values run
/// sequentially.
pub fn mask_fasta_union(
    records: &[(String, Vec<u8>)],
    ks: &[usize],
    threshold: u8,
    window: usize,
    density: f64,
    min_len: usize,
    mismatches: usize,
    cbf_factor: usize,
) -> Vec<MaskedRegion> {
    let mut all = Vec::new();
    for (header, seq) in records {
        let mut covered = BitVec::new(seq.len());
        for k in ks {
            let regions = mask_sequence(
                header, seq, *k, threshold, window, density, min_len, mismatches, cbf_factor,
            );
            for r in regions {
                for p in r.start..r.end.min(seq.len()) {
                    covered.set(p);
                }
            }
        }
        for (start, end) in covered.set_runs() {
            if end - start >= min_len {
                all.push(MaskedRegion {
                    chrom: crate::fasta::seq_id(header).to_string(),
                    start,
                    end,
                name: String::new(),
                });
            }
        }
    }
    all
}

/// Graph-based masking.
///
/// For each k value, high-count k-mer start positions are collected. Positions
/// are linked if they are within `max_gap` bases of another high-count position.
/// Each connected component is expanded by `k - 1` bases at its end, and
/// components spanning at least `min_len` bases are emitted as masked regions.
///
/// This is a graph-theoretic way to bridge small low-count gaps inside repeat
/// copies and merge fragmented k-mer hits into coherent repeat blocks.
pub fn mask_fasta_graph(
    records: &[(String, Vec<u8>)],
    ks: &[usize],
    threshold: u8,
    max_gap: usize,
    min_len: usize,
    mismatches: usize,
    cbf_factor: usize,
) -> Vec<MaskedRegion> {
    let mut all = Vec::new();
    for (header, seq) in records {
        if seq.is_empty() {
            continue;
        }
        let mut covered = BitVec::new(seq.len());

        for k in ks {
            if seq.len() < *k {
                continue;
            }
            let hash_seed = Rng::new(*k as u64).next_u64();
            let cbf = parallel_count_kmers(seq, *k, hash_seed, cbf_factor);
            let high = query_high_positions(seq, *k, &cbf, threshold, mismatches);

            // Sort high positions and merge those within max_gap.
            let mut positions: Vec<usize> = Vec::with_capacity(high.count_ones());
            for p in 0..seq.len() {
                if high.get(p) {
                    positions.push(p);
                }
            }
            if positions.is_empty() {
                continue;
            }
            positions.sort_unstable();

            let mut comp_start = positions[0];
            let mut comp_last = positions[0];
            let k_ext = k.saturating_sub(1);
            for &pos in positions.iter().skip(1) {
                if pos > comp_last + max_gap {
                    // Close current component.
                    let s = comp_start;
                    let e = (comp_last + k_ext).min(seq.len());
                    for p in s..e {
                        covered.set(p);
                    }
                    comp_start = pos;
                }
                comp_last = pos;
            }
            // Close final component.
            let s = comp_start;
            let e = (comp_last + k_ext).min(seq.len());
            for p in s..e {
                covered.set(p);
            }
        }

        for (start, end) in covered.set_runs() {
            if end - start >= min_len {
                all.push(MaskedRegion {
                    chrom: crate::fasta::seq_id(header).to_string(),
                    start,
                    end,
                name: String::new(),
                });
            }
        }
    }
    all
}

/// Build a de Bruijn-style graph over high-count k-mers and return only the
/// k-mers that belong to components whose total CBF abundance is above the
/// threshold. This filters isolated noisy k-mers that happen to be high-count
/// in random background sequence.
fn dbg_retain_kmers(
    seq: &[u8],
    k: usize,
    seed_set: &HashSet<u64>,
    cbf: &CountingBloomFilter,
    abundance_threshold: u64,
) -> HashSet<u64> {
    // Assign a dense index to every seed k-mer so union-find is compact.
    let mut index: HashMap<u64, usize> = HashMap::with_capacity(seed_set.len());
    for &kmer in seed_set {
        let idx = index.len();
        index.insert(kmer, idx);
    }

    let mut uf = UnionFind::new(index.len());

    // Connect consecutive genomic k-mers that are both high-count seeds.
    let mut prev_pos = 0usize;
    let mut prev_kmer: Option<u64> = None;
    for (pos, kmer) in KmerIter::new(seq, k) {
        if let Some(pk) = prev_kmer {
            if pos == prev_pos + 1 {
                if let Some(&ia) = index.get(&pk) {
                    if let Some(&ib) = index.get(&kmer) {
                        uf.union(ia, ib);
                    }
                }
            }
        }
        prev_pos = pos;
        prev_kmer = Some(kmer);
    }

    // Compute total CBF abundance per connected component.
    let mut comp_abundance: HashMap<usize, u64> = HashMap::new();
    for (&kmer, &idx) in &index {
        let root = uf.find(idx);
        let cnt = cbf.count(kmer) as u64;
        *comp_abundance.entry(root).or_insert(0) += cnt;
    }

    // Keep k-mers in components that pass the abundance filter.
    let mut retained = HashSet::new();
    for (&kmer, &idx) in &index {
        let root = uf.find(idx);
        if comp_abundance[&root] >= abundance_threshold {
            retained.insert(kmer);
        }
    }
    retained
}

/// Assembly-based masking.
///
/// High-count k-mers are linked into connected components using a de Bruijn
/// graph of consecutive genomic k-mers. Only components whose total CBF
/// abundance is above `assembly_abundance` are kept. This lets the masker use
/// a lower per-k-mer threshold while still filtering random background noise.
///
/// Positions whose k-mer (or mismatch-1 neighbor, if enabled) is in a retained
/// component are masked, and the full k-mer window is filled in.
#[allow(dead_code)]
pub fn mask_fasta_assembly(
    records: &[(String, Vec<u8>)],
    ks: &[usize],
    threshold: u8,
    min_len: usize,
    mismatches: usize,
    assembly_abundance: u64,
    cbf_factor: usize,
) -> Vec<MaskedRegion> {
    let mut all = Vec::new();
    for (header, seq) in records {
        if seq.is_empty() {
            continue;
        }
        let mut covered = BitVec::new(seq.len());

        for k in ks {
            if seq.len() < *k {
                continue;
            }
            let hash_seed = Rng::new(*k as u64).next_u64();
            let cbf = parallel_count_kmers(seq, *k, hash_seed, cbf_factor);
            let seed_set = build_high_kmer_set(seq, *k, &cbf, threshold);
            if seed_set.is_empty() {
                continue;
            }

            let abundance_threshold = if assembly_abundance > 0 {
                assembly_abundance
            } else {
                (threshold as u64 * 5).max(20)
            };
            let retained = dbg_retain_kmers(seq, *k, &seed_set, &cbf, abundance_threshold);
            if retained.is_empty() {
                continue;
            }

            let high = parallel_query_with_set(seq, *k, &retained, mismatches);
            let k_ext = k.saturating_sub(1) + 1; // cover the whole k-mer window
            for p in 0..seq.len() {
                if high.get(p) {
                    let e = (p + k_ext).min(seq.len());
                    for q in p..e {
                        covered.set(q);
                    }
                }
            }
        }

        for (start, end) in covered.set_runs() {
            if end - start >= min_len {
                all.push(MaskedRegion {
                    chrom: crate::fasta::seq_id(header).to_string(),
                    start,
                    end,
                name: String::new(),
                });
            }
        }
    }
    all
}

/// Convert regions to a simple BED string.
pub fn regions_to_bed(regions: &[MaskedRegion]) -> String {
    let mut s = String::new();
    use std::fmt::Write as _;
    let named = regions.iter().any(|r| !r.name.is_empty());
    for r in regions {
        if named {
            let n = if r.name.is_empty() { "." } else { r.name.as_str() };
            writeln!(s, "{}\t{}\t{}\t{}", r.chrom, r.start, r.end, n).unwrap();
        } else {
            writeln!(s, "{}\t{}\t{}", r.chrom, r.start, r.end).unwrap();
        }
    }
    s
}

/// Convert regions to GTF (1-based inclusive coordinates).
pub fn regions_to_gtf(regions: &[MaskedRegion]) -> String {
    let mut s = String::new();
    use std::fmt::Write as _;
    for (i, r) in regions.iter().enumerate() {
        writeln!(
            s,
            "{}\tkrep\trepeat\t{}\t{}\t.\t.\t.\trepeat_id \"{}\";",
            r.chrom,
            r.start + 1,
            r.end,
            i
        )
        .unwrap();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masks_tandem_repeat() {
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        for _ in 0..50 {
            seq.extend_from_slice(b"ACACACACACACACACACAC");
        }
        seq.extend_from_slice(b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG");

        let regions = mask_sequence("chr1", &seq, 15, 3, 100, 0.5, 50, 0, 8);
        let total_masked: usize = regions.iter().map(|r| r.len()).sum();
        assert!(
            total_masked >= 900,
            "expected at least 900 bp masked, got {}",
            total_masked
        );
    }

    #[test]
    fn empty_sequence() {
        let regions = mask_sequence("chr1", b"", 31, 3, 100, 0.5, 50, 0, 8);
        assert!(regions.is_empty());
    }

    #[test]
    fn too_short_for_window() {
        let regions = mask_sequence("chr1", b"ACGTACGTACGT", 31, 3, 100, 0.5, 50, 0, 8);
        assert!(regions.is_empty());
    }
}

// ---------------------------------------------------------------------------
// Index-driven masking
// ---------------------------------------------------------------------------

/// One seed hit while masking against an index.
#[derive(Clone, Copy, Debug)]
struct SeedHit {
    /// Window start.
    pos: u32,
    /// Window end (exclusive) — `pos + span` of the index that produced it.
    end: u32,
    /// Count reached the nucleation threshold, not just the extension one.
    strong: bool,
}

/// Mask one record using **genome-wide** seed counts from prebuilt indices.
///
/// This is the counterpart to `mask_fasta_graph`, with one crucial difference:
/// the abundance of every seed was measured across the entire genome rather
/// than across whichever record happens to be in memory. A family with a few
/// hundred copies genome-wide is therefore visible even when the target is a
/// single chromosome or a 10 Mb slice.
///
/// Several indices may be supplied (say a contiguous 18-mer index and a
/// spaced-seed index); their hits are unioned before linking.
///
/// **Hysteresis.** A single count threshold cannot serve two jobs at once: low
/// enough to follow a repeat copy through its mutated stretches, yet high
/// enough that random near-unique seeds do not nucleate regions of their own.
/// So there are two. Seeds at or above `t_low` are collected and linked within
/// `max_gap`; a linked component is kept only if it contains at least one seed
/// at or above `t_high`. Weak seeds extend and bridge, never nucleate. With
/// `t_high == t_low` this reduces to the plain single-threshold linker.
///
/// Because indices may store only a hash sample of seed space, hits can be
/// sparse (roughly one per `sample` bases inside a repeat). `max_gap` links
/// them back into contiguous blocks, so it must comfortably exceed the
/// sampling stride.
///
/// **Density.** With dense sampling (every position tested) a genuine repeat
/// copy produces hits at a large fraction of its positions, while a cluster of
/// composition-biased background seeds does not. A component is kept only if
/// it has at least `min_hits` hits and at least `min_density` hits per base;
/// `min_hits = 1, min_density = 0.0` disables the filter.
#[allow(clippy::too_many_arguments)]
pub fn mask_sequence_indexed(
    chrom: &str,
    seq: &[u8],
    indices: &[crate::index::KmerIndex],
    thresholds: &[(u32, u32)],
    max_gap: usize,
    min_len: usize,
    min_hits: usize,
    min_density: f64,
) -> Vec<MaskedRegion> {
    assert_eq!(
        indices.len(),
        thresholds.len(),
        "one (t_high, t_low) pair per index"
    );
    let mut hits: Vec<SeedHit> = Vec::new();

    for (idx, &(t_high, t_low)) in indices.iter().zip(thresholds) {
        // A count of 0 means "absent from the index" (below its min_count),
        // never "occurs zero times" — so neither threshold may reach 0, or
        // every seed outside the table would register as a hit. The extension
        // threshold also cannot exceed the nucleation one. Thresholds are per
        // index because a weight-16 spaced seed and a contiguous 18-mer sit on
        // very different count scales in the same genome.
        let t_high = t_high.max(1);
        let t_low = t_low.clamp(1, t_high);
        let span = idx.span();
        let n_starts = seq.len().saturating_sub(span.saturating_sub(1));
        if n_starts == 0 {
            continue;
        }

        // Disjoint seed start ranges, so every position is tested exactly once.
        let block = (n_starts.div_ceil(rayon::current_num_threads().max(1))).max(1 << 20);
        let mut bounds = Vec::new();
        let mut b = 0usize;
        while b < n_starts {
            bounds.push((b, (b + block).min(n_starts)));
            b += block;
        }

        let parts: Vec<Vec<SeedHit>> = bounds
            .into_par_iter()
            .map(|(bs, be)| {
                let hi = (be + span - 1).min(seq.len());
                let mut local = Vec::new();
                for (pos, kmer) in SeedIter::new(&seq[bs..hi], &idx.seed) {
                    if !idx.sampled(kmer) {
                        continue;
                    }
                    let c = idx.count(kmer);
                    if c >= t_low {
                        let p = bs + pos;
                        local.push(SeedHit {
                            pos: p as u32,
                            end: (p + span).min(seq.len()) as u32,
                            strong: c >= t_high,
                        });
                    }
                }
                local
            })
            .collect();
        for p in parts {
            hits.extend_from_slice(&p);
        }
    }
    if hits.is_empty() {
        return Vec::new();
    }
    // A single index yields hits already in order; the union of several
    // interleaves and must be sorted.
    if indices.len() > 1 {
        hits.par_sort_unstable_by_key(|h| h.pos);
    }
    debug_assert!(hits.windows(2).all(|w| w[0].pos <= w[1].pos));

    let chrom_id = crate::fasta::seq_id(chrom);
    let mut regions = Vec::new();
    let mut comp_start = hits[0].pos as usize;
    let mut comp_end = hits[0].end as usize;
    let mut comp_last = hits[0].pos as usize;
    let mut comp_strong = hits[0].strong;
    let mut comp_hits = 1usize;

    let mut close = |s: usize, e: usize, strong: bool, n: usize, out: &mut Vec<MaskedRegion>| {
        let len = e - s;
        if strong && len >= min_len && n >= min_hits && n as f64 >= min_density * len as f64 {
            out.push(MaskedRegion {
                chrom: chrom_id.to_string(),
                start: s,
                end: e,
                name: String::new(),
            });
        }
    };

    for h in hits.iter().skip(1) {
        let pos = h.pos as usize;
        if pos > comp_last + max_gap {
            close(comp_start, comp_end, comp_strong, comp_hits, &mut regions);
            comp_start = pos;
            comp_end = h.end as usize;
            comp_strong = h.strong;
            comp_hits = 1;
        } else {
            comp_end = comp_end.max(h.end as usize);
            comp_strong |= h.strong;
            comp_hits += 1;
        }
        comp_last = pos;
    }
    close(comp_start, comp_end, comp_strong, comp_hits, &mut regions);

    regions
}

/// Soft-mask a sequence in place from regions belonging to this record.
pub fn apply_soft_mask(seq: &mut [u8], regions: &[MaskedRegion]) {
    let len = seq.len();
    for r in regions {
        for b in &mut seq[r.start.min(len)..r.end.min(len)] {
            b.make_ascii_lowercase();
        }
    }
}

#[cfg(test)]
pub(super) mod hysteresis_tests {
    use super::*;
    use crate::index::{self, KmerIndex};
    use crate::kmer::SpacedSeed;
    use std::fs::{self, File};
    use std::io::Write;

    /// Genome: unique flank, a 300 bp unit x 10 (strong), unique spacer, a
    /// 300 bp unit x 3 (weak), unique flank. Returns (path, seq, offsets).
    pub(super) fn build_fixture(name: &str) -> (KmerIndex, Vec<u8>, (usize, usize, usize, usize)) {
        // Tests run concurrently, so each gets its own scratch directory.
        let dir = std::env::temp_dir().join(format!("krep_hysteresis_{}", name));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        let rnd = |salt: u64, n: usize| -> Vec<u8> {
            (0..n as u64)
                .map(|i| b"ACGT"[(index::splitmix64(i ^ salt) % 4) as usize])
                .collect()
        };
        let strong_unit = rnd(0x1111, 300);
        let weak_unit = rnd(0x2222, 300);

        let mut seq = rnd(0xaaaa, 2000);
        let strong_start = seq.len();
        for _ in 0..10 {
            seq.extend_from_slice(&strong_unit);
        }
        let strong_end = seq.len();
        seq.extend_from_slice(&rnd(0xbbbb, 2000));
        let weak_start = seq.len();
        for _ in 0..3 {
            seq.extend_from_slice(&weak_unit);
        }
        let weak_end = seq.len();
        seq.extend_from_slice(&rnd(0xcccc, 2000));

        let fa = dir.join("h.fa");
        {
            let mut f = File::create(&fa).unwrap();
            writeln!(f, ">chrH").unwrap();
            f.write_all(&seq).unwrap();
            writeln!(f).unwrap();
        }
        let out = dir.join("h.kidx");
        index::build(&fa, &SpacedSeed::contiguous(15), 1, 2, 1 << 16, 1, &dir, &out, false)
            .unwrap();
        let idx = KmerIndex::load(&out).unwrap();
        (idx, seq, (strong_start, strong_end, weak_start, weak_end))
    }

    fn covered(regions: &[MaskedRegion], s: usize, e: usize) -> usize {
        regions
            .iter()
            .map(|r| r.end.min(e).saturating_sub(r.start.max(s)))
            .sum()
    }

    #[test]
    fn weak_seeds_never_nucleate() {
        let (idx, seq, (ss, se, ws, we)) = build_fixture("nucleate");
        let idxs = [idx];
        // t_high=8 (only the x10 array qualifies), t_low=3 (the x3 array has
        // count 3 and would pass a single threshold of 3).
        let regions = mask_sequence_indexed("chrH", &seq, &idxs, &[(8, 3)], 50, 30, 1, 0.0);
        assert!(covered(&regions, ss, se) > (se - ss) * 9 / 10, "strong array masked");
        assert_eq!(covered(&regions, ws, we), 0, "weak-only array must not nucleate");

        // Sanity: a plain low threshold *would* mask the weak array.
        let plain = mask_sequence_indexed("chrH", &seq, &idxs, &[(3, 3)], 50, 30, 1, 0.0);
        assert!(covered(&plain, ws, we) > (we - ws) * 9 / 10);
    }

    #[test]
    fn equal_thresholds_match_plain_linking() {
        let (idx, seq, _) = build_fixture("equal");
        let idxs = [idx];
        let a = mask_sequence_indexed("chrH", &seq, &idxs, &[(5, 5)], 50, 30, 1, 0.0);
        // t_low above t_high is clamped down to it, so 7 behaves as 5.
        let b = mask_sequence_indexed("chrH", &seq, &idxs, &[(5, 7)], 50, 30, 1, 0.0);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!((x.start, x.end), (y.start, y.end));
        }
        // A t_low of 0 must not turn absent seeds (count 0) into hits: the
        // unique flanks stay unmasked.
        let c = mask_sequence_indexed("chrH", &seq, &idxs, &[(5, 0)], 50, 30, 1, 0.0);
        assert!(covered(&c, 0, 1900) == 0, "unique flank must stay unmasked");
    }
}

#[cfg(test)]
mod density_tests {
    use super::*;

    #[test]
    fn density_filter_drops_sparse_components() {
        // Build hits by hand through the public path is awkward; instead reuse
        // the hysteresis fixture and assert the filter's monotonicity: a
        // demanding density must never mask *more* than a permissive one.
        let (idx, seq, (ss, se, _, _)) = super::hysteresis_tests::build_fixture("density");
        let idxs = [idx];
        let loose = mask_sequence_indexed("chrH", &seq, &idxs, &[(5, 5)], 50, 30, 1, 0.0);
        let tight = mask_sequence_indexed("chrH", &seq, &idxs, &[(5, 5)], 50, 30, 50, 0.5);
        let total = |r: &[MaskedRegion]| r.iter().map(|x| x.end - x.start).sum::<usize>();
        assert!(total(&tight) <= total(&loose));
        // The x10 array is hit at every sampled position (sample=1), so it
        // survives even a strict density requirement.
        let cov: usize = tight
            .iter()
            .map(|r| r.end.min(se).saturating_sub(r.start.max(ss)))
            .sum();
        assert!(cov > (se - ss) * 9 / 10, "dense true repeat must survive density filter");
    }
}

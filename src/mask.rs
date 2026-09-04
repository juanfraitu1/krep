//! De novo repeat masking using a counting Bloom filter.
//!
//! The core masking pass is parallelized with rayon: the genome is split into
//! overlapping chunks, each chunk builds a local counting Bloom filter and
//! local coverage bitset, and the results are merged.

use crate::bitvec::BitVec;
use crate::cbf::CountingBloomFilter;
use crate::kmer::{self, KmerIter};
use crate::rng::Rng;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// A predicted repeat region.
#[derive(Debug, Clone)]
pub struct MaskedRegion {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
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
                chrom: chrom.to_string(),
                start,
                end,
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
                    chrom: header.clone(),
                    start,
                    end,
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
                    chrom: header.clone(),
                    start,
                    end,
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
                    chrom: header.clone(),
                    start,
                    end,
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
    for r in regions {
        writeln!(s, "{}\t{}\t{}", r.chrom, r.start, r.end).unwrap();
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

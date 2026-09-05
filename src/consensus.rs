//! De novo repeat consensus building, RepeatScout-style.
//!
//! K-mer abundance finds a repeat copy only if *that copy* still carries
//! high-count k-mers. A copy 30% diverged from its family consensus keeps
//! 0.7^18 ≈ 0.16% of them, so ancient families (MIR, L2, old L1) are
//! invisible to `krep mask` at any threshold. RepeatMasker sidesteps this by
//! aligning every copy to a *consensus*, which is far closer to each copy than
//! copies are to one another. This module builds such consensi without a
//! library:
//!
//! 1. Seeds are high-count k-mers from a `krep index` (most abundant first).
//! 2. One streaming pass finds each seed's genomic occurrences (reservoir
//!    capped) and writes a plain sequence dump so windows around them can be
//!    fetched by random access later.
//! 3. For each seed not already inside a built consensus, ~100 occurrence
//!    windows are oriented by strand and the consensus is grown one base at a
//!    time. For each candidate base, one column of a banded "fit" alignment DP
//!    is advanced against every occurrence; the base maximising the sum of
//!    positive alignment scores wins. Growth stops once the total has not
//!    improved for `lookahead` steps and the consensus is truncated at its
//!    maximum. Left extension mirrors the right.
//! 4. Consensi are filtered for length, support, tandem periodicity and
//!    redundancy against those already built.

use crate::fasta::{self, FastaStream};
use crate::index::KmerIndex;
use crate::kmer::{encode_base, SpacedSeed};
use crate::rng::Rng;
use rayon::prelude::*;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::cmp::Reverse;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

pub struct Params {
    pub min_seed_count: u32,
    pub max_families: usize,
    /// Cap the number of seeds we collect occurrences for. Keeping every seed
    /// above `--min-seed-count` in memory can be huge on a dense spaced-seed
    /// index; limiting to the top `seed_pool` abundant seeds bounds the
    /// occurrence table. Zero means auto (`max_families * 20`).
    pub seed_pool: usize,
    pub max_occ: usize,
    pub flank: usize,
    pub band: usize,
    pub lookahead: usize,
    pub min_len: usize,
    pub min_support: usize,
    pub verbose: bool,
    /// Run a second pass over the same occurrence table with relaxed filters
    /// to recover low-copy or highly diverged families (e.g. CR1, Helitron)
    /// that were rejected in the first pass.
    pub second_pass: bool,
    /// Second-pass minimum consensus length.
    pub second_pass_min_len: usize,
    /// Second-pass minimum support (copies fitting the consensus).
    pub second_pass_min_support: usize,
    /// Maximum additional consensi to build in the second pass.
    pub second_pass_max_families: usize,
}

const NEG_INF: i32 = i32::MIN / 4;
const MATCH: i32 = 1;
const MISMATCH: i32 = -1;
/// Linear gap cost. It must exceed the mismatch penalty by a clear margin:
/// at gap 2 a diagonal switch costs the same as a mismatch, and the maximum
/// over a wide band of near-free paths drifts *upward* through random sequence
/// (the linear phase of alignment statistics), so extension never stops. At 3
/// every off-diagonal move is strictly worse than a mismatch and random
/// sequence drifts at about -0.5 per base, which is what makes "stop when the
/// total stops improving" a real stopping rule.
const GAP: i32 = 3;

// ---------------------------------------------------------------------------
// Sequence dump: plain uppercase bytes per record, for random access.
// ---------------------------------------------------------------------------

pub struct SeqDump {
    file: File,
    records: Vec<(String, u64, u64)>,
}

fn table_path(dump: &Path) -> PathBuf {
    let mut p = dump.as_os_str().to_owned();
    p.push(".tsv");
    PathBuf::from(p)
}

impl SeqDump {
    pub fn exists(path: &Path) -> bool {
        path.exists() && table_path(path).exists()
    }

    pub fn open(path: &Path) -> io::Result<Self> {
        let mut records = Vec::new();
        for line in BufReader::new(File::open(table_path(path))?).lines() {
            let line = line?;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() == 3 {
                records.push((
                    f[0].to_string(),
                    f[1].parse().map_err(|_| io::Error::other("bad dump table"))?,
                    f[2].parse().map_err(|_| io::Error::other("bad dump table"))?,
                ));
            }
        }
        Ok(Self {
            file: File::open(path)?,
            records,
        })
    }

    pub fn record_len(&self, rec: usize) -> usize {
        self.records[rec].2 as usize
    }

    /// Fetch `[start, end)` of record `rec`, clipped to the record.
    pub fn fetch(&mut self, rec: usize, start: usize, end: usize) -> io::Result<Vec<u8>> {
        let (_, off, len) = self.records[rec];
        let start = start.min(len as usize);
        let end = end.min(len as usize);
        let mut buf = vec![0u8; end - start];
        if !buf.is_empty() {
            self.file.seek(SeekFrom::Start(off + start as u64))?;
            self.file.read_exact(&mut buf)?;
        }
        Ok(buf)
    }
}

// ---------------------------------------------------------------------------
// Pass A: seed occurrences + dump.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct Occ {
    pub rec: u32,
    pub pos: u32,
    /// Forward k-mer at `pos` equals the canonical seed (else its reverse
    /// complement does, and the window must be flipped).
    pub forward: bool,
}

struct Reservoir {
    seen: u64,
    occs: Vec<Occ>,
}

/// Stream the genome once: write the sequence dump (unless it exists) and
/// collect up to `max_occ` occurrences per seed by reservoir sampling.
/// `seed` is the spaced (or contiguous) seed used by the index.
pub fn collect_occurrences(
    genome: &Path,
    dump: &Path,
    seeds: &HashSet<u64>,
    seed: &SpacedSeed,
    max_occ: usize,
    verbose: bool,
) -> io::Result<HashMap<u64, Vec<Occ>>> {
    let write_dump = !SeqDump::exists(dump);
    let mut dump_w = if write_dump {
        Some((
            BufWriter::with_capacity(1 << 20, File::create(dump)?),
            BufWriter::new(File::create(table_path(dump))?),
        ))
    } else {
        None
    };
    let mut offset = 0u64;
    let mut table: HashMap<u64, Reservoir> = HashMap::new();
    let mut rng = Rng::new(0xC0FFEE);
    let span = seed.span();
    let mask = if span == 32 { u64::MAX } else { (1u64 << (2 * span)) - 1 };
    let top = 2 * (span as u32 - 1);

    let mut stream = FastaStream::open(genome)?;
    let mut rec_idx = 0u32;
    while let Some(rec) = stream.next_record(true)? {
        if let Some((w, t)) = dump_w.as_mut() {
            w.write_all(&rec.seq)?;
            writeln!(t, "{}\t{}\t{}", fasta::seq_id(&rec.header), offset, rec.seq.len())?;
            offset += rec.seq.len() as u64;
        }
        let seq = &rec.seq;
        let n_starts = seq.len().saturating_sub(span - 1);
        let block = n_starts.div_ceil(rayon::current_num_threads().max(1)).max(1 << 20);
        let mut bounds = Vec::new();
        let mut b = 0;
        while b < n_starts {
            bounds.push((b, (b + block).min(n_starts)));
            b += block;
        }
        let found: Vec<Vec<(u64, Occ)>> = bounds
            .into_par_iter()
            .map(|(bs, be)| {
                let hi = (be + span - 1).min(seq.len());
                let s = &seq[bs..hi];
                let mut out = Vec::new();
                // Cap pushes per seed per block: a satellite array would
                // otherwise emit one hit per base for tens of megabases.
                // The global reservoir downsamples to max_occ later, so a
                // modest per-block cap is enough and keeps memory bounded on
                // dense spaced-seed indices.
                let mut per_seed: HashMap<u64, usize> = HashMap::new();
                let per_block_cap = max_occ.min(20);
                let (mut fwd, mut rev, mut valid) = (0u64, 0u64, 0usize);
                for (i, &base) in s.iter().enumerate() {
                    match encode_base(base) {
                        Some(bits) => {
                            fwd = ((fwd << 2) | bits) & mask;
                            rev = (rev >> 2) | ((3 - bits) << top);
                            valid += 1;
                            if valid >= span {
                                let fwd_seed = seed.extract(fwd);
                                let rev_seed = seed.extract(rev);
                                let canon = fwd_seed.min(rev_seed);
                                if seeds.contains(&canon) {
                                    let n = per_seed.entry(canon).or_insert(0);
                                    if *n >= per_block_cap {
                                        continue;
                                    }
                                    *n += 1;
                                    out.push((
                                        canon,
                                        Occ {
                                            rec: rec_idx,
                                            pos: (bs + i + 1 - span) as u32,
                                            forward: fwd_seed <= rev_seed,
                                        },
                                    ));
                                }
                            }
                        }
                        None => {
                            fwd = 0;
                            rev = 0;
                            valid = 0;
                        }
                    }
                }
                out
            })
            .collect();
        for part in found {
            for (seed, occ) in part {
                let r = table.entry(seed).or_insert_with(|| Reservoir {
                    seen: 0,
                    occs: Vec::new(),
                });
                r.seen += 1;
                if r.occs.len() < max_occ {
                    r.occs.push(occ);
                } else {
                    let j = rng.range_usize(r.seen as usize);
                    if j < max_occ {
                        r.occs[j] = occ;
                    }
                }
            }
        }
        if verbose {
            eprintln!("  scanned {} ({} bp)", fasta::seq_id(&rec.header), rec.seq.len());
        }
        rec_idx += 1;
    }
    if let Some((mut w, mut t)) = dump_w {
        w.flush()?;
        t.flush()?;
    }
    Ok(table.into_iter().map(|(s, r)| (s, r.occs)).collect())
}

// ---------------------------------------------------------------------------
// Greedy consensus extension.
// ---------------------------------------------------------------------------

/// Banded fit-alignment state for one occurrence: `col[d]` holds the score of
/// aligning the whole extension so far against the occurrence prefix of length
/// `j = m - band + d`, where `m` is the extension length.
struct Lane {
    seq: Vec<u8>,
    col: Vec<i32>,
}

/// Advance a lane by one consensus base into `out` and return the fit score.
#[inline]
fn advance(lane: &Lane, m: usize, band: usize, base: u8, out: &mut [i32]) -> i32 {
    let n = lane.seq.len();
    let width = 2 * band + 1;
    let mut best = NEG_INF;
    let new_m = m + 1;
    for d in 0..width {
        // j = new_m - band + d, as a signed quantity.
        let j = new_m as i64 - band as i64 + d as i64;
        if j < 0 || j as usize > n {
            out[d] = NEG_INF;
            continue;
        }
        let j = j as usize;
        let mut v = NEG_INF;
        // Diagonal: old column at j-1 is index d (same d since m shifted by 1).
        if j >= 1 {
            let old = lane.col[d];
            if old > NEG_INF {
                let s = if lane.seq[j - 1] == base { MATCH } else { MISMATCH };
                v = v.max(old + s);
            }
        }
        // Consensus base unmatched (gap in occurrence): old column at j is index d+1.
        if d + 1 < width {
            let old = lane.col[d + 1];
            if old > NEG_INF {
                v = v.max(old - GAP);
            }
        }
        // Extra occurrence base (gap in consensus): new column at j-1 is index d-1.
        if d >= 1 && out[d - 1] > NEG_INF {
            v = v.max(out[d - 1] - GAP);
        }
        out[d] = v;
        if v > best {
            best = v;
        }
    }
    best
}

/// Grow an extension to the right over `seqs` (each already the sequence that
/// follows the anchor, in consensus orientation). Returns the extension and
/// the fit score per occurrence at the chosen length.
fn extend(seqs: Vec<Vec<u8>>, band: usize, lookahead: usize, max_len: usize) -> (Vec<u8>, Vec<i32>) {
    let width = 2 * band + 1;
    let mut lanes: Vec<Lane> = seqs
        .into_iter()
        .map(|seq| {
            // m = 0: aligning nothing against a prefix of length j costs j gaps.
            let mut col = vec![NEG_INF; width];
            for d in 0..width {
                let j = d as i64 - band as i64;
                if j >= 0 && (j as usize) <= seq.len() {
                    col[d] = -(j as i32) * GAP;
                }
            }
            Lane { seq, col }
        })
        .collect();
    let n_lanes = lanes.len();
    // Voters choose each base; held-out judges decide whether the total
    // improved. Choosing and scoring on the same lanes lets the consensus fit
    // noise — with few lanes the plurality base "matches" enough of them that
    // the total creeps upward through random sequence and extension never
    // stops. Judges see a base chosen without them, so their drift is the
    // unbiased -0.5/base in random sequence and positive inside the family.
    let is_judge = |li: usize| n_lanes >= 3 && li % 3 == 2;
    let mut ext: Vec<u8> = Vec::new();
    let mut scores_at_best: Vec<i32> = vec![0; n_lanes];
    let mut best_total = 0i64;
    let mut best_m = 0usize;
    // Candidate columns for every lane and base, flat: [lane][base][width].
    let mut cands: Vec<i32> = vec![NEG_INF; n_lanes * 4 * width];
    let mut cand_scores: Vec<i32> = vec![0; n_lanes * 4];
    let mut m = 0usize;

    while m < max_len {
        let mut totals = [0i64; 4];
        for (li, lane) in lanes.iter().enumerate() {
            for (bi, &base) in b"ACGT".iter().enumerate() {
                let off = (li * 4 + bi) * width;
                let s = advance(lane, m, band, base, &mut cands[off..off + width]);
                cand_scores[li * 4 + bi] = s;
                if !is_judge(li) {
                    totals[bi] += s.max(0) as i64;
                }
            }
        }
        let mut chosen = 0usize;
        for bi in 1..4 {
            if totals[bi] > totals[chosen] {
                chosen = bi;
            }
        }
        let judged: i64 = (0..n_lanes)
            .filter(|&li| is_judge(li))
            .map(|li| cand_scores[li * 4 + chosen].max(0) as i64)
            .sum();
        for (li, lane) in lanes.iter_mut().enumerate() {
            let off = (li * 4 + chosen) * width;
            lane.col.copy_from_slice(&cands[off..off + width]);
        }
        ext.push(b"ACGT"[chosen]);
        m += 1;
        if judged > best_total {
            best_total = judged;
            best_m = m;
            for li in 0..n_lanes {
                scores_at_best[li] = cand_scores[li * 4 + chosen];
            }
        } else if m - best_m >= lookahead {
            break;
        }
    }
    ext.truncate(best_m);
    (ext, scores_at_best)
}

fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

/// True if the sequence is a tandem repeat: some period p (<= 50) explains
/// >= 80% of positions.
fn is_tandem(seq: &[u8]) -> bool {
    let n = seq.len();
    if n < 60 {
        return false;
    }
    for p in 1..=50.min(n / 3) {
        let agree = (0..n - p).filter(|&i| seq[i] == seq[i + p]).count();
        if agree * 5 >= (n - p) * 4 {
            return true;
        }
    }
    false
}

fn trim_runs(mut seq: Vec<u8>) -> Vec<u8> {
    let run_at = |s: &[u8], from_end: bool| -> usize {
        let mut it: Box<dyn Iterator<Item = &u8>> = if from_end {
            Box::new(s.iter().rev())
        } else {
            Box::new(s.iter())
        };
        let first = match it.next() {
            Some(&b) => b,
            None => return 0,
        };
        1 + it.take_while(|&&b| b == first).count()
    };
    let head = run_at(&seq, false);
    if head >= 8 {
        seq.drain(..head);
    }
    let tail = run_at(&seq, true);
    if tail >= 8 {
        let n = seq.len();
        seq.truncate(n - tail);
    }
    seq
}

fn kmer_set(seq: &[u8], k: usize) -> HashSet<u64> {
    crate::kmer::KmerIter::new(seq, k).map(|(_, x)| x).collect()
}

fn spaced_seed_set(seq: &[u8], seed: &SpacedSeed) -> HashSet<u64> {
    crate::kmer::SeedIter::new(seq, seed).map(|(_, x)| x).collect()
}

pub struct Consensus {
    pub id: usize,
    pub seq: Vec<u8>,
    pub support: usize,
    pub occurrences: usize,
    pub seed_count: u32,
}

/// Build consensi from the most abundant seeds down.
pub fn build_library(
    index: &KmerIndex,
    genome: &Path,
    dump_path: &Path,
    p: &Params,
) -> io::Result<Vec<Consensus>> {
    let seed = &index.seed;
    let span = seed.span();
    let weight = seed.weight();
    let seed_pool = if p.seed_pool > 0 {
        p.seed_pool
    } else {
        p.max_families.saturating_mul(20).max(100)
    };

    // Use a bounded min-heap to keep only the top `seed_pool` abundant seeds.
    // For a dense spaced-seed index this avoids materialising and sorting the
    // entire candidate list, which can be tens of millions of entries.
    let mut heap: BinaryHeap<Reverse<(u32, u64)>> = BinaryHeap::new();
    let mut total_candidates = 0u64;
    for (s, c) in index.entries() {
        if c < p.min_seed_count {
            continue;
        }
        total_candidates += 1;
        if heap.len() < seed_pool {
            heap.push(Reverse((c, s)));
        } else if c > heap.peek().unwrap().0 .0 {
            *heap.peek_mut().unwrap() = Reverse((c, s));
        }
    }
    let mut seeds: Vec<(u64, u32)> = heap
        .into_iter()
        .map(|Reverse((c, s))| (s, c))
        .collect();
    seeds.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    if p.verbose {
        eprintln!(
            "{} candidate seeds with count >= {}, using top {} for occurrence collection",
            total_candidates,
            p.min_seed_count,
            seed_pool
        );
    }
    let seed_set: HashSet<u64> = seeds.iter().map(|&(s, _)| s).collect();
    let occs = collect_occurrences(genome, dump_path, &seed_set, seed, p.max_occ, p.verbose)?;
    let mut dump = SeqDump::open(dump_path)?;

    let mut library: Vec<Consensus> = Vec::new();
    // canonical spaced seeds of built consensi
    let mut covered: HashSet<u64> = HashSet::new();
    let mut covered12: HashSet<u64> = HashSet::new(); // 12-mers, for redundancy
    let mut built = 0usize;
    let mut skipped_covered = 0usize;
    let mut rejected = 0usize;
    let mut candidates = 0usize;
    let max_candidates = p.max_families.saturating_mul(20).max(100);

    // -----------------------------------------------------------------------
    // First pass: build high-confidence consensi with the primary filters.
    // -----------------------------------------------------------------------
    let mut pass = 1;
    let mut pass_built = 0usize;
    let mut pass_candidates = 0usize;
    for &(seed, count) in &seeds {
        if library.len() >= p.max_families || candidates >= max_candidates {
            break;
        }
        // A subfamily seed differs from the consensus it belongs to by a
        // mismatch or two, so an exact check lets it through and the family is
        // rebuilt only to be thrown out as redundant. Testing the Hamming-2
        // neighbourhood (~1,400 k-mers for k=18) catches those before the
        // expensive extension. For spaced seeds this is an approximation but
        // still cheap.
        if covered.contains(&seed)
            || crate::kmer::neighbors(seed, weight, 2, true)
                .iter()
                .any(|n| covered.contains(n))
        {
            skipped_covered += 1;
            continue;
        }
        let list = match occs.get(&seed) {
            Some(l) if l.len() >= p.min_support => l,
            _ => continue,
        };
        candidates += 1;
        pass_candidates += 1;
        if p.verbose && pass_candidates % 100 == 0 {
            eprintln!(
                "  progress pass{}: {} candidates, {} built, {} rejected, {} seeds skipped",
                pass,
                pass_candidates,
                pass_built,
                rejected,
                skipped_covered
            );
        }

        // Fetch oriented windows around the seed span. We keep the actual
        // seed-span bases from the first occurrence to insert into the
        // consensus; the extension then grows outward on both sides.
        let mut lefts = Vec::with_capacity(list.len());
        let mut rights = Vec::with_capacity(list.len());
        let mut seed_span: Option<Vec<u8>> = None;
        for o in list {
            let rec = o.rec as usize;
            let pos = o.pos as usize;
            let start = pos.saturating_sub(p.flank);
            let end = (pos + span + p.flank).min(dump.record_len(rec));
            let win = dump.fetch(rec, start, end)?;
            let seed_off = pos - start;
            let (win, seed_off) = if o.forward {
                (win, seed_off)
            } else {
                let n = win.len();
                (revcomp_bytes(&win), n - seed_off - span)
            };
            if seed_span.is_none() {
                seed_span = Some(win[seed_off..seed_off + span].to_vec());
            }
            let mut left = win[..seed_off].to_vec();
            left.reverse();
            lefts.push(left);
            rights.push(win[seed_off + span..].to_vec());
        }

        let (right_ext, right_scores) = extend(rights, p.band, p.lookahead, p.flank);
        let (left_ext, left_scores) = extend(lefts, p.band, p.lookahead, p.flank);

        let mut seq: Vec<u8> = left_ext.iter().rev().copied().collect();
        if let Some(s) = seed_span {
            seq.extend_from_slice(&s);
        }
        seq.extend_from_slice(&right_ext);
        let seq = trim_runs(seq);
        let len = seq.len();

        // Support: copies fitting at roughly >= 62% identity over the consensus.
        let support = right_scores
            .iter()
            .zip(&left_scores)
            .filter(|(r, l)| (**r + **l + span as i32) as f64 >= 0.25 * len as f64)
            .count();

        let k12 = kmer_set(&seq, 12);
        let shared = k12.iter().filter(|x| covered12.contains(x)).count();
        let redundant = !k12.is_empty() && shared * 10 >= k12.len() * 6;

        let tandem = is_tandem(&seq);
        let pass_min_len = if pass == 1 { p.min_len } else { p.second_pass_min_len };
        let pass_min_support = if pass == 1 { p.min_support } else { p.second_pass_min_support };
        if len < pass_min_len || support < pass_min_support || tandem || redundant {
            rejected += 1;
            if p.verbose {
                eprintln!(
                    "  reject{}: seed count {:>6}  occ {:>3}  len {:>5} (L{}+S{}+R{})  support {:>3}  tandem {}  redundant {}",
                    pass,
                    count,
                    list.len(),
                    len,
                    left_ext.len(),
                    span,
                    right_ext.len(),
                    support,
                    tandem,
                    redundant
                );
            }
            continue;
        }

        // Accepted: mark this seed and the consensus spaced seeds as covered so
        // the same family is not rebuilt from another of its seeds.
        covered.insert(seed);
        let ks = spaced_seed_set(&seq, &index.seed);
        covered.extend(ks.iter().copied());
        covered12.extend(k12.iter().copied());

        built += 1;
        pass_built += 1;
        if p.verbose {
            eprintln!(
                "  family {:>5} (pass{}): len {:>5}  support {:>3}/{:<3}  seed count {}",
                built,
                pass,
                len,
                support,
                list.len(),
                count
            );
        }
        library.push(Consensus {
            id: built,
            seq,
            support,
            occurrences: list.len(),
            seed_count: count,
        });
    }

    // -----------------------------------------------------------------------
    // Second pass: revisit the same occurrence table with relaxed filters.
    // The pass-1 consensi already populate `covered` and `covered12`, so their
    // seeds are skipped and new pass-2 consensi are checked against them.
    // -----------------------------------------------------------------------
    if p.second_pass && p.second_pass_max_families > 0 {
        pass = 2;
        pass_built = 0;
        pass_candidates = 0;
        let total_budget = p.max_families + p.second_pass_max_families;
        let candidate_budget = max_candidates + p.second_pass_max_families.saturating_mul(20);
        for &(seed, count) in &seeds {
            if library.len() >= total_budget || candidates >= candidate_budget {
                break;
            }
            if covered.contains(&seed)
                || crate::kmer::neighbors(seed, weight, 2, true)
                    .iter()
                    .any(|n| covered.contains(n))
            {
                skipped_covered += 1;
                continue;
            }
            let list = match occs.get(&seed) {
                Some(l) if l.len() >= p.second_pass_min_support => l,
                _ => continue,
            };
            candidates += 1;
            pass_candidates += 1;
            if p.verbose && pass_candidates % 100 == 0 {
                eprintln!(
                    "  progress pass2: {} candidates, {} built, {} rejected, {} seeds skipped",
                    pass_candidates,
                    pass_built,
                    rejected,
                    skipped_covered
                );
            }

            let mut lefts = Vec::with_capacity(list.len());
            let mut rights = Vec::with_capacity(list.len());
            let mut seed_span: Option<Vec<u8>> = None;
            for o in list {
                let rec = o.rec as usize;
                let pos = o.pos as usize;
                let start = pos.saturating_sub(p.flank);
                let end = (pos + span + p.flank).min(dump.record_len(rec));
                let win = dump.fetch(rec, start, end)?;
                let seed_off = pos - start;
                let (win, seed_off) = if o.forward {
                    (win, seed_off)
                } else {
                    let n = win.len();
                    (revcomp_bytes(&win), n - seed_off - span)
                };
                if seed_span.is_none() {
                    seed_span = Some(win[seed_off..seed_off + span].to_vec());
                }
                let mut left = win[..seed_off].to_vec();
                left.reverse();
                lefts.push(left);
                rights.push(win[seed_off + span..].to_vec());
            }

            let (right_ext, right_scores) = extend(rights, p.band, p.lookahead, p.flank);
            let (left_ext, left_scores) = extend(lefts, p.band, p.lookahead, p.flank);

            let mut seq: Vec<u8> = left_ext.iter().rev().copied().collect();
            if let Some(s) = seed_span {
                seq.extend_from_slice(&s);
            }
            seq.extend_from_slice(&right_ext);
            let seq = trim_runs(seq);
            let len = seq.len();

            let support = right_scores
                .iter()
                .zip(&left_scores)
                .filter(|(r, l)| (**r + **l + span as i32) as f64 >= 0.25 * len as f64)
                .count();

            let k12 = kmer_set(&seq, 12);
            let shared = k12.iter().filter(|x| covered12.contains(x)).count();
            let redundant = !k12.is_empty() && shared * 10 >= k12.len() * 6;

            let tandem = is_tandem(&seq);
            if len < p.second_pass_min_len || support < p.second_pass_min_support || tandem || redundant {
                rejected += 1;
                if p.verbose {
                    eprintln!(
                        "  reject{}: seed count {:>6}  occ {:>3}  len {:>5} (L{}+S{}+R{})  support {:>3}  tandem {}  redundant {}",
                        pass,
                        count,
                        list.len(),
                        len,
                        left_ext.len(),
                        span,
                        right_ext.len(),
                        support,
                        tandem,
                        redundant
                    );
                }
                continue;
            }

            // Accepted in pass 2: mark this seed and its consensus as covered.
            covered.insert(seed);
            let ks = spaced_seed_set(&seq, &index.seed);
            covered.extend(ks.iter().copied());
            covered12.extend(k12.iter().copied());

            built += 1;
            pass_built += 1;
            if p.verbose {
                eprintln!(
                    "  family {:>5} (pass{}): len {:>5}  support {:>3}/{:<3}  seed count {}",
                    built,
                    pass,
                    len,
                    support,
                    list.len(),
                    count
                );
            }
            library.push(Consensus {
                id: built,
                seq,
                support,
                occurrences: list.len(),
                seed_count: count,
            });
        }
    }

    if p.verbose {
        eprintln!(
            "built {} consensi; {} seeds skipped as covered, {} candidates rejected",
            library.len(),
            skipped_covered,
            rejected
        );
    }
    Ok(library)
}

pub fn write_library(lib: &[Consensus], path: &Path) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for c in lib {
        writeln!(
            w,
            ">krep_fam_{:05} len={} support={} occurrences={} seed_count={}",
            c.id,
            c.seq.len(),
            c.support,
            c.occurrences,
            c.seed_count
        )?;
        for chunk in c.seq.chunks(80) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::splitmix64;

    fn rnd(salt: u64, n: usize) -> Vec<u8> {
        (0..n as u64)
            .map(|i| b"ACGT"[(splitmix64(i ^ salt) % 4) as usize])
            .collect()
    }

    /// Mutate a sequence at the given per-base substitution rate.
    fn mutate(seq: &[u8], rate: f64, salt: u64) -> Vec<u8> {
        seq.iter()
            .enumerate()
            .map(|(i, &b)| {
                let h = splitmix64(i as u64 ^ salt);
                if (h % 10_000) as f64 / 10_000.0 < rate {
                    let alt = b"ACGT"[((h >> 20) % 4) as usize];
                    if alt == b { b"ACGT"[((h >> 22) % 4) as usize] } else { alt }
                } else {
                    b
                }
            })
            .collect()
    }

    #[test]
    fn extension_recovers_consensus_from_diverged_copies() {
        // 40 copies of a 300 bp element, each 20% diverged, after an anchor.
        let truth = rnd(0x51, 300);
        let copies: Vec<Vec<u8>> = (0..40u64)
            .map(|i| {
                let mut c = mutate(&truth, 0.20, 0x1000 + i);
                c.extend_from_slice(&rnd(0x9000 + i, 200)); // unique tail
                c
            })
            .collect();
        let (ext, scores) = extend(copies, 16, 60, 600);
        assert!(ext.len() >= 280 && ext.len() <= 320, "length {}", ext.len());
        let n = ext.len().min(300);
        let agree = (0..n).filter(|&i| ext[i] == truth[i]).count();
        assert!(agree * 100 >= n * 97, "consensus identity {}/{}", agree, n);
        assert!(scores.iter().filter(|&&s| s > 0).count() >= 35);
    }

    #[test]
    fn extension_stops_in_random_sequence() {
        // Occurrences that share nothing after the anchor must not grow a
        // consensus: the total score has to fall, not wander upward.
        let copies: Vec<Vec<u8>> = (0..60u64).map(|i| rnd(0x7000 + i, 400)).collect();
        let (ext, _) = extend(copies, 16, 60, 400);
        assert!(ext.len() < 20, "random tails extended to {} bp", ext.len());
    }

    #[test]
    fn tandem_detection() {
        let mut t = Vec::new();
        for _ in 0..40 {
            t.extend_from_slice(b"ACGTG");
        }
        assert!(is_tandem(&t));
        assert!(!is_tandem(&rnd(7, 200)));
    }

    #[test]
    fn trims_long_runs_only() {
        let s = trim_runs(b"AAAAAAAAAACGTACGTTTTTTTTTT".to_vec());
        assert_eq!(s, b"CGTACG");
        let s = trim_runs(b"AAACGTACGTTT".to_vec());
        assert_eq!(s, b"AAACGTACGTTT");
    }
}

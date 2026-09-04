//! Tandem-repeat detection in the spirit of Tandem Repeats Finder.
//!
//! NCBI's soft-mask is WindowMasker plus TRF; krep's k-mer index already plays
//! the WindowMasker role, and this module supplies the other half. It also
//! lifts agreement with RepeatMasker's Simple_repeat / Low_complexity classes,
//! which abundance masking underserves.
//!
//! Method: a table of the last position of every k-mer gives, at each
//! position, the distance `d` back to the previous copy of that k-mer. In a
//! tandem array of period `d` that distance is `d` at most positions, so
//! matches at the same `d` form a run. A run is reported when it is long
//! enough, dense enough in matching k-mers, and — as verification — the
//! sequence agrees with itself shifted by `d` at a high enough fraction of
//! positions. Random sequence produces isolated matches at scattered
//! distances and never forms a dense run at one period.

use crate::fasta;
use crate::kmer::encode_base;
use crate::mask::MaskedRegion;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct TandemParams {
    /// k-mer length used to find recurrences (6 = 4,096-entry table).
    pub k: usize,
    /// Largest period considered.
    pub max_period: usize,
    /// Shortest array reported.
    pub min_len: usize,
    /// Minimum fraction of positions (beyond the first copy) whose k-mer
    /// recurs at the run's period.
    pub min_density: f64,
    /// Minimum fraction of positions agreeing with the base one period back.
    pub min_identity: f64,
}

impl Default for TandemParams {
    fn default() -> Self {
        Self {
            k: 5,
            max_period: 500,
            min_len: 30,
            min_density: 0.25,
            min_identity: 0.7,
        }
    }
}

#[derive(Clone, Copy)]
struct Run {
    start: usize,
    last: usize,
    matches: usize,
}

/// Fraction of positions in `seq[s..e]` equal to the base `d` earlier.
fn periodic_identity(seq: &[u8], s: usize, e: usize, d: usize) -> f64 {
    if e <= s + d {
        return 0.0;
    }
    let n = e - s - d;
    let agree = (s + d..e).filter(|&i| seq[i] == seq[i - d]).count();
    agree as f64 / n as f64
}

fn scan(seq: &[u8], lo: usize, hi: usize, p: &TandemParams, chrom: &str, out: &mut Vec<MaskedRegion>) {
    let k = p.k;
    let table_size = 1usize << (2 * k);
    let mask = (1u64 << (2 * k)) - 1;
    let mut last = vec![u32::MAX; table_size];
    let mut runs: Vec<Option<Run>> = vec![None; p.max_period + 1];
    let mut close = |d: usize, r: Run, out: &mut Vec<MaskedRegion>| {
        let s = r.start;
        let e = (r.last + k).min(seq.len());
        let len = e - s;
        if len < p.min_len || len < 2 * d {
            return;
        }
        // Positions after the first copy that could have matched.
        let candidates = len.saturating_sub(d).max(1);
        if (r.matches as f64) < p.min_density * candidates as f64 {
            return;
        }
        if periodic_identity(seq, s, e, d) < p.min_identity {
            return;
        }
        out.push(MaskedRegion {
            chrom: fasta::seq_id(chrom).to_string(),
            start: s,
            end: e,
            name: format!("tandem:period={}", d),
        });
    };

    let (mut fwd, mut valid) = (0u64, 0usize);
    for i in lo..hi {
        match encode_base(seq[i]) {
            Some(bits) => {
                fwd = ((fwd << 2) | bits) & mask;
                valid += 1;
            }
            None => {
                fwd = 0;
                valid = 0;
                continue;
            }
        }
        if valid < k {
            continue;
        }
        let pos = i + 1 - k;
        let prev = last[fwd as usize];
        last[fwd as usize] = pos as u32;
        if prev == u32::MAX {
            continue;
        }
        let d = pos - prev as usize;
        if d == 0 || d > p.max_period {
            continue;
        }
        // Extend the run for this period, or close it and start afresh.
        let slot = &mut runs[d];
        match slot {
            Some(r) if pos <= r.last + d + k => {
                r.last = pos;
                r.matches += 1;
            }
            _ => {
                if let Some(r) = slot.take() {
                    close(d, r, out);
                }
                *slot = Some(Run {
                    start: prev as usize,
                    last: pos,
                    matches: 1,
                });
            }
        }
    }
    for d in 1..=p.max_period {
        if let Some(r) = runs[d].take() {
            close(d, r, out);
        }
    }
}

/// Find tandem arrays in one record, in parallel blocks. Blocks overlap by
/// more than the longest array a run can span across a boundary, so every
/// array is seen whole by at least one block; duplicates merge downstream.
pub fn find_tandems(chrom: &str, seq: &[u8], p: &TandemParams) -> Vec<MaskedRegion> {
    if seq.len() < p.k {
        return Vec::new();
    }
    let overlap = 4 * p.max_period + 64;
    let block = (seq.len().div_ceil(rayon::current_num_threads().max(1))).max(1 << 20);
    let mut bounds = Vec::new();
    let mut b = 0usize;
    while b < seq.len() {
        bounds.push((b.saturating_sub(overlap), (b + block).min(seq.len())));
        b += block;
    }
    let parts: Vec<Vec<MaskedRegion>> = bounds
        .into_par_iter()
        .map(|(lo, hi)| {
            let mut out = Vec::new();
            scan(seq, lo, hi, p, chrom, &mut out);
            out
        })
        .collect();
    let mut all: Vec<MaskedRegion> = parts.into_iter().flatten().collect();
    crate::align::merge_regions(&mut all);
    all
}

/// DUST-style low-complexity detection.
///
/// RepeatMasker's Simple_repeat and Low_complexity classes are full of
/// sequence whose period wobbles (`CAAACAAAACAAACAAC...`) or that is merely
/// AT/GC-skewed. A fixed-period recurrence test cannot lock onto those; DUST
/// scores triplet-frequency skew over a sliding window instead:
/// `sum_t c_t (c_t - 1) / 2 / (l - 1)` over the 64 triplets, where `l` is the
/// number of triplets in the window. On this scale (measured on chr1): random
/// sequence ~0.5, AT-skewed sequence <1, wobbly period-3/4 microsatellites
/// 3-5, a perfect dinucleotide repeat ~15, a homopolymer ~31. A threshold of
/// 5 catches the low-complexity sequence RepeatMasker calls Simple_repeat /
/// Low_complexity; 3 catches more at a precision cost.
pub fn find_low_complexity(chrom: &str, seq: &[u8], window: usize, threshold: f64) -> Vec<MaskedRegion> {
    if seq.len() < window {
        return Vec::new();
    }
    let block = (seq.len().div_ceil(rayon::current_num_threads().max(1))).max(1 << 20);
    let mut bounds = Vec::new();
    let mut b = 0usize;
    while b < seq.len() {
        bounds.push((b.saturating_sub(window), (b + block).min(seq.len())));
        b += block;
    }
    let chrom_id = fasta::seq_id(chrom).to_string();
    let parts: Vec<Vec<MaskedRegion>> = bounds
        .into_par_iter()
        .map(|(lo, hi)| dust_scan(&chrom_id, seq, lo, hi, window, threshold))
        .collect();
    let mut all: Vec<MaskedRegion> = parts.into_iter().flatten().collect();
    crate::align::merge_regions(&mut all);
    all
}

fn dust_scan(chrom: &str, seq: &[u8], lo: usize, hi: usize, window: usize, threshold: f64) -> Vec<MaskedRegion> {
    let mut out: Vec<MaskedRegion> = Vec::new();
    // Triplet code per position (start of triplet), -1 across non-ACGT.
    let code = |i: usize| -> Option<usize> {
        let a = encode_base(seq[i])?;
        let b = encode_base(seq[i + 1])?;
        let c = encode_base(seq[i + 2])?;
        Some(((a << 4) | (b << 2) | c) as usize)
    };
    let l = window - 2; // triplets per window
    let denom = (l - 1) as f64;
    let mut counts = [0u32; 64];
    let mut pairs = 0u64; // sum of c(c-1)/2
    let mut valid = 0usize; // triplets currently in the window
    let mut ring: Vec<Option<usize>> = vec![None; l];
    let mut cur_start: Option<usize> = None;
    let push = |out: &mut Vec<MaskedRegion>, s: usize, e: usize| {
        out.push(MaskedRegion {
            chrom: chrom.to_string(),
            start: s,
            end: e,
            name: "dust".to_string(),
        });
    };
    let last_start = hi.saturating_sub(3);
    for i in lo..=last_start.min(seq.len() - 3) {
        let slot = i % l;
        if let Some(old) = ring[slot] {
            let c = counts[old];
            pairs -= (c - 1) as u64;
            counts[old] -= 1;
            valid -= 1;
        }
        let t = code(i);
        ring[slot] = t;
        if let Some(t) = t {
            pairs += counts[t] as u64;
            counts[t] += 1;
            valid += 1;
        }
        // Window of triplets [i-l+1, i] covers bases [i-l+1, i+3).
        if valid == l && i + 1 >= l {
            let score = pairs as f64 / denom;
            let ws = i + 1 - l;
            let we = i + 3;
            if score > threshold {
                match cur_start {
                    Some(_) => {}
                    None => cur_start = Some(ws),
                }
                // Extend the open region to the end of this window.
                if let Some(s) = cur_start {
                    if let Some(last) = out.last_mut() {
                        if last.start == s {
                            last.end = we;
                            continue;
                        }
                    }
                    push(&mut out, s, we);
                }
            } else {
                cur_start = None;
            }
        } else {
            cur_start = None;
        }
    }
    out
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

    fn covered(regions: &[MaskedRegion], s: usize, e: usize) -> usize {
        regions.iter().map(|r| r.end.min(e).saturating_sub(r.start.max(s))).sum()
    }

    #[test]
    fn finds_microsatellite_and_diverged_minisatellite() {
        let p = TandemParams::default();
        let mut seq = rnd(1, 2000);
        let ms_start = seq.len();
        for _ in 0..40 {
            seq.extend_from_slice(b"CA");
        }
        let ms_end = seq.len();
        seq.extend_from_slice(&rnd(2, 1500));
        // 37 bp unit x 12, each copy with ~8% substitutions.
        let unit = rnd(3, 37);
        let mini_start = seq.len();
        for c in 0..12u64 {
            for (i, &b) in unit.iter().enumerate() {
                let h = splitmix64(i as u64 ^ (c * 977));
                seq.push(if h % 100 < 8 { b"ACGT"[((h >> 8) % 4) as usize] } else { b });
            }
        }
        let mini_end = seq.len();
        seq.extend_from_slice(&rnd(4, 1500));

        let regions = find_tandems("chrT", &seq, &p);
        assert!(covered(&regions, ms_start, ms_end) * 10 >= (ms_end - ms_start) * 9, "microsat");
        assert!(covered(&regions, mini_start, mini_end) * 10 >= (mini_end - mini_start) * 8, "minisat");
        let total: usize = regions.iter().map(|r| r.end - r.start).sum();
        let truth = (ms_end - ms_start) + (mini_end - mini_start);
        assert!(total <= truth + 100, "over-masked: {} vs {}", total, truth);
    }

    #[test]
    fn random_sequence_is_clean() {
        let seq = rnd(9, 200_000);
        let regions = find_tandems("chrT", &seq, &TandemParams::default());
        let total: usize = regions.iter().map(|r| r.end - r.start).sum();
        assert!(total < 200, "random masked {} bp", total);
    }
}

#[cfg(test)]
mod dust_tests {
    use super::*;
    use crate::index::splitmix64;

    fn rnd(salt: u64, n: usize) -> Vec<u8> {
        (0..n as u64).map(|i| b"ACGT"[(splitmix64(i ^ salt) % 4) as usize]).collect()
    }

    #[test]
    fn dust_flags_wobbly_low_complexity_not_random() {
        let mut seq = rnd(1, 3000);
        let s = seq.len();
        seq.extend_from_slice(b"AATAATCAAACAAACAAACAAAAACAAAAACAAACAAACAACAACAACAAAAACAA");
        let e = seq.len();
        seq.extend_from_slice(&rnd(2, 3000));
        let regions = find_low_complexity("chrT", &seq, 64, 3.0);
        let cov: usize = regions.iter().map(|r| r.end.min(e).saturating_sub(r.start.max(s))).sum();
        assert!(cov * 10 >= (e - s) * 8, "low-complexity covered {} of {}", cov, e - s);
        let total: usize = regions.iter().map(|r| r.end - r.start).sum();
        assert!(total <= (e - s) + 2 * 64, "over-masked {}", total);
        let clean = find_low_complexity("chrT", &rnd(9, 100_000), 64, 3.0);
        assert!(clean.iter().map(|r| r.end - r.start).sum::<usize>() < 300);
    }
}

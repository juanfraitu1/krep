//! Library-based masking: align consensus sequences to the genome.
//!
//! This is the second half of the RepeatMasker idea. `krep consensus` builds
//! family consensi; here every consensus is seeded into a direct-address table
//! of spaced seeds (both strands), the genome is streamed, and two seed hits on
//! one diagonal trigger a banded X-drop alignment in both directions. Hits at
//! or above `min_score` are masked, and seeds inside a freshly masked region
//! are skipped, so an L1 copy costs one extension rather than hundreds.
//!
//! The seed is PatternHunter's weight-11 / span-18 pattern, chosen for
//! sensitivity around 70% identity. With 4^11 possible values the table is
//! addressed directly — no hashing on the 3.1e9-position hot path.

use crate::fasta;
use crate::kmer::{encode_base, SpacedSeed};
use crate::mask::MaskedRegion;
use rayon::prelude::*;
use std::io;
use std::path::Path;

pub const LIB_SEED: &str = "111010010100110111";
const MATCH: i32 = 1;
const MISMATCH: i32 = -1;
const GAP: i32 = 3;
const NEG_INF: i32 = i32::MIN / 4;
/// Two hits must fall within this many genome bases to chain.
const CHAIN_WINDOW: usize = 64;
/// ...and on diagonals this close.
const CHAIN_DIAG: i64 = 3;
/// Recent hits remembered per oriented consensus.
const SLOTS: usize = 2;
/// Single-hit mode: banded DPs attempted per genome position, taken from the
/// entries with the best ungapped gate scores.
const SINGLE_HIT_TOP: usize = 2;

/// One accepted library hit with the features a learned filter can use.
#[derive(Clone, Debug)]
pub struct Candidate {
    pub start: usize,
    pub end: usize,
    /// Consensus index (not oriented).
    pub cid: usize,
    pub strand: bool,
    pub score: i32,
    pub score_fwd: i32,
    pub score_bwd: i32,
    /// Consensus bases consumed forward / backward from the anchor.
    pub cons_fwd: usize,
    pub cons_bwd: usize,
    pub cons_len: usize,
    pub gate_left: i32,
    pub gate_right: i32,
}

pub struct Library {
    pub names: Vec<String>,
    /// Oriented consensi: index 2c is consensus c forward, 2c+1 its reverse
    /// complement.
    seqs: Vec<Vec<u8>>,
    seed: SpacedSeed,
    span: usize,
    /// `entries[table[v]..table[v+1]]` are the library positions with seed
    /// value `v`.
    table: Vec<u32>,
    entries: Vec<(u32, u32)>,
    pub min_score: i32,
    pub band: usize,
    pub xdrop: i32,
    /// Trigger extension on every seed hit (gated by the ungapped check)
    /// instead of requiring two hits on one diagonal. Far more sensitive to
    /// short diverged fragments; costs ~1 µs per genome position.
    pub single_hit: bool,
    /// Single-hit gate: pass when the summed flank score reaches `gate_sum`
    /// or either flank alone reaches `gate_side`. Lower is more sensitive to
    /// short fragments whose flanks run into random sequence, and slower.
    pub gate_sum: i32,
    pub gate_side: i32,
    /// Learned per-consensus minimum score (indexed by consensus id); a hit
    /// must clear both this and `min_score`. Empty means none.
    pub thresholds: Vec<i32>,
    /// Optional logistic filter over hit features with a per-consensus
    /// offset (see `train_logistic.py`).
    pub logistic: Option<Logistic>,
}

/// Logistic hit filter: p = sigmoid(b0 + offset[cid] + w . z(features)),
/// accept when p >= tau. Features, in order: score, ln(length),
/// score/length, consensus coverage, GC, gate sum, gate min, ln(consensus
/// length); `mu`/`sd` standardize them as at training time.
#[derive(Clone, Debug)]
pub struct Logistic {
    pub w: [f64; 8],
    pub b0: f64,
    pub tau: f64,
    pub mu: [f64; 8],
    pub sd: [f64; 8],
    pub offsets: Vec<f64>,
}

impl Logistic {
    #[inline]
    fn prob(&self, cid: usize, x: [f64; 8]) -> f64 {
        let mut z = self.b0 + self.offsets.get(cid).copied().unwrap_or(0.0);
        for j in 0..8 {
            z += self.w[j] * (x[j] - self.mu[j]) / self.sd[j];
        }
        1.0 / (1.0 + (-z.clamp(-30.0, 30.0)).exp())
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
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

/// A seed whose care bases use fewer than three symbols, or are dominated by
/// one, would match every poly-A tail and AT-rich stretch in the genome.
#[inline]
fn complex_enough(value: u64, weight: usize) -> bool {
    let mut counts = [0u8; 4];
    let mut v = value;
    for _ in 0..weight {
        counts[(v & 3) as usize] += 1;
        v >>= 2;
    }
    let distinct = counts.iter().filter(|&&c| c > 0).count();
    let max = *counts.iter().max().unwrap() as usize;
    distinct >= 3 && max + 3 <= weight
}

/// Forward-strand seed values of `seq`: (window start, value).
fn forward_seeds(seq: &[u8], seed: &SpacedSeed) -> Vec<(usize, u64)> {
    let span = seed.span();
    let mask = if span == 32 { u64::MAX } else { (1u64 << (2 * span)) - 1 };
    let mut out = Vec::new();
    let (mut fwd, mut valid) = (0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        match encode_base(b) {
            Some(bits) => {
                fwd = ((fwd << 2) | bits) & mask;
                valid += 1;
                if valid >= span {
                    out.push((i + 1 - span, seed.extract(fwd)));
                }
            }
            None => {
                fwd = 0;
                valid = 0;
            }
        }
    }
    out
}

impl Library {
    pub fn load(path: &Path, min_score: i32, band: usize, xdrop: i32) -> io::Result<Self> {
        Self::load_with_seed(path, LIB_SEED, min_score, band, xdrop)
    }

    pub fn load_with_seed(path: &Path, pattern: &str, min_score: i32, band: usize, xdrop: i32) -> io::Result<Self> {
        let records = fasta::read_uppercase(path)?;
        let seed = SpacedSeed::parse_any(pattern)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        if seed.weight() > 13 || seed.weight() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "library seed weight must be 8..=13 (direct-address table of 4^weight entries)",
            ));
        }
        let weight = seed.weight();
        let mut names = Vec::new();
        let mut seqs = Vec::new();
        for r in records {
            names.push(fasta::seq_id(&r.header).to_string());
            seqs.push(revcomp(&r.seq));
            seqs.push(r.seq);
        }
        // seqs is [c0 rc, c0 fwd, c1 rc, ...]; swap into [fwd, rc] order.
        for c in seqs.chunks_mut(2) {
            c.swap(0, 1);
        }

        let mut raw: Vec<(u64, u32, u32)> = Vec::new();
        for (oi, s) in seqs.iter().enumerate() {
            for (pos, v) in forward_seeds(s, &seed) {
                if complex_enough(v, weight) {
                    raw.push((v, oi as u32, pos as u32));
                }
            }
        }
        raw.sort_unstable();
        let n_values = 1usize << (2 * weight);
        let mut table = vec![0u32; n_values + 1];
        for &(v, _, _) in &raw {
            table[v as usize + 1] += 1;
        }
        for i in 1..=n_values {
            table[i] += table[i - 1];
        }
        let entries = raw.into_iter().map(|(_, oi, pos)| (oi, pos)).collect();
        Ok(Self {
            names,
            seqs,
            span: seed.span(),
            seed,
            table,
            entries,
            min_score,
            band,
            xdrop,
            single_hit: false,
            gate_sum: 4,
            gate_side: 6,
            thresholds: Vec::new(),
            logistic: None,
        })
    }

    /// Load a learned model. Either per-consensus minimum scores
    /// (`name<TAB>min_score`) or a logistic model (`#logistic`, `#mu`, `#sd`
    /// header lines followed by `name<TAB>offset`). Consensi absent from a
    /// threshold file keep `min_score`.
    pub fn load_thresholds(&mut self, path: &Path) -> io::Result<usize> {
        use std::collections::HashMap;
        let text = std::fs::read_to_string(path)?;
        if text.starts_with("#logistic") {
            let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());
            let mut lines = text.lines();
            let nums = |l: &str| -> Vec<f64> { l.split('\t').skip(1).filter_map(|v| v.trim().parse().ok()).collect() };
            let head = nums(lines.next().ok_or_else(|| bad("empty model"))?);
            if head.len() != 10 { return Err(bad("#logistic needs 8 weights, b0, tau")); }
            let mu = nums(lines.next().ok_or_else(|| bad("missing #mu"))?);
            let sd = nums(lines.next().ok_or_else(|| bad("missing #sd"))?);
            if mu.len() != 8 || sd.len() != 8 { return Err(bad("#mu/#sd need 8 values")); }
            let by_name: HashMap<&str, f64> = lines
                .filter_map(|l| { let mut f = l.split('\t'); Some((f.next()?, f.next()?.trim().parse().ok()?)) })
                .collect();
            let mut w = [0.0; 8]; w.copy_from_slice(&head[..8]);
            let mut m = [0.0; 8]; m.copy_from_slice(&mu);
            let mut d = [0.0; 8]; d.copy_from_slice(&sd);
            self.logistic = Some(Logistic {
                w, b0: head[8], tau: head[9], mu: m, sd: d,
                offsets: self.names.iter().map(|n| *by_name.get(n.as_str()).unwrap_or(&0.0)).collect(),
            });
            return Ok(by_name.len());
        }
        let by_name: HashMap<&str, i32> = text
            .lines()
            .filter_map(|l| {
                let mut f = l.split('\t');
                let name = f.next()?;
                let t = f.next()?.trim().parse().ok()?;
                Some((name, t))
            })
            .collect();
        self.thresholds = self
            .names
            .iter()
            .map(|n| *by_name.get(n.as_str()).unwrap_or(&self.min_score))
            .collect();
        Ok(by_name.len())
    }

    #[inline]
    fn accept(&self, cid: usize, total: i32) -> bool {
        total >= self.min_score && self.thresholds.get(cid).is_none_or(|&t| total >= t)
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn total_bases(&self) -> usize {
        self.seqs.iter().map(|s| s.len()).sum::<usize>() / 2
    }

    pub fn seeded_positions(&self) -> usize {
        self.entries.len()
    }

    /// Mask one record. Returns merged regions named by consensus and score.
    pub fn mask_record(&self, chrom: &str, seq: &[u8]) -> Vec<MaskedRegion> {
        self.mask_record_with_candidates(chrom, seq).0
    }

    /// Mask one record and also return every accepted hit with its features
    /// (unmerged), for training or applying a learned filter.
    pub fn mask_record_with_candidates(&self, chrom: &str, seq: &[u8]) -> (Vec<MaskedRegion>, Vec<Candidate>) {
        let n_starts = seq.len().saturating_sub(self.span - 1);
        if n_starts == 0 {
            return (Vec::new(), Vec::new());
        }
        let block = n_starts.div_ceil(rayon::current_num_threads().max(1)).max(1 << 20);
        let mut bounds = Vec::new();
        let mut b = 0usize;
        while b < n_starts {
            bounds.push((b, (b + block).min(n_starts)));
            b += block;
        }
        let chrom_id = fasta::seq_id(chrom).to_string();
        let parts: Vec<Vec<Candidate>> = bounds
            .into_par_iter()
            .map(|(bs, be)| self.scan_block(seq, bs, be))
            .collect();
        let cands: Vec<Candidate> = parts.into_iter().flatten().collect();
        let mut all: Vec<MaskedRegion> = cands
            .iter()
            .map(|c| MaskedRegion {
                chrom: chrom_id.clone(),
                start: c.start,
                end: c.end,
                name: format!("{}:{}", self.names[c.cid], c.score),
            })
            .collect();
        merge_regions(&mut all);
        (all, cands)
    }

    fn scan_block(&self, seq: &[u8], bs: usize, be: usize) -> Vec<Candidate> {
        let span = self.span;
        let weight = self.seed.weight();
        let mask = if span == 32 { u64::MAX } else { (1u64 << (2 * span)) - 1 };
        // Warm the rolling window so the first emitted seed starts at `bs`.
        let warm = bs.saturating_sub(span - 1);
        let mut regions: Vec<Candidate> = Vec::new();
        // Chaining state: the last few hits per oriented consensus, as
        // (genome pos, diagonal, consensus pos). A real copy hits the same
        // consensus on the same diagonal again and again, so remembering two
        // hits per consensus is enough — and it makes the lookup O(1) instead
        // of a scan over every recent hit in the window, which at low seed
        // weight is thousands of comparisons per genome position.
        let mut last: Vec<[(u32, i64, u32); SLOTS]> =
            vec![[(u32::MAX, 0, 0); SLOTS]; self.seqs.len()];
        // Single-hit mode: last failed (genome pos, diagonal) per consensus.
        let mut failed: Vec<(u32, i64)> = vec![(u32::MAX, 0); self.seqs.len()];
        let mut cands: Vec<(i32, u32, u32)> = Vec::new();
        let mut masked_until = 0usize;
        let (mut fwd, mut valid) = (0u64, 0usize);

        for i in warm..(be + span - 1).min(seq.len()) {
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
            if valid < span {
                continue;
            }
            let g = i + 1 - span;
            if g < bs || g >= be || g < masked_until {
                continue;
            }
            let v = self.seed.extract(fwd);
            let lo = self.table[v as usize] as usize;
            let hi = self.table[v as usize + 1] as usize;
            if lo == hi || !complex_enough(v, weight) {
                continue;
            }

            if self.single_hit {
                // Gate every entry cheaply, then spend banded DPs only on the
                // best-gated few. Inside an old L1 fragment a position can
                // pass the gate for a dozen related L1 consensi; running a DP
                // for each is where naive single-hit mode spends its time,
                // while the true consensus is almost always the top gate score.
                cands.clear();
                for &(oi, cpos) in &self.entries[lo..hi] {
                    let diag = g as i64 - cpos as i64;
                    let (fg, fd) = failed[oi as usize];
                    if fg != u32::MAX && (fg as usize) + CHAIN_WINDOW >= g && (fd - diag).abs() <= CHAIN_DIAG {
                        continue;
                    }
                    if let Some(gs) = self.gate(seq, oi as usize, cpos as usize, g) {
                        cands.push((gs, oi, cpos));
                    }
                }
                if cands.len() > SINGLE_HIT_TOP {
                    cands.select_nth_unstable_by(SINGLE_HIT_TOP - 1, |a, b| b.0.cmp(&a.0));
                    cands.truncate(SINGLE_HIT_TOP);
                }
                cands.sort_unstable_by(|a, b| b.0.cmp(&a.0));
                for &(_, oi, cpos) in cands.iter() {
                    if let Some(c) = self.extend_gated(seq, oi as usize, cpos as usize, g) {
                        masked_until = c.end;
                        regions.push(c);
                        break;
                    }
                    failed[oi as usize] = (g as u32, g as i64 - cpos as i64);
                }
                continue;
            }

            for &(oi, cpos) in &self.entries[lo..hi] {
                let diag = g as i64 - cpos as i64;
                let slots = &mut last[oi as usize];
                let anchor = slots.iter().position(|&(g0, d0, _)| {
                    g0 != u32::MAX
                        && (g0 as usize) < g
                        && g0 as usize + CHAIN_WINDOW >= g
                        && (d0 - diag).abs() <= CHAIN_DIAG
                });
                if let Some(ai) = anchor {
                    let (g0, _, cpos0) = slots[ai];
                    if let Some(c) = self.extend_hit(seq, oi as usize, cpos0 as usize, g0 as usize) {
                        masked_until = c.end;
                        regions.push(c);
                        break;
                    }
                    // Failed pair: forget it so the same spurious diagonal does
                    // not re-trigger at every following position.
                    slots[ai] = (u32::MAX, 0, 0);
                    continue;
                }
                // Remember this hit, replacing the oldest slot.
                let oldest = (0..SLOTS)
                    .min_by_key(|&k| if slots[k].0 == u32::MAX { 0 } else { slots[k].0 })
                    .unwrap();
                slots[oldest] = (g as u32, diag, cpos);
            }
        }
        regions
    }

    /// Ungapped look along the anchor diagonal: the cheap test that decides
    /// whether a banded DP is worth running. A spurious chain is random
    /// sequence here (expected -0.5/base); a real copy even at 62% identity
    /// averages +0.25/base. Either side may be broken by an indel, so it
    /// accepts on the sum or on one side alone. Returns the summed score when
    /// the gate passes.
    #[inline]
    fn gate(&self, seq: &[u8], oi: usize, cpos: usize, gpos: usize) -> Option<i32> {
        let cons = &self.seqs[oi];
        // Both windows sit *outside* the seed span. The seed's care positions
        // match by construction, so a window that included them would start
        // ~+6 in the black and the one-sided rule would fire by chance at a
        // large fraction of positions — which is exactly what made naive
        // single-hit mode 10x slower than it should be.
        let span = self.span;
        let left = cpos.min(gpos).min(32);
        let (cr, gr) = (cpos + span, gpos + span);
        let right = cons.len().saturating_sub(cr).min(seq.len().saturating_sub(gr)).min(32);
        let (mut sl, mut sr) = (0i32, 0i32);
        for o in 0..left {
            sl += if cons[cpos - 1 - o] == seq[gpos - 1 - o] { MATCH } else { MISMATCH };
        }
        for o in 0..right {
            sr += if cons[cr + o] == seq[gr + o] { MATCH } else { MISMATCH };
        }
        // Single-hit mode has no second seed hit backing the candidate, so
        // its gate is stricter.
        let pass = if self.single_hit {
            sl + sr >= self.gate_sum || sl >= self.gate_side || sr >= self.gate_side
        } else {
            sl + sr >= -8
        };
        if pass { Some(sl + sr) } else { None }
    }

    /// The two flank scores, for feature output.
    #[inline]
    fn gate_sides(&self, seq: &[u8], oi: usize, cpos: usize, gpos: usize) -> (i32, i32) {
        let cons = &self.seqs[oi];
        let span = self.span;
        let left = cpos.min(gpos).min(32);
        let (cr, gr) = (cpos + span, gpos + span);
        let right = cons.len().saturating_sub(cr).min(seq.len().saturating_sub(gr)).min(32);
        let (mut sl, mut sr) = (0i32, 0i32);
        for o in 0..left {
            sl += if cons[cpos - 1 - o] == seq[gpos - 1 - o] { MATCH } else { MISMATCH };
        }
        for o in 0..right {
            sr += if cons[cr + o] == seq[gr + o] { MATCH } else { MISMATCH };
        }
        (sl, sr)
    }

    /// Extend from an anchored seed in both directions (gate first).
    fn extend_hit(&self, seq: &[u8], oi: usize, cpos: usize, gpos: usize) -> Option<Candidate> {
        self.gate(seq, oi, cpos, gpos)?;
        self.extend_gated(seq, oi, cpos, gpos)
    }

    /// Banded X-drop extension both ways from an anchor that already passed
    /// the gate. Returns the hit if it clears the score floor (and the
    /// learned per-consensus threshold, if any).
    fn extend_gated(&self, seq: &[u8], oi: usize, cpos: usize, gpos: usize) -> Option<Candidate> {
        let cons = &self.seqs[oi];
        let band = self.band;
        // Forward: consensus tail vs a bounded genome slice.
        let c_fwd = &cons[cpos..];
        let g_end = (gpos + c_fwd.len() + band + 64).min(seq.len());
        let (sf, cf, gf) = xdrop_extend(c_fwd, &seq[gpos..g_end], band, self.xdrop, false);
        // Backward: the same DP walked from the ends of the prefixes, so no
        // reversed copies are allocated per extension.
        let c_bwd = &cons[..cpos];
        let g_start = gpos.saturating_sub(c_bwd.len() + band + 64);
        let (sb, cb, gb) = xdrop_extend(c_bwd, &seq[g_start..gpos], band, self.xdrop, true);
        let total = sf + sb;
        if self.accept(oi / 2, total) && gf + gb >= 30 {
            let (gl, gr) = self.gate_sides(seq, oi, cpos, gpos);
            if let Some(m) = &self.logistic {
                let (start, end) = (gpos - gb, gpos + gf);
                let len = (end - start) as f64;
                let gc = seq[start..end].iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 / len;
                let x = [
                    total as f64,
                    len.ln(),
                    total as f64 / len,
                    (cf + cb) as f64 / cons.len().max(1) as f64,
                    gc,
                    (gl + gr) as f64,
                    gl.min(gr) as f64,
                    (cons.len() as f64).ln(),
                ];
                if m.prob(oi / 2, x) < m.tau {
                    return None;
                }
            }
            Some(Candidate {
                start: gpos - gb,
                end: gpos + gf,
                cid: oi / 2,
                strand: oi % 2 == 0,
                score: total,
                score_fwd: sf,
                score_bwd: sb,
                cons_fwd: cf,
                cons_bwd: cb,
                cons_len: cons.len(),
                gate_left: gl,
                gate_right: gr,
            })
        } else {
            None
        }
    }
}

/// Banded X-drop extension of `a` against `b` from their starts (or, with
/// `rev`, from their ends walking backwards). Returns the best score and the
/// (a, b) lengths consumed at that score.
fn xdrop_extend(a: &[u8], b: &[u8], band: usize, xdrop: i32, rev: bool) -> (i32, usize, usize) {
    let (an, bn) = (a.len(), b.len());
    let at = |i: usize| if rev { a[an - i] } else { a[i - 1] };
    let bt = |j: usize| if rev { b[bn - j] } else { b[j - 1] };
    let width = 2 * band + 1;
    let mut prev = vec![NEG_INF; width];
    let mut cur = vec![NEG_INF; width];
    // Row 0: consuming j bases of b with nothing from a costs j gaps.
    for d in 0..width {
        let j = d as i64 - band as i64;
        if j >= 0 && (j as usize) <= b.len() {
            prev[d] = -(j as i32) * GAP;
        }
    }
    let (mut best, mut best_i, mut best_j) = (0i32, 0usize, 0usize);
    for i in 1..=a.len() {
        let mut row_best = NEG_INF;
        for d in 0..width {
            let j = i as i64 - band as i64 + d as i64;
            if j < 0 || j as usize > b.len() {
                cur[d] = NEG_INF;
                continue;
            }
            let j = j as usize;
            let mut v = NEG_INF;
            if j >= 1 && prev[d] > NEG_INF {
                let s = if at(i) == bt(j) { MATCH } else { MISMATCH };
                v = prev[d] + s;
            }
            if d + 1 < width && prev[d + 1] > NEG_INF {
                v = v.max(prev[d + 1] - GAP);
            }
            if d >= 1 && cur[d - 1] > NEG_INF {
                v = v.max(cur[d - 1] - GAP);
            }
            cur[d] = v;
            if v > row_best {
                row_best = v;
            }
            if v > best {
                best = v;
                best_i = i;
                best_j = j;
            }
        }
        if row_best < best - xdrop {
            break;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    (best, best_i, best_j)
}

/// Sort by start and merge overlapping or touching regions in place. A merged
/// region keeps the first non-empty name.
pub fn merge_regions(regions: &mut Vec<MaskedRegion>) {
    if regions.is_empty() {
        return;
    }
    regions.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    let mut out: Vec<MaskedRegion> = Vec::with_capacity(regions.len());
    for r in regions.drain(..) {
        match out.last_mut() {
            Some(last) if last.chrom == r.chrom && r.start <= last.end => {
                last.end = last.end.max(r.end);
                if last.name.is_empty() {
                    last.name = r.name;
                }
            }
            _ => out.push(r),
        }
    }
    *regions = out;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::splitmix64;
    use std::fs::File;
    use std::io::Write;

    fn rnd(salt: u64, n: usize) -> Vec<u8> {
        (0..n as u64)
            .map(|i| b"ACGT"[(splitmix64(i ^ salt) % 4) as usize])
            .collect()
    }
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
    fn xdrop_scores_identity_and_stops_in_noise() {
        let a = rnd(1, 200);
        let (s, i, j) = xdrop_extend(&a, &a, 8, 30, false);
        assert_eq!((s, i, j), (200, 200, 200));
        let (s, i, j) = xdrop_extend(&a, &a, 8, 30, true);
        assert_eq!((s, i, j), (200, 200, 200));
        let b = rnd(2, 200);
        let (s, _, _) = xdrop_extend(&a, &b, 8, 30, false);
        assert!(s < 15, "random score {}", s);
        // Reverse walk must agree with the forward walk on reversed inputs.
        let c = mutate(&a, 0.2, 9);
        let ar: Vec<u8> = a.iter().rev().copied().collect();
        let cr: Vec<u8> = c.iter().rev().copied().collect();
        assert_eq!(xdrop_extend(&a, &c, 8, 30, true), xdrop_extend(&ar, &cr, 8, 30, false));
    }

    #[test]
    fn library_finds_diverged_copies_on_both_strands() {
        let dir = std::env::temp_dir().join("krep_align_test");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let cons = rnd(0x77, 300);
        let lib_path = dir.join("lib.fa");
        {
            let mut f = File::create(&lib_path).unwrap();
            writeln!(f, ">fam1").unwrap();
            f.write_all(&cons).unwrap();
            writeln!(f).unwrap();
        }
        // Genome: 6 copies at 25% divergence, alternating strands, in unique DNA.
        let mut genome = rnd(0xa1, 1000);
        let mut truth = Vec::new();
        for i in 0..6u64 {
            let copy = mutate(&cons, 0.25, 0x500 + i);
            let copy = if i % 2 == 1 { revcomp(&copy) } else { copy };
            truth.push((genome.len(), genome.len() + copy.len()));
            genome.extend_from_slice(&copy);
            genome.extend_from_slice(&rnd(0xb00 + i, 800));
        }
        let lib = Library::load(&lib_path, 50, 16, 40).unwrap();
        let regions = lib.mask_record("chrT", &genome);
        for (s, e) in &truth {
            let cov: usize = regions
                .iter()
                .map(|r| r.end.min(*e).saturating_sub(r.start.max(*s)))
                .sum();
            assert!(cov * 100 >= (e - s) * 80, "copy {}..{} covered {}", s, e, cov);
        }
        let total: usize = regions.iter().map(|r| r.end - r.start).sum();
        let truth_total: usize = truth.iter().map(|(s, e)| e - s).sum();
        assert!(total <= truth_total + 6 * 60, "over-masked: {} vs {}", total, truth_total);
        let _ = std::fs::remove_dir_all(&dir);
    }
}

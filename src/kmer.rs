//! 2-bit DNA k-mer encoding and iteration.
//!
//! Bases are encoded as: A=00, C=01, G=10, T=11.
//! K-mers are stored in the lower 2*k bits of a `u64`, which works for k ≤ 32.

/// Encode a single ASCII base as a 2-bit value.
/// Returns `None` for any base other than A/C/G/T (case-insensitive).
pub fn encode_base(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Decode a 2-bit value back to an ASCII base.
#[allow(dead_code)]
pub fn decode_base(bits: u64) -> u8 {
    match bits & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!(),
    }
}

/// Reverse complement of a k-mer stored in the lower 2*k bits of `kmer`.
pub fn reverse_complement(kmer: u64, k: usize) -> u64 {
    let mut rev = 0u64;
    let mut x = kmer;
    for _ in 0..k {
        rev = (rev << 2) | (3 - (x & 0b11));
        x >>= 2;
    }
    rev
}

/// Canonical form of a k-mer: the smaller of the forward and reverse-complement values.
pub fn canonical(kmer: u64, k: usize) -> u64 {
    let rc = reverse_complement(kmer, k);
    kmer.min(rc)
}

/// Convert a k-mer value into a human-readable string.
#[allow(dead_code)]
pub fn to_string(kmer: u64, k: usize) -> String {
    let mut s = String::with_capacity(k);
    for i in (0..k).rev() {
        s.push(decode_base(kmer >> (2 * i)) as char);
    }
    s
}

/// Generate all canonical k-mers at Hamming distance ≤ `max_mismatches` from
/// `center` (excluding the center itself if `exclude_self` is true).
///
/// This is useful for expanding a set of high-count k-mers to also cover
/// mutated versions, improving recall on diverged repeats.
pub fn neighbors(center: u64, k: usize, max_mismatches: usize, exclude_self: bool) -> Vec<u64> {
    if max_mismatches == 0 {
        if exclude_self {
            return Vec::new();
        }
        return vec![center];
    }

    let mut result = Vec::new();
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    // Depth-first expansion of mismatches.
    // Stack entries: (current value, position being considered next, remaining mismatches).
    let mut stack = vec![(center, 0usize, max_mismatches)];
    while let Some((current, pos, remaining)) = stack.pop() {
        if pos >= k {
            if exclude_self && current == center {
                continue;
            }
            result.push(current);
            continue;
        }
        if remaining == 0 {
            // No more mismatches allowed; copy remaining positions verbatim.
            stack.push((current, k, 0));
            continue;
        }

        let shift = 2 * pos;
        let current_base = (current >> shift) & 0b11;
        for b in 0..4u64 {
            let mut next = current;
            if b != current_base {
                next = (next & !(0b11 << shift)) | (b << shift);
            }
            let used = if b == current_base { 0 } else { 1 };
            stack.push((next, pos + 1, remaining - used));
        }
    }

    // Canonicalize all neighbors and deduplicate.
    result.sort_unstable();
    result.dedup();
    result
        .into_iter()
        .map(|kmer| canonical(kmer & mask, k))
        .collect()
}

/// Iterator over canonical k-mers of a DNA sequence.
/// Emits `(position, canonical_kmer)` for every fully-ACGT window.
pub struct KmerIter<'a> {
    seq: &'a [u8],
    k: usize,
    mask: u64,
    top_shift: u32, // 2 * (k - 1), where the complement base enters `rev`
    pos: usize,
    fwd: u64,
    rev: u64, // reverse complement of the current window, kept incrementally
    valid: usize, // number of consecutive valid bases ending at pos-1
}

impl<'a> KmerIter<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!(k > 0 && k <= 32, "k must be in 1..=32");
        let mask = if k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        Self {
            seq,
            k,
            mask,
            top_shift: 2 * (k as u32 - 1),
            pos: 0,
            fwd: 0,
            rev: 0,
            valid: 0,
        }
    }
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let base = self.seq[self.pos];
            self.pos += 1;

            if let Some(bits) = encode_base(base) {
                // Roll both strands forward. Appending a base on the right of
                // `fwd` prepends its complement on the left of `rev`, so the
                // reverse complement never has to be recomputed from scratch.
                self.fwd = ((self.fwd << 2) | bits) & self.mask;
                self.rev = (self.rev >> 2) | ((3 - bits) << self.top_shift);
                self.valid += 1;
                if self.valid >= self.k {
                    let start = self.pos - self.k;
                    return Some((start, self.fwd.min(self.rev)));
                }
            } else {
                self.fwd = 0;
                self.rev = 0;
                self.valid = 0;
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_roundtrip() {
        assert_eq!(encode_base(b'a'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'g'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b'N'), None);
    }

    #[test]
    fn reverse_complement_small() {
        // k=3, "ACG" = 0b00_01_10 = 6
        let rc = reverse_complement(0b00_01_10, 3);
        // "CGT" = 0b01_10_11 = 27
        assert_eq!(rc, 0b01_10_11);
    }

    #[test]
    fn canonical_palindrome() {
        // "AT" reverse complement is "AT"
        let kmer = 0b00_11; // AT
        assert_eq!(canonical(kmer, 2), kmer);
    }

    #[test]
    fn iterator_skips_n() {
        let seq = b"AACGTNACGTA";
        let kmers: Vec<_> = KmerIter::new(seq, 3).collect();
        // First window AAC, ACG, CGT, then N breaks, then GTA starts after N+2
        assert!(kmers.iter().any(|(_, k)| *k == canonical(0b00_00_01, 3))); // AAC
        assert!(kmers.iter().any(|(_, k)| *k == canonical(0b00_01_10, 3))); // ACG
        assert_eq!(
            kmers.last().map(|(_, k)| *k),
            Some(canonical(0b10_11_00, 3))
        ); // GTA
    }

    #[test]
    fn to_string_matches() {
        let kmer = 0b00_01_10; // ACG
        assert_eq!(to_string(kmer, 3), "ACG");
    }

    #[test]
    fn neighbors_hamming_1() {
        // k=2, center = AA = 0b0000.
        let n = neighbors(0, 2, 1, true);
        // Hamming-1 neighbors of AA are AC, AG, AT, CA, GA, TA (canonical).
        // (CC, GG, TT are 2 substitutions away from AA).
        assert_eq!(n.len(), 6);
    }

    #[test]
    fn neighbors_include_self() {
        let n = neighbors(0, 2, 0, false);
        assert_eq!(n, vec![0]);
    }
}

#[cfg(test)]
mod rolling_tests {
    use super::*;

    #[test]
    fn rolling_canonical_matches_reference() {
        let seq = b"ACGTTGCAANNACGTACGTTTGACCAGTNACGTAGCATCGATCGGATCCAA";
        for k in [3usize, 7, 15, 18, 21, 31] {
            let got: Vec<(usize, u64)> = KmerIter::new(seq, k).collect();
            // Reference: recompute the canonical form from the raw window.
            let mut want = Vec::new();
            for start in 0..seq.len().saturating_sub(k - 1) {
                let win = &seq[start..start + k];
                let mut fwd = 0u64;
                let mut ok = true;
                for &b in win {
                    match encode_base(b) {
                        Some(bits) => fwd = (fwd << 2) | bits,
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    want.push((start, canonical(fwd, k)));
                }
            }
            assert_eq!(got, want, "mismatch at k={}", k);
        }
    }
}

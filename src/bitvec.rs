//! Compact, zero-overhead bit vector.
//!
//! Stores one bit per entry in `u64` blocks for dense memory and fast
//! set/clear/any operations. The API is intentionally small and used by the
//! masker for coverage tracking.

pub struct BitVec {
    blocks: Vec<u64>,
    len: usize,
}

impl BitVec {
    /// Create a `BitVec` with `len` bits, all cleared.
    pub fn new(len: usize) -> Self {
        let n_blocks = len.div_ceil(64);
        Self {
            blocks: vec![0u64; n_blocks],
            len,
        }
    }

    /// Number of bits in the vector.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.len
    }

    #[must_use]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Set bit `i` to 1.
    pub fn set(&mut self, i: usize) {
        debug_assert!(i < self.len, "index {} out of bounds {}", i, self.len);
        self.blocks[i >> 6] |= 1u64 << (i & 63);
    }

    /// Clear bit `i` to 0.
    #[allow(dead_code)]
    pub fn clear(&mut self, i: usize) {
        debug_assert!(i < self.len, "index {} out of bounds {}", i, self.len);
        self.blocks[i >> 6] &= !(1u64 << (i & 63));
    }

    /// True if bit `i` is set.
    #[inline]
    pub fn get(&self, i: usize) -> bool {
        debug_assert!(i < self.len, "index {} out of bounds {}", i, self.len);
        (self.blocks[i >> 6] & (1u64 << (i & 63))) != 0
    }

    /// Clear all bits.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.blocks.fill(0);
    }

    /// Count the number of set bits.
    #[allow(dead_code)]
    pub fn count_ones(&self) -> usize {
        self.blocks.iter().map(|b| b.count_ones() as usize).sum()
    }

    /// Iterator over contiguous runs of set bits. Each run is returned as
    /// `(start, end)` where `end` is exclusive.
    pub fn set_runs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        SetRuns {
            bitvec: self,
            pos: 0,
        }
    }
}

struct SetRuns<'a> {
    bitvec: &'a BitVec,
    pos: usize,
}

impl<'a> Iterator for SetRuns<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        // Skip zeros.
        while self.pos < self.bitvec.len && !self.bitvec.get(self.pos) {
            self.pos += 1;
        }
        if self.pos >= self.bitvec.len {
            return None;
        }
        let start = self.pos;
        while self.pos < self.bitvec.len && self.bitvec.get(self.pos) {
            self.pos += 1;
        }
        Some((start, self.pos))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_and_get() {
        let mut b = BitVec::new(100);
        assert!(!b.get(7));
        b.set(7);
        assert!(b.get(7));
        assert!(!b.get(6));
        assert!(!b.get(8));
    }

    #[test]
    fn count_and_reset() {
        let mut b = BitVec::new(200);
        for i in [0, 1, 63, 64, 65, 199].iter() {
            b.set(*i);
        }
        assert_eq!(b.count_ones(), 6);
        b.reset();
        assert_eq!(b.count_ones(), 0);
    }

    #[test]
    fn set_runs() {
        let mut b = BitVec::new(20);
        b.set(1);
        b.set(2);
        b.set(3);
        b.set(10);
        b.set(18);
        b.set(19);
        let runs: Vec<_> = b.set_runs().collect();
        assert_eq!(runs, vec![(1, 4), (10, 11), (18, 20)]);
    }
}

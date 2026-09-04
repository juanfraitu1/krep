//! Counting Bloom filter with u8 counters and double hashing.

use crate::rng::Rng;

/// Deterministic 64-bit hasher based on splitmix64.
/// Produces two independent-ish hashes from a single `u64` key.
#[derive(Debug, Clone)]
pub struct SplitMix64Hasher {
    seed1: u64,
    seed2: u64,
}

impl SplitMix64Hasher {
    pub fn new(seed: u64) -> Self {
        // Use two different seeds derived from the caller's seed.
        let mut rng = Rng::new(seed);
        Self {
            seed1: rng.next_u64(),
            seed2: rng.next_u64(),
        }
    }

    fn splitmix64(state: u64) -> u64 {
        let mut z = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Produce (h1, h2) for double hashing.
    pub fn hash(&self, key: u64) -> (u64, u64) {
        let h1 = Self::splitmix64(key.wrapping_add(self.seed1));
        let h2 = Self::splitmix64(key.wrapping_add(self.seed2));
        (h1, h2)
    }
}

/// Counting Bloom filter.
///
/// Uses `k` hash functions derived from a single 64-bit key via double hashing:
/// `h_i(key) = (h1 + i * h2) mod m`.
pub struct CountingBloomFilter {
    counters: Vec<u8>,
    m: usize,
    k: usize,
    hasher: SplitMix64Hasher,
}

#[allow(dead_code)]
impl CountingBloomFilter {
    /// Build a new CBF sized for roughly `expected_inserts` items.
    /// `k` is the number of hash functions; `factor` is the multiplier for `m`
    /// (default 8 gives a load factor of ~1/8 when full).
    pub fn new(expected_inserts: usize, k: usize, hash_seed: u64) -> Self {
        Self::with_factor(expected_inserts, k, hash_seed, 8)
    }

    /// Build a new CBF with a configurable size multiplier.
    /// Smaller factors save memory but increase hash collisions.
    pub fn with_factor(expected_inserts: usize, k: usize, hash_seed: u64, factor: usize) -> Self {
        assert!(k > 0, "k must be > 0");
        let m = (expected_inserts.saturating_mul(factor)).max(1024);
        Self {
            counters: vec![0; m],
            m,
            k,
            hasher: SplitMix64Hasher::new(hash_seed),
        }
    }

    /// Return the number of slots.
    pub fn capacity(&self) -> usize {
        self.m
    }

    /// Return the number of hash functions.
    pub fn num_hashes(&self) -> usize {
        self.k
    }

    /// Saturating insert of a key.
    pub fn insert(&mut self, key: u64) {
        let (h1, h2) = self.hasher.hash(key);
        for i in 0..self.k as u64 {
            let idx = ((h1.wrapping_add(i.wrapping_mul(h2))) as usize) % self.m;
            self.counters[idx] = self.counters[idx].saturating_add(1);
        }
    }

    /// Estimated count of a key (minimum over the `k` counters).
    pub fn count(&self, key: u64) -> u8 {
        let (h1, h2) = self.hasher.hash(key);
        let mut min = u8::MAX;
        for i in 0..self.k as u64 {
            let idx = ((h1.wrapping_add(i.wrapping_mul(h2))) as usize) % self.m;
            let c = self.counters[idx];
            if c < min {
                min = c;
            }
        }
        min
    }

    /// True if any counter for the key is non-zero.
    pub fn contains(&self, key: u64) -> bool {
        self.count(key) > 0
    }

    /// Fraction of non-zero counters.
    pub fn load(&self) -> f64 {
        let occupied = self.counters.iter().filter(|c| **c > 0).count();
        occupied as f64 / self.m as f64
    }

    /// Reset all counters to zero.
    pub fn clear(&mut self) {
        self.counters.fill(0);
    }

    /// Merge another CBF into this one by saturating addition of counters.
    /// Both filters must have been created with the same `m` and `k`.
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.m, other.m, "CBF merge requires equal capacity");
        assert_eq!(self.k, other.k, "CBF merge requires equal k");
        for (dst, src) in self.counters.iter_mut().zip(other.counters.iter()) {
            *dst = dst.saturating_add(*src);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_and_count() {
        let mut cbf = CountingBloomFilter::new(100, 4, 1);
        for i in 0u64..10 {
            cbf.insert(i);
        }
        for i in 0u64..10 {
            assert!(cbf.count(i) >= 1);
        }
    }

    #[test]
    fn count_is_min() {
        let mut cbf = CountingBloomFilter::new(100, 4, 2);
        let key = 0xdead_beef_u64;
        for _ in 0..5 {
            cbf.insert(key);
        }
        assert_eq!(cbf.count(key), 5);
    }

    #[test]
    fn saturates_at_u8_max() {
        let mut cbf = CountingBloomFilter::new(10, 4, 3);
        let key = 123u64;
        for _ in 0..300 {
            cbf.insert(key);
        }
        assert_eq!(cbf.count(key), 255);
    }

    #[test]
    fn clear_works() {
        let mut cbf = CountingBloomFilter::new(100, 4, 4);
        cbf.insert(7);
        assert!(cbf.contains(7));
        cbf.clear();
        assert!(!cbf.contains(7));
    }
}

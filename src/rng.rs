//! Seeded splitmix64 pseudo-random number generator.
//!
//! Uses the standard 64-bit splitmix64 algorithm so the PRNG is small,
//! deterministic, and has zero dependencies.

#[derive(Debug, Clone, Copy)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Create a new generator from a 64-bit seed.
    pub fn new(seed: u64) -> Self {
        // Mix the seed once before use to avoid weak low seeds.
        let mut s = Self { state: seed };
        s.next_u64();
        s
    }

    /// Return the next 64-bit value and advance the state.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Return a float in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        // 53 bits of precision via the upper bits of a u64.
        ((self.next_u64() >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// Return an integer in [0, n).
    pub fn range_usize(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        ((self.next_u64() as u128 * n as u128) >> 64) as usize
    }

    /// Return true with probability `p`.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.next_f64() < p
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn range_usize_in_bounds() {
        let mut rng = Rng::new(123);
        for _ in 0..1000 {
            let v = rng.range_usize(10);
            assert!(v < 10);
        }
    }
}

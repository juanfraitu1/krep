//! Mock genome generator with known inserted repetitive elements.

use crate::rng::Rng;
use std::fmt::Write as _;

/// A single placed repeat instance.
#[derive(Debug, Clone)]
pub struct Repeat {
    pub start: usize,
    pub end: usize,
    pub family: String,
}

#[allow(dead_code)]
impl Repeat {
    pub fn len(&self) -> usize {
        self.end - self.start
    }
}

/// Generate a random DNA background of `size` bases with the given GC fraction.
pub fn generate_background(size: usize, gc: f64, rng: &mut Rng) -> Vec<u8> {
    let at = (1.0 - gc) / 2.0;
    let gc_each = gc / 2.0;
    let weights = [(b'A', at), (b'C', gc_each), (b'G', gc_each), (b'T', at)];
    (0..size).map(|_| sample_base(&weights, rng)).collect()
}

fn sample_base(weights: &[(u8, f64)], rng: &mut Rng) -> u8 {
    let r = rng.next_f64();
    let mut cum = 0.0;
    for (base, w) in weights {
        cum += w;
        if r < cum {
            return *base;
        }
    }
    weights.last().unwrap().0
}

/// Mutate a sequence by substituting each base with probability `divergence`.
/// The new base is always different from the original.
pub fn mutate(seq: &[u8], divergence: f64, rng: &mut Rng) -> Vec<u8> {
    let bases = *b"ACGT";
    seq.iter()
        .map(|base| {
            if rng.bernoulli(divergence) {
                let mut new_base = bases[rng.range_usize(4)];
                while new_base == *base {
                    new_base = bases[rng.range_usize(4)];
                }
                new_base
            } else {
                *base
            }
        })
        .collect()
}

/// Repeat family definition.
pub struct Family {
    pub name: String,
    pub min_len: usize,
    pub max_len: usize,
    pub gc: f64,
    pub copies: usize,
    pub kind: FamilyKind,
}

pub enum FamilyKind {
    /// Random sequence with the family's GC content.
    Random,
    /// Tandem repeat of a short motif repeated `motif` times.
    Tandem {
        motif: Vec<u8>,
        repeats: (usize, usize),
    },
}

/// Generate a consensus sequence for a family.
pub fn generate_consensus(family: &Family, rng: &mut Rng) -> Vec<u8> {
    match &family.kind {
        FamilyKind::Random => {
            let len = rng.range_usize(family.max_len - family.min_len + 1) + family.min_len;
            generate_background(len, family.gc, rng)
        }
        FamilyKind::Tandem { motif, repeats } => {
            let min_r = repeats.0;
            let max_r = repeats.1;
            let n = rng.range_usize(max_r - min_r + 1) + min_r;
            motif
                .iter()
                .copied()
                .cycle()
                .take(motif.len() * n)
                .collect()
        }
    }
}

/// Create the default set of repeat families used for benchmarking.
pub fn default_families() -> Vec<Family> {
    vec![
        Family {
            name: "ALU".into(),
            min_len: 280,
            max_len: 320,
            gc: 0.60,
            copies: 100,
            kind: FamilyKind::Random,
        },
        Family {
            name: "LINE1".into(),
            min_len: 5800,
            max_len: 6200,
            gc: 0.40,
            copies: 20,
            kind: FamilyKind::Random,
        },
        Family {
            name: "SINE".into(),
            min_len: 180,
            max_len: 220,
            gc: 0.50,
            copies: 50,
            kind: FamilyKind::Random,
        },
        Family {
            name: "LTR".into(),
            min_len: 1400,
            max_len: 1600,
            gc: 0.50,
            copies: 30,
            kind: FamilyKind::Random,
        },
        Family {
            name: "MICROSAT".into(),
            min_len: 40,
            max_len: 120,
            gc: 0.50,
            copies: 200,
            kind: FamilyKind::Tandem {
                motif: vec![b'C', b'A'],
                repeats: (20, 60),
            },
        },
        Family {
            name: "SAT".into(),
            min_len: 400,
            max_len: 600,
            gc: 0.50,
            copies: 50,
            kind: FamilyKind::Tandem {
                motif: vec![b'A', b'A', b'T', b'C', b'G'],
                repeats: (80, 120),
            },
        },
    ]
}

/// Create the default families plus an optional segmental-duplication family.
pub fn families_with_segdup() -> Vec<Family> {
    let mut families = default_families();
    families.push(Family {
        name: "SEG_DUP".into(),
        min_len: 5000,
        max_len: 10000,
        gc: 0.50,
        copies: 10,
        kind: FamilyKind::Random,
    });
    families
}

#[allow(dead_code)]
/// Build a mock genome.
///
/// Returns `(background_sequence, repeats)` where each repeat has a start, end,
/// and family name. The caller can then write the FASTA and BED files.
pub fn build_mock_genome(
    size: usize,
    gc: f64,
    families: &[Family],
    divergence: f64,
    rng: &mut Rng,
) -> (Vec<u8>, Vec<Repeat>) {
    build_mock_genome_with_segdup(size, gc, families, divergence, divergence, rng)
}

/// Build a mock genome, allowing segmental duplications to have a separate
/// (usually lower) divergence from transposable elements.
pub fn build_mock_genome_with_segdup(
    size: usize,
    gc: f64,
    families: &[Family],
    divergence: f64,
    segdup_divergence: f64,
    rng: &mut Rng,
) -> (Vec<u8>, Vec<Repeat>) {
    let mut genome = generate_background(size, gc, rng);

    // Pre-generate one consensus per family and then create mutated copies.
    let mut copies: Vec<(usize, String, Vec<u8>)> = Vec::new();
    for (fi, family) in families.iter().enumerate() {
        // Use a family-specific seed so the consensus is stable across runs with the same global seed.
        let mut family_rng = Rng::new(rng.next_u64().wrapping_add(fi as u64));
        let consensus = generate_consensus(family, &mut family_rng);
        let is_segdup = family.name == "SEG_DUP";
        let div = if is_segdup {
            segdup_divergence
        } else {
            divergence
        };
        for ci in 0..family.copies {
            let mut copy_rng = Rng::new(rng.next_u64().wrapping_add((fi * 10000 + ci) as u64));
            let seq = mutate(&consensus, div, &mut copy_rng);
            copies.push((seq.len(), family.name.clone(), seq));
        }
    }

    // Shuffle copies so placement order is random.
    for i in (1..copies.len()).rev() {
        let j = rng.range_usize(i + 1);
        copies.swap(i, j);
    }

    // Place copies without overlap using random positions and retries.
    let mut placed: Vec<Repeat> = Vec::with_capacity(copies.len());
    let mut occupied = vec![false; size];
    const MAX_RETRIES: usize = 1000;

    for (len, family, seq) in copies {
        let mut placed_ok = false;
        for _ in 0..MAX_RETRIES {
            if len > size {
                break;
            }
            let start = rng.range_usize(size - len + 1);
            let end = start + len;
            let overlaps = occupied[start..end].iter().any(|x| *x);
            if !overlaps {
                genome[start..end].copy_from_slice(&seq);
                occupied[start..end].fill(true);
                placed.push(Repeat {
                    start,
                    end,
                    family: family.clone(),
                });
                placed_ok = true;
                break;
            }
        }
        if !placed_ok {
            eprintln!("Warning: could not place {} bp {} repeat", len, family);
        }
    }

    placed.sort_by_key(|r| r.start);
    (genome, placed)
}

/// Format repeats as a simple BED-like string.
pub fn repeats_to_bed(chrom: &str, repeats: &[Repeat]) -> String {
    let mut s = String::new();
    for r in repeats {
        writeln!(s, "{}\t{}\t{}\t{}", chrom, r.start, r.end, r.family).unwrap();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn background_gc() {
        let mut rng = Rng::new(1);
        let seq = generate_background(100_000, 0.5, &mut rng);
        let gc =
            seq.iter().filter(|b| **b == b'C' || **b == b'G').count() as f64 / seq.len() as f64;
        assert!((gc - 0.5).abs() < 0.02);
    }

    #[test]
    fn mutation_rate() {
        let mut rng = Rng::new(2);
        let seq = generate_background(10_000, 0.5, &mut rng);
        let mutated = mutate(&seq, 0.10, &mut Rng::new(3));
        let diff = seq
            .iter()
            .zip(mutated.iter())
            .filter(|(a, b)| a != b)
            .count() as f64
            / seq.len() as f64;
        assert!((diff - 0.10).abs() < 0.02);
    }

    #[test]
    fn placed_repeats_no_overlap() {
        let mut rng = Rng::new(4);
        let families = default_families();
        let size = 1_000_000;
        let (_genome, repeats) = build_mock_genome(size, 0.5, &families, 0.05, &mut rng);

        let total_repeat: usize = repeats.iter().map(|r| r.len()).sum();
        let expected: usize = families
            .iter()
            .map(|f| {
                let avg_len = (f.min_len + f.max_len) / 2;
                avg_len * f.copies
            })
            .sum();
        // Allow for a few placement failures, but the vast majority should fit.
        assert!(
            total_repeat > expected * 9 / 10,
            "total repeat {} vs expected {}",
            total_repeat,
            expected
        );

        for window in repeats.windows(2) {
            assert!(window[0].end <= window[1].start, "overlapping repeats");
        }
    }
}

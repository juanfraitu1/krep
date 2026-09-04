//! BED overlap evaluation (precision, recall, F1, per-family recall).

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// A single interval, optionally tagged with a family name.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Interval {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub family: Option<String>,
}

#[allow(dead_code)]
impl Interval {
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}

/// Parse a BED file. Supports 3-column BED or 4-column BED (family/name in column 4).
pub fn read_bed<P: AsRef<Path>>(path: P) -> io::Result<Vec<Interval>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let chrom = cols[0].to_string();
        let start: usize = cols[1]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad start: {e}")))?;
        let end: usize = cols[2]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bad end: {e}")))?;
        let family = cols.get(3).map(|s| s.to_string());
        intervals.push(Interval {
            chrom,
            start,
            end,
            family,
        });
    }

    Ok(intervals)
}

/// Compute the total length covered by a set of intervals (merging overlaps within the set).
pub fn total_covered(intervals: &[Interval]) -> usize {
    if intervals.is_empty() {
        return 0;
    }
    // Group by chromosome first: coordinates from different chromosomes share
    // the same numeric range, so sweeping them on one axis would merge
    // unrelated intervals and undercount.
    let mut by_chrom: HashMap<&str, Vec<(usize, i32)>> = HashMap::new();
    for iv in intervals {
        let e = by_chrom.entry(iv.chrom.as_str()).or_default();
        e.push((iv.start, 1));
        e.push((iv.end, -1));
    }

    let mut covered = 0usize;
    for (_chrom, mut events) in by_chrom {
        events.sort_by_key(|e| e.0);
        let mut depth = 0i32;
        let mut prev = events[0].0;
        for (pos, delta) in events {
            if depth > 0 && pos > prev {
                covered += pos - prev;
            }
            depth += delta;
            prev = pos;
        }
    }
    covered
}

/// Intersection length between two sets of intervals, counting each base once
/// even when intervals overlap within a set.
fn intersection(a: &[Interval], b: &[Interval]) -> usize {
    if a.is_empty() || b.is_empty() {
        return 0;
    }

    // Sweep-line per chromosome with two counters: truth depth and predicted
    // depth. A base is in the intersection when both depths are > 0.
    let mut by_chrom: HashMap<&str, Vec<(usize, i32, i32)>> = HashMap::new();
    for iv in a {
        let e = by_chrom.entry(iv.chrom.as_str()).or_default();
        e.push((iv.start, 1, 0));
        e.push((iv.end, -1, 0));
    }
    for iv in b {
        // Only chromosomes present in the truth set can contribute.
        if let Some(e) = by_chrom.get_mut(iv.chrom.as_str()) {
            e.push((iv.start, 0, 1));
            e.push((iv.end, 0, -1));
        }
    }

    let mut total = 0usize;
    for (_chrom, mut events) in by_chrom {
        events.sort_by_key(|e| e.0);
        let mut truth_depth = 0i32;
        let mut pred_depth = 0i32;
        let mut prev = events[0].0;
        for (pos, d_truth, d_pred) in events {
            if pos > prev && truth_depth > 0 && pred_depth > 0 {
                total += pos - prev;
            }
            truth_depth += d_truth;
            pred_depth += d_pred;
            prev = pos;
        }
    }

    total
}

/// Aggregate evaluation metrics.
#[derive(Debug, Clone, Default)]
pub struct Metrics {
    pub truth_total: usize,
    pub pred_total: usize,
    pub intersection: usize,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
    /// Per truth family recall (family -> recall).
    pub per_family_recall: HashMap<String, f64>,
}

/// Evaluate predicted intervals against ground-truth intervals.
pub fn evaluate(truth: &[Interval], pred: &[Interval]) -> Metrics {
    let truth_total = total_covered(truth);
    let pred_total = total_covered(pred);
    let inter = intersection(truth, pred);

    let precision = if pred_total == 0 {
        0.0
    } else {
        inter as f64 / pred_total as f64
    };
    let recall = if truth_total == 0 {
        0.0
    } else {
        inter as f64 / truth_total as f64
    };
    let f1 = if precision + recall == 0.0 {
        0.0
    } else {
        2.0 * precision * recall / (precision + recall)
    };

    let mut per_family_recall = HashMap::new();
    let mut by_family: HashMap<String, Vec<Interval>> = HashMap::new();
    for iv in truth {
        if let Some(fam) = &iv.family {
            by_family.entry(fam.clone()).or_default().push(iv.clone());
        }
    }
    for (fam, fam_truth) in by_family {
        let fam_total = total_covered(&fam_truth);
        let fam_inter = intersection(&fam_truth, pred);
        let rec = if fam_total == 0 {
            0.0
        } else {
            fam_inter as f64 / fam_total as f64
        };
        per_family_recall.insert(fam, rec);
    }

    Metrics {
        truth_total,
        pred_total,
        intersection: inter,
        precision,
        recall,
        f1,
        per_family_recall,
    }
}

/// Pretty-print metrics.
pub fn format_metrics(m: &Metrics) -> String {
    let mut s = String::new();
    use std::fmt::Write as _;
    writeln!(s, "Truth bases:      {}", m.truth_total).unwrap();
    writeln!(s, "Predicted bases:  {}", m.pred_total).unwrap();
    writeln!(s, "Intersection:     {}", m.intersection).unwrap();
    writeln!(s, "Precision:        {:.4}", m.precision).unwrap();
    writeln!(s, "Recall:           {:.4}", m.recall).unwrap();
    writeln!(s, "F1:               {:.4}", m.f1).unwrap();
    if !m.per_family_recall.is_empty() {
        writeln!(s, "Per-family recall:").unwrap();
        let mut families: Vec<&String> = m.per_family_recall.keys().collect();
        families.sort();
        for fam in families {
            writeln!(s, "  {}: {:.4}", fam, m.per_family_recall[fam]).unwrap();
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: usize, end: usize, family: Option<&str>) -> Interval {
        Interval {
            chrom: chrom.into(),
            start,
            end,
            family: family.map(|s| s.into()),
        }
    }

    #[test]
    fn total_covered_merge() {
        let intervals = vec![
            iv("chr1", 0, 100, None),
            iv("chr1", 50, 150, None),
            iv("chr1", 200, 300, None),
        ];
        assert_eq!(total_covered(&intervals), 250);
    }

    #[test]
    fn intersection_basic() {
        let truth = vec![iv("chr1", 0, 100, None), iv("chr1", 200, 300, None)];
        let pred = vec![iv("chr1", 50, 150, None), iv("chr1", 250, 350, None)];
        let inter = intersection(&truth, &pred);
        assert_eq!(inter, 100); // 50..100 and 250..300
    }

    #[test]
    fn perfect_eval() {
        let truth = vec![iv("chr1", 10, 100, Some("ALU"))];
        let pred = truth.clone();
        let m = evaluate(&truth, &pred);
        assert!((m.precision - 1.0).abs() < 1e-9);
        assert!((m.recall - 1.0).abs() < 1e-9);
        assert!((m.f1 - 1.0).abs() < 1e-9);
    }
}

#[cfg(test)]
mod chrom_tests {
    use super::*;

    fn iv(chrom: &str, start: usize, end: usize) -> Interval {
        Interval {
            chrom: chrom.to_string(),
            start,
            end,
            family: None,
        }
    }

    #[test]
    fn covered_does_not_merge_across_chromosomes() {
        // Same coordinates on two chromosomes are 200 bp in total, not 100.
        let ivs = vec![iv("chr1", 0, 100), iv("chr2", 0, 100)];
        assert_eq!(total_covered(&ivs), 200);
    }

    #[test]
    fn intersection_is_chromosome_aware() {
        // Identical coordinates but different chromosomes must not intersect.
        let truth = vec![iv("chr1", 0, 100)];
        let pred = vec![iv("chr2", 0, 100)];
        assert_eq!(intersection(&truth, &pred), 0);

        let pred_same = vec![iv("chr1", 50, 150)];
        assert_eq!(intersection(&truth, &pred_same), 50);
    }

    #[test]
    fn overlaps_within_a_chromosome_count_once() {
        let ivs = vec![iv("chr1", 0, 100), iv("chr1", 50, 150)];
        assert_eq!(total_covered(&ivs), 150);
    }
}

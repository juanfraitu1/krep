#!/usr/bin/env python3
"""Base-level comparison of two 3-column BEDs (overlap-based P/R/F1)."""
import sys
from collections import defaultdict

def load(path):
    iv = defaultdict(list)
    with open(path) as f:
        for line in f:
            f2 = line.rstrip('\n').split('\t')
            if len(f2) < 3:
                continue
            chrom, b, e = f2[0], int(f2[1]), int(f2[2])
            if e <= b:
                continue
            iv[chrom].append((b, e))
    return iv

def merge(iv):
    merged = {}
    for chrom, lst in iv.items():
        lst.sort()
        out = []
        for b, e in lst:
            if out and b <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((b, e))
        merged[chrom] = out
    return merged

def comp(pred, truth):
    pred_bp = truth_bp = tp = 0
    for chrom in set(pred.keys()) | set(truth.keys()):
        p = pred.get(chrom, [])
        t = truth.get(chrom, [])
        for b, e in p: pred_bp += e - b
        for b, e in t: truth_bp += e - b
        i = j = 0
        while i < len(p) and j < len(t):
            pb, pe = p[i]; tb, te = t[j]
            tp += max(0, min(pe, te) - max(pb, tb))
            if pe < te: i += 1
            else: j += 1
    P = tp/pred_bp if pred_bp else 0
    R = tp/truth_bp if truth_bp else 0
    F = 2*P*R/(P+R) if P+R else 0
    print(f"pred_bp {pred_bp:,} truth_bp {truth_bp:,} tp {tp:,} P {P:.4f} R {R:.4f} F1 {F:.4f}")

if __name__ == '__main__':
    comp(merge(load(sys.argv[1])), merge(load(sys.argv[2])))

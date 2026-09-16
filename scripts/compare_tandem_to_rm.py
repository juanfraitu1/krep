#!/usr/bin/env python3
"""Compare krep --tandem BED to RepeatMasker tandem-related families."""
import sys
from collections import defaultdict

krep_path = sys.argv[1]
rm_path = sys.argv[2]

TANDEM_FAMS = {'Simple_repeat', 'Low_complexity', 'Satellite', 'Satellite/centr', 'Satellite/acro', 'Satellite/subtelo'}

def load_bed(path, fam_filter=None):
    iv = defaultdict(list)
    with open(path) as f:
        for line in f:
            f2 = line.rstrip('\n').split('\t')
            if len(f2) < 3:
                continue
            chrom, b, e = f2[0], int(f2[1]), int(f2[2])
            fam = f2[3] if len(f2) > 3 else 'all'
            if fam_filter and fam not in fam_filter:
                continue
            iv[chrom].append((b, e))
    return iv

def merge_iv(iv):
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

def overlap(pred, truth):
    total_tp = 0
    pred_bp = 0
    truth_bp = 0
    for chrom in set(pred.keys()) | set(truth.keys()):
        p = pred.get(chrom, [])
        t = truth.get(chrom, [])
        for b, e in p:
            pred_bp += e - b
        for b, e in t:
            truth_bp += e - b
        i = j = 0
        while i < len(p) and j < len(t):
            pb, pe = p[i]
            tb, te = t[j]
            tp = max(0, min(pe, te) - max(pb, tb))
            total_tp += tp
            if pe < te:
                i += 1
            else:
                j += 1
    P = total_tp / pred_bp if pred_bp else 0
    R = total_tp / truth_bp if truth_bp else 0
    F = 2 * P * R / (P + R) if P + R else 0
    return P, R, F, pred_bp, truth_bp, total_tp

krep = merge_iv(load_bed(krep_path))
rm_tandem = merge_iv(load_bed(rm_path, TANDEM_FAMS))
P, R, F, pred_bp, truth_bp, tp = overlap(krep, rm_tandem)
print(f"krep tandem vs RM tandem-related families")
print(f"  pred_bp {pred_bp:,}  truth_bp {truth_bp:,}  tp {tp:,}")
print(f"  P {P:.4f}  R {R:.4f}  F1 {F:.4f}")

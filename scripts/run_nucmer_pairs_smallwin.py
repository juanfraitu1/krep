#!/usr/bin/env python3
"""Run nucmer verification on candidate k-mer pairs with small windows.

usage: run_nucmer_pairs_smallwin.py <pairs.tsv> <genome.fa> <sd.bed> <out_prefix> [n_pairs]
"""
import sys, random, subprocess, bisect, re

pairs_path, fa_path, sd_path, out_prefix = sys.argv[1:5]
n_pairs = int(sys.argv[5]) if len(sys.argv) > 5 else 100
window = 2000

sd_iv = []
with open(sd_path) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 3:
            continue
        b, e = int(cols[1]), int(cols[2])
        sd_iv.append((b, e))
sd_iv.sort()
starts = [x[0] for x in sd_iv]
ends = [x[1] for x in sd_iv]

def in_sd(pos):
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return False
    return pos < ends[i]

pos_pairs = []
neg_pairs = []
with open(pairs_path) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        p1, p2 = int(cols[1]), int(cols[2])
        d = abs(p2 - p1)
        if 10000 <= d < 100000:
            a, b = in_sd(p1), in_sd(p2)
            if a and b:
                pos_pairs.append((p1, p2))
            elif not a and not b:
                neg_pairs.append((p1, p2))

random.seed(1)
sample_pos = random.sample(pos_pairs, min(n_pairs, len(pos_pairs)))
sample_neg = random.sample(neg_pairs, min(n_pairs, len(neg_pairs)))
selected = []
for i, (p1, p2) in enumerate(sample_pos):
    selected.append((f"pos_{i}", p1, p2, True))
for i, (p1, p2) in enumerate(sample_neg):
    selected.append((f"neg_{i}", p1, p2, False))

seq = []
name = None
with open(fa_path) as f:
    for line in f:
        if line.startswith('>'):
            name = line[1:].strip().split()[0]
            continue
        seq.append(line.strip().upper())
seq = ''.join(seq)

fa_a = f"{out_prefix}_A.fa"
fa_b = f"{out_prefix}_B.fa"
with open(fa_a, 'w') as out_a, open(fa_b, 'w') as out_b:
    for label, p1, p2, is_pos in selected:
        b = max(0, p1 - window)
        e = min(len(seq), p1 + 50 + window)
        out_a.write(f">{label}\n")
        out_a.write(seq[b:e] + "\n")
        b = max(0, p2 - window)
        e = min(len(seq), p2 + 50 + window)
        out_b.write(f">{label}\n")
        out_b.write(seq[b:e] + "\n")

prefix = f"{out_prefix}_nucmer"
cmd = ["nucmer", "--nosimplify", "-p", prefix, fa_a, fa_b]
print("running", ' '.join(cmd), file=sys.stderr)
subprocess.run(cmd, check=True)

coords_out = f"{prefix}.coords"
cmd = ["show-coords", "-r", "-l", "-T", f"{prefix}.delta"]
print("running", ' '.join(cmd), file=sys.stderr)
with open(coords_out, 'w') as out:
    subprocess.run(cmd, stdout=out, check=True)

verified = []
with open(coords_out) as f:
    for _ in range(4):
        next(f, None)
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 10:
            continue
        tags = cols[9].strip().split()
        if len(tags) != 2:
            continue
        if tags[0] != tags[1]:
            continue
        s1, e1, s2, e2 = int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])
        len1, len2, idy = int(cols[4]), int(cols[5]), float(cols[6])
        if len1 < 200 or idy < 90.0:
            continue
        rec = next((x for x in selected if x[0] == tags[0]), None)
        if rec is None:
            continue
        label, p1, p2, is_pos = rec
        a_start = max(0, p1 - window)
        b_start = max(0, p2 - window)
        chr1_a_start = a_start + min(s1, e1) - 1
        chr1_a_end = a_start + max(s1, e1)
        chr1_b_start = b_start + min(s2, e2) - 1
        chr1_b_end = b_start + max(s2, e2)
        verified.append((label, chr1_a_start, chr1_a_end, chr1_b_start, chr1_b_end, len1, idy, is_pos))

print(f"verified alignments: {len(verified)}", file=sys.stderr)
bed_out = f"{out_prefix}.verified.bed"
with open(bed_out, 'w') as out:
    for v in verified:
        label, a_s, a_e, b_s, b_e, ln, idy, is_pos = v
        out.write(f"NC_060925.1\t{a_s}\t{a_e}\t{label}_A\n")
        out.write(f"NC_060925.1\t{b_s}\t{b_e}\t{label}_B\n")

pos_ver = sum(1 for v in verified if v[7])
neg_ver = sum(1 for v in verified if not v[7])
print(f"positive pairs verified: {pos_ver}/{len(sample_pos)}")
print(f"negative pairs verified: {neg_ver}/{len(sample_neg)}")

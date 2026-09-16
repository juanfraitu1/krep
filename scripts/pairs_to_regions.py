#!/usr/bin/env python3
"""Convert k-mer pair positions to merged BED regions.

usage: pairs_to_regions.py <pairs.tsv> <expand> <min_dist> <max_dist> <out.bed>
"""
import sys

pairs_path, expand, min_d, max_d, out_path = sys.argv[1:6]
expand = int(expand)
min_d = int(min_d)
max_d = int(max_d)

iv = []
with open(pairs_path) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        p1, p2 = int(cols[1]), int(cols[2])
        d = abs(p2 - p1)
        if min_d <= d < max_d:
            for p in (p1, p2):
                b = max(0, p - expand)
                e = p + 64 + expand
                iv.append((b, e))

iv.sort()
merged = []
for b, e in iv:
    if merged and b <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    else:
        merged.append((b, e))

total = 0
with open(out_path, 'w') as out:
    for b, e in merged:
        out.write(f"NC_060925.1\t{b}\t{e}\n")
        total += e - b
print(f"wrote {len(merged)} merged intervals covering {total:,} bp")

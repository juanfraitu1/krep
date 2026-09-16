#!/usr/bin/env python3
import sys, bisect
from collections import defaultdict

pairs_path = sys.argv[1]
sd_path = sys.argv[2]

sd_iv = defaultdict(list)
with open(sd_path) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 3:
            continue
        chrom, b, e = cols[0], int(cols[1]), int(cols[2])
        sd_iv[chrom].append((b, e))

for chrom in sd_iv:
    sd_iv[chrom].sort()
    starts = [x[0] for x in sd_iv[chrom]]
    ends = [x[1] for x in sd_iv[chrom]]
    sd_iv[chrom] = (starts, ends)

def in_sd(chrom, pos):
    if chrom not in sd_iv:
        return False
    starts, ends = sd_iv[chrom]
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return False
    return pos < ends[i]

ranges = [(0, 1000), (1000, 10000), (10000, 100000), (100000, 1000000), (1000000, 10000000), (10000000, 250000000)]
counts = {r: [0, 0, 0, 0] for r in ranges}

with open(pairs_path) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        p1, p2 = int(cols[1]), int(cols[2])
        d = abs(p2 - p1)
        for r in ranges:
            if r[0] <= d < r[1]:
                a = in_sd('NC_060925.1', p1)
                b = in_sd('NC_060925.1', p2)
                counts[r][0] += 1
                if a and b:
                    counts[r][1] += 1
                elif a or b:
                    counts[r][2] += 1
                else:
                    counts[r][3] += 1
                break

print("distance_range\ttotal\tboth_in_SD\tat_least_one\tneither")
for r in ranges:
    total, both, one, none = counts[r]
    if total == 0:
        continue
    print(f"{r[0]}-{r[1]}\t{total}\t{both} ({both/total:.3f})\t{both+one} ({(both+one)/total:.3f})\t{none} ({none/total:.3f})")

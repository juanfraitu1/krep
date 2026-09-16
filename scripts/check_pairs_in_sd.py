#!/usr/bin/env python3
import sys, bisect
from collections import defaultdict

pairs_path = sys.argv[1]
sd_path = sys.argv[2]

# load SD intervals
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

total = both_in = at_least_one = none_in = 0
with open(pairs_path) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        p1, p2 = int(cols[1]), int(cols[2])
        total += 1
        a = in_sd('NC_060925.1', p1)
        b = in_sd('NC_060925.1', p2)
        if a and b:
            both_in += 1
        elif a or b:
            at_least_one += 1
        else:
            none_in += 1

print(f"total pairs: {total:,}")
print(f"both in SD: {both_in:,} ({both_in/total:.4f})")
print(f"at least one in SD: {both_in + at_least_one:,} ({(both_in+at_least_one)/total:.4f})")
print(f"neither in SD: {none_in:,} ({none_in/total:.4f})")

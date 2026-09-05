#!/usr/bin/env python3
"""Label each consensus in a FASTA by the Dfam hit BED produced by krep mask."""
import sys
from collections import defaultdict

fasta = sys.argv[1]
bed = sys.argv[2]

# Parse FASTA headers and lengths
headers = {}
name = None
seq_len = 0
with open(fasta) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                headers[name] = seq_len
            name = line[1:].split()[0]
            seq_len = 0
        else:
            seq_len += len(line.strip())
if name is not None:
    headers[name] = seq_len

# Parse BED: chrom start end label
hits = defaultdict(list)
with open(bed) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 4:
            continue
        chrom, start, end, label = cols[0], int(cols[1]), int(cols[2]), cols[3]
        # label may be name:score; strip score
        family = label.split(":")[0]
        hits[chrom].append((start, end, family))

for name, length in sorted(headers.items(), key=lambda x: int(x[0].split("_")[-1])):
    fam_cov = defaultdict(int)
    for s, e, fam in hits[name]:
        fam_cov[fam] += e - s
    if fam_cov:
        total = sum(fam_cov.values())
        top = sorted(fam_cov.items(), key=lambda x: -x[1])[:3]
        top_str = ", ".join(f"{f}({c/total:.1%})" for f, c in top)
    else:
        top_str = "UNLABELED"
    print(f"{name}\tlen={length}\t{top_str}")

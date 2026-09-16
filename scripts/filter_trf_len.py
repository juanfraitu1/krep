#!/usr/bin/env python3
import sys
in_path, out_path, min_len = sys.argv[1:4]
min_len = int(min_len)
with open(in_path) as f, open(out_path, "w") as out:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 3:
            continue
        chrom, b, e = cols[0], int(cols[1]), int(cols[2])
        if e - b >= min_len:
            out.write(f"NC_060925.1\t{b}\t{e}\n")

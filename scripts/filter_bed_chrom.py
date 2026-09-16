#!/usr/bin/env python3
import sys
chrom = sys.argv[1]
for line in sys.stdin:
    f = line.rstrip('\n').split('\t')
    if f[0] == chrom:
        sys.stdout.write(line)

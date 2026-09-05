#!/usr/bin/env python3
import sys
from collections import defaultdict
import bisect

path = sys.argv[1]
with open(path) as f:
    hdr = None
    for i, line in enumerate(f, 1):
        fcols = line.rstrip('\n').split('\t')
        if hdr is None or fcols[0] == 'chrom':
            hdr = fcols
            continue
        if len(fcols) != len(hdr):
            print(f"line {i}: field count mismatch {len(fcols)} vs {len(hdr)}")
            continue
        r = dict(zip(hdr, fcols))
        if 'consensus' not in r:
            print(f"line {i}: missing consensus: {line[:200]!r}")
            print("  keys:", list(r.keys()))
            break

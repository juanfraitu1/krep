#!/usr/bin/env python3
import sys
path = sys.argv[1]
with open(path) as f:
    hdr = f.readline().rstrip('\n').split('\t')
    print("header fields:", len(hdr), hdr)
    bad = 0
    for i, line in enumerate(f, 2):
        cols = line.rstrip('\n').split('\t')
        if len(cols) != len(hdr):
            print(f"line {i}: {len(cols)} fields: {line[:200]!r}")
            bad += 1
            if bad >= 10:
                break
    print("bad lines:", bad)

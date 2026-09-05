#!/usr/bin/env python3
from collections import Counter
import sys
lengths = []
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith(">"):
            for part in line.split():
                if part.startswith("len="):
                    lengths.append(int(part.split("=")[1]))
c = Counter(lengths)
for l in sorted(c):
    print(l, c[l])

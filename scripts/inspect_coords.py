#!/usr/bin/env python3
import sys
path = sys.argv[1]
same = 0
all_lines = 0
with open(path) as f:
    for _ in range(4):
        next(f, None)
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 10:
            continue
        all_lines += 1
        tags = cols[9].strip().split()
        if len(tags) != 2:
            continue
        if tags[0] == tags[1]:
            same += 1
            print(line.rstrip('\n'))
print(f"same-tag alignments: {same}/{all_lines}", file=sys.stderr)

#!/usr/bin/env python3
import sys
path, tau = sys.argv[1], sys.argv[2]
with open(path) as f:
    lines = f.readlines()
parts = lines[0].rstrip('\n').split('\t')
parts[-1] = tau
lines[0] = '\t'.join(parts) + '\n'
with open(path, 'w') as f:
    f.writelines(lines)

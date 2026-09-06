#!/usr/bin/env python3
import re, sys
pat = re.compile(sys.argv[2], re.I) if len(sys.argv) > 2 else re.compile(r'CR1|Helitron|MIR|L2|Tip100', re.I)
in_rec = False
for l in open(sys.argv[1]):
    if l.startswith('>'):
        in_rec = pat.search(l) is not None
    if in_rec:
        sys.stdout.write(l)

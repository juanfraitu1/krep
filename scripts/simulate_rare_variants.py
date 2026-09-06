#!/usr/bin/env python3
"""Generate mutated variants of a consensus FASTA for augmentation.

Each input sequence is copied N times, with each base mutated at rate `rate`
to a random different base. Output header keeps the original name plus the
variant index.

usage: simulate_rare_variants.py <input.fa> <variants_per_seq> <rate> <output.fa>
"""
import sys, random

in_p, n_var, rate, out_p = sys.argv[1:5]
n_var = int(n_var)
rate = float(rate)
random.seed(7)

ALTS = {b'A': b'CGT', b'C': b'AGT', b'G': b'ACT', b'T': b'ACG'}

def mutate(seq, rate):
    out = bytearray()
    for b in seq:
        if random.random() < rate and b in ALTS:
            out.append(random.choice(ALTS[b]))
        else:
            out.append(b)
    return bytes(out)

records = []
name = None; seq = []
for l in open(in_p, 'rb'):
    if l.startswith(b'>'):
        if name:
            records.append((name, b''.join(seq)))
        name = l[1:].split()[0].decode()
        seq = []
    else:
        seq.append(l.strip().upper())
if name:
    records.append((name, b''.join(seq)))

with open(out_p, 'wb') as f:
    for name, s in records:
        f.write(f">{name}\n".encode())
        for i in range(0, len(s), 80):
            f.write(s[i:i+80] + b'\n')
        for i in range(n_var):
            v = mutate(s, rate)
            f.write(f">{name}_var{i+1} rate={rate}\n".encode())
            for j in range(0, len(v), 80):
                f.write(v[j:j+80] + b'\n')
print(f"wrote {out_p}: {len(records)} originals + {len(records)*n_var} variants")

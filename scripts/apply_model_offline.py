#!/usr/bin/env python3
"""Apply a learned logistic model to a candidate dump offline and write a BED.

This is much faster than re-running `krep mask` because the candidates are
already computed. It is also useful for ensembling multiple models: compute
probabilities for each model, accept if any model accepts, then merge BEDs.

usage: apply_model_offline.py <model.tsv> <cand.tsv> <out.bed> [tau]
"""
import sys, math
from collections import defaultdict

model_p, cand_p, out_p = sys.argv[1:4]
tau_override = float(sys.argv[4]) if len(sys.argv) > 4 else None

def load_model(path):
    lines = open(path).readlines()
    nums = lambda l: [float(x) for x in l.split('\t')[1:]]
    head = nums(lines[0])
    if len(head) != 10:
        raise ValueError("#logistic line needs 8 weights + b0 + tau")
    mu = nums(lines[1]); sd = nums(lines[2])
    if len(mu) != 8 or len(sd) != 8:
        raise ValueError("#mu/#sd need 8 values")
    offsets = {}
    for l in lines[3:]:
        f = l.rstrip('\n').split('\t')
        if len(f) >= 2:
            offsets[f[0]] = float(f[1])
    tau = head[9]
    return head[:8], head[8], mu, sd, offsets, tau

w, b0, mu, sd, offsets, tau = load_model(model_p)
tau = tau_override if tau_override is not None else tau
isd = [1.0 / s if s else 1.0 for s in sd]

def prob(row):
    cons = row['consensus']
    L = int(row['end']) - int(row['start'])
    s = int(row['score']); gl = int(row['gate_left']); gr = int(row['gate_right'])
    cov = (int(row['cons_fwd']) + int(row['cons_bwd'])) / max(1, int(row['cons_len']))
    x = [
        s, math.log(L), s / L, cov, float(row['gc']),
        gl + gr, min(gl, gr), math.log(int(row['cons_len']))
    ]
    zz = b0 + offsets.get(cons, 0.0)
    for j in range(8):
        zz += w[j] * (x[j] - mu[j]) * isd[j]
    zz = max(-30, min(30, zz))
    return 1 / (1 + math.exp(-zz))

# Read candidates, filter by probability, merge intervals per chrom
by_chrom = defaultdict(list)
hdr = None
for l in open(cand_p):
    f = l.rstrip('\n').split('\t')
    if hdr is None or f[0] == 'chrom':
        hdr = f; continue
    r = dict(zip(hdr, f))
    if prob(r) >= tau:
        chrom = r['chrom']; b = int(r['start']); e = int(r['end'])
        by_chrom[chrom].append((b, e, r['consensus']))

with open(out_p, 'w') as out:
    for chrom in sorted(by_chrom.keys()):
        ivs = sorted(by_chrom[chrom])
        cs, ce, cfam = ivs[0]
        for b, e, fam in ivs[1:]:
            if b <= ce:
                if e > ce: ce = e
            else:
                out.write(f"{chrom}\t{cs}\t{ce}\t{cfam}\n")
                cs, ce, cfam = b, e, fam
        out.write(f"{chrom}\t{cs}\t{ce}\t{cfam}\n")
print(f"wrote {out_p} with tau={tau}")

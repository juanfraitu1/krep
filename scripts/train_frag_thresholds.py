#!/usr/bin/env python3
"""Train per-consensus score thresholds for the fragment library.

usage: train_frag_thresholds.py <train.tsv[,...]> <train_rm.bed[,...]> \
                              <labels.tsv> <out_thresholds.tsv> [min_score]
"""
import sys, bisect
from collections import defaultdict

train_p, train_rm, labels_p, out_p = sys.argv[1:5]
min_score = int(sys.argv[5]) if len(sys.argv) > 5 else 8

GROUPS = [
    ('SINE/Alu', {'SINE/Alu'}),
    ('LINE/L1', {'LINE/L1'}),
    ('LTR/ERVK', {'LTR/ERVK'}),
    ('LTR/ERV1', {'LTR/ERV1'}),
    ('LTR/ERVL-MaLR', {'LTR/ERVL-MaLR'}),
    ('LTR/ERVL', {'LTR/ERVL'}),
    ('DNA/TcMar', {'DNA/TcMar-Tigger', 'DNA/TcMar-Tc1', 'DNA/TcMar-Tc2',
                   'DNA/TcMar-Mariner', 'DNA/TcMar-Pogo'}),
    ('DNA/hAT', {'DNA/hAT-Charlie', 'DNA/hAT-Tip100', 'DNA/hAT-Blackjack',
                 'DNA/hAT-Ac', 'DNA/hAT-Tag1'}),
    ('Satellite', {'Satellite', 'Satellite/centr', 'Satellite/acro',
                   'Satellite/subtelo'}),
    ('Simple/Low', {'Simple_repeat', 'Low_complexity'}),
    ('Ancient', {'SINE/MIR', 'LINE/L2', 'LINE/CR1', 'RC/Helitron',
                 'SINE/5S-Deu-L2', 'SINE/tRNA', 'SINE/tRNA-Deu', 'SINE/tRNA-RTE'}),
    ('Other', set()),
]

name_to_fam = {}
for l in open(labels_p):
    if l.startswith('#'):
        continue
    f = l.rstrip('\n').split('\t')
    if len(f) >= 2:
        name_to_fam[f[0]] = f[1]
print(f"loaded {len(name_to_fam)} consensus labels")

def load_rm(paths):
    iv = defaultdict(list)
    for part in paths.split(','):
        for l in open(part):
            f = l.split()
            if len(f) >= 3:
                iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in iv.items():
        v.sort(); st, en = [], []
        for b, e in v:
            if len(en) and b <= en[-1]: en[-1] = max(en[-1], e)
            else: st.append(b); en.append(e)
        merged[c] = (st, en)
    return merged

def overlap(merged, c, b, e):
    if c not in merged:
        return 0
    st, en = merged[c]
    i = max(0, bisect.bisect_left(st, b) - 1); t = 0
    while i < len(st) and st[i] < e:
        t += max(0, min(e, en[i]) - max(b, st[i])); i += 1
    return t

tr_m = load_rm(train_rm)

by_cons = defaultdict(list)  # (score, y, L)
hdr = None
for part in train_p.split(','):
    for l in open(part):
        f = l.rstrip('\n').split('\t')
        if hdr is None or f[0] == 'chrom':
            hdr = f; continue
        r = dict(zip(hdr, f))
        c = r['consensus']
        if c not in name_to_fam:
            continue
        b = int(r['start']); e = int(r['end']); L = e - b
        score = int(r['score'])
        y = 1 if overlap(tr_m, r['chrom'], b, e) * 2 >= L else 0
        by_cons[c].append((score, y, L))

thr = {}
for c, vals in by_cons.items():
    vals.sort(key=lambda x: -x[0])
    best_t = min_score
    best_net = -1e18
    net = 0
    i = 0
    while i < len(vals):
        t = vals[i][0]
        j = i
        while j < len(vals) and vals[j][0] == t:
            # weight by length: TP adds L, FP subtracts L
            net += vals[j][2] if vals[j][1] else -vals[j][2]
            j += 1
        if net > best_net:
            best_net = net; best_t = t
        i = j
    thr[c] = max(min_score, best_t)

for c in name_to_fam:
    thr.setdefault(c, min_score)

with open(out_p, 'w') as f:
    for c in sorted(thr.keys()):
        f.write(f"{c}\t{thr[c]}\n")
print("thresholds written to", out_p)

# quick stats
tvals = list(thr.values())
print(f"thresholds: min={min(tvals)}, max={max(tvals)}, median={sorted(tvals)[len(tvals)//2]}")

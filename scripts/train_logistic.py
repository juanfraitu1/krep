"""Logistic hit filter with per-consensus offsets, trained on one chromosome's
candidate dump and evaluated offline on another's.

usage: train_logistic.py <train.tsv> <train_rm.bed> <test.tsv> <test_rm.bed> <thresholds.tsv> <out_model.tsv>
Offline evaluation scores kept candidates base-by-base against RepeatMasker:
TP = kept bases inside annotation, FP = kept bases outside, FN = annotated bases not kept.
"""
import sys, bisect, math, random
from collections import defaultdict
train_p, train_rm, test_p, test_rm, thr_p, out_p = sys.argv[1:7]

def load_rm(path):
    rm = defaultdict(list)
    for l in open(path):
        f = l.split()
        if len(f) >= 3: rm[f[0]].append((int(f[1]), int(f[2])))
    merged = {}; total = 0
    for c, iv in rm.items():
        iv.sort(); out = []
        for b, e in iv:
            if out and b <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], e))
            else: out.append((b, e))
        merged[c] = (out, [b for b, _ in out]); total += sum(e - b for b, e in out)
    return merged, total
def overlap(merged, c, b, e):
    if c not in merged: return 0
    iv, starts = merged[c]; i = max(0, bisect.bisect_left(starts, b) - 1); t = 0
    while i < len(iv) and iv[i][0] < e:
        t += max(0, min(e, iv[i][1]) - max(b, iv[i][0])); i += 1
    return t

NF = 8
def load(path, merged):
    rows = []; hdr = None
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if hdr is None: hdr = f; continue
        r = dict(zip(hdr, f)); b, e = int(r['start']), int(r['end']); L = e - b
        s = int(r['score']); gl, gr = int(r['gate_left']), int(r['gate_right'])
        cov = (int(r['cons_fwd']) + int(r['cons_bwd'])) / max(1, int(r['cons_len']))
        x = [s, math.log(L), s / L, cov, float(r['gc']), gl + gr, min(gl, gr), math.log(int(r['cons_len']))]
        tp = overlap(merged, r['chrom'], b, e)
        rows.append((r['consensus'], x, tp, L - tp, L, s, r['chrom'], b, e))
    return rows

tr_m, tr_tot = load_rm(train_rm); te_m, te_tot = load_rm(test_rm)
train = load(train_p, tr_m); test = load(test_p, te_m)
print(f"train {len(train):,} rows, test {len(test):,} rows")

# standardize features on training data
mu = [sum(r[1][j] for r in train) / len(train) for j in range(NF)]
sd = [math.sqrt(sum((r[1][j] - mu[j]) ** 2 for r in train) / len(train)) or 1.0 for j in range(NF)]
def z(x): return [(x[j] - mu[j]) / sd[j] for j in range(NF)]
cons = sorted({r[0] for r in train} | {r[0] for r in test}); cidx = {c: i for i, c in enumerate(cons)}

# SGD logistic regression, base-weighted loss, per-consensus offset with L2 shrinkage
w = [0.0] * NF; b0 = 0.0; off = [0.0] * len(cons)
lr, lam, epochs = 0.02, 1e-3, 6
meanL = sum(r[4] for r in train) / len(train)
data = [(cidx[c], z(x), 1.0 if tp * 2 >= L else 0.0, L / meanL) for c, x, tp, fp, L, s, *_ in train]
random.seed(1)
for ep in range(epochs):
    random.shuffle(data); loss = 0.0
    for ci, x, y, wt in data:
        zz = b0 + off[ci] + sum(w[j] * x[j] for j in range(NF)); zz = max(-30, min(30, zz))
        p = 1 / (1 + math.exp(-zz)); g = (p - y) * wt
        loss += -wt * (y * math.log(p + 1e-12) + (1 - y) * math.log(1 - p + 1e-12))
        for j in range(NF): w[j] -= lr * (g * x[j] + lam * w[j])
        b0 -= lr * g; off[ci] -= lr * (g + 10 * lam * off[ci])
    print(f"  epoch {ep+1}: weighted loss {loss/len(data):.4f}")
print("weights:", {n: round(v, 3) for n, v in zip(['score','logL','score/L','cov','gc','gate_sum','gate_min','logclen'], w)}, "b0", round(b0, 3))

def prob(ci, x):
    zz = b0 + off[ci] + sum(w[j] * x[j] for j in range(NF)); return 1 / (1 + math.exp(-max(-30, min(30, zz))))
def evaluate(rows, keep, total, merged):
    # Kept candidates overlap (block boundaries, nested hits): merge first so
    # every base is counted once, exactly as the masker's output is scored.
    by_c = defaultdict(list)
    for r in rows:
        if keep(r): by_c[r[6]].append((r[7], r[8]))
    tp = kept = 0
    for c, iv in by_c.items():
        iv.sort(); cs, ce = iv[0]
        for b, e in iv[1:]:
            if b > ce:
                kept += ce - cs; tp += overlap(merged, c, cs, ce); cs, ce = b, e
            else: ce = max(ce, e)
        kept += ce - cs; tp += overlap(merged, c, cs, ce)
    fp = kept - tp
    P = tp / max(1, tp + fp); R = tp / total; F = 2 * P * R / max(1e-9, P + R); return P, R, F
# pick the probability cutoff on the training set
best = (0, 0.5)
for tau in [i / 20 for i in range(2, 19)]:
    F = evaluate(train, lambda r: prob(cidx[r[0]], z(r[1])) >= tau, tr_tot, tr_m)[2]
    if F > best[0]: best = (F, tau)
tau = best[1]; print(f"cutoff tau={tau} (train F1 {best[0]:.4f})")
th = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(thr_p)}
print("offline on test chromosome (base-level vs RepeatMasker):")
for name, keep in [("floor 30, no model", lambda r: r[5] >= 30),
                   ("per-consensus thresholds", lambda r: r[5] >= th.get(r[0], 30)),
                   ("logistic + per-consensus offset", lambda r: prob(cidx[r[0]], z(r[1])) >= tau)]:
    P, R, F = evaluate(test, keep, te_tot, te_m); print(f"  {name:<34} P {P:.4f}  R {R:.4f}  F1 {F:.4f}")
with open(out_p, 'w') as f:
    f.write("#logistic\t" + "\t".join(f"{v:.6f}" for v in w) + f"\t{b0:.6f}\t{tau}\n#mu\t" + "\t".join(f"{v:.6f}" for v in mu) + "\n#sd\t" + "\t".join(f"{v:.6f}" for v in sd) + "\n")
    for c in cons: f.write(f"{c}\t{off[cidx[c]]:.6f}\n")
print("model written to", out_p)

"""Logistic hit filter with per-consensus offsets, trained on one or more
chromosomes' candidate dumps and evaluated offline on a held-out chromosome.

usage: train_logistic.py <train.tsv[,train2.tsv,...]> <train_rm.bed[,train_rm2.bed,...]> <test.tsv> <test_rm.bed> <thresholds.tsv> <out_model.tsv>
Comma-separated lists are concatenated (train side and/or test side), so the
model can be fit on several chromosomes at once.

Rows are held in flat typed arrays (array module), not object tuples: several
chromosomes of candidate dumps (5M+ hits, ~1 GB of TSV) fit in well under a
gigabyte of RAM, which a tuple-per-row layout would not.

Offline evaluation scores kept candidates base-by-base against RepeatMasker:
TP = kept bases inside annotation, FP = kept bases outside, FN = annotated
bases not kept.
"""
import sys, bisect, math, random
from array import array
from collections import defaultdict
train_p, train_rm, test_p, test_rm, thr_p, out_p = sys.argv[1:7]

def load_rm(path):
    iv = defaultdict(list)
    for part in path.split(','):
        for l in open(part):
            f = l.split()
            if len(f) >= 3: iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}; total = 0
    for c, v in iv.items():
        v.sort(); st, en = array('i'), array('i')
        for b, e in v:
            if len(en) and b <= en[-1]: en[-1] = max(en[-1], e)
            else: st.append(b); en.append(e)
        merged[c] = (st, en)
        total += sum(e - b for b, e in zip(st, en))
    return merged, total

def overlap(merged, c, b, e):
    if c not in merged: return 0
    st, en = merged[c]; i = max(0, bisect.bisect_left(st, b) - 1); t = 0
    while i < len(st) and st[i] < e:
        t += max(0, min(e, en[i]) - max(b, st[i])); i += 1
    return t

NF = 8
def load(path, merged):
    """Candidate hits as parallel typed arrays; one row per accepted hit."""
    D = {'feat': array('d'), 'ci': array('i'), 'chi': array('i'),
         'tp': array('i'), 'L': array('i'), 'sc': array('i'),
         'b': array('i'), 'e': array('i'), 'cons': [], 'chroms': []}
    cmap = {}; chromap = {}; hdr = None
    for part in path.split(','):
        for l in open(part):
            f = l.rstrip('\n').split('\t')
            if hdr is None or f[0] == 'chrom': hdr = f; continue
            r = dict(zip(hdr, f))
            b = int(r['start']); e = int(r['end']); L = e - b
            s = int(r['score']); gl = int(r['gate_left']); gr = int(r['gate_right'])
            cov = (int(r['cons_fwd']) + int(r['cons_bwd'])) / max(1, int(r['cons_len']))
            for v in (s, math.log(L), s / L, cov, float(r['gc']),
                      gl + gr, min(gl, gr), math.log(int(r['cons_len']))):
                D['feat'].append(v)
            c = r['consensus']
            if c not in cmap: cmap[c] = len(D['cons']); D['cons'].append(c)
            ch = r['chrom']
            if ch not in chromap: chromap[ch] = len(D['chroms']); D['chroms'].append(ch)
            D['ci'].append(cmap[c]); D['chi'].append(chromap[ch])
            D['tp'].append(overlap(merged, ch, b, e))
            D['L'].append(L); D['sc'].append(s); D['b'].append(b); D['e'].append(e)
    return D

tr_m, tr_tot = load_rm(train_rm); te_m, te_tot = load_rm(test_rm)
train = load(train_p, tr_m); test = load(test_p, te_m)
ntr, nte = len(train['L']), len(test['L'])
print(f"train {ntr:,} rows, test {nte:,} rows")

# standardize features on training data
mu = [sum(train['feat'][j::NF]) / ntr for j in range(NF)]
sd = [math.sqrt(sum(x * x for x in train['feat'][j::NF]) / ntr - mu[j] ** 2) or 1.0
      for j in range(NF)]
isd = [1.0 / s for s in sd]
cons = sorted(set(train['cons']) | set(test['cons'])); cidx = {c: i for i, c in enumerate(cons)}
tmap_tr = [cidx[c] for c in train['cons']]   # local cons idx -> global idx
tmap_te = [cidx[c] for c in test['cons']]

# SGD logistic regression, base-weighted loss, per-consensus offset with L2 shrinkage
w = [0.0] * NF; b0 = 0.0; off = [0.0] * len(cons)
lr, lam, epochs = 0.02, 1e-3, 6
meanL = sum(train['L']) / ntr
F, TP, Lr, CI = train['feat'], train['tp'], train['L'], train['ci']
order = array('i', range(ntr))
random.seed(1)
for ep in range(epochs):
    random.shuffle(order); loss = 0.0
    for i in order:
        base = i * NF
        x = [(F[base + j] - mu[j]) * isd[j] for j in range(NF)]
        t = tmap_tr[CI[i]]
        zz = b0 + off[t] + sum(w[j] * x[j] for j in range(NF)); zz = max(-30, min(30, zz))
        p = 1 / (1 + math.exp(-zz))
        Li = Lr[i]; y = 1.0 if TP[i] * 2 >= Li else 0.0
        g = (p - y) * (Li / meanL)
        loss += -(Li / meanL) * (y * math.log(p + 1e-12) + (1 - y) * math.log(1 - p + 1e-12))
        for j in range(NF): w[j] -= lr * (g * x[j] + lam * w[j])
        b0 -= lr * g; off[t] -= lr * (g + 10 * lam * off[t])
    print(f"  epoch {ep+1}: weighted loss {loss/ntr:.4f}")
print("weights:", {n: round(v, 3) for n, v in zip(['score','logL','score/L','cov','gc','gate_sum','gate_min','logclen'], w)}, "b0", round(b0, 3))

def prob(i, D, tm):
    base = i * NF
    zz = b0 + off[tm[D['ci'][i]]] + sum(w[j] * (D['feat'][base + j] - mu[j]) * isd[j] for j in range(NF))
    return 1 / (1 + math.exp(-max(-30, min(30, zz))))
def evaluate(D, keep, total, merged, tm):
    # Kept candidates overlap (block boundaries, nested hits): merge first so
    # every base is counted once, exactly as the masker's output is scored.
    by_c = defaultdict(list)
    for i in range(len(D['L'])):
        if keep(i): by_c[D['chroms'][D['chi'][i]]].append((D['b'][i], D['e'][i]))
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
    F = evaluate(train, lambda i, t=tau: prob(i, train, tmap_tr) >= t, tr_tot, tr_m, tmap_tr)[2]
    if F > best[0]: best = (F, tau)
tau = best[1]; print(f"cutoff tau={tau} (train F1 {best[0]:.4f})")
th = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(thr_p)}
print("offline on test chromosome (base-level vs RepeatMasker):")
for name, keep in [("floor 30, no model", lambda i: test['sc'][i] >= 30),
                   ("per-consensus thresholds", lambda i: test['sc'][i] >= th.get(test['cons'][test['ci'][i]], 30)),
                   ("logistic + per-consensus offset", lambda i: prob(i, test, tmap_te) >= tau)]:
    P, R, F = evaluate(test, keep, te_tot, te_m, tmap_te); print(f"  {name:<34} P {P:.4f}  R {R:.4f}  F1 {F:.4f}")
with open(out_p, 'w') as f:
    f.write("#logistic\t" + "\t".join(f"{v:.6f}" for v in w) + f"\t{b0:.6f}\t{tau}\n#mu\t" + "\t".join(f"{v:.6f}" for v in mu) + "\n#sd\t" + "\t".join(f"{v:.6f}" for v in sd) + "\n")
    for c in cons: f.write(f"{c}\t{off[cidx[c]]:.6f}\n")
print("model written to", out_p)
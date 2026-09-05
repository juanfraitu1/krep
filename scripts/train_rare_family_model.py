#!/usr/bin/env python3
"""Train a specialist acceptor model for underrepresented repeat families.

The global logistic model is dominated by Alu/L1/satellite examples and
struggles to learn the very different feature distribution of short, sparse,
low-score ancient hits (CR1, Helitron, Tip100, L2, MIR). This script trains a
separate model whose positive label is "candidate overlaps one of the target
families" and whose negatives include everything else (both non-repeat FP
regions and other repeat families). At test time, a hit is accepted if either
the main model or the specialist model accepts it.

usage: train_rare_family_model.py <train.tsv[,...]> <train_rm.bed[,...]> \
                                  <test.tsv> <test_rm.bed> <out_model.tsv>

Target families are hardcoded below; adjust TARGETS as needed.
"""
import sys, bisect, math, random
from array import array
from collections import defaultdict

train_p, train_rm, test_p, test_rm, out_p = sys.argv[1:6]

# Focus on the biologically important ancient/rare families that the main
# model under-recovers.  Lowering the global cutoff helps these, but at a
# large precision cost.  A specialist model should learn their structure
# without touching Alu/L1 precision.
TARGETS = {'LINE/CR1', 'RC/Helitron', 'DNA/hAT-Tip100', 'LINE/L2', 'SINE/MIR'}

NF = 8
random.seed(1)

def load_rm(paths):
    iv = defaultdict(lambda: defaultdict(list))
    for part in paths.split(','):
        for l in open(part):
            f = l.rstrip('\n').split('\t')
            if len(f) < 4: continue
            chrom, s, e, family = f[0], int(f[1]), int(f[2]), f[3]
            iv[chrom][family].append((s, e))
    merged = {}
    for chrom, fams in iv.items():
        merged[chrom] = {}
        for fam, v in fams.items():
            v.sort(); st, en = array('i'), array('i')
            for b, e in v:
                if len(en) and b <= en[-1]: en[-1] = max(en[-1], e)
                else: st.append(b); en.append(e)
            merged[chrom][fam] = (st, en)
    return merged

def overlap(merged, c, b, e):
    if c not in merged: return 0, set()
    total = 0; fams = set()
    for fam, (st, en) in merged[c].items():
        i = max(0, bisect.bisect_left(st, b) - 1); t = 0
        while i < len(st) and st[i] < e:
            t += max(0, min(e, en[i]) - max(b, st[i])); i += 1
        if t:
            total += t
            fams.add(fam)
    return total, fams

def load(path, merged):
    D = {'feat': array('d'), 'ci': array('i'), 'chi': array('i'),
         'L': array('i'), 'sc': array('i'), 'b': array('i'), 'e': array('i'),
         'cons': [], 'chroms': [], 'target': array('b'), 'any': array('b')}
    cmap = {}; chromap = {}; hdr = None
    for part in path.split(','):
        for l in open(part):
            f = l.rstrip('\n').split('\t')
            if hdr is None or f[0] == 'chrom':
                hdr = f; continue
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
            cov, fams = overlap(merged, ch, b, e)
            D['ci'].append(cmap[c]); D['chi'].append(chromap[ch])
            D['L'].append(L); D['sc'].append(s); D['b'].append(b); D['e'].append(e)
            D['target'].append(1 if fams & TARGETS else 0)
            D['any'].append(1 if fams else 0)
    return D

tr_m = load_rm(train_rm); te_m = load_rm(test_rm)
train = load(train_p, tr_m); test = load(test_p, te_m)
ntr, nte = len(train['L']), len(test['L'])
nt_pos = sum(1 for i in range(ntr) if train['target'][i])
print(f"train {ntr:,} rows, {nt_pos:,} target positives ({100*nt_pos/ntr:.2f}%)")
print(f"test  {nte:,} rows")

# Standardize features on training data
mu = [sum(train['feat'][j::NF]) / ntr for j in range(NF)]
sd = [math.sqrt(sum(x * x for x in train['feat'][j::NF]) / ntr - mu[j] ** 2) or 1.0
      for j in range(NF)]
isd = [1.0 / s for s in sd]

cons = sorted(set(train['cons']) | set(test['cons'])); cidx = {c: i for i, c in enumerate(cons)}
tmap_tr = [cidx[c] for c in train['cons']]
tmap_te = [cidx[c] for c in test['cons']]

# SGD logistic regression, base-weighted loss, per-consensus offset
w = [0.0] * NF; b0 = 0.0; off = [0.0] * len(cons)
lr, lam, epochs = 0.02, 1e-3, 8
meanL = sum(train['L']) / ntr
F, TP, Lr, CI, Y = train['feat'], train['target'], train['L'], train['ci'], train['target']
order = array('i', range(ntr))
for ep in range(epochs):
    random.shuffle(order); loss = 0.0
    for i in order:
        base = i * NF
        x = [(F[base + j] - mu[j]) * isd[j] for j in range(NF)]
        t = tmap_tr[CI[i]]
        zz = b0 + off[t] + sum(w[j] * x[j] for j in range(NF)); zz = max(-30, min(30, zz))
        p = 1 / (1 + math.exp(-zz))
        Li = Lr[i]; y = float(Y[i])
        g = (p - y) * (Li / meanL)
        loss += -(Li / meanL) * (y * math.log(p + 1e-12) + (1 - y) * math.log(1 - p + 1e-12))
        for j in range(NF): w[j] -= lr * (g * x[j] + lam * w[j])
        b0 -= lr * g; off[t] -= lr * (g + 10 * lam * off[t])
    print(f"  epoch {ep+1}: weighted loss {loss/ntr:.4f}")
print("weights:", {n: round(v, 3) for n, v in zip(
    ['score','logL','score/L','cov','gc','gate_sum','gate_min','logclen'], w)}, "b0", round(b0, 3))

def prob(i, D, tm):
    base = i * NF
    zz = b0 + off[tm[D['ci'][i]]] + sum(w[j] * (D['feat'][base + j] - mu[j]) * isd[j] for j in range(NF))
    return 1 / (1 + math.exp(-max(-30, min(30, zz))))

def evaluate(D, keep):
    by_c = defaultdict(list)
    for i in range(len(D['L'])):
        if keep(i): by_c[D['chroms'][D['chi'][i]]].append((D['b'][i], D['e'][i]))
    target_bases = kept = 0
    for c, iv in by_c.items():
        if c not in te_m and c not in tr_m: continue
        m = te_m if c in te_m else tr_m
        iv.sort(); cs, ce = iv[0]
        for b, e in iv[1:]:
            if b > ce:
                seg = ce - cs; kept += seg
                _, fams = overlap(m, c, cs, ce)
                if fams & TARGETS: target_bases += seg
                cs, ce = b, e
            else: ce = max(ce, e)
        seg = ce - cs; kept += seg
        _, fams = overlap(m, c, cs, ce)
        if fams & TARGETS: target_bases += seg
    # Compare to total target truth bases
    truth_target = 0
    for c, fams in (te_m if te_m else tr_m).items():
        for fam, (st, en) in fams.items():
            if fam in TARGETS:
                truth_target += sum(en[i] - st[i] for i in range(len(st)))
    recall = target_bases / truth_target if truth_target else 0
    precision = target_bases / kept if kept else 0
    return precision, recall, kept, target_bases, truth_target

# Pick a cutoff that maximizes target recall while keeping specialist precision high
best = (0, 0.5, 0.0, 0.0)
print("\ntau sweep on training set:")
for tau in [i / 20 for i in range(2, 19)]:
    P, R, kept, tb, tt = evaluate(train, lambda i, t=tau: prob(i, train, tmap_tr) >= t)
    # Objective: recall, but with a precision floor to avoid flooding
    obj = R if P >= 0.70 else R * (P / 0.70)
    print(f"  tau={tau:.2f} P={P:.4f} R={R:.4f} obj={obj:.4f}")
    if obj > best[0]:
        best = (obj, tau, P, R)
tau = best[1]
print(f"chosen tau={tau} (train P={best[2]:.4f} R={best[3]:.4f})")

print("\noffline test chromosome (target-family bases only):")
P, R, kept, tb, tt = evaluate(test, lambda i: prob(i, test, tmap_te) >= tau)
print(f"  specialist only  P {P:.4f}  R {R:.4f}  kept {kept:,} target bases {tb:,}/{tt:,}")

with open(out_p, 'w') as f:
    f.write("#logistic\t" + "\t".join(f"{v:.6f}" for v in w) + f"\t{b0:.6f}\t{tau}\n#mu\t" + "\t".join(f"{v:.6f}" for v in mu) + "\n#sd\t" + "\t".join(f"{v:.6f}" for v in sd) + "\n")
    for c in cons: f.write(f"{c}\t{off[cidx[c]]:.6f}\n")
print("model written to", out_p)

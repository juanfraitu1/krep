"""Gradient-boosted (LightGBM) hit filter over the same features as the
logistic filter, plus consensus identity as a native categorical feature —
the thing a linear model with per-consensus offsets cannot express: family x
divergence interactions.

usage: train_gbm.py <train.tsv[,train2.tsv,...]> <train_rm.bed[,rm2.bed,...]> <test.tsv> <test_rm.bed> <out_report.txt>
Offline evaluation mirrors train_logistic.py exactly: kept candidate intervals
are merged per chromosome, then scored base-by-base against RepeatMasker.
"""
import sys, math
from collections import defaultdict
import numpy as np
import lightgbm as lgb

train_p, train_rm, test_p, test_rm, out_rep = sys.argv[1:6]

def load_rm(path):
    iv = defaultdict(list)
    for part in path.split(','):
        for l in open(part):
            f = l.split()
            if len(f) >= 3: iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}; total = 0
    for c, v in iv.items():
        v.sort(); st, en = [], []
        for b, e in v:
            if en and b <= en[-1]: en[-1] = max(en[-1], e)
            else: st.append(b); en.append(e)
        rs = np.array(st, dtype=np.int64); re = np.array(en, dtype=np.int64)
        lens = re - rs
        cpre = np.concatenate(([0], np.cumsum(lens)[:-1]))
        merged[c] = (rs, re, cpre)
        total += int(lens.sum())
    return merged, total

def coverage(merged, c, x):
    """Bases of RepeatMasker annotation in [0, x) on chromosome c, vectorized."""
    x = np.asarray(x, dtype=np.int64)
    if c not in merged or len(x) == 0: return np.zeros(len(x))
    rs, re, cpre = merged[c]
    k = np.searchsorted(rs, x, side='right') - 1
    kc = np.clip(k, 0, len(rs) - 1)
    f = cpre[kc] + np.maximum(0, np.minimum(re[kc], x) - rs[kc])
    return np.where(k >= 0, f, 0)

# candidate dump columns (fixed order, verified against the header):
# 0 chrom 1 start 2 end 3 consensus 4 strand 5 score 6 score_fwd 7 score_bwd
# 8 cons_fwd 9 cons_bwd 10 cons_len 11 gate_left 12 gate_right 13 gc
COLS = ['chrom', 'start', 'end', 'consensus', 'strand', 'score', 'score_fwd', 'score_bwd',
        'cons_fwd', 'cons_bwd', 'cons_len', 'gate_left', 'gate_right', 'gc']

def load(path):
    D = {k: [] for k in COLS}
    cons_map = {}; chrom_map = {}
    for part in path.split(','):
        hdr = None
        for l in open(part):
            f = l.rstrip('\n').split('\t')
            if hdr is None or f[0] == 'chrom': hdr = f; continue
            if hdr[1] != 'start' or hdr[5] != 'score':
                sys.exit(f"unexpected dump header in {part}: {hdr}")
            ch = f[0]
            if ch not in chrom_map: chrom_map[ch] = len(chrom_map)
            D['chrom'].append(chrom_map[ch]); D['start'].append(int(f[1])); D['end'].append(int(f[2]))
            c = f[3]
            if c not in cons_map: cons_map[c] = len(cons_map)
            D['consensus'].append(cons_map[c])
            D['strand'].append(1 if f[4] == '+' else 0)
            D['score'].append(int(f[5])); D['score_fwd'].append(int(f[6])); D['score_bwd'].append(int(f[7]))
            D['cons_fwd'].append(int(f[8])); D['cons_bwd'].append(int(f[9]))
            D['cons_len'].append(int(f[10])); D['gate_left'].append(int(f[11])); D['gate_right'].append(int(f[12]))
            D['gc'].append(float(f[13]))
    D['cons_names'] = [None] * len(cons_map)
    for c, i in cons_map.items(): D['cons_names'][i] = c
    D['chrom_names'] = [None] * len(chrom_map)
    for c, i in chrom_map.items(): D['chrom_names'][i] = c
    for k in COLS:
        if k in ('chrom', 'consensus'): continue
        D[k] = np.array(D[k], dtype=np.int64 if k != 'gc' else np.float64)
    D['chrom'] = np.array(D['chrom'], dtype=np.int32)
    D['consensus'] = np.array(D['consensus'], dtype=np.int32)
    return D

FEAT_NAMES = ['score', 'score_fwd', 'score_bwd', 'logL', 'L', 'score_per_L', 'cov',
              'cons_fwd', 'cons_bwd', 'gc', 'gate_left', 'gate_right',
              'gate_sum', 'gate_min', 'logclen', 'strand', 'cons']

def features(D):
    L = D['end'] - D['start']
    s = D['score'].astype(np.float32)
    cov = (D['cons_fwd'] + D['cons_bwd']) / np.maximum(1, D['cons_len'])
    X = np.column_stack([
        s, D['score_fwd'].astype(np.float32), D['score_bwd'].astype(np.float32),
        np.log(L).astype(np.float32), L.astype(np.float32), (s / L).astype(np.float32),
        cov.astype(np.float32), D['cons_fwd'].astype(np.float32), D['cons_bwd'].astype(np.float32),
        D['gc'].astype(np.float32), D['gate_left'].astype(np.float32), D['gate_right'].astype(np.float32),
        (D['gate_left'] + D['gate_right']).astype(np.float32),
        np.minimum(D['gate_left'], D['gate_right']).astype(np.float32),
        np.log(D['cons_len']).astype(np.float32), D['strand'].astype(np.float32)])
    return X, L

def true_overlap(merged, D):
    """Per-candidate TP bases against RM annotation, vectorized per chrom."""
    tp = np.zeros(len(D['start']), dtype=np.int64)
    for ci in np.unique(D['chrom']):
        c = D['chrom_names'][ci]
        m = D['chrom'] == ci
        if c in merged:
            tp[m] = coverage(merged, c, D['end'][m]) - coverage(merged, c, D['start'][m])
    return tp

def merged_kept(D, keep):
    """Merge kept intervals per chromosome -> list of (chrom, st, en) arrays."""
    out = []
    for ci in np.unique(D['chrom']):
        m = np.nonzero(keep & (D['chrom'] == ci))[0]
        if len(m) == 0: continue
        st = D['start'][m]; en = D['end'][m]
        o = np.argsort(st, kind='stable'); st, en = st[o], en[o]
        keepf = np.empty(len(st), dtype=bool); keepf[0] = True
        keepf[1:] = st[1:] > np.maximum.accumulate(en[:-1])
        mst = st[keepf]
        block_id = np.cumsum(keepf) - 1
        men = np.zeros(len(mst), dtype=np.int64)
        np.maximum.at(men, block_id, en)
        out.append((D['chrom_names'][ci], mst, men))
    return out

def evaluate(D, keep, merged, total):
    tp = kept = 0
    for c, st, en in merged_kept(D, keep):
        kept += int((en - st).sum())
        tp += int((coverage(merged, c, en) - coverage(merged, c, st)).sum())
    fp = kept - tp
    P = tp / max(1, tp + fp); R = tp / total; F = 2 * P * R / max(1e-9, P + R)
    return P, R, F

tr_m, tr_tot = load_rm(train_rm); te_m, te_tot = load_rm(test_rm)
train = load(train_p); test = load(test_p)
print(f"train {len(train['start']):,} rows, test {len(test['start']):,} rows", flush=True)

Xtr, Ltr = features(train); Xte, Lte = features(test)
tp_tr = true_overlap(tr_m, train)
y = (tp_tr * 2 >= Ltr).astype(np.int8); w = Ltr.astype(np.float32)
print(f"positive rate (base-weighted) {np.average(y, weights=w):.4f}", flush=True)

# per-consensus categorical id
ctr = train['consensus'].astype(np.int32)
cons_names = train['cons_names'] + [c for c in test['cons_names'] if c not in set(train['cons_names'])]
cmap = {c: i for i, c in enumerate(cons_names)}
cte = np.array([cmap[c] for c in test['cons_names']], dtype=np.int32)[test['consensus']]

Xtr2 = np.column_stack([Xtr, ctr]); Xte2 = np.column_stack([Xte, cte])
# Lightweight evaluation copies: keep only what evaluate() needs.
train_eval = {'chrom': train['chrom'], 'start': train['start'], 'end': train['end'],
              'chrom_names': train['chrom_names'], 'sc': train['score'],
              'cons': train['cons_names'], 'ci': train['consensus']}
test_eval  = {'chrom': test['chrom'], 'start': test['start'], 'end': test['end'],
              'chrom_names': test['chrom_names'], 'sc': test['score'],
              'cons': test['cons_names'], 'ci': test['consensus']}
# Drop the heavy original dicts; Xtr2/Xte2 are kept for prediction.
del Xtr, Xte, train, test, ctr, cte

n_val = max(1, len(y) // 20)
val = np.zeros(len(y), dtype=bool); val[-n_val:] = True
params = dict(objective='binary', metric='binary_logloss', learning_rate=0.03,
              num_leaves=127, min_data_in_leaf=100, feature_fraction=0.9,
              bagging_fraction=0.8, bagging_freq=1, seed=1, verbosity=-1,
              max_bin=63)
print("building datasets...", flush=True)
dtr = lgb.Dataset(Xtr2[~val], label=y[~val], weight=w[~val], feature_name=FEAT_NAMES, categorical_feature=['cons'])
dva = lgb.Dataset(Xtr2[val], label=y[val], weight=w[val], feature_name=FEAT_NAMES, categorical_feature=['cons'], reference=dtr)
# We need Xtr2/Xte2 for predictions, so don't delete them yet.  After the
# model is trained we release the Dataset objects (the booster doesn't need them).
model = lgb.train(params, dtr, num_boost_round=3000, valid_sets=[dva],
                  callbacks=[lgb.early_stopping(150), lgb.log_evaluation(100)])
try:
    model.free_dataset()
except Exception:
    pass
del dtr, dva, y, w, val
print("best iteration:", model.best_iteration)

ptr = model.predict(Xtr2, num_iteration=model.best_iteration)
pte = model.predict(Xte2, num_iteration=model.best_iteration)

# keep-threshold sweep on training rows
best = (0, 0.5)
for tau in [0.05 + 0.05 * i for i in range(18)]:
    F = evaluate(train_eval, ptr >= tau, tr_m, tr_tot)[2]
    if F > best[0]: best = (F, tau)
tau = best[1]
print(f"cutoff tau={tau} (train F1 {best[0]:.4f})")

# threshold file for the floor baselines (same as train_logistic.py)
import os
thr_p = sys.argv[6] if len(sys.argv) > 6 else None
if thr_p and os.path.exists(thr_p):
    th = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(thr_p)}
    th_arr = np.array([th.get(c, 30) for c in test_eval['cons']], dtype=np.int64)[test_eval['ci']]
else:
    th_arr = np.full(len(test_eval['start']), 30, dtype=np.int64)

lines = []
lines.append("offline on test chromosome (base-level vs RepeatMasker):")
for name, keep in [("floor 30, no model", test_eval['sc'] >= 30),
                   ("per-consensus thresholds", test_eval['sc'] >= th_arr),
                   ("logistic (chr2 model, ref)", None),
                   ("gbm + categorical consensus", pte >= tau)]:
    if keep is None:
        lines.append(f"  logistic (chr2 model, ref)          P 0.9828  R 0.8314  F1 0.9008  [from train_logistic.py]")
        continue
    P, R, F = evaluate(test_eval, keep, te_m, te_tot)
    lines.append(f"  {name:<34} P {P:.4f}  R {R:.4f}  F1 {F:.4f}")
report = "\n".join(lines)
print(report)
with open(out_rep, 'w') as f:
    f.write(report + "\n")
model.save_model(out_rep + '.lgb')
print("gbm saved to", out_rep + '.lgb')
#!/usr/bin/env python3
"""Apply a learned logistic model to a candidate dump offline and write a BED.

Supports both the single global `#logistic` model and the new per-family-group
`#family_models` format produced by `train_family_models.py`.

usage: apply_model_offline.py <model.tsv> <cand.tsv> <out.bed> [tau]
"""
import sys, math
from collections import defaultdict

model_p, cand_p, out_p = sys.argv[1:4]
tau_override = float(sys.argv[4]) if len(sys.argv) > 4 else None

def load_global(lines):
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
            try: offsets[f[0]] = float(f[1])
            except ValueError: pass
    tau = head[9]
    return (head[:8], head[8], mu, sd, offsets, tau)

def load_family(lines):
    """Returns (models, fallback_thr, cons_to_group).
    models: dict group_name -> (w, b0, mu, sd, offsets_dict, tau)
    """
    it = iter(lines)
    header = next(it).split('\t')
    n_groups = int(header[1]) if len(header) > 1 else 0
    models = {}
    cons_to_group = {}
    fallback_thr = {}
    nums = lambda l: [float(x) for x in l.split('\t')[1:]]
    in_fallback = False
    cur_group = None
    cur_model = None
    cur_cons = []
    for l in it:
        l = l.rstrip('\n')
        if not l:
            continue
        if l.startswith("#fallback_thresholds"):
            in_fallback = True
            cur_group = None
            continue
        if l.startswith("#group"):
            in_fallback = False
            f = l.split('\t')
            cur_group = f[1]
            if len(f) >= 4 and f[2] != "fallback":
                # model group: #group name n_cons tau
                tau = float(f[3])
                # next lines are #logistic, #mu, #sd
                wv = nums(next(it))
                mu = nums(next(it))
                sd = nums(next(it))
                if len(wv) != 9 or len(mu) != 8 or len(sd) != 8:
                    raise ValueError("bad family model header")
                cur_model = [wv[:8], wv[8], mu, sd, {}, tau]
            else:
                cur_model = None
            cur_cons = []
            continue
        f = l.split('\t')
        if len(f) < 2:
            continue
        name, val = f[0], f[1]
        if in_fallback:
            try: fallback_thr[name] = int(val)
            except ValueError: pass
        else:
            cons_to_group[name] = cur_group
            if cur_model is not None:
                try: cur_model[4][name] = float(val)
                except ValueError: pass
    # flush last group
    if cur_group and cur_model:
        models[cur_group] = tuple(cur_model)
    return models, fallback_thr, cons_to_group

lines = open(model_p).readlines()
family_mode = lines[0].startswith("#family_models")
if family_mode:
    models, fallback_thr, cons_to_group = load_family(lines)
    def prob_for(row):
        cons = row['consensus']
        g = cons_to_group.get(cons)
        if g is None or g not in models:
            return None
        w, b0, mu, sd, offsets, tau = models[g]
        isd = [1.0 / s if s else 1.0 for s in sd]
        L = int(row['end']) - int(row['start'])
        s = int(row['score']); gl = int(row['gate_left']); gr = int(row['gate_right'])
        cov = (int(row['cons_fwd']) + int(row['cons_bwd'])) / max(1, int(row['cons_len']))
        x = [s, math.log(L), s / L, cov, float(row['gc']), gl + gr, min(gl, gr), math.log(int(row['cons_len']))]
        zz = b0 + offsets.get(cons, 0.0)
        for j in range(8):
            zz += w[j] * (x[j] - mu[j]) * isd[j]
        zz = max(-30, min(30, zz))
        return 1 / (1 + math.exp(-zz)), tau
else:
    w, b0, mu, sd, offsets, tau = load_global(lines)
    tau = tau_override if tau_override is not None else tau
    isd = [1.0 / s if s else 1.0 for s in sd]
    def prob_for(row):
        cons = row['consensus']
        L = int(row['end']) - int(row['start'])
        s = int(row['score']); gl = int(row['gate_left']); gr = int(row['gate_right'])
        cov = (int(row['cons_fwd']) + int(row['cons_bwd'])) / max(1, int(row['cons_len']))
        x = [s, math.log(L), s / L, cov, float(row['gc']), gl + gr, min(gl, gr), math.log(int(row['cons_len']))]
        zz = b0 + offsets.get(cons, 0.0)
        for j in range(8):
            zz += w[j] * (x[j] - mu[j]) * isd[j]
        zz = max(-30, min(30, zz))
        return 1 / (1 + math.exp(-zz)), tau

by_chrom = defaultdict(list)
hdr = None
for l in open(cand_p):
    f = l.rstrip('\n').split('\t')
    if hdr is None or f[0] == 'chrom':
        hdr = f; continue
    r = dict(zip(hdr, f))
    res = prob_for(r)
    if res is None:
        # family mode: consensus not in any model group; try fallback threshold
        t = fallback_thr.get(r['consensus'], 30)
        if int(r['score']) >= t:
            by_chrom[r['chrom']].append((int(r['start']), int(r['end']), r['consensus']))
    else:
        p, tau = res
        if p >= tau:
            by_chrom[r['chrom']].append((int(r['start']), int(r['end']), r['consensus']))

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
print(f"wrote {out_p}")

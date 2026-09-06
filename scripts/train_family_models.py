#!/usr/bin/env python3
"""Train per-family-group logistic models for library-alignment hits.

Memory-light streaming version: splits the training candidate dump into one
per-group file, then trains each group model independently. Writes a combined
model that krep can apply per consensus.

usage: train_family_models.py <train.tsv[,...]> <train_rm.bed[,...]> \
                              <labels.tsv> <test.tsv> <test_rm.bed> <out_model.tsv>
"""
import sys, os, bisect, math, random, tempfile
from collections import defaultdict

train_p, train_rm, labels_p, test_p, test_rm, out_p = sys.argv[1:7]

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
group_names = [g[0] for g in GROUPS]

def consensus_group(name_to_fam):
    out = {}
    for cons, fam in name_to_fam.items():
        for g, members in GROUPS:
            if fam in members:
                out[cons] = g
                break
        else:
            out[cons] = 'Other'
    return out

name_to_fam = {}
for l in open(labels_p):
    if l.startswith('#'):
        continue
    f = l.rstrip('\n').split('\t')
    if len(f) >= 2:
        name_to_fam[f[0]] = f[1]
cons_to_group = consensus_group(name_to_fam)
print(f"loaded {len(name_to_fam)} consensus labels")
for g in sorted(set(cons_to_group.values())):
    n = sum(1 for c, gr in cons_to_group.items() if gr == g)
    print(f"  {g}: {n} consensi")

NF = 8
random.seed(1)

def load_rm(paths):
    iv = defaultdict(list)
    for part in paths.split(','):
        for l in open(part):
            f = l.split()
            if len(f) >= 3:
                iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}; total = 0
    for c, v in iv.items():
        v.sort(); st, en = [], []
        for b, e in v:
            if len(en) and b <= en[-1]: en[-1] = max(en[-1], e)
            else: st.append(b); en.append(e)
        merged[c] = (st, en)
        total += sum(en[i] - st[i] for i in range(len(st)))
    return merged, total

def overlap(merged, c, b, e):
    if c not in merged:
        return 0
    st, en = merged[c]
    i = max(0, bisect.bisect_left(st, b) - 1); t = 0
    while i < len(st) and st[i] < e:
        t += max(0, min(e, en[i]) - max(b, st[i])); i += 1
    return t

def make_features(r):
    b = int(r['start']); e = int(r['end']); L = e - b
    s = int(r['score']); gl = int(r['gate_left']); gr = int(r['gate_right'])
    cov = (int(r['cons_fwd']) + int(r['cons_bwd'])) / max(1, int(r['cons_len']))
    glen = int(r['cons_len'])
    return [s, math.log(L), s / L, cov, float(r['gc']), gl + gr, min(gl, gr), math.log(glen)], L

# Split train candidates into per-group temp files
print("splitting training candidates by group ...")
tmpdir = tempfile.mkdtemp(prefix="krep_family_")
group_paths = {g: os.path.join(tmpdir, f"{g.replace('/', '_')}.tsv") for g in group_names}
group_files = {g: open(p, 'w') for g, p in group_paths.items()}
group_counts = {g: 0 for g in group_names}
hdr = None
for part in train_p.split(','):
    for l in open(part):
        f = l.rstrip('\n').split('\t')
        if hdr is None or f[0] == 'chrom':
            hdr = f
            for gf in group_files.values(): gf.write(l)
            continue
        r = dict(zip(hdr, f))
        g = cons_to_group.get(r['consensus'], 'Other')
        group_files[g].write(l)
        group_counts[g] += 1
for gf in group_files.values(): gf.close()
print("group counts:", {g: group_counts[g] for g in group_names})

tr_m, tr_tot = load_rm(train_rm)
te_m, te_tot = load_rm(test_rm)

MAX_PER_GROUP = 120000

def load_group_file(path, merged):
    """Load rows for one group, returning list of (consensus, chrom, b, e, L, glen, y, x)."""
    rows = []
    hdr_local = None
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if hdr_local is None or f[0] == 'chrom':
            hdr_local = f; continue
        r = dict(zip(hdr_local, f))
        x, L = make_features(r)
        y = 1 if overlap(merged, r['chrom'], int(r['start']), int(r['end'])) * 2 >= L else 0
        glen = int(r['cons_len'])
        rows.append((r['consensus'], r['chrom'], int(r['start']), int(r['end']), L, glen, y, x))
    return rows

def sample(rows, max_n):
    if len(rows) <= max_n:
        return rows
    random.seed(7)
    return random.sample(rows, max_n)

def train_group(gname, rows, min_rows=200, min_pos_bases=5000):
    rows = sample(rows, MAX_PER_GROUP)
    n = len(rows)
    pos_bases = sum(r[4] for r in rows if r[6])
    if n < min_rows or pos_bases < min_pos_bases:
        print(f"  {gname}: too few data ({n} rows, {pos_bases} pos bases); using fallback")
        return None
    print(f"  {gname}: training on {n:,} rows ({pos_bases:,} pos bases)")
    cons_set = sorted({r[0] for r in rows})
    cidx = {c: i for i, c in enumerate(cons_set)}
    ci = [cidx[r[0]] for r in rows]

    mu = [0.0] * NF; sd = [0.0] * NF
    for j in range(NF):
        s1 = s2 = 0.0
        for r in rows:
            v = r[7][j]; s1 += v; s2 += v * v
        mu[j] = s1 / n
        sd[j] = math.sqrt(max(0.0, s2 / n - mu[j] * mu[j]))
        if sd[j] == 0: sd[j] = 1.0
    isd = [1.0 / s for s in sd]

    Xs = [[(r[7][j] - mu[j]) * isd[j] for j in range(NF)] for r in rows]
    y = [r[6] for r in rows]
    L = [r[4] for r in rows]
    meanL = sum(L) / n
    sw = [li / meanL for li in L]

    n_cons = len(cons_set)
    w = [0.0] * NF; off = [0.0] * n_cons; b0 = 0.0
    lr, lam, epochs = 0.5, 1e-4, 15
    counts = [0] * n_cons
    for c in ci: counts[c] += 1
    counts = [max(1, c) for c in counts]

    order = list(range(n))
    for ep in range(epochs):
        random.shuffle(order)
        dw = [0.0] * NF; doff = [0.0] * n_cons; db0 = 0.0
        loss = 0.0
        for i in order:
            z = b0 + off[ci[i]]
            xi = Xs[i]
            for j in range(NF): z += w[j] * xi[j]
            if z > 30: p, logp, log1p = 1.0, 0.0, -30.0
            elif z < -30: p, logp, log1p = 0.0, -30.0, 0.0
            else:
                ez = math.exp(-z)
                p = 1.0 / (1.0 + ez)
                logp = math.log(p + 1e-12)
                log1p = math.log(1 - p + 1e-12)
            swi = sw[i]
            loss -= swi * (y[i] * logp + (1 - y[i]) * log1p)
            g = (p - y[i]) * swi
            db0 += g
            doff[ci[i]] += g
            for j in range(NF): dw[j] += g * xi[j]
        scale = lr / n
        for j in range(NF): w[j] -= scale * dw[j] + lr * lam * w[j]
        for c in range(n_cons): off[c] -= scale * doff[c] / counts[c] + lr * 10 * lam * off[c]
        b0 -= scale * db0
        if ep % 3 == 0:
            print(f"    epoch {ep+1}: loss {loss / n:.4f}")
    print(f"    final loss {loss / n:.4f}, weights: " + ", ".join(f"{n}={v:.3f}" for n, v in zip(
        ['score','logL','score/L','cov','gc','gate_sum','gate_min','logclen'], w)))
    return {'cons': cons_set, 'w': w, 'b0': b0, 'off': off, 'mu': mu, 'sd': sd,
            'rows': rows, 'ci': ci, 'Xs': Xs, 'y': y, 'L': L}

def pick_tau(model, merged, total_bases, floor=0.85):
    w = model['w']; b0 = model['b0']; off = model['off']; mu = model['mu']; sd = model['sd']
    isd = [1.0 / s for s in sd]
    cons_map = {c: i for i, c in enumerate(model['cons'])}
    probs = []
    for r in model['rows']:
        c = r[0]
        z = b0 + off[cons_map[c]]
        for j in range(NF):
            z += w[j] * (r[7][j] - mu[j]) * isd[j]
        z = max(-30.0, min(30.0, z))
        p = 1.0 / (1.0 + math.exp(-z))
        probs.append(p)
    best = (0, 0.5, 0.0, 0.0)
    for t in [i / 40 for i in range(2, 39)]:
        by_c = defaultdict(list)
        for i, r in enumerate(model['rows']):
            if probs[i] >= t:
                by_c[r[1]].append((r[2], r[3]))
        tp = kept = 0
        for c, iv in by_c.items():
            iv.sort(); cs, ce = iv[0]
            for b, e in iv[1:]:
                if b > ce:
                    kept += ce - cs; tp += overlap(merged, c, cs, ce); cs, ce = b, e
                else: ce = max(ce, e)
            kept += ce - cs; tp += overlap(merged, c, cs, ce)
        P = tp / kept if kept else 0
        R = tp / total_bases
        F = 2 * P * R / max(1e-9, P + R)
        obj = F if P >= floor else F * (P / floor)
        if obj > best[0]:
            best = (obj, t, P, R)
    return best[1], best[2], best[3]

models = {}
for gname in group_names:
    path = group_paths[gname]
    rows = load_group_file(path, tr_m)
    print(f"\nGroup {gname}: {len(rows):,} train rows")
    floor = 0.50 if gname == 'Ancient' else 0.85
    m = train_group(gname, rows)
    if m:
        tau, P, R = pick_tau(m, tr_m, tr_tot, floor=floor)
        m['tau'] = tau
        print(f"  chosen tau={tau:.2f} (train P={P:.4f} R={R:.4f})")
        models[gname] = m
    else:
        models[gname] = None

# Fallback per-consensus thresholds
threshold_path = '/mnt/c/krep_work/model_hybrid_pass2f_thresholds.tsv'
fallback_thr = {}
for l in open(threshold_path):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 2:
        try: fallback_thr[f[0]] = int(f[1])
        except ValueError: pass

# Write combined family model file (format consumed by krep align.rs)
with open(out_p, 'w') as f:
    f.write(f"#family_models\t{len(group_names)}\n")
    for gname in group_names:
        m = models[gname]
        cons_in_group = [c for c, gr in cons_to_group.items() if gr == gname]
        if m:
            f.write(f"#group\t{gname}\t{len(m['cons'])}\t{m['tau']:.4f}\n")
            f.write("#logistic\t" + "\t".join(f"{v:.6f}" for v in m['w']) + f"\t{m['b0']:.6f}\n")
            f.write("#mu\t" + "\t".join(f"{v:.6f}" for v in m['mu']) + "\n")
            f.write("#sd\t" + "\t".join(f"{v:.6f}" for v in m['sd']) + "\n")
            off_map = dict(zip(m['cons'], m['off']))
            for c in m['cons']:
                f.write(f"{c}\t{off_map.get(c, 0.0):.6f}\n")
        else:
            f.write(f"#group\t{gname}\tfallback\n")
            for c in cons_in_group:
                f.write(f"{c}\t0.0\n")
    f.write("#fallback_thresholds\n")
    for c in sorted(fallback_thr.keys()):
        f.write(f"{c}\t{fallback_thr[c]}\n")
print("model written to", out_p)

# Offline evaluation: stream test file
def evaluate(test_path, merged, total_bases, use_family):
    by_c = defaultdict(list)
    hdr_local = None
    for l in open(test_path):
        f = l.rstrip('\n').split('\t')
        if hdr_local is None or f[0] == 'chrom':
            hdr_local = f; continue
        r = dict(zip(hdr_local, f))
        c = r['consensus']
        chrom = r['chrom']
        b = int(r['start']); e = int(r['end']); L = e - b
        accept = False
        if use_family:
            gname = cons_to_group.get(c, 'Other')
            m = models[gname]
            if m and c in m['cons']:
                idx = m['cons'].index(c)
                x, _ = make_features(r)
                z = m['b0'] + m['off'][idx]
                for j in range(NF):
                    z += m['w'][j] * (x[j] - m['mu'][j]) / m['sd'][j]
                z = max(-30.0, min(30.0, z))
                if 1 / (1 + math.exp(-z)) >= m['tau']:
                    accept = True
            if not accept:
                t = fallback_thr.get(c, 30)
                if int(r['score']) >= t:
                    accept = True
        else:
            t = fallback_thr.get(c, 30)
            if L >= t:
                accept = True
        if accept:
            by_c[chrom].append((b, e))
    tp = kept = 0
    for c, iv in by_c.items():
        iv.sort(); cs, ce = iv[0]
        for b, e in iv[1:]:
            if b > ce:
                kept += ce - cs; tp += overlap(merged, c, cs, ce); cs, ce = b, e
            else: ce = max(ce, e)
        kept += ce - cs; tp += overlap(merged, c, cs, ce)
    P = tp / kept if kept else 0
    R = tp / total_bases
    F = 2 * P * R / max(1e-9, P + R)
    return P, R, F

print("\noffline test chromosome evaluation:")
P, R, F = evaluate(test_p, te_m, te_tot, use_family=False)
print(f"  per-consensus thresholds  P {P:.4f}  R {R:.4f}  F1 {F:.4f}")
P, R, F = evaluate(test_p, te_m, te_tot, use_family=True)
print(f"  per-family models         P {P:.4f}  R {R:.4f}  F1 {F:.4f}")

# Cleanup temp files
for p in group_paths.values():
    try: os.remove(p)
    except OSError: pass
try: os.rmdir(tmpdir)
except OSError: pass

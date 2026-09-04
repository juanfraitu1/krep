"""Learn per-consensus minimum scores from a candidate dump labelled by RepeatMasker.

usage: train_thresholds.py <cand.tsv> <rm_family.bed> <out_model.tsv> [min_cands=20]
Labels each candidate by the fraction of its bases inside RepeatMasker annotation.
For each consensus, hits sorted by score descending; the threshold is the score at
which the cumulative (true bases - false bases) is maximal. Consensi with few
candidates fall back to the global optimum.
"""
import sys, bisect
from collections import defaultdict

cand_path, rm_path, out_path = sys.argv[1:4]
min_cands = int(sys.argv[4]) if len(sys.argv) > 4 else 20

# --- RepeatMasker intervals per chrom, merged ---
rm = defaultdict(list)
for l in open(rm_path):
    f = l.split()
    if len(f) >= 3: rm[f[0]].append((int(f[1]), int(f[2])))
merged = {}
for c, iv in rm.items():
    iv.sort(); out = []
    for b, e in iv:
        if out and b <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], e))
        else: out.append((b, e))
    merged[c] = (out, [b for b, _ in out])
def overlap(c, b, e):
    if c not in merged: return 0
    iv, starts = merged[c]
    i = max(0, bisect.bisect_left(starts, b) - 1); t = 0
    while i < len(iv) and iv[i][0] < e:
        t += max(0, min(e, iv[i][1]) - max(b, iv[i][0])); i += 1
    return t

# --- candidates ---
per = defaultdict(list); hdr = None; n = 0
for l in open(cand_path):
    f = l.rstrip('\n').split('\t')
    if hdr is None: hdr = f; continue
    row = dict(zip(hdr, f)); b, e = int(row['start']), int(row['end']); L = e - b
    tp = overlap(row['chrom'], b, e); per[row['consensus']].append((int(row['score']), tp, L - tp)); n += 1
print(f"{n:,} candidates, {len(per)} consensi")

def best_threshold(rows):
    rows = sorted(rows, key=lambda r: -r[0]); cum = 0; best = 0; best_t = None
    i = 0
    while i < len(rows):
        s = rows[i][0]; j = i
        while j < len(rows) and rows[j][0] == s:
            cum += rows[j][1] - rows[j][2]; j += 1
        if cum > best: best, best_t = cum, s
        i = j
    return best_t, best  # None => keep nothing

all_rows = [r for rows in per.values() for r in rows]
g_t, _ = best_threshold(all_rows)
print(f"global optimum threshold: {g_t}")
kept = 0
with open(out_path, 'w') as out:
    for name, rows in per.items():
        if len(rows) < min_cands:
            t = g_t
        else:
            t, _ = best_threshold(rows)
            if t is None: t = 10**6  # nothing from this consensus is worth keeping
        out.write(f"{name}\t{t}\n"); kept += 1
print(f"wrote {kept} thresholds to {out_path}")
# summary of what the thresholds imply on the training set
tp = fp = fn_kept = 0
th = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(out_path)}
for name, rows in per.items():
    for s, t_, f_ in rows:
        if s >= th[name]: tp += t_; fp += f_
print(f"training-set hits kept: true bases {tp:,}, false bases {fp:,}, precision {tp/(tp+fp):.4f}")

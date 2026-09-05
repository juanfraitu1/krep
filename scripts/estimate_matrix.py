"""Estimate a divergence-aware scoring matrix from krep's gate-window
substitution counts (produced by `--lib-subst-dump`).

The counts are 4x4 (consensus x genome) over the ungapped flanking windows
around accepted library hits. From them we compute relative scores that keep
the expected random-match score per base roughly at the current -0.5
(match=+1 / mismatch=-1) while up-weighting transitions and down-weighting
transversions.
"""
import sys, math
path = sys.argv[1]
BASES = 'ACGT'
idx = {b: i for i, b in enumerate(BASES)}
counts = [[0]*4 for _ in range(4)]
n = 0
hdr = None
for l in open(path):
    f = l.rstrip('\n').split('\t')
    if hdr is None or f[0] == 'chrom':
        hdr = f
        continue
    for i, b in enumerate(BASES):
        for j, b2 in enumerate(BASES):
            counts[i][j] += int(f[6 + i*4 + j])
            n += int(f[6 + i*4 + j])
print(f"total counted bases: {n:,}")
print("    A       C       G       T")
for i, b in enumerate(BASES):
    print(f"{b} " + "  ".join(f"{counts[i][j]:7,}" for j in range(4)))

# Symmetrise by grouping complementary substitutions (A<->G and C<->T
# transitions; the four transversion pairs). This makes the matrix strand-
# invariant and less noisy.
trans = counts[idx['A']][idx['G']] + counts[idx['G']][idx['A']] + \
        counts[idx['C']][idx['T']] + counts[idx['T']][idx['C']]
tv = 0
for i in range(4):
    for j in range(4):
        if i != j:
            if not ((i==0 and j==2) or (i==2 and j==0) or (i==1 and j==3) or (i==3 and j==1)):
                tv += counts[i][j]
match = sum(counts[i][i] for i in range(4))
print(f"\nmatches: {match:,}  transitions: {trans:,}  transversions: {tv:,}")
print(f"transition/transversion ratio = {trans/tv:.3f}")

# Derive scores such that E_random_per_base stays at -0.5.
# Under equal base composition:
#   0.25*s_match + 0.5*p_trans*s_trans + 0.25*p_tv*s_tv = -0.5
# where p_trans = fraction of mismatches that are transitions, p_tv = 1-p_trans.
# We pin s_match = 1 and solve for s_trans, s_tv given a chosen ratio of
# penalty increase between tv and trans.
# Simple target: transversions should be penalised ~1.5x as much as transitions.
# So s_tv = s_trans - delta, choose delta so expected score = -0.5.
p_m = 0.25
p_trans = trans / (trans + tv)  # conditional on mismatch
p_tv = 1.0 - p_trans
# We want ratio: |s_tv| / |s_trans| = r (r > 1). With s_match=1:
# s_trans = x (<0), s_tv = r*x (<0)
# E = p_m*1 + (1-p_m)*(p_trans*x + p_tv*r*x) = -0.5
r = 1.5
x = (-0.5 - p_m) / ((1 - p_m) * (p_trans + p_tv * r))
s_trans = x
s_tv = r * x
print(f"\nproposed matrix (match=+1, E_random_per_base = -0.5):")
print(f"  transition score {s_trans:.3f}, transversion score {s_tv:.3f}")

# Print full 4x4 with this rule.
print("\n   A      C      G      T")
for i, b in enumerate(BASES):
    row = []
    for j in range(4):
        if i == j:
            row.append(1.0)
        elif (i==0 and j==2) or (i==2 and j==0) or (i==1 and j==3) or (i==3 and j==1):
            row.append(s_trans)
        else:
            row.append(s_tv)
    print(f"{b} " + "  ".join(f"{v:6.3f}" for v in row))

# Also log-odds style (scaled) for reference.
print("\n--- log-odds (for reference; not used below) ---")
total = sum(sum(r) for r in counts)
q_match = match / total
q_trans = trans / total
q_tv = tv / total
print(f"observed frequencies: match {q_match:.4f}, trans {q_trans:.4f}, tv {q_tv:.4f}")
# arbitrary background: match 0.25, each mismatch type equally likely.
print(f"log-odds transition {math.log(q_trans/0.375):.3f}, transversion {math.log(q_tv/0.375):.3f}")

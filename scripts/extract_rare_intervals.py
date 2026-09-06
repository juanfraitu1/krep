#!/usr/bin/env python3
"""Extract genomic intervals for rare repeat families and build mini-consensi.

For each target family, reads the family-specific RM BED, extracts the
sequences, sorts by length, greedily assembles a representative consensus
using a simple majority-vote multiple alignment of the longest copies.

usage: extract_rare_intervals.py <genome.fa> <rm_family.bed> <out_dir> [families]
"""
import sys, os, re
from collections import defaultdict

genome_p, bed_p, out_dir = sys.argv[1:4]
DEFAULT_FAMILIES = ['LINE/CR1', 'RC/Helitron', 'SINE/MIR', 'LINE/L2', 'DNA/hAT-Tip100']
families = sys.argv[4].split(',') if len(sys.argv) > 4 else DEFAULT_FAMILIES

os.makedirs(out_dir, exist_ok=True)

# load genome
seqs = {}
cur = None
for l in open(genome_p):
    if l.startswith('>'):
        cur = l[1:].split()[0]
        seqs[cur] = []
    else:
        seqs[cur].append(l.strip().upper())
seqs = {k: ''.join(v) for k, v in seqs.items()}

# load BED: chrom start end family
by_fam = defaultdict(list)
for l in open(bed_p):
    if l.startswith('#'):
        continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 4:
        continue
    chrom, b, e, fam = f[0], int(f[1]), int(f[2]), f[3]
    by_fam[fam].append((chrom, b, e))

def consensus_from_copies(copies, min_copies=3, max_copies=30):
    """Greedy majority-vote consensus from sorted copies."""
    copies = sorted(copies, key=lambda s: -len(s))[:max_copies]
    n = len(copies)
    if n < min_copies:
        return None
    # Use the longest as anchor
    anchor = copies[0]
    aln = []
    for s in copies[1:]:
        # simple banded alignment by identity around diagonal; skip for speed
        # just keep sequences that overlap anchor by prefix/suffix roughly
        if len(s) < 50:
            continue
        aln.append(s)
    if not aln:
        return None
    # column-wise majority over copies aligned to anchor by local ungapped match
    # For speed: build a simple profile from copies anchored at best offset to anchor
    offsets = [0]
    for s in aln:
        best_o = 0
        best_id = -1
        max_o = min(20, len(anchor) - 50)
        for o in range(-max_o, max_o + 1):
            idy = 0
            for i in range(max(0, -o), min(len(s), len(anchor) - o)):
                if s[i] == anchor[i + o]:
                    idy += 1
            if idy > best_id:
                best_id = idy
                best_o = o
        offsets.append(best_o)
    # build consensus by majority vote in anchor coordinate
    out = []
    for i in range(len(anchor)):
        counts = defaultdict(int)
        for idx, s in enumerate([anchor] + aln):
            o = offsets[idx]
            j = i - o
            if 0 <= j < len(s):
                counts[s[j]] += 1
        if counts:
            out.append(max(counts, key=counts.get))
    cons = ''.join(out)
    # trim N-heavy edges
    cons = re.sub(r'^N+|N+$', '', cons)
    if len(cons) < 100:
        return None
    return cons

for fam in families:
    ivs = by_fam.get(fam, [])
    print(f"{fam}: {len(ivs)} intervals")
    copies = []
    for chrom, b, e in ivs:
        seq = seqs.get(chrom)
        if seq is None:
            continue
        s = seq[b:e]
        # skip very short or N-heavy
        if len(s) < 100:
            continue
        if s.count('N') / len(s) > 0.3:
            continue
        copies.append(s)
    if not copies:
        continue
    cons = consensus_from_copies(copies)
    if cons is None:
        print(f"  {fam}: not enough/long copies for consensus")
        continue
    safe = fam.replace('/', '_')
    out_p = os.path.join(out_dir, f"{safe}.fa")
    with open(out_p, 'w') as f:
        f.write(f">{fam}_denovo length={len(cons)} copies={len(copies)}\n")
        for i in range(0, len(cons), 80):
            f.write(cons[i:i+80] + '\n')
    print(f"  {fam}: wrote {len(cons)} bp from {len(copies)} copies -> {out_p}")

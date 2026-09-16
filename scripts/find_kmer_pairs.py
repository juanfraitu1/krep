"""Find all positions of a set of k-mers in a FASTA file and emit pairs.

usage: find_kmer_pairs.py <kmers.txt> <genome.fa> <k> <out_pairs.tsv>

kmers.txt format: one canonical k-mer per line (first whitespace-separated token).
out_pairs.tsv format: kmer\tpos1\tpos2 (one line per k-mer with exactly 2 hits)
"""
import sys
from collections import defaultdict

kmers_path, fa_path, k_str, out_path = sys.argv[1:5]
k = int(k_str)

target = set()
with open(kmers_path) as f:
    for line in f:
        seq = line.split()[0].upper()
        if len(seq) == k:
            target.add(seq)
print(f"loaded {len(target)} target k-mers", file=sys.stderr)

hits = defaultdict(list)
seq_parts = []
name = None

def process(name, parts):
    if not name or not parts:
        return
    seq = ''.join(parts).upper()
    n = len(seq)
    for i in range(n - k + 1):
        mer = seq[i:i+k]
        if 'N' in mer:
            continue
        if mer in target:
            hits[mer].append((name, i))
    print(f"scanned {name}: {n} bp, {sum(1 for v in hits.values() if v)}/{len(target)} kmers hit so far", file=sys.stderr)

with open(fa_path) as f:
    for line in f:
        if line.startswith('>'):
            if name:
                process(name, seq_parts)
            name = line[1:].strip().split()[0]
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if name:
        process(name, seq_parts)

print(f"total distinct kmers with hits: {len(hits)}", file=sys.stderr)

with open(out_path, 'w') as out:
    for mer, lst in hits.items():
        if len(lst) == 2:
            out.write(f"{mer}\t{lst[0][1]}\t{lst[1][1]}\n")
print("done", file=sys.stderr)

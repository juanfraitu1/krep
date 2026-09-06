#!/usr/bin/env python3
"""Combine existing de-novo consensus labels with labels for added rare consensi.

usage: make_augmented_labels.py <existing_labels.tsv> <added_fa[,added_fa2,...]> <out.tsv> [label]
"""
import sys, re

existing_labels = sys.argv[1]
added_paths = sys.argv[2]
out = sys.argv[3]
forced_label = sys.argv[4] if len(sys.argv) > 4 else None

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

def family_of_name(name):
    base = name.split('#')[0]
    # explicit de novo headers like LINE/CR1_denovo
    if base.endswith('_denovo'):
        return base.replace('_denovo', '').replace('_', '/')
    # variants keep original family; infer from name
    if '_var' in base:
        base = base.split('_var')[0]
    if re.search(r'\bAlu\b', base): return 'SINE/Alu'
    if re.search(r'\bL1\b', base): return 'LINE/L1'
    if 'ERVK' in base: return 'LTR/ERVK'
    if 'ERV1' in base: return 'LTR/ERV1'
    if 'ERVL-MaLR' in base: return 'LTR/ERVL-MaLR'
    if 'ERVL' in base: return 'LTR/ERVL'
    if re.search(r'TcMar|Tigger|Tc1|Tc2|Mariner|Pogo', base): return 'DNA/TcMar'
    if 'Tip100' in base or 'Charlie' in base or 'Blackjack' in base or 'hAT' in base:
        return 'DNA/hAT'
    if 'Satellite' in base: return 'Satellite'
    if 'Simple_repeat' in base or 'Low_complexity' in base: return 'Simple/Low'
    if 'MIR' in base: return 'SINE/MIR'
    if 'L2' in base: return 'LINE/L2'
    if 'CR1' in base: return 'LINE/CR1'
    if 'Helitron' in base: return 'RC/Helitron'
    return forced_label or 'Unknown'

with open(out, 'w') as f:
    # copy existing
    for l in open(existing_labels):
        if l.startswith('#'):
            f.write(l)
            continue
        f.write(l)
    # add new fasta entries
    for fa_path in added_paths.split(','):
        name = None; seqlen = 0
        for l in open(fa_path):
            if l.startswith('>'):
                if name is not None:
                    fam = family_of_name(name)
                    f.write(f"{name}\t{fam}\t0\t{seqlen}\t0.0000\n")
                name = l[1:].split()[0]
                seqlen = 0
            else:
                seqlen += len(l.strip())
        if name is not None:
            fam = family_of_name(name)
            f.write(f"{name}\t{fam}\t0\t{seqlen}\t0.0000\n")
print("wrote", out)

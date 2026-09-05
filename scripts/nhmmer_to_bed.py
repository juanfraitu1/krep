"""Convert nhmmer --tblout to a 4-column BED for krep evaluate.

Columns in tblout (HMMER 3.4):
  target name, target accession, query name, query accession, hmmfrom,
  hmmto, alifrom, alito, envfrom, envto, sq len, strand, Evalue, score,
  bias, description
"""
import sys, argparse

# map query name to family class used in RM bed (e.g. MIR -> SINE/MIR,
# L2 -> LINE/L2, CR1* -> LINE/CR1)
family_map = {
    'MIR': 'SINE/MIR',
    'L2': 'LINE/L2',
}

def family(name):
    if name.startswith('CR1'):
        return 'LINE/CR1'
    return family_map.get(name, name)

def parse_args():
    p = argparse.ArgumentParser(description='Convert nhmmer tblout to BED')
    p.add_argument('inp')
    p.add_argument('out')
    p.add_argument('--min-score', type=float, default=None,
                   help='minimum nhmmer bit score (inclusive)')
    p.add_argument('--max-evalue', type=float, default=None,
                   help='maximum E-value (inclusive)')
    return p.parse_args()

def main():
    args = parse_args()
    kept = 0
    with open(args.inp) as f, open(args.out, 'w') as w:
        for l in f:
            if l.startswith('#') or not l.strip():
                continue
            # nhmmer tbl is whitespace-delimited; description may contain spaces at end
            p = l.rstrip('\n').split()
            if len(p) < 16:
                continue
            try:
                chrom = p[0]
                qname = p[2]
                alifrom = int(p[6])
                alito = int(p[7])
                evalue = float(p[13])
                score = float(p[14])
            except (ValueError, IndexError):
                continue
            if args.max_evalue is not None and evalue > args.max_evalue:
                continue
            if args.min_score is not None and score < args.min_score:
                continue
            st, en = (alifrom, alito) if alifrom <= alito else (alito, alifrom)
            w.write(f"{chrom}\t{st}\t{en}\t{family(qname)}\n")
            kept += 1
    print(f"Wrote {kept} hits to {args.out}")

if __name__ == '__main__':
    main()

"""Convert nhmmer --tblout to a 4-column BED for krep evaluate.

Columns in tblout (HMMER 3.4):
  target name, target accession, query name, query accession, hmmfrom,
  hmmto, alifrom, alito, envfrom, envto, sq len, strand, Evalue, score,
  bias, description
"""
import sys, csv
inp, out = sys.argv[1:3]
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

with open(inp) as f, open(out, 'w') as w:
    for l in f:
        if l.startswith('#') or not l.strip():
            continue
        # nhmmer tbl is whitespace-delimited; description may contain spaces at end
        p = l.rstrip('\n').split()
        if len(p) < 16:
            continue
        # fields after fixed columns can be variable; parse from end
        # last 3 numeric: Evalue score bias
        # before that strand
        # before that sq len
        # before that envto envfrom alito alifrom hmmto hmmfrom
        try:
            chrom = p[0]
            qname = p[2]
            alifrom = int(p[6])
            alito = int(p[7])
        except (ValueError, IndexError):
            continue
        st, en = (alifrom, alito) if alifrom <= alito else (alito, alifrom)
        w.write(f"{chrom}\t{st}\t{en}\t{family(qname)}\n")

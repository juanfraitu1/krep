import sys
with open('chm13v2.0_SD.bed') as f, open('chr1_SD.bed', 'w') as out:
    total = 0
    for line in f:
        cols = line.rstrip('\n').split('\t')
        if cols[0] == 'chr1':
            out.write(f"NC_060925.1\t{cols[1]}\t{cols[2]}\n")
            total += int(cols[2]) - int(cols[1])
    print(f"wrote chr1_SD.bed with {total:,} bp")

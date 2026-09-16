import sys, statistics
dists = []
for line in sys.stdin:
    f = line.split()
    p1, p2 = int(f[1]), int(f[2])
    d = abs(p2 - p1)
    dists.append(d)
print("n", len(dists))
print("min", min(dists), "max", max(dists), "mean", statistics.mean(dists), "median", statistics.median(dists))
for thresh in [100, 500, 1000, 5000, 10000, 50000, 100000, 1000000, 10000000]:
    print(f">={thresh}", sum(1 for d in dists if d >= thresh))

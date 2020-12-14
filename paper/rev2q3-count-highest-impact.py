from cyvcf2 import VCF
import collections
import sys

order = [line.strip() for line in open(sys.argv[1], 'rt') if line[0] != '#']

counts = collections.defaultdict(int)

for v in VCF(sys.argv[2]):
    hi = v.INFO.get("highest_impact_order")
    counts[order[hi]] += 1

for k, v in counts.items():
    print(f"{k}\t{v}")

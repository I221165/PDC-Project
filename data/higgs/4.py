#!/usr/bin/env python3
import sys
from bisect import bisect_left

if len(sys.argv)!=3:
    print("Usage: remap_ids.py in.edgelist out_prefix")
    sys.exit(1)

infile = sys.argv[1]
out_prefix = sys.argv[2]

# 1) collect all unique IDs
ids = set()
with open(infile) as f:
    for line in f:
        u,v, *rest = line.split()
        ids.add(int(u))
        ids.add(int(v))
ids = sorted(ids)

# 2) write mapping file
with open(out_prefix + "_map.txt","w") as mf:
    for idx,orig in enumerate(ids):
        mf.write(f"{idx} {orig}\n")

# 3) remap edges
with open(infile) as f, open(out_prefix + ".edgelist","w") as out:
    for line in f:
        toks = line.split()
        u = int(toks[0]); v = int(toks[1])
        lu = bisect_left(ids, u)
        lv = bisect_left(ids, v)
        out.write(f"{lu+1} {lv+1} " + " ".join(toks[2:]) + "\n")

print(f"Written {out_prefix}.edgelist and {out_prefix}_map.txt")

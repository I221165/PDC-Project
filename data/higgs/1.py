#!/usr/bin/env python3
import sys
from collections import defaultdict

def merge_snap(retweet, reply, mention, output):
    counts = defaultdict(lambda: [0,0,0])
    for fn, idx in [(retweet,1),(reply,2),(mention,2)]:
        with open(fn) as f:
            for line in f:
                u,v = line.split()
                counts[(u,v)][idx] += 1
    with open(output,'w') as o:
        for (u,v),(lk,sh,cm) in counts.items():
            o.write(f"{u} {v} {lk} {sh} {cm}\n")

if __name__=="__main__":
    if len(sys.argv)!=5:
        print("Usage: merge_higgs.py retweet reply mention out.edgelist")
        sys.exit(1)
    merge_snap(*sys.argv[1:])


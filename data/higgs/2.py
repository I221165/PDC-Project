#!/usr/bin/env python3
import sys
from collections import defaultdict

def merge_snap(retweet_file, reply_file, mention_file, output_file):
    # maps (u,v) -> [likes, shares, comments]
    counts = defaultdict(lambda: [0,0,0])

    # retweets → "shares"
    for fn, idx in [(retweet_file,1), (reply_file,2), (mention_file,2)]:
        with open(fn) as f:
            for line in f:
                toks = line.strip().split()
                if len(toks) < 2: 
                    continue
                u, v = toks[0], toks[1]
                counts[(u,v)][idx] += 1

    with open(output_file, 'w') as out:
        for (u,v),(lk,sh,cm) in counts.items():
            out.write(f"{u} {v} {lk} {sh} {cm}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: merge_higgs.py retweet.edges reply.edges mention.edges out.edgelist")
        sys.exit(1)
    merge_snap(*sys.argv[1:])


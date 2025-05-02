#!/usr/bin/env python3
import csv
import random

# 1) Collect all unique user IDs from your graph.edgelist
uids = set()
with open('out.edgelist') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        u, v = parts[0], parts[1]
        uids.add(int(u))
        uids.add(int(v))
uids = sorted(uids)

# 2) Define number of topics
D = 50

# 3) Write interests.csv
with open('interests.csv', 'w', newline='') as out:
    writer = csv.writer(out)
    # header
    header = ['userID'] + [f'topic{i+1}' for i in range(D)]
    writer.writerow(header)
    # one random 0/1 vector per user
    for u in uids:
        row = [u] + [random.choice([0,1]) for _ in range(D)]
        writer.writerow(row)

print(f"Generated interests.csv with {len(uids)} users and {D} topics.")


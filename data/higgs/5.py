#!/usr/bin/env python3
import sys
import csv

if len(sys.argv) != 4:
    print("Usage: remap_interests.py map.txt old_interests.csv new_interests.csv")
    sys.exit(1)

mapfile, old_csv, new_csv = sys.argv[1:]

# 1) Load the mapping: local_id → original_id
orig2local = {}
with open(mapfile) as f:
    for line in f:
        local, orig = line.split()
        orig2local[int(orig)] = int(local)

# 2) Read old interests and write new ones
with open(old_csv, newline='') as inf, open(new_csv, 'w', newline='') as outf:
    reader = csv.reader(inf)
    writer = csv.writer(outf)
    header = next(reader)
    writer.writerow(header)  # same header

    for row in reader:
        orig_id = int(row[0])
        if orig_id not in orig2local:
            # skip users not in the remapped graph
            continue
        local_id = orig2local[orig_id]
        writer.writerow([local_id] + row[1:])
print(f"Wrote remapped interests to {new_csv}")

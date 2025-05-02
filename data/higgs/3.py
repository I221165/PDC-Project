python3 - <<'PYCODE'
import csv, random
from collections import defaultdict

# first find all user IDs in your merged graph
uids = set()
with open('data/graph.edgelist') as f:
    for u,v,_,_,_ in (line.split() for line in f):
        uids.add(int(u)); uids.add(int(v))
uids = sorted(uids)

D = 50
with open('data/interests.csv','w',newline='') as out:
    writer = csv.writer(out)
    header = ['userID'] + [f'topic{i}' for i in range(1,D+1)]
    writer.writerow(header)
    for u in uids:
        row = [u] + [random.choice([0,1]) for _ in range(D)]
        writer.writerow(row)
PYCODE


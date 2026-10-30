#!/usr/bin/env python3
"""Collect one Polus suite; standard library only."""
import csv
import re
import sys
from pathlib import Path
root = Path(sys.argv[1])
records = {}
for path in root.rglob('runs.csv'):
    match = re.fullmatch(r'(seq|omp)_(\d+)_(\d+)_(\d+)_(.+)', path.parent.name)
    if not match:
        continue
    mode, m, n, p, job = match.groups()
    key = (mode, int(m), int(n), int(p))
    records.setdefault(key, []).extend((row, str(path)) for row in csv.DictReader(path.open()))
best = {}
for key, rows in records.items():
    for row, path in rows:
        if (int(row['M']), int(row['N']), int(row['threads'])) != key[1:]:
            raise ValueError('Grid or actual thread count mismatch: ' + path)
        if float(row['time']) <= 0:
            raise ValueError('Timer resolution too low: ' + path)
    iterations = {int(row['iterations']) for row, _ in rows}
    if len(iterations) != 1:
        raise ValueError('Iteration count varies: ' + str(key))
    best[key] = min(rows, key=lambda item: float(item[0]['time']))
fields = ['mode','M','N','threads','runs','iterations','min_time','max_time','speedup_seq','source']
with (root/'summary.csv').open('w',newline='') as stream:
    writer = csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
    for key, (row, source) in sorted(best.items()):
        mode,m,n,p=key; baseline=best.get(('seq',m,n,1))
        speedup = float(baseline[0]['time'])/float(row['time']) if baseline else ''
        writer.writerow(dict(mode=mode,M=m,N=n,threads=p,runs=len(records[key]),
            iterations=row['iterations'],min_time=row['time'],max_time=max(float(r['time']) for r,_ in records[key]),
            speedup_seq=speedup,source=source))
with (root/'table1.csv').open('w',newline='') as stream:
    writer=csv.writer(stream);writer.writerow(['threads','M','N','iterations','time','speedup_seq'])
    for m,n,ps in [(400,600,[2,4,8,16]),(800,1200,[4,8,16,32])]:
        for p in ps:
            value=best.get(('omp',m,n,p));base=best.get(('seq',m,n,1))
            if value and base:
                row=value[0];writer.writerow([p,m,n,row['iterations'],row['time'],float(base[0]['time'])/float(row['time'])])
            else:writer.writerow([p,m,n,'','',''])
print(root/'summary.csv');print(root/'table1.csv')

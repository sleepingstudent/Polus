#!/usr/bin/env python3
"""Stream saved CSVs and compare OpenMP solutions against sequential solutions."""
import csv, itertools, math, re, sys
from pathlib import Path
root=Path(sys.argv[1]);files={}
for path in root.rglob('solution.csv'):
    m=re.fullmatch(r'(seq|omp)_(\d+)_(\d+)_(\d+)_(.+)',path.parent.name)
    if m:
        mode,M,N,p,job=m.groups();files[(mode,int(M),int(N),int(p))]=path
for (mode,M,N,p),path in sorted(files.items()):
    if mode!='omp':continue
    base=files.get(('seq',M,N,1))
    if base is None:
        print('No sequential baseline for',M,N,p);continue
    maximum=0.;squares=0.;count=0
    with base.open() as a,path.open() as b:
        for x,y in itertools.zip_longest(csv.DictReader(a),csv.DictReader(b)):
            if x is None or y is None or (x['i'],x['j'],x['x'],x['y'])!=(y['i'],y['j'],y['x'],y['y']):
                raise ValueError('CSV grid mismatch')
            d=abs(float(x['value'])-float(y['value']))
            if not math.isfinite(d):raise ValueError('Nonfinite solution')
            maximum=max(maximum,d);squares+=d*d;count+=1
    print(f'{M}x{N}, p={p}: points={count}, max_abs={maximum:.12e}, E_diff={math.sqrt(25/(M*N)*squares):.12e}')
    if maximum>1e-6:raise SystemExit('Solution difference exceeds 1e-6; inspect logs')

#!/usr/bin/env python3
"""Isolated before/after HTSlib feature comparison; full builds, exact output."""
import argparse,gzip,json,os,statistics,subprocess,sys
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--integration',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--parallel',type=Path,required=True);a=p.parse_args();assert not a.work.exists();a.work.mkdir(parents=True)
reports=[];arrays=texts=0
for repeat in range(3):
 variants=['baseline','candidate','parallel'];roots={name:a.work/f'run-{repeat+1}'/name for name in variants}
 for name in variants[repeat:]+variants[:repeat]:
  root=roots[name];root.mkdir(parents=True)
  command=['bash',str(a.integration/'scripts/run_constrained.sh'),sys.executable,str(a.integration/'scripts/benchmark_lifecycle.py'),'--mode','rust','--dataset','airway','--through','build','--cache','fresh','--work',str(root),'--binary',str(getattr(a,name))]
  print('RUN',repeat+1,name,flush=True)
  with (root/'container.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
  report=json.loads((root/'report-fresh.json').read_text());report.update(variant=name,repeat=repeat+1);reports.append(report)
 for variant in ['candidate','parallel']:
  baseline=roots['baseline']/'results';candidate=roots[variant]/'results'
  want={p.relative_to(baseline) for p in baseline.rglob('*') if p.is_file()};got={p.relative_to(candidate) for p in candidate.rglob('*') if p.is_file()};assert want==got
  for path in want:
   if path.suffix=='.hdf5':
    with h5py.File(baseline/path) as x,h5py.File(candidate/path) as y:
     keys=[];x.visit(keys.append);other=[];y.visit(other.append);assert keys==other
     for key in keys:
      if isinstance(x[key],h5py.Dataset):
       assert x[key].dtype==y[key].dtype and x[key].shape==y[key].shape
       np.testing.assert_array_equal(x[key][...],y[key][...]);arrays+=1
   elif path.suffix in ['.gz','.txt','.gff3','.bed']:
    read=lambda p:gzip.open(p,'rb').read() if p.suffix=='.gz' else p.read_bytes()
    assert read(baseline/path)==read(candidate/path);texts+=1
metrics=['wall_seconds','cpu_seconds','mean_cores','max_process_rss_kib','cgroup_peak_bytes','input_blocks','output_blocks']
summary={name:{key:statistics.median(r[key] for r in reports if r['variant']==name) for key in metrics} for name in roots}
report=dict(summary=summary,runs=reports,exact_arrays=arrays,exact_texts=texts,result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)

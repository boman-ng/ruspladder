#!/usr/bin/env python3
"""Run interleaved, isolated 8-CPU/16-GiB lifecycle measurements and parity checks."""
import argparse,json,os,statistics,subprocess,sys
from pathlib import Path
import h5py,numpy,scipy,statsmodels
from compare_lifecycle import compare

def main():
 p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--exporter',type=Path,required=True);p.add_argument('--repeats',type=int,default=3);a=p.parse_args()
 assert a.repeats>=3
 assert not a.work.exists(),'use a fresh benchmark directory'
 a.work.mkdir(parents=True)
 repo=Path(__file__).resolve().parent.parent;reports=[];parity=[]
 task_root=Path(os.environ.get('RUSPLADDER_WORK_ROOT','/home/wubw/data/ruspladder'))
 revision=lambda path:subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip()
 metadata=dict(reference_commit=revision(task_root/'upstream/spladder'),native_commit=revision(repo),python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,h5py=h5py.__version__,filesystem_cache='not flushed; fresh/reused refer to application outputs',measurement='child commands, startup included; staging and container startup excluded; cgroup peak includes worker and staging')
 assert metadata['reference_commit']=='65ceec839b9ff0cf96703c1605ee43667662f410'
 for dataset,through,sparse in [('airway','build',False),('events','test',False),('events','test',True)]:
  label=dataset+('-sparse' if sparse else '-direct')
  for repeat in range(a.repeats):
   order=['reference','rust'] if repeat%2==0 else ['rust','reference']
   roots={mode:a.work/label/f'run-{repeat+1}'/mode for mode in order}
   for cache in ['fresh','reused']:
    for mode in order:
     root=roots[mode];root.mkdir(parents=True,exist_ok=True)
     command=['bash',str(repo/'scripts/run_constrained.sh'),sys.executable,str(repo/'scripts/benchmark_lifecycle.py'),'--mode',mode,'--dataset',dataset,'--cache',cache,'--through',through,'--work',str(root),'--binary',str(a.binary)]
     if sparse:command+=['--sparse']
     print(f'RUN {label} {repeat+1}/{a.repeats} {cache} {mode}',flush=True)
     with (root/f'container-{cache}.log').open('w') as log:subprocess.run(command,cwd=repo,stdout=log,stderr=subprocess.STDOUT,check=True)
     report=json.loads((root/f'report-{cache}.json').read_text());report.update(repeat=repeat+1,label=label);reports.append(report)
    matched=compare(roots['reference']/'results',roots['rust']/'results',a.exporter);matched.update(label=label,repeat=repeat+1,cache=cache);parity.append(matched)
    (a.work/'progress.json').write_text(json.dumps(dict(runs=reports,parity=parity),indent=2)+'\n')
 summaries=[]
 metrics=['wall_seconds','cpu_seconds','mean_cores','utilization_8cpu_percent','max_process_rss_kib','cgroup_peak_bytes','input_blocks','output_blocks']
 for label in sorted({r['label'] for r in reports}):
  for cache in ['fresh','reused']:
   values={mode:{metric:statistics.median(r[metric] for r in reports if r['label']==label and r['cache']==cache and r['mode']==mode) for metric in metrics} for mode in ['reference','rust']}
   summaries.append(dict(label=label,cache=cache,median=values,wall_speedup=values['reference']['wall_seconds']/values['rust']['wall_seconds'],peak_memory_reduction=1-values['rust']['cgroup_peak_bytes']/values['reference']['cgroup_peak_bytes']))
 (a.work/'report.json').write_text(json.dumps(dict(metadata=metadata,repeats=a.repeats,summary=summaries,runs=reports,parity=parity,result='pass'),indent=2)+'\n')
 print(json.dumps(summaries,indent=2),flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Measure fresh/reused build + differential test inside an enforced cgroup.

Input staging precedes timing. Each invocation needs its own container so the
cgroup peak for reuse cannot include the earlier fresh build.
"""
import argparse,hashlib,json,os,resource,signal,subprocess,sys,time
from pathlib import Path
from resource_probe import snapshot

def main():
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['reference','rust'],required=True);p.add_argument('--dataset',choices=['airway','events'],required=True);p.add_argument('--cache',choices=['fresh','reused'],required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--sparse',action='store_true');p.add_argument('--through',choices=['build','test'],default='test')
 a=p.parse_args();task_root=Path(os.environ.get('RUSPLADDER_WORK_ROOT',str(Path.home() / 'data/ruspladder')));a.work.mkdir(parents=True,exist_ok=True)
 before=snapshot();assert len(before['cpu_affinity'])==4 and before['effective_memory_limit']==4*1024**3
 if a.dataset=='airway':
  source=task_root/'fixtures/airway';bams=sorted(source.glob('*.bam'));annotation=source/'Homo_sapiens.GRCh37.75_subset.gtf';reference=None;extra=['--ignore-mismatches','-n','101'];assert len(bams)==8
 else:
  source=task_root/'upstream/spladder/tests/testcase_events/data';bams=[source/f'align/testcase_events_1_sample{i}.bam' for i in range(1,21)];annotation=source/'testcase_events.gtf';reference=source/'genome.fa';extra=[]
 inputs=a.work/'inputs';inputs.mkdir(exist_ok=True);local_annotation=inputs/annotation.name
 if a.cache=='fresh':
  assert not (a.work/'results').exists(),'fresh run needs a new output directory'
  local_annotation.write_bytes(annotation.read_bytes())
  for path in bams:
   (inputs/path.name).symlink_to(path)
   (inputs/(path.name+'.bai')).symlink_to(Path(str(path)+'.bai'))
 else:
  assert (a.work/'report-fresh.json').exists(),'reuse requires completed fresh run'
  assert json.loads((a.work/'report-fresh.json').read_text())['returncode']==0,'fresh workflow did not complete'
 paths=[inputs/b.name for b in bams];names=[b.stem for b in bams];split=len(names)//2
 binary=Path(sys.executable).parent/'spladder' if a.mode=='reference' else a.binary
 build=[str(binary),'build','-a',str(local_annotation),'-b',','.join(map(str,paths)),'-o',str(a.work/'results'),'--parallel','4',*extra]
 if reference:build+=['--reference',str(reference)]
 if a.sparse:build+=['--sparse-bam']
 test=[str(binary),'test','-o',str(a.work/'results'),'-a',','.join(names[:split]),'-b',','.join(names[split:]),'--parallel','4','--event-types','exon_skip','--dpsi','0']
 usage_before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic();stages=[];returncode=0;soft_exceeded=False;timed_out=False
 for name,command in ([('build',build)] if a.through=='build' else [('build',build),('test',test)]):
  stage_start=time.monotonic();cpu_before=resource.getrusage(resource.RUSAGE_CHILDREN)
  with (a.work/f'{a.cache}-{name}.log').open('w') as log:
   process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   while True:
    elapsed=time.monotonic()-start
    if elapsed>=1800 and not soft_exceeded:
     soft_exceeded=True
     (a.work/f'soft-timeout-{a.cache}.json').write_text(json.dumps(dict(wall_seconds=elapsed,hard_timeout_seconds=3600,action='record and continue'))+'\n')
    if elapsed>=3600:
     timed_out=True;os.killpg(process.pid,signal.SIGKILL);process.wait();break
    try:
     process.wait(timeout=min(10,3600-elapsed));break
    except subprocess.TimeoutExpired:pass
  cpu_after=resource.getrusage(resource.RUSAGE_CHILDREN)
  stages.append(dict(stage=name,wall_seconds=time.monotonic()-stage_start,user_seconds=cpu_after.ru_utime-cpu_before.ru_utime,system_seconds=cpu_after.ru_stime-cpu_before.ru_stime,returncode=process.returncode,command=command))
  if process.returncode:returncode=process.returncode;break
 elapsed=time.monotonic()-start;usage=resource.getrusage(resource.RUSAGE_CHILDREN);after=snapshot();cpu=usage.ru_utime-usage_before.ru_utime+usage.ru_stime-usage_before.ru_stime
 with binary.open('rb') as f:binary_hash=hashlib.file_digest(f,'sha256').hexdigest()
 report=dict(mode=a.mode,dataset=a.dataset,cache=a.cache,sparse=a.sparse,through=a.through,binary_sha256=binary_hash,wall_seconds=elapsed,cpu_seconds=cpu,mean_cores=cpu/elapsed,utilization_4cpu_percent=25*cpu/elapsed,max_process_rss_kib=usage.ru_maxrss,cgroup_peak_bytes=after['memory_limits'][0]['peak'],input_blocks=usage.ru_inblock-usage_before.ru_inblock,output_blocks=usage.ru_oublock-usage_before.ru_oublock,stages=stages,resources_before=before,resources_after=after,returncode=returncode,soft_timeout_seconds=1800,hard_timeout_seconds=3600,soft_exceeded=soft_exceeded or elapsed>=1800,timed_out=timed_out)
 (a.work/f'report-{a.cache}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
 if returncode:raise SystemExit(returncode)
if __name__=='__main__':main()

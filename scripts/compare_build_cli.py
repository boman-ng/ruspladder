#!/usr/bin/env python3
"""Compare unmodified build CLI graphs, events, public files, and reused caches."""
import argparse,gzip,json,os,pickle,shutil,subprocess,sys
import h5py
import numpy as np
from pathlib import Path
from compare_annotations import serialize
from compare_collection import event_value
from compare_count_io import compare_hdf5
from compare_merge import load_genes
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('exporter',type=Path);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--airway',type=Path);p.add_argument('--limit',type=int)
a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
basic=a.upstream/'tests/testcase_basic/data';events=a.upstream/'tests/testcase_events/data'
scenarios=[]
for strand in ['pos','neg']:
 for mode in ['merge_graphs','merge_bams','merge_all','single']:
  scenarios.append((f'{strand}-{mode}',basic/f'annotation_{strand}.gtf',[basic/f'align/{strand}_{i}.bam' for i in ([1] if mode=='single' else [1,2])],basic/f'genome_{strand}.fa',['-M',mode]))
for ext in ['bam','cram']:
 scenarios.append((ext+'-outputs',events/'testcase_events.gtf',[events/f'align/testcase_events_1_sample{i}.{ext}' for i in [1,2]],events/'genome.fa',['--output-txt','--output-struc','--output-gff3','--output-bed','--output-struc-conf','--output-conf-bed','--output-conf-tcga','--output-conf-icgc','--labels','ignoredA,ignoredB']))
scenarios += [('validated',events/'testcase_events.gtf',[events/f'align/testcase_events_1_sample{i}.bam' for i in [1,2]],events/'genome.fa',['--validate-sg','--validate-sg-count','1','--use-anno-support','--no-compress-text']),('annotation-only',events/'testcase_events.gtf',[events/'align/testcase_events_1_sample1.bam'],events/'genome.fa',['--no-insert-es','--no-insert-ni','--no-insert-ir','--no-quantify-graph','--no-curate-alt-prime'])]
if a.airway:scenarios.append(('airway',a.airway/'Homo_sapiens.GRCh37.75_subset.gtf',[a.airway/f'SRR103950{i}_subset.bam' for i in [8,9]],None,['--ignore-mismatches','-n','101']))
if a.limit:scenarios=scenarios[:a.limit]
comparisons=dict(cases=0,graphs=0,events=0,hdf5_arrays=0,text_files=0,reused=0)
def run(binary,out,annotation,bams,ref,extra,parallel,label):
 flags=['build','-o',str(out),'-a',str(annotation),'-b',','.join(map(str,bams)),*extra,'--parallel',str(parallel)]
 if ref:flags+=['--reference',str(ref)]
 with open(out.parent/(label+'.log'),'w') as log:subprocess.run([str(binary),*flags],stdout=log,stderr=subprocess.STDOUT,check=True)
def compare(reference,native):
 requests=[];expected=[]
 for path in sorted(reference.rglob('*.pickle')):
  if path.name.endswith('.confirmed.pickle') or path.name.endswith('.quick_ids_segs') or path.name.endswith('.quick_ids_edges'):continue
  if path.parent.name=='spladder':
   requests.append(dict(graph=str(native/path.relative_to(reference).with_suffix('.hdf5'))));expected.append([serialize(g) for g in load_genes(path)]);comparisons['graphs']+=1
  else:
   with open(path,'rb') as f:values=pickle.load(f)
   requests.append(dict(events=str(native/path.relative_to(reference).with_suffix('.events.hdf5'))));expected.append([event_value(e) for e in values]);comparisons['events']+=1
 result=subprocess.run([str(a.exporter)],input=''.join(json.dumps(r)+'\n' for r in requests),stdout=subprocess.PIPE,text=True,check=True)
 actual=[json.loads(line) for line in result.stdout.splitlines()];assert len(actual)==len(expected)
 for i,(want,got) in enumerate(zip(expected,actual)):
  if want!=got:
   (native.parent/'mismatch.json').write_text(json.dumps(dict(request=requests[i],expected=want,actual=got),indent=2));raise AssertionError(f'cache mismatch: {requests[i]}')
 public=[]
 for path in sorted(reference.rglob('*')):
  if not path.is_file():continue
  rel=path.relative_to(reference);target=native/rel
  if path.suffix=='.hdf5':comparisons['hdf5_arrays']+=compare_hdf5(path,target);public.append(rel)
  elif path.suffix in ['.txt','.gz','.bed','.gff3']:
   read=lambda x:gzip.open(x,'rb').read() if x.suffix=='.gz' else x.read_bytes()
   assert read(path)==read(target),(path,target);comparisons['text_files']+=1;public.append(rel)
 expected_public=set(public)
 actual_public={p.relative_to(native) for p in native.rglob('*') if p.is_file() and (p.suffix in ['.txt','.gz','.bed','.gff3'] or p.name.endswith(('.counts.hdf5','.count.hdf5','.gene_exp.hdf5')))}
 assert expected_public==actual_public,(expected_public-actual_public,actual_public-expected_public)
for name,annotation,bams,ref,extra in scenarios:
 root=a.work/name;root.mkdir(exist_ok=True);local=root/annotation.name;shutil.copyfile(annotation,local)
 reference=root/'reference';reference.mkdir(exist_ok=True)
 run(Path(sys.executable).parent/'spladder',reference,local,bams,ref,extra,1,'reference')
 for threads in [1,4,8]:
  native=root/f'native-{threads}';native.mkdir(exist_ok=True)
  run(a.binary,native,local,bams,ref,extra,threads,f'native-{threads}')
  compare(reference,native)
  before={p.relative_to(native):p.stat().st_mtime_ns for p in native.rglob('*') if p.is_file()}
  run(a.binary,native,local,bams,ref,extra,threads,f'native-{threads}-reuse')
  assert before=={p.relative_to(native):p.stat().st_mtime_ns for p in native.rglob('*') if p.is_file()},'reused cache was rewritten'
  compare(reference,native);comparisons['cases']+=1;comparisons['reused']+=1
  if threads>1:
   for path in (root/'native-1').rglob('*.hdf5'):
    if not path.name.endswith(('.count.hdf5','.counts.hdf5','.gene_exp.hdf5')):continue
    with h5py.File(path) as one,h5py.File(native/path.relative_to(root/'native-1')) as parallel:
     for key in one:
      x,y=one[key][...],parallel[key][...]
      np.testing.assert_array_equal(x,y,err_msg=f'thread determinism {path} {key}')
      if x.dtype.kind=='f':assert x.tobytes()==y.tobytes(),f'float bits differ: {path} {key}'
 print('PASS:',name,flush=True)
report=dict(**comparisons,threads=[1,4,8],result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

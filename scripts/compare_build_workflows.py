#!/usr/bin/env python3
"""Exercise externally chunked merging, per-sample quantification and collection."""
import argparse,json,os,pickle,shutil,subprocess,sys
from pathlib import Path
from compare_annotations import serialize
from compare_collection import event_value
from compare_count_io import compare_hdf5
from compare_merge import load_genes
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('exporter',type=Path);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
source=a.upstream/'tests/testcase_basic/data';bams=[source/f'align/pos_{i}.bam' for i in [2,1]];genome=source/'genome_pos.fa'
annotation=a.work/'annotation.gtf';shutil.copyfile(source/'annotation_pos.gtf',annotation)
arrays=0;graphs=0;failures=[]
def run(binary,dest,selected,extra,label,expect=None):
 dest.mkdir(parents=True,exist_ok=True)
 command=[str(binary),'build','-a',str(annotation),'-b',','.join(map(str,selected)),'-o',str(dest),'--reference',str(genome),*extra]
 result=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (dest.parent/(label+'.log')).write_text(result.stdout)
 if expect is None:assert result.returncode==0,(command,result.stdout[-2000:])
 else:assert result.returncode!=0 and expect in result.stdout,(command,result.stdout[-2000:])
def compare(reference,native):
 global arrays,graphs
 for path in reference.rglob('*.hdf5'):arrays+=compare_hdf5(path,native/path.relative_to(reference))
 req=[];want=[]
 for path in (reference/'spladder').glob('*.pickle'):
  req.append(dict(graph=str(native/'spladder'/path.with_suffix('.hdf5').name)));want.append([serialize(g) for g in load_genes(path)]);graphs+=1
 result=subprocess.run([str(a.exporter)],input=''.join(json.dumps(x)+'\n' for x in req),text=True,stdout=subprocess.PIPE,check=True)
 assert [json.loads(x) for x in result.stdout.splitlines()]==want
for workflow in ['collect','chunks']:
 root=a.work/workflow;reference=root/'reference';native=root/'native'
 for binary,dest in [(Path(sys.executable).parent/'spladder',reference),(a.binary,native)]:
  if workflow=='collect':
   run(binary,dest,bams,['--no-quantify-graph','--no-extract-ase'],dest.name+'-graph')
   for i,bam in enumerate(bams):run(binary,dest,[bam],['--qmode','single','--no-extract-ase'],dest.name+f'-single-{i}')
   run(binary,dest,bams,['--qmode','collect'],dest.name+'-collect')
  else:
   for start in [0,1]:run(binary,dest,bams,['--chunked-merge','1','2',str(start),str(start+1),'--no-quantify-graph','--no-extract-ase'],dest.name+f'-chunk-{start}')
   run(binary,dest,bams,['--chunked-merge','2','2','0','2'],dest.name+'-final')
 compare(reference,native)
for name,extra,source_error,native_error in [
 ('reinfer',['--re-infer-sg'],"name 'infer_splice_graph' is not defined",'undefined infer_splice_graph'),
 ('single-multiple',['-M','single'],'No such file or directory','No such file'),
 ('qmode-single-multiple',['--qmode','single'],'IndexError','invalid expression sample index'),
 ('nonfinal-default',['--chunked-merge','1','2','0','1'],'No such file or directory','No such file'),
 ('validate-merge-bams',['-M','merge_bams','--validate-sg'],'No such file or directory','No such file'),
]:
 root=a.work/name
 for binary,dest,error in [(Path(sys.executable).parent/'spladder',root/'reference',source_error),(a.binary,root/'native',native_error)]:run(binary,dest,bams,extra,dest.name,expect=error)
 failures.append(name)
report=dict(workflows=2,graphs=graphs,hdf5_arrays=arrays,upstream_failures=failures,result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

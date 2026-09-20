#!/usr/bin/env python3
"""Compare annotation and alignment prep CLI caches against unmodified SplAdder."""
import argparse,json,pickle,shutil,subprocess,sys
from pathlib import Path
from compare_annotations import serialize
from compare_count_io import compare_hdf5
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('exporter',type=Path);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
source=a.upstream/'tests/testcase_basic/data';cases=arrays=0
for confidence in range(4):
 for parallel in [1,4]:
  root=a.work/f'C{confidence}-threads{parallel}';outputs=[]
  for binary,name in [(Path(sys.executable).parent/'spladder','reference'),(a.binary,'native')]:
   dest=root/name;dest.mkdir(parents=True,exist_ok=True);(dest/'tmp').mkdir(exist_ok=True);gtf=dest/'annotation.gtf';shutil.copyfile(source/'annotation_pos.gtf',gtf);bams=[]
   for sample in [1,2]:
    original=source/f'align/pos_{sample}.bam';bam=dest/original.name;shutil.copyfile(original,bam);shutil.copyfile(Path(str(original)+'.bai'),Path(str(bam)+'.bai'));bams.append(bam)
   command=[str(binary),'prep','-a',str(gtf),'-b',','.join(map(str,bams)),'--sparse-bam','-c',str(confidence),'--parallel',str(parallel),'--tmp-dir',str(dest/'tmp')]
   with (dest/'command.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
   before={p.name:p.stat().st_mtime_ns for p in dest.glob('*.hdf5')}
   with (dest/'reuse.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
   assert before=={p.name:p.stat().st_mtime_ns for p in dest.glob('*.hdf5')}
   outputs.append(dest)
  reference,native=outputs
  for path in reference.glob('*.hdf5'):arrays+=compare_hdf5(path,native/path.name)
  with (reference/'annotation.gtf.pickle').open('rb') as f:genes=pickle.load(f)
  for g in genes:
   g.from_sparse()
   if not hasattr(g,'introns_anno'):g.populate_annotated_introns()
  result=subprocess.run([str(a.exporter)],input=json.dumps(dict(graph=str(native/'annotation.gtf.ruspladder.hdf5')))+'\n',text=True,stdout=subprocess.PIPE,check=True)
  assert json.loads(result.stdout)==[serialize(g) for g in genes]
  cases+=1
report=dict(cases=cases,hdf5_arrays=arrays,threads=[1,4],result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

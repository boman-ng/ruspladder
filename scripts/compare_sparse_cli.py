#!/usr/bin/env python3
"""Compare complete sparse-input builds and reuse, including public COO files."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
import h5py
import numpy as np
from compare_lifecycle import compare
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('exporter',type=Path);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--limit',type=int);a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
basic=a.upstream/'tests/testcase_basic/data';events=a.upstream/'tests/testcase_events/data'
cases=[(f'{s}-merge',basic/f'annotation_{s}.gtf',[basic/f'align/{s}_{i}.bam' for i in [1,2]],[]) for s in ['pos','neg']]
cases += [('events-C3',events/'testcase_events.gtf',[events/f'align/testcase_events_1_sample{i}.bam' for i in [1,2]],[]),('events-C0',events/'testcase_events.gtf',[events/f'align/testcase_events_1_sample{i}.bam' for i in [1,2]],['-c','0'])]
cases += [(mode,basic/'annotation_pos.gtf',[basic/f'align/pos_{i}.bam' for i in ([1] if mode=='single' else [1,2])],['-M',mode]) for mode in ['single','merge_bams','merge_all']]
if a.limit:cases=cases[:a.limit]
reports=[];arrays=0
for name,annotation,bams,extra in cases:
 root=a.work/name;outputs=[]
 for binary,tag,parallel in [(Path(sys.executable).parent/'spladder','reference',1),(a.binary,'native-1',1),(a.binary,'native-4',4)]:
  dest=root/tag;inputs=dest/'inputs';inputs.mkdir(parents=True,exist_ok=True);gtf=inputs/annotation.name;shutil.copyfile(annotation,gtf);local=[]
  for path in bams:
   bam=inputs/path.name;shutil.copyfile(path,bam);shutil.copyfile(Path(str(path)+'.bai'),Path(str(bam)+'.bai'));local.append(bam)
  command=[str(binary),'build','-a',str(gtf),'-b',','.join(map(str,local)),'-o',str(dest/'results'),'--sparse-bam','--parallel',str(parallel),*extra]
  with (dest/'command.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
  before={p.relative_to(dest):p.stat().st_mtime_ns for p in dest.rglob('*.hdf5')}
  with (dest/'reuse.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
  assert before=={p.relative_to(dest):p.stat().st_mtime_ns for p in dest.rglob('*.hdf5')}
  outputs.append(dest)
 reference=outputs[0]
 for native in outputs[1:]:
  reports.append(compare(reference/'results',native/'results',a.exporter))
  for path in (reference/'inputs').glob('*.hdf5'):
   with h5py.File(path) as want,h5py.File(native/'inputs'/path.name) as got:
    assert sorted(want)==sorted(got)
    for key in want:
     assert want[key].dtype==got[key].dtype and want[key].shape==got[key].shape
     np.testing.assert_array_equal(want[key][...],got[key][...]);arrays+=1
     expected_compression='gzip' if key.endswith('_reads_shp') and native.name=='native-4' else want[key].compression
     assert got[key].compression==expected_compression
 # Native 1/4-thread outputs must also have exact floating-point values.
 for path in outputs[1].rglob('*.hdf5'):
  if not path.name.endswith(('.count.hdf5','.counts.hdf5','.gene_exp.hdf5')):continue
  with h5py.File(path) as one,h5py.File(outputs[2]/path.relative_to(outputs[1])) as four:
   for key in one:
    x,y=one[key][...],four[key][...]
    np.testing.assert_array_equal(x,y)
    if x.dtype.kind=='f':assert x.tobytes()==y.tobytes(),f'float bits differ: {path} {key}'
 print('PASS:',name,flush=True)
report=dict(cases=len(reports),summary_arrays=arrays,outputs=reports,result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

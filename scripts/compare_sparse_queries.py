#!/usr/bin/env python3
"""Compare bounded COO queries with upstream add_reads_from_sparse_bam."""
import argparse,json,subprocess
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
from spladder.reads import add_reads_from_sparse_bam
p=argparse.ArgumentParser();p.add_argument('probe',type=Path);p.add_argument('--summaries',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
requests=[];expected=[]
for path in sorted(a.summaries.glob('case-*.hdf5'))[::4]:
 with h5py.File(path) as f:chromosomes=[k[:-10] for k in f if k.endswith('_reads_shp')]
 for chromosome in chromosomes:
  for start,stop in [(0,100),(99,1100),(500,500),(2000,15000)]:
   for strand in ['+','-']:
    for unstranded in [False,True]:
     for minimum in [None,3]:
      filters=None if minimum is None else dict(intron=350000,exon_len=100,mismatch=0,mincount=minimum)
      gene=SimpleNamespace(chr=chromosome,start=start,stop=stop,strand=strand);cache={}
      tracks=add_reads_from_sparse_bam(gene,str(path),chromosome,3,types=['exon_track'],filter=filters,cache=cache,unstranded=unstranded)[0]
      introns=[]
      for intron_strand in ['+','-']:
       gene.strand=intron_strand
       introns.append(add_reads_from_sparse_bam(gene,str(path),chromosome,3,types=['intron_list'],filter=filters,cache=cache,unstranded=False)[0])
      requests.append(dict(path=str(path),chromosome=chromosome,start=start,stop=stop,options=dict(strand=strand,filter=filters),unstranded=unstranded))
      expected.append(dict(coverage=np.asarray(tracks.sum(axis=0)).ravel().astype(int).tolist(),introns_plus=introns[0].tolist(),introns_minus=introns[1].tolist()))
assert requests, 'no sparse summaries were exercised'
result=subprocess.run([str(a.probe)],input=''.join(json.dumps(r)+'\n' for r in requests),text=True,stdout=subprocess.PIPE,check=True)
actual=[json.loads(x) for x in result.stdout.splitlines()];assert len(actual)==len(expected)
for request,want,got in zip(requests,expected,actual):
 if want!=got:
  (a.work/'mismatch.json').write_text(json.dumps(dict(request=request,expected=want,actual=got),indent=2));raise AssertionError(request)
report=dict(cases=len(requests),result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

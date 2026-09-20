#!/usr/bin/env python3
"""Check graphs, six event sets, quantification and differential results end to end."""
import argparse,gzip,json,pickle,subprocess
from pathlib import Path
import numpy as np
from compare_annotations import serialize
from compare_collection import event_value
from compare_count_io import compare_hdf5
from compare_merge import load_genes

def compare(reference,native,exporter):
 requests=[];expected=[];arrays=texts=tables=0
 for path in sorted(reference.glob('spladder/*.pickle')):
  requests.append(dict(graph=str(native/'spladder'/path.with_suffix('.hdf5').name)));expected.append([serialize(g) for g in load_genes(path)])
 graphs=len(requests)
 for path in sorted(reference.glob('*.pickle')):
  if path.name.endswith('.confirmed.pickle'):continue
  with path.open('rb') as f:values=pickle.load(f)
  requests.append(dict(events=str(native/path.with_suffix('.events.hdf5').name)));expected.append([event_value(e) for e in values])
 result=subprocess.run([str(exporter)],input=''.join(json.dumps(r)+'\n' for r in requests),text=True,stdout=subprocess.PIPE,check=True)
 actual=[json.loads(x) for x in result.stdout.splitlines()];assert len(actual)==len(expected)
 for req,want,got in zip(requests,expected,actual):assert want==got,req
 for path in sorted(reference.rglob('*.hdf5')):arrays+=compare_hdf5(path,native/path.relative_to(reference))
 for path in sorted(reference.rglob('*')):
  if path.suffix in ['.gz','.txt','.gff3','.bed']:
   read=lambda p:gzip.open(p,'rb').read() if p.suffix=='.gz' else p.read_bytes()
   assert read(path)==read(native/path.relative_to(reference)),path;texts+=1
  elif path.suffix=='.tsv':
   want=np.loadtxt(path,delimiter='\t',dtype=str,ndmin=2);got=np.loadtxt(native/path.relative_to(reference),delimiter='\t',dtype=str,ndmin=2)
   np.testing.assert_array_equal(want[0],got[0]);np.testing.assert_array_equal(want[1:,:6],got[1:,:6])
   want=want[1:,6:].astype(float);got=got[1:,6:].astype(float)
   np.testing.assert_allclose(want,got,atol=1e-8,rtol=1e-6,equal_nan=True);np.testing.assert_array_equal(np.round(want,6),np.round(got,6))
   for alpha in [0.01,0.05,0.1]:np.testing.assert_array_equal(want[:,:2]<=alpha,got[:,:2]<=alpha)
   tables+=1
 return dict(graphs=graphs,event_sets=len(requests)-graphs,hdf5_arrays=arrays,text_files=texts,test_tables=tables,result='pass')

def main():
 p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('native',type=Path);p.add_argument('exporter',type=Path);p.add_argument('--report',type=Path,required=True);a=p.parse_args()
 report=compare(a.reference,a.native,a.exporter);a.report.write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)
if __name__=='__main__':main()

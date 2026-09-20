#!/usr/bin/env python3
"""Complete test CLI for all six event types, with heterogeneous counted fixtures."""
import argparse,copy,json,pickle,shutil,subprocess,sys
from pathlib import Path
import h5py
import numpy as np
from compare_annotations import serialize
from compare_collection import event_value,KINDS
from compare_merge import load_genes
from compare_verification import sources
parser=argparse.ArgumentParser();parser.add_argument('binary',type=Path);parser.add_argument('importer',type=Path);parser.add_argument('--upstream',type=Path,required=True);parser.add_argument('--work',type=Path,required=True)
a=parser.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
source=a.upstream/'tests/testcase_events/results_merged';reference=a.work/'reference';native=a.work/'native'
for root in [reference,native]:
    (root/'spladder').mkdir(parents=True,exist_ok=True)
    shutil.copyfile(source/'spladder/genes_graph_conf3.merge_graphs.gene_exp.hdf5',root/'spladder/genes_graph_conf3.merge_graphs.gene_exp.hdf5')
graph=source/'spladder/genes_graph_conf3.merge_graphs.pickle';shutil.copyfile(graph,reference/'spladder'/graph.name)
with h5py.File(source/'spladder/genes_graph_conf3.merge_graphs.gene_exp.hdf5') as f:labels=f['samples'][:];n=len(labels);gene_count=len(f['gene_ids'])
representatives={}
for gene,event,*_ in sources(a.upstream):representatives.setdefault(event.event_type,event)
rng=np.random.default_rng(1717);count=64;cache_events=[]
for kind in KINDS:
    with h5py.File(source/f'merge_graphs_{kind}_C3.counts.hdf5') as f:features=f['event_features'][:]
    means=np.exp(rng.uniform(2,5,count));disp=rng.uniform(.01,.8,count)
    values=np.stack([rng.negative_binomial(1/d,(1/d)/(1/d+mu),(n,len(features))) for mu,d in zip(means,disp)],axis=2).astype(float)
    path=reference/f'merge_graphs_{kind}_C3.counts.hdf5'
    with h5py.File(path,'w') as f:
        f['samples']=labels;f['event_features']=features;f['event_counts']=values;f['conf_idx']=np.arange(count);f['gene_idx']=rng.integers(0,gene_count,count);f['psi']=rng.uniform(0,1,(n,count))
    shutil.copyfile(path,native/path.name)
    events=[copy.deepcopy(representatives[kind]) for _ in range(count)]
    for i,e in enumerate(events):e.id=i+1
    with open(reference/f'merge_graphs_{kind}_C3.pickle','wb') as f:pickle.dump(np.array(events,dtype=object),f)
    cache_events.append(dict(path=str(native/f'merge_graphs_{kind}_C3.events.hdf5'),values=[event_value(e) for e in events]))
request=dict(genes=[serialize(g) for g in load_genes(graph)],graph=str(native/'spladder/genes_graph_conf3.merge_graphs.hdf5'),events=cache_events)
subprocess.run([str(a.importer)],input=json.dumps(request)+'\n',text=True,stdout=subprocess.PIPE,check=True)
flags=['--dpsi','0','-a',','.join(x.decode() for x in labels[:n//2]),'-b',','.join(x.decode() for x in labels[n//2:])]
with open(a.work/'reference.log','w') as log:subprocess.run([str(Path(sys.executable).parent/'spladder'),'test','-o',str(reference),*flags],stdout=log,stderr=subprocess.STDOUT,check=True)
files=sorted((reference/'testing').glob('*.tsv'));assert len(files)==18
for threads in [1,4]:
    with open(a.work/f'native-{threads}.log','w') as log:subprocess.run([str(a.binary),'test','-o',str(native),*flags,'--parallel',str(threads)],stdout=log,stderr=subprocess.STDOUT,check=True)
    for path in files:
        want=np.loadtxt(path,dtype=str,delimiter='\t',ndmin=2);got=np.loadtxt(native/'testing'/path.name,dtype=str,delimiter='\t',ndmin=2)
        np.testing.assert_array_equal(got[0],want[0]);np.testing.assert_array_equal(got[1:,:6],want[1:,:6])
        np.testing.assert_allclose(got[1:,6:].astype(float),want[1:,6:].astype(float),atol=1e-8,rtol=1e-6,equal_nan=True,err_msg=path.name)
        np.testing.assert_array_equal(np.round(got[1:,6:].astype(float),6),np.round(want[1:,6:].astype(float),6))
        if threads==1:shutil.copyfile(native/'testing'/path.name,a.work/('threads1-'+path.name))
        else:assert (a.work/('threads1-'+path.name)).read_bytes()==(native/'testing'/path.name).read_bytes()
report=dict(event_types=KINDS,events_per_type=count,threads=[1,4],text_files=len(files)*2,result='pass')
(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

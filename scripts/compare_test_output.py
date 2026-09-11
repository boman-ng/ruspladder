#!/usr/bin/env python3
"""Exercise the upstream test writer with fixed statistical results and compare files."""
import argparse,contextlib,copy,importlib,inspect,io,json,pickle,subprocess,warnings
from pathlib import Path
from unittest.mock import patch
import h5py
import numpy as np
from spladder.spladder import parse_options
from compare_verification import sources
from compare_collection import event_value,KINDS
module=importlib.import_module('spladder.spladder_test')
parser=argparse.ArgumentParser()
parser.add_argument('probe',type=Path);parser.add_argument('--upstream',type=Path,required=True);parser.add_argument('--work',type=Path,required=True)
args=parser.parse_args();args.work.mkdir(parents=True,exist_ok=True);(args.work/'report.json').unlink(missing_ok=True)
representatives={}
for gene,event,*_ in sources(args.upstream):representatives.setdefault(event.event_type,event)
assert all(k in representatives for k in KINDS)
rng=np.random.default_rng(2346);requests=[];roots=[]
for kind,event in representatives.items():
    for count in [1,17,65,261]:
        for correction in ['BH','BY']:
            root=args.work/str(len(requests));roots.append(root);(root/'spladder').mkdir(parents=True,exist_ok=True)
            (root/'spladder/genes_graph_conf3.merge_graphs.pickle').touch()
            with h5py.File(root/f'merge_graphs_{kind}_C3.counts.hdf5','w') as f:f['conf_idx']=np.arange(count)
            events=np.array([copy.deepcopy(event) for _ in range(count)])
            with open(root/f'merge_graphs_{kind}_C3.pickle','wb') as f:pickle.dump(events,f)
            samples=8;group_a=3;labels=np.array([f'sample_{s}' for s in range(samples)])
            genes=rng.integers(5,1000,(11,samples)).astype(float);gene_ids=np.array([f'g{i}' for i in range(11)]);gene_symbols=np.array([f'G{i}' for i in range(11)])
            coverage=rng.integers(2,100,(2,count,samples)).astype(float);coverage[0,0,:group_a]=0
            psi=rng.uniform(0,1,(count,samples));gene_idx=rng.integers(0,11,count);event_idx=rng.permutation(count)
            _,options=parse_options(['spladder','test','-o',str(root),'-a',','.join(labels[:group_a]),'-b',','.join(labels[group_a:]),'--event-types',kind,'-C',correction])
            options.min_dpsi=0.05
            captured={}
            def get_expression(*args,**kwargs):return genes.copy(),labels.copy(),np.arange(samples),gene_ids.copy(),gene_symbols.copy()
            def quantify(*args,**kwargs):return [coverage[0].copy(),coverage[1].copy()],psi.copy(),gene_idx.copy(),event_idx.copy(),None,labels.copy()
            def result(cov,null,alternative,sf,options,event_type,selected):
                frame=inspect.currentframe().f_back.f_locals
                captured['prepared']=dict(counts=cov.copy(),size_factors=sf.copy(),event_size_factors=frame['sf_ev'].copy(),null=null.copy(),alternative=alternative.copy(),selected=selected.copy(),delta_psi=frame['delta_psi'].copy(),event_idx=frame['event_idx'].copy(),gene_idx=frame['gene_idx'].copy())
                pvalues=np.resize(np.array([1.,0.,0.05,0.1,1e-20,0.5,0.5]),count)
                raw=np.linspace(0.001,1,count)[:,None];adjusted=np.linspace(.1,.5,count)[:,None];covused=cov[:count].copy()
                captured['result']=dict(pvalues=pvalues,coverage=covused,dispersion_raw=raw.ravel(),dispersion_adjusted=adjusted.ravel())
                return pvalues,covused,raw,adjusted
            with patch.object(module,'_get_gene_expression',get_expression),patch.object(module.quantify,'quantify_from_counted_events',quantify),patch.object(module,'run_testing',result),contextlib.redirect_stdout(io.StringIO()),warnings.catch_warnings():
                warnings.simplefilter('ignore');module.spladder_test(options)
            requests.append(dict(output=str(root/'native'),options=dict(confidence=3,kind=kind,correction=correction,group_a=group_a,labels=['condA','condB']),metadata=dict(samples=labels,gene_ids=gene_ids,gene_symbols=gene_symbols),events=[event_value(e) for e in events],**captured))
run=subprocess.run([str(args.probe)],input=''.join(json.dumps(r,default=lambda x:x.tolist())+'\n' for r in requests),text=True,stdout=subprocess.PIPE,check=True)
assert len(run.stdout.splitlines())==len(requests)
files=0
for i,root in enumerate(roots):
    for expected in (root/'testing').glob('*.tsv'):
        want=np.loadtxt(expected,delimiter='\t',dtype=str,ndmin=2);got=np.loadtxt(root/'native'/expected.name,delimiter='\t',dtype=str,ndmin=2)
        np.testing.assert_array_equal(got[0],want[0]);np.testing.assert_array_equal(got[1:,:6],want[1:,:6])
        np.testing.assert_allclose(got[1:,6:].astype(float),want[1:,6:].astype(float),atol=1e-8,rtol=1e-6,equal_nan=True,err_msg=f'{i} {expected.name}')
        np.testing.assert_array_equal(np.round(got[1:,6:8].astype(float),6),np.round(want[1:,6:8].astype(float),6))
        files+=1
    for expected in (root/'testing').glob('test_setup*.pickle'):
        with open(expected,'rb') as f:want=pickle.load(f)
        with h5py.File(root/'native'/expected.name.replace('.pickle','.hdf5')) as got:
            assert sorted(got)==sorted(want)
            for name,value in want.items():
                actual=got[name][()]
                if got[name].dtype.kind=='S':actual=actual.astype(str)
                np.testing.assert_array_equal(actual,value,err_msg=f'{i} {name}')
report=dict(cases=len(requests),text_files=files,setup_files=len(requests),result='pass')
(args.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

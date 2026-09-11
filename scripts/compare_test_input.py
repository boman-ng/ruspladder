#!/usr/bin/env python3
"""Capture the real upstream test CLI's input boundary before statistical fitting."""
import argparse
import contextlib
import importlib
import inspect
import io
import json
import pickle
import subprocess
import warnings
from pathlib import Path
from unittest.mock import patch
import h5py
import numpy as np
from spladder.spladder import parse_options

module = importlib.import_module('spladder.spladder_test')
parser = argparse.ArgumentParser()
parser.add_argument('probe', type=Path)
parser.add_argument('--work', type=Path, required=True)
args = parser.parse_args()
args.work.mkdir(parents=True,exist_ok=True)
(args.work/'report.json').unlink(missing_ok=True)
root=args.work/'reference'
(root/'spladder').mkdir(parents=True,exist_ok=True)
(root/'spladder/genes_graph_conf3.merge_graphs.pickle').touch()
with h5py.File(root/'merge_graphs_exon_skip_C3.counts.hdf5','w') as f: f['conf_idx']=np.arange(60)
with open(root/'merge_graphs_exon_skip_C3.pickle','wb') as f: pickle.dump(np.arange(60),f)
rng=np.random.default_rng(1368)
requests,expected=[],[]
class Captured(Exception): pass
for samples in [4,8,16,33]:
    group_a=samples//2
    labels=np.array([f'sample_{i}' for i in range(samples)])
    for profile in ['regular','outliers','filtered']:
        genes=rng.negative_binomial(3,0.04,(30,samples)).astype(float)
        genes[0,0]=1e8; genes[1,-1]=5e7
        coverage=rng.negative_binomial(2,0.1,(2,60,samples)).astype(float)
        coverage[0,0]=0;coverage[1,1]=0;coverage[:,2,:group_a]=0
        if profile=='outliers': coverage[:,3,0]=1e7;coverage[:,4,-1]=1e9
        elif profile=='filtered': coverage[:]=1
        psi=rng.uniform(0,1,(60,samples));psi[0,:group_a]=np.nan;psi[1,group_a:]=np.nan;psi[2,0]=np.nan
        gene_idx=rng.integers(0,30,60);event_idx=rng.permutation(60)
        for cap_exp in [False,True]:
            for cap_ev in [False,True]:
                for min_dpsi in [0.,0.05,0.6]:
                    _,options=parse_options(['spladder','test','-o',str(root),'-a',','.join(labels[:group_a]),'-b',','.join(labels[group_a:]),'--event-types','exon_skip'])
                    options.cap_exp_outliers=cap_exp;options.cap_outliers=cap_ev;options.min_dpsi=min_dpsi
                    captured={}
                    def get_expression(*args,**kwargs): return genes.copy(),labels.copy(),np.arange(samples),np.array([f'g{i}' for i in range(30)]),np.array([f'G{i}' for i in range(30)])
                    def quantify(*args,**kwargs):
                        frame=inspect.currentframe().f_back.f_locals
                        captured['expression']=dict(counts=frame['gene_counts'].copy(),size_factors=frame['sf_ge'].copy(),capped=frame.get('outlier_cnt',0))
                        return [coverage[0].copy(),coverage[1].copy()],psi.copy(),gene_idx.copy(),event_idx.copy(),None,labels.copy()
                    def capture(cov,null,alternative,sf,options,event_type,selected):
                        frame=inspect.currentframe().f_back.f_locals
                        captured['prepared']=dict(counts=cov.copy(),size_factors=sf.copy(),event_size_factors=frame['sf_ev'].copy(),null=null.copy(),alternative=alternative.copy(),selected=selected.copy(),delta_psi=frame['delta_psi'].copy(),event_idx=frame['event_idx'].copy(),gene_idx=frame['gene_idx'].copy())
                        a=np.nanmean(cov[:,:group_a]/sf[:group_a],axis=1);b=np.nanmean(cov[:,group_a:samples]/sf[group_a:samples],axis=1)
                        ga=np.nanmean(cov[:,samples:samples+group_a]/sf[samples:samples+group_a],axis=1);gb=np.nanmean(cov[:,samples+group_a:]/sf[samples+group_a:],axis=1)
                        captured['means']=np.c_[a,b,np.log2(a)-np.log2(b),ga,gb,np.log2(ga)-np.log2(gb)]
                        raise Captured()
                    with patch.object(module,'_get_gene_expression',get_expression),patch.object(module.quantify,'quantify_from_counted_events',quantify),patch.object(module,'run_testing',capture),contextlib.redirect_stdout(io.StringIO()),warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        try: module.spladder_test(options)
                        except Captured: pass
                    captured.setdefault('prepared',None);captured.setdefault('means',None)
                    expected.append(captured)
                    requests.append(dict(counts=genes.tolist(),samples=samples,cap_expression=cap_exp,group_a=group_a,options=dict(cap_outliers=cap_ev,max_zero_fraction=0.5,min_dpsi=min_dpsi),
                        events=dict(coverage=coverage.tolist(),psi=[[None if np.isnan(v) else v for v in row] for row in psi],gene_idx=gene_idx.tolist(),event_idx=event_idx.tolist(),samples=labels.tolist())))
run=subprocess.run([str(args.probe)],input=''.join(json.dumps(r)+'\n' for r in requests),text=True,stdout=subprocess.PIPE,check=True)
actual=[json.loads(s) for s in run.stdout.splitlines()]
assert len(actual)==len(expected)
def compare(want,got,key=''):
    if want is None: assert got is None,key
    elif isinstance(want,dict):
        for k in want:compare(want[k],got[k],key+'.'+k)
    elif key.endswith(('counts','capped','null','alternative','selected','event_idx','gene_idx')):
        if key=='.expression.counts':np.testing.assert_allclose(got,want,atol=1e-10,rtol=1e-8,err_msg=key)
        else:np.testing.assert_array_equal(got,want,err_msg=key)
    else:np.testing.assert_allclose(np.asarray(got,dtype=float),want,atol=1e-10,rtol=1e-8,equal_nan=True,err_msg=key)
for i,(want,got) in enumerate(zip(expected,actual)):
    try:compare(want,got)
    except AssertionError:
        (args.work/'mismatch.json').write_text(json.dumps(dict(case=i,request=requests[i],expected=want,actual=got),default=lambda x:x.tolist(),indent=2)+'\n');raise
report=dict(cases=len(actual),skipped=sum(w['prepared'] is None for w in expected),result='pass')
(args.work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS:',report)

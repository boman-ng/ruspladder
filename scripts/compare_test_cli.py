#!/usr/bin/env python3
"""Run unmodified SplAdder test and the native test CLI on counted BAM/CRAM fixtures."""
import argparse,json,os,pickle,shutil,subprocess,sys
from pathlib import Path
import h5py
import numpy as np
from compare_annotations import serialize
from compare_collection import event_value
from compare_merge import load_genes
parser=argparse.ArgumentParser();parser.add_argument('binary',type=Path);parser.add_argument('importer',type=Path);parser.add_argument('--upstream',type=Path,required=True);parser.add_argument('--work',type=Path,required=True)
args=parser.parse_args();args.work.mkdir(parents=True,exist_ok=True);(args.work/'report.json').unlink(missing_ok=True)
files=0;cases=0
for suffix,order_a,order_b in [('',list(range(1,11)),list(range(11,21))),('_cram',[10,8,2,1,7,6,5,3,9,4],[20,13,17,11,12,19,15,14,16,18])]:
    source=args.upstream/'tests/testcase_events'/('results_merged'+suffix)
    root=args.work/('cram' if suffix else 'bam');reference=root/'reference';native=root/'native'
    for dest in [reference,native]:(dest/'spladder').mkdir(parents=True,exist_ok=True)
    for f in (source/'spladder').glob('*gene_exp.hdf5'):
        for dest in [reference,native]:shutil.copyfile(f,dest/'spladder'/f.name)
    graph=source/'spladder/genes_graph_conf3.merge_graphs.pickle';shutil.copyfile(graph,reference/'spladder'/graph.name)
    cache_events=[]
    for f in source.glob('*.counts.hdf5'):
        for dest in [reference,native]:shutil.copyfile(f,dest/f.name)
        event_path=source/f.name.replace('.counts.hdf5','.pickle')
        shutil.copyfile(event_path,reference/event_path.name)
        with open(event_path,'rb') as stream:events=pickle.load(stream,encoding='latin1')
        cache_events.append(dict(path=str(native/event_path.name.replace('.pickle','.events.hdf5')),values=[event_value(e) for e in events]))
    request=dict(genes=[serialize(g) for g in load_genes(graph)],graph=str(native/'spladder/genes_graph_conf3.merge_graphs.hdf5'),events=cache_events)
    subprocess.run([str(args.importer)],input=json.dumps(request)+'\n',text=True,stdout=subprocess.PIPE,check=True)
    a=root/'conditionA.txt';b=root/'conditionB.txt'
    for path,order in [(a,order_a),(b,order_b)]:path.write_text(''.join(f'/align/testcase_events_1_sample{i}'+('.cram' if suffix else '.bam')+'\n' for i in order))
    variants=[
        ('baseline','testing',['--event-types','exon_skip','--dpsi','0']),
        ('four-kinds','testing_four',['--event-types','exon_skip,intron_retention,alt_3prime,alt_5prime','--dpsi','0','--out-tag','four']),
        ('options','testing.non_alt_A_vs_B_opts',['--event-types','exon_skip,intron_retention,alt_3prime,alt_5prime','--non-alt-norm','--cap-outliers','--no-cap-exp-outliers','--high-memory','--timestamp','--labelA','A','--labelB','B','--out-tag','opts']),
    ]
    for profile,output_name,extra in variants:
        flags=['-a',str(a),'-b',str(b),*extra]
        with open(root/f'reference-{profile}.log','w') as log:
            subprocess.run([str(Path(sys.executable).parent/'spladder'),'test','-o',str(reference),*flags],stdout=log,stderr=subprocess.STDOUT,check=True)
        expected_files=sorted((reference/output_name).glob('*.tsv'));assert expected_files
        for parallel in [1,4]:
            with open(root/f'native-{profile}-{parallel}.log','w') as log:
                subprocess.run([str(args.binary),'test','-o',str(native),*flags,'--parallel',str(parallel)],stdout=log,stderr=subprocess.STDOUT,check=True)
            for expected in expected_files:
                want=np.loadtxt(expected,delimiter='\t',dtype=str,ndmin=2);got=np.loadtxt(native/output_name/expected.name,delimiter='\t',dtype=str,ndmin=2)
                np.testing.assert_array_equal(got[0],want[0]);np.testing.assert_array_equal(got[1:,:6],want[1:,:6])
                np.testing.assert_allclose(got[1:,6:].astype(float),want[1:,6:].astype(float),atol=1e-8,rtol=1e-6,equal_nan=True,err_msg=f'{suffix} {parallel} {expected.name}')
                np.testing.assert_array_equal(np.round(got[1:,6:].astype(float),6),np.round(want[1:,6:].astype(float),6))
                for alpha in [0.01,0.05,0.1]:np.testing.assert_array_equal(got[1:,6:8].astype(float)<=alpha,want[1:,6:8].astype(float)<=alpha)
                snapshot=root/f'{profile}-threads-1-{expected.name}'
                if parallel==1:shutil.copyfile(native/output_name/expected.name,snapshot)
                else:assert (native/output_name/expected.name).read_bytes()==snapshot.read_bytes()
                files+=1
            cases+=1
# These are observed upstream failures/skip paths, not successful model fits.
failures=[]
for name,extra,reference_error,native_error in [
    ('singleton',['--event-types','exon_skip'], 'iteration over a 0-d array', 'singleton condition list'),
    ('empty-trend',['--event-types','mutex_exons','--dpsi','0'],'zero-size array','GLM observation shape mismatch'),
]:
    condition_a=a
    if name=='singleton':
        condition_a=root/'singleton.txt';condition_a.write_text('testcase_events_1_sample1.cram\n')
    shared=['-a',str(condition_a),'-b',str(b),'--out-tag','failure_'+name,*extra]
    ref=subprocess.run([str(Path(sys.executable).parent/'spladder'),'test','-o',str(reference),*shared],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    got=subprocess.run([str(args.binary),'test','-o',str(native),*shared],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (root/f'failure-{name}-reference.log').write_text(ref.stdout);(root/f'failure-{name}-native.log').write_text(got.stdout)
    assert ref.returncode!=0 and reference_error in ref.stdout
    assert got.returncode!=0 and native_error in got.stdout
    failures.append(name)
limit=subprocess.run([str(args.binary),'test','-o',str(native),'-a','one','-b','two','--parallel','65'],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
assert limit.returncode==2 and '1..=64' in limit.stdout
report=dict(cases=cases,text_files=files,threads=[1,4],upstream_failures=failures,cpu_limit=64,result='pass')
(args.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

#!/usr/bin/env python3
"""Compare sparse-input graph generation, including upstream multi-BAM failures."""
import argparse,contextlib,copy,io,json,pickle,shutil,subprocess
from pathlib import Path
from types import SimpleNamespace
import pysam
from spladder import settings
from spladder.init import init_genes_gtf
from spladder.spladder_prep import prep_sparse_bam_filtered
from spladder.core.gen_graphs import gen_graphs
from compare_annotations import serialize
from compare_build_graphs import config
p=argparse.ArgumentParser();p.add_argument('probe',type=Path);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
basic=a.upstream/'tests/testcase_basic/data';events=a.upstream/'tests/testcase_events/data'
sources=[(basic/f'annotation_{s}.gtf',[basic/f'align/{s}_{i}.bam' for i in [1,2]]) for s in ['pos','neg']]+[(events/'testcase_events.gtf',[events/f'align/testcase_events_1_sample{i}.bam' for i in [1,2]])]
requests=[];expected=[];failures=[]
for source_index,(gtf,original_bams) in enumerate(sources):
 root=a.work/str(source_index);root.mkdir(exist_ok=True);annotation=root/gtf.name;shutil.copyfile(gtf,annotation);bams=[]
 for original in original_bams:
  bam=root/original.name;shutil.copyfile(original,bam);pysam.index(str(bam));bams.append(bam)
 for confidence in [0,3]:
  for multiple in [False,True]:
   options=SimpleNamespace(annotation=str(annotation),verbose=False,filter_overlap_genes=False,filter_overlap_exons=False,filter_overlap_transcripts=False,confidence=confidence,readlen=50,primary_only=True,var_aware=False,ignore_mismatches=False,mm_tag='NM',ref_genome=None,logfile='-',infer_sg=False,sparse_bam=True,insert_es=True,insert_ir=True,insert_ni=True,remove_se=False,insert_intron_iterations=5,filter_consensus='',parallel=1,tmpdir=str(root),bam_fnames=list(map(str,bams)))
   settings.default_settings(options);settings.set_confidence_level(options)
   options.bam_fnames=list(map(str,bams))
   with contextlib.redirect_stdout(io.StringIO()):
    genes,options=init_genes_gtf(options);prep_sparse_bam_filtered(options)
   paths=bams if multiple else bams[:1]
   request=dict(genes=[serialize(g) for g in genes],bams=list(map(str,paths)),reference=None,options=config(options));request['options']['reads']['sparse_confidence']=confidence
   log=io.StringIO()
   try:
    with contextlib.redirect_stdout(log):result,inserted=gen_graphs(copy.deepcopy(genes),list(map(str,paths)) if multiple else str(paths[0]),options)
   except Exception as error:
    failures.append(dict(source=source_index,confidence=confidence,multiple=multiple,error=str(error),exception=type(error).__name__,request=request));continue
   (root/f'reference-C{confidence}-multi{multiple}.log').write_text(log.getvalue())
   requests.append(request);expected.append(dict(genes=[serialize(g) for g in result],inserted=inserted,read_filter=options.read_filter))
(a.work/'source-failures.json').write_text(json.dumps(failures,indent=2)+'\n')
assert requests, 'no successful upstream graph builds were exercised'
result=subprocess.run([str(a.probe)],input=''.join(json.dumps(r)+'\n' for r in requests),text=True,stdout=subprocess.PIPE,check=True)
actual=[json.loads(x) for x in result.stdout.splitlines()];assert len(actual)==len(expected)
for i,(want,got) in enumerate(zip(expected,actual)):
 if want!=got:
  (a.work/'mismatch.json').write_text(json.dumps(dict(case=i,request=requests[i],expected=want,actual=got),indent=2));raise AssertionError(f'sparse graph case {i}')
for failure in failures:
 assert failure['exception']=='IndexError'
 result=subprocess.run([str(a.probe)],input=json.dumps(failure['request'])+'\n',text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 assert result.returncode!=0 and 'sparse multi-BAM intron index' in result.stderr,result
report=dict(cases=len(requests),upstream_errors=[{k:v for k,v in x.items() if k!='request'} for x in failures],result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

#!/usr/bin/env python3
"""Positive cassette insertion from spliced BAMs, across independent genes."""
import argparse,contextlib,copy,json,os,subprocess
from pathlib import Path
from types import SimpleNamespace
import pysam
from spladder import settings
from spladder.init import init_genes_gtf
from spladder.core.gen_graphs import gen_graphs
from compare_annotations import serialize
from compare_build_graphs import config
p=argparse.ArgumentParser();p.add_argument('probe',type=Path);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
annotation=a.work/'cassette.gtf';bam=a.work/'cassette.bam';lines=[]
with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':12000}]}) as out:
 for i in range(8):
  start=100+i*1000;strand='+' if i%2==0 else '-';tags=f'gene_id "g{i}"; transcript_id "t{i}";'
  lines.append(f'chr1\tfixture\tgene\t{start+1}\t{start+250}\t.\t{strand}\t.\t{tags}')
  for lo,hi in [(start,start+50),(start+200,start+250)]:lines.append(f'chr1\tfixture\texon\t{lo+1}\t{hi}\t.\t{strand}\t.\t{tags}')
  for j in range(20):
   read=pysam.AlignedSegment();read.query_name=f'g{i}-r{j}';read.query_sequence='A'*150;read.flag=0;read.reference_id=0;read.reference_start=start;read.mapping_quality=60;read.cigartuples=[(0,50),(3,50),(0,50),(3,50),(0,50)];read.set_tag('NM',0);read.set_tag('XS',strand);out.write(read)
annotation.write_text('\n'.join(lines)+'\n');pysam.index(str(bam));requests=[];expected=[]
for confidence in range(4):
 o=SimpleNamespace(annotation=str(annotation),verbose=False,filter_overlap_genes=False,filter_overlap_exons=False,filter_overlap_transcripts=False,confidence=confidence,readlen=150,primary_only=True,var_aware=False,ignore_mismatches=False,mm_tag='NM',ref_genome=None,logfile='-',infer_sg=False,sparse_bam=False,insert_es=True,insert_ir=True,insert_ni=True,remove_se=False,insert_intron_iterations=5,filter_consensus='')
 settings.default_settings(o);settings.set_confidence_level(o)
 with (a.work/f'reference-C{confidence}.log').open('w') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
  genes,o=init_genes_gtf(o);requests.append(dict(genes=[serialize(g) for g in genes],bams=[str(bam)],reference=None,options=copy.deepcopy(config(o))));result,inserted=gen_graphs(genes,[str(bam)],o)
 assert inserted['cassette_exon']==8,inserted
 expected.append(dict(genes=[serialize(g) for g in result],inserted=inserted,read_filter=o.read_filter))
for threads in [1,4]:
 run=subprocess.run([str(a.probe)],input=''.join(json.dumps(r)+'\n' for r in requests),text=True,stdout=subprocess.PIPE,check=True,env={**os.environ,'RAYON_NUM_THREADS':str(threads)})
 actual=[json.loads(line) for line in run.stdout.splitlines()]
 if actual!=expected:
  (a.work/'mismatch.json').write_text(json.dumps(dict(expected=expected,actual=actual),indent=2));raise AssertionError(f'cassette build mismatch ({threads} threads)')
report=dict(cases=8,genes_per_case=8,inserted_cassettes=64,threads=[1,4],result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

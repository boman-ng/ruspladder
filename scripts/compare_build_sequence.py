#!/usr/bin/env python3
"""Exercise >200-intron feasibility and shared options across sample graph builds."""
import argparse,contextlib,copy,io,json,subprocess
from pathlib import Path
from types import SimpleNamespace
import pysam
from spladder import settings
from spladder.init import init_genes_gtf
from spladder.core.gen_graphs import gen_graphs
from compare_annotations import serialize
from compare_build_graphs import config
parser=argparse.ArgumentParser();parser.add_argument('probe',type=Path);parser.add_argument('--work',type=Path,required=True)
a=parser.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
gtf=a.work/'annotation.gtf'
gtf.write_text('chr1\tfixture\tgene\t101\t400\t.\t+\t.\tgene_id "g"; gene_type "protein_coding";\n'+''.join(f'chr1\tfixture\texon\t{s}\t{e}\t.\t+\t.\tgene_id "g"; transcript_id "t";\n' for s,e in [(101,200),(301,400)]))
bams=[]
for sample,support in enumerate([2,4]):
    rows=[]
    for donor in range(200,401):
        for _ in range(support):rows.append((donor-20,[(0,20),(3,800-donor),(0,20)]))
    for start in [100,300]:
        for _ in range(60):rows.append((start,[(0,100)]))
    path=a.work/f'sample{sample}.bam';bams.append(path)
    with pysam.AlignmentFile(path,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000}]}) as out:
        for i,(start,cigar) in enumerate(sorted(rows)):
            read=pysam.AlignedSegment();read.query_name=f'r{i}';read.flag=0;read.reference_id=0;read.reference_start=start;read.mapping_quality=60;read.cigartuples=cigar
            read.query_sequence='A'*sum(length for op,length in cigar if op==0);read.query_qualities=pysam.qualitystring_to_array('I'*len(read.query_sequence));read.set_tag('NM',0);read.set_tag('XS','+');out.write(read)
    pysam.index(str(path))
o=SimpleNamespace(annotation=str(gtf),verbose=False,filter_overlap_genes=False,filter_overlap_exons=False,filter_overlap_transcripts=False,confidence=3,readlen=50,primary_only=True,var_aware=False,ignore_mismatches=False,mm_tag='NM',ref_genome=None,logfile='-',infer_sg=False,sparse_bam=False,insert_es=False,insert_ir=True,insert_ni=False,remove_se=False,insert_intron_iterations=5,filter_consensus='')
settings.default_settings(o);settings.set_confidence_level(o)
with contextlib.redirect_stdout(io.StringIO()):genes,o=init_genes_gtf(o)
request=dict(genes=[serialize(g) for g in genes],samples=[[str(p)] for p in bams],options=config(o))
expected=[];filters=[];log=io.StringIO()
with contextlib.redirect_stdout(log):
    for path in bams:
        generated,inserted=gen_graphs(copy.deepcopy(genes),str(path),o)
        expected.append(dict(genes=[serialize(g) for g in generated],inserted=inserted,read_filter=o.read_filter.copy()))
        filters.append(dict(current=o.read_filter.copy(),retention=o.intron_retention['read_filter'].copy()))
(a.work/'reference.log').write_text(log.getvalue());(a.work/'filters.json').write_text(json.dumps(filters,indent=2)+'\n')
run=subprocess.run([str(a.probe)],input=json.dumps(request)+'\n',text=True,stdout=subprocess.PIPE,check=True);actual=json.loads(run.stdout)
if actual!=expected:
    (a.work/'mismatch.json').write_text(json.dumps(dict(request=request,expected=expected,actual=actual),indent=2)+'\n');raise AssertionError('sequential graph generation mismatch')
assert [x['current']['exon_len'] for x in filters]==[17,21]
assert [x['retention']['exon_len'] for x in filters]==[17,17]
report=dict(samples=len(bams),introns=201,filter_states=filters,result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

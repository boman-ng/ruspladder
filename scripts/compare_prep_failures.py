#!/usr/bin/env python3
"""Record prep CLI defects and verify CRAM output against its working core API."""
import argparse,contextlib,io,json,pickle,shutil,subprocess,sys
from pathlib import Path
from types import SimpleNamespace
import h5py,pysam,numpy as np
from spladder import settings
from spladder.reads import summarize_chr
from compare_count_io import compare_hdf5
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('--upstream',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
source=a.upstream/'tests/testcase_events/data';ref=source/'genome.fa';results={}
for binary,name in [(Path(sys.executable).parent/'spladder','reference'),(a.binary,'native')]:
 root=a.work/name;root.mkdir(exist_ok=True);annotation=root/'annotation.gtf';shutil.copyfile(source/'testcase_events.gtf',annotation);bam=root/'sample.cram';shutil.copyfile(source/'align/testcase_events_1_sample1.cram',bam);shutil.copyfile(source/'align/testcase_events_1_sample1.cram.crai',Path(str(bam)+'.crai'))
 command=[str(binary),'prep','-a',str(annotation),'-b',str(bam),'--sparse-bam','--reference',str(ref),'--parallel','1']
 result=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(root/'command.log').write_text(result.stdout);results[name]=result
assert results['reference'].returncode!=0 and "has no attribute 'cram_ref'" in results['reference'].stdout
assert results['native'].returncode==0,results['native'].stdout
with (a.work/'reference/annotation.gtf.pickle').open('rb') as f:genes=pickle.load(f)
chromosomes=sorted({g.chr for g in genes});o=SimpleNamespace(verbose=False,primary_only=True,var_aware=False,mm_tag='NM',cram_ref=str(ref),confidence=3,readlen=50);settings.default_settings(o);settings.set_confidence_level(o)
arrays=0;oracle=a.work/'core-oracle';oracle.mkdir(exist_ok=True)
for filtered in [True,False]:
 output=oracle/('sample.conf_3.filt.hdf5' if filtered else 'sample.hdf5')
 with h5py.File(output,'w') as out:
  for chromosome in chromosomes:
   with contextlib.redirect_stdout(io.StringIO()):_,coo,minus,plus=summarize_chr(str(a.work/'reference/sample.cram'),chromosome,o,filter=o.read_filter if filtered else None)
   values=dict(reads_row=coo.row.astype('uint8'),reads_col=coo.col,reads_dat=coo.data,reads_shp=np.array(coo.shape),introns_m=minus,introns_p=plus)
   for key,value in values.items():out.create_dataset(chromosome+'_'+key,data=value,**({} if key=='reads_shp' else dict(compression='gzip')))
 arrays+=compare_hdf5(output,a.work/'native'/output.name)
# summarize_chr ignores --ignore-mismatches even though direct BAM reads honor it.
basic=a.upstream/'tests/testcase_basic/data'
for binary,name in [(Path(sys.executable).parent/'spladder','reference'),(a.binary,'native')]:
 root=a.work/(name+'-no-nm');root.mkdir(exist_ok=True);annotation=root/'annotation.gtf';shutil.copyfile(basic/'annotation_pos.gtf',annotation);bam=root/'without-nm.bam'
 with pysam.AlignmentFile(str(basic/'align/pos_1.bam'),'rb') as src,pysam.AlignmentFile(str(bam),'wb',template=src) as out:
  for read in src:read.set_tag('NM',None);out.write(read)
 pysam.index(str(bam))
 result=subprocess.run([str(binary),'prep','-a',str(annotation),'-b',str(bam),'--sparse-bam','--ignore-mismatches'],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (root/'command.log').write_text(result.stdout);assert result.returncode!=0 and 'NM' in result.stdout
report=dict(cram_core_arrays=arrays,upstream_cli_errors=['CRAM uses missing cram_ref attribute','sparse prep ignores ignore-mismatches'],native_cram_reference_forwarded=True,result='pass');(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

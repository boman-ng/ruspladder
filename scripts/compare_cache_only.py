#!/usr/bin/env python3
"""Re-quantify from completed sparse caches after the original BAMs are absent."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
from compare_lifecycle import compare
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('exporter',type=Path);p.add_argument('--source',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--regenerate',action='store_true');a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True);(a.work/'report.json').unlink(missing_ok=True)
for binary,name,original in [(Path(sys.executable).parent/'spladder','reference','reference'),(a.binary,'native','native-1')]:
 dest=a.work/name
 assert not dest.exists(),'use a fresh fixture directory'
 shutil.copytree(a.source/original,dest)
 bams=sorted((dest/'inputs').glob('*.bam'))
 for path in bams:
  path.unlink();Path(str(path)+'.bai').unlink()
 if a.regenerate:shutil.rmtree(dest/'results')
 else:
  for suffix in ['.count.hdf5','.gene_exp.hdf5']:(dest/'results/spladder'/('genes_graph_conf3.merge_graphs'+suffix)).unlink()
 annotation=next((dest/'inputs').glob('*.gtf'))
 command=[str(binary),'build','-a',str(annotation),'-b',','.join(map(str,bams)),'-o',str(dest/'results'),'--sparse-bam','--parallel','4']
 result=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(dest/'cache-only.log').write_text(result.stdout)
 assert result.returncode==0,(name,result.stdout[-1600:])
report=compare(a.work/'reference/results',a.work/'native/results',a.exporter);(a.work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS:',report)

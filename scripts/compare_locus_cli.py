#!/usr/bin/env python3
"""Check locus normalization, complete Python parity, reuse and mode isolation."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
from audit_locus_gtf import audit
from compare_result_dirs import compare

p = argparse.ArgumentParser()
p.add_argument('binary', type=Path)
p.add_argument('exporter', type=Path)
p.add_argument('--upstream', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
assert not (a.work/'native').exists(), 'use a fresh test directory'
source = a.upstream/'tests/testcase_events/data'
# Reuse real upstream read evidence on both strands. Only the IDs of gene7
# are changed to collide with gene1; its positions and transcripts stay intact.
text = (source/'testcase_events.gtf').read_text()
text = text.replace('gene_id "gene7";', 'gene_id "gene1";').replace('transcript_id "gene7.', 'transcript_id "gene1.')
annotation = a.work/'colliding.gtf'
annotation.write_text(text)
normalized = Path(str(annotation)+'.locus.gtf')
common = ['build','--bams',','.join(str(source/f'align/testcase_events_1_sample{i}.bam') for i in [1,2]),
          '--parallel','4','--reference',str(source/'genome.fa'),'--output-txt']
native = [str(a.binary),*common,'-a',str(annotation),'-o',str(a.work/'native'),'--annotation-mode','locus']
with (a.work/'native.log').open('w') as log:
    subprocess.run(native, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=3600)
audit_report = audit(annotation, normalized)
reference = [str(Path(sys.executable).parent/'spladder'),*common,'-a',str(normalized),'-o',str(a.work/'reference')]
with (a.work/'reference.log').open('w') as log:
    subprocess.run(reference, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=3600)
comparison = compare(a.work/'reference', a.work/'native', a.exporter)
before = {str(f):f.stat().st_mtime_ns for f in a.work.rglob('*') if f.is_file()}
result = subprocess.run(native, capture_output=True, check=True, timeout=3600)
assert before == {str(f):f.stat().st_mtime_ns for f in a.work.rglob('*') if f.is_file()}, 'cache reuse changed artifacts'
result = subprocess.run(native[:-2], capture_output=True, text=True)
assert result.returncode and 'different annotation mode' in result.stderr, 'mode isolation'
report = dict(result='pass', annotation=audit_report, comparison=comparison, reuse='pass', mode_isolation='pass')
(a.work/'report.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)

#!/usr/bin/env python3
"""Compare complete direct-BAM gen_graphs calls, from annotation to augmented graph."""
import argparse
import contextlib
import copy
import io
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
import pysam

from spladder import settings
from spladder.init import init_genes_gtf
from spladder.core.gen_graphs import gen_graphs
from compare_annotations import serialize


def config(options):
    return dict(reads=dict(filter=options.read_filter.copy(), primary_only=options.primary_only,
                           var_aware=options.var_aware, no_mm=options.ignore_mismatches, mm_tag=options.mm_tag),
                cassette=options.cassette_exon,
                retention={k: v for k, v in options.intron_retention.items() if k != "read_filter"},
                intron_edges={k: bool(v) if k in ["insert_intron_retention", "gene_merges", "append_new_terminal_exons"] else v
                              for k, v in options.intron_edges.items()},
                remove_exons=options.remove_exons,
                insert_es=options.insert_es, insert_ir=options.insert_ir, insert_ni=options.insert_ni,
                remove_se=options.remove_se, introns_unstranded=options.introns_unstranded,
                insert_intron_iterations=options.insert_intron_iterations,
                consensus=None if not options.filter_consensus else options.filter_consensus == "lenient")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--airway", type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "report.json").unlink(missing_ok=True)
    basic = args.upstream / "tests/testcase_basic/data"
    events = args.upstream / "tests/testcase_events/data"
    sources = [(basic / f"annotation_{s}.gtf", [basic / f"align/{s}_{i}.bam" for i in [1, 2]], basic / f"genome_{s}.fa") for s in ["pos", "neg"]]
    sources.append((events / "testcase_events.gtf", [events / f"align/testcase_events_1_sample{i}.bam" for i in [1, 2]], events / "genome.fa"))
    if args.airway:
        sources.append((args.airway / "Homo_sapiens.GRCh37.75_subset.gtf",
                        [args.airway / f"SRR103950{i}_subset.bam" for i in [8, 9]], None))
    scenarios = [dict(confidence=i) for i in range(4)]
    scenarios += [dict(insert_es=False), dict(insert_ir=False), dict(insert_ni=False),
                  dict(insert_es=False, insert_ir=False, insert_ni=False), dict(remove_se=True),
                  dict(introns_unstranded=True), dict(filter_consensus="strict"), dict(filter_consensus="lenient")]
    requests, expected = [], []
    totals = {}
    for source_index, (annotation, bams, reference) in enumerate(sources):
        for bam in bams:
            if not Path(str(bam) + ".bai").exists(): pysam.index(str(bam))
        local = args.work / f"annotation-{source_index}.gtf"
        shutil.copyfile(annotation, local)
        source_scenarios = scenarios + [dict(multiple=True), dict(multiple=True, confidence=0)] if reference else [dict(confidence=0), dict(confidence=3), dict(multiple=True)]
        for changes in source_scenarios:
            options = SimpleNamespace(annotation=str(local), verbose=False, filter_overlap_genes=False,
                filter_overlap_exons=False, filter_overlap_transcripts=False, confidence=3, readlen=50,
                primary_only=True, var_aware=False, ignore_mismatches=False, mm_tag="NM", ref_genome=str(reference) if reference else None,
                logfile="-", infer_sg=False, sparse_bam=False, insert_es=True, insert_ir=True, insert_ni=True,
                remove_se=False, insert_intron_iterations=5, filter_consensus="")
            settings.default_settings(options)
            if reference is None:
                options.readlen = 101
                # The distributed airway alignments have no NM tags. Use the
                # documented upstream option on BOTH sides; never invent tags.
                options.ignore_mismatches = True
            for k, v in changes.items(): setattr(options, k, v)
            settings.set_confidence_level(options)
            selected = bams if changes.get("multiple") else bams[:1]
            with contextlib.redirect_stdout(io.StringIO()):
                genes, options = init_genes_gtf(options)
            requests.append(dict(genes=[serialize(g) for g in genes], bams=list(map(str, selected)),
                                 reference=str(reference) if reference else None, options=copy.deepcopy(config(options))))
            with (args.work / f"reference-{len(requests) - 1}.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                result, inserted = gen_graphs(genes, list(map(str, selected)), options)
            expected.append(dict(genes=[serialize(g) for g in result], inserted=inserted, read_filter=options.read_filter))
            for k, v in inserted.items(): totals[k] = totals.get(k, 0) + int(v)
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"complete graph generation mismatch: {path}")
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), inserted=totals,
        airway_ignore_mismatches=bool(args.airway), result="pass"), indent=2) + "\n")
    print(f"PASS: {len(actual)} complete direct-BAM graph generation comparisons; {totals}")


if __name__ == "__main__": main()

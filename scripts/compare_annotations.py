#!/usr/bin/env python3
"""Compare annotations, filters, gene ordering and graph cache roundtrips."""
import argparse
import contextlib
import io
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from spladder.init import init_genes_gtf, init_genes_gff3


def serialize(gene):
    result = {key: getattr(gene, key) for key in
              ("name", "start", "stop", "chr", "strand", "source", "gene_type", "symbol",
               "transcripts", "is_alt", "is_alt_spliced")}
    result["start"] = int(result["start"])
    result["stop"] = int(result["stop"])
    for key in ["is_alt", "is_alt_spliced"]:
        if result[key] is not None: result[key] = bool(result[key])
    result["exons"] = [e.tolist() for e in gene.exons]
    result["introns_anno"] = [[int(v) for v in x] for x in sorted(gene.introns_anno)]
    sg = gene.splicegraph
    result["splicegraph"] = {"vertices": sg.vertices.T.tolist(),
                             "edges": [np.where(row)[0].tolist() for row in sg.edges],
                             "terminals": sg.terminals.T.astype(bool).tolist()}
    sg = gene.segmentgraph
    result["segmentgraph"] = {"segments": sg.segments.T.tolist(),
                              "matches": [np.where(row)[0].tolist() for row in sg.seg_match],
                              "edges": np.array(np.where(sg.seg_edges)).T.tolist()}
    return result


def make_gtf(path):
    # Keep a separate non-overlapping gene so all filter combinations are valid.
    path.write_text("""chr1\tsource\tgene\t1\t100\t.\t+\t.\tgene_id "a"; gene_name "A";
chr1\tsource\tgene\t80\t200\t.\t+\t.\tgene_id "b";
chr1\tsource\tgene\t300\t400\t.\t-\t.\tgene_id "empty";
chr1\tsource\texon\t1\t20\t.\t+\t.\tgene_id "a"; transcript_id "a1";
chr1\tsource\texon\t60\t100\t.\t+\t.\tgene_id "a"; transcript_id "a1";
chr1\tsource\texon\t80\t110\t.\t+\t.\tgene_id "b"; transcript_id "b1";
chr1\tsource\texon\t100\t200\t.\t+\t.\tgene_id "b"; transcript_id "b1";
chr1\tsource\texon\t1\t20\t.\t-\t.\tgene_id "shared"; transcript_id "s1"; gene_name "not_inferred";
chr1\tsource\texon\t501\t520\t.\t+\t.\tgene_id "safe"; transcript_id "safe1";
chr1\tsource\texon\t551\t600\t.\t+\t.\tgene_id "safe"; transcript_id "safe1";
""")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--extra", type=Path, action="append", default=[])
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    fixtures = []
    for source in sorted((args.upstream / "tests").glob("*/data/*.gtf")):
        target = args.work / f"{source.parent.parent.name}-{source.name}"
        shutil.copyfile(source, target)
        fixtures.append(target)
    for source in args.extra:
        target = args.work / source.name
        shutil.copyfile(source, target)
        fixtures.append(target)
    special = args.work / "filters.gtf"
    make_gtf(special)
    gff = args.work / "filters.gff3"
    gff.write_text("""##gff-version 3
chr2\tsource\tgene\t101\t300\t.\t-\t.\tID=g;gene_name=G
chr2\tsource\tmRNA\t101\t300\t.\t-\t.\tID=t;Parent=g
chr2\tsource\texon\t101\t150\t.\t-\t.\tID=e1;Parent=t
chr2\tsource\texon\t201\t300\t.\t-\t.\tID=e2;Parent=t
##FASTA
>chr2
ACGT
""")
    tasks = [(p, 0) for p in fixtures] + [(special, flags) for flags in range(8)] + [(gff, 0)]
    for index, (path, flags) in enumerate(tasks):
        # Each run has its own annotation, pickle, and exclusion reports.
        local = args.work / f"case-{index}{path.suffix}"
        # Reusing this test directory must not inherit reports from a former
        # case occupying the same index (e.g. when adding the real-data GTF).
        for suffix in ["no_exons", "gene_overlap", "exon_shared", "exon_overlap"]:
            Path(str(local) + ".genes_excluded_" + suffix).unlink(missing_ok=True)
        shutil.copyfile(path, local)
        options = SimpleNamespace(annotation=str(local), verbose=False,
                                  filter_overlap_genes=bool(flags & 1),
                                  filter_overlap_exons=bool(flags & 2),
                                  filter_overlap_transcripts=bool(flags & 4))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            genes, _ = (init_genes_gff3 if local.suffix == ".gff3" else init_genes_gtf)(options)
        excluded = []
        for suffix in ["no_exons", "gene_overlap", "exon_shared", "exon_overlap"]:
            report = Path(str(local) + ".genes_excluded_" + suffix)
            if report.exists(): excluded.append([suffix, report.read_text().splitlines()])
        expected = {"genes": [serialize(g) for g in genes], "excluded": excluded}
        process = subprocess.run([str(args.probe), str(local), str(args.work / f"case-{index}.hdf5"), str(flags)],
                                 text=True, capture_output=True, check=True)
        actual = json.loads(process.stdout)
        if actual != expected:
            failure = args.work / "mismatch.json"
            failure.write_text(json.dumps({"path": str(local), "flags": flags, "expected": expected, "actual": actual}, indent=2))
            raise AssertionError(f"annotation mismatch: {failure}")
    print(f"PASS: {len(tasks)} annotation/filter cases and HDF5 roundtrips agree exactly")


if __name__ == "__main__":
    main()

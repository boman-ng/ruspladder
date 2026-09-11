#!/usr/bin/env python3
"""Compare ambiguity and FASTA splice-site consensus filtering."""
import argparse
import contextlib
import io
import json
import random
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from spladder.classes.gene import Gene
from spladder.helpers import filter_introns, filter_introns_consensus
from compare_annotations import serialize


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("probe", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    reference = args.work / "reference.fa"
    sequence = bytearray(b"ACGT" * 100)
    for start, stop, donor, acceptor in [(20, 60, b"GT", b"AG"), (80, 140, b"GC", b"AG"),
                                       (160, 220, b"CT", b"AC"), (240, 280, b"CT", b"GC"),
                                       (300, 360, b"gt", b"ag")]:
        sequence[start:start + 2], sequence[stop - 2:stop] = donor, acceptor
    reference.write_text("".join(f">chr{i}\n{sequence.decode()}\n" for i in [1, 2]))
    Path(str(reference) + ".fai").unlink(missing_ok=True)
    rng = random.Random(761)
    requests, expected = [], []
    for trial in range(20):
        genes = []
        for chromosome in ["chr1", "chr2"]:
            for strand in ["+", "-"]:
                for start, stop in [[10, 100], [100, 180], [90, 250], [280, 380]]:
                    genes.append(Gene(name=f"g{len(genes)}", start=start, stop=stop, chr=chromosome, strand=strand))
        pairs = [[20, 60], [80, 140], [160, 220], [240, 280], [300, 360], [100, 180]]
        pairs += [sorted(rng.sample(range(2, 395), 2)) for _ in range(20)]
        lists = [[[a, b, rng.randrange(1, 20)] for a, b in pairs] for _ in range(2)]
        raw = [lists for _ in genes]
        for consensus, lenient, offset in [(False, False, 0), (False, False, 30), (True, False, None), (True, True, None), (True, True, 5)]:
            request = dict(genes=[serialize(g) for g in genes], introns=raw, offset=offset,
                           reference=str(reference) if consensus else None, lenient=lenient)
            requests.append(request)
            introns = np.empty((len(genes), 2), dtype=object)
            for i in range(len(genes)):
                for si in range(2): introns[i, si] = np.array(raw[i][si], dtype=np.int64)
            options = SimpleNamespace(ref_genome=str(reference), filter_consensus="lenient" if lenient else "strict",
                                      intron_edges=dict(append_new_terminal_exons_len=offset))
            with contextlib.redirect_stdout(io.StringIO()):
                if consensus: introns = filter_introns_consensus(introns, genes, options)
                if offset is not None: introns = filter_introns(introns, genes, options)
            expected.append([[row.tolist() for row in lists] for lists in introns])
    run = subprocess.run([str(args.probe)], input="".join(json.dumps(r) + "\n" for r in requests), text=True, stdout=subprocess.PIPE, check=True)
    actual = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(actual) == len(expected)
    for i, (want, got) in enumerate(zip(expected, actual)):
        if want != got:
            path = args.work / "mismatch.json"
            path.write_text(json.dumps(dict(case=i, request=requests[i], expected=want, actual=got), indent=2))
            raise AssertionError(f"intron-filter mismatch: {path}")
    (args.work / "report.json").write_text(json.dumps(dict(cases=len(actual), result="pass"), indent=2) + "\n")
    print(f"PASS: {len(actual)} intron ambiguity and reference-consensus filter cases")


if __name__ == "__main__": main()

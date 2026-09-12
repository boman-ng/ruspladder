# Performance evidence

The acceptance workload uses four logical CPUs and a hard 8 GiB cgroup memory
limit, with swap disabled. The local Docker image is pinned by digest; it runs
the same reference environment and native libraries used for correctness checks.
Every measured command records its actual affinity and cgroup limit. Up to 64
threads are accepted by the CLI; the performance baseline stays at four.

`scripts/run_lifecycle_benchmarks.py` runs three interleaved repetitions of:

- Eight real airway BAMs: the complete build, including annotation, all sample
  graphs, merge, quantification, six event types and public outputs.
- Twenty upstream event BAMs: complete build followed by exon-skip differential
  testing, with ten samples per condition, from direct alignment input.
- The same twenty-sample build/test workflow with sparse alignment preparation
  and sparse input enabled.

Fresh means new annotation/output caches. Reused means rerunning the identical
successful command with its completed caches. Each measurement gets a separate
container so reused-run memory peaks cannot include the earlier fresh build.
Filesystem caches are not flushed; this is not a cold-disk benchmark. Input
staging and container startup precede timing. Child command startup, all their
I/O and completion are included. Cgroup peaks include the lightweight measuring
worker and staging. BLAS and OpenMP use one thread; Rayon and NumPy/Numba's
upstream worker setting use four. NUMBA_CACHE_DIR is shared, but the pinned
event detectors use `@jit(nopython=True)` without persistent caching. Their
compilation cost remains in fresh reference builds; the source is unmodified.
This startup cost has a large effect on these small subsets. Each build and test
stage has its own wall/CPU timing.

CPU seconds include child workers; mean cores = CPU seconds / wall seconds.
Percentage utilization is mean cores / 4 * 100. RSS is the largest individual
child process, while the cgroup peak measures charged memory across the container,
including charged file cache. Shared pages are charged to the cgroup that first
touches them, so process RSS can exceed a container's charged-memory counter;
these are distinct metrics. See the [Linux memory-controller documentation](https://docs.kernel.org/admin-guide/cgroup-v1/memory.html#shared-page-accounting).
`input_blocks` and `output_blocks` are Linux getrusage I/O
counts (512-byte blocks). Each summary metric is independently reduced to the median of three runs.
Raw commands, binary SHA-256, versions, output parity
and every repetition are retained in the report.

The airway subset has too few confirmed events for the source Gamma trend fit;
both test implementations fail there. That failed build/test pilot is excluded
from successful timing claims. The twenty-sample event fixture completes the
full build/test workflow and supplies the statistical end-to-end comparison.
These small subsets do not establish whole-genome or 64-core scaling.

## Adopted optimization investigations

- [Bounded chromosome coverage](RESEARCH.md): replaces chromosome-wide dense
  allocation with indexed windows and difference accumulation. The source's
  four-contig, 294-byte BAM reproducer OOMs and leaves a hung worker pool under
  8 GiB. Native completes it, and 96 summaries / 1,008 arrays compare exactly.
  Three real-chromosome runs reduce roughly 1.93 GiB cgroup peak to 19–20 MiB.
- Explicit HDF5 chunks: separate `perf/p0-hdf5-chunks` branch reproduces 64 MiB
  automatic chunks for small unlimited sample matrices. Bounded chunks retain
  logical outputs and reduce its local 7.7-second pathology to about 9 ms.
  That microbenchmark uses an address-space cap, not the lifecycle cgroup.
- [BAM decompression and independent gene stages](DECOMPRESSION.md): reuse
  HTSlib's libdeflate feature and Rayon for cassette/retention stages. Relative
  to the earlier native implementation, median real-airway build time falls
  38.2%, CPU seconds fall 29.6%, and mean cores rise from 0.93 to 1.05. This
  carries a measured 2.9 MiB peak-memory cost. The whole-pipeline comparison
  against SplAdder, including overall memory improvement, is reported separately.

## Final lifecycle results

Measured 2026-09-12 (Asia/Shanghai), native source commit `2049531c2bcc4a79d61e27489155127fff7344a5`.
Release built on Linux x86_64 with rustc 1.94.1 (`e408947bf`, 2026-03-25).
Native binary SHA-256: `5250cc8f6afeaa05e5e3d4a1b2adde1958a4da154d3a6fd1683c2d5a7f3b7c41`.

All 36 measured runs complete. Eighteen source/native comparisons pass: 306 graph
comparisons, 108 event collections, 2,100 public HDF5 arrays, 192 scientific text
files and 36 differential-test TSV comparisons. These totals include reused-cache
rechecks. Discrete values/order/decisions are exact; quantification uses
atol=1e-10/rtol=1e-8, statistics atol=1e-8/rtol=1e-6 and exact six-decimal results.
The separate release CLI suite verifies raw floating-point bits at 1/4 threads.

| Workflow | Cache | SplAdder wall s | Rust wall s | Speedup | SplAdder peak MiB | Rust peak MiB | Peak reduction |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Airway build | fresh | 36.1212 | 1.0735 | 33.65× | 830.59 | 27.43 | 96.70% |
| Airway build | reused | 1.4320 | 0.0492 | 29.10× | 149.02 | 17.57 | 88.21% |
| 20-sample direct build/test | fresh | 44.8325 | 1.8410 | 24.35× | 830.77 | 27.06 | 96.74% |
| 20-sample direct build/test | reused | 3.2656 | 0.0553 | 59.03× | 170.96 | 23.95 | 85.99% |
| 20-sample sparse build/test | fresh | 35.3547 | 2.5144 | 14.06× | 838.38 | 29.35 | 96.50% |
| 20-sample sparse build/test | reused | 3.2202 | 0.0555 | 58.01× | 170.89 | 23.88 | 86.02% |

| Workflow | Cache | CPU s source / Rust | Mean cores source / Rust | Four-CPU utilization source / Rust | RSS MiB source / Rust | I/O blocks in,out source / Rust |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Airway build | fresh | 111.9703 / 0.9954 | 3.100 / 0.927 | 77.5% / 23.2% | 299.02 / 26.13 | 128,2960 / 32,4920 |
| Airway build | reused | 1.4015 / 0.0494 | 0.972 / 0.994 | 24.3% / 24.9% | 245.81 / 17.86 | 48,24 / 0,0 |
| 20-sample direct build/test | fresh | 120.4555 / 1.8577 | 2.687 / 1.009 | 67.2% / 25.2% | 298.48 / 25.80 | 472,2000 / 0,7112 |
| 20-sample direct build/test | reused | 3.2884 / 0.0675 | 1.007 / 1.221 | 25.2% / 30.5% | 249.41 / 25.64 | 0,96 / 0,56 |
| 20-sample sparse build/test | fresh | 108.8417 / 4.1134 | 3.079 / 1.627 | 77.0% / 40.7% | 300.97 / 26.15 | 0,7200 / 0,9888 |
| 20-sample sparse build/test | reused | 3.2464 / 0.0681 | 1.007 / 1.210 | 25.2% / 30.3% | 250.32 / 25.57 | 0,96 / 0,56 |

Native CPU time is much lower, but fresh-run average CPU occupancy is also lower
than the source. Do not label the source-relative occupancy numbers an improvement:
upstream JIT compilation is included, while the native executable is precompiled.
The controlled native-before/after experiment in DECOMPRESSION.md measures the
12.3% mean-core increase from independent gene stages. These small workflows do
not saturate four cores throughout their lifetime. Native HDF5 caches also write
more blocks in fresh runs than source pickle caches; no I/O-volume reduction is
claimed. Runtime and peak-memory improvements pass on the small real-data
airway workload above; this does not establish full-size sample completion.

| Workflow | Cache | Build wall s source / Rust | Test wall s source / Rust |
| --- | --- | ---: | ---: |
| Airway build | fresh | 36.1212 / 1.0735 | — / — |
| Airway build | reused | 1.4320 / 0.0492 | — / — |
| 20-sample direct build/test | fresh | 42.8632 / 1.8253 | 1.8858 / 0.0156 |
| 20-sample direct build/test | reused | 1.4027 / 0.0414 | 1.8629 / 0.0139 |
| 20-sample sparse build/test | fresh | 33.5063 / 2.4996 | 1.8483 / 0.0152 |
| 20-sample sparse build/test | reused | 1.3809 / 0.0414 | 1.8398 / 0.0144 |

Raw evidence: [lifecycle report](../runs/lifecycle-final/report.json). The per-run
reports include all stage CPU times, limits, I/O, hashes and commands. Earlier
failed pilots remain under ../runs/lifecycle-pilot-* and are not included here.

Reproduce from the prepared environment (use a new output directory):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /home/wubw/data/ruspladder/envs/reference/bin/python \
  scripts/run_lifecycle_benchmarks.py \
  --work /home/wubw/data/ruspladder/runs/lifecycle-repeat \
  --binary /home/wubw/data/ruspladder/target/release/ruspladder \
  --exporter /home/wubw/data/ruspladder/target/release/examples/cache_export_probe
```


## Full-size lncRNA follow-up, 2026-09-12

The three BAMs from the previously timed-out Python workflow were run separately
with the same scientific parameters: RefSeq hg19, confidence 3, read length 150,
`merge_graphs`, `exon_skip,mult_exon_skip`, and text outputs. Each native process
used four disjoint CPUs and a hard 8 GiB container limit. Shared annotation prep
took 26.235 seconds; the build durations below exclude this shared step. The user
set a 20-minute deadline while builds were running, applied from each build's
original start time. The supervisor sent SIGTERM at that deadline and retained
all files. The small scheduling/signal delay is included in measured wall time.

| Sample | BAM GB (decimal) | Build wall seconds | Mean cores | Peak process RSS GiB | Outcome at deadline |
| --- | ---: | ---: | ---: | ---: | --- |
| S025 | 11.249 | 1200.893 | 1.295 | 1.876 | Timed out during graph generation |
| S026 | 9.193 | 1201.373 | 1.252 | 2.198 | Timed out during graph generation |
| S027 | 10.674 | 1201.884 | 2.009 | 2.105 | Graph saved; timed out during graph quantification |

All three containers reached 8 GiB charged peak, including file cache. None had
an OOM kill. Charged peak is not process RSS. None produced complete event count
or event text outputs. S027 saved its sample and merged graph caches; its
unfinished temporary count file is not a usable completed count result.

The S027 merged graph contains 31,130 genes, all matching the earlier completed
Python graph under the existing graph serializer, including graph structure,
ordering, annotation fields, and segment graph. The earlier and current BAMs
have identical per-contig index counts but different bytes; this is a comparison
of observed outputs, not a claim of input byte identity. No complete event or
quantification parity result is available for these timed-out runs.

These are execution observations, not an isolated Python/Rust speedup benchmark:
the samples shared storage, and separate CPU sets also ran the bounded P0
investigation. The Python retries timed out around 59 minutes; older completed
Python runs took 9–10.5 hours under their earlier resource configuration. Neither
is used to calculate a speedup ratio here.

The serial initial intron-query loop was evaluated on branch
`perf/p0-intron-extraction`. The candidate passed the existing source and thread
comparisons and reduced median chr22-subset runtime by 18.6%, but its full S026
run also timed out at 20 minutes after saving only the sample graph. It was not
merged; the production binary remains the tested baseline. See the branch
research note for the bounded experiment and its limitations.
Evidence and commands: [run manifest](../runs/lncrna-three-20260912/manifest.json),
[baseline outcomes](../runs/lncrna-three-20260912/baseline-summary.json), and
[S027 graph comparison](../runs/lncrna-three-20260912/S027_HK20260811047RNA-5_rna_LncRNA_2556568/graph-comparison-previous.json).

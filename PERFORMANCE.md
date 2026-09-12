# Performance evidence

The current acceptance baseline is **8 CPUs / 16 GiB**, with no additional
swap, a 1,800-second soft notification and a 3,600-second hard SIGKILL boundary.
The [HPC investigation](HPC_LIFECYCLE.md) records the current experiment. The
measurements below retain their original resource limits and provenance.

## Eight-CPU lifecycle, 2026-09-12

The main HPC change `1a7e920` was compared with the already-tested Rust
`a8cf33c` on the three full lncRNA BAMs. These are native/native comparisons,
with identical scientific settings and annotation mode on both sides. All 36
builds completed before the 30-minute soft boundary without OOM. Each metric
below is independently reduced to the median of three interleaved paired runs.
Each run has eight affinity CPUs, an eight-CPU quota, 16 GiB memory, no additional
swap and a one-hour hard kill. Annotation preparation is included; staging and
container startup are excluded. OS cache/storage are shared and not flushed.

| Mode / sample | Baseline s | HPC candidate s | Wall reduction | CPU s before / after | Mean cores before / after | RSS MiB before / after |
| --- | ---: | ---: | ---: | --- | --- | --- |
| locus / S025 | 271.966 | 203.866 | 25.0% | 815.1 / 655.8 | 2.989 / 3.267 | 514.0 / 825.1 |
| locus / S026 | 255.308 | 186.754 | 26.9% | 739.0 / 570.5 | 2.893 / 3.113 | 510.8 / 815.1 |
| locus / S027 | 258.735 | 189.654 | 26.7% | 765.0 / 615.6 | 2.949 / 3.226 | 511.5 / 818.6 |
| spladder / S025 | 531.247 | 433.737 | 18.4% | 2887.2 / 2568.9 | 5.435 / 5.930 | 1326.6 / 4418.4 |
| spladder / S026 | 487.138 | 400.628 | 17.8% | 2521.2 / 2244.2 | 5.206 / 5.728 | 1291.4 / 4498.1 |
| spladder / S027 | 517.366 | 405.628 | 21.6% | 2615.4 / 2368.7 | 5.237 / 5.829 | 1365.7 / 4535.7 |

The whole-lifecycle mean-core increase accompanies a 19.5–22.8% reduction in
child CPU seconds for locus mode and 9.4–11.0% for default compatibility mode.
Eight-CPU occupancy reaches 38.9–40.8% and 71.6–74.1%, respectively. Serial
annotation and graph-cache object processing still limit overall parallelism;
the HDF5 memory driver reduces small I/O operations but does not parallelize
the library. Default-mode RSS increases substantially and remains within the
16 GiB budget; this is a measured tradeoff, not an across-the-board memory saving.
Cgroup charged peaks, including file cache, are retained in the raw reports.

Six first-round full comparisons passed: 48 HDF5 files, 6,854,364 datasets and
36 scientific text files. Native float bits, object inventories, types, shapes,
attributes, links, compression and values match. All 24 later repeated outputs
also pass six non-graph HDF5 files and six texts each against their own first
repetition; their two private graph caches are not repeated in this later gate.
The parent's full regression suite passed, including Python fixtures and failure
paths. Full Python output parity on these three real BAMs remains unverified:
the previous Python runs timed out at one hour.

Median kernel output-block accounting decreases by 1.9–2.9%. Live ten-second
I/O samples also show fewer logical reads and small write calls. For S026 locus,
the last samples contain 4,411,937 versus 28,372 write calls, taken at 260.0 of
269.4 seconds and 190.0 of 194.0 seconds, respectively. These omit each run's
remaining tail and are not exact terminal counters. Physical reads are largely
served by shared OS cache, so no cold-disk read reduction is established.

The adopted [conditional-store investigation](ZERO_DEPTH_STORES.md) follows the
RSS increase; its immutable binary and additional measurements are kept distinct
from the completed 36-run experiment above. Its all-sample and repeated paired
gates passed, as recorded below.

Evidence: [acceptance report](../runs/hpc-lifecycle-20260912/acceptance/README.md),
[resource/input audit](../runs/hpc-lifecycle-20260912/acceptance/audit.json),
[three-repeat summary](../runs/hpc-lifecycle-20260912/acceptance/parent-summary.json),
[sampled I/O](../runs/hpc-lifecycle-20260912/acceptance/sampled-io.json).

## Adopted final candidate and worker scaling

The final scientific code is `53f8c71`, binary SHA256
`0b88334e1c8d4872dd0112ff1e93894c2114666e13f88d225fcf4fe2a007944d`.
It adds a store only when the prefix-sum destination needs to change, preserving
the parent's exact result. In three simultaneous S026/default parent-versus-final
pairs, median wall time decreases 3.8%, CPU time 4.2%, and peak RSS 30.2%
(4,310.8 to 3,008.9 MiB); mean occupied cores remain about 5.75. CPU-set placement
alternates across rounds. This local component measurement is separate from the
36-run parent comparison and is not multiplied into an inferred global speedup.

Final binary, one fresh complete build per sample/mode, including preparation:

| Mode / sample | Wall s | CPU s | Mean cores | Peak RSS MiB | Cgroup peak MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| locus / S025 | 188.640 | 627.1 | 3.325 | 812.5 | 2003.2 |
| locus / S026 | 173.602 | 548.6 | 3.161 | 807.9 | 2007.4 |
| locus / S027 | 183.417 | 588.8 | 3.211 | 815.2 | 2014.0 |
| spladder / S025 | 429.722 | 2533.9 | 5.897 | 3711.1 | 3945.8 |
| spladder / S026 | 381.469 | 2169.6 | 5.688 | 3969.4 | 4204.0 |
| spladder / S027 | 400.005 | 2309.6 | 5.775 | 3966.7 | 4201.5 |

All six final output comparisons pass against the parent, including both private
graph caches, all HDF5 data/structure and scientific texts with exact native float
bits. The two stages together comprise twelve full directory comparisons. All
54 measured acceptance builds completed before the soft limit without OOM;
resource/annotation audits and 36 later non-graph repeat comparisons also pass.
Full Python parity on the three real BAMs is still unverified after the earlier
timeouts. Fixture Python comparisons and the final 39-case CLI/reuse/one-four-eight
worker raw-bit gate pass. No additional dependency or schema/version change is made.

S026 worker scaling, one observation at each count, same eight-CPU/16-GiB quota
and fixed affinity 56–63; all associated non-graph outputs match exactly:

| Annotation mode | Four-worker s | Eight-worker s | Speedup | Mean cores four / eight | RSS MiB four / eight |
| --- | ---: | ---: | ---: | --- | --- |
| locus | 237.457 | 185.511 | 1.28× | 2.376 / 3.111 | 681.2 / 815.4 |
| spladder | 637.453 | 385.454 | 1.65× | 3.441 / 5.692 | 2387.5 / 4531.8 |

The serial portions of the shorter locus lifecycle limit its scaling gain.
Four-worker runs reserve the same eight-CPU envelope for this comparison; their
mean-core figures are not percentages of a four-CPU quota. The observations
share OS cache/storage and do not establish 64-core scaling or cold-disk behavior.
All raw repeats, exclusions, commands and source hashes remain in the
[acceptance record](../runs/hpc-lifecycle-20260912/acceptance/README.md).

## Earlier four-CPU benchmark methodology

The earlier acceptance workload used four logical CPUs and a hard 8 GiB cgroup memory
limit, with swap disabled. The local Docker image is pinned by digest; it runs
the same reference environment and native libraries used for correctness checks.
Every measured command records its actual affinity and cgroup limit. Up to 64
threads are accepted by the CLI; the performance baseline at that revision was four.

At the recorded source revision, `scripts/run_lifecycle_benchmarks.py` ran
three interleaved repetitions of:

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

For the historical four-CPU reproduction, check out the recorded source
revision first. The current runner enforces eight CPUs / 16 GiB. Use a new
output directory:

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

## Timeout root cause: annotation scope and skewed work

A subsequent bounded S026 investigation found 714 of the 31,130 GTF gene IDs
on more than one chromosome/strand, including 641 on multiple contigs. The
pinned Python parser keys genes only by gene_id and combines these exon
coordinates under the first chromosome/strand. Rust preserves that behavior.
For example, ZNF84 becomes a 133.60 Mb interval on chromosome 12 after exons
from Un_gl000223 are added. These mixed-locus models account for 85.8% of the
summed gene spans. Historical Python graphs contain the same abnormal bounds.

Instrumentation on branch `perf/p0-timeout-root-cause` (code commit `7bf2d6c`)
measured 1,104,073,305 alignment record returns during the serial initial
intron stage (382.306 seconds). The cassette stage then took 405.006 seconds;
its last worker continued querying for 346.585 seconds after the other three
workers' final evidence calls. Costly chromosome 6 genes were concentrated in
a large Rayon work item. Initial and cassette-tail profiles attributed 59.08%
and 43.18% of sampled user cycles to BAM decompression. The diagnostic build
also timed out at 1,200.189 seconds, during graph quantification, without OOM.

As a causal counterfactual, a separate GTF copy split only conflicting gene IDs
by chromosome and strand. All 1,490,438 data lines, coordinates, transcripts
and the full BAM were retained. Using the same instrumented binary and
4-CPU / 8-GiB constraints, this build completed in 351.802 seconds with
438.1 MiB peak process RSS; separate annotation preparation took 28.288 seconds.
Its initial intron stage returned 227,109,460 records in 96.244 seconds. The
instrumented binary was checked against production on chromosome 22:
23,658 HDF5 datasets and 6 text files matched, including raw float bits.

**The counterfactual changes gene identities and scientific outputs. It is not
an output-equivalent optimization or a successful rerun of all three original
samples.** These are single diagnostic measurements on a shared warm-cache
host, not a production speedup benchmark. Correcting annotation scope needs
an explicit compatibility decision. Output-preserving work should investigate
repeated decode work, task granularity and oversized coverage allocations.
After the interval problem is removed, graph-cache reads/writes are another
measured serial cost (156.855 seconds in the counterfactual).

All diagnostic jobs are stopped; the production binary is unchanged. The
[full causal report and raw evidence](/home/wubw/data/ruspladder/runs/lncrna-three-20260912/root-cause/README.md)
include the GTF audit, source links, per-region traces, counterfactual rules,
resource records, exclusions and limitations.

## Locus identity and lifecycle optimization, 2026-09-12

Implemented on `perf/p0-locus-lifecycle` (scientific code `bca30f7`, regression
integration `d608b26`). This follows the bounded root-cause study above. The
[annotation proposal](ANNOTATION_LOCUS_PROPOSAL.md) records the primary sources
and compatibility decision; the [complete experiment report](../runs/locus-lifecycle-20260912/README.md)
contains commands, binary hashes, original logs, resource records and comparisons.

`--annotation-mode locus` on `prep`/`build` normalizes conflicting GTF IDs by
`(gene_id, seqname, strand)` and scopes transcript identities to their parent
placement. It writes an auditable `.locus.gtf` plus `.loci.tsv`, then uses the
existing annotation importer and scientific algorithms. All original rows,
coordinates and non-ID attributes are retained. Only conflicting IDs change;
there is no distance-based splitting or removal of alternate contigs. Generated
names are deterministic for the same annotation, independent of record order,
and avoid collisions with existing names. They are not cross-release stable IDs.
This mode currently accepts GTF; the existing GFF3 ID/Parent importer is unchanged.
The tuple separates coordinate frames; it does not infer distinct same-strand
loci on one contig or reconstruct cross-locus transcripts.

The default `spladder` mode preserves the earlier interpretation. Output folders
record their mode to prevent cross-mode reuse; annotation companions follow the
existing immutable-input cache contract. Correcting mixed-locus graphs can also
change adaptive filters and outputs outside the affected genes. Locus comparisons
therefore use the same canonical GTF on both sides. They do not establish parity
with the old mixed-coordinate graph.

The output-preserving performance changes reuse the previously tested parallel
initial intron extraction, bound Rayon jobs to at most 64 genes, move the first
BAM coverage vector instead of allocating an additional full-span zero vector,
and retain one graph between counting, event collection and reporting. Counting
can temporarily update gene bounds; the retained graph restores those bounds
before later stages, matching the former disk-read behavior. These changes add
no dependency, HDF5 schema, read filter or event algorithm. Initial queries still
retain each gene's evidence and the upstream ordering/missing-contig behavior.
The aggregate benchmark does not attribute a separate speedup to each change.

All three full BAMs use confidence 3, read length 150, `merge_graphs`,
`exon_skip,mult_exon_skip` and text outputs. Each build has a hard four-CPU quota,
four disjoint affinity CPUs, an 8 GiB memory limit and no additional swap. Fresh
input/output caches are used per run, including annotation preparation in the
measured build duration. At 30 minutes the supervisor records a soft timeout;
at one hour it sends SIGKILL without a grace period. The timeout smoke test
terminated a sleeping child at 2.002 seconds for a two-second hard limit.
Benchmark containers shared storage and the OS page cache and ran concurrently.
These are one observation per sample/condition, not replicated cold-cache trials.
Peak process RSS is distinct from cgroup memory including file cache; complete
reports contain both, CPU seconds and the effective limits.

Same canonical annotation: old Rust versus optimized locus mode:

| Sample | Old seconds | New seconds | Wall reduction | Mean cores old → new | Peak RSS MiB old → new | Complete output |
| --- | ---: | ---: | ---: | --- | --- | --- |
| S025 | 458.011 | 320.333 | 30.1% | 1.644 → 2.253 | 434.9 → 391.1 | exact |
| S026 | 426.533 | 297.339 | 30.3% | 1.588 → 2.164 | 455.5 → 384.1 | exact |
| S027 | 438.690 | 304.924 | 30.5% | 1.629 → 2.209 | 422.1 → 379.1 | exact |

Original annotation: old Rust versus optimized default compatibility mode:

| Sample | Old seconds | New seconds | Wall reduction | Mean cores old → new | Peak RSS MiB old → new | Complete output |
| --- | ---: | ---: | ---: | --- | --- | --- |
| S025 | 1630.492 | 817.638 | 49.9% | 1.780 → 3.308 | 2485.3 → 939.9 | exact |
| S026 | 1423.627 | 708.765 | 50.2% | 1.790 → 3.262 | 2150.1 → 865.1 | exact |
| S027 | 1493.363 | 756.165 | 49.4% | 1.778 → 3.263 | 2651.8 → 881.1 | exact |

Both native comparisons passed for every sample: 8 HDF5 and 6 text files per
comparison, with 1,164,030 datasets per locus sample and 1,120,758 per default
sample. All twelve native builds completed before the 30-minute soft boundary
without OOM. The optimized default runs used 81.5–82.7% of their four-CPU
quota on average; optimized locus runs used 54.1–56.3% and finished sooner.
These percentages cover the entire build, including serial annotation and output.


The whole annotation audit retained 1,490,438 rows and changed IDs on 49,969
rows. There were 714 conflicting gene IDs: 31,130 original IDs become 32,332
placements. This file's 81,407 transcript IDs require no disambiguation.
The source SHA256 is `a9250d41a6601d40546cd0fa825c772c4591ca52475b1bb015e22495bb1258c1`;
all three canonical annotations have SHA256
`dcdcf934254a8483e7a5039ca302de76f419e08190fba556fbdfb2671765f3eb`.
The independent complete annotation-model comparison matched all 32,332 genes
against Python, including metadata, transcripts, exons, introns and graphs.

Validation passed formatting, release clippy with warnings denied, three unit
cases, 13 annotation/filter comparisons, 26 build CLI configurations (one/four
threads, BAM/CRAM, all six event types, all supported merge strategies, validation
and cache reuse), chunked/per-sample workflows, sparse build/CLI and prep CLI.
A two-BAM end-to-end fixture reuses upstream reads and changes only annotation
IDs to induce both gene and transcript placement conflicts. Its full Python
comparison passed, with maximum absolute float difference 5.684e-14 under the
project's existing atol=1e-10/rtol=1e-8. Cache reuse and mode isolation also passed.
Native/native comparisons require exact float bits, HDF5 objects, shape, type,
attributes, links, compression and data, plus decompressed text equality.

Full Python comparisons on the same canonical annotation:

| Sample | Wall seconds | Build outcome | Full output comparison | Last logged progress |
| --- | ---: | --- | --- | --- |
| S025 | 3600.060 | hard timeout | not completed | 9800(32332) genes done (found 3509 new retentions in 80884 tested eligible introns, 4.3%) |
| S026 | 3600.067 | hard timeout | not completed | 18800(32332) genes done (found 6876 new retentions in 156406 tested eligible introns, 4.4%) |
| S027 | 3600.064 | hard timeout | not completed | 9800(32332) genes done (found 2504 new retentions in 80885 tested eligible introns, 3.1%) |

Python builds keep their partial files and terminal records. The whole-annotation
and small end-to-end comparisons described above remain valid; a timed-out full
Python build does not establish full-sample Python output parity. The recorded
hard-timeout duration includes signal delivery and process reaping; no additional
processing grace period or resumed run is used. The last log line is buffered
progress, not proof of the exact instruction executing when the process stopped.

To use the normalized mode, keep the BAM and annotation at their existing paths
and choose a new result directory:

```sh
OPENBLAS_NUM_THREADS=1 /home/wubw/data/ruspladder/target/release/ruspladder build \
  --bams sample.bam --annotation RefSeq_hg19.gtf --outdir results-locus \
  --annotation-mode locus --parallel 4 --readlen 150 --confidence 3 \
  --merge-strat merge_graphs --event-types exon_skip,mult_exon_skip --output-txt
```

`--parallel 4` selects worker count; enforce memory/CPU limits and the deadline
with the recorded container supervisor when reproducing the benchmark. Removing
`--annotation-mode locus` selects default compatibility behavior and needs a
separate result directory. The original GTF and previous results are preserved.

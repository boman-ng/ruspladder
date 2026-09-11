# Migration contract and evidence

## Fixed scope

- Reference: SplAdder v3.1.1, commit 65ceec839b9ff0cf96703c1605ee43667662f410.
- Complete prep/build/test algorithms and non-visual CLI behavior, all six events,
  confidence 0–3, four merge strategies, all qmode and chunked-merge paths.
- Independent Rust production runtime; HTSlib/HDF5 native libraries allowed.
- Public HDF5, TXT/TSV, GFF3, BED, structured, TCGA, and ICGC outputs compatible.
- Rust-owned HDF5 graph caches; no Python pickle compatibility.
- Exact discrete values, IDs/order, missingness and statistical decisions.
  Quantification atol=1e-10/rtol=1e-8; statistics atol=1e-8/rtol=1e-6;
  also run upstream's six-decimal statistics comparison. Do not relax tolerances
  to pass failures. Rust 1/4-thread logical outputs must be bit-identical.
- Benchmark complete fresh and reused-cache workflows, 4 CPU / 8 GiB actual
  enforced limits, >=3 interleaved runs. Report wall/CPU time, full-life CPU
  utilization, RSS, job peak memory, I/O and stage timings. Require measured
  runtime and peak-memory improvement on real subsets; no invented speedup target.
- P0 candidates get separate perf/p0-* worktrees. Minimal reproduction, research,
  alternatives, parity and measured benefit precede adoption; user authorized
  autonomous decisions. Never change scientific thresholds for performance.
- Work and caches stay in ~/data/ruspladder. Do not release/publish or change
  established project versions without authorization.
- User's execution ceiling: at most 64 CPU cores.
  Build defaults to 16 jobs; scientific benchmark remains 4 threads / 8 GiB.

## Implementation checkpoints (not a reduced completion definition)

- [x] Reproducible reference environment, complete CLI inventory, fixture manifests.
- [x] Annotation parsing/filtering, splice/segment graph primitives and caches.
- [ ] BAM/CRAM filtering, coverage, sparse prep, reference lookup.
- [ ] Graph augmentation, re-inference, pruning, merge strategies.
- [x] Six event detectors, collection, sorting, curation, feature verification.
- [x] Segment/junction counts, gene expression, PSI, public build result formats (CLI wiring pending).
- [ ] Normalization, NB/Gamma GLM, dispersion fitting/shrinkage, LRT and correction.
- [ ] Full prep/build/test orchestration and cache reuse.
- [ ] Upstream fixtures and airway real-data differential verification.
- [ ] 1/4-thread determinism and 4 CPU / 8 GiB lifecycle benchmarks.
- [ ] P0 investigation, all required checks, final compatibility and performance reports.

## Sources and licensing

SplAdder paper and supplement: https://www.ong-home.my/papers/kahles16spladder.pdf
Upstream classes/init/reads/editgraph/detect/count/test and their called modules
are the behavior reference. Preserve code-level details even when unusual;
document upstream failures instead of calling skipped comparisons successful.

SplAdder's COPYRIGHT states BSD-3-Clause, but some individual modules carry
GPL-3.0-or-later notices. The combined migration is conservatively GPL-3.0-or-later,
with upstream BSD and file-specific notices retained in licenses/ and source
attribution. Third-party libraries retain their own licenses.

Reuse rust-htslib, Rayon, hdf5-metno, nalgebra and statrs where applicable.
Port the called statsmodels/SciPy numerical paths when required for parity;
do not replace the statistical method with another model.

## Current environment / continuation

Development worktree: /home/wubw/data/ruspladder/worktree (branch implement).
Reference checkout: ../upstream/spladder; reference Python: ../envs/reference/bin/python.
Existing Git entrypoint: /home/wubw/ruspladder (empty initial master).
Slurm partition soma enforces 4 logical CPUs and 8,589,934,592 bytes memory,
verified in runs/resource-probe.json. Slurm accounting is disabled, so benchmarks
must gather their own CPU, memory and I/O measurements.
No full migration or benchmark success has been established.

P0 candidate identified in reads.summarize_chr: dense (3, chromosome_length)
uint32 coverage allocation before sparse conversion, ~3 GB per 250 Mb contig.
The separate perf/p0-chromosome-coverage branch contains its reproduction.
One worker reached 1,998,336 KiB RSS with only two short reads per contig;
four workers on four 300 Mb contigs caused two OOM kills under the enforced
8 GiB limit, followed by a hung upstream pool until the job timeout.
The bounded-window native prototype is committed separately at 9955eb2;
96 sparse summaries / 1,008 HDF5 arrays pass exact parity. It remains unmerged
pending performance evidence. Slurmctld is currently DOWN; the benchmark failed
before starting with an allocation/connect error. Ordinary local timings do not
substitute for the required enforced 4 CPU / 8 GiB run.

A second candidate, perf/p0-hdf5-chunks at 7e2f370, reproduced 64 MiB automatic
chunks on an unlimited sample axis. Its release 8-by-2 matrix microbenchmark
used four-CPU affinity and an 8 GiB address-space cap (not a cgroup memory limit).
Three default runs took 7.69–7.71 s; bounded chunks took 0.0088–0.0091 s with exact
logical output. The explicit chunk layout is adopted for count collection and
event batches; public file comparisons pass. This measures a local pathology,
not overall pipeline speed. See ../p0-hdf5-chunks/RESEARCH.md and
../runs/p0-chunks/report.json for the method, alternatives and evidence.

## Verified implementation boundaries

Reference NumPy 2.2.6 / SciPy 1.13.1 / statsmodels 0.14.4 is locked;
all six upstream tests pass without source changes. The latest statsmodels
removed an imported upstream symbol; NumPy 1.x cannot read the committed
NumPy 2 pickle fixtures. See runs/upstream-locked.log for the accepted baseline.

Current differential checks against that reference:

| Boundary | Comparisons | Evidence under ../runs |
| --- | ---: | --- |
| Splice/segment graph primitives | 614 | compare_graphs.py |
| Annotation/filter/cache roundtrip (including airway GTF) | 14 | annotation-parity-real.log |
| BAM/CRAM coverage and filtering | 106 | reads-parity.log |
| Six raw event detectors | 208 | detector-parity/report.json |
| Collection, curation, coordinates and IDs, 1/4 threads | 20 | collection-parity/report.json |
| Gene labels, short exons, duplicate merge | 832 | edit-parity/report.json |
| Cassette insertion and intron retention from matched coverage | 1278 | augmentation-parity/report.json |
| Intron edge insertion from matched coverage | 840 | intron-parity/report.json |
| Intron ambiguity / FASTA consensus filtering | 100 | intron-filter-parity/report.json |
| Complete direct-BAM graph generation, including airway | 45 | build-graph-parity/report.json |
| Sample/chunk merge, support filtering, cache roundtrip | 42 | merge-parity/report.json |
| Graph segment/junction counting, 1/4 threads | 24 | count-parity/report.json |
| Six event feature vectors, flags, PSI | 1330 | verify-parity/report.json |
| Public graph counts and gene expression | 60 files / 560 arrays | check-20260911T205705/count-io/report.json |
| Geometric-mean and total-count normalization | 160 | check-20260911T205705/count-io/report.json |
| Single-quantification collection, first-file types and unlimited axes | 12 / 108 arrays | check-20260911T205705/count-collect/report.json |
| Six event analysis HDF5, 1/4 threads and 128-event batches | 130 / 1580 arrays | check-20260911T205705/analysis/report.json |
| TXT, structured, BED, GFF3/GTF, TCGA, ICGC and gzip text | 990 files | check-20260911T205705/outputs/report.json |
| Isoform counts for testing, sample and PSI ordering | 324 | quantify-parity/report.json |
| Six working multiple-testing corrections and decisions | 174 | correction-parity/report.json |
| NB/log and Gamma/identity GLM, including exact iteration counts | 263 valid / 57 reference errors | check-20260911T213736/glm/report.json |
| Bounded scalar optimizer, statuses and evaluation counts | 168 | check-20260911T213736/optimize/report.json |
| NB likelihood, Cox–Reid, shrinkage, trigamma and chi-square kernels | 45 | check-20260911T213736/likelihood/report.json |
| Raw/trend/shrunken dispersion, LRT and final isoform selection | 12 × 60 rows, 1/4 threads | check-20260911T213736/statistics/report.json |

The augmentation comparison injects identical coverage at the upstream I/O
boundary; it verifies graph algorithms, not the full build workflow. Raw detector
comparisons separately record 24 upstream empty-graph exceptions; those do not
count as passing comparisons. Cargo clippy --all-targets -- -D warnings passes.

Graph generation includes confidence 0–3, individual augmentation switches,
short-exon removal, strand and consensus settings, and multiple direct BAMs.
The airway BAMs contain no NM tags; both implementations explicitly use the
upstream --ignore-mismatches option (recorded in the graph-generation report).
Event verification flags are exact. Feature vectors and PSI use the agreed
quantification tolerance; observed maximum absolute difference is 1.78e-15.
Sample merge streams one sample at a time and retains sparse edge support.
The reference's metadata ownership, terminal behavior and validation side
effects are preserved. Advanced paths that fail upstream need separate
failure-contract documentation; no invented successful results are substituted.

Public array comparisons include shapes, dtypes, links, compression and values;
collection also checks unlimited axes. Chunk geometry is intentionally changed
for bounded access. Text comparisons are exact bytes after gzip decompression.
Upstream's unimplemented structured multi-exon output and BED multi-exon/mutex
output are reproduced and documented in COMPATIBILITY.md, not invented.
The complete component regression is recorded in check-20260911T213736 and
check-statistics-current.log. Its earlier quantification fixture enumeration
missed a synthetic multi-exon fixture; that guard was corrected before this run.
The new run also confirms exact event feature/PSI values after enabling exact
JSON float round trips. See NUMERICS.md for the native kernel decisions and
the retained failed numerical comparisons that led to them.

The CLI currently implements annotation prep only. Sparse-input integration,
test input preparation/final text output, complete CLI orchestration/cache reuse
and lifecycle benchmarking remain open. Passing the statistical numerical
sequence does not constitute a completed differential test workflow.
No end-to-end replacement or full-lifecycle runtime/memory improvement is claimed yet.

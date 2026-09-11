# P0 candidate: BAM decompression dominates real-data CPU

A user-space `perf record -e cpu-clock:u -F 199` sample of the complete native
8-sample airway build (before this change) captured 265 samples, no lost samples:
56.98% in `inflate_fast`, 5.28% in `crc32_z`, and a further 3.40% in inflate and
its table setup. Evidence: ../runs/profile-airway/{perf.data,report.txt}.
This is hotspot evidence, not a constrained benchmark or claimed speedup.

The [HTSlib paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7931820/) identifies
libdeflate as a contributor to BAM performance. HTSlib's own
[compression benchmarks](https://www.htslib.org/benchmarks/zlib.html) recommend
it. These published speedups are not ruspladder measurements. The installed
rust-htslib 1.0.1 already exposes the `libdeflate` feature, forwarded to
hts-sys 2.2.1. Its build script defines HAVE_LIBDEFLATE and links the library.
No BAM filtering, region selection or scientific algorithm needs to change.

Candidate: enable that existing feature. Baseline and candidate use identical
source at 1a44852 and identical release settings; copy both finished binaries
before measurements. Compare full airway builds in three interleaved, fresh
4-CPU/8-GiB containers per variant, including all graph and output I/O. Compare
all native HDF5 datasets and decompressed scientific texts exactly. Separately
rerun the 106-case upstream BAM/CRAM reader comparison with the candidate.

Alternative: add HTSlib decompression worker pools to every gene reader. This
would interact with existing Rayon workers and multiply per-handle threads and
buffers. Do not add that orchestration before measuring the supported library
switch. Replacing zlib globally would affect HDF5 too; the HTSlib feature is the
smaller change for the observed BAM hotspot. No custom codec is warranted.

Adoption is pending parity, provenance/license review and measured lifecycle
benefit without peak-memory regression. The final source-vs-Rust lifecycle
benchmark remains separate from this isolated candidate comparison.

The lockfile resolves libdeflate-sys 1.26.0, containing upstream libdeflate 1.26.
Bindings are Apache-2.0; the C library is MIT (Eric Biggers / Google LLC).
Both notices are retained in licenses/. The Cargo registry checksum pins the
source archive. The dependency uses its existing static build; production gains
no service, Python dependency, custom codec or runtime library search path.

A second candidate combines libdeflate with the existing Rayon gene loop pattern
already used for graph counting. Cassette insertion and retention insertion each
modify only their own gene after the shared feasibility/filter step. Parallelize
those two loops with worker-local indexed readers and integer count reduction;
keep samples, mutable filter transitions and neighbor-gene intron insertion
sequential. Compare this combination as a third variant in a rotating three-run
order. This directly tests whether independent gene work improves full-life CPU
occupancy beyond the codec switch; do not assume that more threads are faster.

## Constrained results and decision

Three interleaved fresh airway builds per variant, all in separate Docker
4-CPU/8-GiB cgroups. Median results:

| Variant | Wall s | CPU s | Mean cores | Cgroup peak bytes |
| --- | ---: | ---: | ---: | ---: |
| Baseline | 1.51368 | 1.39215 | 0.93183 | 27,267,072 |
| libdeflate | 1.17217 | 0.91055 | 0.77309 | 28,041,216 |
| libdeflate + independent gene stages | 0.93612 | 0.97958 | 1.04643 | 30,269,440 |

The combined candidate reduces median wall time by 38.2% and CPU seconds by
29.6%, while increasing median mean-core occupancy by 12.3%. Peak memory rises
by about 2.9 MiB (11.0% of this very small native baseline). The initial assumption
of no peak-memory regression is not supported; report this tradeoff explicitly.
All runs remain below 31.1 MB, far below the enforced 8 GiB ceiling. Adopt the
combined candidate for the measured time/CPU benefit; require the final complete
source-vs-native benchmark to establish overall time and memory improvement.

21,948 datasets across internal graph/event caches and public HDF5, and 48
scientific text files, compare exactly between the baseline and both candidates.
The candidate reader passes 106 upstream BAM/CRAM cases. Parallel graph generation
passes 45 upstream cases and the 201-intron sequential filter-state regression.
An additional positive fixture inserts eight new cassette exons at each of four
confidence levels: 1/4-thread outputs agree exactly with upstream (64 insertions).
Reports: ../runs/p0-decompression/{benchmark/report.json,read-parity.log,
graph-parity/report.json,sequence-parity/report.json,cassette-parity/report.json}.
Raw per-run RSS, cgroup limits/peaks, I/O and commands are retained in the benchmark.
No published HTSlib speedup is substituted for these measurements.

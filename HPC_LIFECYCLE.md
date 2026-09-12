# Eight-CPU lifecycle investigation

Work is isolated on `perf/p0-hpc-lifecycle`, based on `a8cf33c`. The acceptance
configuration is now 8 CPUs / 16 GiB, no additional swap, with a 1,800-second
soft notification and 3,600-second hard SIGKILL boundary. Each run uses fresh
annotation/output caches and the full S025/S026/S027 BAM. Both annotation modes
are compared only against the same mode. Cache files and public schemas stay
unchanged. Existing artifacts and the frozen baseline are preserved.

## Sources and decisions

- Williams, Waterman and Patterson, *Roofline* (CACM 2009,
  [paper](https://aiichironakano.github.io/cs596/Williams-Roofline-CACM09.pdf),
  DOI 10.1145/1498765.1498785): distinguish work, memory traffic and attainable
  parallel performance. We report wall time, CPU seconds, mean occupied cores,
  process RSS and cgroup peak separately. This irregular integer/I/O workload
  is not a measured floating-point Roofline plot, and occupancy alone is not
  evidence of efficient execution.
- Bonfield et al., *HTSlib* (GigaScience 2021, DOI 10.1093/gigascience/giab007),
  [upstream library](https://github.com/samtools/htslib) and
  [Tabix documentation](https://www.htslib.org/doc/tabix.html): the existing
  BGZF cache avoids rereading/reinflating compressed blocks across overlapping
  queries. The locked HTSlib source has BGZF_CACHE enabled. Its public
  `hts_set_cache_size` API was tested in the pilot; the final candidate does not
  enable this additional cache. CRAM is unchanged. No codec,
  alignment cache format or new dependency is introduced.
- Pedersen and Quinlan, *mosdepth* (Bioinformatics 2018, DOI
  10.1093/bioinformatics/btx699,
  [author-hosted paper](https://ucgd.genetics.utah.edu/wp-content/uploads/2017/11/btx699.pdf)):
  start/end differences and one prefix sum replace per-base updates for each
  CIGAR block. This is already used by the sparse preparation code. Direct
  queries now reuse the output buffer for differences, with modular uint64
  arithmetic for negative deltas. SplAdder's deletion coverage, overlapping
  mates, filters, strands and junction/read counts are preserved. Mosdepth's
  biological filtering semantics are not substituted.
- Rayon 1.12 [ParallelBridge](https://docs.rs/rayon/1.12.0/rayon/iter/trait.ParallelBridge.html)
  and the locked `iter/par_bridge.rs`: dynamically pull independent genes while
  retaining the initializer's reader across work. Existing `with_max_len(64)`
  produced many short-lived readers and repeated index loading. Bridge's
  synchronized `next()` may cost throughput on cheap tasks; the gene kernels
  do substantial independent work. Count slots are paired with genes before
  scheduling because bridge completion order is unspecified. Mutable filter
  transitions and neighbor-dependent graph updates remain ordered.
- HDF Group [memory file driver](https://portal.hdfgroup.org/documentation/hdf5/latest/_h5_f__u_g.html)
  and [metadata caching](https://portal.hdfgroup.org/documentation/hdf5/latest/_t_n_m_d_c.html):
  test the supported core driver for existing graph files, with backing-store
  enabled for atomic writes and disabled for read-only loads. This trades a
  whole graph-file image in memory for fewer small file operations. It does
  not change HDF5 object names, types, data, compression or cache contracts.
  The driver does not make the HDF5 library parallel. Its object-creation CPU
  cost remains. Increasing cache sizes alone is not presumed beneficial.

The measured P0 mechanisms are repetitive BGZF decode/index work and serial
cache I/O. The selected candidate retains the single most recent graph across
build, merge, validation and quantification, extending the existing retention
through event reporting. It consumes a matching in-memory graph directly and
reads other graphs in the established order. There is no unbounded collection
of sample graphs. The BGZF cache pilot is not adopted: with the same reuse/core-driver changes,
64 MiB caching improved locus wall time by only 6 seconds (174.8 vs 180.9),
but made default-mode wall time worse (397.3 vs 382.1). Both CPU times and
RSS increased. The final candidate retains readers/indexes without enabling
an additional BGZF cache. This is a single paired pilot, so small wall-time
differences are not treated as statistically established benefits.

## Pilot evidence (not final acceptance)

The first S026 locus baseline at eight CPUs completed in 239.55 seconds with
2.84 mean cores. Two 20/30-second user CPU profiles sampled decode + CRC at
43.8–48.0%, and index decoding around 5%. They are interval profiles rather
than whole-lifecycle attribution; this profiled baseline will be excluded
from final performance estimates. The separately frozen read-cache/difference
candidate completed 106 exact Python BAM/CRAM filtering/region cases.

A full 32,332-gene cache roundtrip took 70.2 seconds with the default driver
and 62.1 seconds with the core driver, with exact in-memory graph equality.
RSS increased from 397 to 681 MiB. These single pilot runs share filesystem
cache; the core pilot overlapped one baseline CPU set briefly, so final
adoption requires full lifecycle measurements under disjoint resource sets.

Raw pilots, failed launcher attempts, immutable binaries, source diffs, hashes,
commands, profiles and resource samples are under
`../runs/hpc-lifecycle-20260912/`. The corrected parent candidate passed the
acceptance gate recorded below. No publication or version change is involved.

A further I/O check found that an untracked core driver increased kernel
accounted write blocks despite identical final file sizes. The pinned
`H5FD__core_flush` implementation rewrites its whole image on a dirty flush
unless write tracking is enabled. The selected candidate reuses
[H5Pset_core_write_tracking](https://portal.hdfgroup.org/documentation/hdf5/latest/group___f_a_p_l.html)
with 64 KiB aggregation pages and closes the root group before the file.
No file schema or data layout changes are involved. This is measured again
in a cache roundtrip and the complete paired benchmark, not assumed from the
API documentation.

## Compatibility gate and measurement exclusions

The first candidate failed an existing non-final chunk workflow assertion: the
core driver omitted the OS `ENOENT` reason. Graph-cache opening now translates
failed opens to the corresponding OS error when present, retaining the HDF5
error for readable but invalid files. Successful opens do not perform extra
file operations. The unchanged workflow test passed on the corrected binary
(two workflows, eight graphs, 203 arrays and five required upstream failures).
The first paired-run controller was stopped and all twelve launched builds
were allowed to finish. They are retained as predecessor evidence and excluded
from acceptance; the corrected binary has its own fresh three-round experiment
in `../runs/hpc-lifecycle-20260912/acceptance/`.

An experimental attempt to collect `/proc/PID/io` after exit using WNOWAIT was
rejected: the unprivileged container denies access to those counters after a
process exits. It was tested only with a two-second timeout smoke command;
that child was killed on time, while its experimental measuring helper needed
to be stopped. The existing proven deadline runner was restored before the
acceptance children started. All six initial acceptance measuring processes
were inspected live and had only the standard main thread. A second smoke
command with the restored runner exited at 2.0018 seconds with both soft/hard
flags and the expected SIGKILL return. No privilege escalation or instrumentation
was added to the Rust executable. I/O evidence is getrusage block accounting
and live ten-second samples, not unavailable exact terminal /proc counters.

## Completed parent acceptance

The immutable `candidate-final` matches scientific source `1a7e920` and SHA256
`ce2c4b547e540a51497007a4fac6d2c6f556431aad474bc21115a10acdffed08`.
All 36 actual containers passed the CPU/memory/input/output audit and completed
before the soft boundary. Three-run median wall times fall 25.0–26.9% in locus
mode and 17.8–21.6% in compatibility mode. CPU seconds decrease while mean
occupied cores increase in both modes. Peak RSS rises from about 0.5 to 0.8 GiB
in locus mode and about 1.3 to 4.3–4.4 GiB in compatibility mode.

Six unchanged full native comparisons passed, covering 48 HDF5 files,
6,854,364 datasets and 36 texts with exact native float bits. Every later
repetition also passed the established comparator on the six non-graph HDF5
files and six texts. The entire regression suite completed successfully at
`../runs/check-20260912T134231`, with its log and source hashes preserved in
`acceptance/validation/`. The existing CLI suite was subsequently corrected to
apply its raw-bit one-worker comparison to eight as well as four workers; the
frozen conditional-store candidate passed all 39 corrected CLI cases.

[PERFORMANCE.md](PERFORMANCE.md) reports the complete parent table, resource
tradeoffs and I/O interpretation. The [raw acceptance record](../runs/hpc-lifecycle-20260912/acceptance/README.md)
links immutable inputs, commands, limits, source hashes and exclusions. The
[zero-store branch](ZERO_DEPTH_STORES.md) investigates the measured RSS increase
without changing the parent's completed experiment. Its all-sample and repeated
simultaneous-pair gates passed; the final scientific source is `53f8c71`. The
conditional store reduced S026/default median wall time by 3.8%, CPU time by
4.2% and RSS by 30.2% in three simultaneous pairs, with exact output checks.
The twelve full native comparisons and all 54 measured builds passed their
respective output and resource gates. See the final tables in PERFORMANCE.md.

# P0 investigation: redundant stores in coverage prefix sums

Parent: `1a7e920` on `perf/p0-hpc-lifecycle`. The parent replaces repeated
per-read base increments with an in-place difference/prefix sum. Its full
S026 default-mode pilot uses about 4.3 GiB RSS, compared with about 1.3 GiB in
the earlier binary, while the locus mode remains below 1 GiB. The parent was
tested separately; its frozen candidate and three-repeat results are preserved.

Hypothesis: unconditional writes of zero prefix sums materialize otherwise
untouched pages in large, mostly uncovered gene intervals. The candidate only
stores the new prefix sum if it differs from the value already in that slot.
For every possible uint64 input state the resulting buffer is identical: when
the store is skipped its target already equals the intended value. This does
not alter read filtering, coverage arithmetic, regions, junctions or graphs.
It adds a branch, which could offset any memory-traffic benefit on dense data.
The optimization is tested rather than presumed beneficial.

The experiment uses existing Python read/filter comparisons, then the entire
S026 BAM and original GTF, eight CPUs / 16 GiB, fresh caches, 30-minute soft and
one-hour hard limits. The initial pilot runs on CPUs 56–63 while the parent's
paired experiment uses 0–47. The source diff, frozen binary/hash, command,
telemetry and outputs are preserved under
`../runs/hpc-lifecycle-20260912/acceptance/candidate-spladder-zero/`.
The completed adoption decision is recorded below; original pilot evidence is retained.

The kernel's [memory mapping documentation](https://cdn.kernel.org/doc/html/latest/admin-guide/mm/nommu-mmap.html#further-notes-on-no-mmu-mmap)
describes the relevant MMU case: anonymous zero-filled regions can initially
read from shared zero pages; writing a page materializes its physical storage.
This explains why an arithmetically redundant store can matter. Allocator reuse
and actual read coverage also affect RSS, so the full observed memory change
is not attributed solely to this mechanism without a page-level profile.

The first full S026 candidate completed in 391.559 seconds, using 2,213.615 CPU
seconds and 2,793.6 MiB peak RSS. The parent's first acceptance observation was
430.840 seconds, 2,467.373 CPU seconds and 4,136.0 MiB. These used different
CPU sets and times on the shared host; they justify a fixed-CPU follow-up rather
than establish a paired performance estimate. The 106 existing BAM/CRAM
region/filter comparisons passed exactly. The fixed-CPU parent/child pair is
reported below; all three samples in both annotation modes subsequently passed
the completed gate with immutable binaries.

Review found that the CLI suite exercised eight workers but its raw float-bit
thread comparison ran only for four workers. The existing comparison now runs
for both four and eight against one worker, without changing its criteria.
The frozen zero-store binary passed all 39 CLI scenarios, including cache reuse,
Python comparisons and one/four/eight-worker raw float bits. The parent's full
regression suite already passed; no unrelated algorithm or schema is changed here.

The same-CPU S026 follow-up (56–63, sequential parent then zero candidate) gave:

| Metric | Parent | Conditional stores |
| --- | ---: | ---: |
| Wall seconds | 409.722 | 385.454 |
| Child CPU seconds | 2208.327 | 2193.920 |
| Mean occupied cores | 5.390 | 5.692 |
| Peak RSS MiB | 4687.7 | 4531.8 |
| Cgroup peak MiB | 5436.7 | 5272.0 |

Wall time fell 5.9%, CPU time 0.7% and RSS 3.3% in this pair. The large RSS
reduction in the initial observation did not repeat in this sequential pair; task timing and overlapping
large intervals affect the peak. This supports a modest local benefit, not a
stable multi-gigabyte saving. The subsequent all-sample output/resource gate passed.

## All-sample and scaling observations

Scientific code `53f8c71`, frozen binary SHA256
`0b88334e1c8d4872dd0112ff1e93894c2114666e13f88d225fcf4fe2a007944d`,
completed all three full samples in both modes before the soft timeout, without
OOM. Each row is one fresh build, with preparation included. These additional
cohorts differ from the parent's paired schedule; no speedup is inferred by
comparing their absolute times with the parent's first run or median.

| Mode / sample | Wall s | CPU s | Mean cores | Peak RSS MiB | Cgroup peak MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| locus / S025 | 188.640 | 627.1 | 3.325 | 812.5 | 2003.2 |
| locus / S026 | 173.602 | 548.6 | 3.161 | 807.9 | 2007.4 |
| locus / S027 | 183.417 | 588.8 | 3.211 | 815.2 | 2014.0 |
| spladder / S025 | 429.722 | 2533.9 | 5.897 | 3711.1 | 3945.8 |
| spladder / S026 | 381.469 | 2169.6 | 5.688 | 3969.4 | 4204.0 |
| spladder / S027 | 400.005 | 2309.6 | 5.775 | 3966.7 | 4201.5 |

S026 worker-count observations use the same 8-CPU/16-GiB envelope and affinity
56–63. They are sequential single observations, not replicated scaling estimates.
All associated non-graph HDF5 and scientific text comparisons passed exactly.

| Annotation mode | Four-worker s | Eight-worker s | Speedup | Mean cores four / eight | RSS MiB four / eight |
| --- | ---: | ---: | ---: | --- | --- |
| locus | 237.457 | 185.511 | 1.28× | 2.376 / 3.111 | 681.2 / 815.4 |
| spladder | 637.453 | 385.454 | 1.65× | 3.441 / 5.692 | 2387.5 / 4531.8 |

The smaller locus scaling gain is consistent with substantial serial annotation
and graph-cache work in its shorter lifecycle. Occupancy is measured over the
whole build; CPU time and wall time are considered separately.

The final gate passed six full comparisons, including private graph caches,
and three rounds of simultaneous parent/zero S026 runs. Within each round the
two variants shared background conditions, alternating CPU-set assignment
(0–7 and 8–15). This addresses the earlier sequential pilot's changing background
workload. The paired runs began after the all-sample cohort finished; full
comparisons used 48–55.

## Decision: adopt the conditional store

All three simultaneous pairs completed successfully; the same-mode non-graph
outputs in all six repeated builds match exactly. Each metric below is reduced
independently to the median of three, not reconstructed from another median.

| Metric, independent median of three | Parent | Conditional stores |
| --- | ---: | ---: |
| Wall seconds | 402.018 | 386.720 |
| Child CPU seconds | 2320.802 | 2222.497 |
| Mean occupied cores | 5.751 | 5.754 |
| Peak RSS MiB | 4310.8 | 3008.9 |
| Cgroup peak MiB | 4977.9 | 3354.4 |

On this S026/default workload, median wall time decreases 3.8%, child CPU time
4.2%, and RSS 30.2%; median occupied cores stay about 5.75. All three individual
pairs improve wall time and RSS. The original large single-pilot RSS decrease
was variable; the repeated evidence supports the modest runtime gain and a
material memory benefit on this case. This is a local three-repeat measurement,
not a universal memory bound or a statistical population estimate.

The six final all-sample native comparisons passed 48 HDF5 files, 6,854,364
datasets and 36 texts, including raw float bits. Together with the parent's
comparisons, twelve full directory comparisons pass. All 54 acceptance builds
(parent experiment plus follow-ups) completed before the 30-minute soft limit
without OOM under the eight-CPU/16-GiB envelope. Two worker-scaling builds use
four workers inside that same quota; all other builds use eight. The source
binary, 39-case one/four/eight-worker CLI validation, 106 read comparisons,
format/release-clippy check, all limits and exact commands are retained in the
[acceptance report](../runs/hpc-lifecycle-20260912/acceptance/README.md).

The implementation remains the three-line conditional in `src/reads.rs`.
No new dependency, coverage representation, tuning option, schema or biological
filter is introduced. Python full-sample parity remains limited by the earlier
one-hour timeouts; this decision does not claim otherwise.

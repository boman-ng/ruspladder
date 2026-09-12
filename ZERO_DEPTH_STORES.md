# P0 investigation: redundant stores in coverage prefix sums

Parent: `1a7e920` on `perf/p0-hpc-lifecycle`. The parent replaces repeated
per-read base increments with an in-place difference/prefix sum. Its full
S026 default-mode pilot uses about 4.3 GiB RSS, compared with about 1.3 GiB in
the earlier binary, while the locus mode remains below 1 GiB. The parent is
being tested separately and its frozen candidate is unchanged.

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
Further adoption depends on complete runtime/memory evidence and output checks.

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
reported below; all three samples in both annotation modes remain pending,
with immutable binaries.

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
reduction in the initial observation did not repeat; task timing and overlapping
large intervals affect the peak. This supports a modest local benefit, not a
stable multi-gigabyte saving. The all-sample output/resource gate remains pending.

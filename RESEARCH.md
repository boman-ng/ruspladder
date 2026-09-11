# P0 candidate: chromosome coverage allocation

Reference: SplAdder v3.1.1 reads.py summarize_chr, lines 589 and 638.
It creates a three-row dense uint32 chromosome array and later reduces it into
an unstranded sum (NumPy widens the sum) before conversion to COO.

Hypothesis: chromosome-length allocation/reduction dominates memory even for
tiny sparse BAMs, and concurrently summarizing long contigs can exceed 8 GiB.
This is a hypothesis until recorded resource measurements confirm it.

Minimal reproducer: indexed BAM with four 250 Mb contigs, two short alignments
per contig; call the unmodified upstream summarizer in four workers under an
actual 4 CPU / 8 GiB Slurm cgroup. No production input or external data is used.

Candidate methods to compare after reproduction: bounded coordinate windows
and coverage difference accumulation, using the same CIGAR/filter semantics.
Do not introduce reduced resolution or scientific filtering changes.

## Reproduction results (2026-09-11)

- A real Slurm step confirmed affinity to four logical CPUs and an effective
  hard cgroup memory limit of 8,589,934,592 bytes (job 13932).
- One worker, four 250 Mb contigs: 8.84 seconds inside the parent; peak child
  RSS 1,998,336 KiB despite just 200 covered bases and two introns per contig.
  Evidence: ../runs/p0-coverage-single.json and .log.
- Four workers, four 300 Mb contigs: the input BAM is 294 bytes. Slurm reported
  two cgroup OOM kills (job 13934); Python's multiprocessing pool did not finish
  after its workers died, and the two-minute job limit terminated the step.
  Evidence: ../runs/p0-coverage-parallel.json and .log. This is an observed
  failure, not a measured successful runtime. Slurm accounting is disabled;
  the step log is the available authoritative termination evidence.

The memory and pool-failure hypotheses are now reproduced. Preserve output
semantics while eliminating dense chromosome-length storage. Do not reproduce
the worker-death hang in the Rust runtime.

## Candidate and decision

Reference methods: [mosdepth's implementation description](https://github.com/brentp/mosdepth#how-it-works)
and [its paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6030888/) explain start/end
difference accumulation followed by cumulative sums. Its full-chromosome array
retains chromosome-length memory cost, so that storage scheme alone does not
resolve this reproduction under concurrent work.

Candidate implemented here: indexed fixed-size windows with difference arrays;
reuse one HTSlib handle per worker, calculate batches on Rayon, and append their
COO output to compressed HDF5 in coordinate order. Each worker's coverage memory
is bounded by window size, independent of chromosome length. Intron counts are
owned by the window containing the read start, preventing duplicates when a read
spans multiple windows. Keep SplAdder deletion, overlap, XS, filtering, wrapping
uint32 arithmetic and missing-contig dtypes. Do not import mosdepth's mate
deduplication or its different CIGAR counting policy.

Alternative considered: retain chromosome-wide sparse endpoint maps. Their
memory scales with all distinct alignment endpoints and requires a second
conversion step; bounded windows make the 8 GiB requirement easier to verify
with less implementation complexity. Window re-fetch and serial HDF5 compression
are candidate costs to measure; no speed benefit is assumed.

The prototype is in src/sparse.rs and reads::coverage_window. Tests compare
window sizes 73/1024, 1/4 workers, BAM/CRAM, filtering and stranded/unstranded
arrays against the unmodified upstream summarizer. Larger performance runs use
1 MiB windows. No optimization has yet been adopted into the implementation
branch; parity and resource measurements are required first.

## Validation / pending benchmark

96 valid BAM/CRAM summaries passed (1,008 arrays with exact values, shapes and
dtypes), at 1/4 workers and two window sizes, including the lazy-allocation
optimization for empty windows. Four reference missing-XM failures were checked
as explicit Rust failures rather than counted as successful summaries.
Evidence: ../runs/sparse-parity/report.json. Clippy passes with warnings denied.

The first interleaved benchmark submission failed before any workload ran:
`srun: Unable to contact slurm controller (connect failure)`.
Evidence: ../runs/p0-benchmark-20260911T195314/reference-0.log.
This is an infrastructure failure, not a native performance result. The 4 CPU /
8 GiB benchmark remains required; do not substitute an unenforced local run.

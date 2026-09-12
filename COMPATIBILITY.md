# Compatibility, consistency and differences

Ruspladder v0.1.0 targets **SplAdder v3.1.1**, commit
`65ceec839b9ff0cf96703c1605ee43667662f410`, with the Python dependencies in
`reference-requirements.lock`. The supported release is Linux amd64 with glibc
2.28 or newer. Compatibility means the tested scientific contracts below;
it does **not** mean every input, parameter combination or file byte is identical.

| Area | Preserved behavior | Differences and limits |
| --- | --- | --- |
| Input | GTF/GFF3, indexed BAM/CRAM, reference FASTA, public sparse summaries | No Python pickle-cache import; use fresh output directories when switching programs. |
| Commands | `prep`, `build`, differential `test`; six event types; tested merge, count and output options | No `viz` or diagnostic plots; native help/errors/progress differ. `--annotation-mode` is a Rust extension. Consult each command's help rather than assuming every Python option/abbreviation exists. |
| Core algorithms | Annotation graphs, evidence augmentation, event detection, verification, counts, PSI and statistical testing | Parallel scheduling, storage and native kernels differ; numerical implementation details are in [NUMERICS.md](NUMERICS.md). |
| Annotation | Default grouping and adaptive-filter behavior follow upstream | Explicit `locus` normalization changes affected gene/transcript identities and graphs. Both programs must receive the same normalized GTF for a meaningful parity comparison. |
| Public results | Tested HDF5 names, shapes, dtypes, links, event identities, text fields, coordinates and statistical decisions | Floating-point comparisons use stated tolerances; HDF5 layout and gzip headers are not byte-identity contracts. Some upstream output limitations are retained. |
| Private caches | Native graph/event caches support reuse, recount and downstream testing | HDF5 replaces pickle; filenames, sizes and bytes differ. Rust caches can occupy more disk. |
| Intentional repair | CRAM preparation forwards the reference FASTA to HTSlib | Upstream CRAM prep CLI has missing-attribute wiring; parity there is against the working upstream API. |
| Resources | `--parallel 1..64`; 4 CPU / 4 GiB benchmark, swap disabled | Workload-dependent memory, I/O and scaling. No guarantee of full CPU occupancy, 64-core scaling or universal 4 GiB sufficiency. |

## How consistency is checked

`make test` uses the public fixtures and six scenarios from upstream's `make test`
(`pytest`): positive/negative strand merged and single builds, and 20-sample BAM
and CRAM builds followed by exon-skipping differential tests. Ruspladder replays
these without visualization, uses four workers, and compares against fresh,
unmodified Python CLI runs. The CRAM differential case retains reversed sample
order. The single-sample sparse rerun checks cache reuse; it is not a fresh
sparse-preparation test. The broader `scripts/check.sh` suite separately covers
sparse prep/build, other options, failure cases and all six statistical event types.

Comparators check full graph/event contents, public HDF5 datasets, decompressed
text and differential TSVs. Integer/string values, shapes and dtypes match
exactly. HDF5 floats use `atol=1e-10, rtol=1e-8`, with matching NaNs.
Differential TSV numeric values use `atol=1e-8, rtol=1e-6`, also require equal
six-decimal rounding (as upstream tests do), and preserve rejection decisions at
0.01, 0.05 and 0.1. Passing the fixtures establishes these tested cases, not
universal equivalence or validation on an independent biological cohort.

## Resource baseline

The benchmark runner enforces CPU affinity to four logical CPUs, a four-CPU
quota and **4,294,967,296 bytes (4 GiB)** of memory, with swap disabled. It measures
wall time, child CPU time, average cores, utilization relative to four cores,
maximum process RSS, cgroup peak memory (including charged file cache), and I/O.
Three repetitions interleave Python/Rust runs using fresh and reused application
outputs; filesystem cache is not flushed. Input staging is outside timed stages,
while startup and annotation parsing are inside. Staging and the benchmark driver
can contribute to cgroup memory peak. This is distinct from a cold-disk benchmark.

After building the image, preparing the local reference environment and fetching
the small airway fixture, run:

```sh
export RUSPLADDER_WORK_ROOT="$HOME/data/ruspladder"
"$RUSPLADDER_WORK_ROOT/envs/reference/bin/python" scripts/run_lifecycle_benchmarks.py \
  --work "$RUSPLADDER_WORK_ROOT/runs/baseline-4cpu-4g" \
  --binary "$RUSPLADDER_WORK_ROOT/target/release/ruspladder" \
  --exporter "$RUSPLADDER_WORK_ROOT/target/release/examples/cache_export_probe"
```

Choose a new output directory. `RUSPLADDER_BENCH_CPUS` can select a different set
of four available logical CPUs. The matrix includes the public airway subset and
SplAdder's 20-sample event fixture (direct and sparse build plus differential test).
GitHub-hosted runner allocation is separate: private-repository runners may have
only two available logical CPUs even when four workers are requested. CI is a
correctness gate; the enforced local benchmark supplies resource measurements.

Earlier 8 CPU / 16 GiB measurements on a private 10-million-read subset are not
4 CPU / 4 GiB evidence and are not used to promise resource sufficiency here.

### v0.1.0 measurements

Measured on 2026-09-12 with the packaged glibc 2.28 binary and the enforced
4 CPU / 4 GiB configuration above. All 36 runs completed, and all 18 paired
scientific-output comparisons passed. Values below are medians of three runs;
CPU time includes child processes and threads. Memory is the charged cgroup
peak, including file cache and the benchmark driver.

| Fresh workflow | Python wall (s) | Rust wall (s) | Speedup | Python / Rust CPU time (s) | Python / Rust peak (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Airway, 8 public subset BAMs, build | 29.796 | 0.822 | 36.25× | 97.003 / 0.789 | 835.56 / 90.09 |
| SplAdder, 20 synthetic BAMs, build + test | 39.708 | 1.441 | 27.55× | 105.500 / 1.615 | 839.86 / 85.14 |
| Same event fixture, sparse build + test | 31.289 | 2.336 | 13.40× | 97.230 / 3.133 | 838.53 / 92.28 |

For reused outputs, Python/Rust wall times were 1.317/0.064 s (airway),
3.234/0.065 s (events, direct), and 2.935/0.050 s (events, sparse).
These sub-second measurements include the runner's process-wait polling overhead;
they are observed workflow times, not kernel microbenchmarks.
Fresh Rust runs averaged 0.95–1.33 CPU cores (24–33% of the four-CPU quota),
with median maximum process RSS of 83–89 MiB. Native cgroup peak never exceeded
93.3 MiB across this small-fixture matrix. These inputs emphasize startup,
annotation and orchestration costs; the figures do not predict whole-transcriptome
or 10-million-read workloads. Rust also wrote more data on fresh runs
(3.05/4.82/6.18 MiB versus Python's 1.44/0.98/3.51 MiB for the three rows).
Reduced elapsed time and CPU work do not imply higher CPU occupancy or less I/O.

The release package also passed all six public make-test scenarios: 56 graph
cache files, 36 event collections, 558 HDF5 datasets, 152 decompressed text files
and six differential TSVs. Relocation was exercised by actual builds in a
Debian 12 container without Python and directly on a glibc 2.28 host.

## Detailed observed behavior

This migration targets SplAdder v3.1.1, commit
`65ceec839b9ff0cf96703c1605ee43667662f410`, with dependencies pinned in
`reference-requirements.lock`. These are observed contracts and failures.

- Text starts use the source's mixed coordinate conventions in each format.
  Event IDs, annotation bits and isoform order retain upstream behavior.
- Structured output prints Python byte-string representations for sample tags.
  It stops after an unterminated prefix for the first multi-exon-skip event.
  BED output stops without records for multi-exon-skip and mutex events.
  These limitations are exercised in the 990 exact text comparisons.
- Event HDF5 positions for mutex events have shape (events, 2, 4), unlike other
  types. Empty event sets use the source's one-element event_counts sentinel
  and metadata, without invented empty PSI/confirmation datasets.
- Single-quantification collection copies metadata and numeric dtypes from the
  first file, extends the sample axes, truncates labels to S255, and has no
  strains soft link. Ordinary graph/event/expression files do have that link.
- The differential-test quantifier restores requested coverage/sample order
  after sorting within each group, but leaves PSI in sorted group order.
  The 324-case comparison includes reversed groups, filter_idx, fractional
  features, empty groups and the source's high-memory mode.
- Normalization `uq` raises NameError because count.py does not import
  scoreatpercentile. Eighty reference failures are recorded separately from
  160 successful normalization comparisons. Native code reports this error.
- Correction `TSBH` passes `tsbh` to statsmodels 0.14.4, which rejects that
  method name. Correction on no non-NaN hypotheses also fails upstream
  by division by zero. Native code reports errors; 29 unrecognized-method and
  91 empty-input reference cases are separate from 174 successful comparisons.
  The six working methods retain missingness, six-decimal results, and tested
  rejection decisions at alpha 0.01, 0.05 and 0.1.
- NB/log and Gamma/identity IRLS reuse statsmodels' iteration, scale and
  stopping rules, NumPy's LAPACK least-squares and pseudoinverse kernels,
  and SciPy's bounded scalar optimizer. The Cox–Reid term retains the LU
  determinant followed by log; replacing it with a stable log-determinant
  changes observed estimates for zero-heavy, ill-conditioned designs.
- SciPy 1.13.1's `xlog1py` calls its renamed Cephes `log1p`, despite the
  Cython source importing `libc.math.log1p`. The wheel's disassembly confirms
  the call to `cephes_log1p`. Reusing unity.c's approximation removes small
  likelihood differences that otherwise perturb dispersion optimization.
  Source and SHA-256 are retained in vendor/scipy-special/SOURCES.json.
- The final test chooses the larger available isoform p-value, retaining
  the first on ties; it replaces missing p-values with 1. LRT uses 1-CDF
  with one degree of freedom, preserving cancellation at very small values.
- JSON comparison probes enable exact float round trips. The default parser
  changed individual inputs by one ULP, enough to amplify differences in
  the ill-conditioned dispersion fixture. HDF5 production values do not
  pass through JSON.
- Test input processing retains the source's two different outlier caps
  (expression: 1.5 IQR after normalization; events: 3 IQR before normalization),
  low-coverage union across isoforms, all-missing group PSI replacement,
  ties-to-even count rounding, design columns and dPSI mask. The 144-case
  comparison captures the actual upstream test function at the GLM boundary.
- Test output retains NumPy's unstable floating-point argsort for equal
  p-values, first-row-per-gene selection, coordinate strings and all TSV
  columns. Internal test_setup pickle dictionaries become HDF5 with the same
  values. The 48 writer comparisons inject fixed statistical results into
  the actual upstream test function; they do not count as full CLI tests.
- The real gene-expression column selection and quantified PSI are Fortran-order
  arrays in NumPy. Their row reductions use sequential addition; count matrices
  use pairwise addition. The complete CLI option check exposed a 5.6e-6 p-value
  error before this distinction was reproduced. The input mock now preserves
  the actual array order, and the CLI checks include non-alt normalization,
  both outlier switches, high-memory input, labels, tags and reversed CRAM groups.
- The upstream boolean `--timestamp` is compared to the string `y`; the CLI
  consequently never adds a timestamp. Native naming preserves this behavior.
  A one-entry condition text file fails upstream with `iteration over a 0-d
  array`; native parsing reports the corresponding singleton-list error.
- The committed mutex fixture has too few distinct raw dispersions after
  percentile trimming, leaving an empty Gamma trend design. Upstream and native
  test commands both fail explicitly. A separate heterogeneous 64-event
  fixture successfully tests mutex and all five other event types through
  the full command. These upstream failures are not counted as successful fits.

The CLI test path explicitly disables exon-count augmentation of isoform counts
and event-ID construction. Legacy direct graph quantifiers and unused experimental
helper paths are not substituted for that production path. The reproduced advanced graph failures are listed below.

Primary implementation references: SplAdder v3.1.1 count.py,
alt_splice/{analyze,quantify,write}.py and spladder_test.py; statsmodels 0.14.4
stats/multitest.py. The HDF5 chunk decision follows the
[HDF Group partial-I/O guidance](https://support.hdfgroup.org/documentation/hdf5/latest/hdf5_chunking.html).
Gzip output uses [flate2's streaming encoder](https://docs.rs/flate2/1.1.10/flate2/write/struct.GzEncoder.html)
over the already linked zlib; gzip headers are not a byte-identity interface.

Sequential graph generation shares read-filter state across samples. The source
feasibility loop raises the minimum exon length/support when a gene has more
than 200 introns. The first retention stage saves a copy into `read_filter`,
detaching it from the filter dictionary retained by `intron_retention`.
Later samples use the current filter for intron selection/cassette coverage and
the retained first-stage filter for retention coverage. The regression with
201 distinct introns changes the current exon minimum from 17 to 21, while
retention stays at 17; both samples insert one retention. Treating these filters
as one value lost the second retention and changed the graph and alternative
label. The native options now retain the same state across sequential builds.

Build orchestration preserves cache-existence reuse and the source's file naming
for public counts, expression and event outputs. Full graph/event caches use
`.hdf5` / `.events.hdf5` in place of internal pickle. Collected count labels retain
S255 width in downstream gene-expression HDF5. Diagnostic progress/log formatting
is native; scientific text files are compared byte-for-byte after decompression.

The following build paths are separately reproduced as upstream failures:
`--re-infer-sg` calls an undefined function; `single` with multiple samples lacks
an event cache for the second sample; `--qmode single` with multiple supplied
samples later indexes beyond its one count column; a nonfinal chunk with default
quantification tries to load the absent final graph; and graph validation with
`merge_bams` requests a validated graph that the source never creates. These are
explicit native errors. Successful external chunk workflows disable extraction
and quantification until the final level; per-sample quantification supplies one
sample per call and disables event extraction until collection.


Sparse alignment prep uses bounded windows and preserves the public COO names,
shapes, dtypes and values. The serial reference leaves `_reads_shp` uncompressed;
its parallel collector compresses that dataset too, which native prep preserves.
Chunks and unlimited append extents are storage choices, not value interfaces.
The source ignores `--ignore-mismatches` in sparse prep; missing NM tags still
fail there even when direct-alignment build succeeds with that option.

The upstream CRAM prep CLI refers to an absent `cram_ref` attribute. Native prep
forwards `--reference` to HTSlib. Twelve native arrays are compared against the
unmodified working upstream `summarize_chr` API supplied with its expected
attribute. This is an intentional CLI wiring repair, not a successful upstream
CRAM CLI comparison. The direct-CRAM build CLI remains separately verified.

Sparse multi-BAM graph generation reuses the first file's chromosome cache.
Its duplicate-intron loop also overwrites the gene index, sometimes assigning
introns to another gene or failing out of bounds. The migration retains those
observed results and explicit failures; eight valid core builds and four source
IndexErrors are checked separately. Normal `merge_graphs` builds each sample
before merging, and the complete sparse CLI suite covers all four merge modes.
Sparse counting retains its strand behavior, which differs from the direct-BAM
helper's default treatment of untagged junctions.

When a BAM is absent, `--sparse-bam` accepts its complete `.hdf5` summary without
requiring the original BAM index. Existing graphs can be recounted; rebuilding
from annotation also preserves `init_regions` skipping absent BAMs. Both paths
are compared against source output using isolated copies with BAMs removed.

## Annotation

The default `spladder` mode groups annotation records by gene ID, retaining
upstream behavior even when an ID spans reference sequences or strands.
Explicit `locus` mode normalizes GTF IDs by `(gene_id, seqname, strand)` before
import and writes an identity mapping. Coordinates and non-ID attributes are
preserved. Same-reference, same-strand records are not split by distance.
Compare both programs on the same normalized GTF: normalization can change
graphs and adaptive filters. GFF3 uses feature IDs and Parent relationships.

NCBI biological GeneIDs can occur on multiple features; they are not equivalent
to unique GFF3 feature IDs. See the [NCBI annotation format documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).

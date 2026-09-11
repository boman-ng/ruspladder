# Observed SplAdder v3.1.1 behavior

This migration targets the pinned source and reference environment specified in
MIGRATION.md. These are observed contracts and failures, not proposed fixes.

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

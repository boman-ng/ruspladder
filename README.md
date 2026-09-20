# Ruspladder

Ruspladder is an independent Rust implementation of the non-visual analysis in
[SplAdder](https://github.com/ratschlab/spladder) v3.1.1. It takes gene annotation
and RNA-seq alignments, builds and augments splice graphs, detects alternative
splicing events, quantifies their support and PSI, and tests differences between
conditions. It implements `prep`, `build` and differential `test`, including exon
skipping, intron retention, alternative 3′/5′ splice sites, multiple-exon skipping
and mutually exclusive exons.

The scientific methods originate in SplAdder. Ruspladder is maintained by
[boman-ng](https://github.com/boman-ng); it is not an official release of SplAdder,
ratschlab or OpenGene. Its focus is preserving tested scientific behavior while
reducing runtime and memory use through Rust, bounded data processing and
multithreading. See [compatibility and differences](COMPATIBILITY.md) before
substituting it in an existing workflow.

## Install the Linux binary

Download `ruspladder-v0.1.0-linux-amd64.tar.gz` and its `.sha256` file from
[Releases](https://github.com/boman-ng/ruspladder/releases/tag/v0.1.0).
The release is public; no GitHub account is required:

```sh
curl -fLO https://github.com/boman-ng/ruspladder/releases/download/v0.1.0/ruspladder-v0.1.0-linux-amd64.tar.gz
curl -fLO https://github.com/boman-ng/ruspladder/releases/download/v0.1.0/ruspladder-v0.1.0-linux-amd64.tar.gz.sha256
sha256sum -c ruspladder-v0.1.0-linux-amd64.tar.gz.sha256
tar -xzf ruspladder-v0.1.0-linux-amd64.tar.gz
cd ruspladder-v0.1.0-linux-amd64
export OPENBLAS_NUM_THREADS=1
./ruspladder --version
```

Requires **Linux amd64 (x86_64), glibc 2.28 or newer**, including RHEL/AlmaLinux 8, Debian 10+ and
Ubuntu 20.04+. Keep the executable and adjacent `lib/` directory together; the
whole directory can be moved. Python, Rust, a compiler and a system BLAS/HDF5
installation are unnecessary at runtime. Alpine/musl and ARM are outside this
release's supported platform. `OPENBLAS_NUM_THREADS=1` prevents the bundled BLAS
loader from creating an idle thread pool; `--parallel` controls analysis workers.

## Quick start

Use coordinate-sorted, indexed BAM files and a matching GTF/GFF3 annotation.
For CRAM, also supply `--reference genome.fa`.

```sh
./ruspladder build -a annotation.gtf \
  -b control1.bam,control2.bam,treated1.bam,treated2.bam \
  -o results --parallel 4 --readlen 150 --output-txt

./ruspladder test -o results \
  -a control1,control2 -b treated1,treated2 --parallel 4 --readlen 150
```

Sample names for `test` correspond to the alignment basenames without extensions.
Set `--readlen` to your data's read length. Public results include graph and gene
expression counts, event counts/PSI in HDF5, optional event text files, and
non-visual differential-testing TSV files. Use `./ruspladder COMMAND --help` for
parameters. Internal Rust caches cannot be exchanged with Python pickle caches.

For reusable sparse alignment summaries:

```sh
./ruspladder prep -a annotation.gtf -b sample.bam --sparse-bam --parallel 4
./ruspladder build -a annotation.gtf -b sample.bam -o results-sparse \
  --sparse-bam --parallel 4
```

For GTFs reusing a gene ID across reference sequences or strands, explicit
`--annotation-mode locus` writes normalized `.locus.gtf` and `.loci.tsv`
companions beside the annotation. It preserves coordinates and records original
identities. This changes affected gene models: compare Python and Rust on the
**same normalized GTF**. Default `--annotation-mode spladder` retains upstream
import behavior. Use separate result directories for the two modes. Annotation
files and their caches are treated as immutable; after an edit, use a new path or
remove the generated companions. See [annotation details](COMPATIBILITY.md#annotation).

## Resources and verification

The resource baseline is **4 logical CPUs / 4 GiB RAM**, with a hard memory limit
and swap disabled. Commands accept 1–64 threads. Memory requirements depend on
annotation spans, sample count and event density; 4 GiB is the benchmark
configuration, not a cap enforced by the executable or a promise for every input.

`make test` replays the six non-visual scenarios in SplAdder's public
[`tests/test_end_to_end.py`](https://github.com/ratschlab/spladder/blob/v3.1.1/tests/test_end_to_end.py):
positive/negative strand, merged/single analysis, and 20-sample BAM/CRAM build
plus differential testing. Both programs run from fresh annotations; graphs,
events, HDF5 values and text/TSV results are compared. CI tests the release package
and its relocation. [COMPATIBILITY.md](COMPATIBILITY.md) defines comparison
precision and the limits of this evidence.

## Build and test from source

Build/reference caches and data default to `~/data/ruspladder`; override with
`RUSPLADDER_WORK_ROOT`. The supported build uses the pinned Linux toolchain:

```sh
export RUSPLADDER_WORK_ROOT="$HOME/data/ruspladder"
mkdir -p "$RUSPLADDER_WORK_ROOT"
docker build -f Dockerfile.build -t ruspladder-build:0.1.0 .
docker run --rm --cpus 4 --memory 4g --memory-swap 4g \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/src" -v "$RUSPLADDER_WORK_ROOT:/work" \
  ruspladder-build:0.1.0 bash scripts/ci.sh
```

This builds, checks, packages and compares public data. The runtime archive,
corresponding native-source archive, checksums and verification JSON are written
to `$RUSPLADDER_WORK_ROOT/dist/`. See [redistribution materials](licenses/README.md). Use a fresh work directory when rerunning this
complete packaging/CI command. Python is used only for dependency extraction,
packaging and reference comparisons. [NUMERICS.md](NUMERICS.md) describes why the
native numerical dependencies are pinned.

With Rust ≥1.89, C/C++17, Clang/libclang, CMake ≥3.26, zlib/bzip2 development
headers and Python ≥3.11 installed, `make build` builds locally.
The maintenance scripts use Python 3.11 standard-library hashing and TOML APIs. To run comparisons,
install `uv`, run `bash scripts/bootstrap.sh`, then `make test`. The broader
comparison suite is available through `bash scripts/check.sh`
after `python3 scripts/fetch_fixtures.py --work-root "$RUSPLADDER_WORK_ROOT"`.

## Attribution, citation and license

Thank you to André Kahles, Cheng Soon Ong, Yi Zhong and Gunnar Rätsch for the
SplAdder methods and publication, and to the SplAdder authors and contributors
for making their implementation and test data available. Ruspladder derives its
algorithms from SplAdder v3.1.1, commit
`65ceec839b9ff0cf96703c1605ee43667662f410`. It also reuses work from NumPy,
SciPy, statsmodels, OpenBLAS/LAPACK, HTSlib, HDF5 and x86-simd-sort; source notices
and provenance are retained in `licenses/` and `vendor/`.

When using Ruspladder in research, cite the original method and record the
Ruspladder version:

Kahles A, Ong CS, Zhong Y, Rätsch G. **SplAdder: identification, quantification and
testing of alternative splicing events from RNA-Seq data.** *Bioinformatics*
32(12), 1840–1847 (2016). [doi:10.1093/bioinformatics/btw076](https://doi.org/10.1093/bioinformatics/btw076).
Machine-readable citation metadata is in [CITATION.cff](CITATION.cff).

Ruspladder inherits SplAdder's **BSD 3-Clause** license, retaining the original
copyright and adding attribution for the Rust implementation. See [LICENSE](LICENSE).
Third-party components retain their own licenses; binary archives include their
notices. Upstream attribution does not imply endorsement.

## Contributing and maintenance

See [CONTRIBUTING.md](CONTRIBUTING.md) for development checks and scientific
compatibility requirements, and [SECURITY.md](SECURITY.md) for private reports.

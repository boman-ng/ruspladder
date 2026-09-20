# v0.1.0

Initial Linux amd64 release of the independent Rust implementation of
SplAdder v3.1.1's non-visual analysis: annotation/alignment preparation,
splice-graph construction and augmentation, all six event types, quantification,
and differential testing.

The archive includes the executable, pinned native libraries and license notices.
Extract the complete directory; run `OPENBLAS_NUM_THREADS=1 ./ruspladder --help`.
Requires Linux amd64 with glibc 2.28 or newer. Python is not needed at runtime.

Public test scenarios follow SplAdder's `make test`. The resource benchmark is
4 logical CPUs and a 4 GiB hard memory limit with swap disabled. See
`COMPATIBILITY.md` for tested contracts, numerical tolerances, cache differences
and annotation behavior; this release does not claim universal byte identity
or provide visualization.

Method credit belongs to the SplAdder authors; please cite Kahles et al.,
Bioinformatics (2016), doi:10.1093/bioinformatics/btw076. Ruspladder is maintained
by boman-ng and distributed under BSD-3-Clause, retaining upstream notices.

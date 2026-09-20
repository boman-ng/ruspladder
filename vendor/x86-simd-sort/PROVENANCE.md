Unmodified src/ headers and license from intel/x86-simd-sort commit
9a1b616d5cd4eaf49f7664fb86ccc1d18bad2b8d, the submodule revision vendored by
NumPy 2.2.6. Download: https://github.com/intel/x86-simd-sort/tree/9a1b616d5cd4eaf49f7664fb86ccc1d18bad2b8d

Used for NumPy-compatible argsort order, including equal keys. Rust stable sort
does not preserve NumPy's equal-key ordering. The native wrapper follows
numpy/_core/src/npysort/x86_simd_argsort.dispatch.cpp.

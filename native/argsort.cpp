// Use the same library revision and call as NumPy 2.2.6, to preserve
// observable ordering of equal keys. The scientific algorithms remain Rust.
#include <cstdint>
#include "x86simdsort-static-incl.h"

extern "C" void RUSPLADDER_ARGSORT(int64_t *values, size_t *indices, size_t size) {
    x86simdsortStatic::argsort(values, indices, size, true);
}
extern "C" void RUSPLADDER_ARGSORT_F64(double *values, size_t *indices, size_t size) {
    x86simdsortStatic::argsort(values, indices, size, true);
}

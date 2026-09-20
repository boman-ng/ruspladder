// Reuse the pinned SciPy/Cephes functions for dispersion likelihood parity.
#include "cephes/gamma.h"
#include "cephes/zeta.h"

extern "C" double ruspladder_lgamma(double x) { return special::cephes::lgam(x); }
extern "C" double ruspladder_trigamma(double x) { return special::cephes::zeta(2.0, x); }
// SciPy 1.13.1 _xlogy.pxd's log1p call is renamed to cephes_log1p by
// cephes_names.h. Preserve unity.c's approximation, not the system log1p.
extern "C" double ruspladder_xlog1py(double y, double x) {
    if (y == 0.0 && !std::isnan(x)) return 0.0;
    static const double LP[] = {
        4.5270000862445199635215E-5, 4.9854102823193375972212E-1,
        6.5787325942061044846969E0, 2.9911919328553073277375E1,
        6.0949667980987787057556E1, 5.7112963590585538103336E1,
        2.0039553499201281259648E1,
    };
    static const double LQ[] = {
        1.5062909083469192043167E1, 8.3047565967967209469434E1,
        2.2176239823732856465394E2, 3.0909872225312059774938E2,
        2.1642788614495947685003E2, 6.0118660497603843919306E1,
    };
    double z = 1.0 + x;
    if ((z < M_SQRT1_2) || (z > M_SQRT2)) return y * std::log(z);
    z = x * x;
    z = -0.5 * z + x * (z * special::cephes::polevl(x, LP, 6) /
                        special::cephes::p1evl(x, LQ, 6));
    return y * (x + z);
}
// The production LRT always has one degree of freedom. DLMF 8.4.1 gives
// P(1/2, x/2) = erf(sqrt(x/2)); use the platform math-library implementation.
extern "C" double ruspladder_chi2_cdf1(double x) {
    return x < 0.0 ? 0.0 : std::erf(std::sqrt(x / 2.0));
}

#include "fastmath.hpp"

#include <boost/math/special_functions/log1p.hpp>

#include "jacoblog.hpp"
#include "icsilogw.hpp"
#include "PowFast.hpp"


/** Switch to disable fast math routines for debugging. */
static const bool DisableFastMath = true;


/** Allocate precomputed tables. */
static ICSILog flog(14);
static PowFast fpow(14);
static JacobianLogTable fjacoblog(1 << 11);


/** Fast, approxiate jacobian logarithm: log(1+exp(x)) */
double fast_jacoblog(double x) {

    // Avoid numerical issues arising from extreme values 
    if (x >= 60.0) return x;
    else if (x < -60.0) return 0.0;
    else {
        if (DisableFastMath) 
            return boost::math::log1p(std::exp(x));
        else
            return fjacoblog.jacobianLog(x);
    }
}


/** Fast, approximate natural logarithm. */
double fast_log(double x) {

    if (DisableFastMath) return std::log(x);

    return flog.log(static_cast<float>(x));
}


/** Fast, approximate exponentiation. */
double fast_exp(double x) {

    if (DisableFastMath) return std::exp(x);

    return fpow.e(static_cast<float>(x));
}

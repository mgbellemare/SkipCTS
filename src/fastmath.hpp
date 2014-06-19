#ifndef __FASTMATH_HPP__
#define __FASTMATH_HPP__

/** A fast, approxiate jacobian logarithm function that computes log(1+exp(x)). */
extern double fast_jacoblog(double x);

/** A fast, approximate natural logarithm function. */
extern double fast_log(double x);

/** A fast, approximate exponentiation function. */
extern double fast_exp(double x);

#endif // __FASTMATH_HPP__


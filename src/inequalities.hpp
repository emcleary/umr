#ifndef INEQUALITIES_H
#define INEQUALITIES_H

#include <cmath> // abs
#include <functional> // less, hash
#include <cassert>
#include <string>


void debug_init();


namespace umr {

namespace inequalities {


static const double RTOL = 1e-6;
static const double ATOL = 1e-8;

/**
 * Check if \f$ a \approx b \f$.
 *
 * @param a Value
 * @param b Value
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 */
inline bool is_close(double a, double b, double rtol = RTOL, double atol = ATOL) {
    return std::abs(a - b) <= atol + std::abs(b) * rtol;
}

/**
 * Check if \f$ a < b \f$.
 *
 * @param a Value
 * @param b Value
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 */
inline bool is_lt(double a, double b, double rtol = RTOL, double atol = ATOL) {
    return a < b && !is_close(a, b, rtol, atol);
}

/**
 * Check if \f$ a \leq b \f$.
 *
 * @param a Value
 * @param b Value
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 */
inline bool is_le(double a, double b, double rtol = RTOL, double atol = ATOL) {
    return a < b || is_close(a, b, rtol, atol);
}

/**
 * Check if \f$ a > b \f$.
 *
 * @param a Value
 * @param b Value
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 */
inline bool is_gt(double a, double b, double rtol = RTOL, double atol = ATOL) {
    return a > b && !is_close(a, b, rtol, atol);
}

/**
 * Check if \f$ a \geq b \f$.
 *
 * @param a Value
 * @param b Value
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 */
inline bool is_ge(double a, double b, double rtol = RTOL, double atol = ATOL) {
    return a > b || is_close(a, b, rtol, atol);
}


} // namespace inequalities

} // namespace umr

#endif // INEQUALITIES_H

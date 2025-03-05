#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "integration_autogk.hpp"


namespace umr {

namespace math {

namespace integrate {


/**
 * Acts as an integrand for calculating integrals.
 *
 * @param t A parameter
 * @param t_minus_tmin Difference between t and its lower bound
 * @param t_minus_tmin Difference between t's upper bound and itself
 * @param result The result
 * @param ptr Assumed to be an IParametric pointer
 */
void segment_integrand(double t, double t_minus_tmin, double tmax_minus_t, double &result, void *ptr);

/**
 * @class SegmentLength
 * @brief Calculates the line integral
 *
 * Calculates the line integral of an IParametric object using the
 * AlgLib library.
 */
class SegmentLength : public AutoGK {
public:
    /** Constructor for SegmentLength */
    SegmentLength() : AutoGK(segment_integrand) {}

    /**
     * Initializes the object for integration.
     *
     * @param ptr A pointer to an IParametric object
     * @param tmin The lower bound of integration
     * @param tmax The upper bound of integration
     */
    virtual void initialize(void* ptr, double& tmin, double& tmax) final;
};


} // namespace integrate

} // namespace math

} // namespace umr

#endif // INTEGRATORS_H

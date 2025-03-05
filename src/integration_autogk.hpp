#ifndef INTEGRATION_AUTOGK_H
#define INTEGRATION_AUTOGK_H

#include "../ext/alglib-cpp/src/integration.h"


namespace umr {

namespace math {

namespace integrate {


/** A type for a function to be integrated */
using integrand_func = void(*)(double t, double t_minus_tmin, double tmax_minus_t,
        double &result, void *ptr);


/**
 * @class AutoGK
 * @brief Calculates an integral
 *
 * A wrapper for calculating integrals using the AlgLib library.
 */
class AutoGK {
public:
    /**
     * Constructor for AutoGK
     *
     * @param integrand A function representing the integrand
     */
    AutoGK(integrand_func integrand) : m_integrand_function(integrand) {};

    /** Destructor for AutoGK */
    ~AutoGK() {};

    /**
     * Evaluates the integral
     *
     * @param ptr A pointer passed to the integrand function
     */
    double run(void* ptr);

private:
    /** Set upper/lower bounds of the parameter for integration */
    virtual void initialize(void* ptr, double& tmin, double& tmax) = 0;

    integrand_func m_integrand_function;
};


} // namespace integrate

} // namespace math

} // namespace umr

#endif // INTEGRATION_AUTOGK_H

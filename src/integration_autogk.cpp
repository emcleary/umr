#include "integration_autogk.hpp"

#include "../ext/alglib-cpp/src/integration.h"


namespace umr {

namespace math {

namespace integrate {


double AutoGK::run(void* ptr) {
    double result;

    try {
        double tmin;
        double tmax;
        initialize(ptr, tmin, tmax);

        alglib::autogkstate state;
        alglib::autogkreport report;

        alglib::autogksmooth(tmin, tmax, state);
        alglib::autogkintegrate(state, m_integrand_function, ptr);
        alglib::autogkresults(state, result, report);

    } catch(alglib::ap_error alglib_exception) {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        exit(-1);
    }

    return result;
}


} // namespace integrate

} // namespace math

} // namespace umr

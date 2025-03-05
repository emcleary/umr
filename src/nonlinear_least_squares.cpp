#include "nonlinear_least_squares.hpp"

#include <cassert>

#include "../ext/alglib-cpp/src/optimization.h"


namespace umr {

namespace math {

namespace optimize {


std::vector<double> NonlinearLeastSquares::run(void* ptr) {
    assert(m_num_params >= 0 && "NonlinearLeastSquares::run: Number of parameters must be set!");
    assert(m_num_params <= m_num_resid
            && "NonlinearLeastSquares::run: Infinite solutions, more parameters than residuals!");

    m_params.setlength(m_num_params);
    m_scale.setlength(m_num_params);
    m_lower_bound.setlength(m_num_params);
    m_upper_bound.setlength(m_num_params);

    initialize(ptr);

    try {
        // Create optimizer

        // derivative-free nonlinear least squares; first arg is size of fi, NOT size of x
        alglib::nlscreatedfo(m_num_resid, m_params, m_state);
        // stopping conditions; maxits=0 --> unlimited, epsx ->
        // stopping trust region radius; both 0 -> automatic stopping
        // criteria (small EpsX)
        alglib::nlssetcond(m_state, m_epsilon, m_maxits);
        // scale s; for my application setting all to 1 should be
        // fine; TODO: test with high changes in curvature
        alglib::nlssetscale(m_state, m_scale);
        alglib::nlssetbc(m_state, m_lower_bound, m_upper_bound);

        // Select algorithm
        // alglib::nlssetalgodfolsa(m_state); // fewer target evaluations with no parallelism support
        alglib::nlssetalgo2ps(m_state); // more target evaluations but can support parallelism

        // Optimize (serially; parallelization requires additional alglib::parallelcallbacks parameter)
        alglib::nlsoptimize(m_state, m_func_resid, nullptr, ptr);

        // Test optimization results
        alglib::nlsresults(m_state, m_params, m_report);

    } catch(alglib::ap_error alglib_exception) {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        exit(-1);
    }

    return finalize(ptr);
}


} // namespace optimize

} // namespace math

} // namespace umr

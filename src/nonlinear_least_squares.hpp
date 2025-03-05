#ifndef NONLINEAR_LEAST_SQUARES_H
#define NONLINEAR_LEAST_SQUARES_H

#include "../ext/alglib-cpp/src/optimization.h"


namespace umr {

namespace math {

namespace optimize {


/** A type for a residual function used while optimizing */
using residual_func = void(*)(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr);


/**
 * @class NonlinearLeastSquares
 * @brief Optimizes parameters
 *
 * A wrapper for optimizing parameters using the nonlinear least squares
 * algorithm of the AlgLib library.
 */
class NonlinearLeastSquares {
public:
    /**
     * Constructor for NonlinearLeastSquares
     *
     * @param func_resid A residual function
     */
    NonlinearLeastSquares(residual_func func_resid) : m_func_resid(func_resid) {}

    /** Destructor for NonlinearLeastSquares */
    ~NonlinearLeastSquares() {}

    /**
     * Sets the number of parameters to be optimized. The user must
     * ensure that the number of parameters is less than the number of
     * residuals.
     */
    void set_num_params(size_t n) { m_num_params = n; }

    /**
     * Sets the number of residuals to be evaluated. The user must
     * ensure that the number of parameters is less than the number of
     * residuals.
     */
    void set_num_resid(size_t n) { m_num_resid = n; }

    /** Set the max number of iterations */
    void set_num_iter(alglib::ae_int_t num_iter) { m_maxits = num_iter; }

    /** Set the terminating trust region radius */
    void set_epsilon(double eps) { m_epsilon = eps; }

    // passing model in ensures initialize and finalize use the same model

    /**
     * Runs the optimization algorithm
     *
     * @param ptr A pointer passed to initialize, finalize, and the
     * optimization algorithm
     * @return A vector of optimized parameter values
     */
    std::vector<double> run(void* ptr);

protected:
    /**
     * This function initializes arrays and variables used by the optimizer.
     * The user is responsible to define this function. It must set:
     *
     * - Parameter bounds m_lower_bound, m_upper_bound
     * - Initial parameters m_params
     * - Scale m_scale
     *
     * @param ptr The same pointer passed to run
     */
    virtual void initialize(void* ptr) = 0;

    /**
     * This function is intended to postprocess the parameters,
     * converting them from alglib::read_1d_array to vector type.
     *
     * @param ptr The same pointer passed to run
     */
    virtual std::vector<double> finalize(void* ptr) = 0;

    int m_num_params, m_num_resid;
    alglib::ae_int_t m_maxits = 0;
    double m_epsilon = 0;
    alglib::nlsstate m_state;
    alglib::nlsreport m_report;
    alglib::real_1d_array m_scale;
    alglib::real_1d_array m_params;
    alglib::real_1d_array m_lower_bound, m_upper_bound;
    residual_func m_func_resid;
};


} // namespace optimize

} // namespace math

} // namespace umr

#endif // NONLINEAR_LEAST_SQUARES_H

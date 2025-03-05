#ifndef OPTIMIZERS_H
#define OPTIMIZERS_H

#include <cassert>

#include "nonlinear_least_squares.hpp"
#include "point.hpp"
#include "parametrics.hpp"
#include "quadedge.hpp"


/** A pair of doubles */
using Real2 = std::pair<double, double>;


namespace umr {

namespace math {

namespace optimize {


////////////////////////
// Residual Functions //
////////////////////////

/**
 * A residual function used in SplitSegmentUniformly
 *
 * @param t An array of parameters
 * @param fi An array of residuals
 * @param ptr A pointer passed to the function; used as needed
 */
void segment_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr);

/**
 * A residual function used in ShortestDistanceParametricParametric
 *
 * @param t An array of parameters
 * @param fi An array of residuals
 * @param ptr A pointer passed to the function; used as needed
 */
void shortest_distance_parametric_parametric_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr);

/**
 * A residual function used in ShortestDistanceParametricPoint
 *
 * @param t An array of parameters
 * @param fi An array of residuals
 * @param ptr A pointer passed to the function; used as needed
 */
void shortest_distance_parametric_point_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr);

/**
 * A residual function used in SplitSegmentAtLength
 *
 * @param t An array of parameters
 * @param fi An array of residuals
 * @param ptr A pointer passed to the function; used as needed
 */
void split_segment_at_length_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr);



///////////////////////
// Optimizer Classes //
///////////////////////

/**
 * @class SplitSegmentUniformly
 * @brief Splits a segment into uniform lengths
 *
 * This class uses an the nonlinear least squares algorithm of the
 * AlgLib library to find the best parameters of an IParametric object
 * that split its parametric function into uniform lengths. Lengths
 * are measured by the Euclidean distance between each pair of
 * consecutive points.
 */
class SplitSegmentUniformly : public NonlinearLeastSquares {
public:

    /** Constructor for SplitSegmentUniformly */
    SplitSegmentUniformly() : NonlinearLeastSquares(segment_residual) {}

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
    virtual void initialize(void* ptr) final;


    /**
     * This function is intended to postprocess the parameters,
     * converting them from alglib::read_1d_array to vector type.
     *
     * @param ptr The same pointer passed to run
     * @return A vector of parameters
     */
    virtual std::vector<double> finalize(void* ptr) final;

    /**
     * Set the number of subsegments to split the parametric into.
     * This method sets the number of parameters and number of
     * residuals for the optimizer.
     *
     * @param n_subsegments Number of subsegments
     */
    void set_num_subsegments(size_t n_subsegments);
};


/**
 * @class ShortestDistanceParametricParametric
 * @brief Find the shortest distance between two parametrics
 *
 * This class uses an the nonlinear least squares algorithm of the
 * AlgLib library to find the shortest distance between parametrics of
 * two IParametric objects.
 */
class ShortestDistanceParametricParametric : public NonlinearLeastSquares {
public:

    /** Constructor for ShortestDistanceParametricParametric */
    ShortestDistanceParametricParametric()
            : NonlinearLeastSquares(shortest_distance_parametric_parametric_residual) {
        m_num_params = 2;
        m_num_resid = 2;
    }

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
    virtual void initialize(void* ptr) final;

    /**
     * This function is intended to postprocess the parameters,
     * converting them from alglib::read_1d_array to vector type.
     *
     * @param ptr The same pointer passed to run
     * @return A vector of parameters
     */
    virtual std::vector<double> finalize(void* ptr) final;

    /** Sets a parametric */
    void set_parametric_0(parametric::IParametric* p) { m_param0 = p; }

    /** Sets a parametric */
    void set_parametric_1(parametric::IParametric* p) { m_param1 = p; }

    /** Gets a parametric */
    const parametric::IParametric* const get_parametric_0() { return m_param0; }

    /** Gets a parametric */
    const parametric::IParametric* const get_parametric_1() { return m_param1; }

private:
    parametric::IParametric* m_param0;
    parametric::IParametric* m_param1;
};


/**
 * @class ShortestDistanceParametricPoint
 * @brief Find the shortest distance between a parametrics and a point
 *
 * This class uses an the nonlinear least squares algorithm of the
 * AlgLib library to find the shortest distance between the parametric
 * of an IParametric object and a Point object.
 */
class ShortestDistanceParametricPoint : public NonlinearLeastSquares {
public:
    /** Constructor for ShortestDistanceParametricPoint */
    ShortestDistanceParametricPoint()
            : NonlinearLeastSquares(shortest_distance_parametric_point_residual) {
        m_num_params = 1;
        m_num_resid = 2;
    }

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
    virtual void initialize(void* ptr) final;

    /**
     * This function is intended to postprocess the parameters,
     * converting them from alglib::read_1d_array to vector type.
     *
     * @param ptr The same pointer passed to run
     * @return A vector of parameters
     */
    virtual std::vector<double> finalize(void* ptr) final;

    /** Set the parametric */
    void set_parametric(parametric::IParametric* p) { m_param = p; }

    /** Set the point */
    void set_point(Real2 p) { m_point = p; }

    /** Set the point */
    void set_point(mesh::Point* p) { m_point = {p->x, p->y}; }

    /** Get the parametric */
    const parametric::IParametric* const get_parametric() { return m_param; }

    /** Get the point */
    const Real2& get_point() { return m_point; }

private:
    parametric::IParametric* m_param;
    Real2 m_point;
};


/**
 * @class SplitSegmentAtLength
 * @brief Splits a segment at a specified length from an endpoint
 *
 * This class uses an the nonlinear least squares algorithm of the
 * AlgLib library to split a parametric segments at a specified length
 * from one of its endpoints.
 */
class SplitSegmentAtLength : public NonlinearLeastSquares {
public:
    /** Constructor for SplitSegmentAtLength */
    SplitSegmentAtLength()
            : NonlinearLeastSquares(split_segment_at_length_residual) {
        m_num_params = 1;
        m_num_resid = 1;
    }

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
    virtual void initialize(void* ptr) final;

    /**
     * This function is intended to postprocess the parameters,
     * converting them from alglib::read_1d_array to vector type.
     *
     * @param ptr The same pointer passed to run
     * @return A vector of parameters
     */
    virtual std::vector<double> finalize(void* ptr) final;

    /** Set initial parameter */
    void set_initial_parameter(double t) { m_t = t; }

    /** Set the parametric */
    void set_parametric(parametric::IParametric* p) { m_param = p; }

    /** Set the endpoint used */
    void set_point(Real2 p) { m_point = p; }

    /** Set the endpoint used */
    void set_point(mesh::Point* p) { m_point = {p->x, p->y}; }

    /** Set the distance from the initial point */
    void set_length(double len) { m_length = len; }

    /** Get the initial parameter */
    double get_initial_parameter() { return m_t; }

    /** Get the parametric */
    const parametric::IParametric* const get_parametric() { return m_param; }

    /** Get the endpoint used */
    const Real2& get_point() { return m_point; }

    /** Get the length */
    const double get_length() { return m_length; }

private:
    parametric::IParametric* m_param;
    Real2 m_point;
    double m_length;
    double m_t;
};


} // namespace optimize

} // namespace math

} // namespace umr

#endif // OPTIMIZERS_H

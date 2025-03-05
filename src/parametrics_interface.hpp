#ifndef PARAMETERICS_INTERFACE_H
#define PARAMETERICS_INTERFACE_H

#include <memory>
#include <functional>


/** Represents a point in 2d space */
using Real2 = std::pair<double, double>;


namespace umr {

namespace parametric {

// f : t -> (x, y)
/** A type for function used in the IParametric framework */
using ParametricFunction = std::function<Real2(double)>;


/**
 * @class IParametric
 * @brief An interface to handle parametric functions
 *
 * This class is designed to handle parametric functions \f$
 * \mathbf{x} = f(t) \f$ where \f$ f : R \rightarrow R^2 \f$. The
 * values of the parameter \f$ t \f$ are considered
 * bounded. Evaluating the parameteric at the lower bound returns and
 * initial point, while evaluating at the upper bound returns a final
 * point. This is done to represent an input segment in a mesh.
 *
 * This class is also responsible for splitting parametrics. Its
 * purpose is to facilitate splitting input segments in the mesh.
 */
class IParametric {
public:
    /** Constructor for IParametric */
    IParametric();

    /** Constructor for IParametric */
    IParametric(const double tmin, const double tmax);

    /** Constructor for IParametric */
    IParametric(IParametric& ip);

    /** Constructor for IParametric */
    IParametric(IParametric& ip, const double tmin, const double tmax);

    /** Destructor for IParametric */
    virtual ~IParametric() {};

    /** Copy IParametric */
    virtual IParametric* clone() = 0;

    /**
     * Split IParametric
     *
     * @param tmin Initial parameter value
     * @param tmax Final parameter value
     * @return A new IParametric object from tmin to tmax
     */
    virtual IParametric* split(double tmin, double tmax) = 0;

    /**
     * Evaluate IParametric at a given parameter within the
     * parametric's domain
     *
     * @param t Parameter value
     * @return An (x, y) pair
     */
    Real2 evaluate(const double t) const;

    /**
     * Evaluate IParametric at a given parameter
     *
     * @param t Parameter value
     * @return An (x, y) pair
     */
    virtual Real2 evaluate_model(const double t) const = 0;

    /**
     * Evaluate a numerical derivative at a given parameter within the
     * parametric's domain
     *
     * @param t Parameter value
     * @return An (x, y) pair
     */
    virtual Real2 differentiate_model(double t) const;

    /**
     * Calculate a line integral of the parametric
     *
     * @return Line integral of the parametric
     */
    virtual double get_length() = 0;

    /**
     * Evaluate at the lower bound of the parameter
     *
     * @return An (x, y) pair
     */
    Real2 evaluate_tmin() const { return {m_x0_b, m_y0_b}; };

    /**
     * Evaluate at the lower bound of the parameter
     *
     * @return An (x, y) pair
     */
    Real2 evaluate_tmax() const { return {m_x1_b, m_y1_b}; };

    /** Get the lower bound of the parameter */
    double get_tmin() const { return m_parameters.front(); }

    /** Get the upper bound of the parameter */
    double get_tmax() const { return m_parameters.back(); }

    /** Get the number of subsegments for splitting */
    size_t get_num_subsegments() { return m_num_subsegments; }

    /** Set the number of subsegments for splitting */
    void set_num_subsegments(size_t n);

    /** Get the subsegments after splitting */
    std::vector<std::shared_ptr<IParametric>> get_subsegments();

    void set_parameters(std::vector<double>& params);
    void set_parameters(std::vector<double>&& params);

    /**
     * Set the bounds of the parameter
     *
     * @param tmin Lower bound
     * @param tmax Upper bound
     */
    void set_parameter_bounds(double tmin, double tmax);

    /**
     * Get a vector of parameters with a size equal to the number of
     * subsegments + 1. Defaults to a size of 2 representing the
     * parameter's lower and upper bounds.  It changes when running
     * set_num_subsegments. Parameters will be correct after running
     * optimize_parameters.
     */
    std::vector<double>& get_parameters() { return m_parameters; }

    /** Sets the initial point corresponding to the minimum parameter */
    void set_initial_point(double x, double y);

    /** Sets the final point corresponding to the maximum parameter */
    void set_final_point(double x, double y);

    /**
     * Finds parameters bounds for the parameteric subsegments,
     * splitting into subsegments of uniform length.
     */
    void optimize_parameters();

    /**
     * Finds a point in parameter space such that the Euclidean
     * distance between the point and the initial point equals this
     * function's input length.
     */
    Real2 point_from_start(double length);

    /**
     * Finds a point in parameter space such that the Euclidean
     * distance between the point and the final point equals this
     * function's input length.
     */
    Real2 point_from_end(double length);

protected:

    /** Confirms the input value is within the parameter's bounds */
    bool validate_parameter_in_bounds(double t) const;

    /**
     * Initializes the parameteric, overriding input points for
     * precision purposes if necessary (e.g. periodic function)
     */
    void initialize();

    /**
     * Initializes the parameteric, overriding input points for
     * precision purposes if necessary (e.g. matching final and
     * initial points of consecutive subsegments).
     */
    void initialize(IParametric& ip);

    std::vector<double> m_parameters;

    // Equals model evaluate at current parametric parameters tmin and tmax
    double m_x0_b, m_y0_b, m_x1_b, m_y1_b;
    size_t m_num_subsegments = 1;
};


} // namespace parametric

} // namespace umr

#endif // PARAMETERICS_INTERFACE_H

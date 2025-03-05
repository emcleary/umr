#ifndef DENSITY_H
#define DENSITY_H

#include <functional>
#include <vector>
#include <string>

#include "../ext/alglib-cpp/src/interpolation.h"

#include "utilities.hpp"


namespace umr {

namespace density {


/** A function mapping an (x, y) pair to lengthscale */
using DensityFunction = std::function<double(double, double)>;


/**
 * @class DensityInterpolator
 * @brief Calculates density using interpolation
 *
 * This class load a data file and calculates density by
 * interpolating. Interpolation uses the RBF model from AlgLib.
 * Data files must be of CSV format, each line being \f$ x, y, \rho \f$
 * where \f$ x \f$ and \f$ y \f$ are coordinates of the point,
 * and \f$ \rho \f$ is the density at that point.
 */
class DensityInterpolator {
public:

    /** Constructor for DensityInterpolator */
    DensityInterpolator(std::string& filename);

    /** Constructor for DensityInterpolator */
    DensityInterpolator(std::string&& filename);

    /** Calculates density at point x, y */
    double interpolate(double x, double y);

    /** Gets the filename of the interpolation data used */
    std::string& get_filename() const;

private:
    /** Loads the data and builds a model */
    void initialize();

    /** Loads the CSV file provided at construction */
    void load_data();

    /** Builds the AlgLib model used for interpolation */
    void build_model();

    std::string& m_filename;
    alglib::real_2d_array m_data;
    alglib::rbfmodel m_model;
    alglib::rbfreport m_rep;
};


/**
 * @class DensityManager
 * @brief Manages density functions and interpolators
 *
 * This class manages all density functions and interpolators added to
 * it. When calculating density, it returns the minimum calculated
 * from all functions and interpolators present.
 */
class DensityManager {
public:

    /** Constructor for DensityManager */
    DensityManager() {}

    /** Destructor for DensityManager */
    ~DensityManager();

    /**
     * Sets a minimum value for the density. It is useful when
     * calculating densities to ensure that no function or
     * interpolator ever returns a negavite value. Setting it is
     * optional.
     */
    void set_minimum_density(double d);

    /** Sets no minimum value for the density. */
    void set_no_minimum_density();

    /**
     * Returns the minimum value for the density. If none are set it
     * will return the max value of double.
     */
    double get_minimum_density() const;

    /** Add an arbitrary function for calculating density */
    void add_functor(DensityFunction function);

    /** Build and add an interpolator calculating density */
    void add_interpolator(std::string& filename);

    /** Build and add an interpolator calculating density */
    void add_interpolator(std::string&& filename);

    /** Evaluates the density at point (x, y) */
    double evaluate(double x, double y) const;

private:
    bool m_min_set = false;
    double m_min_density = std::numeric_limits<double>::max();
    std::vector<DensityInterpolator*> m_interpolators;
    std::vector<DensityFunction> m_functors;
};


} // density

} // umr

#endif // DENSITY_H

#ifndef PARAMETRICS_GENERAL_H
#define PARAMETRICS_GENERAL_H

#include "parametrics_interface.hpp"


namespace umr {

namespace parametric {


/**
 * @class Parametric
 * @brief A general parameteric function
 *
 * This parameteric function is supplied by the user.
 * The function must take in a double prepresenting the parameter
 * and return a pair of doubles representing a point at the coordinate
 * (x, y). For simplicity, use the type ParametricFunction.
 */
class Parametric : public IParametric {
public:

    /**
     * Constructor for the parametric. Without supplying an additional
     * function for the derivative, the function's derivative will be
     * calculated numerically as needed.
     *
     * @param function The input function
     * @param tmin Lower bound of the parameter
     * @param tmax Upper bound of the parameter
     */
    Parametric(ParametricFunction function, const double tmin, const double tmax);

    /**
     * Constructor for the parametric
     *
     * @param function The input function
     * @param derivative The input function's derivative
     * @param tmin Lower bound of the parameter
     * @param tmax Upper bound of the parameter
     */
    Parametric(ParametricFunction function, ParametricFunction derivative,
            const double tmin, const double tmax);

    /** Constructor for Parametric */
    Parametric(Parametric& parametric);

    /** Constructor for Parametric */
    Parametric(Parametric& parametric, double tmin, double tmax);

    /** Copy Parametric */
    virtual IParametric* clone() override { return new Parametric(*this); }

    /**
     * Split Parametric
     *
     * @param tmin Initial parameter value
     * @param tmax Final parameter value
     * @return A new IParametric object from tmin to tmax
     */
    virtual IParametric* split(double tmin, double tmax)
        override { return new Parametric(*this, tmin, tmax); }

    /**
     * Evaluate IParametric at a given parameter
     *
     * @param t Parameter value
     * @return An (x, y) pair
     */
    virtual Real2 evaluate_model(const double t) const;

    /**
     * Evaluate a derivative at a given parameter within the
     * parametric's domain
     *
     * @param t Parameter value
     * @return An (x, y) pair
     */
    virtual Real2 differentiate_model(double t) const;

    /**
     * Calculate the line integral of the function from the
     * parameter's lower bound to its upper bound.
     *
     * @return Line integral of the Parametric
     */
    virtual double get_length();

private:
    ParametricFunction m_function;
    ParametricFunction m_derivative;
};


} // namespace parametric

} // namespace umr

#endif // PARAMETRICS_GENERAL_H

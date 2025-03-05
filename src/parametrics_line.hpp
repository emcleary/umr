#ifndef PARAMETRICS_LINE_H
#define PARAMETRICS_LINE_H

#include "parametrics_interface.hpp"

#include "point.hpp"


namespace umr {

namespace parametric {


/**
 * @class Line
 * @brief A parameteric function representing an line
 *
 * This parameteric function is for an line.
 * The model used takes the form
 * \f[ x = x_0 + m * t \f]
 * \f[ y = y_0 + m * t \f]
 *
 * where \f$ (x_0, y_0) \f$ represents the initial point and
 * \f$ m \f$ represents the lines slope. Bounds of parameter \f$ t \f$
 * must be between 0 and 1.
 */
class Line : public IParametric {
public:
    /**
     * Constructor for the parametric
     *
     * @param p0 The initial point
     * @param p1 The final point
     */
    Line(const mesh::Point& p0, const mesh::Point& p1);

    /**
     * Constructor for the parametric. Parameter bounds must be in \f$
     * [0, 1] \f$.
     *
     * @param p0 The initial point
     * @param p1 The final point
     * @param tmin Lower bound of the parameter
     * @param tmax Upper bound of the parameter
     */
    Line(const mesh::Point& p0, const mesh::Point& p1, const double tmin, const double tmax);

    /** Constructor for Line */
    Line(Line& line);

    /** Constructor for Line */
    Line(Line& line, double tmin, double tmax);

    /** Copy Line */
    virtual IParametric* clone() override { return new Line(*this); }

    /**
     * Split Line
     *
     * @param tmin Initial parameter value
     * @param tmax Final parameter value
     * @return A new IParametric object from tmin to tmax
     */
    virtual IParametric* split(double tmin, double tmax)
        override { return new Line(*this, tmin, tmax); }

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
     * Calculate the length of the Line from the parameter's lower
     * bound to its upper bound.
     *
     * @return Line integral of the Line
     */
    virtual double get_length();

private:
    double m_x0, m_y0, m_x1, m_y1;
};


} // namespace parametric

} // namespace umr

#endif // PARAMETRICS_LINE_H

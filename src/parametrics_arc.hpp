#ifndef PARAMETRICS_ARC_H
#define PARAMETRICS_ARC_H

#include "parametrics_interface.hpp"

#include "point.hpp"


namespace umr {

namespace parametric {

/**
 * @class Arc
 * @brief A parameteric function representing an arc
 *
 * This parameteric function is for an arc, i.e. a subsegment of a circle.
 * The model used takes the form
 * \f[ x = x_c + r * \cos(t) \f]
 * \f[ y = y_c + r * \sin(t) \f]
 *
 * where \f$ (x_c, y_c) \f$ represents the center point of the circle and
 * \f$ r \f$ represents the circle's radius. Bounds of parameter \f$ t \f$
 * default to 0 and \f$ 2\pi \f$.
 */
class Arc : public IParametric {
public:
    /**
     * Constructor for the parametric
     *
     * @param pc A center point
     * @param radius The radius of the arc
     * @param tmin Lower bound of the parameter
     * @param tmax Upper bound of the parameter
     */
    Arc(const mesh::Point& pc, const double radius, const double tmin, const double tmax);

    /** Constructor for Arc */
    Arc(Arc& arc);

    /** Constructor for Arc */
    Arc(Arc& arc, double tmin, double tmax);

    /** Copy Arc */
    virtual IParametric* clone() override { return new Arc(*this); }

    /**
     * Split Arc
     *
     * @param tmin Initial parameter value
     * @param tmax Final parameter value
     * @return A new IParametric object from tmin to tmax
     */
    virtual IParametric* split(double tmin, double tmax)
        override { return new Arc(*this, tmin, tmax); }

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
     * Calculate a line integral of the Arc from the parameter's lower
     * bound to its upper bound.
     *
     * @return Line integral of the Arc
     */
    virtual double get_length();

private:
    const double m_xc, m_yc, m_radius;
};


} // namespace parametric

} // namespace umr

#endif // PARAMETRICS_ARC_H

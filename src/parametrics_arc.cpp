#include "parametrics_arc.hpp"

#include <iostream>

#include "inequalities.hpp"


namespace umr {

namespace parametric {


Arc::Arc(Arc& arc) : IParametric(arc), m_xc(arc.m_xc), m_yc(arc.m_yc),
                     m_radius(arc.m_radius) {}


Arc::Arc(Arc& arc, double tmin, double tmax)
        : IParametric(arc, tmin, tmax),
          m_xc(arc.m_xc), m_yc(arc.m_yc), m_radius(arc.m_radius) {
    initialize(arc);
}


Arc::Arc(const mesh::Point& pc, const double radius, const double tmin,
        const double tmax)
        : IParametric(tmin, tmax),
          m_xc(pc.x), m_yc(pc.y), m_radius(radius) {
    double dtheta = tmax - tmin;
    if (inequalities::is_gt(dtheta, 2 * std::numbers::pi)) {
        std::cerr << "Arc::Arc: parameters must have a range less than 2*PI\n";
        exit(-1);
    }
    initialize();
    set_num_subsegments(3);
}


Real2 Arc::evaluate_model(const double t) const {
    double x = m_xc + m_radius * cos(t);
    double y = m_yc + m_radius * sin(t);
    return {x, y};
}


Real2 Arc::differentiate_model(double t) const {
    double dxdt = m_radius * sin(t);
    double dydt = m_radius * cos(t);
    return {dxdt, dydt};
}


double Arc::get_length() {
    double dtheta = get_tmax() - get_tmin();
    return m_radius * dtheta;
}


} // namespace parametric

} // namespace umr

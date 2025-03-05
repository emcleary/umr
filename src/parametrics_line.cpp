#include "parametrics_line.hpp"

#include <cassert>
#include <iostream>

#include "inequalities.hpp"


namespace umr {

namespace parametric {


Line::Line(const mesh::Point& p0, const mesh::Point& p1)
        : IParametric(),
          m_x0(p0.x), m_y0(p0.y),
          m_x1(p1.x), m_y1(p1.y) {
    initialize();
}


Line::Line(const mesh::Point& p0, const mesh::Point& p1,
        const double tmin, const double tmax)
        : IParametric(tmin, tmax),
          m_x0(p0.x), m_y0(p0.y),
          m_x1(p1.x), m_y1(p1.y) {
    if (inequalities::is_lt(tmin, 0) || inequalities::is_gt(tmax, 1)
            || inequalities::is_ge(tmin, tmax)) {
        std::cerr << "ERROR: Line::Line initial tmin, tmax error; "
            "need 0 <= tmin < tmax <= 1\n";
        exit(-1);
    }
    initialize();
}


Line::Line(Line& line) : IParametric(line),
			 m_x0(line.m_x0), m_y0(line.m_y0),
			 m_x1(line.m_x1), m_y1(line.m_y1) {}


Line::Line(Line& line, double tmin, double tmax)
        : IParametric(line, tmin, tmax),
          m_x0(line.m_x0), m_y0(line.m_y0),
          m_x1(line.m_x1), m_y1(line.m_y1) {
    assert(inequalities::is_ge(tmin, 0));
    assert(inequalities::is_le(tmax, 1));
    assert(inequalities::is_lt(tmin, tmax));
    initialize(line);
}


Real2 Line::evaluate_model(const double t) const {
    double dx = m_x1 - m_x0;
    double dy = m_y1 - m_y0;
    double x = m_x0 + dx * t;
    double y = m_y0 + dy * t;
    return {x, y};
}


Real2 Line::differentiate_model(double t) const {
    double dxdt = m_x1 - m_x0;
    double dydt = m_y1 - m_y0;
    return {dxdt, dydt};
}


double Line::get_length() {
    double dx = m_x1_b - m_x0_b;
    double dy = m_y1_b - m_y0_b;
    return sqrt(dx * dx + dy * dy);
}


} // namespace parametric

} // namespace umr

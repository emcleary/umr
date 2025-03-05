#include <cmath> // abs
#include <cassert>

#include "triangle.hpp"
#include "utilities.hpp"
#include "math.hpp"


namespace umr {

namespace mesh {


Triangle::Triangle(const Point& p0, const Point& p1, const Point& p2)
        : m_p0(p0), m_p1(p1), m_p2(p2),
          m_center(make_circumcenter(p0, p1, p2)) {

    double length01 = p0.distance_to(p1);
    double length12 = p1.distance_to(p2);
    double length20 = p2.distance_to(p0);
    m_min_length = std::min({length01, length12, length20});
    m_max_length = std::max({length01, length12, length20});

    double angle0 = math::loc_angle(length12, length20, length01);
    double angle1 = math::loc_angle(length20, length01, length12);
    double angle2 = math::loc_angle(length01, length12, length20);

    static const double to_deg = 180 / std::numbers::pi;
    m_min_angle = to_deg * std::min({angle0, angle1, angle2});
    m_max_angle = to_deg * std::max({angle0, angle1, angle2});

    m_radius = m_center.distance_to(p0);
}


Point make_circumcenter(const Point& p0,
        const Point& p1, const Point& p2) {
    auto [xc, yc] = math::calculate_circumcenter(p0, p1, p2);
    return Point(xc, yc);
}


} // namespace mesh

} // namespace umr
